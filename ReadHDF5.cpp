/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

//-------------------------------------------------------------------------
// Filename      : ReadHDF5.cpp
//
// Purpose       : HDF5 Writer 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/18/04
//-------------------------------------------------------------------------

#include <assert.h>
/* include our MPI header before any HDF5 because otherwise
   it will get included indirectly by HDF5 */
#ifdef USE_MPI
#  include "moab_mpi.h"
#  include "moab/ParallelComm.hpp"
#endif 
#include <H5Tpublic.h>
#include <H5Ppublic.h>
#include "moab/Interface.hpp"
#include "Internals.hpp"
#include "MBTagConventions.hpp"
#include "ReadHDF5.hpp"
#include "moab/CN.hpp"
#include "FileOptions.hpp"
#ifdef HDF5_PARALLEL
#  include <H5FDmpi.h>
#  include <H5FDmpio.h>
#endif
//#include "WriteHDF5.hpp"

#include <stdlib.h>
#include <string.h>
#include <limits>
#include <functional>

#include "IODebugTrack.hpp"
#include "ReadHDF5Dataset.hpp"

namespace moab {

#undef BLOCKED_COORD_IO

#define READ_HDF5_BUFFER_SIZE (40*1024*1024)

#define assert_range( PTR, CNT ) \
  assert( (PTR) >= (void*)dataBuffer ); assert( ((PTR)+(CNT)) <= (void*)(dataBuffer + bufferSize) );


// This function doesn't do anything useful.  It's just a nice
// place to set a break point to determine why the reader fails.
static inline ErrorCode error( ErrorCode rval )
  { return rval; }

static void copy_sorted_file_ids( const EntityHandle* sorted_ids, 
                                  long num_ids,
                                  Range& results )
{
  Range::iterator hint = results.begin();
  long i = 0;
  while (i < num_ids) {
    EntityHandle start = sorted_ids[i];
    for (++i; i < num_ids && sorted_ids[i] == 1+sorted_ids[i-1]; ++i);
    hint = results.insert( hint, start, sorted_ids[i-1] );
  }
}

static void intersect( const mhdf_EntDesc& group, const Range& range, Range& result )
{
  Range::const_iterator s, e;
  s = Range::lower_bound( range.begin(), range.end(), group.start_id );
  e = Range::lower_bound( s, range.end(), group.start_id + group.count );
  result.merge( s, e );
}

#define debug_barrier() debug_barrier_line(__LINE__)
void ReadHDF5::debug_barrier_line(int lineno)
{
#ifdef USE_MPI
  const unsigned threshold = 2;
  static unsigned long count = 0;
  if (dbgOut.get_verbosity() >= threshold) {
    dbgOut.printf( threshold, "*********** Debug Barrier %lu (@%d)***********\n", ++count, lineno);
    MPI_Barrier( myPcomm->proc_config().proc_comm() );
  }
#endif
}

ReaderIface* ReadHDF5::factory( Interface* iface )
  { return new ReadHDF5( iface ); }

ReadHDF5::ReadHDF5( Interface* iface )
  : bufferSize( READ_HDF5_BUFFER_SIZE ),
    dataBuffer( 0 ),
    iFace( iface ), 
    filePtr( 0 ), 
    fileInfo( 0 ), 
    readUtil( 0 ),
    handleType( 0 ),
    indepIO( H5P_DEFAULT ),
    collIO( H5P_DEFAULT ),
    debugTrack( false ),
    dbgOut(stderr)
{
}

ErrorCode ReadHDF5::init()
{
  ErrorCode rval;

  if (readUtil) 
    return MB_SUCCESS;
  
  indepIO = collIO = H5P_DEFAULT;
  //WriteHDF5::register_known_tag_types( iFace );
  
  handleType = H5Tcopy( H5T_NATIVE_ULONG );
  if (handleType < 0)
    return error(MB_FAILURE);
  
  if (H5Tset_size( handleType, sizeof(EntityHandle)) < 0)
  {
    H5Tclose( handleType );
    return error(MB_FAILURE);
  }
  
  void* ptr = 0;
  rval = iFace->query_interface( "ReadUtilIface", &ptr );
  if (MB_SUCCESS != rval)
  {
    H5Tclose( handleType );
    return error(rval);
  }
  readUtil = reinterpret_cast<ReadUtilIface*>(ptr);
  
  idMap.clear();
  fileInfo = 0;
  debugTrack = false;
  myPcomm = 0;
  
  return MB_SUCCESS;
}
  

ReadHDF5::~ReadHDF5()
{
  if (!readUtil) // init() failed.
    return;

  iFace->release_interface( "ReadUtilIface", readUtil );
  H5Tclose( handleType );
}

ErrorCode ReadHDF5::set_up_read( const char* filename,
                                 const FileOptions& opts )
{
  ErrorCode rval;
  mhdf_Status status;
  indepIO = collIO = H5P_DEFAULT;
  mpiComm = 0;

  if (MB_SUCCESS != init())
    return error(MB_FAILURE);
  
  // Set up debug output
  int tmpval;
  if (MB_SUCCESS == opts.get_int_option("DEBUG_IO", 1, tmpval)) {
    dbgOut.set_verbosity(tmpval);
    dbgOut.set_prefix("H5M ");
  }
  
  // Enable some extra checks for reads.  Note: amongst other things this
  // will print errors if the entire file is not read, so if doing a 
  // partial read that is not a parallel read, this should be disabled.
  debugTrack = (MB_SUCCESS == opts.get_null_option("DEBUG_BINIO"));
    
    // Handle parallel options
  std::string junk;
  bool use_mpio = (MB_SUCCESS == opts.get_null_option("USE_MPIO"));
  rval = opts.match_option("PARALLEL", "READ_PART");
  bool parallel = (rval != MB_ENTITY_NOT_FOUND);
  nativeParallel = (rval == MB_SUCCESS);
  if (use_mpio && !parallel) {
    readUtil->report_error( "'USE_MPIO' option specified w/out 'PARALLEL' option" );
    return MB_NOT_IMPLEMENTED;
  }
  
  // This option is intended for testing purposes only, and thus
  // is not documented anywhere.  Decreasing the buffer size can
  // expose bugs that would otherwise only be seen when reading
  // very large files.
  rval = opts.get_int_option( "BUFFER_SIZE", bufferSize );
  if (MB_SUCCESS != rval) {
    bufferSize = READ_HDF5_BUFFER_SIZE;
  }
  else if (bufferSize < (int)std::max( sizeof(EntityHandle), sizeof(void*) )) {
    return error(MB_INVALID_SIZE);
  }
  
  dataBuffer = (char*)malloc( bufferSize );
  if (!dataBuffer)
    return error(MB_MEMORY_ALLOCATION_FAILED);
  
  if (use_mpio || nativeParallel) {
#ifndef HDF5_PARALLEL
    readUtil->report_error("MOAB not configured with parallel HDF5 support");
    free(dataBuffer);
    return MB_NOT_IMPLEMENTED;
#else
    int pcomm_no = 0;
    rval = opts.get_int_option("PARALLEL_COMM", pcomm_no);
    if (rval == MB_TYPE_OUT_OF_RANGE) {
      readUtil->report_error("Invalid value for PARALLEL_COMM option");
      return rval;
    }
    myPcomm = ParallelComm::get_pcomm(iFace, pcomm_no);
    if (0 == myPcomm) {
      myPcomm = new ParallelComm(iFace);
    }
    const int rank = myPcomm->proc_config().proc_rank();
    dbgOut.set_rank(rank);
    mpiComm = new MPI_Comm(myPcomm->proc_config().proc_comm());

      // Open the file in serial on root to read summary
    dbgOut.tprint( 1, "Getting file summary\n" );
    fileInfo = 0;
    unsigned long size = 0;
    if (rank == 0) {
      filePtr = mhdf_openFile( filename, 0, NULL, &status );
      if (filePtr) {  
        fileInfo = mhdf_getFileSummary( filePtr, handleType, &status );
        if (!is_error(status)) {
          size = fileInfo->total_size;
          fileInfo->offset = (unsigned char*)fileInfo;
        }
      }
      mhdf_closeFile( filePtr, &status );
      if (fileInfo && mhdf_isError(&status)) {
        free(fileInfo);
        fileInfo = 0;
      }
    }
      // Broadcast the size of the struct (zero indicates an error)
    dbgOut.tprint( 1, "Communicating file summary\n" );
    int err = MPI_Bcast( &size, 1, MPI_UNSIGNED_LONG, 0, myPcomm->proc_config().proc_comm() );
    if (err || !size)
      return MB_FAILURE;
    
      // allocate structure
    if (rank != 0) 
      fileInfo = reinterpret_cast<mhdf_FileDesc*>( malloc( size ) );
      // bcast file summary
    MPI_Bcast( fileInfo, size, MPI_BYTE, 0, myPcomm->proc_config().proc_comm() );
      // fix up internal pointers in file summary struct
    if (rank != 0)
      mhdf_fixFileDesc( fileInfo, reinterpret_cast<mhdf_FileDesc*>(fileInfo->offset) );
  
      // configure HDF5 properties  
    hid_t file_prop = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(file_prop, myPcomm->proc_config().proc_comm(), MPI_INFO_NULL);
    collIO = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(collIO, H5FD_MPIO_COLLECTIVE);
    indepIO = nativeParallel ? H5P_DEFAULT : collIO;

      // re-open file in parallel
    dbgOut.tprintf( 1, "Re-opening \"%s\" for parallel IO\n", filename );
    filePtr = mhdf_openFileWithOpt( filename, 0, NULL, file_prop, &status );
    H5Pclose( file_prop );
    if (!filePtr)
    {
      readUtil->report_error( mhdf_message( &status ));
      free( dataBuffer );
      H5Pclose( indepIO ); 
      if (collIO != indepIO)
        H5Pclose( collIO );
      collIO = indepIO = H5P_DEFAULT;
      return error(MB_FAILURE);
    }
#endif
  }
  else {
  
      // Open the file
    filePtr = mhdf_openFile( filename, 0, NULL, &status );
    if (!filePtr)
    {
      readUtil->report_error( "%s", mhdf_message( &status ));
      free( dataBuffer );
      return error(MB_FAILURE);
    }

      // get file info
    fileInfo = mhdf_getFileSummary( filePtr, handleType, &status );
    if (is_error(status)) {
      free( dataBuffer );
      mhdf_closeFile( filePtr, &status );
      return error(MB_FAILURE);
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::clean_up_read( const FileOptions& )
{
  free( dataBuffer );
  free( fileInfo );
  delete mpiComm;
  mpiComm = 0;

  if (indepIO != H5P_DEFAULT)
    H5Pclose( indepIO );
  if (collIO != indepIO)
    H5Pclose( collIO );
  collIO = indepIO = H5P_DEFAULT;

  mhdf_Status status;
  mhdf_closeFile( filePtr, &status );
  filePtr = 0;
  return is_error(status) ? MB_FAILURE : MB_SUCCESS;
}

ErrorCode ReadHDF5::load_file( const char* filename, 
                               const EntityHandle* file_set, 
                               const FileOptions& opts,
                               const ReaderIface::SubsetList* subset_list,
                               const Tag* file_id_tag )
{
  ErrorCode rval;
 
  rval = set_up_read( filename, opts );
  if (MB_SUCCESS != rval)
    return rval;
 
  if (subset_list) 
    rval = load_file_partial( subset_list->tag_list, 
                              subset_list->tag_list_length, 
                              subset_list->num_parts,
                              subset_list->part_number,
                              opts );
  else
    rval = load_file_impl( opts );
    
  if (MB_SUCCESS == rval && file_id_tag) {
    dbgOut.tprint( 1, "Storing file IDs in tag\n" );
    rval = store_file_ids( *file_id_tag );
  }
  
  if (MB_SUCCESS == rval && 0 != file_set) {
    dbgOut.tprint( 1, "Reading QA records\n" );
    rval = read_qa( *file_set );
  }
  
  
  dbgOut.tprint( 1, "Cleaining up\n" );
  ErrorCode rval2 = clean_up_read( opts );
  if (rval == MB_SUCCESS && rval2 != MB_SUCCESS)
    rval = rval2;
  
  dbgOut.tprint(1, "Read finished.\n");
  
  if (H5P_DEFAULT != collIO)
    H5Pclose( collIO );
  if (H5P_DEFAULT != indepIO)
    H5Pclose( indepIO );
  collIO = indepIO = H5P_DEFAULT;
  
  return rval;
}
  


ErrorCode ReadHDF5::load_file_impl( const FileOptions& )
{
  ErrorCode rval;
  mhdf_Status status;
  std::string tagname;
  int i;

  dbgOut.tprint(1, "Reading all nodes...\n");
  Range ids;
  if (fileInfo->nodes.count) {
    ids.insert( fileInfo->nodes.start_id,
                fileInfo->nodes.start_id + fileInfo->nodes.count - 1);
    rval = read_nodes( ids );
    if (MB_SUCCESS != rval)
      return error(rval);
  }


  dbgOut.tprint(1, "Reading all element connectivity...\n");
  std::vector<int> polyhedra; // need to do these last so that faces are loaded
  for (i = 0; i < fileInfo->num_elem_desc; ++i) {
    if (CN::EntityTypeFromName(fileInfo->elems[i].type) == MBPOLYHEDRON) {
      polyhedra.push_back(i);
      continue;
    }
    
    rval = read_elems( i );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  for (std::vector<int>::iterator it = polyhedra.begin();
       it != polyhedra.end(); ++it) {
    rval = read_elems( *it );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
  dbgOut.tprint(1, "Reading all sets...\n");
  ids.clear();
  if (fileInfo->sets.count) {
    ids.insert( fileInfo->sets.start_id,
                fileInfo->sets.start_id + fileInfo->sets.count - 1);
    rval = read_sets( ids );
    if (rval != MB_SUCCESS) {
      return error(rval);
    }
  }
  
  dbgOut.tprint(1, "Reading all adjacencies...\n");
  for (i = 0; i < fileInfo->num_elem_desc; ++i) {
    if (!fileInfo->elems[i].have_adj)
      continue;
    
    long table_len;
    hid_t table = mhdf_openAdjacency( filePtr, 
                                      fileInfo->elems[i].handle,
                                      &table_len,
                                      &status );
    if (is_error(status))
      return error(MB_FAILURE);
      
    rval = read_adjacencies( table, table_len );
    mhdf_closeData( filePtr, table, &status );
    if (MB_SUCCESS != rval)
      return error(rval);
    if (is_error(status))
      return error(MB_FAILURE);
  }

  dbgOut.tprint(1, "Reading all tags...\n");
  for (i = 0; i < fileInfo->num_tag_desc; ++i) {
    rval = read_tag( i );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
  dbgOut.tprint(1, "Core read finished.  Cleaning up...\n");
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::find_int_tag( const char* name, int& index )
{
  for (index = 0; index < fileInfo->num_tag_desc; ++index) 
    if (!strcmp( name, fileInfo->tags[index].name))
      break;

  if (index == fileInfo->num_tag_desc) {
    readUtil->report_error( "File does not contain subset tag '%s'", name );
    return error(MB_TAG_NOT_FOUND);
  }

  if (fileInfo->tags[index].type != mhdf_INTEGER ||
      fileInfo->tags[index].size != 1) {
    readUtil->report_error( "Tag '%s' does not containa single integer value", name );
    return error(MB_TYPE_OUT_OF_RANGE);
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::get_subset_ids( const ReaderIface::IDTag* subset_list,
                                      int subset_list_length,
                                      Range& file_ids )
{
  ErrorCode rval;
  
  for (int i = 0; i < subset_list_length; ++i) {  
    
    int tag_index;
    rval = find_int_tag( subset_list[i].tag_name, tag_index );
    if (MB_SUCCESS != rval)
      return error(rval);
  
    Range tmp_file_ids;
    if (!subset_list[i].num_tag_values) {
      rval = get_tagged_entities( tag_index, tmp_file_ids );
    }
    else {
      std::vector<int> ids( subset_list[i].tag_values, 
                            subset_list[i].tag_values + subset_list[i].num_tag_values );
      std::sort( ids.begin(), ids.end() );
      rval = search_tag_values( tag_index, ids, tmp_file_ids );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
    
    if (tmp_file_ids.empty())
      return error(MB_ENTITY_NOT_FOUND);
    
    if (i == 0) 
      file_ids.swap( tmp_file_ids );
    else 
      file_ids = intersect( tmp_file_ids, file_ids );
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::get_partition( Range& tmp_file_ids, int num_parts, int part_number )
{    
     // check that the tag only identified sets
   if ((unsigned long)fileInfo->sets.start_id > tmp_file_ids.front()) {
     dbgOut.print(1,"Ignoreing non-set entities with partition set tag\n");
     tmp_file_ids.erase( tmp_file_ids.begin(), 
                         tmp_file_ids.lower_bound( 
                           (EntityHandle)fileInfo->sets.start_id ) );
   }
   unsigned long set_end = (unsigned long)fileInfo->sets.start_id + fileInfo->sets.count;
   if (tmp_file_ids.back() >= set_end) {
     dbgOut.print(1,"Ignoreing non-set entities with partition set tag\n");
     tmp_file_ids.erase( tmp_file_ids.upper_bound( (EntityHandle)set_end ),
                         tmp_file_ids.end() );
   }
      
  Range::iterator s = tmp_file_ids.begin();
  size_t num_per_proc = tmp_file_ids.size() / num_parts;
  size_t num_extra = tmp_file_ids.size() % num_parts;
  Range::iterator e;
  if (part_number < (long)num_extra) {
    s += (num_per_proc+1) * part_number;
    e = s;
    e += (num_per_proc+1);
  }
  else {
    s += num_per_proc * part_number + num_extra;
    e = s;
    e += num_per_proc;
  }
  tmp_file_ids.erase(e, tmp_file_ids.end());
  tmp_file_ids.erase(tmp_file_ids.begin(), s);

  return MB_SUCCESS;
}


ErrorCode ReadHDF5::load_file_partial( const ReaderIface::IDTag* subset_list,
                                       int subset_list_length,
                                       int num_parts,
                                       int part_number,
                                       const FileOptions& opts )
{
  mhdf_Status status;
  
  for (int i = 0; i < subset_list_length; ++i) {
    dbgOut.printf( 1, "Select by \"%s\" with num_tag_values = %d\n",
                   subset_list[i].tag_name, subset_list[i].num_tag_values );
    if (subset_list[i].num_tag_values) {
      assert(0 != subset_list[i].tag_values);
      dbgOut.printf( 1, "  \"%s\" values = { %d",
        subset_list[i].tag_name, subset_list[i].tag_values[0] );
      for (int j = 1; j < subset_list[i].num_tag_values; ++j)
        dbgOut.printf( 1, ", %d", subset_list[i].tag_values[j] );
      dbgOut.printf(1," }\n");
    }
  }
  if (num_parts) 
    dbgOut.printf( 1, "Partition with num_parts = %d and part_number = %d\n", 
                   num_parts, part_number );
  
  dbgOut.tprint( 1, "RETREIVING TAGGED ENTITIES\n" );
    
  Range file_ids;
  ErrorCode rval = get_subset_ids( subset_list, subset_list_length, file_ids );
  if (MB_SUCCESS != rval)
    return error(rval);
  
  if (num_parts) {
    rval = get_partition( file_ids, num_parts, part_number );
    if (MB_SUCCESS != rval)
      return error(rval);
  }

  dbgOut.print_ints( 3, "Set file IDs for partial read: ", file_ids );
  
  dbgOut.tprint( 1, "GATHERING ADDITIONAL ENTITIES\n" );
  
  const char* const set_opts[] = { "NONE", "SETS", "CONTENTS" };
  int child_mode;
  rval = opts.match_option( "CHILDREN", set_opts, child_mode );
  if (MB_ENTITY_NOT_FOUND == rval)
    child_mode = 2;
  else if (MB_SUCCESS != rval) {
    readUtil->report_error( "Invalid value for 'CHILDREN' option" );
    return error(rval);
  }
  int content_mode;
  rval = opts.match_option( "SETS", set_opts, content_mode );
  if (MB_ENTITY_NOT_FOUND == rval)
    content_mode = 2;
  else if (MB_SUCCESS != rval) {
    readUtil->report_error( "Invalid value for 'SETS' option" );
    return error(rval);
  }
  
    // If we want the contents of contained/child sets, 
    // search for them now (before gathering the non-set contents
    // of the sets.)
  Range sets;
  intersect( fileInfo->sets, file_ids, sets );
  if (content_mode == 2 || child_mode == 2) {
    rval = read_set_ids_recursive( sets, content_mode == 2, child_mode == 2 );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
    // get elements and vertices contained in sets
  rval = get_set_contents( sets, file_ids );
  if (MB_SUCCESS != rval)
    return error(rval);

  dbgOut.print_ints( 3, "File IDs for partial read: ", file_ids );
  debug_barrier();
    
  dbgOut.tprint( 1, "GATHERING NODE IDS\n" );
  
    // if input contained any polyhedra, need to get faces
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    EntityType type = CN::EntityTypeFromName( fileInfo->elems[i].type );
    if (type != MBPOLYHEDRON)
      continue;
    
    debug_barrier();
    dbgOut.print( 2, "    Getting polyhedra faces\n" );
    
    Range polyhedra;
    intersect( fileInfo->elems[i].desc, file_ids, polyhedra );
    rval = read_elems( i, polyhedra, file_ids );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
    // get node file ids for all elements
  Range nodes;
  intersect( fileInfo->nodes, file_ids, nodes );
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    EntityType type = CN::EntityTypeFromName( fileInfo->elems[i].type );
    if (type <= MBVERTEX || type >= MBENTITYSET) {
      assert( false ); // for debug code die for unknown element tyoes
      continue; // for release code, skip unknown element types
    }
    if (MBPOLYHEDRON == type)
      continue;
    
    debug_barrier();
    dbgOut.printf( 2, "    Getting element node IDs for: %s\n", fileInfo->elems[i].handle );
    
    Range subset;
    intersect( fileInfo->elems[i].desc, file_ids, subset );
    rval = read_elems( i, subset, nodes );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
    
  debug_barrier();
  dbgOut.tprint( 1, "READING NODE COORDINATES\n" );
  
    // Read node coordinates and create vertices in MOAB
    // NOTE:  This populates the RangeMap with node file ids,
    //        which is expected by read_node_adj_elems.
  rval = read_nodes( nodes );
  if (MB_SUCCESS != rval)
    return error(rval);

 
  debug_barrier();
  dbgOut.tprint( 1, "READING ELEMENTS\n" );
 
    // decide if we need to read additional elements
  int side_mode;
  const char* const options[] = { "EXPLICIT", "NODES", "SIDES", 0 };
  rval = opts.match_option( "ELEMENTS", options, side_mode );
  if (MB_ENTITY_NOT_FOUND == rval) {
      // Chose default based on whether or not any elements have been
      // specified.
    Range tmp;
    intersect( fileInfo->nodes, file_ids, tmp );
      // If only nodes were specified, then default to "NODES", otherwise
      // default to "SIDES".
    if (!tmp.empty() && (tmp.size() + sets.size()) == file_ids.size())
      side_mode = 1;
    else
      side_mode = 2;
  }
  else if (MB_SUCCESS != rval) {
    readUtil->report_error( "Invalid value for 'ELEMENTS' option" );
    return error(rval);
  }
  
  switch (side_mode) {
    case 0: // ELEMENTS=EXPLICIT : read only specified element IDS
    for (int dim = 1; dim <= 3; ++dim) {
      for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
        EntityType type = CN::EntityTypeFromName( fileInfo->elems[i].type );
        if (CN::Dimension(type) == dim) {
          debug_barrier();
          dbgOut.tprintf( 2, "    Reading element connectivity for: %s\n", fileInfo->elems[i].handle );
          Range subset;
          intersect( fileInfo->elems[i].desc, file_ids, subset );
          rval = read_elems( fileInfo->elems[i],  subset );
          if (MB_SUCCESS != rval)
            return error(rval);
        }
      }
    }
    break;
    
    case 1: // ELEMENTS=NODES : read all elements for which all nodes have been read
    for (int dim = 1; dim <= 3; ++dim) {
      for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
        EntityType type = CN::EntityTypeFromName( fileInfo->elems[i].type );
        if (CN::Dimension(type) == dim) {
          debug_barrier();
          dbgOut.tprintf( 2, "    Reading node-adjacent elements for: %s\n", fileInfo->elems[i].handle );
          rval = read_node_adj_elems( fileInfo->elems[i] );
          if (MB_SUCCESS != rval)
            return error(rval);
        }
      }
    }
    break;
    
    case 2: // ELEMENTS=SIDES : read explicitly specified elems and any sides of those elems
    debug_barrier();
    dbgOut.tprint( 2, "    Reading elements and sides\n");
    rval = read_elements_and_sides( file_ids );
    if (MB_SUCCESS != rval)
      return error(rval);
    break;
  }
  
  debug_barrier();
  dbgOut.tprint( 1, "READING SETS\n" );
    
    // If reading contained/child sets but not their contents then find
    // them now. If we were also reading their contents we would
    // have found them already.
  if (content_mode == 1 || child_mode == 1) {
    rval = read_set_ids_recursive( sets, content_mode != 0, child_mode != 0 );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
    // Append file IDs of sets containing any of the nodes or elements
    // we've read up to this point.
  rval = find_sets_containing( sets );
  if (MB_SUCCESS != rval)
    return error(rval);
    // Now actually read all set data and instantiate sets in MOAB.
    // Get any contained sets out of file_ids.
  EntityHandle first_set = fileInfo->sets.start_id;
  sets.merge( file_ids.lower_bound( first_set ),
              file_ids.lower_bound( first_set + fileInfo->sets.count ) );
  rval = read_sets( sets );
  if (MB_SUCCESS != rval)
    return error(rval);

  
  dbgOut.tprint( 1, "READING ADJACENCIES\n" );
    
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    if (fileInfo->elems[i].have_adj &&
        idMap.intersects( fileInfo->elems[i].desc.start_id, fileInfo->elems[i].desc.count )) {
      long len;
      hid_t th = mhdf_openAdjacency( filePtr, fileInfo->elems[i].handle, &len, &status );
      if (is_error(status))
        return error(MB_FAILURE);

      rval = read_adjacencies( th, len );
      mhdf_closeData( filePtr, th, &status );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
  }
  
  dbgOut.tprint( 1, "READING TAGS\n" );
  
  for (int i = 0; i < fileInfo->num_tag_desc; ++i) {
    rval = read_tag( i );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
  dbgOut.tprint( 1, "PARTIAL READ COMPLETE.\n" );
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::search_tag_values( int tag_index,
                                       const std::vector<int>& sorted_values,
                                       Range& file_ids,
                                       bool sets_only )
{
  ErrorCode rval;
  mhdf_Status status;
  std::vector<EntityHandle>::iterator iter;
  const mhdf_TagDesc& tag = fileInfo->tags[tag_index];
  long size;
  long start_id;

  debug_barrier();
   
    // do dense data
    
  hid_t table;
  const char* name;
  std::vector<EntityHandle> indices;
    // These are probably in order of dimension, so iterate
    // in reverse order to make Range insertions more efficient.
  std::vector<int> grp_indices( tag.dense_elem_indices, tag.dense_elem_indices+tag.num_dense_indices );
  for (std::vector<int>::reverse_iterator i = grp_indices.rbegin(); i != grp_indices.rend(); ++i)
  {
    int idx = *i;
    if (idx == -2) {
      name = mhdf_set_type_handle();
      start_id = fileInfo->sets.start_id;
    }
    else if (sets_only) {
      continue;
    }
    else if (idx == -1) {
      name = mhdf_node_type_handle();
     start_id = fileInfo->nodes.start_id;
    }
    else {
      if (idx < 0 || idx >= fileInfo->num_elem_desc) 
        return error(MB_FAILURE);
      name = fileInfo->elems[idx].handle;
      start_id = fileInfo->elems[idx].desc.start_id;
    }
    table = mhdf_openDenseTagData( filePtr, tag.name, name, &size, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    rval = search_tag_values( table, size, sorted_values, indices );
    mhdf_closeData( filePtr, table, &status );
    if (MB_SUCCESS != rval || is_error(status))
      return error(MB_FAILURE);
      // Convert from table indices to file IDs and add to result list
    std::sort( indices.begin(), indices.end(), std::greater<EntityHandle>() );
    std::transform( indices.begin(), indices.end(), range_inserter(file_ids),
                    std::bind1st( std::plus<long>(), start_id ) );
    indices.clear();
  }
  
  if (!tag.have_sparse)
    return MB_SUCCESS;
  
    // do sparse data
    
  hid_t tables[2]; 
  long junk; // redundant value for non-variable-length tags
  mhdf_openSparseTagData( filePtr, tag.name, &size, &junk, tables, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  rval = search_tag_values( tables[1], size, sorted_values, indices );
  mhdf_closeData( filePtr, tables[1], &status );
  if (MB_SUCCESS != rval || is_error(status)) {
    mhdf_closeData( filePtr, tables[0], &status );
    return error(MB_FAILURE);
  }
    // convert to ranges
  std::sort( indices.begin(), indices.end() );
  std::vector<EntityHandle> ranges;
  iter = indices.begin();
  while (iter != indices.end()) {
    ranges.push_back( *iter );
    EntityHandle last = *iter;
    for (++iter; iter != indices.end() && (last + 1) == *iter; ++iter, ++last);
    ranges.push_back( last );
  }
    // read file ids
  iter = ranges.begin();
  unsigned long offset = 0;
  while (iter != ranges.end()) {
    long begin = *iter; ++iter;
    long end   = *iter; ++iter;
    mhdf_readSparseTagEntitiesWithOpt( tables[0], begin, end - begin + 1, 
                                handleType, &indices[offset], indepIO, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, tables[0], &status );
      return error(MB_FAILURE);
    }
    offset += end - begin + 1;
  }
  mhdf_closeData( filePtr, tables[0], &status );
  if (is_error(status))
    return error(MB_FAILURE);
  assert( offset == indices.size() );
  std::sort( indices.begin(), indices.end() );
  
  if (sets_only) {
    iter = std::lower_bound( indices.begin(), indices.end(), 
                             fileInfo->sets.start_id + fileInfo->sets.count );
    indices.erase( iter, indices.end() );
    iter = std::lower_bound( indices.begin(), indices.end(), 
                             fileInfo->sets.start_id );
    indices.erase( indices.begin(), iter );
  }
  copy_sorted_file_ids( &indices[0], indices.size(), file_ids );
  
  return MB_SUCCESS;  
}

ErrorCode ReadHDF5::get_tagged_entities( int tag_index, Range& file_ids )
{
  const mhdf_TagDesc& tag = fileInfo->tags[tag_index];
   
    // do dense data
  Range::iterator hint = file_ids.begin();
  for (int i = 0; i < tag.num_dense_indices; ++i)
  {
    int idx = tag.dense_elem_indices[i];
    mhdf_EntDesc* ents;
    if (idx == -2)
      ents = &fileInfo->sets;
    else if (idx == -1) 
      ents = &fileInfo->nodes;
    else {
      if (idx < 0 || idx >= fileInfo->num_elem_desc) 
        return error(MB_FAILURE);
      ents = &(fileInfo->elems[idx].desc);
    }
    
    EntityHandle h = (EntityHandle)ents->start_id;
    hint = file_ids.insert( hint, h, h + ents->count - 1 );
  }
  
  if (!tag.have_sparse)
    return MB_SUCCESS;
  
    // do sparse data
    
  mhdf_Status status;
  hid_t tables[2]; 
  long size, junk; 
  mhdf_openSparseTagData( filePtr, tag.name, &size, &junk, tables, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  mhdf_closeData( filePtr, tables[1], &status );
  if (is_error(status)) {
    mhdf_closeData( filePtr, tables[0], &status );
    return error(MB_FAILURE);
  }
  
  hint = file_ids.begin();
  EntityHandle* buffer = reinterpret_cast<EntityHandle*>(dataBuffer);
  const long buffer_size = bufferSize / sizeof(EntityHandle);
  long remaining = size, offset = 0;
  while (remaining) {
    long count = std::min( buffer_size, remaining );
    assert_range( buffer, count );
    mhdf_readSparseTagEntitiesWithOpt( *tables, offset, count, 
                                handleType, buffer, collIO, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, *tables, &status );
      return error(MB_FAILURE);
    }
    
    std::sort( buffer, buffer + count );
    for (long i = 0; i < count; ++i)
      hint = file_ids.insert( hint, buffer[i], buffer[i] );
    
    remaining -= count;
    offset += count;
  }

  mhdf_closeData( filePtr, *tables, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  
  return MB_SUCCESS;  
}

ErrorCode ReadHDF5::search_tag_values( hid_t tag_table, 
                                       unsigned long table_size,
                                       const std::vector<int>& sorted_values,
                                       std::vector<EntityHandle>& value_indices )
{

  debug_barrier();

  mhdf_Status status;
  size_t chunk_size = bufferSize / sizeof(unsigned);
  unsigned * buffer = reinterpret_cast<unsigned*>(dataBuffer);
  size_t remaining = table_size, offset = 0;
  while (remaining) {
      // Get a block of tag values
    size_t count = std::min( chunk_size, remaining );
    assert_range( buffer, count );
    mhdf_readTagValuesWithOpt( tag_table, offset, count, H5T_NATIVE_UINT, buffer, collIO, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    
      // search tag values
    for (size_t i = 0; i < count; ++i)
      if (std::binary_search( sorted_values.begin(), sorted_values.end(), buffer[i] ))
        value_indices.push_back( i + offset );
    
    offset += count;
    remaining -= count;
  }

  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_nodes( const Range& node_file_ids )
{
  ErrorCode rval;
  mhdf_Status status;
  const int dim = fileInfo->nodes.vals_per_ent;
  Range range;
  
  if (node_file_ids.empty())
    return MB_SUCCESS;
  
  int cdim;
  rval = iFace->get_dimension( cdim );
  if (MB_SUCCESS != rval)
    return error(rval);
  
  if (cdim < dim)
  {
    rval = iFace->set_dimension( dim );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
  hid_t data_id = mhdf_openNodeCoordsSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);

  EntityHandle handle;
  std::vector<double*> arrays(dim);
  const size_t num_nodes = node_file_ids.size();
  rval = readUtil->get_node_coords( dim, (int)num_nodes, 0, handle, arrays );
  if (MB_SUCCESS != rval)
  {
    mhdf_closeData( filePtr, data_id, &status );
    return error(rval);
  }

#ifdef BLOCKED_COORD_IO
  try {
    for (int d = 0; d < dim; ++d) {
      ReadHDF5Dataset reader( data_id, H5T_NATIVE_DOUBLE, d );
      reader.set_file_ids( node_file_ids, fileInfo->nodes.start_id, collIO, mpiComm );
      dbgOut.tprintf(1,"Reading Node %c Coords (mode=%s)\n",(char)('X' + d),reader.get_mode_str()); 
      // should normally only have one read call, unless sparse nature
      // of file_ids caused reader to do something strange
      size_t remaining = num_nodes;
      while (!reader.done()) {
        size_t count;
        reader.read( arrays[d], remaining, count );
        arrays[d] += count;
        remaining -= count;
      }
      assert(!remaining);
    }
  }
  catch (ReadHDF5Dataset::Exception e) {
    mhdf_closeData( filePtr, data_id, &status );
    return error(MB_FAILURE);
  }
  mhdf_closeData( filePtr, data_id, &status );
#else
  double* buffer = (double*)dataBuffer;
  long chunk_size = bufferSize / (3*sizeof(double));
  long coffset = 0;
  try {
    ReadHDF5Dataset reader( data_id, H5T_NATIVE_DOUBLE );
    reader.set_file_ids( node_file_ids, fileInfo->nodes.start_id, indepIO, mpiComm );
    dbgOut.tprintf(1,"Reading Node Coords with mode = %s\n", reader.get_mode_str() );
    while (!reader.done()) {
      dbgOut.tprintf(2,"Reading chunk %ld of node coords\n", coffset/chunk_size+1);
      
      size_t count;
      reader.read( buffer, chunk_size, count );
      
      for (size_t i = 0; i < count; ++i)
        for (int d = 0; d < dim; ++d) 
          arrays[d][coffset+i] = buffer[dim*i+d];
      coffset += count;
    }
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
#endif
  for (int d = dim; d < cdim; ++d)
    memset( arrays[d], 0, num_nodes*sizeof(double) );
    
  return insert_in_id_map( node_file_ids, handle );
}

ErrorCode ReadHDF5::read_elems( int i )
{
  Range ids;
  ids.insert( fileInfo->elems[i].desc.start_id,
              fileInfo->elems[i].desc.start_id + fileInfo->elems[i].desc.count - 1);
  return read_elems( i, ids );
}

ErrorCode ReadHDF5::read_elems( int i, const Range& file_ids )
{
  if (fileInfo->elems[i].desc.vals_per_ent < 0)
    return read_poly( fileInfo->elems[i], file_ids );
  else
    return read_elems( fileInfo->elems[i], file_ids );
}

ErrorCode ReadHDF5::read_elems( const mhdf_ElemDesc& elems, const Range& file_ids )
{

  debug_barrier();

  ErrorCode rval = MB_SUCCESS;
  mhdf_Status status;
  
  EntityType type = CN::EntityTypeFromName( elems.type );
  if (type == MBMAXTYPE)
  {
    readUtil->report_error( "Unknown element type: \"%s\".\n", elems.type );
    return error(MB_FAILURE);
  }
  
  const int nodes_per_elem = elems.desc.vals_per_ent;
  const size_t count = file_ids.size();
  hid_t data_id = mhdf_openConnectivitySimple( filePtr, elems.handle, &status );
  if (is_error(status))
    return error(MB_FAILURE);

  EntityHandle handle;
  EntityHandle* array = 0;
  rval = readUtil->get_element_connect( count, nodes_per_elem, type,
                                        0, handle, array );
  if (MB_SUCCESS != rval)
    return error(rval);
  
  try {
    ReadHDF5Dataset reader( data_id, handleType );
    reader.set_file_ids( file_ids, elems.desc.start_id, collIO, mpiComm );
    dbgOut.tprintf(1,"Reading elements in '%s' using mode '%s'\n", 
                   elems.handle, reader.get_mode_str());
    // should normally only have one read call, unless sparse nature
    // of file_ids caused reader to do something strange
    EntityHandle* iter = array;
    size_t remaining = count;
    while (!reader.done()) {
      size_t num_read;
      reader.read( iter, remaining, num_read );
      iter += num_read * nodes_per_elem;
      remaining -= num_read;
    }
    assert(!remaining);
    assert(iter - array == (ptrdiff_t)count * nodes_per_elem);
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
  
  rval = convert_id_to_handle( array, count*nodes_per_elem );
  if (MB_SUCCESS != rval)
    return error(rval);
  
  rval = readUtil->update_adjacencies( handle, count, nodes_per_elem, array );
  if (MB_SUCCESS != rval)
    return error(rval);
 
  return insert_in_id_map( file_ids, handle );
}

ErrorCode ReadHDF5::read_node_adj_elems( const mhdf_ElemDesc& group )
{
  mhdf_Status status;
  ErrorCode rval;
  
  hid_t table = mhdf_openConnectivitySimple( filePtr, group.handle, &status );
  if (is_error(status))
    return error(MB_FAILURE);
    
  rval = read_node_adj_elems( group, table );
  
  mhdf_closeData( filePtr, table, &status );
  if (MB_SUCCESS == rval && is_error(status))
    return error(rval = MB_FAILURE);
  return rval;
}

ErrorCode ReadHDF5::read_node_adj_elems( const mhdf_ElemDesc& group, 
                                         hid_t table_handle )
{

  debug_barrier();

  mhdf_Status status;
  ErrorCode rval;
  IODebugTrack debug_track( debugTrack, std::string(group.handle) );

    // copy data to local variables (makes other code clearer)
  const int node_per_elem = group.desc.vals_per_ent;
  long start_id = group.desc.start_id;
  long remaining = group.desc.count;
  const EntityType type = CN::EntityTypeFromName( group.type );
  
    // figure out how many elements we can read in each pass
  long* const buffer = reinterpret_cast<long*>( dataBuffer );
  const long buffer_size = bufferSize / (node_per_elem * sizeof(buffer[0]));
    // read all element connectivity in buffer_size blocks
  long offset = 0;
  while (remaining) {
      // read a block of connectivity data
    const long count = std::min( remaining, buffer_size );
    debug_track.record_io( offset, count );
    assert_range( buffer, count*node_per_elem );
    mhdf_readConnectivityWithOpt( table_handle, offset, count, H5T_NATIVE_LONG, buffer, collIO, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    offset += count;
    remaining -= count;
    
      // count the number of elements in the block that we want,
      // zero connectivity for other elements
    long num_elem = 0;
    long* iter = buffer;
    for (long i = 0; i < count; ++i) {
      for (int j = 0; j < node_per_elem; ++j) {
        iter[j] = (long)idMap.find( iter[j] );
        if (!iter[j]) {
          iter[0] = 0;
          break;
        }
      }
      if (iter[0])
        ++num_elem;
      iter += node_per_elem;
    }
    
    if (!num_elem) {
      start_id += count;
      continue;
    }
    
      // create elements
    EntityHandle handle;
    EntityHandle* array;
    rval = readUtil->get_element_connect( (int)num_elem,
                                         node_per_elem,
                                         type,
                                         0,
                                         handle, 
                                         array );
    if (MB_SUCCESS != rval)
      return error(rval);
    
      // copy all non-zero connectivity values
    iter = buffer;
    EntityHandle* iter2 = array;
    EntityHandle h = handle;
    for (long i = 0; i < count; ++i) {
      if (!*iter) {
        iter += node_per_elem;
        continue;
      }
      if (!idMap.insert( start_id + i, h++, 1 ).second) 
        return error(MB_FAILURE);
        
      long* const end = iter + node_per_elem;
      for (; iter != end; ++iter, ++iter2)
        *iter2 = (EntityHandle)*iter;
    }
    assert( iter2 - array == num_elem * node_per_elem );
    start_id += count;
  }
  
  debug_track.all_reduce();
  return MB_SUCCESS;
}
  

ErrorCode ReadHDF5::read_elems( int i, const Range& elems_in, Range& nodes )
{
  EntityHandle* const buffer = reinterpret_cast<EntityHandle*>(dataBuffer);
  const int node_per_elem = fileInfo->elems[i].desc.vals_per_ent;
  const size_t buffer_size = bufferSize / (node_per_elem*sizeof(EntityHandle));
  
  if (elems_in.empty())
    return MB_SUCCESS;
    
  assert( (long)elems_in.front() >= fileInfo->elems[i].desc.start_id );
  assert( (long)elems_in.back() - fileInfo->elems[i].desc.start_id < fileInfo->elems[i].desc.count );
  
    // we don't support version 3 style poly element data
  if (fileInfo->elems[i].desc.vals_per_ent <= 0)
    return error(MB_TYPE_OUT_OF_RANGE);
  
  mhdf_Status status;
  hid_t table = mhdf_openConnectivitySimple( filePtr, fileInfo->elems[i].handle, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  
  try {
    ReadHDF5Dataset reader( table, handleType );
    reader.set_file_ids( elems_in, fileInfo->elems[i].desc.start_id, indepIO, mpiComm );
    dbgOut.tprintf(1,"Reading node ids for elements in '%s' using mode '%s'\n",
                   fileInfo->elems[i].handle, reader.get_mode_str() );
    while (!reader.done()) {
      size_t num_read;
      reader.read( buffer, buffer_size, num_read );
      std::sort( buffer, buffer + num_read*node_per_elem );
      num_read = std::unique( buffer, buffer + num_read*node_per_elem ) - buffer;
      copy_sorted_file_ids( buffer, num_read, nodes );
    }
  } 
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_poly( const mhdf_ElemDesc& elems, const Range& file_ids )
{
  class PolyReader : public ContentReader {
    private:
      const EntityType type;
      Interface *const mb;
      IDMap& idMap;
    public:
    PolyReader( EntityType elem_type, Interface* iface, IDMap& id_map )
               : type(elem_type), mb(iface), idMap(id_map) 
               {}
    ErrorCode store_data( EntityHandle, long file_id, EntityHandle* conn, long len, bool )
    {
      size_t valid;
      convert_id_to_handle( conn, len, valid, idMap );
      if (valid != (size_t)len)
        return error(MB_ENTITY_NOT_FOUND);
      EntityHandle handle;
      ErrorCode rval = mb->create_element( type, conn, len, handle );
      if (MB_SUCCESS != rval)
        return error(rval);
      if (!idMap.insert( file_id, handle, 1 ).second) 
        return error(MB_FAILURE);
      return MB_SUCCESS;
    }
  };

  debug_barrier();
  
  EntityType type = CN::EntityTypeFromName( elems.type );
  if (type == MBMAXTYPE)
  {
    readUtil->report_error( "Unknown element type: \"%s\".\n", elems.type );
    return error(MB_FAILURE);
  }
  
  hid_t handles[2];
  mhdf_Status status;
  long num_poly, num_conn, first_id;
  mhdf_openPolyConnectivity( filePtr, elems.handle, &num_poly, &num_conn, &first_id, 
                             handles, &status );
  if (is_error(status))
    return error(MB_FAILURE);

  ReadHDF5Dataset offset_reader( handles[0], handleType, true );
  ReadHDF5Dataset connect_reader( handles[1], handleType, true );
  

  PolyReader tool( type, iFace, idMap );
  Range empty;
  return read_contents( tool, offset_reader, connect_reader,
                        file_ids, first_id, 0, num_poly, num_conn, empty );
}
/*
ErrorCode ReadHDF5::read_poly( const char* elem_group )
{
  ErrorCode rval;
  mhdf_Status status;
  char name[64];
   
  mhdf_getElemTypeName( filePtr, elem_group, name, sizeof(name), &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return error(MB_FAILURE);
  }
  EntityType type = CN::EntityTypeFromName( name );

  long count, first_id, data_len;
  hid_t handles[2];
  mhdf_openPolyConnectivity( filePtr, elem_group, &count, &data_len,
                             &first_id, handles, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return error(MB_FAILURE);
  }

  ElemSet empty_set;
  empty_set.type = CN::EntityTypeFromName( name );
  empty_set.type2 = elem_group;
  
  EntityHandle h;
  bool first = true;
  long connend = -1;
  std::vector<EntityHandle> connectivity; 
  for (long i = 0; i < count; ++i) {
    long prevend = connend;
    mhdf_readPolyConnIndicesWithOpt( handles[0], i, 1, H5T_NATIVE_LONG, &connend, ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(MB_FAILURE);
    }
    
    connectivity.resize( connend - prevend );
    mhdf_readPolyConnIDsWithOpt( handles[1], prevend+1, connectivity.size(), handleType,
                          &connectivity[0], ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(MB_FAILURE);
    }
    
    rval= convert_id_to_handle( &connectivity[0], connectivity.size() );
    if (MB_SUCCESS != rval)     {
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(rval);
    }
    
    rval = iFace->create_element( type, &connectivity[0], connectivity.size(), h );
    if (MB_SUCCESS != rval) 
    {
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(rval);
    }
    
    if (first || elemList.back().range.back() + 1 >= h) {
      elemList.push_back( empty_set );
      elemList.back().first_id = first_id + i;
      first = false;
    }
    elemList.back().range.insert( h );
  }
 
  ErrorCode result = MB_SUCCESS;
  mhdf_closeData( filePtr, handles[0], &status );
  if (mhdf_isError( &status )) {
    readUtil->report_error( mhdf_message( &status ));
    result = error(MB_FAILURE);
  }
  mhdf_closeData( filePtr, handles[1], &status );
  if (mhdf_isError( &status )) {
    readUtil->report_error( mhdf_message( &status ));
    result = error(MB_FAILURE);
  }
  return result;
}
*/

ErrorCode ReadHDF5::read_elements_and_sides( const Range& file_ids )
{
  ErrorCode rval;
  EntityType type;
    
    // determine largest element dimension
  int max_dim = 0;
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    type = CN::EntityTypeFromName( fileInfo->elems[i].type );
    int dim = CN::Dimension(type);
    if (dim > max_dim) {
      EntityHandle start = (EntityHandle)fileInfo->elems[i].desc.start_id;
      Range::iterator it = file_ids.lower_bound( start );
      if (it != file_ids.end() && (long)(*it - start ) < fileInfo->elems[i].desc.count)
        max_dim = dim;
    }
  }
  
    // Need to do some communication here to make sure that
    // things work correctly if some proc ended up w/ no entities
#ifdef USE_MPI
  if (nativeParallel) {
    int send_val = max_dim;
    MPI_Allreduce( &send_val, &max_dim, 1, MPI_INT, MPI_MAX, 
                   myPcomm->proc_config().proc_comm() );
  }
#endif

  dbgOut.printf(3, "read_elements_and_sides: max_dim = %d\n", max_dim);
  
    // Get Range of element IDs only
  Range elem_ids( file_ids );
  Range::iterator s, e;
  if (fileInfo->nodes.count) {
    EntityHandle first = (EntityHandle)fileInfo->nodes.start_id;
    s = elem_ids.lower_bound( first );
    e = Range::lower_bound( s, elem_ids.end(), first + fileInfo->nodes.count );
    elem_ids.erase( s, e );
  }
  if (fileInfo->sets.count) {
    EntityHandle first = (EntityHandle)fileInfo->sets.start_id;
    s = elem_ids.lower_bound( first );
    e = Range::lower_bound( s, elem_ids.end(), first + fileInfo->sets.count );
    elem_ids.erase( s, e );
  }
  
    // read all node-adjacent elements of smaller dimensions
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    type = CN::EntityTypeFromName( fileInfo->elems[i].type );
    if (CN::Dimension(type) < max_dim) {
      dbgOut.tprintf(3, "   Reading node-adjacent elements for: %s\n", fileInfo->elems[i].handle);
      rval = read_node_adj_elems( fileInfo->elems[i] );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
  }
  
    // Read the subset of the explicitly sepecified elements that are of the
    // largest dimension.   We could read all explicitly specified elements,
    // but then the lower-dimension elements will be read a second time when 
    // reading node-adjacent elements.
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    type = CN::EntityTypeFromName( fileInfo->elems[i].type );
    if (CN::Dimension(type) == max_dim) {
      Range subset;
      intersect( fileInfo->elems[i].desc, elem_ids, subset );
      if (!subset.empty() || nativeParallel) {
        subtract( elem_ids,  subset );
        dbgOut.tprintf(3, "   Reading connectivity for: %s\n", fileInfo->elems[i].handle);
        rval = read_elems( fileInfo->elems[i],  subset );
        if (MB_SUCCESS != rval)
          return error(rval); 
      } 
    }
  } 

  debug_barrier();
  dbgOut.tprint(3, "  Deleting entities\n");
      
    // delete anything we read in that we should not have
    // (e.g. an edge spanning two disjoint blocks of elements)
    
    // get Range of all explicitly specified elements, grouped by dimension
  Range explicit_elems[4];
  Range::iterator hints[4] = { explicit_elems[0].begin(),
                               explicit_elems[1].begin(),
                               explicit_elems[2].begin(),
                               explicit_elems[3].begin() };
  RangeMap<long, EntityHandle>::iterator rit;
  while (!elem_ids.empty()) {
    long start = elem_ids.front();
    long count = elem_ids.const_pair_begin()->second - start + 1;
    rit = idMap.lower_bound( start );
    assert( rit != idMap.end() );
    assert( rit->begin <= start && rit->begin + rit->count > start );
    long offset = start - rit->begin;
    long avail = rit->count - offset;
    if (avail > count)
      count = avail;
    EntityHandle start_h = rit->value + offset;
    int d = CN::Dimension( TYPE_FROM_HANDLE( start_h ) );
    if (CN::Dimension( TYPE_FROM_HANDLE( start_h + count - 1 ) ) != d) 
      count = start - LAST_HANDLE( TYPE_FROM_HANDLE(start_h) ) + 1;
    
    elem_ids.erase( elem_ids.begin(), elem_ids.begin() + count );
    hints[d] = explicit_elems[d].insert( hints[d], rit->value + offset, rit->value + offset + count - 1 );
  }
  
    // get Range of all read elements, and remove any that were explicity specified
    // get handles for everything we've read
  Range all_elems;
  Range::iterator hint = all_elems.begin();
  for (rit = idMap.begin(); rit != idMap.end(); ++rit)
    hint = all_elems.insert( hint, rit->value, rit->value + rit->count - 1 );
    // remove any vertex handles
  all_elems.erase( all_elems.begin(), all_elems.upper_bound( MBVERTEX ) );
    // remove any set handles and any handles for elements >= max_dim
  all_elems.erase( all_elems.lower_bound( CN::TypeDimensionMap[max_dim].first ), all_elems.end() );
    // remove explicit elements < max_dim
  if (!explicit_elems[1].empty() && max_dim > 1)
    all_elems = subtract( all_elems,  explicit_elems[1] );
  if (!explicit_elems[2].empty() && max_dim > 2)
    all_elems = subtract( all_elems,  explicit_elems[2] );
 
    // remove any elements that are adjacent to some explicity specified element.
  if (max_dim > 1 && !explicit_elems[2].empty()) {
    Range adj;
    rval = iFace->get_adjacencies( explicit_elems[2], 1, false, adj, Interface::UNION );
    if (MB_SUCCESS != rval)
      return error(rval);
    if (!adj.empty())
      all_elems = subtract( all_elems,  adj );
  }
  if (max_dim == 3) {
    Range adj;
    rval = iFace->get_adjacencies( explicit_elems[3], 1, false, adj, Interface::UNION );
    if (MB_SUCCESS != rval)
      return error(rval);
    if (!adj.empty())
      all_elems = subtract( all_elems,  adj );
    adj.clear();
    rval = iFace->get_adjacencies( explicit_elems[3], 2, false, adj, Interface::UNION );
    if (MB_SUCCESS != rval)
      return error(rval);
    if (!adj.empty())
      all_elems = subtract( all_elems,  adj );
  }
  
    // now delete anything remaining in all_elems
  rval = iFace->delete_entities( all_elems );
  if (MB_SUCCESS != rval)
    return error(rval);
  
    // remove dead entities from ID map
  while (!all_elems.empty()) {
    EntityHandle start = all_elems.front();
    EntityID count = all_elems.const_pair_begin()->second - start + 1;
    for (rit = idMap.begin(); rit != idMap.end(); ++rit) 
      if (rit->value <= start && (long)(start - rit->value) < rit->count)
        break;
    if (rit == idMap.end())
      return error(MB_FAILURE);
  
    EntityID offset = start - rit->value;
    EntityID avail = rit->count - offset;
    if (avail < count)
      count = avail;
    
    all_elems.erase( all_elems.begin(), all_elems.begin() + count );
    idMap.erase( rit->begin + offset, count );
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_sets( const Range& file_ids )
{

  debug_barrier();

  mhdf_Status status;
  ErrorCode rval;
  if (fileInfo->sets.count == 0 || file_ids.empty()) {
    assert(file_ids.empty());
    return MB_SUCCESS;
  }
  
  hid_t meta_handle = mhdf_openSetMetaSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);


    // create sets 
  Range ranged_set_ids;
  EntityHandle start_handle;
  rval = read_sets( file_ids, meta_handle, ranged_set_ids, start_handle );
  if (MB_SUCCESS != rval) {
    mhdf_closeData( filePtr, meta_handle, &status );
    return error(rval);
  }
    
    // read contents
  if (fileInfo->have_set_contents) {
    long len = 0;
    hid_t handle = mhdf_openSetData( filePtr, &len, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
    
    rval = read_contents( file_ids, start_handle, meta_handle, handle, len,
                          ranged_set_ids );
    mhdf_closeData( filePtr, handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    if (MB_SUCCESS != rval) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(rval);
    }
  }
  
    // read set child lists
  if (fileInfo->have_set_children) {
    long len = 0;
    hid_t handle = mhdf_openSetChildren( filePtr, &len, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
    
    rval = read_children( file_ids, start_handle, meta_handle, handle, len );
    mhdf_closeData( filePtr, handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    if (MB_SUCCESS != rval) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(rval);
    }
  }
  
    // read set parent lists
  if (fileInfo->have_set_parents) {
    long len = 0;
    hid_t handle = mhdf_openSetParents( filePtr, &len, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
    
    rval = read_parents( file_ids, start_handle, meta_handle, handle, len );
    mhdf_closeData( filePtr, handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    if (MB_SUCCESS != rval) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(rval);
    }
  }
    
  mhdf_closeData( filePtr, meta_handle, &status );
  return is_error(status) ? error(MB_FAILURE) : MB_SUCCESS;
}

ErrorCode ReadHDF5::read_set_ids_recursive( Range& sets_in_out,
                                            bool contained_sets,
                                            bool child_sets )
{
  if (!fileInfo->have_set_children)
    child_sets = false;
  if (!fileInfo->have_set_contents)
    contained_sets = false;
  if (!child_sets && !contained_sets)
    return MB_SUCCESS;

    // open data tables
  if (fileInfo->sets.count == 0) {
    assert( sets_in_out.empty() );
    return MB_SUCCESS;
  }
  
  if (!contained_sets && !child_sets)
    return MB_SUCCESS;
  
  hid_t meta_handle, content_handle = 0, child_handle = 0;
  
  mhdf_Status status;
  meta_handle = mhdf_openSetMetaSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  
  if (contained_sets) {
    long content_len = 0;
    content_handle = mhdf_openSetData( filePtr, &content_len, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
  }
  
  if (child_sets) {
    long child_len = 0;
    child_handle = mhdf_openSetChildren( filePtr, &child_len, &status );
    if (is_error(status)) {
      if (contained_sets)
        mhdf_closeData( filePtr, content_handle, &status );
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
  }
  
  ErrorCode rval = MB_SUCCESS;
  Range children, new_children(sets_in_out);
  do {
    children.clear();
    if (child_sets) {
      rval = read_child_ids( new_children, meta_handle, child_handle, children );
      if (MB_SUCCESS != rval)
        break;
    }
    if (contained_sets) {
      rval = read_contained_set_ids( new_children, meta_handle, content_handle, children );
      if (MB_SUCCESS != rval)
        break;
    }
    new_children = subtract( children,  sets_in_out );
    sets_in_out.merge( new_children );
  } while (!new_children.empty());
  
  if (child_sets) {
    mhdf_closeData( filePtr, child_handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = error(MB_FAILURE);
  }
  if (contained_sets) {
    mhdf_closeData( filePtr, content_handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = error(MB_FAILURE);
  }
  mhdf_closeData( filePtr, meta_handle, &status );
  if (MB_SUCCESS == rval && is_error(status))
    rval = error(MB_FAILURE);
  
  return rval;
}

ErrorCode ReadHDF5::find_sets_containing( Range& sets_out )
{
  ErrorCode rval;
  mhdf_Status status;

  if (!fileInfo->have_set_contents)
    return MB_SUCCESS;
  assert( fileInfo->sets.count );

    // open data tables
  hid_t meta_handle = mhdf_openSetMetaSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  long content_len = 0;
  hid_t content_handle = mhdf_openSetData( filePtr, &content_len, &status );
  if (is_error(status)) {
    mhdf_closeData( filePtr, meta_handle, &status );
    return error(MB_FAILURE);
  }

  rval = find_sets_containing( meta_handle, content_handle, content_len, sets_out );

  mhdf_closeData( filePtr, content_handle, &status );
  if(MB_SUCCESS == rval && is_error(status))
    rval = error(MB_FAILURE);
  mhdf_closeData( filePtr, meta_handle, &status );
  if(MB_SUCCESS == rval && is_error(status))
    rval = error(MB_FAILURE);
    
  return rval;
}

static bool set_map_intersect( unsigned short flags,
                               const long* contents,
                               int content_len,
                               const RangeMap<long,EntityHandle>& id_map  )
{
  if (flags & mhdf_SET_RANGE_BIT) {
    if (!content_len || id_map.empty())
      return false;
      
    const long* j = contents;
    const long* const end = contents + content_len;
    assert(content_len % 2 == 0);
    while (j != end) {
      long start = *(j++);
      long count = *(j++);
      if (id_map.intersects( start, count ))
        return true;
    }
  }
  else {
    const long* const end = contents + content_len;
    for (const long* i = contents; i != end; ++i)
      if (id_map.exists( *i ))
        return true;
  }
  return false;
}

ErrorCode ReadHDF5::find_sets_containing( hid_t meta_handle,
                                          hid_t contents_handle, 
                                          long contents_len,
                                          Range& file_ids )
{
  const long avg_set_len = contents_len / fileInfo->sets.count;
  long sets_per_buffer = bufferSize / (sizeof(short) + sizeof(long) * (2+avg_set_len));
    // round to down multiple of 8 to avoid alignment issues
  sets_per_buffer = 8 * (sets_per_buffer / 8);
  if (sets_per_buffer < 10) // just in case there's one huge set
    sets_per_buffer = 10;  
  unsigned short* flag_buffer = (unsigned short*)dataBuffer;
  long* offset_buffer = (long*)(flag_buffer + sets_per_buffer);
  long* content_buffer = offset_buffer + sets_per_buffer;
  assert(bufferSize % sizeof(long) == 0);
  long content_len = (long*)(dataBuffer + bufferSize) - content_buffer;
  assert(dataBuffer + bufferSize >= (char*)(content_buffer + content_len));
    // scan set table  
  mhdf_Status status;
  Range::iterator hint = file_ids.begin();
  long remaining = fileInfo->sets.count;
  long offset = 0;
  long prev_idx = -1;
  while (remaining) {
    long count = std::min( remaining, sets_per_buffer );
    assert_range( flag_buffer, count );
    mhdf_readSetFlagsWithOpt( meta_handle, offset, count, H5T_NATIVE_USHORT, flag_buffer, collIO, &status );
    if (is_error(status)) 
      return error(MB_FAILURE);
    assert_range( offset_buffer, count );
    mhdf_readSetContentEndIndicesWithOpt( meta_handle, offset, count, H5T_NATIVE_LONG, offset_buffer, collIO, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    
    long sets_remaining = count;
    long sets_offset = 0;
    while (sets_remaining) {
        // figure how many of the remaining sets are required to 
        // fill the set contents buffer.
      long sets_count = std::lower_bound( offset_buffer + sets_offset, 
                          offset_buffer + count, content_len + prev_idx )
                          - offset_buffer - sets_offset;
      if (!sets_count) { // contents of single set don't fit in buffer
        long content_remaining = offset_buffer[sets_offset] - prev_idx;
        long content_offset = prev_idx+1;
        while (content_remaining) {
          long content_count = content_len < content_remaining ?
                               2*(content_len/2) : content_remaining;
          assert_range( content_buffer, content_count );
          mhdf_readSetDataWithOpt( contents_handle, content_offset,
                                   content_count, H5T_NATIVE_LONG, 
                                   content_buffer, collIO, &status );
          if (is_error(status))
            return error(MB_FAILURE);
          if (set_map_intersect( flag_buffer[sets_offset],
                                 content_buffer, content_count, idMap )) {
            long id = fileInfo->sets.start_id + offset + sets_offset;
            hint = file_ids.insert( hint, id, id );
            break;
          }
          content_remaining -= content_count;
          content_offset += content_count;
        }
        prev_idx = offset_buffer[sets_offset];
        sets_count = 1;
      }
      else if (long read_num = offset_buffer[sets_offset + sets_count - 1] - prev_idx) {
        assert(sets_count > 0);
        assert_range( content_buffer, read_num );
        mhdf_readSetDataWithOpt( contents_handle, prev_idx+1, read_num, 
                                 H5T_NATIVE_LONG, content_buffer, collIO, &status );
        if (is_error(status))
          return error(MB_FAILURE);
        
        long* buff_iter = content_buffer;
        for (long i = 0; i < sets_count; ++i) {
          long set_size = offset_buffer[i+sets_offset] - prev_idx;
          prev_idx += set_size;
          if (set_map_intersect( flag_buffer[sets_offset+i],
                                 buff_iter, set_size, idMap )) {
            long id = fileInfo->sets.start_id + offset + sets_offset + i;
            hint = file_ids.insert( hint, id, id );
          }
          buff_iter += set_size;
        }
      }
    
      sets_offset += sets_count;
      sets_remaining -= sets_count;
    }
    
    offset += count;
    remaining -= count;
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_child_ids( const Range& input_file_ids,
                                    hid_t meta_handle,
                                    hid_t child_handle,
                                    Range& child_file_ids )
{
  mhdf_Status status;
  long* buffer = reinterpret_cast<long*>(dataBuffer);
  long buffer_size = bufferSize / sizeof(long);
  long first, range[2], count, remaining;
  Range sets(input_file_ids);
  Range::iterator hint;
  while (!sets.empty()) {
    count = (long)sets.const_pair_begin()->second - sets.front() + 1;
    first = (long)sets.front() - fileInfo->sets.start_id;
    sets.erase( sets.begin(), sets.begin() + count );
    
    if (!first) {
      range[0] = -1;
      mhdf_readSetChildEndIndicesWithOpt( meta_handle, first+count-1, 1, 
                                          H5T_NATIVE_LONG, range+1, 
                                          indepIO, &status );
      if (is_error(status))
        return error(MB_FAILURE);
    }
    else if (count == 1) {
      mhdf_readSetChildEndIndicesWithOpt( meta_handle, first-1, 2, 
                                          H5T_NATIVE_LONG, range, 
                                          indepIO, &status );
      if (is_error(status))
        return error(MB_FAILURE);
    }
    else {
      mhdf_readSetChildEndIndicesWithOpt( meta_handle, first-1, 1, 
                                          H5T_NATIVE_LONG, range, 
                                          indepIO, &status );
      if (is_error(status))
        return error(MB_FAILURE);
      mhdf_readSetChildEndIndicesWithOpt( meta_handle, first+count-1, 1, 
                                          H5T_NATIVE_LONG, range+1, 
                                          indepIO, &status );
      if (is_error(status))
        return error(MB_FAILURE);
    }
    
    if (range[0] > range[1]) 
      return error(MB_FAILURE);
    remaining = range[1] - range[0];
    long offset = range[0] + 1;
    while (remaining) {
      count = std::min( buffer_size, remaining );
      remaining -= count;
      assert_range( buffer, count );
      mhdf_readSetParentsChildrenWithOpt( child_handle, offset, count, 
                                          H5T_NATIVE_LONG, buffer, 
                                          indepIO, &status );
  
      std::sort( buffer, buffer + count );
      count = std::unique( buffer, buffer + count ) - buffer;
      hint = child_file_ids.begin();
      for (long i = 0; i < count; ++i) {
        EntityHandle h = (EntityHandle)buffer[i];
        hint = child_file_ids.insert( hint, h, h );
      }
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_contained_set_ids( const Range& input_file_ids,
                                            hid_t meta_handle,
                                            hid_t content_handle,
                                            Range& contained_set_file_ids )
{
  mhdf_Status status;
  long buffer_size = bufferSize / (sizeof(long) + sizeof(short));
    // don't want to worry about reading half of a range pair later
  if (buffer_size % 2) --buffer_size;
  long* content_buffer = reinterpret_cast<long*>(dataBuffer);
  unsigned short* flag_buffer = reinterpret_cast<unsigned short*>(content_buffer + buffer_size);
  long first, range[2], count, remaining, sets_offset;

  Range sets(input_file_ids);
  Range::iterator hint;
  while (!sets.empty()) {
    count = (long)sets.const_pair_begin()->second - sets.front() + 1;
    first = (long)sets.front() - fileInfo->sets.start_id;
    if (count > buffer_size)
      count = buffer_size;
    sets.erase( sets.begin(), sets.begin() + count );
    
    assert_range( flag_buffer, count );
    mhdf_readSetFlags( meta_handle, first, count, H5T_NATIVE_USHORT, flag_buffer, &status );
    if (is_error(status))
      return MB_FAILURE;
    
    sets_offset = 0;
    while (sets_offset < count) {
        // Find block of sets with same value for ranged flag
      long start_idx = sets_offset;
      unsigned short ranged = flag_buffer[start_idx] & mhdf_SET_RANGE_BIT;
      for (++sets_offset; sets_offset < count; ++sets_offset)
        if ((flag_buffer[sets_offset] & mhdf_SET_RANGE_BIT) != ranged)
          break;
          
      if (!first && !start_idx) { // first set
        range[0] = -1;
        mhdf_readSetContentEndIndicesWithOpt( meta_handle, first+sets_offset-1, 
                                              1, H5T_NATIVE_LONG, range+1, 
                                              indepIO, &status );
        if (is_error(status))
          return error(MB_FAILURE);
      }
      else if (count == 1) {
        mhdf_readSetContentEndIndicesWithOpt( meta_handle, first+start_idx-1, 
                                              2, H5T_NATIVE_LONG, range, 
                                              indepIO, &status );
        if (is_error(status))
          return error(MB_FAILURE);
      }
      else {
        mhdf_readSetContentEndIndicesWithOpt( meta_handle, first+start_idx-1, 
                                              1, H5T_NATIVE_LONG, range, 
                                              indepIO, &status );
        if (is_error(status))
          return error(MB_FAILURE);
        mhdf_readSetContentEndIndicesWithOpt( meta_handle, first+sets_offset-1, 
                                              1, H5T_NATIVE_LONG, range+1, 
                                              indepIO, &status );
        if (is_error(status))
          return error(MB_FAILURE);
      }
    
      remaining = range[1] - range[0];
      long offset = range[0] + 1;
      while (remaining) {
        assert( !ranged || !(remaining % 2) );
        long content_count = std::min( buffer_size, remaining );
        remaining -= content_count;
        assert_range( content_buffer, content_count );
        mhdf_readSetDataWithOpt( content_handle, offset, content_count, 
                                 H5T_NATIVE_LONG, content_buffer, indepIO, 
                                 &status );
  
        if (ranged) {
          hint = contained_set_file_ids.begin();
          for (long i = 0; i < content_count; i += 2) {
            EntityHandle s = (EntityHandle)content_buffer[i];
            EntityHandle e = s + content_buffer[i+1];
            if ((long)s < fileInfo->sets.start_id)
              s = fileInfo->sets.start_id;
            if ((long)e > fileInfo->sets.start_id + fileInfo->sets.count)
              e = fileInfo->sets.start_id + fileInfo->sets.count;
            if (s < e) 
              hint = contained_set_file_ids.insert( hint, s, e - 1 );
          }
        }
        else {
          std::sort( content_buffer, content_buffer + content_count );
          long* s = std::lower_bound( content_buffer, content_buffer + content_count,
                                      fileInfo->sets.start_id );
          long* e = std::lower_bound( s, content_buffer + content_count, 
                                      fileInfo->sets.start_id + fileInfo->sets.count );
          e = std::unique( s, e );
          hint = contained_set_file_ids.begin();
          for ( ; s != e; ++s) {
            EntityHandle h = *s;
            hint = contained_set_file_ids.insert( hint, h, h );
          }
        }
      }
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_sets( const Range& file_ids,
                               hid_t meta_handle, 
                               Range& ranged_file_ids,
                               EntityHandle& start_handle,
                               bool create )
{
  ErrorCode rval;

  size_t num_sets = file_ids.size();
  if (!num_sets)
    return MB_SUCCESS;

  std::vector<unsigned> tmp_buffer;
  unsigned* buffer;
  if (num_sets > bufferSize/sizeof(unsigned)) {
    tmp_buffer.resize( num_sets );
    buffer = &tmp_buffer[0];
  }
  else {
    buffer = reinterpret_cast<unsigned*>(dataBuffer);
  }
  
  try {
    ReadHDF5Dataset reader( meta_handle, H5T_NATIVE_UINT, false );
    reader.set_column( reader.columns() - 1 );
    reader.set_file_ids( file_ids, fileInfo->sets.start_id, collIO, mpiComm );
    dbgOut.tprintf(1,"Reading set descriptions with create = %s and mode = %s\n",
      create ? "true" : "false", reader.get_mode_str() );
    // should normally only have one read call, unless sparse nature
    // of file_ids caused reader to do something strange
    size_t remaining = num_sets;
    unsigned* iter = buffer;
    while (!reader.done()) {
      size_t count = 0;
      reader.read( iter, remaining, count );
      iter += count;
      remaining -= count;
    }
    assert(!remaining);
    assert(iter - buffer == (ptrdiff_t)num_sets);
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
  
  unsigned* buff_iter = buffer;
  Range::iterator hint = ranged_file_ids.begin();
  for (Range::iterator i = file_ids.begin(); i != file_ids.end(); ++i, ++buff_iter) {
    if ((*buff_iter) & mhdf_SET_RANGE_BIT) {
      *buff_iter &= ~(unsigned)mhdf_SET_RANGE_BIT;
      hint = ranged_file_ids.insert( hint, *i, *i );
    }
  }
  
  if (create) {
    rval = readUtil->create_entity_sets( num_sets, buffer, 0, start_handle );
    if (MB_SUCCESS != rval)
      return error(rval);
    
    rval = insert_in_id_map( file_ids, start_handle );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
  return MB_SUCCESS;
}


ErrorCode ReadHDF5::read_contents( ContentReader& tool,
                                   ReadHDF5Dataset& idx_reader,
                                   ReadHDF5Dataset& dat_reader,
                                   const Range& file_ids,
                                   const long start_id,
                                   const EntityHandle start_handle,
                                   const long entity_count,
                                   const long content_len,
                                   const Range& ranged_ids_in )
{
  ErrorCode rval;
  if (file_ids.empty())
    return MB_SUCCESS;
    
    // things will get messed up if this isn't true
  assert( subtract( ranged_ids_in, file_ids ).empty() );
  assert( file_ids.front() >= (EntityHandle)start_id );
  assert( file_ids.back() - start_id < (EntityHandle)entity_count );
  
  // Set up buffers
  const long avg_set_len = content_len / entity_count;
  long sets_per_buffer = bufferSize / (sizeof(long) + (avg_set_len+1)*sizeof(EntityHandle));
    // round to down multiple of 8 to avoid alignment issues
  sets_per_buffer = 8 * (sets_per_buffer / 8);
  if (sets_per_buffer < 10) { // just in case there's one huge set
    sets_per_buffer = 10;  
    if (sets_per_buffer * (long)sizeof(long) > bufferSize) 
      sets_per_buffer = bufferSize/(sizeof(long)+2*sizeof(EntityHandle));
  }
  long* offset_buffer = (long*)dataBuffer;
  EntityHandle* content_buffer = (EntityHandle*)(offset_buffer + sets_per_buffer);
  assert(bufferSize % sizeof(long) == 0);
  long content_size = (EntityHandle*)(dataBuffer + bufferSize) - content_buffer;
  assert(dataBuffer + bufferSize >= (char*)(content_buffer + content_size));
    // content_size must be an even number or we might end up with half of a 
    // range pair in the buffer, which will break code below
  if (content_size % 2)
    --content_size;
    
  if (!sets_per_buffer || !content_size) {
      // buffer is too small to be usable
    return MB_FAILURE;
  }

  // build list of offsets into end index table from file id list
  // Note: must make offsets one-based because Range can't hold zeros
  const EntityHandle one = 1;
  const EntityHandle off = start_id - one;
  Range::const_pair_iterator pi = file_ids.begin();
  Range offsets;
  Range::iterator ins = offsets.begin();
  Range::const_iterator range_iter = ranged_ids_in.begin();
  size_t prev_count = 0;
  if ((long)pi->first == start_id) {
    prev_count = 1;
    offset_buffer[0] = -1;
    ins = offsets.insert( ins, pi->first - off, pi->second - off );
    ++pi;
  }
  for ( ; pi != file_ids.const_pair_end(); ++pi) 
    ins = offsets.insert( ins, pi->first - off - 1, pi->second - off );

    
  assert( offsets.size() + prev_count == file_ids.size() + file_ids.psize() );
  try {
    idx_reader.set_file_ids( offsets, one, indepIO, mpiComm );
    idx_reader.set_data_type( H5T_NATIVE_LONG );
    dat_reader.set_data_type( handleType );
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
      
  // Do actual read
  EntityHandle h = start_handle;
  Range::iterator fid_iter = file_ids.begin();
  std::vector<EntityHandle> tmp_buffer;
  dbgOut.tprintf(1,"Reading set data/children/parents/polyconn with idx_reader.ioMode = %s\n", idx_reader.get_mode_str());
  while (!idx_reader.done()) {
    size_t count;
    try {
      idx_reader.read( offset_buffer+prev_count, sets_per_buffer-prev_count, count ); 
      count += prev_count;
    }
    catch (ReadHDF5Dataset::Exception e) {
      return error(MB_FAILURE);
    }
    
    Range rows;
    ins = rows.begin();
    long prev_end = offset_buffer[0];
    size_t r = 1, w = 0;
    Range::iterator tmp_iter = fid_iter;
    while (r < count) {
        // iterate until taking any more sets would overflow the buffer
      long s = offset_buffer[r] - prev_end;

        // build list of values from data table to read
      ins = rows.insert( ins, prev_end + 1 + one, offset_buffer[r] + one );

        // If this entity is immediately after the previous, then the
        // end if the previous data is the start of this data.  Otherwise
        // we will have read an extra value so we know the start of the
        // next chunk of data.  If the later is the case, just increment
        // r so that we take that extra value.
      EntityHandle curr = *tmp_iter;
      ++tmp_iter;
      if (tmp_iter != file_ids.end() && *tmp_iter != (curr+1))
        ++r;

        // overwrite end offsets with counts
      prev_end = offset_buffer[r];
      offset_buffer[w++] = s;
      ++r;
    }

      // we now have:
      // w: number of sets for which contents/whatever are to be read
      // offset_buffer[0..w) : number of values for each set
      // fid_iter: file id of first set to read
      // tmp_tier: one past last file id to read
      // rows    : offsets (plus one) of set values to read

      // read set contents
    size_t count2;
    dat_reader.set_file_ids( rows, one, indepIO, mpiComm );
    size_t i = 0;
    size_t prev = 0;
    dbgOut.tprintf(1,"Reading set data/children/parents/polyconn with dat_reader.ioMode = %s\n", dat_reader.get_mode_str());
    while (!dat_reader.done()) {
      try {
        dat_reader.read( content_buffer+prev, content_size-prev, count2 );
        count2 += prev;
      }
      catch (ReadHDF5Dataset::Exception e) {
        return error(MB_FAILURE);
      }

        // now process data for each set
      EntityHandle* cont_iter = content_buffer;
      for (; i < w; ++i) {
          // stop at first set for which not all content data has been read
        if (cont_iter + offset_buffer[i] > content_buffer + content_size)
          break;
      
        bool ranged = false;
        if (range_iter != ranged_ids_in.end() && *range_iter == *fid_iter) {
          ranged = true;
          ++range_iter;
        }

        rval = tool.store_data( h++,
                                *fid_iter,
                                cont_iter,
                                offset_buffer[i],
                                ranged );
        if (MB_SUCCESS != rval)
          return error(rval);

        cont_iter += offset_buffer[i];
        ++fid_iter;
      }
      
        // number of content entries not processed in this iteration
      prev = content_buffer + count2 - cont_iter;
      
        // if unread set has contents larger then entire buffer,
        // then we need a bigger buffer
      if (i < w && offset_buffer[i] > content_size) {
        content_size = offset_buffer[i];
        tmp_buffer.resize( content_size );
        content_buffer = &tmp_buffer[0];
      }
      
        // shift unprocessed entries to start of buffer
      memmove( content_buffer, cont_iter, sizeof(EntityHandle) * prev );
    }

      // set up next iteration:
      // do we need to keep end index from previous set?
    if (fid_iter != file_ids.end() && fid_iter.start_of_block() != fid_iter)
      --r;
      // move unconsumed portion of offset buffer to beginning
    prev_count = count - r;
    memmove( offset_buffer, offset_buffer+r, sizeof(EntityHandle)*prev_count );
  }
  
  return MB_SUCCESS;
}


ErrorCode ReadHDF5::read_contents( const Range& set_file_ids,
                                   EntityHandle start_handle,
                                   hid_t set_meta_data_table,
                                   hid_t set_contents_table,
                                   long set_contents_length,
                                   const Range& ranged_set_file_ids )
{

  class ReadSetContents : public ReadHDF5::ContentReader {
    Interface *const mb;
    const IDMap& idMap;
  public:
    ReadSetContents( Interface* iface, const IDMap& id_map )
                    : mb(iface), idMap(id_map) 
                    {}
    ErrorCode store_data( EntityHandle set, long, EntityHandle* array, long len, bool ranged ) 
    {
      if (ranged) {
        if (len % 2) 
          return error(MB_INDEX_OUT_OF_RANGE);
        Range range;
        convert_range_to_handle( array, len/2, idMap, range );
        return mb->add_entities( set, range );
      }
      else {
        if (!len)
          return MB_SUCCESS;
        size_t valid;
        convert_id_to_handle( array, len, valid, idMap );
        return mb->add_entities( set, array, valid );
      }
    }
  };
  
  dbgOut.print(1,"Reading set contents\n");
  try {
    ReadHDF5Dataset met_dat( set_meta_data_table, H5T_NATIVE_LONG, false );
    met_dat.set_column( 0 );
    ReadHDF5Dataset cnt_dat( set_contents_table, handleType, false );
    ReadSetContents tool( iFace, idMap );
    return read_contents( tool, met_dat, cnt_dat, 
                          set_file_ids, fileInfo->sets.start_id, start_handle, 
                          fileInfo->sets.count, set_contents_length, ranged_set_file_ids );
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
}


ErrorCode ReadHDF5::read_children( const Range& set_file_ids,
                                   EntityHandle start_handle,
                                   hid_t set_meta_data_table,
                                   hid_t set_contents_table,
                                   long set_contents_length )
{
  class ReadSetChildren : public ReadHDF5::ContentReader {
    Interface *const mb;
    const IDMap& idMap;
  public:
    ReadSetChildren( Interface* iface, const IDMap& id_map )
                    : mb(iface), idMap(id_map) 
                    {}
    ErrorCode store_data( EntityHandle set, long, EntityHandle* array, long len, bool ranged ) 
    {
      assert(!ranged);
      size_t valid;
      convert_id_to_handle( array, len, valid, idMap );
      return mb->add_child_meshsets( set, array, valid );
    }
  };

  dbgOut.print(1,"Reading set children\n");
  try {
    Range empty;
    ReadHDF5Dataset met_dat( set_meta_data_table, H5T_NATIVE_LONG, false );
    met_dat.set_column( 1 );
    ReadHDF5Dataset cnt_dat( set_contents_table, handleType, false );
    ReadSetChildren tool( iFace, idMap );
    return read_contents( tool, met_dat, cnt_dat,
                          set_file_ids, fileInfo->sets.start_id, start_handle, 
                          fileInfo->sets.count, set_contents_length, empty );
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
}

ErrorCode ReadHDF5::read_parents( const Range& set_file_ids,
                                  EntityHandle start_handle,
                                  hid_t set_meta_data_table,
                                  hid_t set_contents_table,
                                  long set_contents_length )
{
  class ReadSetParents : public ReadHDF5::ContentReader {
    Interface *const mb;
    const IDMap& idMap;
  public:
    ReadSetParents( Interface* iface, const IDMap& id_map )
                    : mb(iface), idMap(id_map) 
                    {}
    ErrorCode store_data( EntityHandle set, long, EntityHandle* array, long len, bool ranged ) 
    {
      assert(!ranged);
      size_t valid;
      convert_id_to_handle( array, len, valid, idMap );
      return mb->add_parent_meshsets( set, array, valid );
    }
  };

  dbgOut.print(1,"Reading set parents\n");
  try {
    Range empty;
    ReadHDF5Dataset met_dat( set_meta_data_table, H5T_NATIVE_LONG, false );
    met_dat.set_column( 2 );
    ReadHDF5Dataset cnt_dat( set_contents_table, handleType, false );
    ReadSetParents tool( iFace, idMap );
    return read_contents( tool, met_dat, cnt_dat,
                          set_file_ids, fileInfo->sets.start_id, start_handle, 
                          fileInfo->sets.count, set_contents_length, empty );
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
}

static void copy_set_contents( int ranged,
                               const EntityHandle* contents,
                               long length,
                               Range& results )
{
  if (ranged) {
    assert( length%2 == 0 );
    Range::iterator hint = results.begin();
    for (long i = 0; i < length; i += 2)
      hint = results.insert( hint, contents[i], contents[i] + contents[i+1] - 1 );
  }
  else {
    for (long i = 0; i < length; ++i)
      results.insert( contents[i] );
  }
}


ErrorCode ReadHDF5::get_set_contents( const Range& sets, Range& file_ids )
{
  class GetContentList : public ReadHDF5::ContentReader {
    Range *const resultList;
  public:
    GetContentList( Range* result_set )
                    : resultList(result_set) 
                    {}

    ErrorCode store_data( EntityHandle, long, EntityHandle* array, long len, bool ranged ) 
    {
      if (ranged) {
        if (len % 2)
          return error(MB_INDEX_OUT_OF_RANGE);
        copy_set_contents( 1, array, len, *resultList );
      }
      else {
        std::sort( array, array+len );
        copy_sorted_file_ids( array, len, *resultList );
      }
      return MB_SUCCESS;
    }
  };

  ErrorCode rval;
  if (sets.empty())
    return MB_SUCCESS;

  dbgOut.print(1,"Reading set contained file IDs\n");
  try {
    mhdf_Status status;
    long content_len;
    hid_t meta, contents;
    meta = mhdf_openSetMetaSimple( filePtr, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    ReadHDF5Dataset meta_reader( meta, H5T_NATIVE_LONG, true );
    meta_reader.set_column( 0 );

    contents = mhdf_openSetData( filePtr, &content_len, &status );
    if (is_error(status)) 
      return error(MB_FAILURE);
    ReadHDF5Dataset data_reader( contents, handleType, true );

    EntityHandle junk;
    Range ranged;
    rval = read_sets( sets, meta, ranged, junk, false );
    if (MB_SUCCESS != rval)
      return error(rval);

    GetContentList tool( &file_ids );
    rval = read_contents( tool, meta_reader, data_reader,
                          sets, fileInfo->sets.start_id, junk, 
                          fileInfo->sets.count, content_len, ranged );
    return rval;
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
}


ErrorCode ReadHDF5::read_adjacencies( hid_t table, long table_len )
{
  ErrorCode rval;
  mhdf_Status status;

  debug_barrier();

  
  EntityHandle* buffer = (EntityHandle*)dataBuffer;
  size_t chunk_size = bufferSize / sizeof(EntityHandle);
  size_t remaining = table_len;
  size_t left_over = 0;
  size_t offset = 0;
  while (remaining)
  {
    size_t count = std::min( chunk_size, remaining );
    count -= left_over;
    remaining -= count;
    
    assert_range( buffer + left_over, count );
    mhdf_readAdjacencyWithOpt( table, offset, count, handleType, buffer + left_over,
                               collIO, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    
    EntityHandle* iter = buffer;
    EntityHandle* end = buffer + count + left_over;
    while (end - iter >= 3)
    {
      EntityHandle h = idMap.find( *iter++ );
      EntityHandle count2 = *iter++;
      if (!h) {
        iter += count2;
        continue;
      }

      if (count2 < 1)
        return error(MB_FAILURE);

      if (end < count2 + iter)
      {
        iter -= 2;
        break;
      }
      
      size_t valid;
      convert_id_to_handle( iter, count2, valid, idMap );
      rval = iFace->add_adjacencies( h, iter, valid, false );
      if (MB_SUCCESS != rval)
        return error(rval);
     
      iter += count2;
    }
    
    left_over = end - iter;
    assert_range( (char*)buffer, left_over );
    assert_range( (char*)iter, left_over );
    memmove( buffer, iter, left_over );
  }
  
  assert(!left_over);  // unexpected truncation of data
  
  return MB_SUCCESS;  
}


ErrorCode ReadHDF5::read_tag( int tag_index )
{
  dbgOut.tprintf(2, "Reading tag \"%s\"\n", fileInfo->tags[tag_index].name );

  debug_barrier();


  ErrorCode rval;
  mhdf_Status status;
  Tag tag = 0;
  hid_t read_type;
  bool table_type;
  rval = create_tag( fileInfo->tags[tag_index], tag, read_type ); 
  if (MB_SUCCESS != rval)
    return error(rval);

  if (fileInfo->tags[tag_index].have_sparse) {
    hid_t handles[3];
    long num_ent, num_val;
    mhdf_openSparseTagData( filePtr, 
                            fileInfo->tags[tag_index].name,
                            &num_ent, &num_val,
                            handles, &status );
    if (is_error(status)) {
      if (read_type) H5Tclose( read_type );
      return error(MB_FAILURE);
    }
    
    table_type = false;
    if (read_type == 0) {
      read_type = H5Dget_type( handles[1] );
      if (read_type == 0) {
        mhdf_closeData( filePtr, handles[0], &status );
        mhdf_closeData( filePtr, handles[0], &status );
        if (fileInfo->tags[tag_index].size <= 0) 
          mhdf_closeData( filePtr, handles[2], &status );
        return error(MB_FAILURE);
      }
      table_type = true;
    }

    if (fileInfo->tags[tag_index].size > 0) {
      dbgOut.printf(1, "Reading sparse data for tag \"%s\"\n", fileInfo->tags[tag_index].name );
      rval = read_sparse_tag( tag, read_type, handles[0], handles[1], num_ent );
    }
    else {
      dbgOut.printf(1, "Reading var-len sparse data for tag \"%s\"\n", fileInfo->tags[tag_index].name );
      rval = read_var_len_tag( tag, read_type, handles[0], handles[1], handles[2], num_ent, num_val );
    }

    if (table_type) {
      H5Tclose(read_type);
      read_type = 0;
    }
    
    mhdf_closeData( filePtr, handles[0], &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    mhdf_closeData( filePtr, handles[1], &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    if (fileInfo->tags[tag_index].size <= 0) {
      mhdf_closeData( filePtr, handles[2], &status );
      if (MB_SUCCESS == rval && is_error(status))
        rval = MB_FAILURE;
    }
    if (MB_SUCCESS != rval) {
      if (read_type) H5Tclose( read_type );
      return error(rval);
    }
  }
  
  for (int j = 0; j < fileInfo->tags[tag_index].num_dense_indices; ++j) {
    long count;
    const char* name = 0;
    mhdf_EntDesc* desc;
    int elem_idx = fileInfo->tags[tag_index].dense_elem_indices[j];
    if (elem_idx == -2) {
      desc = &fileInfo->sets;
      name = mhdf_set_type_handle();
    }
    else if (elem_idx == -1) {
      desc = &fileInfo->nodes;
      name = mhdf_node_type_handle();
    }
    else if (elem_idx >= 0 && elem_idx < fileInfo->num_elem_desc) {
      desc = &fileInfo->elems[elem_idx].desc;
      name = fileInfo->elems[elem_idx].handle;
    }
    else {
      return error(MB_FAILURE);
    }
    
    dbgOut.printf(1, "Read dense data block for tag \"%s\" on \"%s\"\n", fileInfo->tags[tag_index].name, name );
    
    hid_t handle = mhdf_openDenseTagData( filePtr, 
                                          fileInfo->tags[tag_index].name,
                                          name,
                                          &count, &status );
    if (is_error(status)) {
      rval = error(MB_FAILURE);
      break;
    }
    
    if (count > desc->count) {
      readUtil->report_error( "Invalid data length for dense tag data: %s/%s\n",
                              name, fileInfo->tags[tag_index].name );
      mhdf_closeData( filePtr, handle, &status );
      rval = error(MB_FAILURE);
      break;
    }
    
    table_type = false;
    if (read_type == 0) {
      read_type = H5Dget_type( handle );
      if (read_type == 0) {
        mhdf_closeData( filePtr, handle, &status );
        return error(MB_FAILURE);
      }
      table_type = true;
    }

    rval = read_dense_tag( tag, read_type, handle, desc->start_id, count );
    
    if (table_type) {
      H5Tclose( read_type );
      read_type = 0;
    }
    
    mhdf_closeData( filePtr, handle, &status );
    if (MB_SUCCESS != rval)
      break;
    if (is_error(status)) {
      rval = error(MB_FAILURE);
      break;
    }
  }
  
  if (read_type) 
    H5Tclose( read_type );
  return rval;
}
                              
                                            
    


ErrorCode ReadHDF5::create_tag( const mhdf_TagDesc& info,
                                Tag& handle,
                                hid_t& hdf_type )
{
  ErrorCode rval;
  mhdf_Status status;
  TagType storage;
  DataType mb_type;
  bool re_read_default = false;

  switch (info.storage) {
    case mhdf_DENSE_TYPE : storage = MB_TAG_DENSE ; break;
    case mhdf_SPARSE_TYPE: storage = MB_TAG_SPARSE; break;
    case mhdf_BIT_TYPE   : storage = MB_TAG_BIT;    break;
    case mhdf_MESH_TYPE  : storage = MB_TAG_MESH;   break;
    default:
      readUtil->report_error( "Invalid storage type for tag '%s': %d\n", info.name, info.storage );
      return error(MB_FAILURE);
  }

    // Type-specific stuff
  if (info.type == mhdf_BITFIELD) {
    if (info.size < 1 || info.size > 8)
    {
      readUtil->report_error( "Invalid bit tag:  class is MB_TAG_BIT, num bits = %d\n", info.size );
      return error(MB_FAILURE);
    }
    hdf_type = H5Tcopy(H5T_NATIVE_B8);
    mb_type = MB_TYPE_BIT;
    if (hdf_type < 0)
      return error(MB_FAILURE);
  }
  else if (info.type == mhdf_OPAQUE) {
    mb_type = MB_TYPE_OPAQUE;

      // Check for user-provided type
    Tag type_handle;
    std::string tag_type_name = "__hdf5_tag_type_";
    tag_type_name += info.name;
    rval = iFace->tag_get_handle( tag_type_name.c_str(), type_handle );
    if (MB_SUCCESS == rval) {
      rval = iFace->tag_get_data( type_handle, 0, 0, &hdf_type );
      if (MB_SUCCESS != rval)
        return error(rval);
      hdf_type = H5Tcopy( hdf_type );
      re_read_default = true;
    }
    else if (MB_TAG_NOT_FOUND == rval) {
      hdf_type = 0;
    }
    else
      return error(rval);
      
    if (hdf_type < 0)
      return error(MB_FAILURE);
  }
  else {
    switch (info.type)
    {
      case mhdf_INTEGER:
        hdf_type = H5T_NATIVE_INT;
        mb_type = MB_TYPE_INTEGER;
        break;

      case mhdf_FLOAT:
        hdf_type = H5T_NATIVE_DOUBLE;
        mb_type = MB_TYPE_DOUBLE;
        break;

      case mhdf_BOOLEAN:
        hdf_type = H5T_NATIVE_UINT;
        mb_type = MB_TYPE_INTEGER;
        break;

      case mhdf_ENTITY_ID:
        hdf_type = handleType;
        mb_type = MB_TYPE_HANDLE;
        break;

      default:
        return error(MB_FAILURE);
    }
    
    if (info.size > 1) { // array
        hsize_t tmpsize = info.size;
#if defined(H5Tarray_create_vers) && H5Tarray_create_vers > 1  
        hdf_type = H5Tarray_create2( hdf_type, 1, &tmpsize );
#else
        hdf_type = H5Tarray_create( hdf_type, 1, &tmpsize, NULL );
#endif
    }
    else {
      hdf_type = H5Tcopy( hdf_type );
    }
    if (hdf_type < 0)
      return error(MB_FAILURE);
  }

  
    // If default or global/mesh value in file, read it.
  if (info.default_value || info.global_value)
  {
    if (re_read_default) {
      mhdf_getTagValues( filePtr, info.name, hdf_type, info.default_value, info.global_value, &status );
      if (mhdf_isError( &status ))
      {
        readUtil->report_error( "%s", mhdf_message( &status ) );
        if (hdf_type) H5Tclose( hdf_type );
        return error(MB_FAILURE);
      }
    }
    
    if (MB_TYPE_HANDLE == mb_type) {
      if (info.default_value) {
        rval = convert_id_to_handle( (EntityHandle*)info.default_value, info.default_value_size );
        if (MB_SUCCESS != rval) {
          if (hdf_type) H5Tclose( hdf_type );
          return error(rval);
        }
      }
      if (info.global_value) {
        rval = convert_id_to_handle( (EntityHandle*)info.global_value, info.global_value_size );
        if (MB_SUCCESS != rval) {
          if (hdf_type) H5Tclose( hdf_type );
          return error(rval);
        }
      }
    }
  }
  
  
    // Check if tag already exists
  rval = iFace->tag_get_handle( info.name, handle );
  if (MB_SUCCESS == rval) {
    // If tag exists, make sure it is consistant with the type in the file
    int curr_size;
    DataType curr_type;
    TagType curr_store;
    
    rval = iFace->tag_get_size( handle, curr_size );
    if (MB_VARIABLE_DATA_LENGTH == rval)
      curr_size = -1;
    else if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
    
    rval = iFace->tag_get_data_type( handle, curr_type );
    if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
    
    rval = iFace->tag_get_type( handle, curr_store );
    if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
    
    if ((curr_store != MB_TAG_BIT && curr_size != info.bytes) || curr_type != mb_type ||
        ((curr_store == MB_TAG_BIT || storage == MB_TAG_BIT) && 
          curr_store != storage))
    {
      readUtil->report_error( "Tag type in file does not match type in "
                              "database for \"%s\"\n", info.name );
      if (hdf_type) H5Tclose( hdf_type );
      return error(MB_FAILURE);
    }
  }
    // Create the tag if it doesn't exist
  else if (MB_TAG_NOT_FOUND == rval)
  {
    if (info.size < 0) {
      size_t size = hdf_type ? H5Tget_size(hdf_type) : 1;
      rval = iFace->tag_create_variable_length( info.name, storage, mb_type,
                                                handle, info.default_value, 
                                                info.default_value_size * size );
    }
    else
      rval = iFace->tag_create( info.name, info.bytes, storage, mb_type,
                                handle, info.default_value );
    if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
  }
    // error
  else {
    if (hdf_type) H5Tclose( hdf_type );
    return error(rval);
  }
    
  if (info.global_value) {
    int type_size = hdf_type ? H5Tget_size(hdf_type) : 1;
    int tag_size = info.global_value_size * type_size;
    rval = iFace->tag_set_data( handle, 0, 0, &info.global_value, &tag_size );
    if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_dense_tag( Tag tag_handle,
                                    hid_t hdf_read_type,
                                    hid_t data,
                                    long start_id,
                                    long num_values )
{
  ErrorCode rval;
  DataType mb_type;
  
  rval = iFace->tag_get_data_type( tag_handle, mb_type );
  if (MB_SUCCESS != rval) 
    return error(rval);

  
  int read_size;
  rval = iFace->tag_get_size( tag_handle, read_size );
  if (MB_SUCCESS != rval) // wrong function for variable-length tags
    return error(rval);
  if (MB_TYPE_BIT == mb_type) 
    read_size = (read_size + 7)/8; // convert bits to bytes, plus 7 for ceiling
    
  if (hdf_read_type) { // if not opaque
    hsize_t hdf_size = H5Tget_size( hdf_read_type );
    if (hdf_size != (hsize_t)read_size) 
      return error(MB_FAILURE);
  }
  
    // get actual entities read from file
  Range file_ids, handles;
  Range::iterator f_ins = file_ids.begin(), h_ins = handles.begin();
  IDMap::iterator l, u;
  l = idMap.lower_bound( start_id );
  u = idMap.lower_bound( start_id + num_values - 1 );
  if (l != idMap.end() && start_id + num_values > l->begin) {
    if (l == u) {
      size_t beg = std::max(start_id, l->begin);
      size_t end = std::min(start_id + num_values, u->begin + u->count) - 1;
      f_ins = file_ids.insert( f_ins, beg, end );
      h_ins = handles.insert( h_ins, l->value + (beg - l->begin),
                                     l->value + (end - l->begin) );
    }
    else {
      size_t beg = std::max(start_id, l->begin);
      f_ins = file_ids.insert( f_ins, beg, l->begin + l->count - 1 );
      h_ins = handles.insert( h_ins, l->value + (beg - l->begin), l->value + l->count - 1 );
      for (++l; l != u; ++l) {
        f_ins = file_ids.insert( f_ins, l->begin, l->begin + l->count - 1 );
        h_ins = handles.insert( h_ins, l->value, l->value + l->count - 1 );
      }
      if (u != idMap.end() && u->begin < start_id + num_values) {
        size_t end = std::min( start_id + num_values, u->begin + u->count - 1 );
        f_ins = file_ids.insert( f_ins, u->begin, end );
        h_ins = handles.insert( h_ins, u->value, u->value + end - u->begin );
      }
    }
  }
  
    // Given that all of the entities for this dense tag data should
    // have been created as a single contiguous block, the resulting
    // MOAB handle range should be contiguous. 
  assert( handles.empty() || handles.size() == (handles.back() - handles.front() + 1));
  
  try {
    h_ins = handles.begin();
    ReadHDF5Dataset reader( data, hdf_read_type, false );
    reader.set_file_ids( file_ids, start_id, indepIO, mpiComm );
    dbgOut.tprintf(1,"Reading dense tag data with mode = %s\n", reader.get_mode_str());
    long buffer_size = bufferSize / read_size;
    while (!reader.done()) {
      size_t count;
      Range::const_iterator iter = reader.next_file_id();
      reader.read( dataBuffer, buffer_size, count );

      if (MB_TYPE_HANDLE == mb_type) {
        rval = convert_id_to_handle( (EntityHandle*)dataBuffer, count * read_size / sizeof(EntityHandle) );
        if (MB_SUCCESS != rval)
          return error(rval);
      }

      Range ents;
      Range::iterator end = h_ins;
      end += count;
      ents.insert( h_ins, end );
      h_ins = end;

      rval = iFace->tag_set_data( tag_handle, ents, dataBuffer );
      if (MB_SUCCESS != rval)
        return error(MB_FAILURE);
    }
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
    
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_sparse_tag( Tag tag_handle,
                                     hid_t hdf_read_type,
                                     hid_t id_table,
                                     hid_t value_table,
                                     long num_values )
{
  ErrorCode rval;
  DataType mbtype;
  
  rval = iFace->tag_get_data_type( tag_handle, mbtype );
  if (MB_SUCCESS != rval) 
    return error(rval);
  
  int read_size;
  rval = iFace->tag_get_size( tag_handle, read_size );
  if (MB_SUCCESS != rval) // wrong function for variable-length tags
    return error(rval);
  if (MB_TYPE_BIT == mbtype) 
    read_size = (read_size + 7)/8; // convert bits to bytes, plus 7 for ceiling
    
  if (hdf_read_type) { // if not opaque
    hsize_t hdf_size = H5Tget_size( hdf_read_type );
    if (hdf_size != (hsize_t)read_size) 
      return error(MB_FAILURE);
  }

  const int handles_per_tag = read_size/sizeof(EntityHandle);
  
  
    // Split buffer into two portions: one for handles and one for data.
    
    // We want to read the same number of handles as data values for each
    // iteration.  Calculate the total number of entries to read in each
    // pass as the size of the buffer over the sum of the size of a handle
    // and value.  Subtract off the size of one value so we reserve space
    // for adjusting for data alignment.
  size_t chunk_size = (bufferSize - read_size) / (2*sizeof(EntityHandle) + read_size);
  
    // Allocate a portion of the buffer to hold IDs read from the file
  EntityHandle* idbuf = (EntityHandle*)dataBuffer;
    // Allocate another portion of the buffer for handles of entities for which 
    // data will be read
  EntityHandle* handles = idbuf + chunk_size;
    // Use the latter portion of the buffer for data
  char* databuf = dataBuffer + (2 * chunk_size * sizeof(EntityHandle));
    // To be safe, align tag data to the size of an entire tag value
  if ((size_t)databuf % read_size)
    databuf += read_size - ((size_t)databuf % read_size);
      // Make sure the above calculations are correct
  assert( databuf + chunk_size*read_size < dataBuffer + bufferSize );

  try {
    ReadHDF5Dataset id_reader( id_table, handleType, false );
    id_reader.set_all_file_ids( collIO, mpiComm );
    ReadHDF5Dataset val_reader( value_table, hdf_read_type, false );
    dbgOut.tprintf(1,"Reading sparse tag data with ID mode = %s\n", id_reader.get_mode_str());
    ReadHDF5Dataset::Mode valmode = (ReadHDF5Dataset::Mode)-1;

      // loop until we've read the entire ID table
    size_t offset = 0;
    size_t handle_count = 0;
    Range offsets; // NOTE: must store 1-based offsets because can't put zero in moab::range
    Range::iterator ins = offsets.begin();
    while (!id_reader.done())
    {
      size_t count;
      id_reader.read( idbuf, chunk_size, count ); 

      rval = convert_id_to_handle( idbuf, count );
      if (MB_SUCCESS != rval)
        return error(rval);

        // idbuf will now contain zero-valued handles for those
        // tag values that correspond to entities we are not reading
        // from the file.
      for (size_t i = 0; i < count; ++i) {
        if (idbuf[i]) {
            // if we've filled the handle buffer, flush it
          if (handle_count == chunk_size) {
            assert( offsets.size() == handle_count );

            val_reader.set_file_ids( offsets, 1, indepIO, mpiComm );
            if (val_reader.get_mode() != valmode) {
              valmode = val_reader.get_mode();
              dbgOut.tprintf(1,"Sparse tag value reader mode set to '%s'\n", val_reader.get_mode_str());
            }
            // should normally only have one read call, unless sparse nature
            // of file_ids caused reader to do something strange
            char* iter = databuf;
            size_t remaining = handle_count;
            while (!val_reader.done()) {
              size_t count2;
              val_reader.read( iter, remaining, count2 );
              iter += count2 * read_size;
              remaining -= count2;
            }
            assert(!remaining);
            assert(iter - databuf == read_size * (ptrdiff_t)handle_count);

            if (MB_TYPE_HANDLE == mbtype) {
              rval = convert_id_to_handle( (EntityHandle*)databuf, handle_count*handles_per_tag );
              if (MB_SUCCESS != rval)
                return error(rval);
            }

            rval = iFace->tag_set_data( tag_handle, handles, handle_count, databuf );
            if (MB_SUCCESS != rval)
              return error(rval);

            handle_count = 0;
            offsets.clear();
            ins = offsets.begin();
          }

          ins = offsets.insert( ins, i+offset+1 );
          handles[handle_count++] = idbuf[i];
        }
      }

      offset += count;
    }

      // flush remaining buffered handles
    if (handle_count) {
      assert( offsets.size() == handle_count );
      val_reader.set_file_ids( offsets, 1, indepIO, mpiComm );
      // should normally only have one read call, unless sparse nature
      // of file_ids caused reader to do something strange
      char* iter = databuf;
      size_t remaining = handle_count;
      while (!val_reader.done()) {
        size_t count2;
        val_reader.read( databuf, handle_count, count2 );
        iter += read_size * count2;
        remaining -= count2;
      }
      assert(!remaining);
      assert(iter - databuf == read_size * (ptrdiff_t)handle_count);

      if (MB_TYPE_HANDLE == mbtype) {
        rval = convert_id_to_handle( (EntityHandle*)databuf, handle_count*handles_per_tag );
        if (MB_SUCCESS != rval)
          return error(rval);
      }

      rval = iFace->tag_set_data( tag_handle, handles, handle_count, databuf );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_var_len_tag( Tag tag_handle,
                                      hid_t hdf_read_type,
                                      hid_t ent_table,
                                      hid_t val_table,
                                      hid_t off_table,
                                      long num_entities,
                                      long num_values )
{
  ErrorCode rval;
  DataType mbtype;
  
  rval = iFace->tag_get_data_type( tag_handle, mbtype );
  if (MB_SUCCESS != rval) 
    return error(rval);
    
    // can't do variable-length bit tags
  if (MB_TYPE_BIT == mbtype)
    return error(MB_VARIABLE_DATA_LENGTH);

    // if here, MOAB tag must be variable-length
  int mbsize;
  if (MB_VARIABLE_DATA_LENGTH != iFace->tag_get_size( tag_handle, mbsize )) {
    assert(false);
    return error(MB_VARIABLE_DATA_LENGTH);
  }
  
  int read_size;
  if (hdf_read_type) {
    hsize_t hdf_size = H5Tget_size( hdf_read_type );
    if (hdf_size < 1)
      return error(MB_FAILURE);
    read_size = hdf_size;
  }
  else {
    // opaque
    read_size = 1;
  }
  
  try {
    ReadHDF5Dataset ents( ent_table, handleType, false );
    ents.set_all_file_ids( collIO, mpiComm );
    ReadHDF5Dataset vals( val_table, hdf_read_type, false );
    ReadHDF5Dataset offs( off_table, H5T_NATIVE_LONG, false );

      // Fill rows with the ranges to be read from the offsets table.
      // Fill handles with the corresponding entity handles.
      // Note: for each contiguous range of handles, the read from 
      //  the offsets table must also include the previous value.  So
      //  include those in 'rows' also, and zero values in handles at
      //  the corresponding positions.
    Range rows; // must be 1-based because can't insert zero in moab::Range
    Range::iterator ins = rows.begin();
    std::vector<EntityHandle> handles;
    std::vector<long> offsets;

    size_t count;
    size_t buffer_size = bufferSize/sizeof(EntityHandle);
    EntityHandle* buffer = reinterpret_cast<EntityHandle*>(dataBuffer);
    size_t offset = 0;
    while (!ents.done()) {
      ents.read( buffer, buffer_size, count );

      rval = convert_id_to_handle( buffer, count );
      if (MB_SUCCESS != rval)
        return error(rval);

      for (size_t i = 0; i < count; ++i) {
        if (!buffer[i])
          continue;

        if (rows.empty() || rows.back() != offset+i) {
          handles.push_back(0);
          if (offset+i == 0)
            offsets.push_back(-1);
          else
            ins = rows.insert( ins, offset+i );
        }
        handles.push_back( buffer[i] );
        ins = rows.insert( ins, offset+i+1 );
      }

      offset += count;
    }

      // read offsets 
    size_t nrows = rows.size();
    offsets.resize( offsets.size() + nrows );
    offs.set_file_ids( rows, 1, collIO, mpiComm );
    // should normally only have one read call, unless sparse nature
    // of file_ids caused reader to do something strange
    size_t remaining = nrows;
    long* iter = &offsets[offsets.size()-nrows];
    while (!offs.done()) {
      offs.read( iter, remaining, count );
      iter += count;
      remaining -= count;
    }
    assert(!remaining);
    assert(handles.size() == offsets.size());

      // build ranges to read from data table and 
      // convert values in offsets to per-entity data lengths
    rows.clear();
    ins = rows.begin();
    size_t i = 0;
    long max_count = 0;
    assert(!handles.front());
    while (i != handles.size()) {
      size_t j;
      for (j = i+1; j != handles.size() && handles[j]; ++j);
      assert(j - i >= 2);
      ins = rows.insert( ins, offsets[i]+2, offsets[j-1]+1 );
      for (size_t k = j-1; k > i; --k) {
        offsets[k] -= offsets[k-1];
        if (offsets[k] > max_count)
          max_count = offsets[k];
      }
      i = j;
    }

      // remove zero values from handles array
    size_t total_ents = 0;
    for (i = 0; i < handles.size(); ++i) {
      if (handles[i]) {
        handles[total_ents] = handles[i];
        offsets[total_ents] = offsets[i];
        ++total_ents;
      }
    }
    handles.resize( total_ents );
    
      // convert from long to int for offsets
    assert(sizeof(long)>=sizeof(int)); // required by ANSI C
    int* counts = reinterpret_cast<int*>(&offsets[0]);
    if (sizeof(int) < sizeof(long)) { // compact
      for (size_t i = 0; i < total_ents; ++i)
        counts[i] = offsets[i];
    }        

      // if any entity has a tag value larger than the buffer, then we
      // need a bigger buffer
    std::vector<char> space;
    buffer_size = bufferSize / read_size;
    char* data_buffer = dataBuffer;
    if (buffer_size < (size_t)max_count) {
      space.resize( max_count*read_size );
      buffer_size = max_count;
      data_buffer = &space[0];
    }

      // read tag data
    vals.set_file_ids( rows, 1, indepIO, mpiComm );
    offset = 0;
    std::vector<const void*> ptrs;
    while (!vals.done()) { 
        // figure out how many tag values we can read completely
      size_t sum = 0;
      for (i = offset; i < total_ents; ++i) {
        sum += counts[i];
        if (sum > buffer_size) {
          sum -= counts[i];
          break;
        }
      }
      const size_t n = i - offset; // number of entities for which all tag data was read
      assert(n && sum); // should be at least one whole tag value
      
        // read 'sum' tag values for 'n' entities
      for (size_t j = 0; j < sum; j += count)
        vals.read( data_buffer + j*read_size, sum-j, count );

      if (MB_TYPE_HANDLE == mbtype) {
        rval = convert_id_to_handle( (EntityHandle*)data_buffer, sum );
        if (MB_SUCCESS != rval)
          return error(rval);
      }

        // get per-value pointers
        // convert counts from num values to num bytes
      ptrs.resize(n);
      ptrs[0] = data_buffer;
      counts[offset] *= read_size;
      for (i = 1; i < n; ++i) {
        ptrs[i] = (char*)ptrs[i-1] + counts[offset+i-1];
        counts[offset+i] *= read_size;
      }

        // hand off tag data to moab
      rval = iFace->tag_set_data( tag_handle, &handles[offset], n, &ptrs[0], &counts[offset] );
      if (MB_SUCCESS != rval)
        return error(rval);

      offset += n;
    }
  }
  catch (ReadHDF5Dataset::Exception e) {
    return error(MB_FAILURE);
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::convert_id_to_handle( EntityHandle* array, 
                                            size_t size )
{
  convert_id_to_handle( array, size, idMap );
  return MB_SUCCESS;
}

void ReadHDF5::convert_id_to_handle( EntityHandle* array, 
                                     size_t size,
                                     const RangeMap<long,EntityHandle>& id_map )
{
  for (EntityHandle* const end = array + size; array != end; ++array)
    *array = id_map.find( *array );
}

void ReadHDF5::convert_id_to_handle( EntityHandle* array, 
                                     size_t size, size_t& new_size,
                                     const RangeMap<long,EntityHandle>& id_map )
{
  RangeMap<long,EntityHandle>::const_iterator it;
  new_size = 0;
  for (size_t i = 0; i < size; ++i) {
    it = id_map.lower_bound( array[i] );
    if (it != id_map.end() && it->begin <= (long)array[i])
      array[new_size++] = it->value + (array[i] - it->begin);
  }
}

void ReadHDF5::convert_range_to_handle( const EntityHandle* ranges,
                                        size_t num_ranges,
                                        const RangeMap<long,EntityHandle>& id_map,
                                        Range& merge )
{
  RangeMap<long,EntityHandle>::iterator it = id_map.begin();
  for (size_t i = 0; i < num_ranges; ++i) {
    long id = ranges[2*i];
    const long end = id + ranges[2*i+1];
      // we assume that 'ranges' is sorted, but check just in case it isn't.
    if (it == id_map.end() || it->begin > id)
      it = id_map.begin();
    it = id_map.lower_bound( it, id_map.end(), id );
    if (it == id_map.end())
      continue;
    if (id < it->begin)
      id = it->begin;
    while (id < end) {
      const long off = id - it->begin;
      long count = std::min( it->count - off,  end - id );
      merge.insert( it->value + off, it->value + off + count - 1 );
      id += count;
      if (id < end)
        if (++it == id_map.end())
          break;
    }
  }
}

ErrorCode ReadHDF5::convert_range_to_handle( const EntityHandle* array,
                                             size_t num_ranges,
                                             Range& range )
{
  convert_range_to_handle( array, num_ranges, idMap, range );
  return MB_SUCCESS;
}
  

ErrorCode ReadHDF5::insert_in_id_map( const Range& file_ids,
                                      EntityHandle start_id )
{
  Range::const_pair_iterator p;
  for (p = file_ids.const_pair_begin(); p != file_ids.const_pair_end(); ++p) {
    size_t count = p->second - p->first + 1;
    if (!idMap.insert( p->first, start_id, count ).second) 
      return error(MB_FAILURE);
    start_id += count;
  }
  return MB_SUCCESS;
}


ErrorCode ReadHDF5::read_qa( EntityHandle  )
{
  mhdf_Status status;
  std::vector<std::string> qa_list;
  
  int qa_len;
  char** qa = mhdf_readHistory( filePtr, &qa_len, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( "%s", mhdf_message( &status ) );
    return error(MB_FAILURE);
  }
  qa_list.resize(qa_len);
  for (int i = 0; i < qa_len; i++)
  {
    qa_list[i] = qa[i];
    free( qa[i] );
  }
  free( qa );
  
  /** FIX ME - how to put QA list on set?? */

  return MB_SUCCESS;
}

ErrorCode ReadHDF5::store_file_ids( Tag tag )
{
  typedef int tag_type;
  tag_type* buffer = reinterpret_cast<tag_type*>(dataBuffer);
  const long buffer_size = bufferSize / sizeof(tag_type);
  for (IDMap::iterator i = idMap.begin(); i != idMap.end(); ++i) {
    IDMap::Range range = *i;
    
      // make sure the values will fit in the tag type
    IDMap::key_type rv = range.begin + (range.count - 1);
    tag_type tv = (tag_type)rv;
    if ((IDMap::key_type)tv != rv) {
      assert(false);
      return MB_INDEX_OUT_OF_RANGE;
    }
    
    while (range.count) {
      long count = buffer_size < range.count ? buffer_size : range.count;

      Range handles;
      handles.insert( range.value, range.value + count - 1 );
      range.value += count;
      range.count -= count;
      for (long j = 0; j < count; ++j) 
        buffer[j] = (tag_type)range.begin++;

      ErrorCode rval = iFace->tag_set_data( tag, handles, buffer );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_tag_values( const char* file_name,
                                     const char* tag_name,
                                     const FileOptions& opts,
                                     std::vector<int>& tag_values_out,
                                     const SubsetList* subset_list )
{
  ErrorCode rval;
  
  rval = set_up_read( file_name, opts );
  if (MB_SUCCESS != rval)
    return error(rval);
  
  int tag_index;
  rval = find_int_tag( tag_name, tag_index );
  if (MB_SUCCESS != rval) {
    clean_up_read( opts );
    return error(rval);
  }
  
  if (subset_list) {
    Range file_ids;
    rval = get_subset_ids( subset_list->tag_list, subset_list->tag_list_length, file_ids );
    if (MB_SUCCESS != rval) {
      clean_up_read( opts );
      return error(rval);
    }
    
    rval = read_tag_values_partial( tag_index, file_ids, tag_values_out );
    if (MB_SUCCESS != rval) {
      clean_up_read( opts );
      return error(rval);
    }
  }
  else {
    rval = read_tag_values_all( tag_index, tag_values_out );
    if (MB_SUCCESS != rval) {
      clean_up_read( opts );
      return error(rval);
    }
  }
    
  return clean_up_read( opts );
}

ErrorCode ReadHDF5::read_tag_values_partial( int tag_index,
                                             const Range& file_ids,
                                             std::vector<int>& tag_values )
{
  mhdf_Status status;
  const mhdf_TagDesc& tag = fileInfo->tags[tag_index];
  long num_ent, num_val;
  size_t count;
  
    // read sparse values
  if (tag.have_sparse) {
    hid_t handles[3];
    mhdf_openSparseTagData( filePtr, tag.name, &num_ent, &num_val, handles, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      return error(MB_FAILURE);
    }
    
    try {
      ReadHDF5Dataset ids( handles[0], H5T_NATIVE_LONG );
      ids.set_all_file_ids( collIO, mpiComm );
      ReadHDF5Dataset vals( handles[1], H5T_NATIVE_INT );

        // read all entity handles and fill 'offsets' with ranges of
        // offsets into the data table for entities that we want.
      Range offsets;
      long* buffer = reinterpret_cast<long*>(dataBuffer);
      const long buffer_size = bufferSize/sizeof(long);
      size_t offset = 0;
      while (!ids.done()) {
        ids.read( buffer, buffer_size, count );

        std::sort( buffer, buffer+count );
        Range::iterator ins = offsets.begin();
        Range::const_iterator i = file_ids.begin();
        for (size_t j = 0; j < count; ++j) {
          while (i != file_ids.end() && (long)*i < buffer[j])
            ++i;
          if (i == file_ids.end())
            break;
          if ((long)*i == buffer[j]) {
            ins = offsets.insert( ins, j+offset, j+offset );
          }
        }
        
        offset += count;
      }

      vals.set_file_ids( offsets, 0, collIO, mpiComm );
      tag_values.resize( offsets.size() );
      // should normally only have one read call, unless sparse nature
      // of file_ids caused reader to do something strange
      size_t off = 0;
      while (!vals.done()) {
        vals.read( &tag_values[off], tag_values.size() - off, count );
        off += count;
      }
      assert(off == tag_values.size());
    }
    catch (ReadHDF5Dataset::Exception e) {
      return error(MB_FAILURE);
    }
  }
  
  std::sort( tag_values.begin(), tag_values.end() );
  tag_values.erase( std::unique(tag_values.begin(), tag_values.end()), tag_values.end() );
  
    // read dense values
  std::vector<int> prev_data, curr_data;
  for (int i = 0; i < tag.num_dense_indices; ++i) {
    int grp = tag.dense_elem_indices[i];
    const char* gname = 0;
    mhdf_EntDesc* desc = 0;
    if (grp == -1) {
      gname = mhdf_node_type_handle();
      desc = &fileInfo->nodes;
    }
    else if (grp == -2) {
      gname = mhdf_set_type_handle();
      desc = &fileInfo->sets;
    }
    else {
      assert(grp >= 0 && grp < fileInfo->num_elem_desc);
      gname = fileInfo->elems[grp].handle;
      desc = &fileInfo->elems[grp].desc;
    }
    
    Range::iterator s = file_ids.lower_bound( (EntityHandle)(desc->start_id) );
    Range::iterator e = Range::lower_bound( s, file_ids.end(),  
                                   (EntityHandle)(desc->start_id) + desc->count );
    Range subset;
    subset.merge( s, e );
    
    hid_t handle = mhdf_openDenseTagData( filePtr, tag.name, gname, &num_val, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      return error(MB_FAILURE);
    }
    
    try {
      ReadHDF5Dataset reader( handle, H5T_NATIVE_INT );
      reader.set_file_ids( subset, desc->start_id, indepIO, mpiComm );
      curr_data.resize( subset.size() );
      // should normally only have one read call, unless sparse nature
      // of file_ids caused reader to do something strange
      size_t off = 0;
      while (!reader.done()) {
        reader.read( &curr_data[off], curr_data.size() - off, count );
        off += count;
      }
      assert( off == curr_data.size() );
    }
    catch (ReadHDF5Dataset::Exception e) {
      return error(MB_FAILURE);
    }
    
    std::sort( curr_data.begin(), curr_data.end() );
    curr_data.erase( std::unique(curr_data.begin(), curr_data.end()), curr_data.end() );
    prev_data.clear();
    tag_values.swap( prev_data );
    std::set_union( prev_data.begin(), prev_data.end(),
                    curr_data.begin(), curr_data.end(),
                    std::back_inserter( tag_values ) );
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadHDF5::read_tag_values_all( int tag_index,
                                         std::vector<int>& tag_values )
{
  mhdf_Status status;
  const mhdf_TagDesc& tag = fileInfo->tags[tag_index];
  long junk, num_val;
  
    // read sparse values
  if (tag.have_sparse) {
    hid_t handles[3];
    mhdf_openSparseTagData( filePtr, tag.name, &junk, &num_val, handles, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      return error(MB_FAILURE);
    }
    
    mhdf_closeData( filePtr, handles[0], &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(MB_FAILURE);
    }
    
    tag_values.resize( num_val );
    mhdf_readTagValuesWithOpt( handles[1], 0, num_val, H5T_NATIVE_INT,
                               &tag_values[0], collIO, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(MB_FAILURE);
    }
    
    mhdf_closeData( filePtr, handles[1], &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      return error(MB_FAILURE);
    }
  }
  
  std::sort( tag_values.begin(), tag_values.end() );
  tag_values.erase( std::unique(tag_values.begin(), tag_values.end()), tag_values.end() );
  
    // read dense values
  std::vector<int> prev_data, curr_data;
  for (int i = 0; i < tag.num_dense_indices; ++i) {
    int grp = tag.dense_elem_indices[i];
    const char* gname = 0;
    if (grp == -1)
      gname = mhdf_node_type_handle();
    else if (grp == -2)
      gname = mhdf_set_type_handle();
    else
      gname = fileInfo->elems[grp].handle;
    hid_t handle = mhdf_openDenseTagData( filePtr, tag.name, gname, &num_val, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      return error(MB_FAILURE);
    }
    
    curr_data.resize( num_val );
    mhdf_readTagValuesWithOpt( handle, 0, num_val, H5T_NATIVE_INT, &curr_data[0], collIO, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      mhdf_closeData( filePtr, handle, &status );
      return error(MB_FAILURE);
    }
    
    mhdf_closeData( filePtr, handle, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( "%s", mhdf_message( &status ) );
      return error(MB_FAILURE);
    }
 
    std::sort( curr_data.begin(), curr_data.end() );
    curr_data.erase( std::unique(curr_data.begin(), curr_data.end()), curr_data.end() );
    
    prev_data.clear();
    tag_values.swap( prev_data );
    std::set_union( prev_data.begin(), prev_data.end(),
                    curr_data.begin(), curr_data.end(),
                    std::back_inserter( tag_values ) );
  }
  
  return MB_SUCCESS;
}

} // namespace moab
