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
// Filename      : WriteHDF5.cpp
//
// Purpose       : TSTT HDF5 Writer 
//
// Special Notes : WriteSLAC used as template for this
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/01/04
//-------------------------------------------------------------------------

#ifndef HDF5_FILE
#  error Attempt to compile WriteHDF5 with HDF5 support disabled
#endif

#include <assert.h>
#if defined(_MSC_VER) || defined(__MINGW32__)
#include <sys/time.h>
#endif
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits>
#include <cstdio>
#include <iostream>
/* include our MPI header before any HDF5 because otherwise
   it will get included indirectly by HDF5 */
#ifdef USE_MPI
#  include "moab_mpi.h"
#endif 
#include <H5Tpublic.h>
#include <H5Ppublic.h>
#include "moab/Interface.hpp"
#include "Internals.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"
#include "WriteHDF5.hpp"
#include "moab/WriteUtilIface.hpp"
#include "FileOptions.hpp"
#include "moab/Version.h"
#include "IODebugTrack.hpp"
#include "mhdf.h"
/* Access HDF5 file handle for debugging
#include <H5Fpublic.h>
struct file { uint32_t magic; hid_t handle; };
*/
#undef DEBUG

#ifdef DEBUG
/*
# include <H5Epublic.h>
  extern "C" herr_t hdf_error_handler( void*  )
  {
    H5Eprint( stderr );
    assert( 0 );
  }
*/
# define myassert(A) assert(A)
#else
# define myassert(A)
#endif


#ifdef VALGRIND
#  include <valgrind/memcheck.h>
#else
#  ifndef VALGRIND_CHECK_MEM_IS_DEFINED
#    define VALGRIND_CHECK_MEM_IS_DEFINED(a, b)
#  endif
#  ifndef VALGRIND_CHECK_MEM_IS_ADDRESSABLE
#    define VALGRIND_CHECK_MEM_IS_ADDRESSABLE(a, b)
#  endif
#  ifndef VALGRIND_MAKE_MEM_UNDEFINED
#    define VALGRIND_MAKE_MEM_UNDEFINED(a, b)
#  endif
#endif

namespace moab {

template <typename T> inline 
void VALGRIND_MAKE_VEC_UNDEFINED( std::vector<T>& v ) {
    VALGRIND_MAKE_MEM_UNDEFINED( &v[0], v.size() * sizeof(T) );
}

#define WRITE_HDF5_BUFFER_SIZE (40*1024*1024)

static hid_t get_id_type()
{
  if (8 == sizeof(WriteHDF5::id_t)) {
    if (8 == sizeof(long))  
      return H5T_NATIVE_ULONG;
    else 
      return H5T_NATIVE_UINT64;
  }
  else if (4 == sizeof(WriteHDF5::id_t)) {
    if (4 == sizeof(int))
      return H5T_NATIVE_UINT;
    else
      return H5T_NATIVE_UINT32;
  }
  else {
    assert(0);
    return (hid_t)-1;
  }
}

  // This is the HDF5 type used to store file IDs
const hid_t WriteHDF5::id_type = get_id_type();


// This function doesn't do anything useful.  It's just a nice
// place to set a break point to determine why the reader fails.
static inline ErrorCode error( ErrorCode rval )
  { return rval; }


  // Some macros to handle error checking.  The
  // CHK_MHDF__ERR* macros check the value of an mhdf_Status 
  // object.  The CHK_MB_ERR_* check the value of an ErrorCode.
  // The *_0 macros accept no other arguments. The *_1
  // macros accept a single hdf5 handle to close on error.
  // The *_2 macros accept an array of two hdf5 handles to
  // close on error.  The _*2C macros accept one hdf5 handle
  // to close on error and a bool and an hdf5 handle where
  // the latter handle is conditionally closed depending on
  // the value of the bool.  All macros contain a "return"
  // statement.
#define CHK_MHDF_ERR_0( A )                                 \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    return error(MB_FAILURE);                               \
} while(false)                                               

#define CHK_MHDF_ERR_1( A, B )                              \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    mhdf_closeData( filePtr, (B), &(A) );                   \
    return error(MB_FAILURE);                               \
} while(false)                                               

#define CHK_MHDF_ERR_2( A, B )                              \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    mhdf_closeData( filePtr, (B)[0], &(A) );                \
    mhdf_closeData( filePtr, (B)[1], &(A) );                \
    return error(MB_FAILURE);                               \
} while(false)                                               

#define CHK_MHDF_ERR_3( A, B )                              \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    mhdf_closeData( filePtr, (B)[0], &(A) );                \
    mhdf_closeData( filePtr, (B)[1], &(A) );                \
    mhdf_closeData( filePtr, (B)[2], &(A) );                \
    return error(MB_FAILURE);                               \
} while(false)                                               

#define CHK_MHDF_ERR_2C( A, B, C, D )                       \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    mhdf_closeData( filePtr, (B), &(A) );                   \
    if (C) mhdf_closeData( filePtr, (D), &(A) );            \
    return error(MB_FAILURE);                               \
} while(false)                                               


#define CHK_MB_ERR_0( A ) \
do if (MB_SUCCESS != (A)) return error(A); while(false)

#define CHK_MB_ERR_1( A, B, C )         \
do if (MB_SUCCESS != (A)) {             \
  mhdf_closeData( filePtr, (B), &(C) ); \
  myassert(0);                          \
  return error(A);                      \
} while(false)

#define CHK_MB_ERR_2( A, B, C )            \
do if (MB_SUCCESS != (A)) {                \
  mhdf_closeData( filePtr, (B)[0], &(C) ); \
  mhdf_closeData( filePtr, (B)[1], &(C) ); \
  write_finished();                        \
  myassert(0);                             \
  return error(A);                         \
} while(false)

#define CHK_MB_ERR_3( A, B, C )            \
do if (MB_SUCCESS != (A)) {                \
  mhdf_closeData( filePtr, (B)[0], &(C) ); \
  mhdf_closeData( filePtr, (B)[1], &(C) ); \
  mhdf_closeData( filePtr, (B)[2], &(C) ); \
  write_finished();                        \
  myassert(0);                             \
  return error(A);                         \
} while(false)

#define CHK_MB_ERR_2C( A, B, C, D, E )          \
do if (MB_SUCCESS != (A)) {                     \
  mhdf_closeData( filePtr, (B), &(E) );         \
  if (C) mhdf_closeData( filePtr, (D), &(E) );  \
  write_finished();                             \
  myassert(0);                                  \
  return error(A);                              \
} while(false)


#define debug_barrier() debug_barrier_line(__LINE__)
void WriteHDF5::debug_barrier_line(int )
{
}

bool WriteHDF5::convert_handle_tag( const EntityHandle* source,
                                    EntityHandle* dest, size_t count ) const
{
  bool some_valid = false;
  for (size_t i = 0; i < count; ++i) {
    if (!source[i])
      dest[i] = 0;
    else {
      dest[i] = idMap.find( source[i] );
      if (dest[i])
        some_valid = true;
    }
  }
  return some_valid;
}

bool WriteHDF5::convert_handle_tag( EntityHandle* data, size_t count ) const
{
  assert( sizeof(EntityHandle) == sizeof(id_t) );
  return convert_handle_tag( data, data, count );
}

ErrorCode WriteHDF5::assign_ids( const Range& entities, id_t id )
{
  Range::const_pair_iterator pi;
  for (pi = entities.const_pair_begin(); pi != entities.const_pair_end(); ++pi) {
    const EntityHandle n = pi->second - pi->first + 1;
    dbgOut.printf( 3, "Assigning %s %lu to %lu to file IDs [%lu,%lu]\n",
      CN::EntityTypeName(TYPE_FROM_HANDLE(pi->first)),
      (unsigned long)(ID_FROM_HANDLE(pi->first)),
      (unsigned long)(ID_FROM_HANDLE(pi->first)+n-1),
      (unsigned long)id,
      (unsigned long)(id+n-1));
    if (!idMap.insert( pi->first, id, n ).second)
      return error(MB_FAILURE);
    id += n;
  }
  return MB_SUCCESS;
}

const char* WriteHDF5::ExportSet::name() const
{
  static char buffer[32];
  switch (type) {
    case MBVERTEX:
      return mhdf_node_type_handle();
    case MBENTITYSET:
      return mhdf_set_type_handle();
    default:
      sprintf( buffer, "%s%d", CN::EntityTypeName( type ), num_nodes );
      return buffer;
  }
}
  

WriterIface* WriteHDF5::factory( Interface* iface )
  { return new WriteHDF5( iface ); }

WriteHDF5::WriteHDF5( Interface* iface )
  : bufferSize( WRITE_HDF5_BUFFER_SIZE ),
    dataBuffer( 0 ),
    iFace( iface ), 
    writeUtil( 0 ), 
    filePtr( 0 ), 
    setContentsOffset( 0 ),
    setChildrenOffset( 0 ),
    setParentsOffset( 0 ),
    maxNumSetContent( 0 ),
    maxNumSetChildren( 0 ),
    maxMumSetParents( 0 ),
    writeSets(false),
    writeSetContents(false),
    writeSetChildren(false),
    writeSetParents(false),
    parallelWrite(false),
    collectiveIO(false),
    writeTagDense(false),
    writeProp( H5P_DEFAULT ),
    dbgOut("H5M ", stderr)
{
}

ErrorCode WriteHDF5::init()
{
  ErrorCode rval;

  if (writeUtil) // init has already been called
    return MB_SUCCESS;
/* 
#ifdef DEBUG
  H5Eset_auto( &hdf_error_handler, writeUtil );  // HDF5 callback for errors
#endif
*/ 
    // For known tag types, store the corresponding HDF5 in which
    // the tag data is to be written in the file.
  //register_known_tag_types( iFace ); 
 
    // Get the util interface
  void* ptr = 0;
  rval = iFace->query_interface( "WriteUtilIface", &ptr );
  CHK_MB_ERR_0(rval);

  idMap.clear();
  
  if (MB_SUCCESS != rval)
  {
    writeUtil = 0;
    return error(rval);
  }

  writeUtil = reinterpret_cast<WriteUtilIface*>(ptr);
  return MB_SUCCESS;
}
  
ErrorCode WriteHDF5::write_finished()
{
    // release memory allocated in lists
  exportList.clear();
  nodeSet.range.clear();
  setSet.range.clear();
  tagList.clear();
  idMap.clear();
  return MB_SUCCESS;
}

WriteHDF5::~WriteHDF5()
{
  if (!writeUtil) // init() failed.
    return;

  iFace->release_interface( "WriteUtilIface", writeUtil );
}


ErrorCode WriteHDF5::write_file( const char* filename,
                                   bool overwrite,
                                   const FileOptions& opts,
                                   const EntityHandle* set_array,
                                   const int num_sets,
                                   const std::vector<std::string>& qa_records,
                                   const Tag* tag_list,
                                   int num_tags,
                                   int user_dimension )
{
  mhdf_Status status;
  
  parallelWrite = false;
  collectiveIO = false;

  // Enable debug output
  int tmpval = 0;
  if (MB_SUCCESS == opts.get_int_option("DEBUG_IO", 1, tmpval)) {
    dbgOut.set_verbosity(tmpval);
  }
  
  writeTagDense = (MB_SUCCESS == opts.get_null_option("DENSE_TAGS"));
    

  // Enable some extra checks for reads.  Note: amongst other things this
  // will print errors if the entire file is not read, so if doing a 
  // partial read that is not a parallel read, this should be disabled.
  debugTrack = (MB_SUCCESS == opts.get_null_option("DEBUG_BINIO"));
    
  bufferSize = WRITE_HDF5_BUFFER_SIZE;
  int buf_size;
  ErrorCode rval = opts.get_int_option( "BUFFER_SIZE", buf_size );
  if (MB_SUCCESS == rval && buf_size >= 24)
    bufferSize = buf_size;

    // Allocate internal buffer to use when gathering data to write.
  dataBuffer = (char*)malloc( bufferSize );
  if (!dataBuffer)
    return error(MB_MEMORY_ALLOCATION_FAILED);

    // Clear filePtr so we know if it is open upon failure
  filePtr = 0;

    // Do actual write.
  writeProp = H5P_DEFAULT;
  ErrorCode result = write_file_impl( filename, overwrite, opts, 
                                        set_array, num_sets, 
                                        qa_records, 
                                        tag_list, num_tags,
                                        user_dimension );
    // close writeProp if it was opened
  if (writeProp != H5P_DEFAULT)
    H5Pclose(writeProp);
  
    // Free memory buffer
  free( dataBuffer );
  dataBuffer = 0;
  
    // Close file
  bool created_file = false;
  if (filePtr) {
    created_file = true;
    mhdf_closeFile( filePtr, &status );
    filePtr = 0;
    if (mhdf_isError( &status )) {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      if (MB_SUCCESS == result)
        result = MB_FAILURE;
    }
  }

    // Release other resources
  if (MB_SUCCESS == result)
    result = write_finished();
  else
    write_finished();
  
    // If write failed, remove file unless KEEP option was specified
  if (MB_SUCCESS != result && created_file && 
      MB_ENTITY_NOT_FOUND == opts.get_null_option( "KEEP" ))
    remove( filename );
  
  return result;
}  


ErrorCode WriteHDF5::write_file_impl( const char* filename,
                                        bool overwrite,
                                        const FileOptions& opts,
                                        const EntityHandle* set_array,
                                        const int num_sets,
                                        const std::vector<std::string>& qa_records,
                                        const Tag* tag_list, 
                                        int num_tags,
                                        int user_dimension )
{
  ErrorCode result;
  std::list<SparseTag>::const_iterator t_itor;
  std::list<ExportSet>::iterator ex_itor;
  EntityHandle elem_count, max_id;
  
  if (MB_SUCCESS != init())
    return error(MB_FAILURE);

  dbgOut.tprint(1,"Gathering Mesh\n");
  
    // Gather mesh to export
  exportList.clear();
  if (0 == num_sets || (1 == num_sets && set_array[0] == 0))
  {
    result = gather_all_mesh( );
    CHK_MB_ERR_0(result);
  }
  else
  {
    std::vector<EntityHandle> passed_export_list(set_array, set_array+num_sets);
    result = gather_mesh_info( passed_export_list );
    if (MB_SUCCESS != result) 
      return error(result);
  }
  
  //if (nodeSet.range.size() == 0)
  //  return error(MB_ENTITY_NOT_FOUND);
  
  dbgOut.tprint(1,"Checking ID space\n");

    // Make sure ID space is sufficient
  elem_count = nodeSet.range.size() + setSet.range.size();
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
    elem_count += ex_itor->range.size();
  max_id = (EntityHandle)1 << (8*sizeof(id_t)-1);
  if (elem_count > max_id)
  {
    writeUtil->report_error("ID space insufficient for mesh size.\n");
    return error(result);
  }

  dbgOut.tprint(1, "Creating File\n" );  

    // Figure out the dimension in which to write the mesh.  
  int mesh_dim;
  result = iFace->get_dimension( mesh_dim );
  CHK_MB_ERR_0(result);
  
  if (user_dimension < 1) 
    user_dimension = mesh_dim;
  user_dimension = user_dimension > mesh_dim ? mesh_dim : user_dimension;
  
    // Create the file layout, including all tables (zero-ed) and
    // all structure and meta information.
  const char* optnames[] = { "WRITE_PART", "FORMAT", 0 };
  int junk;
  parallelWrite = (MB_SUCCESS == opts.match_option( "PARALLEL", optnames, junk ));
  if (parallelWrite) {
    int pcomm_no = 0;
    opts.get_int_option("PARALLEL_COMM", pcomm_no);
      // Just store Boolean value based on string option here.
      // parallel_create_file will set writeProp accordingly.
    collectiveIO =  (MB_SUCCESS == opts.get_null_option("COLLECTIVE"));
    dbgOut.printf(2,"'COLLECTIVE' option = %s\n", collectiveIO ? "YES" : "NO" );
    result = parallel_create_file( filename, overwrite, qa_records, tag_list, num_tags, user_dimension, pcomm_no );
  }
  else {
    result = serial_create_file( filename, overwrite, qa_records, tag_list, num_tags, user_dimension );
  }
  if (MB_SUCCESS != result)
    return error(result);

  dbgOut.tprint(1,"Writing Nodes.\n");
  
    // Write node coordinates
  if (!nodeSet.range.empty() || parallelWrite) {
    result = write_nodes();
    if (MB_SUCCESS != result)
      return error(result);
  }

  dbgOut.tprint(1,"Writing connectivity.\n");
  
    // Write element connectivity
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor) {
    result = write_elems( *ex_itor );
    if (MB_SUCCESS != result)
      return error(result);
  }

  dbgOut.tprint(1,"Writing sets.\n");
  
    // Write meshsets
  result = write_sets();
  if (MB_SUCCESS != result)
    return error(result);

  dbgOut.tprint(1,"Writing adjacencies.\n");
  
    // Write adjacencies
  // Tim says don't save node adjacencies!
#ifdef MB_H5M_WRITE_NODE_ADJACENCIES
  result = write_adjacencies( nodeSet );
  if (MB_SUCCESS != result)
    return error(result);
#endif
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor) {
    result = write_adjacencies( *ex_itor );
    if (MB_SUCCESS != result)
      return error(result);
  }

  dbgOut.tprint(1,"Writing tags.\n");
  

    // Write tags
  for (t_itor = tagList.begin(); t_itor != tagList.end(); ++t_itor) {
    result = write_tag( *t_itor );
    if (MB_SUCCESS != result)
      return error(result);
  }
  
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::initialize_mesh( const Range ranges[5] )
{
  ErrorCode rval;
  
  if (!ranges[0].all_of_type(MBVERTEX))
    return error(MB_FAILURE);
  nodeSet.range = ranges[0];
  nodeSet.type = MBVERTEX;
  nodeSet.num_nodes = 1;
  nodeSet.max_num_ents = nodeSet.max_num_adjs = 0;
  
  if (!ranges[4].all_of_type(MBENTITYSET))
    return error(MB_FAILURE);
  setSet.range = ranges[4];
  setSet.type = MBENTITYSET;
  setSet.num_nodes = 0;
  setSet.max_num_ents = setSet.max_num_adjs = 0;
  maxNumSetContent = maxNumSetChildren = maxMumSetParents = 0;

  exportList.clear();
  std::vector<Range> bins(1024); // sort entities by connectivity length
                                   // resize is expensive due to Range copy, so start big
  for (EntityType type = MBEDGE; type < MBENTITYSET; ++type)
  {
    ExportSet set;
    set.max_num_ents = set.max_num_adjs = 0;
    const int dim = CN::Dimension(type);

      // Group entities by connectivity length
    bins.clear();
    assert(dim >= 0 && dim <= 4);
    std::pair<Range::const_iterator,Range::const_iterator> p = ranges[dim].equal_range(type);
    Range::const_iterator i = p.first;
    while (i != p.second) {
      Range::const_iterator first = i;
      EntityHandle const* conn;
      int len, firstlen;
      rval = iFace->get_connectivity( *i, conn, firstlen );
      if (MB_SUCCESS != rval)
        return error(rval);
      
      for (++i; i != p.second; ++i) {
        rval = iFace->get_connectivity( *i, conn, len );
        if (MB_SUCCESS != rval)
          return error(rval);
        
        if (len != firstlen)
          break;
      }
      
      if (firstlen >= (int)bins.size())
        bins.resize(firstlen+1);
      bins[firstlen].merge( first, i );
    }

      // Create ExportSet for each group
    for (std::vector<Range>::iterator j = bins.begin(); j != bins.end(); ++j) {
      if (j->empty())
        continue;
        
      set.range.clear();
      set.type = type;
      set.num_nodes = j - bins.begin();
      exportList.push_back( set );
      exportList.back().range.swap( *j );
    }
  }
    
  return MB_SUCCESS;  
}

                                         
  // Gather the mesh to be written from a list of owning meshsets.
ErrorCode WriteHDF5::gather_mesh_info( 
                           const std::vector<EntityHandle>& export_sets )
{
  ErrorCode rval;
  
  int dim;
  Range range;      // temporary storage
  Range ranges[5];  // lists of entities to export, grouped by dimension
  
    // Gather list of all related sets
  std::vector<EntityHandle> stack(export_sets);
  std::copy( export_sets.begin(), export_sets.end(), stack.begin() );
  std::vector<EntityHandle> set_children;
  while( !stack.empty() )
  {
    EntityHandle meshset = stack.back(); stack.pop_back();
    ranges[4].insert( meshset );
  
      // Get contained sets
    range.clear();
    rval = iFace->get_entities_by_type( meshset, MBENTITYSET, range );
    CHK_MB_ERR_0(rval);
    for (Range::iterator ritor = range.begin(); ritor != range.end(); ++ritor)
      if (ranges[4].find( *ritor ) == ranges[4].end())
        stack.push_back( *ritor );
    
      // Get child sets
    set_children.clear();
    rval = iFace->get_child_meshsets( meshset, set_children, 1 );
    CHK_MB_ERR_0(rval);
    for (std::vector<EntityHandle>::iterator vitor = set_children.begin();
         vitor != set_children.end(); ++vitor )
      if (ranges[4].find( *vitor ) == ranges[4].end())
        stack.push_back( *vitor );
  }
  
    // Gather list of all mesh entities from list of sets,
    // grouped by dimension.
  for (Range::iterator setitor = ranges[4].begin();
       setitor != ranges[4].end(); ++setitor)
  {
    for (dim = 0; dim < 4; ++dim)
    {
      range.clear();
      rval = iFace->get_entities_by_dimension( *setitor, dim, range, false );
      CHK_MB_ERR_0(rval);

      ranges[dim].merge(range);
    }
  }
  
    // For each list of elements, append adjacent children and
    // nodes to lists.
  for (dim = 3; dim > 0; --dim)
  {
    for (int cdim = 1; cdim < dim; ++cdim)
    {
      range.clear();
      rval = iFace->get_adjacencies( ranges[dim], cdim, false, range );
      CHK_MB_ERR_0(rval);
      ranges[cdim].merge( range );
    }  
    range.clear();
    rval = writeUtil->gather_nodes_from_elements( ranges[dim], 0, range );
    CHK_MB_ERR_0(rval);
    ranges[0].merge( range );      
  }
  
  return initialize_mesh( ranges );
}

  // Gather all the mesh and related information to be written.
ErrorCode WriteHDF5::gather_all_mesh( )
{
  ErrorCode rval;
  Range ranges[5];

  rval = iFace->get_entities_by_type( 0, MBVERTEX, ranges[0] );
  if (MB_SUCCESS != rval)
    return error(rval);

  rval = iFace->get_entities_by_dimension( 0, 1, ranges[1] );
  if (MB_SUCCESS != rval)
    return error(rval);

  rval = iFace->get_entities_by_dimension( 0, 2, ranges[2] );
  if (MB_SUCCESS != rval)
    return error(rval);

  rval = iFace->get_entities_by_dimension( 0, 3, ranges[3] );
  if (MB_SUCCESS != rval)
    return error(rval);

  rval = iFace->get_entities_by_type( 0, MBENTITYSET, ranges[4] );
  if (MB_SUCCESS != rval)
    return error(rval);

  return initialize_mesh( ranges );
}
  
ErrorCode WriteHDF5::write_nodes( )
{
  mhdf_Status status;
  int dim, mesh_dim;
  ErrorCode rval;
  hid_t node_table;
  long first_id, num_nodes;
  
  rval = iFace->get_dimension( mesh_dim );
  CHK_MB_ERR_0(rval);
  
  debug_barrier();
  dbgOut.print(3, "Opening Node Coords\n");
  node_table = mhdf_openNodeCoords( filePtr, &num_nodes, &dim, &first_id, &status );
  CHK_MHDF_ERR_0(status);
  IODebugTrack track( debugTrack, "nodes", num_nodes );
  
  double* buffer = (double*)dataBuffer;
  int chunk_size = bufferSize / sizeof(double);
  
  long remaining = nodeSet.range.size();
  long num_writes = (remaining+chunk_size-1) / chunk_size;
  if (nodeSet.max_num_ents) {
    assert( nodeSet.max_num_ents >= remaining );
    num_writes = (nodeSet.max_num_ents+chunk_size-1) / chunk_size;
  }
  long remaining_writes = num_writes;

  long offset = nodeSet.offset;
  Range::const_iterator iter = nodeSet.range.begin();
  dbgOut.printf(3, "Writing %ld nodes in %ld blocks of %d\n", remaining, (remaining+chunk_size-1)/chunk_size, chunk_size);
  while (remaining)
  {
    VALGRIND_MAKE_MEM_UNDEFINED( dataBuffer, bufferSize );
    long count = chunk_size < remaining ? chunk_size : remaining;
    remaining -= count;
    Range::const_iterator end = iter;
    end += count;
    
    for (int d = 0; d < dim; d++)
    {
      if (d < mesh_dim)
      {
        rval = writeUtil->get_node_coords( d, iter, end, count, buffer );
        CHK_MB_ERR_1(rval, node_table, status);
      }
      else
      {
        memset( buffer, 0, count * sizeof(double) );
      }
    
      dbgOut.printf(3,"  writing %c node chunk %ld of %ld, %ld values at %ld\n",
             (char)('X'+d), num_writes - remaining_writes + 1, num_writes, count, offset );
      if (d == 0) track.record_io( offset, count );
      mhdf_writeNodeCoordWithOpt( node_table, offset, count, d, buffer, writeProp, &status );
      CHK_MHDF_ERR_1(status, node_table);
    }
    
    iter = end;
    offset += count;
    --remaining_writes;
  }
  
  // Do empty writes if necessary for parallel collective IO
  if (collectiveIO) {
    while (remaining_writes--) {
      assert(writeProp != H5P_DEFAULT);
      for (int d = 0; d < dim; ++d) {
        dbgOut.printf(3,"  writing (empty) %c node chunk %ld of %ld.\n",
               (char)('X'+d), num_writes - remaining_writes, num_writes );
        mhdf_writeNodeCoordWithOpt( node_table, offset, 0, d, 0, writeProp, &status );
        CHK_MHDF_ERR_1(status, node_table);
      }
    }
  }
  
  mhdf_closeData( filePtr, node_table, &status );
  CHK_MHDF_ERR_0(status);
 
  track.all_reduce();
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::write_elems( ExportSet& elems )
{
  mhdf_Status status;
  ErrorCode rval;
  long first_id;
  int nodes_per_elem;
  long table_size;

  debug_barrier();
  dbgOut.printf(2,"Writing %lu elements of type %s%d\n",
    (unsigned long)elems.range.size(),
    CN::EntityTypeName(elems.type), elems.num_nodes );
  dbgOut.print(3,"Writing elements",elems.range);

  hid_t elem_table = mhdf_openConnectivity( filePtr, 
                                            elems.name(), 
                                            &nodes_per_elem,
                                            &table_size,
                                            &first_id,
                                            &status );
  IODebugTrack track( debugTrack, elems.name() && strlen(elems.name()) 
    ? elems.name() : "<ANONYMOUS ELEM SET?>", table_size );
                                            
  CHK_MHDF_ERR_0(status);
  assert ((unsigned long)first_id <= elems.first_id);
  assert ((unsigned long)table_size >= elems.offset + elems.range.size());
  
  
  EntityHandle* buffer = (EntityHandle*)dataBuffer;
  int chunk_size = bufferSize / (elems.num_nodes * sizeof(id_t));
  long offset = elems.offset;
  long remaining = elems.range.size();
  long num_writes = (remaining+chunk_size-1) / chunk_size;
  if (elems.max_num_ents) {
    assert( elems.max_num_ents >= remaining );
    num_writes = (elems.max_num_ents+chunk_size-1) / chunk_size;
  }
  long remaining_writes = num_writes;
  Range::iterator iter = elems.range.begin();
  
  while (remaining)
  {
    VALGRIND_MAKE_MEM_UNDEFINED( dataBuffer, bufferSize );
    long count = chunk_size < remaining ? chunk_size : remaining;
    remaining -= count;
  
    Range::iterator next = iter;
    next += count;
    rval = writeUtil->get_element_connect( iter, next, elems.num_nodes, 
                                         count * elems.num_nodes, buffer );
    CHK_MB_ERR_1(rval, elem_table, status);
    iter = next;
    
    for (long i = 0; i < count*nodes_per_elem; ++i)
      if (0 == (buffer[i] = idMap.find( buffer[i] )))
        return error(MB_FAILURE);
    
    dbgOut.printf(3,"  writing node connectivity %ld of %ld, %ld values at %ld\n",
           num_writes - remaining_writes + 1, num_writes, count, offset );
    track.record_io( offset, count );
    mhdf_writeConnectivityWithOpt( elem_table, offset, count, 
                                   id_type, buffer, writeProp, &status );
    CHK_MHDF_ERR_1(status, elem_table);
    
    offset += count;
    --remaining_writes;
  }
  
  // Do empty writes if necessary for parallel collective IO
  if (collectiveIO) {
    while (remaining_writes--) {
      assert(writeProp != H5P_DEFAULT);
      dbgOut.printf(3,"  writing (empty) connectivity chunk %ld of %ld.\n",
             num_writes - remaining_writes + 1, num_writes );
      mhdf_writeConnectivityWithOpt( elem_table, offset, 0, id_type, 0, writeProp, &status );
      CHK_MHDF_ERR_1(status, elem_table);
    }
  }
  
  mhdf_closeData( filePtr, elem_table, &status );
  CHK_MHDF_ERR_0(status);
 
  track.all_reduce();
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::get_set_info( EntityHandle set,
                                     long& num_entities,
                                     long& num_children,
                                     long& num_parents,
                                     unsigned long& flags )
{
  ErrorCode rval;
  int i;
  unsigned int u;
  
  rval = iFace->get_number_entities_by_handle( set, i, false );
  CHK_MB_ERR_0(rval);
  num_entities = i;

  rval = iFace->num_child_meshsets( set, &i );
  CHK_MB_ERR_0(rval);
  num_children = i;

  rval = iFace->num_parent_meshsets( set, &i );
  CHK_MB_ERR_0(rval);
  num_parents = i;

  rval = iFace->get_meshset_options( set, u );
  CHK_MB_ERR_0(rval);
  flags = u;
  
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::write_parents_children( bool children )
{
  mhdf_Status status;
  ErrorCode rval = MB_SUCCESS;
  long table_size;
  hid_t table;
  Range::const_iterator iter = setSet.range.begin();
  const Range::const_iterator end = setSet.range.end();
  std::vector<id_t> id_list;
  std::vector<EntityHandle> handle_list;
  
  assert(writeSets);
  
  debug_barrier();

  if (children)
    table = mhdf_openSetChildren( filePtr, &table_size, &status );
  else
    table = mhdf_openSetParents( filePtr, &table_size, &status );
  CHK_MHDF_ERR_0(status);
  IODebugTrack track( debugTrack, children ? "SetChildren" : "SetParents", table_size );
    
  id_t* buffer = reinterpret_cast<id_t*>(dataBuffer);
  const unsigned long buffer_size = bufferSize / sizeof(id_t);
  unsigned long offset = children ? setChildrenOffset : setParentsOffset;
  unsigned long count = 0;
  VALGRIND_MAKE_MEM_UNDEFINED( dataBuffer, bufferSize );
  for (iter = setSet.range.begin(); iter != end; ++iter)
  {
    handle_list.clear();
    if (children)
      rval = iFace->get_child_meshsets( *iter, handle_list, 1 );
    else
      rval = iFace->get_parent_meshsets( *iter, handle_list, 1 );
    CHK_MB_ERR_1(rval, table, status);

    if (handle_list.size() == 0)
      continue;

    id_list.clear();
    rval = vector_to_id_list( handle_list, id_list );
    CHK_MB_ERR_1(rval, table, status);

    if (id_list.size() + count > buffer_size) {
        // buffer is full, flush it
      dbgOut.print(3,"  writing parent/child link chunk\n");
      track.record_io( offset, count );
      mhdf_writeSetParentsChildren( table, offset, count, id_type, buffer, &status );
      CHK_MHDF_ERR_1(status, table);
      offset += count;
      count = 0;
      VALGRIND_MAKE_MEM_UNDEFINED( dataBuffer, bufferSize );

        // If id_list still doesn't it in empty buffer, write it
        // directly rather than trying to buffer it
      if (id_list.size() > buffer_size) {
        dbgOut.print(3,"  writing parent/child link chunk\n");
        track.record_io( offset, id_list.size() );
        mhdf_writeSetParentsChildren( table, offset, id_list.size(), id_type, &id_list[0], &status );
        CHK_MHDF_ERR_1(status, table);
        offset += id_list.size();
        id_list.clear();
      }
    }

    std::copy( id_list.begin(), id_list.end(), buffer + count );
    count += id_list.size();
  }
    
  if (count) {
    dbgOut.print(3,"  writing final parent/child link chunk\n");
    track.record_io( offset, count );
    mhdf_writeSetParentsChildren( table, offset, count, id_type, buffer, &status );
    CHK_MHDF_ERR_1(status, table);
  }

  if (parallelWrite)
    rval = write_shared_set_children( table, &track );
  mhdf_closeData( filePtr, table, &status );
  CHK_MB_ERR_0(rval);

  track.all_reduce();
  return rval;
}


ErrorCode WriteHDF5::write_sets( )
{
  mhdf_Status status;
  Range& sets = setSet.range;
  ErrorCode rval;
  long first_id, meta_size, table_size, content_size, parent_size, child_size;
  hid_t set_table = 0, content_table = 0;
  
  /* If no sets, just return success */
  if (!writeSets)
    return MB_SUCCESS;
  
  /* Write set description table and set contents table */
  
  debug_barrier();
  dbgOut.printf(2,"Writing %lu non-shared sets\n", (unsigned long)setSet.range.size() );
  dbgOut.print(3,"Non-shared sets", setSet.range );
  
  /* Create the table */
  set_table = mhdf_openSetMeta( filePtr, &meta_size, &first_id, &status );
  CHK_MHDF_ERR_0(status);
  IODebugTrack track_meta( debugTrack, "SetMeta", meta_size );
  
  long* buffer = reinterpret_cast<long*>(dataBuffer);
  int chunk_size = bufferSize / (4*sizeof(long));
  long remaining = sets.size();

  id_t* content_buffer = 0;
  unsigned long content_chunk_size = 0;
  unsigned long data_count = 0;
  long content_buffer_offset = setContentsOffset;
  if (writeSetContents)
  {
    content_table = mhdf_openSetData( filePtr, &table_size, &status );
    CHK_MHDF_ERR_1(status, set_table);

    long avg_set_size = (table_size + meta_size - 1) / meta_size;
    if (!avg_set_size)
      avg_set_size = 1;
    chunk_size = bufferSize / (4*sizeof(long) + avg_set_size*sizeof(id_t));
    if (!chunk_size)
      ++chunk_size;
    content_chunk_size = (bufferSize - 4*sizeof(long)*chunk_size)/sizeof(id_t);
    assert(content_chunk_size>0);
    content_buffer = reinterpret_cast<id_t*>(buffer+4*chunk_size);
    VALGRIND_MAKE_MEM_UNDEFINED( content_buffer, content_chunk_size*sizeof(content_buffer[0]) );
  }
  IODebugTrack track_contents( debugTrack, "SetContents", writeSetContents ? table_size : 0  );
    
  Range set_contents;
  Range::const_iterator iter = sets.begin();
  long set_offset = setSet.offset;
  long content_offset = setContentsOffset;
  long child_offset = setChildrenOffset;
  long parent_offset = setParentsOffset;
  unsigned long flags;
  std::vector<id_t> id_list;
  std::vector<EntityHandle> handle_list;
  while (remaining) {
    long* set_data = buffer;
    long count = remaining < chunk_size ? remaining : chunk_size;
    remaining -= count;
      // tell valgrind that buffer portion used for set descriptions
      // is uninitialized (as it is garbage data from the last iteration)
    VALGRIND_MAKE_MEM_UNDEFINED( buffer, 4*sizeof(buffer[0])*chunk_size );

    for (long i = 0; i < count; ++i, ++iter, set_data += 4) {
    
      rval = get_set_info( *iter, content_size, child_size, parent_size, flags );
      CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);

      id_list.clear();
      if (flags & MESHSET_SET) {
        set_contents.clear();

        rval = iFace->get_entities_by_handle( *iter, set_contents, false );
        CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);

        bool blocked_list;
        rval = range_to_blocked_list( set_contents, id_list, blocked_list );
        CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);

        assert (id_list.size() <= (unsigned long)content_size);
        if (blocked_list) {
          assert (id_list.size() % 2 == 0);
          flags |= mhdf_SET_RANGE_BIT;
        }
      }
      else
      {
        handle_list.clear();

        rval = iFace->get_entities_by_handle( *iter, handle_list, false );
        CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);

        rval = vector_to_id_list( handle_list, id_list );
        CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);
      }

      child_offset += child_size;
      parent_offset += parent_size;
      set_data[0] = content_offset + id_list.size() - 1;
      set_data[1] = child_offset - 1;
      set_data[2] = parent_offset - 1;
      set_data[3] = flags;
    
      if (id_list.size())
      {
        if (data_count + id_list.size() > content_chunk_size) {
          dbgOut.print(3,"  writing set content chunk\n");
            // If there isn't enough space remaining in the buffer,
            // flush the buffer.
          track_contents.record_io( content_buffer_offset, data_count );
          mhdf_writeSetData( content_table, 
                             content_buffer_offset,
                             data_count,
                             id_type,
                             content_buffer,
                             &status );
          CHK_MHDF_ERR_2C(status, set_table, writeSetContents, content_table );
          content_buffer_offset += data_count;
          data_count = 0;
          VALGRIND_MAKE_MEM_UNDEFINED( content_buffer, content_chunk_size*sizeof(content_buffer[0]) );
        
            // If there still isn't enough space in the buffer because
            // the size of id_list is bigger than the entire buffer,
            // write id_list directly.
          if (id_list.size() > content_chunk_size) {
            dbgOut.print(3,"  writing set content chunk\n");
            track_contents.record_io( content_buffer_offset, id_list.size() );
            mhdf_writeSetData( content_table, 
                               content_buffer_offset,
                               id_list.size(),
                               id_type,
                               &id_list[0],
                               &status );
            CHK_MHDF_ERR_2C(status, set_table, writeSetContents, content_table );
            content_buffer_offset += id_list.size();
            content_offset += id_list.size();
            id_list.clear();
          }
        }
        
        std::copy( id_list.begin(), id_list.end(), content_buffer+data_count );
        data_count += id_list.size();
        content_offset += id_list.size();
      }
    }

    dbgOut.print(3,"  writing set description chunk.\n");
    track_meta.record_io( set_offset, count );
    mhdf_writeSetMeta( set_table, set_offset, count, H5T_NATIVE_LONG, 
                       buffer, &status );
    CHK_MHDF_ERR_2C(status, set_table, writeSetContents, content_table );
    set_offset += count;
  }
  
  if (data_count) {
    dbgOut.print(3,"  writing final set content chunk\n");
    track_contents.record_io( content_buffer_offset, data_count );
    mhdf_writeSetData( content_table, 
                       content_buffer_offset,
                       data_count,
                       id_type,
                       content_buffer,
                       &status );
    CHK_MHDF_ERR_2C(status, set_table, writeSetContents, content_table );
  }    
  
  if (parallelWrite) {
    rval = write_shared_set_descriptions( set_table, &track_meta );
    CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);
  }
  mhdf_closeData( filePtr, set_table, &status );
  
  rval = MB_SUCCESS;
  if (writeSetContents && parallelWrite) 
    rval = write_shared_set_contents( content_table, &track_contents );
  if (writeSetContents)
    mhdf_closeData( filePtr, content_table, &status );
  CHK_MB_ERR_0( rval );
  
  track_meta.all_reduce();
  track_contents.all_reduce();
  
    /* Write set children */
  if (writeSetChildren)
  {
    rval = write_parents_children( true );
    CHK_MB_ERR_0(rval);
  }
  
    /* Write set parents */
  if (writeSetParents)
  {
    rval = write_parents_children( false );
    CHK_MB_ERR_0(rval);
  }

  return MB_SUCCESS;
}

/*
ErrorCode WriteHDF5::range_to_blocked_list( const Range& input_range,
                                              std::vector<id_t>& output_id_list )
{
  Range::const_iterator r_iter;
  Range::const_iterator const r_end = input_range.end();
  std::vector<id_t>::iterator i_iter, w_iter;
  ErrorCode rval;
  
    // Get file IDs from handles
  output_id_list.resize( input_range.size() );
  w_iter = output_id_list.begin();
  for (r_iter = input_range.begin(); r_iter != r_end; ++r_iter, ++w_iter) 
    *w_iter = idMap.find( *r_iter );
  std::sort( output_id_list.begin(), output_id_list.end() );
  
    // Count the number of ranges in the id list
  unsigned long count = 0;
  bool need_to_copy = false;
  std::vector<id_t>::iterator const i_end = output_id_list.end();
  i_iter = output_id_list.begin();
  while (i_iter != i_end)
  {
    ++count;
    id_t prev = *i_iter;
    for (++i_iter; (i_iter != i_end) && (++prev == *i_iter); ++i_iter);
    if (i_iter - output_id_list.begin() < (long)(2*count))
      need_to_copy = true;
  }
  
    // If the range format is larger than half the size of the
    // the simple list format, just keep the list format
  if (4*count >= output_id_list.size())
    return MB_SUCCESS;
  
    // Convert to ranged format
  std::vector<id_t>* range_list = &output_id_list;
  if (need_to_copy)
    range_list = new std::vector<id_t>( 2*count );

  w_iter = range_list->begin();
  i_iter = output_id_list.begin();
  while (i_iter != i_end)
  {
    unsigned long range_size = 1;
    id_t prev = *w_iter = *i_iter;
    w_iter++;
    for (++i_iter; (i_iter != i_end) && (++prev == *i_iter); ++i_iter)
      ++range_size;
    *w_iter = range_size;
    ++w_iter;
  }

  if (need_to_copy)
  {
    std::swap( *range_list, output_id_list );
    delete range_list;
  }
  else
  {
    assert( w_iter - output_id_list.begin() == (long)(2*count) );
    output_id_list.resize( 2*count );
  }
  
  return MB_SUCCESS;
}
*/
ErrorCode WriteHDF5::range_to_blocked_list( const Range& input_range,
                                              std::vector<id_t>& output_id_list, 
                                              bool& ranged_list )
{
  output_id_list.clear();
  ranged_list = false;
  if (input_range.empty()) {
    return MB_SUCCESS;
  }

    // first try ranged format, but give up if we reach the 
    // non-range format size.
  RangeMap<EntityHandle,id_t>::iterator ri = idMap.begin();
  Range::const_pair_iterator pi;
    // if we end up with more than this many range blocks, then
    // we're better off just writing the set as a simple list
  size_t pairs_remaining = input_range.size() / 2; 
  for (pi = input_range.const_pair_begin(); pi != input_range.const_pair_end(); ++pi) {
    EntityHandle h = pi->first;
    while (h <= pi->second) {
      ri = idMap.lower_bound( ri, idMap.end(), h );
      if (ri == idMap.end() || ri->begin > h) {
        ++h;
        continue;
      }

      id_t n = pi->second - pi->first + 1;
      if (n > ri->count)
        n = ri->count;
  
        // if we ran out of space, (or set is empty) just do list format
      if (!pairs_remaining) {
        output_id_list.resize( input_range.size() );
        range_to_id_list( input_range, &output_id_list[0] );
        output_id_list.erase( std::remove( output_id_list.begin(), 
                                           output_id_list.end(), 
                                           0 ), 
                              output_id_list.end() );
        return MB_SUCCESS;
      }

      --pairs_remaining;
      id_t id = ri->value + (h - ri->begin);
      output_id_list.push_back(id);
      output_id_list.push_back(n);
      h += n;
    }
  }
  
    // if we aren't writing anything (no entities in Range are
    // being written to to file), clean up and return
  if (output_id_list.empty())
    return MB_SUCCESS;
  
    // otherwise check if we can compact the list further
  ranged_list = true;
  size_t r, w = 2;
  const size_t e = output_id_list.size() - 1;
  for (r = 2; r < e; r += 2) {
    if (output_id_list[w-2] + output_id_list[w-1] == output_id_list[r])
      output_id_list[w-1] += output_id_list[r+1];
    else {
      if (w != r) {
        output_id_list[w] = output_id_list[r];
        output_id_list[w+1] = output_id_list[r+1];
      }
      w += 2;
    }
  }
  output_id_list.resize( w );
  assert(output_id_list.size() % 2 == 0);
  
  return MB_SUCCESS;
}
  

ErrorCode WriteHDF5::range_to_id_list( const Range& range,
                                         id_t* array )
{
  VALGRIND_MAKE_MEM_UNDEFINED( array, sizeof(id_t)*range.size() );
  ErrorCode rval = MB_SUCCESS;
  RangeMap<EntityHandle,id_t>::iterator ri = idMap.begin();
  Range::const_pair_iterator pi;
  id_t* i = array;
  for (pi = range.const_pair_begin(); pi != range.const_pair_end(); ++pi) {
    EntityHandle h = pi->first;
    while (h <= pi->second) {
      ri = idMap.lower_bound( ri, idMap.end(), h );
      if (ri == idMap.end() || ri->begin > h) {
        rval = MB_ENTITY_NOT_FOUND;
        *i = 0; 
        ++i;
        ++h;
        continue;
      }

      id_t n = pi->second - h + 1;
      if (n > ri->count)
        n = ri->count;

      id_t id = ri->value + (h - ri->begin);
      for (id_t j = 0; j < n; ++i, ++j)
        *i = id + j;
      h += n;
    }
  }
  assert( i == array + range.size() );
  return rval;
}
 
ErrorCode WriteHDF5::vector_to_id_list( 
                                 const std::vector<EntityHandle>& input,
                                 std::vector<id_t>& output,
                                 bool remove_zeros )
{
  std::vector<EntityHandle>::const_iterator i_iter = input.begin();
  const std::vector<EntityHandle>::const_iterator i_end = input.end();
  output.resize(input.size());
  VALGRIND_MAKE_VEC_UNDEFINED( output );
  std::vector<id_t>::iterator o_iter = output.begin();
  for (; i_iter != i_end; ++i_iter) {
    id_t id = idMap.find( *i_iter );
    if (!remove_zeros || id != 0) {
      *o_iter = id;
      ++o_iter;
    }
  }
  output.erase( o_iter, output.end() );
    
  return MB_SUCCESS;
}



inline ErrorCode WriteHDF5::get_adjacencies( EntityHandle entity,
                                        std::vector<id_t>& adj )
{
  const EntityHandle* adj_array;
  int num_adj;
  ErrorCode rval = writeUtil->get_adjacencies( entity, adj_array, num_adj );
  if (MB_SUCCESS != rval)
    return error(rval);
  
  size_t j = 0;
  adj.resize( num_adj );
  for (int i = 0; i < num_adj; ++i) 
    if (id_t id = idMap.find( adj_array[i] ))
      adj[j++] = id;
  adj.resize( j );
  return MB_SUCCESS;
}


ErrorCode WriteHDF5::write_adjacencies( const ExportSet& elements )
{
  ErrorCode rval;
  mhdf_Status status;
  Range::const_iterator iter;
  const Range::const_iterator end = elements.range.end();
  std::vector<id_t> adj_list;
  
  debug_barrier();
  
  /* Count Adjacencies */
  long count = 0;
  //for (iter = elements.range.begin(); iter != end; ++iter)
  //{
  //  adj_list.clear();
  //  rval = get_adjacencies( *iter, adj_list);
  //  CHK_MB_ERR_0(rval);
  //
  //  if (adj_list.size() > 0)
  //    count += adj_list.size() + 2;
  //}
  
  //if (count == 0)
  //  return MB_SUCCESS;

  long offset = elements.adj_offset;
  if (offset < 0)
    return MB_SUCCESS;
  
  /* Create data list */
  hid_t table = mhdf_openAdjacency( filePtr, elements.name(), &count, &status );
  CHK_MHDF_ERR_0(status);
  IODebugTrack track( debugTrack, "Adjacencies", count );
  
  /* Write data */
  id_t* buffer = (id_t*)dataBuffer;
  long chunk_size = bufferSize / sizeof(id_t); 
  long num_writes = (elements.max_num_adjs + chunk_size - 1)/chunk_size;
  VALGRIND_MAKE_MEM_UNDEFINED( dataBuffer, bufferSize );
  count = 0;
  for (iter = elements.range.begin(); iter != end; ++iter)
  {
    adj_list.clear();
    rval = get_adjacencies( *iter, adj_list );
    CHK_MB_ERR_1(rval, table, status);
    if (adj_list.size() == 0)
      continue;
    
      // If buffer is full, flush it
    if (count + adj_list.size() + 2 > (unsigned long)chunk_size)
    {
      dbgOut.print(3,"  writing adjacency chunk.\n");
      track.record_io( offset, count );
      mhdf_writeAdjacencyWithOpt( table, offset, count, id_type, buffer, writeProp, &status );
      CHK_MHDF_ERR_1(status, table);
      VALGRIND_MAKE_MEM_UNDEFINED( dataBuffer, bufferSize );
      
      offset += count;
      count = 0;
    }
    
    buffer[count++] = idMap.find( *iter );
    buffer[count++] = adj_list.size();
    
    assert (adj_list.size()+2 < (unsigned long)chunk_size);
    memcpy( buffer + count, &adj_list[0], adj_list.size() * sizeof(id_t) );
    count += adj_list.size();
  }
  
  if (count)
  {
    dbgOut.print(2,"  writing final adjacency chunk.\n");
    mhdf_writeAdjacencyWithOpt( table, offset, count, id_type, buffer, writeProp, &status );
    CHK_MHDF_ERR_1(status, table);

    offset += count;
    count = 0;
    --num_writes;
  }

  // Do empty writes if necessary for parallel collective IO
  if (collectiveIO) {
    while (num_writes > 0) {
      --num_writes;
      assert(writeProp != H5P_DEFAULT);
      dbgOut.print(2,"  writing empty adjacency chunk.\n");
      mhdf_writeAdjacencyWithOpt( table, offset, 0, id_type, 0, writeProp, &status );
      CHK_MHDF_ERR_1(status, table );
    }
  }
  
  mhdf_closeData( filePtr, table, &status );
  CHK_MHDF_ERR_0(status);
  
  track.all_reduce();
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::write_tag( const SparseTag& tag_data )
{
  std::string name;
  ErrorCode rval = iFace->tag_get_name( tag_data.tag_id, name );
  if (MB_SUCCESS != rval)
    return error(rval);

  debug_barrier();
  dbgOut.tprintf( 1, "Writing tag: \"%s\"\n", name.c_str() );
 
  int moab_size, elem_size, array_len;
  DataType moab_type;
  mhdf_TagDataType mhdf_type;
  hid_t hdf5_type;
  rval = get_tag_size( tag_data.tag_id, moab_type, moab_size, elem_size,
                       array_len, mhdf_type, hdf5_type );
  if (MB_SUCCESS != rval)
    return error(rval);

  if (array_len == MB_VARIABLE_LENGTH && tag_data.write) {
    dbgOut.printf( 2, "Writing sparse data for var-len tag: \"%s\"\n", name.c_str() );
    rval = write_var_len_tag( tag_data, name, moab_type, hdf5_type, elem_size );
  }
  else {
    int data_len = elem_size;
    if (moab_type != MB_TYPE_BIT)
      data_len *= array_len;
    if (tag_data.write) {
      dbgOut.printf( 2, "Writing sparse data for tag: \"%s\"\n", name.c_str() );
      rval = write_sparse_tag( tag_data, name, moab_type, hdf5_type, data_len );
    }
    for (size_t i = 0; MB_SUCCESS == rval && i < tag_data.denseList.size(); ++i) {
      const ExportSet* set = find( tag_data.denseList[i] );
      assert(0 != set);
      debug_barrier();
      dbgOut.printf( 2, "Writing dense data for tag: \"%s\" on group \"%s\"\n", name.c_str(), set->name() );
      rval = write_dense_tag( tag_data, *set, name, moab_type, hdf5_type, data_len );
    }
  }
 
  H5Tclose( hdf5_type );
  return MB_SUCCESS == rval ? MB_SUCCESS : error(rval);
}

ErrorCode WriteHDF5::write_sparse_ids( const SparseTag& tag_data,
                                       const Range& range,
                                       hid_t id_table,
                                       size_t table_size,
                                       const char* name )
{
  ErrorCode rval;
  mhdf_Status status;

  std::string tname(name ? name : "<UNKNOWN TAG?>");
  tname += " - Ids";
  IODebugTrack track( debugTrack, tname, table_size );

    // Set up data buffer for writing IDs
  size_t chunk_size = bufferSize / sizeof(id_t);
  id_t* id_buffer = (id_t*)dataBuffer;
  
    // Write IDs of tagged entities.
  long remaining = range.size();
  long offset = tag_data.offset;
  long num_writes = (remaining + chunk_size - 1)/chunk_size;
  if (tag_data.max_num_ents) {
    assert(tag_data.max_num_ents >= (unsigned long)remaining);
    num_writes = (tag_data.max_num_ents + chunk_size - 1)/chunk_size;
  }
  Range::const_iterator iter = range.begin();
  while (remaining)
  {
    VALGRIND_MAKE_MEM_UNDEFINED( dataBuffer, bufferSize );

      // write "chunk_size" blocks of data
    long count = (unsigned long)remaining > chunk_size ? chunk_size : remaining;
    remaining -= count;
    Range::const_iterator stop = iter;
    stop += count;
    Range tmp;;
    tmp.merge( iter, stop );
    iter = stop;
    assert(tmp.size() == (unsigned)count);
    
    rval = range_to_id_list( tmp, id_buffer );
    CHK_MB_ERR_0( rval );
    
      // write the data
    dbgOut.print(3,"  writing sparse tag entity chunk.\n");
    track.record_io( offset, count );
    mhdf_writeSparseTagEntitiesWithOpt( id_table, offset, count, id_type, 
                                        id_buffer, writeProp, &status );
    CHK_MHDF_ERR_0( status );
   
    offset += count;
    --num_writes;
  } // while (remaining)

  // Do empty writes if necessary for parallel collective IO
  if (collectiveIO) {
    while (num_writes--) {
      assert(writeProp != H5P_DEFAULT);
      dbgOut.print(3,"  writing empty sparse tag entity chunk.\n");
      mhdf_writeSparseTagEntitiesWithOpt( id_table, offset, 0, id_type, 
                                          0, writeProp, &status );
      CHK_MHDF_ERR_0( status );
    }
  }
  
  track.all_reduce();
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::write_sparse_tag( const SparseTag& tag_data,
                                       const std::string& name,
                                       DataType mb_data_type,
                                       hid_t value_type,
                                       int   value_type_size )
{
  ErrorCode rval;
  mhdf_Status status;
  hid_t tables[3];
  long table_size, data_size;
  
    // get entities for which to write tag values
  Range range;
  rval = get_sparse_tagged_entities( tag_data, range );
  
    //open tables to write info
  mhdf_openSparseTagData( filePtr,
                          name.c_str(),
                          &table_size,
                          &data_size,
                          tables,
                          &status);
  CHK_MHDF_ERR_0(status);
  assert( range.size() + tag_data.offset <= (unsigned long)table_size );
    // fixed-length tag
  assert( table_size == data_size );

    // Write IDs for tagged entities
  rval = write_sparse_ids( tag_data, range, tables[0], table_size, name.c_str() );
  CHK_MB_ERR_2( rval, tables, status );
  mhdf_closeData( filePtr, tables[0], &status );
  CHK_MHDF_ERR_1(status, tables[1]);
  
    // Set up data buffer for writing tag values
  IODebugTrack track( debugTrack, name + " Data", data_size );
  rval = write_tag_values( tag_data.tag_id,
                           tables[1], 
                           tag_data.offset,
                           range,
                           mb_data_type,
                           value_type,
                           value_type_size,
                           tag_data.max_num_ents,
                           track );
  
  mhdf_closeData( filePtr, tables[1], &status );
  CHK_MB_ERR_0(rval);
  CHK_MHDF_ERR_0(status);
  
  track.all_reduce();
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::write_var_len_tag( const SparseTag& tag_data,
                                        const std::string& name,
                                        DataType mb_data_type,
                                        hid_t hdf_type,
                                        int type_size )
{
  ErrorCode rval;
  mhdf_Status status;
  hid_t tables[3];
  long table_size;
  long data_table_size;
  
    // get entities for which to write tag values
  Range range;
  rval = get_sparse_tagged_entities( tag_data, range );
  
    //open tables to write info
  mhdf_openSparseTagData( filePtr,
                          name.c_str(),
                          &table_size,
                          &data_table_size,
                          tables,
                          &status);
  CHK_MHDF_ERR_0(status);
  assert( range.size() + tag_data.offset <= (unsigned long)table_size );

    // Write IDs for tagged entities
  rval = write_sparse_ids( tag_data, range, tables[0], table_size, name.c_str() );
  CHK_MB_ERR_2( rval, tables, status );
  mhdf_closeData( filePtr, tables[0], &status );
  CHK_MHDF_ERR_2(status, tables + 1);


    // Split buffer into four chunks ordered such that there are no 
    // alignment issues:
    //  1) tag data pointer buffer (data from MOAB)
    //  2) tag offset buffer (to be written)
    //  3) tag size buffer (data from MOAB)
    //  4) concatenated tag data buffer (to be written)
  const size_t quarter_buffer_size = bufferSize / 4;
  const size_t num_entities = quarter_buffer_size / sizeof(void*);
  assert( num_entities > 0 );
  const void** const pointer_buffer = reinterpret_cast<const void**>(dataBuffer);
  long* const offset_buffer = reinterpret_cast<long*>(pointer_buffer + num_entities);
  int* const size_buffer = reinterpret_cast<int*>(offset_buffer + num_entities);
  char* const data_buffer = reinterpret_cast<char*>(size_buffer + num_entities);
  assert( data_buffer < bufferSize + dataBuffer );
  const size_t data_buffer_size = dataBuffer + bufferSize - data_buffer;
  VALGRIND_MAKE_MEM_UNDEFINED( data_buffer, data_buffer_size );
  
  IODebugTrack track_dat( debugTrack, name + " Data", data_table_size );
  IODebugTrack track_idx( debugTrack, name + " Indices", table_size );
  
    // offsets into tables
  long offset_offset = tag_data.offset;      // offset at which to write indices
  long data_offset = tag_data.varDataOffset; // offset at which to write data buffer
  long offset = tag_data.varDataOffset;      // used to calculate indices
  
    // iterate in blocks of num_entities entities
  size_t bytes = 0; // occupied size of data buffer
  size_t remaining = range.size();
  Range::const_iterator i = range.begin();
  while (remaining) {
    VALGRIND_MAKE_MEM_UNDEFINED( pointer_buffer, num_entities*sizeof(pointer_buffer[0]) );
    VALGRIND_MAKE_MEM_UNDEFINED( offset_buffer,  num_entities*sizeof(offset_buffer[0]) );
    VALGRIND_MAKE_MEM_UNDEFINED( size_buffer   , num_entities*sizeof(size_buffer[0]) );
  
    const size_t count = remaining < num_entities ? remaining : num_entities;
    remaining -= count;
    
      // get subset of entity handles
    Range::const_iterator e = i; e += count;
    Range subrange; subrange.merge( i, e );
    i = e;
  
      // get pointers and sizes for entities from MOAB
    rval = iFace->tag_get_data( tag_data.tag_id, subrange, pointer_buffer, size_buffer );
    CHK_MB_ERR_2(rval, tables + 1, status);
    
      // calculate end indices from sizes, and process tag data
    for (size_t j = 0; j < count; ++j) {
      const size_t size = size_buffer[j];
      offset += size / type_size;
      offset_buffer[j] = offset - 1;
      
        // if space in data buffer, add current tag value and continue
      assert(size_buffer[j] >= 0);
      const void* ptr = pointer_buffer[j];
      
        // flush buffer if need more room
      if (bytes + size > data_buffer_size) {
          // write out tag data buffer
        if (bytes) { // bytes might be zero if tag value is larger than buffer
          dbgOut.print(2,"  writing var-length sparse tag value chunk.\n");
          track_dat.record_io( data_offset, bytes/type_size );
          assert(hdf_type > 0);
          mhdf_writeTagValues( tables[1], data_offset, 
                               bytes / type_size, 
                               hdf_type, data_buffer, 
                               &status );
          CHK_MHDF_ERR_2(status, tables + 1);
          data_offset += bytes / type_size;
          bytes = 0;
          VALGRIND_MAKE_MEM_UNDEFINED( data_buffer, data_buffer_size );
        }
      }
      
        // special case: if tag data is larger than buffer write it w/out buffering
      if (size > data_buffer_size) {
        if (mb_data_type == MB_TYPE_HANDLE) {
          std::vector<EntityHandle> tmp_storage(size/sizeof(EntityHandle));
          VALGRIND_MAKE_VEC_UNDEFINED( tmp_storage );
          convert_handle_tag( reinterpret_cast<const EntityHandle*>(ptr),
                              &tmp_storage[0], tmp_storage.size() );
          ptr = &tmp_storage[0];
        }
        dbgOut.print(2,"  writing var-length sparse tag value chunk.\n");
        track_dat.record_io( data_offset, size / type_size );
        assert(hdf_type > 0);
        mhdf_writeTagValues( tables[1], data_offset, 
                             size / type_size, hdf_type, ptr,
                             &status );
        CHK_MHDF_ERR_2(status, tables + 1);
        data_offset += size / type_size;
      }
        // otherwise copy data into buffer to be written during a later ieration
      else {
        if (mb_data_type == MB_TYPE_HANDLE) 
          convert_handle_tag( reinterpret_cast<const EntityHandle*>(ptr),
                              reinterpret_cast<EntityHandle*>(data_buffer + bytes), 
                              size/sizeof(EntityHandle) );
        else
          memcpy( data_buffer + bytes, pointer_buffer[j], size );
        bytes += size;
      }
    }
    
      // write offsets
    dbgOut.print(2,"  writing var-length sparse tag index chunk.\n");
    track_idx.record_io( offset_offset, count );
    mhdf_writeSparseTagIndices( tables[2], offset_offset, count, 
                                H5T_NATIVE_LONG, offset_buffer, 
                                &status );
    CHK_MHDF_ERR_2(status, tables + 1);
    offset_offset += count;
  }
  assert( (unsigned long)offset_offset == tag_data.offset + range.size() );
  
    // flush data buffer
  if (bytes) {
      // write out tag data buffer
    dbgOut.print(2,"  writing final var-length sparse tag value chunk.\n");
    track_dat.record_io( data_offset, bytes/type_size );
    assert(hdf_type > 0);
    mhdf_writeTagValues( tables[1], data_offset, bytes / type_size,
                         hdf_type, data_buffer, &status );
    CHK_MHDF_ERR_2(status, tables + 1);
    data_offset += bytes / type_size;
  }
  assert( offset == data_offset );
  
  mhdf_closeData( filePtr, tables[1], &status );
  CHK_MHDF_ERR_1(status, tables[2]);
  mhdf_closeData( filePtr, tables[2], &status );
  CHK_MHDF_ERR_0(status);
  
  track_idx.all_reduce();
  track_dat.all_reduce();
  
  return MB_SUCCESS;
}


ErrorCode WriteHDF5::write_dense_tag( const SparseTag& tag_data,
                                      const ExportSet& elem_data,
                                      const std::string& name,
                                      DataType mb_data_type,
                                      hid_t value_type,
                                      int value_type_size )
{
    //open tables to write info
  mhdf_Status status;
  long table_size;
  hid_t table = mhdf_openDenseTagData( filePtr,
                                       name.c_str(),
                                       elem_data.name(),
                                       &table_size,
                                       &status);
  CHK_MHDF_ERR_0(status);
  assert( elem_data.range.size() + elem_data.offset <= (unsigned long)table_size );
 
  IODebugTrack track( debugTrack, name + " " + elem_data.name() + " Data", table_size );
  ErrorCode rval = write_tag_values( tag_data.tag_id, 
                                     table, 
                                     elem_data.offset,
                                     elem_data.range,
                                     mb_data_type,
                                     value_type,
                                     value_type_size,
                                     elem_data.max_num_ents,
                                     track );
  
  mhdf_closeData( filePtr, table, &status );
  CHK_MB_ERR_0(rval);
  CHK_MHDF_ERR_0(status);
  
  return MB_SUCCESS;
}
  
ErrorCode WriteHDF5::write_tag_values( Tag tag_id,
                                       hid_t data_table,
                                       unsigned long offset_in,
                                       const Range& range,
                                       DataType mb_data_type,
                                       hid_t value_type,
                                       int   value_type_size,
                                       unsigned long max_num_ents,
                                       IODebugTrack& track )
{
  mhdf_Status status;
 
    // Set up data buffer for writing tag values
  size_t chunk_size = bufferSize / value_type_size;
  assert( chunk_size > 0 );
  char* tag_buffer = (char*)dataBuffer;
  
    // Write the tag values
  size_t remaining = range.size();
  size_t offset = offset_in;
  Range::const_iterator iter = range.begin();
  long num_writes = (remaining + chunk_size - 1)/chunk_size;
  if (max_num_ents) {
    assert( max_num_ents >= remaining );
    num_writes = (max_num_ents + chunk_size - 1)/chunk_size;
  }
  while (remaining)
  {
    VALGRIND_MAKE_MEM_UNDEFINED( dataBuffer, bufferSize );
 
      // write "chunk_size" blocks of data
    long count = (unsigned long)remaining > chunk_size ? chunk_size : remaining;
    remaining -= count;
    memset( tag_buffer, 0, count * value_type_size );
    Range::const_iterator stop = iter;
    stop += count;
    Range range;
    range.merge( iter, stop );
    iter = stop;
    assert(range.size() == (unsigned)count);
 
    ErrorCode rval = iFace->tag_get_data( tag_id, range, tag_buffer );
    CHK_MB_ERR_0(rval);
    
      // Convert EntityHandles to file ids
    if (mb_data_type == MB_TYPE_HANDLE)
      convert_handle_tag( reinterpret_cast<EntityHandle*>(tag_buffer), 
                          count * value_type_size / sizeof(EntityHandle) );
    
      // write the data
    dbgOut.print(2,"  writing tag value chunk.\n");
    track.record_io( offset, count );
    assert(value_type > 0);
    mhdf_writeTagValuesWithOpt( data_table, offset, count,
                                value_type, tag_buffer, writeProp, &status );
    CHK_MHDF_ERR_0(status);
   
    offset += count;
    --num_writes;
  } // while (remaining)

  // Do empty writes if necessary for parallel collective IO
  if (collectiveIO) {
    while (num_writes--) {
      assert(writeProp != H5P_DEFAULT);
      dbgOut.print(2,"  writing empty tag value chunk.\n");
      assert(value_type > 0);
      mhdf_writeTagValuesWithOpt( data_table, offset, 0,
                                  value_type, 0, writeProp, &status );
      CHK_MHDF_ERR_0( status );
    }
  }
  
  track.all_reduce();
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::write_qa( const std::vector<std::string>& list )
{
  const char* app = "MOAB";
  const char* vers = MB_VERSION;
  char date_str[64];
  char time_str[64];
  
  std::vector<const char*> strs(list.size() ? list.size() : 4);
  if (list.size() == 0)
  {
    time_t t = time(NULL);
    tm* lt = localtime( &t );
    strftime( date_str, sizeof(date_str), "%D", lt );
    strftime( time_str, sizeof(time_str), "%T", lt );
    
    strs[0] = app;
    strs[1] = vers;
    strs[2] = date_str;
    strs[3] = time_str;
  }
  else
  {
    for (unsigned int i = 0; i < list.size(); ++i)
      strs[i] = list[i].c_str();
  }
  
  mhdf_Status status;
  dbgOut.print(2,"  writing QA history.\n");
  mhdf_writeHistory( filePtr, &strs[0], strs.size(), &status );
  CHK_MHDF_ERR_0(status);
  
  return MB_SUCCESS;
}

/*
ErrorCode WriteHDF5::register_known_tag_types( Interface* iface )
{
  hid_t int4, double16;
  hsize_t dim[1];
  int error = 0;
  ErrorCode rval;
  
  dim[0] = 4;
  int4 = H5Tarray_create( H5T_NATIVE_INT, 1, dim, NULL );
  
  dim[0] = 16;
  double16 = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, dim, NULL );
  
  if (int4 < 0 || double16 < 0)
    error = 1;
  
  struct { const char* name; hid_t type; } list[] = {
    { GLOBAL_ID_TAG_NAME, H5T_NATIVE_INT } ,
    { MATERIAL_SET_TAG_NAME, H5T_NATIVE_INT },
    { DIRICHLET_SET_TAG_NAME, H5T_NATIVE_INT },
    { NEUMANN_SET_TAG_NAME, H5T_NATIVE_INT },
    { HAS_MID_NODES_TAG_NAME, int4 },
    { GEOM_DIMENSION_TAG_NAME, H5T_NATIVE_INT },
    { MESH_TRANSFORM_TAG_NAME, double16 },
    { 0, 0 } };
  
  for (int i = 0; list[i].name; ++i)
  {
    if (list[i].type < 1)
      { ++error; continue; }
    
    Tag handle;
    
    std::string name("__hdf5_tag_type_");
    name += list[i].name;
    
    rval = iface->tag_get_handle( name.c_str(), handle );
    if (MB_TAG_NOT_FOUND == rval)
    {
      rval = iface->tag_create( name.c_str(), sizeof(hid_t), MB_TAG_SPARSE, handle, NULL );
      if (MB_SUCCESS != rval)
        { ++error; continue; }
      
      hid_t copy_id = H5Tcopy( list[i].type );
      rval = iface->tag_set_data( handle, 0, 0, &copy_id );
      if (MB_SUCCESS != rval)
        { ++error; continue; }
    }
  }
  
  H5Tclose( int4 );
  H5Tclose( double16 );
  return error ? MB_FAILURE : MB_SUCCESS;
}
*/

ErrorCode WriteHDF5::gather_tags( const Tag* user_tag_list, int num_tags )
{
  ErrorCode result;
  std::string tagname;
  std::vector<Tag> tag_list;
  std::vector<Tag>::iterator t_itor;
  Range range;
    
    // Get list of Tags to write
  result = writeUtil->get_tag_list( tag_list, user_tag_list, num_tags );
  CHK_MB_ERR_0(result);

    // Get list of tags
  for (t_itor = tag_list.begin(); t_itor != tag_list.end(); ++t_itor)
  {
      // Add tag to export list
    SparseTag tag_data;
    tag_data.tag_id = *t_itor;
    tag_data.offset = 0;
    tag_data.varDataOffset = 0;
    tag_data.max_num_ents = 0;
    tagList.push_back( tag_data );
  }

  return MB_SUCCESS;
}

  // If we support paralle, then this function will have been
  // overridden with an alternate version in WriteHDF5Parallel
  // that supports parallel I/O.  If we're here 
  // then MOAB was not built with support for parallel HDF5 I/O.
ErrorCode WriteHDF5::parallel_create_file( const char* ,
                                    bool ,
                                    const std::vector<std::string>& ,
                                    const Tag*,
                                    int ,
                                    int,
                                    int  )
{
  return error(MB_NOT_IMPLEMENTED);
}

ErrorCode WriteHDF5::serial_create_file( const char* filename,
                                    bool overwrite,
                                    const std::vector<std::string>& qa_records,
                                    const Tag* user_tag_list,
                                    int num_user_tags,
                                    int dimension )
{
  long first_id;
  mhdf_Status status;
  hid_t handle;
  std::list<ExportSet>::iterator ex_itor;
  ErrorCode rval;
  
  const char* type_names[MBMAXTYPE];
  memset( type_names, 0, MBMAXTYPE * sizeof(char*) );
  for (EntityType i = MBEDGE; i < MBENTITYSET; ++i)
    type_names[i] = CN::EntityTypeName( i );
 
    // Create the file
  filePtr = mhdf_createFile( filename, overwrite, type_names, MBMAXTYPE, &status );
  CHK_MHDF_ERR_0(status);
  assert(!!filePtr);

  rval = write_qa( qa_records );
  CHK_MB_ERR_0(rval);
  
    // Create node table
  if (nodeSet.range.size()) {
    nodeSet.total_num_ents = nodeSet.range.size();
    handle = mhdf_createNodeCoords( filePtr, dimension, nodeSet.total_num_ents,
                                    &first_id, &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
    CHK_MHDF_ERR_0(status);
    nodeSet.first_id = (id_t)first_id;
    rval = assign_ids( nodeSet.range, nodeSet.first_id );
    CHK_MB_ERR_0(rval);
  }
  else {
    nodeSet.first_id = std::numeric_limits<id_t>::max();
  } 
  nodeSet.offset = 0;

    // Create element tables
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
  {
    ex_itor->total_num_ents = ex_itor->range.size();
    rval = create_elem_tables( *ex_itor, first_id );
    CHK_MB_ERR_0(rval);
      
    ex_itor->first_id = (id_t)first_id;
    ex_itor->offset = 0;
    rval = assign_ids( ex_itor->range, ex_itor->first_id );
    CHK_MB_ERR_0(rval);
  }

    // create node adjacency table
  id_t num_adjacencies;
#ifdef MB_H5M_WRITE_NODE_ADJACENCIES  
  rval = count_adjacencies( nodeSet.range, num_adjacencies );
  CHK_MB_ERR_0(rval);
  if (num_adjacencies > 0)
  {
    handle = mhdf_createAdjacency( filePtr,
                                   mhdf_node_type_handle(),
                                   num_adjacencies,
                                   &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
    nodeSet.adj_offset = 0;
  }
  else
    nodeSet.adj_offset = -1;
#endif
  
    // create element adjacency tables
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
  {
    rval = count_adjacencies( ex_itor->range, num_adjacencies );
    CHK_MB_ERR_0(rval);
    
    if (num_adjacencies > 0)
    {
      handle = mhdf_createAdjacency( filePtr,
                                     ex_itor->name(),
                                     num_adjacencies,
                                     &status );
      CHK_MHDF_ERR_0(status);
      mhdf_closeData( filePtr, handle, &status );
      ex_itor->adj_offset = 0;
    }
    else
      ex_itor->adj_offset = -1;
  }
  
    // create set tables
  writeSets = !setSet.range.empty();
  if (writeSets)
  {
    long contents_len, children_len, parents_len;
    writeSets = true;
    
    setSet.total_num_ents = setSet.range.size();
    rval = create_set_meta( first_id );
    CHK_MB_ERR_0(rval);

    setSet.first_id = (id_t)first_id;
    rval = assign_ids( setSet.range, setSet.first_id );
    CHK_MB_ERR_0(rval);
    
    rval = count_set_size( setSet.range, contents_len, children_len, parents_len );
    CHK_MB_ERR_0(rval);
    
    rval = create_set_tables( contents_len, children_len, parents_len );
    CHK_MB_ERR_0(rval);
   
    setSet.offset = 0;
    setContentsOffset = 0;
    setChildrenOffset = 0;
    setParentsOffset = 0;
    writeSetContents = !!contents_len;
    writeSetChildren = !!children_len;
    writeSetParents = !!parents_len;
  } // if(!setSet.range.empty())
  
  
  dbgOut.tprint( 1, "Gathering Tags\n" );
  
  rval = gather_tags( user_tag_list, num_user_tags );
  CHK_MB_ERR_0(rval);

    // Create the tags and tag data tables
  std::list<SparseTag>::iterator tag_iter = tagList.begin();
  for ( ; tag_iter != tagList.end(); ++tag_iter)
  {
      // As we haven't yet added any ExportSets for which to write
      // dense tag data to the SparseTag struct pointed to by
      // tag_iter, this call will initially return all tagged entities
      // in the set of entities to be written.
    Range range;
    rval = get_sparse_tagged_entities( *tag_iter, range );
    CHK_MB_ERR_0(rval);
    
    int s;
    bool var_len = (MB_VARIABLE_DATA_LENGTH == iFace->tag_get_size( tag_iter->tag_id, s ));
    
      // Determine which ExportSets we want to write dense
      // data for. We never write dense data for variable-length
      // tag data.
    if (!var_len && writeTagDense) {
      if (!nodeSet.range.empty() && range.contains(nodeSet.range)) {
        range -= nodeSet.range;
        tag_iter->denseList.push_back( nodeSet );
      }

      std::list<ExportSet>::const_iterator ex = exportList.begin();
      for ( ; ex != exportList.end(); ++ex) {
        if (!ex->range.empty() && range.contains( ex->range )) {
          range -= ex->range;
          tag_iter->denseList.push_back( *ex );
        }
      }

      if (!setSet.range.empty() && range.contains(setSet.range)) {
        range -= setSet.range;
        tag_iter->denseList.push_back( setSet );
      }
    }
    
    tag_iter->write = !range.empty();
  
    unsigned long var_len_total = 0;
    if (var_len) {
      rval = get_tag_data_length( *tag_iter, range, var_len_total ); 
      CHK_MB_ERR_0(rval);
    }
  
    rval = create_tag( *tag_iter, range.size(), var_len_total );
    CHK_MB_ERR_0(rval);
  } // for(tags)
  
  return MB_SUCCESS;
}


ErrorCode WriteHDF5::count_adjacencies( const Range& set, id_t& result )
{
  ErrorCode rval;
  std::vector<id_t> adj_list;
  Range::const_iterator iter = set.begin();
  const Range::const_iterator end = set.end();
  result = 0;
  for ( ; iter != end; ++iter )
  {
    adj_list.clear();
    rval = get_adjacencies( *iter, adj_list );
    CHK_MB_ERR_0(rval);
    
    if (adj_list.size() > 0)
      result += 2 + adj_list.size();
  }
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::create_elem_tables( const ExportSet& block,
                                         long& first_id_out )
{
  char name[64];
  mhdf_Status status;
  hid_t handle;
  
  sprintf( name, "%s%d", CN::EntityTypeName(block.type), block.num_nodes );
  mhdf_addElement( filePtr, name, block.type, &status );
  CHK_MHDF_ERR_0(status);
  
  handle = mhdf_createConnectivity( filePtr, 
                                    name,
                                    block.num_nodes,
                                    block.total_num_ents,
                                    &first_id_out,
                                    &status );
  CHK_MHDF_ERR_0(status);
  mhdf_closeData( filePtr, handle, &status );
  CHK_MHDF_ERR_0(status);
  
  return MB_SUCCESS;
}


ErrorCode WriteHDF5::count_set_size( const Range& sets, 
                                     long& contents_length_out,
                                     long& children_length_out,
                                     long& parents_length_out )
{
  ErrorCode rval;
  Range set_contents;
  Range::const_iterator iter = sets.begin();
  const Range::const_iterator end = sets.end();
  long contents_length_set, children_length_set, parents_length_set;
  unsigned long flags;
  std::vector<id_t> set_contents_ids;
  
  contents_length_out = 0;
  children_length_out = 0;
  parents_length_out = 0;
  
  for (; iter != end; ++iter)
  {
    rval = get_set_info( *iter, contents_length_set, children_length_set,
                         parents_length_set, flags );
    CHK_MB_ERR_0(rval);
    
      // check if can and should compress as ranges
    if ((flags&MESHSET_SET) && !(flags&MESHSET_ORDERED) && contents_length_set > 2)
    {
      set_contents.clear();
      rval = iFace->get_entities_by_handle( *iter, set_contents, false );
      CHK_MB_ERR_0(rval);
      
      bool blocked_list;
      rval = range_to_blocked_list( set_contents, set_contents_ids, blocked_list );
      CHK_MB_ERR_0(rval);
      
      if (blocked_list)
      {
        assert (set_contents_ids.size() % 2 == 0);
        contents_length_set = set_contents_ids.size();
      }
    }

    contents_length_out += contents_length_set;
    children_length_out += children_length_set;
    parents_length_out += parents_length_set;
  }
  
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::create_set_meta( long& first_id_out )
{
  hid_t handle;
  mhdf_Status status;
  
  handle = mhdf_createSetMeta( filePtr, setSet.total_num_ents, &first_id_out, &status );
  CHK_MHDF_ERR_0(status);
  mhdf_closeData( filePtr, handle, &status );
  
  return MB_SUCCESS;
}


ErrorCode WriteHDF5::create_set_tables( long num_set_contents,
                                        long num_set_children,
                                        long num_set_parents )
{
  hid_t handle;
  mhdf_Status status;
  
  if (num_set_contents > 0)
  {
    handle = mhdf_createSetData( filePtr, num_set_contents, &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
  }
  
  if (num_set_children > 0)
  {
    handle = mhdf_createSetChildren( filePtr, num_set_children, &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
  }
  
  if (num_set_parents > 0)
  {
    handle = mhdf_createSetParents( filePtr, num_set_parents, &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
  }
  
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::get_tag_size( Tag tag,
                                   DataType& moab_type,
                                   int& num_bytes,
                                   int& type_size,
                                   int& array_length,
                                   mhdf_TagDataType& file_type,
                                   hid_t& hdf_type )
{
  ErrorCode rval;
  Tag type_handle;
  std::string tag_name, tag_type_name;
   
    // We return NULL for hdf_type if it can be determined from
    // the file_type.  The only case where it is non-zero is
    // if the user specified a specific type via a mesh tag.
  hdf_type = (hid_t)0;
  bool close_hdf_type;
  
  rval = iFace->tag_get_data_type( tag, moab_type ); CHK_MB_ERR_0(rval);
  rval = iFace->tag_get_size( tag, num_bytes );     
  if (MB_VARIABLE_DATA_LENGTH == rval)
    num_bytes = MB_VARIABLE_LENGTH;
  else if (MB_SUCCESS != rval)
    return error(rval);

  switch (moab_type)
  {
  case MB_TYPE_INTEGER:
    type_size = sizeof(int);
    file_type = mhdf_INTEGER;
    hdf_type = H5T_NATIVE_INT;
    close_hdf_type = false;
    array_length = num_bytes/type_size;
    break;
  case MB_TYPE_DOUBLE:
    type_size = sizeof(double);
    file_type = mhdf_FLOAT;
    hdf_type = H5T_NATIVE_DOUBLE;
    close_hdf_type = false;
    array_length = num_bytes/type_size;
    break;
  case MB_TYPE_BIT:
    type_size = sizeof(bool);
    file_type = mhdf_BITFIELD;
    if (num_bytes <= 8)
      hdf_type = H5Tcopy( H5T_NATIVE_B8 );
    else if (num_bytes <= 16)
      hdf_type = H5Tcopy( H5T_NATIVE_B16 );
    else if (num_bytes <= 32)
      hdf_type = H5Tcopy( H5T_NATIVE_B32 );
    else if (num_bytes <= 64)
      hdf_type = H5Tcopy( H5T_NATIVE_B64 );
    else
      return error(MB_FAILURE);
    H5Tset_precision( hdf_type, num_bytes );
    close_hdf_type = true;
    array_length = num_bytes;
    num_bytes = (num_bytes+7)/8;
    break;
  case MB_TYPE_HANDLE:
    type_size = sizeof(EntityHandle);
    file_type = mhdf_ENTITY_ID;
    hdf_type = id_type;
    close_hdf_type = false;
    array_length = num_bytes/type_size;
    break;
  case MB_TYPE_OPAQUE:
    file_type = mhdf_OPAQUE;

    rval = iFace->tag_get_name( tag, tag_name ); CHK_MB_ERR_0(rval);
    tag_type_name = "__hdf5_tag_type_";
    tag_type_name += tag_name;
    rval = iFace->tag_get_handle( tag_type_name.c_str(), type_handle );
    if (MB_TAG_NOT_FOUND == rval) {
      if (num_bytes == MB_VARIABLE_LENGTH)
        type_size = 1;
      else
        type_size = num_bytes;
      hdf_type = H5Tcreate( H5T_OPAQUE, type_size );
      close_hdf_type = true;
    }
    else if (MB_SUCCESS == rval) {
      int hsize;
      rval = iFace->tag_get_size( type_handle, hsize );
      if (hsize != sizeof(hid_t))
        return error(MB_FAILURE);
      
      rval = iFace->tag_get_data( type_handle, 0, 0, &hdf_type );
      if (rval != MB_SUCCESS)
        return error(rval);
        
      type_size = H5Tget_size(hdf_type);
      if (type_size != num_bytes)
        return error(MB_FAILURE);
        
      close_hdf_type = false;
    }
    else {
      return error(rval);
    }
    array_length = 1;
  }
  
  if (num_bytes == MB_VARIABLE_LENGTH) {
    array_length = MB_VARIABLE_LENGTH;
    if (!close_hdf_type) {
      hdf_type = H5Tcopy( hdf_type );
      close_hdf_type = true;
    }
  }
  else if (array_length > 1 && moab_type != MB_TYPE_BIT) {
    hsize_t len = array_length;
#if defined(H5Tarray_create_vers) && (H5Tarray_create_vers > 1)
    hid_t temp_id = H5Tarray_create2( hdf_type, 1, &len);
#else
    hid_t temp_id = H5Tarray_create( hdf_type, 1, &len, NULL );
#endif
    if (close_hdf_type)
      H5Tclose( hdf_type );
    hdf_type = temp_id;
  }
  else if (!close_hdf_type) {
    hdf_type = H5Tcopy( hdf_type );
    close_hdf_type = true;
  }
  
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::get_tag_data_length( const SparseTag& tag_info, 
                                          const Range& range,
                                          unsigned long& result )
{
  ErrorCode rval;
  result = 0;
  
    // split buffer into two pieces, one for pointers and one for sizes
  size_t step, remaining;
  step = bufferSize / (sizeof(int) + sizeof(void*));
  const void** ptr_buffer = reinterpret_cast<const void**>(dataBuffer);
  int* size_buffer = reinterpret_cast<int*>(ptr_buffer + step); 
  Range subrange;
  Range::const_iterator iter = range.begin();
  for (remaining = range.size(); remaining >= step; remaining -= step) {
      // get subset of range containing 'count' entities
    Range::const_iterator end = iter; end += step;
    subrange.clear();
    subrange.merge( iter, end );
    iter = end;
      // get tag sizes for entities
    rval = iFace->tag_get_data( tag_info.tag_id, subrange, ptr_buffer, size_buffer );
    if (MB_SUCCESS != rval)
      return error(rval);
      // sum lengths
    for (size_t i = 0; i < step; ++i)
      result += size_buffer[i];
  }
    // process remaining
  subrange.clear();
  subrange.merge( iter, range.end() );
  assert( subrange.size() == remaining );
  rval = iFace->tag_get_data( tag_info.tag_id, subrange, ptr_buffer, size_buffer );
  if (MB_SUCCESS != rval)
    return error(rval);
  for (size_t i= 0; i < remaining; ++i)
    result += size_buffer[i];
    
  DataType type;
  rval = iFace->tag_get_data_type( tag_info.tag_id, type );
  if (MB_SUCCESS != rval)
    return error(rval);
  switch (type) {
    case MB_TYPE_INTEGER: result /= sizeof(int);            break;
    case MB_TYPE_DOUBLE:  result /= sizeof(double);         break;
    case MB_TYPE_HANDLE:  result /= sizeof(EntityHandle); break;
    case MB_TYPE_OPAQUE:                                    break;
      // We fail for MB_TYPE_BIT because MOAB currently does
      // not support variable-length bit tags.
    default:          return error(MB_FAILURE);
  }
    
  return MB_SUCCESS;
}
    
                                     

ErrorCode WriteHDF5::create_tag( const SparseTag& tag_data,
                                 unsigned long num_sparse_entities,
                                 unsigned long data_table_size )
{
  TagType mb_storage;
  DataType mb_type;
  mhdf_TagDataType mhdf_type;
  int tag_size, elem_size, mhdf_size, storage;
  hid_t hdf_type = (hid_t)0;
  hid_t handles[3];
  std::string tag_name;
  ErrorCode rval;
  mhdf_Status status;
  

    // get tag properties
  rval = iFace->tag_get_type( tag_data.tag_id, mb_storage  ); CHK_MB_ERR_0(rval);
  switch (mb_storage) {
    case MB_TAG_DENSE :  storage = mhdf_DENSE_TYPE ; break;
    case MB_TAG_SPARSE:  storage = mhdf_SPARSE_TYPE; break;
    case MB_TAG_BIT:     storage = mhdf_BIT_TYPE;    break;
    case MB_TAG_MESH:    storage = mhdf_MESH_TYPE;   break;
    default: return error(MB_FAILURE);
  }
  rval = iFace->tag_get_name( tag_data.tag_id, tag_name ); CHK_MB_ERR_0(rval);
  rval = get_tag_size( tag_data.tag_id, mb_type, tag_size, elem_size, mhdf_size, mhdf_type, hdf_type );
  CHK_MB_ERR_0(rval);
  
    // get default value
  const void *def_value, *mesh_value;
  int def_val_len, mesh_val_len;
  rval = iFace->tag_get_default_value( tag_data.tag_id, def_value, def_val_len );
  if (MB_ENTITY_NOT_FOUND == rval) {
    def_value = 0;
    def_val_len = 0;
  }
  else if (MB_SUCCESS != rval) {
    H5Tclose( hdf_type );
    return error(rval);
  }
    
    // get mesh value
  rval = iFace->tag_get_data( tag_data.tag_id, 0, 0, &mesh_value, &mesh_val_len );
  if (MB_TAG_NOT_FOUND == rval) {
    mesh_value = 0;
    mesh_val_len = 0;
  }
  else if (MB_SUCCESS != rval) {
    H5Tclose( hdf_type );
    return error(rval);
  }
  
    // for handle-type tags, need to convert from handles to file ids
  if (MB_TYPE_HANDLE == mb_type) {
      // make sure there's room in the buffer for both
    assert( (def_val_len + mesh_val_len) * sizeof(long) < (size_t)bufferSize );

      // convert default value
    if (def_value) {
      memcpy( dataBuffer, def_value, def_val_len );
      if (convert_handle_tag( reinterpret_cast<EntityHandle*>(dataBuffer), 
                              def_val_len / sizeof(EntityHandle) ))
        def_value = dataBuffer;
      else
        def_value = 0;
    }
    
      // convert mesh value
    if (mesh_value) {
      EntityHandle* ptr = reinterpret_cast<EntityHandle*>(dataBuffer + def_val_len);
      memcpy( ptr, mesh_value, mesh_val_len );
      if (convert_handle_tag( ptr, mesh_val_len / sizeof(EntityHandle) ))
        mesh_value = ptr;
      else
        mesh_value = 0;
    }
  }
     
 
  if (MB_VARIABLE_LENGTH != tag_size) {
      // write the tag description to the file
    mhdf_createTag( filePtr,
                    tag_name.c_str(),
                    mhdf_type,
                    mhdf_size,
                    storage,
                    def_value,
                    mesh_value,
                    hdf_type,
                    mb_type == MB_TYPE_HANDLE ? id_type : 0,
                    &status );
    H5Tclose(hdf_type);
    CHK_MHDF_ERR_0(status);


      // create empty table for tag data
    if (num_sparse_entities)
    {
      mhdf_createSparseTagData( filePtr, 
                                tag_name.c_str(), 
                                num_sparse_entities,
                                handles,
                                &status );
      CHK_MHDF_ERR_0(status);
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
    }
    
    for (size_t i = 0; i < tag_data.denseList.size(); ++i) {
      const ExportSet* ex = find( tag_data.denseList[i] );
      assert(0 != ex);
      handles[0] = mhdf_createDenseTagData( filePtr,
                                            tag_name.c_str(),
                                            ex->name(),
                                            ex->total_num_ents,
                                            &status );
      CHK_MHDF_ERR_0(status);
      mhdf_closeData( filePtr, handles[0], &status );
    }
  }
  else {
    mhdf_createVarLenTag( filePtr,
                          tag_name.c_str(),
                          mhdf_type,
                          storage,
                          def_value, def_val_len / elem_size,
                          mesh_value, mesh_val_len / elem_size,
                          hdf_type, mb_type == MB_TYPE_HANDLE ? id_type : 0,
                          &status );
    H5Tclose(hdf_type);
    CHK_MHDF_ERR_0(status);
    
      // create empty table for tag data
    if (num_sparse_entities) {
      mhdf_createVarLenTagData( filePtr, 
                                tag_name.c_str(),
                                num_sparse_entities,
                                data_table_size,
                                handles,
                                &status );
      CHK_MHDF_ERR_0(status);
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      mhdf_closeData( filePtr, handles[2], &status );
    }
  }
    
  return MB_SUCCESS;
}

ErrorCode WriteHDF5::get_num_sparse_tagged_entities( const SparseTag& tag, 
                                                     size_t& count )
{
  Range tmp;
  ErrorCode rval = get_sparse_tagged_entities( tag, tmp );
  count = tmp.size();
  return rval;
}

ErrorCode WriteHDF5::get_sparse_tagged_entities( const SparseTag& tag,
                                                    Range& results )
{
  results.clear();
  if (!tag.have_dense(setSet))
    results.merge( setSet.range );
  std::list<ExportSet>::const_reverse_iterator e;
  for (e = exportList.rbegin(); e != exportList.rend(); ++e) 
    if (!tag.have_dense(*e))
      results.merge( e->range );
  if (!tag.have_dense(nodeSet))
    results.merge( nodeSet.range );
  if (results.empty())
    return MB_SUCCESS;  
    
  return iFace->get_entities_by_type_and_tag( 0, MBMAXTYPE, 
                                              &tag.tag_id, 0, 1, 
                                              results, Interface::INTERSECT );
}

void WriteHDF5::get_write_entities( Range& range )
{
  range.clear();
  range.merge( setSet.range );
  std::list<ExportSet>::const_reverse_iterator e;
  for (e = exportList.rbegin(); e != exportList.rend(); ++e) 
    range.merge( e->range );
  range.merge( nodeSet.range );
}

void WriteHDF5::print_id_map( ) const
{
  print_id_map( std::cout, "" ) ;
}

void WriteHDF5::print_id_map( std::ostream& s, const char* pfx ) const
{
  RangeMap<EntityHandle,id_t>::const_iterator i;
  for (i = idMap.begin(); i != idMap.end(); ++i) {
    const char* n1 = CN::EntityTypeName(TYPE_FROM_HANDLE(i->begin));
    EntityID id = ID_FROM_HANDLE(i->begin);
    if (i->count == 1) {
      s << pfx << n1 << " " << id << " -> " << i->value << std::endl;
    }
    else {
      const char* n2 = CN::EntityTypeName(TYPE_FROM_HANDLE(i->begin + i->count - 1));
      if (n1 == n2) {
        s << pfx << n1 << " " << id << "-" << id + i->count-1
          << " -> " << i->value << "-" << i->value + i->count-1 << std::endl;
      }
      else {
        s << pfx << n1 << " " << id << "-" 
          << n1 << " " << ID_FROM_HANDLE(i->begin + i->count-1)
          << " -> " << i->value << "-" << i->value + i->count-1 << std::endl;
      }
    }
  }
}

} // namespace moab

