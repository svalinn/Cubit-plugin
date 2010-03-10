#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "assert.h"

#include "ReadIDEAS.hpp"
#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBReadUtilIface.hpp"
#include "FileTokenizer.hpp"
#include "MBRange.hpp"

MBReaderIface* ReadIDEAS::factory( MBInterface* iface )
  { return new ReadIDEAS( iface ); }

ReadIDEAS::ReadIDEAS(MBInterface* impl)
    : MBI(impl)
{
  void* ptr = 0;
  impl->query_interface("MBReadUtilIface", &ptr);
  readMeshIface = reinterpret_cast<MBReadUtilIface*>(ptr);
}


MBErrorCode ReadIDEAS::read_tag_values( const char* /* file_name */,
                                        const char* /* tag_name */,
                                        const FileOptions& /* opts */,
                                        std::vector<int>& /* tag_values_out */,
                                        const IDTag* /* subset_list */,
                                        int /* subset_list_length */ )
{
  return MB_NOT_IMPLEMENTED;
}


MBErrorCode ReadIDEAS::load_file(const char* fname, 
                                 const MBEntityHandle* , 
                                 const FileOptions& options,
                                 const MBReaderIface::IDTag* subset_list,
                                 int subset_list_length,
                                 const MBTag* file_id_tag ) {

  if (subset_list && subset_list_length) {
    readMeshIface->report_error( "Reading subset of files not supported for IDEAS." );
    return MB_UNSUPPORTED_OPERATION;
  }

  file.open( fname );
  if (!file.good()) {
    readMeshIface->report_error("Failed to open file: %s", fname);
    return MB_FILE_DOES_NOT_EXIST;
  }

  MBErrorCode rval;

  char line[10000];
  file.getline(line, 10000);
  std::string s = line;
  if (s.find("-1") > s.length()) return MB_FAILURE;

  MBEntityHandle first_vertex = 0;

  while (! file.eof() ) {
    file.getline(line, 10000);
    s = line;

    unsigned int header_id = (unsigned int) strtol(line, NULL, 10);
    switch (header_id) {
      case VERTEX_LIST :
        if (first_vertex) // multiple vertex blocks?
          return MB_FAILURE;
        rval = create_vertices( first_vertex, file_id_tag ); 
      break;
      case MAKE_TETRAHEDRA :
        if (!first_vertex) // need to read vertices first
          return MB_FAILURE;
        rval = create_tetrahedral_elements( first_vertex, file_id_tag );
      break;
      default:
        rval = skip_header(); 
      break;
    }
  }

  file.close();
  return MB_SUCCESS;

}

MBErrorCode ReadIDEAS::skip_header() {

  // Go until finding a pair of -1 lines
  char *ctmp;
  char line[10000];
  std::string s;

  int end_of_block = 0;

  long int il;

  while (file.getline(line, 10000)) {
 
    il = std::strtol(line, &ctmp, 10);
    if (il == -1) {
      s = ctmp;
      if (s.empty()) end_of_block++;
    }
    else end_of_block = 0;

    if (end_of_block >= 2)
      return MB_SUCCESS;

  }

  return MB_FAILURE;
}



MBErrorCode ReadIDEAS::create_vertices(MBEntityHandle& first_vertex,
                                       const MBTag* file_id_tag) {

  // Read two lines: first has some data, second has coordinates
  char line1[10000], line2[10000];
  int il1, il2;
  char *ctmp1, *ctmp2;
  std::string s1, s2;

  MBErrorCode rval;

  int top_of_block = file.tellg();
  unsigned int num_verts = 0;

  for (;;) {

    if (!file.getline(line1, 10000))
      return MB_FAILURE;
    if (!file.getline(line2, 10000))
      return MB_FAILURE;

    // Check if we are at the end of the block
    il1 = std::strtol(line1, &ctmp1, 10);
    il2 = std::strtol(line2, &ctmp2, 10);
    if ((il1 == -1) && (il2 == -1)) {
      s1 = ctmp1;
      s2 = ctmp2;
      if ((s1.empty()) && (s2.empty())) break;     
    }
    num_verts++;
  }

  file.seekg( top_of_block );
  
  std::vector<double*> arrays;
  rval = readMeshIface->get_node_arrays( 3, num_verts, 0, first_vertex, arrays );
  if (MB_SUCCESS != rval)
    return rval;

  MBRange verts;
  verts.insert( first_vertex, first_vertex + num_verts - 1 );
  
  double *x = arrays[0];
  double *y = arrays[1];
  double *z = arrays[2];
  for (unsigned int i = 0; i < num_verts; i++) {

    if (!file.getline(line1, 10000))
      return MB_FAILURE;
    if (!file.getline(line2, 10000))
      return MB_FAILURE;

    // Get the doubles out of the 2nd line
    x[i] = std::strtod(line2, &ctmp2);
    y[i] = std::strtod(ctmp2+1, &ctmp2);
    z[i] = std::strtod(ctmp2+1, NULL);

  }

  if (!file.getline(line1, 10000))
    return MB_FAILURE;
  if (!file.getline(line2, 10000))
    return MB_FAILURE;

  if (file_id_tag) {
    rval = readMeshIface->assign_ids( *file_id_tag, verts, 1 );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  return MB_SUCCESS;
}


MBErrorCode ReadIDEAS::create_tetrahedral_elements(MBEntityHandle vstart,
                                                   const MBTag* file_id_tag) {

  MBEntityHandle connect[4];
  char line1[10000], line2[10000];
  int il1, il2, id = 0;
  char *ctmp1, *ctmp2;
  std::string s1, s2;
  long verts[4];

  MBErrorCode rval;
  MBEntityHandle handle;

  MBTag mat_prop_tag, phys_prop_tag;
  rval = MBI->tag_create( MAT_PROP_TABLE_TAG  , sizeof(int), MB_TAG_DENSE, mat_prop_tag, 0); 
  if (MB_SUCCESS != rval && MB_ALREADY_ALLOCATED != rval) return rval;
  rval = MBI->tag_create( PHYS_PROP_TABLE_TAG , sizeof(int), MB_TAG_DENSE, phys_prop_tag, 0); 
  if (MB_SUCCESS != rval && MB_ALREADY_ALLOCATED != rval) return rval;
 
  for (;;) {

    if (!file.getline(line1, 10000) || !file.getline(line2, 10000))
      return MB_FAILURE;

    // Check if we are at the end of the block
    il1 = std::strtol(line1, &ctmp1, 10);
    il2 = std::strtol(line2, &ctmp2, 10);
    if ((il1 == -1) && (il2 == -1)) {
      s1 = ctmp1;
      s2 = ctmp2;
      if ((s1.empty()) && (s2.empty())) 
        return MB_SUCCESS;     
    }

    // Get property tables out of 1st line
    int phys_table = strtol(line1+21, &ctmp1, 10);
    int mat_table  = strtol(line1+31, &ctmp1, 10);

    // Get the connectivity out of the 2nd line
    if (4 != sscanf( line2, "%ld %ld %ld %ld", verts,verts+1,verts+2,verts+3))
      return MB_FAILURE;
    connect[0] = vstart + verts[0] - 1;
    connect[1] = vstart + verts[1] - 1;
    connect[2] = vstart + verts[2] - 1;
    connect[3] = vstart + verts[3] - 1;

    // Make the element
    rval = MBI->create_element(MBTET, connect, 4, handle);
    assert( MB_SUCCESS == rval );
    
    rval = MBI->tag_set_data(phys_prop_tag,&handle,1,&phys_table);
    assert( MB_SUCCESS == rval);
    
    rval = MBI->tag_set_data(mat_prop_tag,&handle,1,&mat_table);
    assert( MB_SUCCESS == rval);
    
    if (file_id_tag) {
      rval = MBI->tag_set_data( *file_id_tag, &handle, 1, &id );
      ++id;
    }
  }
  
  return MB_FAILURE;
}
