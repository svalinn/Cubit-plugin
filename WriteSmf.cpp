/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */


#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif


#include "WriteSmf.hpp"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>

#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include "MBTagConventions.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBInternals.hpp"
#include "FileOptions.hpp"
#include "MBVersion.h"

const int DEFAULT_PRECISION = 10;
const bool DEFAULT_STRICT = true;

MBWriterIface *WriteSmf::factory( MBInterface* iface )
  { return new WriteSmf( iface ); }

WriteSmf::WriteSmf(MBInterface *impl) 
    : mbImpl(impl), writeTool(0)
{
  assert(impl != NULL);

  void* ptr = 0;
  impl->query_interface( "MBWriteUtilIface", &ptr );
  writeTool = reinterpret_cast<MBWriteUtilIface*>(ptr);
}

WriteSmf::~WriteSmf() 
{
  mbImpl->release_interface("MBWriteUtilIface", writeTool);
}

MBErrorCode WriteSmf::write_file(const char *file_name, 
                                 const bool overwrite,
                                 const FileOptions& opts,
                                 const MBEntityHandle *output_list,
                                 const int num_sets,
                                 const std::vector<std::string>& ,
                                 const MBTag* tag_list,
                                 int num_tags,
                                 int )
{
  MBErrorCode rval;

    // Get precision for node coordinates
  int precision;
  if (MB_SUCCESS != opts.get_int_option( "PRECISION", precision ))
    precision = DEFAULT_PRECISION;
  
   // Honor overwrite flag
  if (!overwrite)
  {
    rval = writeTool->check_doesnt_exist( file_name );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
    // Create file
  std::ofstream file( file_name );
  if (!file)
  {
    writeTool->report_error("Could not open file: %s\n", file_name );
    return MB_FILE_WRITE_ERROR;
  }
  file.precision( precision );
    // Get entities to write
  
  MBRange triangles;
  if (!output_list || !num_sets)
  {
    rval = mbImpl->get_entities_by_type( 0, MBTRI, triangles, false);
    if (MB_SUCCESS != rval) return rval;

    // somehow get all the nodes from this range, order them, uniquify, then use binary search
  }
  else
  {
    // not implemented yet, get out
    // support export only of all triangles from the mesh
    return  MB_NOT_IMPLEMENTED;
  }
  // use an array with all the connectivities in the triangles; it will be converted later to ints
  int numTriangles = triangles.size();
  int array_alloc = 3*numTriangles;       // allocated size of 'array'
  MBEntityHandle* array = new MBEntityHandle[array_alloc]; // ptr to working array of result handles
  // fill up array with node handles; reorder and uniquify 
  if (!array)
     return MB_MEMORY_ALLOCATION_FAILED;
  int fillA = 0;
  for (MBRange::const_iterator e = triangles.begin(); e != triangles.end(); ++e)
    {
      const MBEntityHandle* conn;
      int conn_len;
      rval = mbImpl->get_connectivity( *e, conn, conn_len );
      if (MB_SUCCESS != rval  ) 
	return rval;
      if ( 3!=conn_len) 
        return MB_INVALID_SIZE;

      for (int i = 0; i < conn_len; ++i)
        array[fillA++] = conn[i];
    }
  if (fillA != array_alloc)
	 return MB_INVALID_SIZE;
    
  std::sort( array, array + array_alloc);
  int numNodes = std::unique(array, array + array_alloc ) - array;
  
  file << "#$SMF 1.0\n";
  file << "#$vertices " << numNodes << std::endl;
  file << "#$faces " << numTriangles << std::endl;
  file << "# \n";
  file << "# output from MOAB \n";
  file << "# \n";
  
  // output first the nodes
  // num nodes??
  // write the nodes 
  double coord[3];
  for(int i=0; i<numNodes; i++)
    {
      MBEntityHandle node_handle = array[i];
     
      rval = mbImpl->get_coords(&node_handle,1, coord);
      if(rval !=MB_SUCCESS) return rval;
      
      file << "v " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl; 
    }
  // write faces now
  // leave a blank line for cosmetics
  file << " \n";
  for (MBRange::const_iterator e = triangles.begin(); e != triangles.end(); ++e)
    {
      const MBEntityHandle* conn;
      int conn_len;
      rval = mbImpl->get_connectivity( *e, conn, conn_len );
      if (MB_SUCCESS != rval  ) 
	return rval;
      if ( 3!=conn_len) 
        return MB_INVALID_SIZE;
      file << "f ";
      for (int i = 0; i < conn_len; ++i)
      {
	int indexInArray = std::lower_bound( array, array + numNodes, conn[i] ) - array;
        file << indexInArray + 1 << " " ;
      }
      file << std::endl;
    }

  file.close();
  delete [] array;
  return MB_SUCCESS;
}

  
