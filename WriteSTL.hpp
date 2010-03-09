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

/**
 * \class WriteSTL
 * \brief ASCII and Binary Stereo Lithography File writers.
 * \author Jason Kraftcheck
 *
 * This writer will write only the MBTRI elements in the mesh.  It
 * will not decompose other 2-D elements into triangles, nor will
 * it skin the mesh or do any other high-level operation to generate
 * triangles from 3-D elements.  
 *
 * Binary files will be written with a little-endian byte order by
 * default.  The byte order can be controlled by creating an integer
 * tag named "__STL_BYTE_ORDER" and setting the global/mesh value to
 * 0 for little endian or 1 for big endian.
 */

#ifndef WRITE_STL_HPP
#define WRITE_STL_HPP

#include "MBForward.hpp"
#include "MBWriterIface.hpp"

#include <stdio.h>

class MBWriteUtilIface;

class WriteSTL : public MBWriterIface
{
 
public:
  
    //! factory method forSTL writer
  static MBWriterIface* factory( MBInterface* );

   //! Constructor
  WriteSTL(MBInterface *impl);

   //! Destructor
  virtual ~WriteSTL();
  
    //! writes out a file
  MBErrorCode write_file(const char *file_name,
                         const bool overwrite,
                         const FileOptions& opts,
                         const MBEntityHandle *output_list,
                         const int num_sets,
                         const std::vector<std::string>& qa_list,
                         const MBTag* tag_list,
                         int num_tags,
                         int export_dimension);  

protected:
  
  enum ByteOrder { STL_BIG_ENDIAN, STL_LITTLE_ENDIAN, STL_UNKNOWN_BYTE_ORDER };
  
    //! Write list of triangles to an STL file.  
  MBErrorCode ascii_write_triangles( FILE* file,
                                     const char header[82],
                                     const MBRange& triangles,
                                     int precision );
    //! Write list of triangles to an STL file.  
  MBErrorCode binary_write_triangles( FILE* file,
                                      const char header[82],
                                      ByteOrder byte_order,
                                      const MBRange& triangles );

    //! Given an array of vertex coordinates for a triangle,
    //! pass back individual point coordinates as floats and 
    //! calculate triangle normal.
  MBErrorCode get_triangle_data( const double vtx_coords[9],
                                 float v1[3],
                                 float v2[3],
                                 float v3[3],
                                 float n[3] );
                                       
    //! interface instance
  MBInterface *mbImpl;
  MBWriteUtilIface* mWriteIface;
  
private:

    //! Construct 80-byte, null-terminated description string from
    //! qa_list.  Unused space in header will be null-char padded.
  MBErrorCode make_header( char header[82], const std::vector<std::string>& qa_list );
  
    //! Get triangles to write from input array of entity sets.  If
    //! no sets, gets all triangles.
  MBErrorCode get_triangles( const MBEntityHandle* set_array,
                             int set_array_length,
                             MBRange& triangles );  
  
    //! Open a file, respecting passed overwrite value and
    //! subclass-specified value for need_binary_io().
  FILE* open_file( const char* name, bool overwrite, bool binary );
};


#endif
