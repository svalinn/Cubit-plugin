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
 * \class WriteGmsh
 * \brief Export Gmsh files.
 * \author Jason Kraftcheck
 *
 * Known limitations:
 *  - Tag data and most sets cannot be saved.
 *  - For blocks, geometry sets, and parallel partitions, the
 *     sets are saved but the IDs may be lost.
 *  - File format does not support general polygons or polyhedra.
 *  - File format does not support higher-order volume elements 
 *     other than TET10 and HEX27.
 */

#ifndef WRITE_GMSH_HPP
#define WRITE_GMSH_HPP

#include "MBForward.hpp"
#include "MBWriterIface.hpp"
class MBWriteUtilIface;

#include <stdio.h>

class WriteGmsh : public MBWriterIface
{
 
public:
  
    //! factory method
  static MBWriterIface* factory( MBInterface* );

   //! Constructor
  WriteGmsh(MBInterface *impl);

   //! Destructor
  virtual ~WriteGmsh();
  
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

private:
                                       
    //! interface instance
  MBInterface *mbImpl;
  MBWriteUtilIface* mWriteIface;

};

#endif
