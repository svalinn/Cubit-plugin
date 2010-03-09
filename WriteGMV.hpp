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
// Filename      : WriteGMV.hpp
//
// Purpose       : Writer template
//
// Special Notes : 
//
// Creator       : Tim Tautges
//
// Date          : 2/04
//
//-------------------------------------------------------------------------

#ifndef WRITEGMV_HPP
#define WRITEGMV_HPP

#include "MBForward.hpp"
#include "MBWriterIface.hpp"
class MBWriteUtilIface;

//! Output Exodus File for VERDE
class MB_DLL_EXPORT WriteGMV : public MBWriterIface
{
 
public:

   //! Constructor
   WriteGMV(MBInterface *impl);

   //! Destructor
  virtual ~WriteGMV();

  static MBWriterIface* factory( MBInterface* );

  MBErrorCode write_file( const char* filename,
                          const bool overwite,
                          const FileOptions& opts,
                          const MBEntityHandle* output_sets,
                          const int num_output_sets,
                          const std::vector<std::string>& qa_list,
                          const MBTag* tag_list,
                          int num_tags,
                          int requested_dimension );

    //! writes out a mesh file
  MBErrorCode write_file(const char *file_name,
                         const MBEntityHandle output_set,
                         const int user_dimension = 3,
                         const bool mesh = true,
                         const bool poly_mesh = true);
  
protected:

private:

    //! interface instance
  MBInterface *mbImpl;
  MBWriteUtilIface* mWriteIface;
  
    //! Meshset Handle for the mesh that is currently being written
  MBEntityHandle mCurrentMeshHandle;

  //! Cached tags for reading.  Note that all these tags are defined when the
  //! core is initialized.
  MBTag mMaterialSetTag;
  MBTag mDirichletSetTag;
  MBTag mNeumannSetTag;
  MBTag mHasMidNodesTag;
  MBTag mGeomDimensionTag;
  MBTag mGlobalIdTag;

  static const char *gmvTypeNames[MBMAXTYPE];
  
  MBErrorCode local_write_mesh(const char *file_name,
                               const MBEntityHandle output_set,
                               const int user_dimension,
                               const bool mesh,
                               const bool poly_mesh);
};

#endif
