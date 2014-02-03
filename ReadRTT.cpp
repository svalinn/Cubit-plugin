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
 * \class ReadRTT
 * \brief ReadRTT based on ReadNASTRAN
 *
 * See: 
 *
 * \author Andrew Davis
 */



#include "ReadRTT.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <assert.h>
#include <cmath>

#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "Internals.hpp" // for MB_START_ID
#include "moab/Range.hpp"
#include "moab/FileOptions.hpp"
#include "FileTokenizer.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"

namespace moab {

ReaderIface* ReadRTT::factory( Interface* iface ) { 
  return new ReadRTT( iface );
}

// constructor
ReadRTT::ReadRTT(Interface* impl)
  : MBI(impl) {
    assert(NULL != impl);
    MBI->query_interface(readMeshIface);
    assert(NULL != readMeshIface);
}

// destructor
ReadRTT::~ReadRTT() {
  if (readMeshIface) {
    MBI->release_interface(readMeshIface);
    readMeshIface = 0;
  }
}

ErrorCode ReadRTT::read_tag_values( const char*        /*file_name*/,
                                    const char*        /*tag_name*/,
                                    const FileOptions& /*opts*/,
                                    std::vector<int>&  /*tag_values_out*/,
                                    const SubsetList*  /*subset_list*/ )
{
  return MB_NOT_IMPLEMENTED;
}

// load the file as called by the Interface function
ErrorCode ReadRTT::load_file(const char                      *filename, 
                             const EntityHandle            *, 
                             const FileOptions             &,
                             const ReaderIface::SubsetList *subset_list,
                             const Tag*                     file_id_tag) {
  // at this time there is no support for reading a subset of the file
  if (subset_list) {
    readMeshIface->report_error( "Reading subset of files not supported for RTT." );
    return MB_UNSUPPORTED_OPERATION;
  }

  nodeIdMap.clear();
  elemIdMap.clear();

  return MB_SUCCESS;
}


} // namespace moab
