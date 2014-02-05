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
// Filename      : ReadRTT.hpp                        
//                                
// Purpose       : RTT file reader
//                                             
// Creator       : Andrew Davis
//                                   
// Date          : 08/2009                
//                                                  
//-------------------------------------------------------------------------     
                                    
#ifndef READRTT_HPP                     
#define READRTT_HPP              
                                     
#ifndef IS_BUILDING_MB                   
  #error "ReadRTT.hpp isn't supposed to be included into an application"
#endif   

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "moab/Interface.hpp"
#include "moab/ReaderIface.hpp"
#include "FileTokenizer.hpp"
#include "moab/RangeMap.hpp"

namespace moab {

class ReadUtilIface;

class ReadRTT : public ReaderIface
{

public:
  // factory method
  static ReaderIface* factory( Interface* );
  
  ErrorCode load_file( const char* file_name,
                       const EntityHandle* file_set,
                       const FileOptions& opts,
                       const SubsetList* subset_list = 0,
                       const Tag* file_id_tag = 0 );
  // constructor
  ReadRTT(Interface* impl = NULL);

  // destructor
  virtual ~ReadRTT();

  ErrorCode read_tag_values( const char* file_name,
                             const char* tag_name,
                             const FileOptions& opts,
                             std::vector<int>& tag_values_out,
                             const SubsetList* subset_list = 0 );

protected:
  
private:  
  // read mesh interface
  ReadUtilIface* readMeshIface;
  // Moab Interface
  Interface* MBI;
  
};

} // namespace moab

#endif
