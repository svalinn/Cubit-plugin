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
// Filename      : ReadNASTRAN.hpp                        
//                                
// Purpose       : NASTRAN file reader
//                                             
// Creator       : Brandon Smith             
//                                   
// Date          : 08/2009                
//                                                  
//-------------------------------------------------------------------------     
                                    
#ifndef READNASTRAN_HPP                     
#define READNASTRAN_HPP              
                                     
#ifndef IS_BUILDING_MB                   
  #error "ReadNASTRAN.hpp isn't supposed to be included into an application"
#endif   

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "MBInterface.hpp"
#include "MBReaderIface.hpp"
#include "FileTokenizer.hpp"
#include "RangeMap.hpp"

class MBReadUtilIface;

class ReadNASTRAN : public MBReaderIface
{

public:
  // factory method
  static MBReaderIface* factory( MBInterface* );
  
  MBErrorCode load_file( const char                  *filename,
                         const MBEntityHandle        *file_set,
                         const FileOptions           &options,
                         const MBReaderIface::IDTag  *subset_list = 0,
                         int                         subset_list_length = 0,
                         const MBTag                 *file_id_tag = 0 );
  // constructor
  ReadNASTRAN(MBInterface* impl = NULL);

  // destructor
  virtual ~ReadNASTRAN();

  MBErrorCode read_tag_values( const char         *file_name,
			       const char         *tag_name,
			       const FileOptions  &opts,
			       std::vector<int>   &tag_values_out,
                               const IDTag        *subset_list,
		      	       int                subset_list_length );

protected:
  
private:  
  // read mesh interface
  MBReadUtilIface* readMeshIface;
  
  // MOAB Interface
  MBInterface* MBI;
  
  RangeMap<int, MBEntityHandle> nodeIdMap, elemIdMap;

  enum line_format { SMALL_FIELD,                     
                     LARGE_FIELD,                 
                     FREE_FIELD }; 

  MBErrorCode determine_line_format( const std::string line, 
                                     line_format &format );
  
  MBErrorCode tokenize_line( const std::string line, 
                             const line_format format,
                             std::vector<std::string> &tokens );  

  MBErrorCode determine_entity_type( const std::string token, MBEntityType &type); 

  MBErrorCode get_real( const std::string, double &real );

  MBErrorCode read_node(const std::vector<std::string> tokens, 
                        const bool           debug, 
                        double*              coord_arrays[3], 
                        int                  &node_id);

  MBErrorCode read_element(const std::vector<std::string> tokens, 
                           std::vector<MBRange>           &materials,
                           const MBEntityType             element_type,
                           const bool                     debug );

  MBErrorCode create_materials( const std::vector<MBRange> &materials );

  MBErrorCode assign_ids( const MBTag* file_id_tag );
};
#endif
