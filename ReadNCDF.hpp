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
// Filename      : ReadNCDF.hpp
//
// Purpose       : ExodusII reader
//
// Special Notes : Lots of code taken from verde implementation
//
// Creator       : Tim Tautges & Corey Ernst
//
// Date          : 3/02
//
// Owner         : Tim Tautges & Corey Ernst
//-------------------------------------------------------------------------

#ifndef READNCDF_HPP
#define READNCDF_HPP

#ifndef IS_BUILDING_MB
//#error "ReadNCDF.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <string>

#include "MBForward.hpp"
#include "MBReaderIface.hpp"
#include "ExoIIInterface.hpp"
#include "MBRange.hpp"

class MBReadUtilIface;

struct ReadBlockData
{
  int blockId;
  int startExoId; 
  MBEntityHandle startMBId; 
  int numElements;
  bool reading_in;
  ExoIIElementType elemType;
};

class NcFile;

//! Output Exodus File for VERDE
class ReadNCDF : public MBReaderIface
{
   
public:
  
  static MBReaderIface* factory( MBInterface* );
  
  void tokenize( const std::string& str,
                 std::vector<std::string>& tokens,
                 const char* delimiters );

    //! load an ExoII file
  MBErrorCode load_file( const char *exodus_file_name,
                         const MBEntityHandle* file_set,
                         const FileOptions& opts,
                         const MBReaderIface::IDTag* subset_list = 0,
                         int subset_list_length = 0,
                         const MBTag* file_id_tag = 0 );

  MBErrorCode read_tag_values( const char* file_name,
                               const char* tag_name,
                               const FileOptions& opts,
                               std::vector<int>& tag_values_out,
                               const IDTag* subset_list = 0,
                               int subset_list_length = 0 );
  
   //! Constructor
   ReadNCDF(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadNCDF();

  //update the coords for deformed mesh according to FileOptions
  MBErrorCode update(const char *exodus_file_name, const FileOptions& opts,
                     const int num_blocks, const int *blocks_to_load,
                     const MBEntityHandle file_set );

private:

  MBReadUtilIface* readMeshIface;

  bool dimension_exists(const char *attrib_name);
  
  void reset();

    //! read the header from the ExoII file
  MBErrorCode read_exodus_header();
  
    //! read the nodes
  MBErrorCode read_nodes(const MBTag* file_id_tag);
  
    //! read block headers, containing info about element type, number, etc.
  MBErrorCode read_block_headers(const int *blocks_to_load,
                                  const int num_blocks);
  
    //! read the element blocks
  MBErrorCode read_elements(const MBTag* file_id_tag);
  
    //! read in the global element ids
  MBErrorCode read_global_ids();

    //! read the nodesets into meshsets
  MBErrorCode read_nodesets();
  
    //! read the sidesets (does nothing for now)
  MBErrorCode read_sidesets();

    //! exodus file bound to this object
  int exodus_file();

    //! number of dimensions in this exo file
  int number_dimensions();

  //! map a character exodusII element type to a TSTT type & topology
  MBErrorCode get_type(char *exo_element_type,
                        MBEntityType &elem_type);
 
  MBErrorCode get_type(MBEntityType &elem_type,
                        std::string &exo_element_type);

  /* 
  int get_int_tag(const MB_MeshSet *this_ms,
                  const TagHandle tag_id);
 */

  //qa record stuff 
  MBErrorCode read_qa_records(MBEntityHandle file_set);
  MBErrorCode read_qa_information( std::vector<std::string> &qa_record_list);

  MBErrorCode read_qa_string(char *string,
                              int record_number,
                              int record_position); 

  MBErrorCode create_ss_elements( int *element_ids, int *side_list,
                                   int num_sides, int num_dist_factors,
                                   std::vector<MBEntityHandle> &entities_to_add,
                                   std::vector<MBEntityHandle> &reverse_entities,
                                   std::vector<double> &dist_factor_vector,
                                   int ss_seq_id);

  MBErrorCode find_side_element_type( const int element_id, ExoIIElementType &type, 
                                       ReadBlockData &block_data, int &df_index, int side_id );
  
 /* MBErrorCode assign_block_ids_to_ssets(MBEntityHandle ss_handle,
                                         MB_MeshSet *ss_mesh_set);
                                         */

  //! creates an element with the given connectivity
  MBErrorCode create_sideset_element( const std::vector<MBEntityHandle>&, MBEntityType, MBEntityHandle&);

  int get_number_nodes( MBEntityHandle handle );



  //------------member variables ------------//

    //! interface instance
  MBInterface* mdbImpl;
  
  NcFile *ncFile;        // netcdf/exodus file

  int CPU_WORD_SIZE;
  int IO_WORD_SIZE;

    //! int to offset vertex ids with
  int vertexOffset;

    //! file name
  std::string exodusFile;

    //! number of nodes in the current exoII file
  int numberNodes_loading;

    //! number of elements in the current exoII file
  int numberElements_loading;

    //! number of blocks in the current exoII file
  int numberElementBlocks_loading; 

    //! number of nodesets in the current exoII file
  int numberNodeSets_loading; 

    //! number of sidesets in the current exoII file
  int numberSideSets_loading; 

  int numberDimensions_loading;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;

  //keeps track of the exodus ids of the elements and nodes just loaded
  std::vector<char> nodesInLoadedBlocks;
  //note- vector<bool> has limited capabilities

  //vector of blocks that are loading 
  std::vector< ReadBlockData > blocksLoading;

  //! Cached tags for reading.  Note that all these tags are defined when the
  //! core is initialized.
  MBTag mMaterialSetTag;
  MBTag mDirichletSetTag;
  MBTag mNeumannSetTag;
  MBTag mHasMidNodesTag;
  MBTag mDistFactorTag;
  MBTag mGlobalIdTag;
  MBTag mQaRecordTag;

  int max_line_length, max_str_length;

    //! range of entities in initial mesh, before this read
  MBRange initRange;
};

// inline functions
inline int ReadNCDF::number_dimensions() 
{
   return numberDimensions_loading;
}


#endif




