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
// Filename      : WriteNCDF.hpp
//
// Purpose       : ExodusII writer
//
// Special Notes : Lots of code taken from verde implementation
//
// Creator       : Corey Ernst 
//
// Date          : 8/02
//
// Owner         : Corey Ernst 
//-------------------------------------------------------------------------

#ifndef WRITENCDF_HPP
#define WRITENCDF_HPP

#ifndef IS_BUILDING_MB
#error "WriteNCDF.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <string>

#include "netcdf.hh"
#include "MBForward.hpp"
#include "MBRange.hpp"
#include "ExoIIInterface.hpp"
#include "MBWriterIface.hpp"

class MBWriteUtilIface;

//! struct used to hold data for each block to be output in Exodus; used by
//! initialize_exodus_file to initialize the file header for increased speed
struct MaterialSetData
{
  int id;
  int number_elements;
  int number_nodes_per_element;
  int number_attributes;
  ExoIIElementType element_type;
  MBRange elements;
};

//! struct used to hold data for each nodeset to be output in Exodus; used by
//! initialize_exodus_file to initialize the file header for increased speed
struct DirichletSetData
{
  int id;
  int number_nodes;
  std::vector< MBEntityHandle > nodes;
  std::vector< double > node_dist_factors;
  
};

//! struct used to hold data for each sideset to be output in Exodus; used by
//! initialize_exodus_file to initialize the file header for increased speed
struct NeumannSetData
{
  int id;
  int number_elements;
  std::vector<MBEntityHandle> elements;
  std::vector<int> side_numbers;
  MBEntityHandle mesh_set_handle;
  std::vector< double > ss_dist_factors;
};


//! Output Exodus File for VERDE
class MB_DLL_EXPORT WriteNCDF : public MBWriterIface
{
 
public:
  
  static MBWriterIface* factory( MBInterface* );

   //! Constructor
   WriteNCDF(MBInterface *impl);

   //! Destructor
  virtual ~WriteNCDF();

    //! writes out an ExoII file
  MBErrorCode write_file(const char *exodus_file_name,
                         const bool overwrite,
                         const FileOptions& opts,
                          const MBEntityHandle *output_list,
                          const int num_sets,
                          const std::vector<std::string> &qa_records, 
                          const MBTag*,
                          int,
                          int user_dimension);
  
protected:

    //! number of dimensions in this exo file
  //int number_dimensions();

    //! open an ExoII file for writing
  MBErrorCode open_file(const char *file_name);

  //! contains the general information about a mesh
  struct ExodusMeshInfo
  {
    unsigned int num_dim;
    unsigned int num_nodes;
    unsigned int num_elements;
    unsigned int num_elementblocks;
    std::vector<std::string> qaRecords;
    MBRange nodes;
  };
  
private:

    //! interface instance
  MBInterface *mdbImpl;
  MBWriteUtilIface* mWriteIface;
  
    //! file name
  std::string exodusFile;
  NcFile *ncFile;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;

  //! Cached tags for reading.  Note that all these tags are defined when the
  //! core is initialized.
  MBTag mMaterialSetTag;
  MBTag mDirichletSetTag;
  MBTag mNeumannSetTag;
  MBTag mHasMidNodesTag;
  MBTag mGeomDimensionTag;
  MBTag mDistFactorTag;
  MBTag mGlobalIdTag;
  MBTag mQaRecordTag;

  MBTag mEntityMark;   //used to say whether an entity will be exported

  MBErrorCode gather_mesh_information(ExodusMeshInfo &mesh_info,
                                       std::vector<MaterialSetData> &block_info,
                                       std::vector<NeumannSetData> &sideset_info,
                                       std::vector<DirichletSetData> &nodeset_info,
                                       std::vector<MBEntityHandle> &blocks,
                                       std::vector<MBEntityHandle> &sidesets,
                                       std::vector<MBEntityHandle> &nodesets);

  MBErrorCode write_header(ExodusMeshInfo& mesh_info,
                            std::vector<MaterialSetData> &block_info,
                            std::vector<NeumannSetData> &sideset_info,
                            std::vector<DirichletSetData> &nodeset_info,
                            const char *filename);

  MBErrorCode initialize_exodus_file(ExodusMeshInfo &mesh_info,
                                      std::vector<MaterialSetData> &block_data,
                                      std::vector<NeumannSetData> & sideset_data,
                                      std::vector<DirichletSetData> & nodeset_data,
                                      const char* title_string,
                                      bool write_maps = true,
                                      bool write_sideset_distribution_factors = true);

  MBErrorCode write_qa_string(const char *string,
                               int record_number,
                               int record_position);

  MBErrorCode write_qa_records( std::vector<std::string> &qa_record_list);

  MBErrorCode write_nodes(int num_nodes, MBRange& nodes, int dimension );

  MBErrorCode write_elementblocks(std::vector<MaterialSetData> &block_data );

  MBErrorCode write_exodus_integer_variable(const char* variable_name,
                                                      int *variable_array,
                                                      int start_position,
                                                      int number_values);

  MBErrorCode write_global_node_order_map(int num_nodes, MBRange& nodes);
  MBErrorCode write_global_element_order_map(int num_elements);
  MBErrorCode write_element_order_map(int num_elements);

  MBErrorCode write_BCs(std::vector<NeumannSetData> &sidesets,
                         std::vector<DirichletSetData> &nodesets);

  MBErrorCode find_side_element_type( const int element_id, ExoIIElementType &type );

 /* MBErrorCode assign_block_ids_to_ssets(MBEntityHandle ss_handle,
                                         MB_MeshSet *ss_mesh_set);
                                         */
  //! free up allocated MBRanges
  void reset_block(std::vector<MaterialSetData> &block_info);

    //! recursive function; given a meshset handle, get the entities and put them
    //! on the right list, then call recursively for any contained meshsets, first
    //! checking for sense reversals
  MBErrorCode get_sideset_elems(MBEntityHandle sideset, int current_sense,
                                 MBRange &forward_elems, MBRange &reverse_elems);
  

  MBErrorCode get_valid_sides(MBRange &elems, ExodusMeshInfo &mesh_info, 
                               const int sense,
                               NeumannSetData &sideset_data);
  
    //! get the time and date in strings
  static void time_and_date(char* time_string, char* date_string);
};

#endif
