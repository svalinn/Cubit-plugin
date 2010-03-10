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
// Filename      : WriteCCMIO.hpp
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

#ifndef WRITECCMIO_HPP
#define WRITECCMIO_HPP

#ifndef IS_BUILDING_MB
#error "WriteCCMIO.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <string>

#include "MBForward.hpp"
#include "MBRange.hpp"
#include "ExoIIInterface.hpp"
#include "MBWriterIface.hpp"
#include "ccmio.h"

class MBWriteUtilIface;

class MB_DLL_EXPORT WriteCCMIO : public MBWriterIface
{
 
public:

   //! Constructor
   WriteCCMIO(MBInterface *impl);

   //! Destructor
  virtual ~WriteCCMIO();
  
  static MBWriterIface* factory( MBInterface* );

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
  
//! struct used to hold data for each block to be output; used by
//! initialize_file to initialize the file header for increased speed
  struct MaterialSetData
  {
    int id;
    int number_elements;
    int number_nodes_per_element;
    int number_attributes;
    ExoIIElementType element_type;
    MBEntityType moab_type;
    MBRange *elements;
  };

//! struct used to hold data for each nodeset to be output; used by
//! initialize_file to initialize the file header for increased speed
  struct DirichletSetData
  {
    int id;
    int number_nodes;
    std::vector< MBEntityHandle > nodes;
    std::vector< double > node_dist_factors;
  
  };

//! struct used to hold data for each sideset to be output; used by
//! initialize_file to initialize the file header for increased speed
  struct NeumannSetData
  {
    int id;
    int number_elements;
    std::vector<MBEntityHandle> elements;
    std::vector<int> side_numbers;
    MBEntityHandle mesh_set_handle;
  };


protected:

    //! number of dimensions in this file
  //int number_dimensions();

    //! open a file for writing
  MBErrorCode open_file(const char *filename);

  //! contains the general information about a mesh
  class MeshInfo
  {
  public:
    unsigned int num_dim;
    unsigned int num_nodes;
    unsigned int num_elements;
    unsigned int num_matsets;
    unsigned int num_dirsets;
    unsigned int num_neusets;
    MBRange nodes;

    MeshInfo() 
        : num_dim(0), num_nodes(0), num_elements(0), num_matsets(0), 
          num_dirsets(0), num_neusets(0)
      {}
    
  };
  
private:

    //! interface instance
  MBInterface *mbImpl;
  MBWriteUtilIface* mWriteIface;
  
    //! file name
  std::string fileName;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;

  //! Cached tags for reading.  Note that all these tags are defined when the
  //! core is initialized.
  MBTag mMaterialSetTag;
  MBTag mDirichletSetTag;
  MBTag mNeumannSetTag;
  MBTag mPartitionSetTag;
  MBTag mHasMidNodesTag;
  MBTag mGlobalIdTag;
  MBTag mMatSetIdTag;

  MBTag mEntityMark;   //used to say whether an entity will be exported

  MBErrorCode gather_mesh_information(MeshInfo &mesh_info,
                                      std::vector<MaterialSetData> &matset_info,
                                      std::vector<NeumannSetData> &neuset_info,
                                      std::vector<DirichletSetData> &dirset_info,
                                      std::vector<MBEntityHandle> &matsets,
                                      std::vector<MBEntityHandle> &neusets,
                                      std::vector<MBEntityHandle> &dirsets);
  
  MBErrorCode initialize_file(MeshInfo &mesh_info);

  MBErrorCode write_nodes(CCMIOID rootID, const MBRange& nodes, 
                          const int dimension, int *&vgids);
  
    // get global ids for these entities; allocates gids and passes back,
    // caller is responsible for deleting
  MBErrorCode get_gids(const MBRange &ents, int *&gids,
                       int &minid, int &maxid);
  
  MBErrorCode write_matsets(MeshInfo &mesh_info, 
                            std::vector<MaterialSetData> &matset_data,
                            std::vector<NeumannSetData> &neuset_data,
                            MBRange &verts,
                            const int *vgids);
  
  MBErrorCode get_valid_sides(MBRange &elems, const int sense,
                              WriteCCMIO::NeumannSetData &neuset_data);
  
  void reset_matset(std::vector<MaterialSetData> &matset_info);
  
  MBErrorCode get_neuset_elems(MBEntityHandle neuset, int current_sense,
                               MBRange &forward_elems, MBRange &reverse_elems);
  
  MBErrorCode transform_coords(const int dimension, const int num_nodes, double *coords);

  MBErrorCode write_problem_description(CCMIOID rootID, CCMIOID stateID);
};

#endif
