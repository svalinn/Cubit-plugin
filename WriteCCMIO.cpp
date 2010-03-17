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


#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif


#include "WriteCCMIO.hpp"
#include "ccmio.h"
#include "ccmioutility.h"
#include "ccmiocore.h"
#include <utility>
#include <algorithm>
#include <time.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <algorithm>

#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include "assert.h"
#include "Internals.hpp"
#include "ExoIIUtil.hpp"
#include "moab/MBTagConventions.hpp"
#include "moab/WriteUtilIface.hpp"

namespace moab {

static char const kStateName[] = "default";

static const int ccm_types[] = {
    1,   // MBVERTEX
    2,   // MBEDGE      
    -1,  // MBTRI
    -1,  // MBQUAD
    -1,  // MBPOLYGON
    13,  // MBTET
    14,  // MBPYRAMID
    12,  // MBPRISM
    -1,  // MBKNIFE
    11,  // MBHEX
    255  // MBPOLYHEDRON
};

#define INS_ID(stringvar, prefix, id)           \
    sprintf(stringvar, prefix, id)

WriterIface* WriteCCMIO::factory( Interface* iface )
{ return new WriteCCMIO( iface ); }

WriteCCMIO::WriteCCMIO(Interface *impl) 
        : mbImpl(impl), mCurrentMeshHandle(0)
{
  assert(impl != NULL);

  void* ptr = 0;
  impl->query_interface( "WriteUtilIface", &ptr );
  mWriteIface = reinterpret_cast<WriteUtilIface*>(ptr);

    // initialize in case tag_get_handle fails below
    //! get and cache predefined tag handles
  int dum_val = 0;
  ErrorCode result = impl->tag_get_handle(MATERIAL_SET_TAG_NAME,  mMaterialSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(MATERIAL_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mMaterialSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(DIRICHLET_SET_TAG_NAME, mDirichletSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(DIRICHLET_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mDirichletSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(NEUMANN_SET_TAG_NAME, mNeumannSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(NEUMANN_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mNeumannSetTag,
                              &dum_val);

#ifdef USE_MPI  
  result = impl->tag_get_handle(PARALLEL_PARTITION_TAG_NAME, mPartitionSetTag);
    // no need to check result, if it's not there, we don't create one
#endif
  
  result = impl->tag_get_handle(HAS_MID_NODES_TAG_NAME, mHasMidNodesTag);
  if (MB_TAG_NOT_FOUND == result) {
    int dum_val_array[] = {0, 0, 0, 0};
    result = impl->tag_create(HAS_MID_NODES_TAG_NAME, 4*sizeof(int), MB_TAG_SPARSE, mHasMidNodesTag,
                              dum_val_array);
  }
  
  result = impl->tag_get_handle(GLOBAL_ID_TAG_NAME, mGlobalIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mGlobalIdTag,
                              &dum_val);
  
  dum_val = -1;
  result = impl->tag_get_handle("__matSetIdTag", mMatSetIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create("__matSetIdTag", sizeof(int), MB_TAG_DENSE, mMatSetIdTag,
                              &dum_val);
  

  impl->tag_create("WriteCCMIO element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

}

WriteCCMIO::~WriteCCMIO() 
{
  std::string iface_name = "WriteUtilIface";
  mbImpl->release_interface(iface_name, mWriteIface);

  mbImpl->tag_delete(mEntityMark);

}

void WriteCCMIO::reset_matset(std::vector<WriteCCMIO::MaterialSetData> &matset_info)
{
  std::vector<WriteCCMIO::MaterialSetData>::iterator iter;
  
  for (iter = matset_info.begin(); iter != matset_info.end(); iter++)
  {
    delete (*iter).elements;
  }
}

ErrorCode WriteCCMIO::write_file(const char *file_name, 
                                   const bool /* overwrite (commented out to remove warning) */,
                                   const FileOptions& opts,
                                   const EntityHandle *ent_handles,
                                   const int num_sets,
                                   const std::vector<std::string>&,
                                   const Tag* ,
                                   int ,
                                   int )
{
  assert(0 != mMaterialSetTag &&
         0 != mNeumannSetTag &&
         0 != mDirichletSetTag);

    // check the file name
  if (NULL == strstr(file_name, ".ccmio"))
    return MB_FAILURE;

  std::vector<EntityHandle> matsets, dirsets, neusets, partsets, entities;

  fileName = file_name;
  
    // separate into material sets, dirichlet sets, neumann sets

  if (num_sets == 0) {
      // default to all defined sets
    Range this_range;
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mMaterialSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(matsets));
    this_range.clear();
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mDirichletSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(dirsets));
    this_range.clear();
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mNeumannSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(neusets));
    if (mPartitionSetTag) {
      this_range.clear();
      mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mPartitionSetTag, NULL, 1, this_range);
      std::copy(this_range.begin(), this_range.end(), std::back_inserter(partsets));
    }
  }
  
  else {
    int dummy;
    for (const EntityHandle *iter = ent_handles; iter < ent_handles+num_sets; iter++) 
    {
      if (MB_SUCCESS == mbImpl->tag_get_data(mMaterialSetTag, &(*iter), 1, &dummy))
        matsets.push_back(*iter);
      else if (MB_SUCCESS == mbImpl->tag_get_data(mDirichletSetTag, &(*iter), 1, &dummy))
        dirsets.push_back(*iter);
      else if (MB_SUCCESS == mbImpl->tag_get_data(mNeumannSetTag, &(*iter), 1, &dummy))
        neusets.push_back(*iter);
      else if (mPartitionSetTag &&
               MB_SUCCESS == mbImpl->tag_get_data(mPartitionSetTag, &(*iter), 1, &dummy))
        partsets.push_back(*iter);
    }
  }

  ErrorCode result = mbImpl->tag_get_handle(HAS_MID_NODES_TAG_NAME, mHasMidNodesTag);
  if (MB_TAG_NOT_FOUND == result) {
    int dum_val_array[] = {0, 0, 0, 0};
    result = mbImpl->tag_create(HAS_MID_NODES_TAG_NAME, 4*sizeof(int), MB_TAG_SPARSE, mHasMidNodesTag,
                                dum_val_array);
  }
  
    // if there is nothing to write just return.
  if (matsets.empty() && dirsets.empty() && neusets.empty() && partsets.empty())
    return MB_FILE_WRITE_ERROR;

  std::vector<WriteCCMIO::MaterialSetData> matset_info;
  std::vector<WriteCCMIO::DirichletSetData> dirset_info;
  std::vector<WriteCCMIO::NeumannSetData> neuset_info;

  MeshInfo mesh_info;
  
  matset_info.clear();
  if(gather_mesh_information(mesh_info, matset_info, neuset_info, dirset_info,
                             matsets, neusets, dirsets) != MB_SUCCESS)
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }


    // try to open the file after gather mesh info succeeds
  CCMIOSize_t i = CCMIOSIZEC(0);
  CCMIOID stateID, processorID, rootID;
  CCMIOError error = kCCMIONoErr;

  CCMIOOpenFile(&error, file_name, kCCMIOWrite, &rootID);
  if(kCCMIONoErr != error)
  {
    mWriteIface->report_error("Cannot open %s", file_name);
    return MB_FAILURE;
  }

    // Create a new state (or re-use an existing one).
  CCMIONewState(&error, rootID, kStateName, NULL, NULL, &stateID);
  if (kCCMIONoErr != error) {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

//  for (; i < CCMIOSIZEC(partsets.size()); i++) {
  if (CCMIONextEntity(NULL, stateID, kCCMIOProcessor, &i, &processorID) != kCCMIONoErr)
    CCMIONewEntity(&error, stateID, kCCMIOProcessor, NULL, &processorID);

    // Get rid of any data that may be in this processor (if the state was
    // not new).
  else
    CCMIOClearProcessor(&error, stateID, processorID, TRUE, TRUE, TRUE, TRUE,
                        TRUE);
//  }

  int *vgids;
  if( write_nodes(rootID, mesh_info.nodes, mesh_info.num_dim, vgids) 
      != MB_SUCCESS )
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  if( write_matsets(mesh_info, matset_info, neuset_info, mesh_info.nodes, vgids) )
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  if (write_problem_description(rootID, stateID)) {
    return MB_FAILURE;
  }
  
  CCMIOCloseFile(&error, rootID);

  if (error != kCCMIONoErr)
  {
    return MB_FAILURE;
  }

  // The CCMIO library uses ADF to store the actual data.  Unfortunately,
  // ADF leaks disk space;  deleting a node does not recover all the disk
  // space.  Now that everything is successfully written it might be useful
  // to call CCMIOCompress() here to ensure that the file is as small as
  // possible.  Please see the Core API documentation for caveats on its
  // usage.
  if (CCMIOCompress(NULL, const_cast<char*>(file_name)) != kCCMIONoErr)
  {
    std::cout << "Error compressing file.  Check that you have "
              << "adequate disk space " << std::endl << "and that you have write "
              << "permission to the current directory." << std::endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode WriteCCMIO::write_problem_description(CCMIOID rootID, CCMIOID stateID) 
{
  // Write out a dummy problem description.  If we happen to know that
  // there already is a problem description previously recorded that
  // is valid we could skip this step.
  CCMIOID problem, constants, id;
  CCMIOError error = kCCMIONoErr;

  CCMIONewEntity(&error, rootID, kCCMIOProblemDescription, "Dummy description",
                 &problem);
  CCMIONewIndexedEntity(&error, problem, kCCMIOCellType, 1, "Dummy celltypes", &id);
  CCMIOWriteOptstr(&error, id, "MaterialType", "solid");
  CCMIONewIndexedEntity(&error, problem, kCCMIOCellType, 2, "Dummy celltypes", &id);
  CCMIOWriteOptstr(&error, id, "MaterialType", "solid");

  CCMIONewEntity(&error, problem, kCCMIOModelConstants, "Constant values",
                 &constants);
  CCMIOWriteOptf(&error, constants, "Gravity", 9.82);
  CCMIOWriteOptf(&error, constants, "B.P. of water", 373);

  // We have problem description recorded but our state does not know
  // about it.  So tell the state that it has a problem description.
  CCMIOWriteState(&error, stateID, problem, "Example state");

  return MB_SUCCESS;
}


ErrorCode WriteCCMIO::gather_mesh_information(MeshInfo &mesh_info,
                                                std::vector<WriteCCMIO::MaterialSetData> &matset_info,
                                                std::vector<WriteCCMIO::NeumannSetData> &neuset_info,
                                                std::vector<WriteCCMIO::DirichletSetData> &dirset_info,
                                                std::vector<EntityHandle> &matsets,
                                                std::vector<EntityHandle> &neusets,
                                                std::vector<EntityHandle> &dirsets)
{

  std::vector<EntityHandle>::iterator vector_iter, end_vector_iter;

  mesh_info.num_nodes = 0;
  mesh_info.num_elements = 0;
  mesh_info.num_matsets = 0;
  
  int id = 0;

  vector_iter= matsets.begin();
  end_vector_iter = matsets.end();

  mesh_info.num_matsets = matsets.size();

  std::vector<EntityHandle> parent_meshsets;

    // clean out the bits for the element mark
  mbImpl->tag_delete(mEntityMark);
  mbImpl->tag_create("WriteCCMIO element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

  int highest_dimension_of_element_matsets = 0;

  for(vector_iter = matsets.begin(); vector_iter != matsets.end(); vector_iter++)
  {
       
    WriteCCMIO::MaterialSetData matset_data;
    matset_data.elements = new Range;

      //for the purpose of qa records, get the parents of these matsets 
    if( mbImpl->get_parent_meshsets( *vector_iter, parent_meshsets ) != MB_SUCCESS )
      return MB_FAILURE;

      // get all Entity Handles in the mesh set
    Range dummy_range;
    mbImpl->get_entities_by_handle(*vector_iter, dummy_range, true );

      // find the dimension of the last entity in this range
    Range::iterator entity_iter = dummy_range.end();
    entity_iter = dummy_range.end();
    entity_iter--;
    int this_dim = CN::Dimension(TYPE_FROM_HANDLE(*entity_iter));
    entity_iter = dummy_range.begin();
    while (entity_iter != dummy_range.end() &&
           CN::Dimension(TYPE_FROM_HANDLE(*entity_iter)) != this_dim)
      entity_iter++;
    
    if (entity_iter != dummy_range.end())
      std::copy(entity_iter, dummy_range.end(), range_inserter(*(matset_data.elements)));

    assert(matset_data.elements->begin() == matset_data.elements->end() ||
           CN::Dimension(TYPE_FROM_HANDLE(*(matset_data.elements->begin()))) == this_dim);
    
      // get the matset's id
    if(mbImpl->tag_get_data(mMaterialSetTag, &(*vector_iter), 1, &id) != MB_SUCCESS ) {
      mWriteIface->report_error("Couldn't get matset id from a tag for an element matset.");
      return MB_FAILURE;
    }
    
    matset_data.id = id; 
    matset_data.number_attributes = 0;
 
      // iterate through all the elements in the meshset
    Range::iterator elem_range_iter, end_elem_range_iter;
    elem_range_iter = matset_data.elements->begin();
    end_elem_range_iter = matset_data.elements->end();

      // get the entity type for this matset, verifying that it's the same for all elements
      // THIS ASSUMES HANDLES SORT BY TYPE!!!
    EntityType entity_type = TYPE_FROM_HANDLE(*elem_range_iter);
    end_elem_range_iter--;
    if (entity_type != TYPE_FROM_HANDLE(*(end_elem_range_iter++))) {
      mWriteIface->report_error("Entities in matset %i not of common type", id);
      return MB_FAILURE;
    }

    int dimension = CN::Dimension(entity_type);

    if( dimension > highest_dimension_of_element_matsets )
      highest_dimension_of_element_matsets = dimension;

    matset_data.moab_type = mbImpl->type_from_handle(*(matset_data.elements->begin()));
    if (MBMAXTYPE == matset_data.moab_type) return MB_FAILURE;
    
    std::vector<EntityHandle> tmp_conn;
    mbImpl->get_connectivity(&(*(matset_data.elements->begin())), 1, tmp_conn);
    matset_data.element_type = 
        ExoIIUtil::get_element_type_from_num_verts(tmp_conn.size(), entity_type, dimension);
    
    if (matset_data.element_type == EXOII_MAX_ELEM_TYPE) {
      mWriteIface->report_error("Element type in matset %i didn't get set correctly", id);
      return MB_FAILURE;
    }
    
    matset_data.number_nodes_per_element = ExoIIUtil::VerticesPerElement[matset_data.element_type];

      // number of nodes for this matset
    matset_data.number_elements = matset_data.elements->size();

      // total number of elements
    mesh_info.num_elements += matset_data.number_elements;

      // get the nodes for the elements
    mWriteIface->gather_nodes_from_elements(*matset_data.elements, mEntityMark, mesh_info.nodes);

    if(!neusets.empty())
    {
        // if there are neusets, keep track of which elements are being written out
      for(Range::iterator iter = matset_data.elements->begin(); 
          iter != matset_data.elements->end(); ++iter)
      {
        unsigned char bit = 0x1;
        mbImpl->tag_set_data(mEntityMark, &(*iter), 1, &bit);
      }
    }

    matset_info.push_back( matset_data );
  
  }
 

    //if user hasn't entered dimension, we figure it out
  if( mesh_info.num_dim == 0 )
  {
      //never want 1 or zero dimensions
    if( highest_dimension_of_element_matsets < 2 )
      mesh_info.num_dim = 3;
    else
      mesh_info.num_dim = highest_dimension_of_element_matsets;
  }

  Range::iterator range_iter, end_range_iter;
  range_iter = mesh_info.nodes.begin();
  end_range_iter = mesh_info.nodes.end();

  mesh_info.num_nodes = mesh_info.nodes.size(); 

    //------dirsets--------
  
  vector_iter= dirsets.begin();
  end_vector_iter = dirsets.end();

  for(; vector_iter != end_vector_iter; vector_iter++)
  {
    
    WriteCCMIO::DirichletSetData dirset_data;
    dirset_data.id = 0;
    dirset_data.number_nodes = 0;

      // get the dirset's id
    if(mbImpl->tag_get_data(mDirichletSetTag,&(*vector_iter), 1,&id) != MB_SUCCESS) {
      mWriteIface->report_error("Couldn't get id tag for dirset %i", id);
      return MB_FAILURE;
    }
    
    dirset_data.id = id; 

    std::vector<EntityHandle> node_vector;
      //get the nodes of the dirset that are in mesh_info.nodes
    if( mbImpl->get_entities_by_handle(*vector_iter, node_vector, true) != MB_SUCCESS ) {
      mWriteIface->report_error("Couldn't get nodes in dirset %i", id);
      return MB_FAILURE;
    }

    std::vector<EntityHandle>::iterator iter, end_iter;
    iter = node_vector.begin();
    end_iter= node_vector.end();
 
    int j=0; 
    unsigned char node_marked = 0;
    ErrorCode result;
    for(; iter != end_iter; iter++)
    {
      if (TYPE_FROM_HANDLE(*iter) != MBVERTEX) continue;
      result = mbImpl->tag_get_data(mEntityMark, &(*iter), 1, &node_marked);
      if (MB_SUCCESS != result) {
        mWriteIface->report_error("Couldn't get mark data.");
        return result;
      }
      
      if(node_marked == 0x1) dirset_data.nodes.push_back( *iter );    
      j++;
    } 
    
    dirset_data.number_nodes = dirset_data.nodes.size(); 
    dirset_info.push_back( dirset_data );
  }

    //------neusets--------
  vector_iter= neusets.begin();
  end_vector_iter = neusets.end();

  for(; vector_iter != end_vector_iter; vector_iter++)
  {
    WriteCCMIO::NeumannSetData neuset_data;

      // get the neuset's id
    if(mbImpl->tag_get_data(mNeumannSetTag,&(*vector_iter), 1,&id) != MB_SUCCESS)
      return MB_FAILURE;

    neuset_data.id = id; 
    neuset_data.mesh_set_handle = *vector_iter; 
 
      //get the sides in two lists, one forward the other reverse; starts with forward sense
      // by convention
    Range forward_elems, reverse_elems;
    if(get_neuset_elems(*vector_iter, 0, forward_elems, reverse_elems) == MB_FAILURE)
      return MB_FAILURE;

    ErrorCode result = get_valid_sides(forward_elems, 1, neuset_data);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get valid sides data.");
      return result;
    }
    result = get_valid_sides(reverse_elems, -1, neuset_data);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get valid sides data.");
      return result;
    }
    
    neuset_data.number_elements = neuset_data.elements.size(); 
    neuset_info.push_back( neuset_data );
  }

  return MB_SUCCESS;
}

ErrorCode WriteCCMIO::get_valid_sides(Range &elems, const int sense,
                                        WriteCCMIO::NeumannSetData &neuset_data) 
{
    // this is where we see if underlying element of side set element is included in output 

  unsigned char element_marked = 0;
  ErrorCode result;
  for(Range::iterator iter = elems.begin(); iter != elems.end(); iter++)
  {
      // should insert here if "side" is a quad/tri on a quad/tri mesh
    result = mbImpl->tag_get_data(mEntityMark, &(*iter), 1, &element_marked);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get mark data.");
      return result;
    }
    
    if(element_marked == 0x1)
    {
      neuset_data.elements.push_back( *iter );

        // TJT TODO: the sense should really be # edges + 1or2
      neuset_data.side_numbers.push_back((sense == 1 ? 1 : 2));
    }
    else //then "side" is probably a quad/tri on a hex/tet mesh
    {
      std::vector<EntityHandle> parents;
      int dimension = CN::Dimension( TYPE_FROM_HANDLE(*iter));

        //get the adjacent parent element of "side"
      if( mbImpl->get_adjacencies( &(*iter), 1, dimension+1, false, parents) != MB_SUCCESS ) {
        mWriteIface->report_error("Couldn't get adjacencies for neuset.");
        return MB_FAILURE;
      }
       
      if(!parents.empty())     
      {
          //make sure the adjacent parent element will be output
        for(unsigned int k=0; k<parents.size(); k++)
        {
          result = mbImpl->tag_get_data(mEntityMark, &(parents[k]), 1, &element_marked);
          if (MB_SUCCESS != result) {
            mWriteIface->report_error("Couldn't get mark data.");
            return result;
          }
        
          int side_no, this_sense, this_offset;
          if(element_marked == 0x1 &&
             mbImpl->side_number(parents[k], *iter, side_no, 
                                 this_sense, this_offset) == MB_SUCCESS &&
             this_sense == sense) {
            neuset_data.elements.push_back(parents[k]);
            neuset_data.side_numbers.push_back(side_no+1);
            break;
          }
        }
      }
      else
      {
        mWriteIface->report_error("No parent element exists for element in neuset %i", neuset_data.id);
        return MB_FAILURE;
      }
    }
  }

  return MB_SUCCESS;
}

ErrorCode WriteCCMIO::get_gids(const Range &ents, int *&gids,
                                 int &minid, int &maxid) 
{
  int num_ents = ents.size();
  gids = new int[num_ents];
  ErrorCode result = mbImpl->tag_get_data(mGlobalIdTag, ents, &gids[0]);
  if (MB_SUCCESS != result) {
    mWriteIface->report_error("Couldn't get global id data.");
    return result;
  }
  minid = *std::min_element(gids, gids+num_ents);
  maxid = *std::max_element(gids, gids+num_ents);
  if (0 == minid) {
      // gids need to be assigned
    for (int i = 1; i <= num_ents; i++) gids[i] = i;
    result = mbImpl->tag_set_data(mGlobalIdTag, ents, &gids[0]);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't set global id data.");
      return result;
    }
    maxid = num_ents;
  }
  return MB_SUCCESS;
}

ErrorCode WriteCCMIO::write_nodes(CCMIOID rootID, const Range& nodes, 
                                    const int dimension, int *&vgids)
{
    // get/write map (global ids) first
  const int num_nodes = nodes.size();
  int minid, maxid;
  ErrorCode result = get_gids(nodes, vgids, minid, maxid);
  if (MB_SUCCESS != result) return result;
  
  CCMIOID mapID;
  CCMIOError error;
  CCMIONewEntity(&error, rootID, kCCMIOMap, "Vertex map", &mapID);
  CCMIOWriteMap(&error, mapID, CCMIOSIZEC(nodes.size()),
                CCMIOSIZEC(maxid), vgids,
                CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));

    // get the vertex locations
  double *coords = new double[dimension*num_nodes];
  std::vector<double*> coord_arrays(3);
  coord_arrays[0] = coords;
  coord_arrays[1] = coords+num_nodes;
  coord_arrays[2] = (dimension == 3 ? coords+2*num_nodes : NULL);
  result = mWriteIface->get_node_arrays(dimension, num_nodes, nodes,
                                        mGlobalIdTag, 0, coord_arrays);
  if(result != MB_SUCCESS)
  {
    delete [] coords;
    return result;
  }
  
    // transform coordinates, if necessary
  result = transform_coords(dimension, num_nodes, coords);
  if(result != MB_SUCCESS)
  {
    delete [] coords;
    return result;
  }
  

    // write the vertices
  CCMIOID verticesID;
  CCMIONewEntity(&error, rootID, kCCMIOVertices, "Vertices", &verticesID);
  CCMIOWriteVerticesd(&error, verticesID,
                     CCMIOSIZEC(dimension*num_nodes), 1.0, mapID, coords,
                     CCMIOINDEXC(1), CCMIOINDEXC(dimension));
  if (error != kCCMIONoErr) {
    mWriteIface->report_error("CCMIOWriteVertices failed.");
    return result;
  }
    
    // clean up
  delete [] coords;

  return MB_SUCCESS;
}

ErrorCode WriteCCMIO::transform_coords(const int dimension, const int num_nodes, double *coords) 
{
  Tag trans_tag;
  ErrorCode result = mbImpl->tag_get_handle( MESH_TRANSFORM_TAG_NAME, trans_tag);
  if( result == MB_TAG_NOT_FOUND ) return MB_SUCCESS;
  double trans_matrix[16]; 
  result = mbImpl->tag_get_data( trans_tag, NULL, 0, trans_matrix ); 
  if (MB_SUCCESS != result) {
    mWriteIface->report_error("Couldn't get transform data.");
    return result;
  }
      
  double *tmp_coords = coords;
  for( int i=0; i<num_nodes; i++, tmp_coords += 1) {
    double vec1[3] = {0.0, 0.0, 0.0};
    for( int row=0; row<3; row++ ) {
      vec1[row] += ( trans_matrix[ (row*4)+0 ] * coords[0]);
      vec1[row] += ( trans_matrix[ (row*4)+1 ] * coords[num_nodes]);
      if (3 == dimension) vec1[row] += ( trans_matrix[ (row*4)+2 ] * coords[2*num_nodes]);
    }

    coords[0] = vec1[0];
    coords[num_nodes] = vec1[1];
    coords[2*num_nodes] = vec1[2];
  }

  return MB_SUCCESS;
}

ErrorCode WriteCCMIO::write_matsets(MeshInfo & /* mesh_info (commented out to remove warning) */,
                                      std::vector<WriteCCMIO::MaterialSetData> &matset_data,
                                      std::vector<WriteCCMIO::NeumannSetData> &/* neuset_data  */,
                                        // (commented out to remove warning)
                                      Range &verts,
                                      const int *vgids)
{
  std::vector<int> connect;
  ErrorCode result;
  CCMIOID rootID, cellMapID, topologyID, id;
  
    // don't usually have anywhere near 31 nodes per element
  connect.reserve(31);
  Range::const_iterator rit;

  Range all_elems;
  for (unsigned int i = 0; i < matset_data.size(); i++)
    all_elems.merge(*(matset_data[i].elements));

  
  const int num_elems = all_elems.size();
  int *egids, minid, maxid;
  result = get_gids(all_elems, egids, minid, maxid);
  if (MB_SUCCESS != result) return result;

    // Write the cells
  CCMIOError error;
  CCMIONewEntity(&error, rootID, kCCMIOMap, "Cell map", &cellMapID);
  CCMIOWriteMap(&error, cellMapID, CCMIOSIZEC(all_elems.size()),
                CCMIOSIZEC(maxid), egids,
                CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CCMIONewEntity(&error, rootID, kCCMIOTopology, "Topology", &topologyID);
  CCMIONewEntity(&error, topologyID, kCCMIOCells, "Cells", &id);

    // get cell types
  int *ctypes = new int[all_elems.size()];
  int i = 0;
  rit = all_elems.begin();
  for (; i < num_elems; i++, rit++) {
    ctypes[i] = ccm_types[mbImpl->type_from_handle(*rit)];
    assert(-1 != ctypes[i]);
  }

  CCMIOWriteCells(&error, id, cellMapID, ctypes,
                  CCMIOINDEXC(1), CCMIOINDEXC(num_elems));
  delete [] ctypes;

    // Write the faces
    // first, allocate a tag of length 6 (max # faces per region, except
    // for polyhedra)
  Tag mark_tag;
  short int def_val = 0;
  result = mbImpl->tag_create("__mark", 1, MB_TAG_DENSE, MB_TYPE_OPAQUE, 
                              mark_tag, &def_val);
  if (MB_SUCCESS != result) {
    mWriteIface->report_error("Couldn't create mark tag.");
    return result;
  }
  
    // now faces
  unsigned char markt;
  std::vector<EntityHandle> tmp_face_cells, storage;
  std::vector<int> iface_connect, iface_cells;
  std::vector<int> eface_connect, eface_cells;
  EntityHandle tmp_connect[CN::MAX_NODES_PER_ELEMENT]; // tmp connect vector
  const EntityHandle *connectc; int num_connectc; // cell connectivity
  const EntityHandle *connectf; int num_connectf; // face connectivity
  i = 0;
  rit = all_elems.begin();
  for (; i < num_elems; i++, rit++) {
      // if not polyh, get mark
    EntityType etype = TYPE_FROM_HANDLE(*rit);
    if (MBPOLYHEDRON != etype && MBPOLYGON != etype) {
      result = mbImpl->tag_get_data(mark_tag, &(*rit), 1, &markt);
      if (MB_SUCCESS != result) {
        mWriteIface->report_error("Couldn't get mark data.");
        return result;
      }
    }

      // get connect
    result = mbImpl->get_connectivity(*rit, connectc, num_connectc, false, &storage);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get entity connectivity.");
      return result;
    }

      // if polyh, write faces directly
    bool is_polyh = (MBPOLYHEDRON == etype);

    int num_faces = CN::NumSubEntities(etype, 2);
      // for each face (from CN)
    for (int f = 0; f < num_faces; f++) {
        // if this face marked, skip
      if (!is_polyh && ((markt >> f) & 0x1)) continue;
        
        // get face connect
      if (!is_polyh) {
          // (from CN)
        CN::SubEntityConn(connectc, etype, 2, f, tmp_connect, num_connectf);
        connectf = tmp_connect;
      }
      else {
          // get face connectivity
        result = mbImpl->get_connectivity(connectc[f], connectf, num_connectf, false);
        if (MB_SUCCESS != result) {
          mWriteIface->report_error("Couldn't get polyhedron connectivity.");
          return result;
        }
      }
        
        // get adj cells from face connect
      tmp_face_cells.clear();
      result = mbImpl->get_adjacencies(connectf, num_connectf, 3, false, tmp_face_cells);
      if (MB_SUCCESS != result || tmp_face_cells.empty()) {
        mWriteIface->report_error("Error getting adj hexes.");
        return result;
      }

      bool is_internal = (tmp_face_cells.size() == 2);
      if (!is_polyh && is_internal) {
          // make sure 1st is forward sense
        int side_num, sense, offset;
        CN::SideNumber(etype, connectc, connectf, num_connectf,
                         2, side_num, sense, offset);
        if (sense == 1 && tmp_face_cells[0] != *rit) {
          assert(2 == tmp_face_cells.size());
          EntityHandle tmph = tmp_face_cells[0]; 
          tmp_face_cells[1] = tmp_face_cells[0]; 
          tmp_face_cells[0] = tmph;
        }
      }
        
        // get ids of cells in all_elems
      std::vector<int> *fcells_ptr, *fconnect_ptr;
      fcells_ptr = (is_internal ? &iface_cells : &eface_cells);
      fconnect_ptr = (is_internal ? &iface_connect : &eface_connect);
      fcells_ptr->push_back(egids[all_elems.index(tmp_face_cells[0])]);
      if (is_internal) fcells_ptr->push_back(egids[all_elems.index(tmp_face_cells[1])]);
      fconnect_ptr->push_back(num_connectf);

        // get indices of face vertices, add one
      for (int fv = 0; fv < num_connectf; fv++)
        fconnect_ptr->push_back(vgids[verts.index(connectf[fv])]);

      if (!is_polyh && is_internal) {
          // mark other cell for this face, if there is another cell
        EntityHandle other_cell = tmp_face_cells[0];
        const EntityHandle *connecto; int num_connecto;
        if (other_cell == *rit) other_cell = tmp_face_cells[1];
        result = mbImpl->get_connectivity(other_cell, connecto, num_connecto, 
                                          false, &storage);
        if (MB_SUCCESS != result) {
          mWriteIface->report_error("Couldn't get other entity connectivity.");
          return result;
        }
          // get side number
        int side_num, sense, offset;
        CN::SideNumber(TYPE_FROM_HANDLE(other_cell), connecto, connectf, num_connectf,
                         2, side_num, sense, offset);
          // set mark for this face
        short int tmp_mark, tmp_mark2;
        result = mbImpl->tag_get_data(mark_tag, &other_cell, 1, &tmp_mark);
        if (MB_SUCCESS != result) {
          mWriteIface->report_error("Couldn't get mark data for other cell.");
          return result;
        }
        tmp_mark2 = 0x1 << (unsigned int)side_num;
        assert("mark for this side on other entity shouldn't be set already" &&
               !(tmp_mark & tmp_mark2));
        tmp_mark |= tmp_mark2;
        result = mbImpl->tag_set_data(mark_tag, &other_cell, 1, &tmp_mark);
        if (MB_SUCCESS != result) {
          mWriteIface->report_error("Couldn't set mark data for other cell.");
          return result;
        }
      } // !is_polyh
    } // loop over faces in elem
  } // loop over elems

  int num_ifaces = iface_cells.size()/2,
      num_efaces = eface_cells.size();

    // write internal faces
  CCMIOID mapID;
  CCMIONewEntity(&error, rootID, kCCMIOMap, NULL, &mapID);
    // set gids for internal faces
  if ((int)all_elems.size() < (num_ifaces + num_efaces)) {
    delete [] egids;
    egids = new int[num_ifaces + num_efaces];
  }
  for (int i = 1; i <= (num_ifaces+num_efaces); i++) egids[i-1] = i;
  CCMIOWriteMap(&error, mapID, CCMIOSIZEC(num_ifaces),
                CCMIOSIZEC(num_ifaces+num_efaces),
                &egids[0],
                CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CCMIONewEntity(&error, topologyID, kCCMIOInternalFaces, "Internal faces", &id);
  CCMIOWriteFaces(&error, id, kCCMIOInternalFaces, mapID,
                  CCMIOSIZEC(iface_connect.size()), &iface_connect[0],
                  CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CCMIOWriteFaceCells(&error, id, kCCMIOInternalFaces, mapID, &iface_cells[0],
                      CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));

    // now external faces
  CCMIONewEntity(&error, rootID, kCCMIOMap, NULL, &mapID);
  CCMIOWriteMap(&error, mapID, CCMIOSIZEC(num_efaces),
                CCMIOSIZEC(num_efaces), &egids[num_ifaces],
                CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CCMIONewIndexedEntity(&error, topologyID, kCCMIOBoundaryFaces, 0,
                        "Boundary faces", &id);
  CCMIOWriteFaces(&error, id, kCCMIOBoundaryFaces, mapID,
                  CCMIOSIZEC(eface_connect.size()), &eface_connect[0],
                  CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
  CCMIOWriteFaceCells(&error, id, kCCMIOBoundaryFaces, mapID, &eface_cells[0],
                      CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));

  return MB_SUCCESS;
}

ErrorCode WriteCCMIO::open_file(const char* filename)
{
    // not a valid filname
  if(strlen((const char*)filename) == 0)
  {
    mWriteIface->report_error("Output filename not specified");
    return MB_FAILURE;
  }

    /* template - open file & store somewhere */

    // file couldn't be opened
  return MB_SUCCESS;
}

ErrorCode WriteCCMIO::get_neuset_elems(EntityHandle neuset, int current_sense,
                                         Range &forward_elems, Range &reverse_elems) 
{
  Range neuset_elems, neuset_meshsets;

    // get the sense tag; don't need to check return, might be an error if the tag
    // hasn't been created yet
  Tag sense_tag = 0;
  mbImpl->tag_get_handle("SENSE", sense_tag);

    // get the entities in this set
  ErrorCode result = mbImpl->get_entities_by_handle(neuset, neuset_elems, true);
  if (MB_FAILURE == result) return result;
  
    // now remove the meshsets into the neuset_meshsets; first find the first meshset,
  Range::iterator range_iter = neuset_elems.begin();
  while (TYPE_FROM_HANDLE(*range_iter) != MBENTITYSET && range_iter != neuset_elems.end())
    range_iter++;
  
    // then, if there are some, copy them into neuset_meshsets and erase from neuset_elems
  if (range_iter != neuset_elems.end()) {
    std::copy(range_iter, neuset_elems.end(), range_inserter(neuset_meshsets));
    neuset_elems.erase(range_iter, neuset_elems.end());
  }
  

    // ok, for the elements, check the sense of this set and copy into the right range
    // (if the sense is 0, copy into both ranges)

    // need to step forward on list until we reach the right dimension
  Range::iterator dum_it = neuset_elems.end();
  dum_it--;
  int target_dim = CN::Dimension(TYPE_FROM_HANDLE(*dum_it));
  dum_it = neuset_elems.begin();
  while (target_dim != CN::Dimension(TYPE_FROM_HANDLE(*dum_it)) &&
         dum_it != neuset_elems.end()) 
    dum_it++;

  if (current_sense == 1 || current_sense == 0)
    std::copy(dum_it, neuset_elems.end(), range_inserter(forward_elems));
  if (current_sense == -1 || current_sense == 0)
    std::copy(dum_it, neuset_elems.end(), range_inserter(reverse_elems));
  
    // now loop over the contained meshsets, getting the sense of those and calling this
    // function recursively
  for (range_iter = neuset_meshsets.begin(); range_iter != neuset_meshsets.end(); range_iter++) {

      // first get the sense; if it's not there, by convention it's forward
    int this_sense;
    if (0 == sense_tag ||
        MB_FAILURE == mbImpl->tag_get_data(sense_tag, &(*range_iter), 1, &this_sense))
      this_sense = 1;
      
      // now get all the entities on this meshset, with the proper (possibly reversed) sense
    get_neuset_elems(*range_iter, this_sense*current_sense,
                     forward_elems, reverse_elems);
  }
  
  return result;
}


} // namespace moab

  
