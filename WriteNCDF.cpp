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

#ifndef NETCDF_FILE
#  error Attempt to compile WriteNCDF with NetCDF support disabled
#endif

#include "WriteNCDF.hpp"

#include <utility>
#include <algorithm>
#include <time.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/MBCN.hpp"
#include "moab/MBTagConventions.hpp"
#include "Internals.hpp"
#include "ExoIIUtil.hpp"
#include "moab/WriteUtilIface.hpp"
#include "exodus_order.h"

namespace moab {

const int TIME_STR_LEN = 11;

#define INS_ID(stringvar, prefix, id) \
          sprintf(stringvar, prefix, id)

WriterIface* WriteNCDF::factory( Interface* iface )
  { return new WriteNCDF( iface ); }

WriteNCDF::WriteNCDF(Interface *impl) 
    : mdbImpl(impl), ncFile(0), mCurrentMeshHandle(0)
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
    result = impl->tag_create(MATERIAL_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER, mMaterialSetTag, &dum_val);
  
  result = impl->tag_get_handle(DIRICHLET_SET_TAG_NAME, mDirichletSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(DIRICHLET_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE,
                              MB_TYPE_INTEGER, mDirichletSetTag, &dum_val);
  
  result = impl->tag_get_handle(NEUMANN_SET_TAG_NAME,   mNeumannSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(NEUMANN_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE,
                              MB_TYPE_INTEGER, mNeumannSetTag, &dum_val);
  
  result = impl->tag_get_handle(HAS_MID_NODES_TAG_NAME, mHasMidNodesTag);
  if (MB_TAG_NOT_FOUND == result) {
    int dum_val_array[] = {0, 0, 0, 0};
    result = impl->tag_create(HAS_MID_NODES_TAG_NAME, 4*sizeof(int), MB_TAG_SPARSE,
                              MB_TYPE_INTEGER, mHasMidNodesTag, dum_val_array);
  }
  
  result = impl->tag_get_handle("distFactor",           mDistFactorTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create_variable_length( "distFactor", 
                                               MB_TAG_SPARSE, 
                                               MB_TYPE_DOUBLE,
                                               mDistFactorTag );
  
  result = impl->tag_get_handle("qaRecord",             mQaRecordTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create_variable_length( "qaRecord", 
                                               MB_TAG_SPARSE, 
                                               MB_TYPE_OPAQUE,
                                               mQaRecordTag );
  
  result = impl->tag_get_handle(GLOBAL_ID_TAG_NAME,             mGlobalIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER, mGlobalIdTag, &dum_val);
  

  impl->tag_create("WriteNCDF element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

}

WriteNCDF::~WriteNCDF() 
{
  std::string iface_name = "WriteUtilIface";
  mdbImpl->release_interface(iface_name, mWriteIface);

  mdbImpl->tag_delete(mEntityMark);

  if (NULL != ncFile)
    delete ncFile;
}

void WriteNCDF::reset_block(std::vector<MaterialSetData> &block_info)
{
  std::vector<MaterialSetData>::iterator iter;
  
  for (iter = block_info.begin(); iter != block_info.end(); iter++)
  {
    iter->elements.clear();
  }
}

void WriteNCDF::time_and_date(char* time_string, char* date_string)
{
   struct tm* local_time;
   time_t calendar_time;

   calendar_time = time(NULL);
   local_time = localtime(&calendar_time);

   assert(NULL != time_string && NULL != date_string);
   
   strftime(time_string, TIME_STR_LEN, "%H:%M:%S", local_time);
   strftime(date_string, TIME_STR_LEN, "%m/%d/%Y", local_time);

     // terminate with NULL character
   time_string[10] = (char)NULL;
   date_string[10] = (char)NULL;
}

ErrorCode WriteNCDF::write_file(const char *exodus_file_name, 
                                    const bool overwrite,
                                    const FileOptions&,
                                    const EntityHandle *ent_handles,
                                    const int num_sets,
                                    const std::vector<std::string> &qa_records,
                                    const Tag*,
                                    int,
                                    int user_dimension)
{
  assert(0 != mMaterialSetTag &&
         0 != mNeumannSetTag &&
         0 != mDirichletSetTag);

  if (user_dimension == 0)
   mdbImpl->get_dimension( user_dimension );
  

  std::vector<EntityHandle> blocks, nodesets, sidesets, entities;

    // separate into blocks, nodesets, sidesets

  if (num_sets == 0) {
      // default to all defined block, nodeset and sideset-type sets
    Range this_range;
    mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mMaterialSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(blocks));
    this_range.clear();
    mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mDirichletSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(nodesets));
    this_range.clear();
    mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mNeumannSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(sidesets));
  
    // If there is nothing to write, write everything as one block.
    if (blocks.empty() && nodesets.empty() && sidesets.empty())
    {
      this_range.clear();
      for (int d = user_dimension; d > 0 && this_range.empty(); --d)
        mdbImpl->get_entities_by_dimension( 0, d, this_range, false );
      if (this_range.empty())
        return MB_FILE_WRITE_ERROR;

      EntityHandle block_handle;
      int block_id = 1;
      mdbImpl->create_meshset( MESHSET_SET, block_handle );
      mdbImpl->tag_set_data( mMaterialSetTag, &block_handle, 1, &block_id );
      mdbImpl->add_entities( block_handle, this_range );
      blocks.push_back( block_handle );
    }
  }
  else {
    int dummy;
    for (const EntityHandle *iter = ent_handles; iter < ent_handles+num_sets; iter++) 
    {
      if (MB_SUCCESS == mdbImpl->tag_get_data(mMaterialSetTag, &(*iter), 1, &dummy))
        blocks.push_back(*iter);
      else if (MB_SUCCESS == mdbImpl->tag_get_data(mDirichletSetTag, &(*iter), 1, &dummy))
        nodesets.push_back(*iter);
      else if (MB_SUCCESS == mdbImpl->tag_get_data(mNeumannSetTag, &(*iter), 1, &dummy))
        sidesets.push_back(*iter);
    }
  }
  
  
    // if there is nothing to write just return.
  if (blocks.empty() && nodesets.empty() && sidesets.empty())
    return MB_FILE_WRITE_ERROR;

  // try to get mesh information
  ExodusMeshInfo mesh_info;

  std::vector<MaterialSetData> block_info;
  std::vector<NeumannSetData> sideset_info;
  std::vector<DirichletSetData> nodeset_info;

  mesh_info.num_dim = user_dimension;

  if (qa_records.empty()) {
      // qa records are empty - initialize some MB-standard ones
    mesh_info.qaRecords.push_back("MB");
    mesh_info.qaRecords.push_back("0.99");
    char string1[80], string2[80];
    time_and_date(string2, string1);
    mesh_info.qaRecords.push_back(string2);
    mesh_info.qaRecords.push_back(string1);
  }
  else {
        // constrained to multiples of 4 qa records
    assert(qa_records.size()%4 == 0);
    
    std::copy(qa_records.begin(), qa_records.end(), 
              std::back_inserter(mesh_info.qaRecords));
  }

  block_info.clear();
  if(gather_mesh_information(mesh_info, block_info, sideset_info, nodeset_info,
                             blocks, sidesets, nodesets) != MB_SUCCESS)
  {
    reset_block(block_info);
    return MB_FAILURE;
  }


  // try to open the file after gather mesh info succeeds
  ncFile = new NcFile(exodus_file_name, overwrite ? NcFile::Replace : NcFile::New );
  if (NULL == ncFile) {
    reset_block(block_info);
    return MB_FAILURE;
  }

  if( write_header(mesh_info, block_info, sideset_info,
                   nodeset_info, exodus_file_name) != MB_SUCCESS)
  {
    reset_block(block_info);
    return MB_FAILURE;
  }

  if( write_nodes(mesh_info.num_nodes, mesh_info.nodes, mesh_info.num_dim) != MB_SUCCESS )
  {
    reset_block(block_info);
    return MB_FAILURE;
  }

  if( write_elementblocks(block_info) )
  {
    reset_block(block_info);
    return MB_FAILURE;
  }

  // write the three maps
  if( write_global_node_order_map(mesh_info.num_nodes, mesh_info.nodes) != MB_SUCCESS )
  {
    reset_block(block_info);
    return MB_FAILURE;
  }

  if(write_global_element_order_map(mesh_info.num_elements) != MB_SUCCESS )
  {
    reset_block(block_info);
    return MB_FAILURE;
  }

  if(write_element_order_map(mesh_info.num_elements) != MB_SUCCESS )
  {
    reset_block(block_info);
    return MB_FAILURE;
  }

/*
 if(write_elementmap(mesh_info) != MB_SUCCESS)
   return MB_FAILURE;
*/

  if(write_BCs(sideset_info, nodeset_info) != MB_SUCCESS)
  {
    reset_block(block_info);
    return MB_FAILURE;
  }

  if( write_qa_records( mesh_info.qaRecords) != MB_SUCCESS )
    return MB_FAILURE;

  // copy the qa records into the argument
  //mesh_info.qaRecords.swap(qa_records);
  
  return MB_SUCCESS;
}

ErrorCode WriteNCDF::gather_mesh_information(
                        ExodusMeshInfo &mesh_info,
                        std::vector<MaterialSetData> &block_info,
                        std::vector<NeumannSetData> &sideset_info,
                        std::vector<DirichletSetData> &nodeset_info,
                        std::vector<EntityHandle> &blocks,
                        std::vector<EntityHandle> &sidesets,
                        std::vector<EntityHandle> &nodesets)
{
  ErrorCode rval;
  std::vector<EntityHandle>::iterator vector_iter, end_vector_iter;

  mesh_info.num_nodes = 0;
  mesh_info.num_elements = 0;
  mesh_info.num_elementblocks = 0;
  
  int id = 0;

  vector_iter= blocks.begin();
  end_vector_iter = blocks.end();

  std::vector<EntityHandle> parent_meshsets;

  // clean out the bits for the element mark
  rval = mdbImpl->tag_delete(mEntityMark);
  if (MB_SUCCESS != rval)
    return rval;
  rval = mdbImpl->tag_create("WriteNCDF element mark", 1, MB_TAG_BIT, mEntityMark, NULL);
  if (MB_SUCCESS != rval)
    return rval;

  int highest_dimension_of_element_blocks = 0;

  for(vector_iter = blocks.begin(); vector_iter != blocks.end(); vector_iter++)
  {
       
    MaterialSetData block_data;

    //for the purpose of qa records, get the parents of these blocks 
    if( mdbImpl->get_parent_meshsets( *vector_iter, parent_meshsets ) != MB_SUCCESS )
      return MB_FAILURE;

    // get all Entity Handles in the mesh set
    Range dummy_range;
    rval = mdbImpl->get_entities_by_handle(*vector_iter, dummy_range, true );
    if (MB_SUCCESS != rval)
      return rval;

    // skip empty blocks
    if (dummy_range.empty())
      continue;

    // wait a minute, we are doing some filtering here that doesn't make sense at this level  CJS

      // find the dimension of the last entity in this range
    Range::iterator entity_iter = dummy_range.end();
    entity_iter = dummy_range.end();
    entity_iter--;
    int this_dim = MBCN::Dimension(TYPE_FROM_HANDLE(*entity_iter));
    entity_iter = dummy_range.begin();
    while (entity_iter != dummy_range.end() &&
           MBCN::Dimension(TYPE_FROM_HANDLE(*entity_iter)) != this_dim)
      entity_iter++;
    
    if (entity_iter != dummy_range.end())
      std::copy(entity_iter, dummy_range.end(), range_inserter(block_data.elements));

    assert(block_data.elements.begin() == block_data.elements.end() ||
           MBCN::Dimension(TYPE_FROM_HANDLE(*(block_data.elements.begin()))) == this_dim);
    
    // end of -- wait a minute, we are doing some filtering here that doesn't make sense at this level CJS
   

    // get the block's id
    if(mdbImpl->tag_get_data(mMaterialSetTag, &(*vector_iter), 1, &id) != MB_SUCCESS ) {
      mWriteIface->report_error("Couldn't get block id from a tag for an element block.");
      return MB_FAILURE;
    }
    
    block_data.id = id; 
    block_data.number_attributes = 0;
 
     // iterate through all the elements in the meshset
    Range::iterator elem_range_iter, end_elem_range_iter;
    elem_range_iter = block_data.elements.begin();
    end_elem_range_iter = block_data.elements.end();

      // get the entity type for this block, verifying that it's the same for all elements
      // THIS ASSUMES HANDLES SORT BY TYPE!!!
    EntityType entity_type = TYPE_FROM_HANDLE(*elem_range_iter);
    end_elem_range_iter--;
    if (entity_type != TYPE_FROM_HANDLE(*(end_elem_range_iter++))) {
      mWriteIface->report_error("Entities in block %i not of common type", id);
      return MB_FAILURE;
    }

    int dimension = -1;
    if(entity_type == MBQUAD || entity_type == MBTRI)
      dimension = 3;   // output shells by default
    else if(entity_type == MBEDGE)
      dimension = 2;
    else
      dimension = MBCN::Dimension(entity_type);

    if( dimension > highest_dimension_of_element_blocks )
      highest_dimension_of_element_blocks = dimension;

    std::vector<EntityHandle> tmp_conn;
    rval = mdbImpl->get_connectivity(&(*(block_data.elements.begin())), 1, tmp_conn);
    if (MB_SUCCESS != rval)
      return rval;
    block_data.element_type = ExoIIUtil::get_element_type_from_num_verts(tmp_conn.size(), entity_type, dimension);
    
    if (block_data.element_type == EXOII_MAX_ELEM_TYPE) {
      mWriteIface->report_error("Element type in block %i didn't get set correctly", id);
      return MB_FAILURE;
    }
    
    block_data.number_nodes_per_element = ExoIIUtil::VerticesPerElement[block_data.element_type];

    // number of nodes for this block
    block_data.number_elements = block_data.elements.size();

    // total number of elements
    mesh_info.num_elements += block_data.number_elements;

    // get the nodes for the elements
    rval = mWriteIface->gather_nodes_from_elements(block_data.elements, mEntityMark, mesh_info.nodes);
    if (MB_SUCCESS != rval)
      return rval;

    if(!sidesets.empty())
    {
      // if there are sidesets, keep track of which elements are being written out
      for(Range::iterator iter = block_data.elements.begin(); 
          iter != block_data.elements.end(); ++iter)
      {
        unsigned char bit = 0x1;
        rval = mdbImpl->tag_set_data(mEntityMark, &(*iter), 1, &bit);
        if (MB_SUCCESS != rval)
          return rval;
      }
    }

    block_info.push_back( block_data );
  
    std::vector<char*> *qa_rec = NULL;

    if (mdbImpl->tag_get_data(mQaRecordTag, &(*vector_iter), 1, &qa_rec) == MB_SUCCESS &&
        NULL != qa_rec) 
    {
        // there are qa records on this block - copy over to mesh qa records

        // constrained to multiples of 4 qa records
      assert(qa_rec->size()%4 == 0);
      
      for(unsigned int k=0; k<qa_rec->size(); k++)
        mesh_info.qaRecords.push_back( (*qa_rec)[k] );
    }
  }
 

  mesh_info.num_elementblocks = block_info.size();

  //if user hasn't entered dimension, we figure it out
  if( mesh_info.num_dim == 0 )
  {
    //never want 1 or zero dimensions
    if( highest_dimension_of_element_blocks < 2 )
      mesh_info.num_dim = 3;
    else
      mesh_info.num_dim = highest_dimension_of_element_blocks;
  }


  Range::iterator range_iter, end_range_iter;
  range_iter = mesh_info.nodes.begin();
  end_range_iter = mesh_info.nodes.end();

  mesh_info.num_nodes = mesh_info.nodes.size(); 

  //------nodesets--------
  
  vector_iter= nodesets.begin();
  end_vector_iter = nodesets.end();

  for(; vector_iter != end_vector_iter; vector_iter++)
  {
    
    DirichletSetData nodeset_data;
    nodeset_data.id = 0;
    nodeset_data.number_nodes = 0;

    // get the nodeset's id
    if(mdbImpl->tag_get_data(mDirichletSetTag,&(*vector_iter), 1,&id) != MB_SUCCESS) {
      mWriteIface->report_error("Couldn't get id tag for nodeset %i", id);
      return MB_FAILURE;
    }
    
    nodeset_data.id = id; 

    std::vector<EntityHandle> node_vector;
    //get the nodes of the nodeset that are in mesh_info.nodes
    if( mdbImpl->get_entities_by_handle(*vector_iter, node_vector, true) != MB_SUCCESS ) {
      mWriteIface->report_error("Couldn't get nodes in nodeset %i", id);
      return MB_FAILURE;
    }



    //get the tag for distribution factors
    const double *dist_factor_vector;
    int dist_factor_size;
    const void* ptr = 0;
   
    int has_dist_factors = 0; 
    if(mdbImpl->tag_get_data(mDistFactorTag,&(*vector_iter), 1, &ptr, &dist_factor_size) == MB_SUCCESS &&
       dist_factor_size)
      has_dist_factors = 1;
    dist_factor_size /= sizeof(double);
    dist_factor_vector = reinterpret_cast<const double*>(ptr);
    std::vector<EntityHandle>::iterator iter, end_iter;
    iter = node_vector.begin();
    end_iter= node_vector.end();
 
    int j=0; 
    unsigned char node_marked = 0;
    ErrorCode result;
    for(; iter != end_iter; iter++)
    {
      if (TYPE_FROM_HANDLE(*iter) != MBVERTEX) continue;
      result = mdbImpl->tag_get_data(mEntityMark, &(*iter), 1, &node_marked);
      if (MB_SUCCESS != result) {
        mWriteIface->report_error("Couldn't get mark data.");
        return result;
      }
      
      if(node_marked == 0x1)
      {
        nodeset_data.nodes.push_back( *iter );    
        if( has_dist_factors != 0)
          nodeset_data.node_dist_factors.push_back( dist_factor_vector[j] );
        else
          nodeset_data.node_dist_factors.push_back( 1.0 );
      }
      j++;
    } 
    
    nodeset_data.number_nodes = nodeset_data.nodes.size(); 
    nodeset_info.push_back( nodeset_data );

  }

  //------sidesets--------
  vector_iter= sidesets.begin();
  end_vector_iter = sidesets.end();

  for(; vector_iter != end_vector_iter; vector_iter++)
  {

    NeumannSetData sideset_data;

    // get the sideset's id
    if(mdbImpl->tag_get_data(mNeumannSetTag,&(*vector_iter), 1,&id) != MB_SUCCESS)
      return MB_FAILURE;

    sideset_data.id = id; 
    sideset_data.mesh_set_handle = *vector_iter; 
 
    //get the sides in two lists, one forward the other reverse; starts with forward sense
      // by convention
    Range forward_elems, reverse_elems;
    if(get_sideset_elems(*vector_iter, 0, forward_elems, reverse_elems) == MB_FAILURE)
      return MB_FAILURE;

    ErrorCode result = get_valid_sides(forward_elems, mesh_info, 1, sideset_data);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get valid sides data.");
      return result;
    }
    result = get_valid_sides(reverse_elems, mesh_info, -1, sideset_data);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get valid sides data.");
      return result;
    }
    
    sideset_data.number_elements = sideset_data.elements.size(); 
    sideset_info.push_back( sideset_data );

  }

  return MB_SUCCESS;

}

ErrorCode WriteNCDF::get_valid_sides(Range &elems, ExodusMeshInfo& /*mesh_info*/, 
                                         const int sense,
                                         NeumannSetData &sideset_data) 
{
    // this is where we see if underlying element of side set element is included in output 

    //get the sideset-based info for distribution factors
  const double *dist_factor_vector = 0;
  int dist_factor_size = 0;

  // initialize dist_fac_iter to get rid of compiler warning
  const double* dist_fac_iter = 0;
  const void* ptr = 0;
  bool has_dist_factors = false; 
  if(mdbImpl->tag_get_data(mDistFactorTag,
                           &(sideset_data.mesh_set_handle), 1, &ptr, &dist_factor_size) == MB_SUCCESS &&
     dist_factor_size)
  {
    has_dist_factors = true;
    dist_factor_vector = reinterpret_cast<const double*>(ptr);
    dist_fac_iter = dist_factor_vector;
    dist_factor_size /= sizeof(double);
  }

  unsigned char element_marked = 0;
  ErrorCode result;
  for(Range::iterator iter = elems.begin(); iter != elems.end(); iter++)
  {
      // should insert here if "side" is a quad/tri on a quad/tri mesh
    result = mdbImpl->tag_get_data(mEntityMark, &(*iter), 1, &element_marked);
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get mark data.");
      return result;
    }
    
    if(element_marked == 0x1)
    {
      sideset_data.elements.push_back( *iter );

        // TJT TODO: the sense should really be # edges + 1or2
      sideset_data.side_numbers.push_back((sense == 1 ? 1 : 2));
    }
    else //then "side" is probably a quad/tri on a hex/tet mesh
    {
      std::vector<EntityHandle> parents;
      int dimension = MBCN::Dimension( TYPE_FROM_HANDLE(*iter));

        //get the adjacent parent element of "side"
      if( mdbImpl->get_adjacencies( &(*iter), 1, dimension+1, false, parents) != MB_SUCCESS ) {
        mWriteIface->report_error("Couldn't get adjacencies for sideset.");
        return MB_FAILURE;
      }
       
      if(!parents.empty())     
      {
          //make sure the adjacent parent element will be output
        for(unsigned int k=0; k<parents.size(); k++)
        {
          result = mdbImpl->tag_get_data(mEntityMark, &(parents[k]), 1, &element_marked);
          if (MB_SUCCESS != result) {
            mWriteIface->report_error("Couldn't get mark data.");
            return result;
          }
        
          int side_no, this_sense, this_offset;
          if(element_marked == 0x1 &&
             mdbImpl->side_number(parents[k], *iter, side_no, 
                                  this_sense, this_offset) == MB_SUCCESS &&
             this_sense == sense) {
            sideset_data.elements.push_back(parents[k]);
            sideset_data.side_numbers.push_back(side_no+1);
            break;
          }
        }
      }
      else
      {
        mWriteIface->report_error("No parent element exists for element in sideset %i", sideset_data.id);
        return MB_FAILURE;
      }
    }

    if( sideset_data.elements.size() != 0 )
    {
        // distribution factors
      int num_nodes = MBCN::VerticesPerEntity(TYPE_FROM_HANDLE(*iter));
      if( has_dist_factors )
      {
        std::copy(dist_fac_iter, dist_fac_iter + num_nodes, 
                  std::back_inserter(sideset_data.ss_dist_factors) );
        dist_fac_iter += num_nodes;
      }
      else
      {
        for(int j=0; j<num_nodes; j++)
          sideset_data.ss_dist_factors.push_back( 1.0 );
      }
    }
  }

  return MB_SUCCESS;
}

ErrorCode WriteNCDF::write_qa_records(std::vector<std::string> &qa_record_list)
{
  int i = 0;
  
  for(std::vector<std::string>::iterator string_it = qa_record_list.begin();
      string_it != qa_record_list.end(); )
  {
    for (int j = 0; j < 4; j++)
    {
      write_qa_string((*string_it++).c_str(), i, j);
    }
    i++;
  }

  return MB_SUCCESS;
}

ErrorCode WriteNCDF::write_qa_string(const char *string,
                                         int record_number,
                                         int record_position)
{
      // get the variable id in the exodus file

  NcVar *qavar = ncFile->get_var("qa_records");
  if (NULL == qavar) {
    mWriteIface->report_error("Failed to find qa record array in file");
    return MB_FAILURE;
  }

//  write out the record

  long start[3], count[3];

  start[0] = record_number;
  start[1] = record_position;
  start[2] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = (long)strlen(string) +1;
  NcBool status = qavar->set_cur(record_number, record_position, 0);
  if (!status) {
    mWriteIface->report_error("Failed to position qa string variable.");
    return MB_FAILURE;
  }
  status = qavar->put(string, 1, 1, strlen(string)+1);
  if (!status) {
    mWriteIface->report_error("Failed to write qa string variable.");
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}


ErrorCode WriteNCDF::write_nodes(int num_nodes, Range& nodes, int dimension)
{
  // write coordinates names
  NcVar *coord_var = ncFile->get_var("coor_names");
  if (NULL == coord_var) return MB_FAILURE;
  NcBool status = coord_var->put("x", 1, 1);
  if (!status) return MB_FAILURE;
  status = coord_var->set_cur(1, 0);
  if (!status) return MB_FAILURE;
  status = coord_var->put("y", 1, 1);
  if (!status) return MB_FAILURE;
  status = coord_var->set_cur(2, 0);
  if (!status) return MB_FAILURE;
  status = coord_var->put("z", 1, 1);
  if (!status) return MB_FAILURE;
  
  //see if should transform coordinates
  ErrorCode result;
  Tag trans_tag;
  result = mdbImpl->tag_get_handle( MESH_TRANSFORM_TAG_NAME, trans_tag);
  bool transform_needed = true;
  if( result == MB_TAG_NOT_FOUND )
    transform_needed = false;

  int num_coords_to_fill = transform_needed ? 3 : dimension;

  std::vector<double*> coord_arrays(3);
  coord_arrays[0] = new double[num_nodes];
  coord_arrays[1] = new double[num_nodes];
  coord_arrays[2] = NULL;

  if( num_coords_to_fill == 3 ) 
    coord_arrays[2] = new double[num_nodes];
 
  result = mWriteIface->get_node_arrays(dimension, num_nodes, nodes, mGlobalIdTag, 1, coord_arrays);
  if(result != MB_SUCCESS)
  {
    delete [] coord_arrays[0];
    delete [] coord_arrays[1];
    if(coord_arrays[2]) delete [] coord_arrays[2];
    return result;
  }

  if( transform_needed )
  {
    double trans_matrix[16]; 
    result = mdbImpl->tag_get_data( trans_tag, NULL, 0, trans_matrix ); 
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get transform data.");
      return result;
    }
      
    for( int i=0; i<num_nodes; i++)
    {

      double vec1[3];
      double vec2[3];

      vec2[0] =  coord_arrays[0][i];
      vec2[1] =  coord_arrays[1][i];
      vec2[2] =  coord_arrays[2][i];

      for( int row=0; row<3; row++ )
      {
        vec1[row] = 0.0;
        for( int col = 0; col<3; col++ )
        {
          vec1[row] += ( trans_matrix[ (row*4)+col ] * vec2[col] );
        }
      }

      coord_arrays[0][i] = vec1[0];
      coord_arrays[1][i] = vec1[1];
      coord_arrays[2][i] = vec1[2];

    }
  }


  // write the nodes 
  NcVar *coord = ncFile->get_var("coord");
  if (NULL == coord) return MB_FAILURE;
  status = coord->put(coord_arrays[0], 1, num_nodes);
  if (!status) return MB_FAILURE;
  status = coord->set_cur(1, 0);
  if (!status) return MB_FAILURE;
  status = coord->put(&(coord_arrays[1][0]), 1, num_nodes);
  if (!status) return MB_FAILURE;
  status = coord->set_cur(2, 0);
  if (!status) return MB_FAILURE;
  status = coord->put(&(coord_arrays[2][0]), 1, num_nodes);
  if (!status) return MB_FAILURE;
  
  delete [] coord_arrays[0];
  delete [] coord_arrays[1];
  if(coord_arrays[2]) delete [] coord_arrays[2];

  return MB_SUCCESS;

}


ErrorCode WriteNCDF::write_header(ExodusMeshInfo& mesh_info,
                                     std::vector<MaterialSetData> &block_info,
                                     std::vector<NeumannSetData> &sideset_info,
                                     std::vector<DirichletSetData> &nodeset_info,
                                     const char *filename)
{

  // get the date and time
  char time[TIME_STR_LEN];
  char date[TIME_STR_LEN];
  time_and_date(time,date);

  std::string title_string = "MOAB"; 
  title_string.append( "(" );
  title_string.append( filename );
  title_string.append( "): ");
  title_string.append( date );
  title_string.append( ": " );
  title_string.append( "time " );

  if(title_string.length() > ExoIIInterface::MAX_LINE_LENGTH)
    title_string.resize( ExoIIInterface::MAX_LINE_LENGTH );

  // initialize the exodus file

  int result = initialize_exodus_file(mesh_info,
                                      block_info, sideset_info,
                                      nodeset_info, title_string.c_str());

  if(result == MB_FAILURE)
    return MB_FAILURE;

  return MB_SUCCESS;
}

ErrorCode WriteNCDF::write_elementblocks(std::vector<MaterialSetData> &block_data )
{

  unsigned int i;
  int block_index = 0;  // index into block list, may != 1 if there are inactive blocks
  int exodus_id = 1;

  for(i=0; i< block_data.size(); i++)
  {
    MaterialSetData& block = block_data[i];

    unsigned int num_nodes_per_elem = block.number_nodes_per_element; 
                                     
    // write out the id for the block

    int id = block.id;
    int num_values = 1;

    if( write_exodus_integer_variable("eb_prop1", &id, block_index, num_values) != MB_SUCCESS )
    {
      mWriteIface->report_error("Problem writing element block id %i", id);
    }

    // write out the block status

    int status = 1;
    if (!block.number_elements)
    {
      mWriteIface->report_error("Warning: No elements in block %i", id); 
    }

    if( write_exodus_integer_variable("eb_status", &status, block_index, num_values) != MB_SUCCESS )
    {
      mWriteIface->report_error("Problem writing element block status");
    }

    // map the connectivity to the new nodes
    const unsigned int num_elem = block.number_elements;
    const unsigned int num_nodes = num_nodes_per_elem * num_elem;
    int* connectivity = new int[num_nodes];

    ErrorCode result = mWriteIface->get_element_array(
        num_elem, num_nodes_per_elem, mGlobalIdTag, block.elements, mGlobalIdTag, exodus_id ,connectivity);

    if(result != MB_SUCCESS) {
      mWriteIface->report_error("Couldn't get element array to write from.");
      delete [] connectivity;
      return result;
    }
    
    // if necessary, convert from EXODUS to MBCN node order
    const EntityType elem_type = ExoIIUtil::ExoIIElementMBEntity[block.element_type];
    assert( block.elements.all_of_type( elem_type ) );
    const int* reorder = exodus_elem_order_map[elem_type][block.number_nodes_per_element];
    if (reorder)
      WriteUtilIface::reorder( reorder, connectivity, 
                                 block.number_elements,
                                 block.number_nodes_per_element );

    exodus_id += num_elem;

    char wname[80];
    INS_ID(wname, "connect%d", i+1);
    NcVar *conn = ncFile->get_var(wname);
    if (NULL == conn) {
      delete [] connectivity;
      return MB_FAILURE;
    }
    NcBool err = conn->put(connectivity, num_elem, num_nodes_per_elem);
    delete [] connectivity;
    if(!err)
      return MB_FAILURE;

    block_index++;

  }

  return MB_SUCCESS;
}


ErrorCode WriteNCDF::write_global_node_order_map(int num_nodes, Range& nodes)
{
  // note: this routine bypasses the standard exodusII interface for efficiency!

  // node order map
  int* map = new int[num_nodes];

  // for now, output a dummy map!

   Range::iterator range_iter, end_iter;
   range_iter = nodes.begin();
   end_iter = nodes.end();

   int i=0;

   for(; range_iter != end_iter; range_iter++)
   {
      // TODO -- do we really want to cast this to an int?
      map[i++] = (int)ID_FROM_HANDLE(*range_iter);
   }

   // output array and cleanup

   int error = write_exodus_integer_variable("node_num_map",
                                             map,
                                             0,
                                             num_nodes);


  if(map)
    delete [] map;

  if( error < 0 )
  {
    mWriteIface->report_error("Failed writing global node order map");
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode WriteNCDF::write_global_element_order_map(int num_elements)
{

  // allocate map array
  int* map = new int[num_elements];


  // Many Sandia codes assume this map is unique, and CUBIT does not currently
  // have unique ids for all elements.  Therefore, to make sure nothing crashes,
  // insert a dummy map...

   for(int i=0; i<num_elements; i++)
   {
      map[i] = i+1;
   }


   // output array and cleanup

   int error = write_exodus_integer_variable("elem_num_map",
                                             map,
                                             0,
                                             num_elements);


  if(map)
    delete [] map;

  if( error < 0 )
  {
    mWriteIface->report_error("Failed writing global element order map");
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}


ErrorCode WriteNCDF::write_element_order_map(int num_elements)
{
  // note: this routine bypasses the standard exodusII interface for efficiency!

  // element order map
  int* map = new int[num_elements];

  // for now, output a dummy map!

   for(int i=0; i<num_elements; i++)
   {
      map[i] = i+1;
   }

   // output array and cleanup

   int error = write_exodus_integer_variable("elem_map",
                                             map,
                                             0,
                                             num_elements);


  if(map)
    delete [] map;

  if( error < 0 )
  {
    mWriteIface->report_error("Failed writing element map");
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}




ErrorCode WriteNCDF::write_exodus_integer_variable(const char* variable_name,
                                                      int *variable_array,
                                                      int start_position,
                                                      int number_values)
{

  // note: this routine bypasses the standard exodusII interface for efficiency!

   // write directly to netcdf interface for efficiency

   // get the variable id of the element map

  NcVar *ivar = ncFile->get_var(variable_name);
  if (NULL == ivar)
  {
    mWriteIface->report_error("WriteNCDF: failed to locate variable %s in file.", variable_name);
    return MB_FAILURE;
  }
   // this contortion is necessary because netCDF is expecting nclongs;
   // fortunately it's necessary only when ints and nclongs aren't the same size

   long start[1], count[1];
   NcBool iresult = true;
   start[0] = start_position;
   count[0] = number_values;

   if (0 != start_position) iresult = ivar->set_cur(start_position);
   
   if (sizeof(int) == sizeof(ncLong)) {
     iresult = ivar->put(variable_array, number_values);
   } else {
     long *lptr = new long[number_values];
     for (int jj = 0; jj < number_values; jj++)
       lptr[jj] = variable_array[jj];
     iresult = ivar->put(lptr, number_values);
     delete [] lptr;
   }

   //int error = 0;

   if (!iresult)
   {
     mWriteIface->report_error("Failed to store variable %s", variable_name);
     return MB_FAILURE;
   }

   return MB_SUCCESS;
}



ErrorCode WriteNCDF::write_BCs(std::vector<NeumannSetData> &sidesets,
                                       std::vector<DirichletSetData> &nodesets)
{
  unsigned int i,j;
  int id; 
  int ns_index = -1;

  for(std::vector<DirichletSetData>::iterator ns_it = nodesets.begin();
      ns_it != nodesets.end(); ns_it++)
  { 

    //get number of nodes in set
    int number_nodes = (*ns_it).number_nodes;
    if (0 == number_nodes) continue;

    // if we're here, we have a non-empty nodeset; increment the index
    ns_index++;

    //get the node set id 
    id = (*ns_it).id;

    //build new array to old exodus ids
    int * exodus_id_array = new int[number_nodes];
    double * dist_factor_array = new double[number_nodes];

    std::vector<EntityHandle>::iterator begin_iter, end_iter;
    std::vector<double>::iterator other_iter;
    begin_iter = (*ns_it).nodes.begin();
    end_iter = (*ns_it).nodes.end();
    other_iter = (*ns_it).node_dist_factors.begin();

    j=0;
    int exodus_id;
    ErrorCode result;
    //fill up node array and dist. factor array at the same time
    for(; begin_iter != end_iter; begin_iter++)
    {
      result = mdbImpl->tag_get_data(mGlobalIdTag, &(*begin_iter), 1, &exodus_id);
      if (MB_SUCCESS != result) {
        mWriteIface->report_error("Problem getting id tag data.");
        return result;
      }
      
      exodus_id_array[j] = exodus_id;
      dist_factor_array[j] = *(other_iter);
      other_iter++;
      j++;
    }

    // write out the id for the nodeset

    int num_values = 1;
    
    result = write_exodus_integer_variable("ns_prop1",
                                           &id, ns_index, num_values);

    if (result != MB_SUCCESS)
    {
      mWriteIface->report_error("Problem writing node set id %d", id);
      return MB_FAILURE; 
    }
    
      // write out the nodeset status

    int status = 1;
    if (!number_nodes)
      status = 0;

    result = write_exodus_integer_variable("ns_status",
                                           &status, ns_index, num_values);

    if (result != MB_SUCCESS)
    {
      mWriteIface->report_error("Problem writing node set status");
      return MB_FAILURE; 
    }

      //write it out
    char wname[80];
    INS_ID(wname, "node_ns%d", ns_index+1);
    NcVar *nsvar = ncFile->get_var(wname);
    if (NULL == nsvar) return MB_FAILURE;
    
    NcBool error = nsvar->put(exodus_id_array, number_nodes);
    if(!error)
    {
      mWriteIface->report_error("Failed writing exodus id array");
      return MB_FAILURE;
    }

      // write out nodeset distribution factors
    INS_ID(wname, "dist_fact_ns%d", ns_index+1);
    nsvar = ncFile->get_var(wname);
    if (NULL == nsvar) return MB_FAILURE;
    
    error = nsvar->put(dist_factor_array, number_nodes);
    if(!error)
    {
      mWriteIface->report_error("Failed writing nodeset distribution factors");
      return MB_FAILURE;
    }

    delete [] dist_factor_array;
    delete [] exodus_id_array;
  }


    //now do sidesets 
  int ss_index = 0; // index of sideset - not the same as 'i' because 
                    // only writing non-empty side sets
  for( i=0; i<sidesets.size(); i++)
  { 

    NeumannSetData sideset_data = sidesets[i]; 

      //get the side set id 
    int side_set_id = sideset_data.id;

      //get number of elements in set
    int number_elements = sideset_data.number_elements;
    if( number_elements == 0 )
      continue;

      //build new array to old exodus ids
    int * output_element_ids = new int[number_elements];
    int * output_element_side_numbers = new int[number_elements]; 

    std::vector<EntityHandle>::iterator begin_iter, end_iter;
    begin_iter = sideset_data.elements.begin();
    end_iter = sideset_data.elements.end();
    std::vector<int>::iterator side_iter = sideset_data.side_numbers.begin();

      //get the tag handle
    j=0;
    int exodus_id;

      //for each "side"
    for(; begin_iter != end_iter; begin_iter++, side_iter++)
    {
      ErrorCode result = mdbImpl->tag_get_data(mGlobalIdTag,
                                                  &(*begin_iter), 
                                                  1, &exodus_id);
      if (MB_FAILURE == result) {
        mWriteIface->report_error("Problem getting exodus id for sideset element %lu", 
                                  (long unsigned int) ID_FROM_HANDLE(*begin_iter));
        return result;
      }
      
      output_element_ids[j] = exodus_id;
      output_element_side_numbers[j++] = *side_iter;
    }

    if(number_elements)
    {
        // write out the id for the nodeset

      int num_values = 1;
    
        // ss_prop1[ss_index] = side_set_id
      ErrorCode result = write_exodus_integer_variable("ss_prop1",
                                                          &side_set_id, 
                                                          ss_index,num_values);

      if (result != MB_SUCCESS)
      {
        mWriteIface->report_error("Problem writing node set id %d", id);
        return MB_FAILURE; 
      }

        // FIXME : Something seems wrong here.  The we are within a block
        // started with if(number_elements), so this condition is always
        // false.  This code seems to indicate that we want to write all
        // sidesets, not just empty ones.  But the code both here and in
        // initialize_exodus_file() skip empty side sets.
        //  - j.k. 2007-03-09
      int status = 1;
      if (!number_elements)
        status = 0;

        // ss_status[ss_index] = status
      result = write_exodus_integer_variable("ss_status",
                                             &status, 
                                             ss_index, num_values);
      if (result != MB_SUCCESS)
      {
        mWriteIface->report_error("Problem writing side set status");
        return MB_FAILURE; 
      }
      
        // Increment ss_index now because we want a) we need to 
        // increment it somewhere within the if(number_elements)
        // block and b) the above calls need a zero-based index
        // while the following use a one-based index.
      ++ss_index;

      char wname[80];
      INS_ID(wname, "elem_ss%d", ss_index);
      NcVar *ssvar = ncFile->get_var(wname);
      if (NULL == ssvar) return MB_FAILURE;
      NcBool error = ssvar->put(output_element_ids, number_elements);
      if(!error)
      {
        mWriteIface->report_error("Failed writing sideset element array");
        return MB_FAILURE;
      }

      INS_ID(wname, "side_ss%d", ss_index);
      ssvar = ncFile->get_var(wname);
      if (NULL == ssvar) return MB_FAILURE;
      error = ssvar->put(output_element_side_numbers, number_elements);
      if(!error)
      {
        mWriteIface->report_error("Failed writing sideset side array");
        return MB_FAILURE;
      }

      INS_ID(wname, "dist_fact_ss%d", ss_index);
      ssvar = ncFile->get_var(wname);
      if (NULL == ssvar) return MB_FAILURE;
      error = ssvar->put(&(sideset_data.ss_dist_factors[0]), sideset_data.ss_dist_factors.size());
      if(!error)
      {
        mWriteIface->report_error("Failed writing sideset dist factors array");
        return MB_FAILURE;
      }
    } 
    delete [] output_element_ids;
    delete [] output_element_side_numbers;
    
  }

  return MB_SUCCESS;
}

ErrorCode WriteNCDF::initialize_exodus_file(ExodusMeshInfo &mesh_info,
                                                std::vector<MaterialSetData> &block_data,
                                                std::vector<NeumannSetData> & sideset_data,
                                                std::vector<DirichletSetData> & nodeset_data,
                                                const char* title_string,
                                                bool write_maps,
                                                bool /* write_sideset_distribution_factors */)
{
    // This routine takes the place of the exodusii routine ex_put_init,
    // and additionally pre-defines variables such as qa, element blocks,
    // sidesets and nodesets in a single pass.
    //
    // This is necessary because of the way exodusII works.  Each time the
    // netcdf routine endef is called to take the file out of define mode,
    // the entire file is copied before the new information is added.
    //
    // With very large files, this is simply not workable.  This routine takes
    // the definition portions of all applicable exodus routines and puts them
    // in a single definition, preventing repeated copying of the file.
    //
    // Most of the code is copied directly from the applicable exodus routine,
    // and thus this routine may not seem very consistent in usage or variable
    // naming, etc.

    // perform the initializations

  int element_block_index;

    // inquire on defined string dimension and general dimension for qa

  NcDim *dim_str, *dim_four, *dim_line, *dim_time;
  if ((dim_str = ncFile->add_dim("len_string", ExoIIInterface::MAX_STR_LENGTH)) == NULL)
  {
    mWriteIface->report_error("WriteNCDF: failed to get string length in file");
    return (MB_FAILURE);
  }
  if ((dim_line = ncFile->add_dim("len_line", ExoIIInterface::MAX_LINE_LENGTH)) == NULL)
  {
    mWriteIface->report_error("WriteNCDF: failed to get line length in file");
    return (MB_FAILURE);
  }
  if ((dim_four = ncFile->add_dim("four", 4)) == NULL)
  {
    mWriteIface->report_error("WriteNCDF: failed to locate four in file");
    return (MB_FAILURE);
  }
  if ((dim_time = ncFile->add_dim("time_step")) == NULL)
  {
    mWriteIface->report_error("WriteNCDF: failed to locate time step in file");
    return (MB_FAILURE);
  }

  

/* put file into define mode */

    // it is possible that an NT filename using backslashes is in the title string
    // this can give fits to unix codes where the backslash is assumed to be an escape
    // sequence.  For the exodus file, just flip the backslash to a slash to prevent
    // this problem

    // get a working copy of the title_string;

  char working_title[80];
  strncpy(working_title,title_string,79);

  int length = strlen(working_title);
  for(int pos = 0; pos < length; pos++)
  {
    if (working_title[pos] == '\\')
      strncpy(&working_title[pos],"/",1);
  }

  if (!ncFile->add_att("title", working_title))
  {
    mWriteIface->report_error("WriteNCDF: failed to define title attribute");
    return (MB_FAILURE);
  }

    // add other attributes while we're at it
  float dum_vers = 3.04F;
  if (!ncFile->add_att("api_version", dum_vers)) {
    mWriteIface->report_error("WriteNCDF: failed to define api_version attribute");
    return (MB_FAILURE);
  }
  dum_vers = 2.05F;
  if (!ncFile->add_att("version", dum_vers)) {
    mWriteIface->report_error("WriteNCDF: failed to define version attribute");
    return (MB_FAILURE);
  }
  int dum_siz = sizeof(double);
  if (!ncFile->add_att("floating_point_word_size", dum_siz)) {
    mWriteIface->report_error("WriteNCDF: failed to define floating pt word size attribute");
    return (MB_FAILURE);
  }

    // set up number of dimensions

  NcDim *num_el_blk, *num_elem, *num_nodes, *num_dim;
  if (!(num_dim = ncFile->add_dim("num_dim", (long)mesh_info.num_dim)))
  {
    mWriteIface->report_error("WriteNCDF: failed to define number of dimensions");
    return (MB_FAILURE);
  }

  if (!(num_nodes = ncFile->add_dim("num_nodes", (long)mesh_info.num_nodes)))
  {
    mWriteIface->report_error("WriteNCDF: failed to define number of nodes");
    return (MB_FAILURE);
  }

  if (!(num_elem = ncFile->add_dim("num_elem", (long)mesh_info.num_elements)))
  {
    mWriteIface->report_error("WriteNCDF: failed to define number of elements");
    return (MB_FAILURE);
  }

  if (!(num_el_blk = ncFile->add_dim("num_el_blk", (long)mesh_info.num_elementblocks)))
  {
    mWriteIface->report_error("WriteNCDF: failed to define number of element blocks");
    return (MB_FAILURE);
  }

/* ...and some variables */

    /* element block id status array */
  NcVar *idstat;
  if ((idstat = ncFile->add_var("eb_status", ncLong, num_el_blk)) == NULL)
  {
    mWriteIface->report_error("WriteNCDF: failed to define element block status array");
    return (MB_FAILURE);
  }

    /* element block id array */

  NcVar *idarr;
  if ((idarr = ncFile->add_var("eb_prop1", ncLong, num_el_blk)) == NULL)
  {
    mWriteIface->report_error("WriteNCDF: failed to define element block id array");
    return (MB_FAILURE);
  }

    /*   store property name as attribute of property array variable */
  if (!idarr->add_att("name", "ID"))
  {
    mWriteIface->report_error("WriteNCDF: failed to store element block property name ID");
    return (MB_FAILURE);
  }

    // define element blocks


  char wname[80];
  for (unsigned int i = 0; i < block_data.size(); i++)
  {
    MaterialSetData block = block_data[i];

    element_block_index = i+1;
    NcDim *num_el_in_blk, *num_att_in_blk = NULL;
    NcVar *blk_attrib, *connect;

      /* define number of elements in this block */
    
    INS_ID(wname, "num_el_in_blk%d", element_block_index);
    if ((num_el_in_blk = ncFile->add_dim(wname, (long)block.number_elements))
        == NULL)
    {
      mWriteIface->report_error("WriteNCDF: failed to define number of elements/block for block %d", 
                                i+1);
      return (MB_FAILURE);
    }
    

      /* define number of nodes per element for this block */
    INS_ID(wname, "num_nod_per_el%d", element_block_index);
    NcDim *num_nod_per_el;
    if ((num_nod_per_el = ncFile->add_dim(wname, (long)block.number_nodes_per_element))
         == NULL)
    {
      mWriteIface->report_error("WriteNCDF: failed to define number of nodes/element for block %d", 
                                block.id);
      return (MB_FAILURE);
    }

      /* define element attribute array for this block */
    if (block.number_attributes > 0)
    {
      INS_ID(wname, "num_att_in_blk%d", element_block_index);
      if ((num_att_in_blk = ncFile->add_dim(wname, (long)block.number_attributes)) == NULL)
      {
        mWriteIface->report_error("WriteNCDF: failed to define number of attributes in block %d", block.id);
        return (MB_FAILURE);
      }
      
      INS_ID(wname, "attrib%d", element_block_index);
      if ((blk_attrib = ncFile->add_var(wname, ncDouble, 
                                        num_el_in_blk, num_att_in_blk)) == NULL)
      {
        mWriteIface->report_error("WriteNCDF: failed to define attributes for element block %d", block.id);
        return (MB_FAILURE);
      }
    }


      /* define element connectivity array for this block */

    INS_ID(wname, "connect%d", element_block_index);
    if ((connect = ncFile->add_var(wname, ncLong, 
                                   num_el_in_blk, num_nod_per_el)) == NULL)
    {
      mWriteIface->report_error("WriteNCDF: failed to create connectivity array for block %d", i+1);
      return (MB_FAILURE);
    }

      /* store element type as attribute of connectivity variable */

    char element_type_string[ExoIIInterface::MAX_STR_LENGTH+1];
    strcpy(element_type_string, ExoIIUtil::ElementTypeNames[ block.element_type ]);

    if (!connect->add_att("elem_type", element_type_string))
    {
      mWriteIface->report_error("WriteNCDF: failed to store element type name %d", (int)block.element_type);
      return (MB_FAILURE);
    }
  }


/* node set id array: */

  int non_empty_nss = 0;
    // need to go through nodesets to compute # nodesets, some might be empty
  std::vector<DirichletSetData>::iterator ns_it;
  for( ns_it = nodeset_data.begin();
       ns_it != nodeset_data.end(); ns_it++) {
    if (0 != (*ns_it).number_nodes) non_empty_nss++;
  }

  NcDim *num_ns = NULL;
  NcVar *ns_idstat = NULL, *ns_idarr = NULL;
  if (non_empty_nss > 0)
  {
    if ((num_ns = ncFile->add_dim("num_node_sets", (long)(non_empty_nss))) == NULL)
    {
      mWriteIface->report_error("WriteNCDF: failed to define number of node sets");
      return (MB_FAILURE);
    }

      /* node set id status array: */

    if ((ns_idstat = ncFile->add_var("ns_status", ncLong, num_ns)) == NULL)
    {
      mWriteIface->report_error("WriteNCDF: failed to create node sets status array");
      return (MB_FAILURE);
    }

      /* node set id array: */
    if ((ns_idarr = ncFile->add_var("ns_prop1", ncLong, num_ns)) == NULL)
    {
      mWriteIface->report_error("WriteNCDF: failed to create node sets property array");
      return (MB_FAILURE);
    }


/*   store property name as attribute of property array variable */
    if (!ncFile->add_att("name", "ID"))
    {
      mWriteIface->report_error("WriteNCDF: failed to store node set property name ID");
      return (MB_FAILURE);
    }

      // now, define the arrays needed for each node set


    int index = 0;

    for(unsigned i = 0; i <nodeset_data.size(); i++)
    {
      DirichletSetData node_set = nodeset_data[i];

      if (node_set.number_nodes == 0) {
        mWriteIface->report_error("WriteNCDF: empty nodeset %d", node_set.id);
         continue;
       }
       index++;

       NcDim *num_nod_ns;
       INS_ID(wname, "num_nod_ns%d", index);
       if ((num_nod_ns = ncFile->add_dim(wname, (long)node_set.number_nodes)) == NULL)
       {
         mWriteIface->report_error("WriteNCDF: failed to define number of nodes for set %d", node_set.id);
         return (MB_FAILURE);
       }

/* create variable array in which to store the node set node list */
       NcVar *node_ns;
       INS_ID(wname, "node_ns%d", index);
       if ((node_ns = ncFile->add_var(wname, ncLong, num_nod_ns)) == NULL)
       {
         mWriteIface->report_error("WriteNCDF: failed to create node set %d node list", node_set.id);
         return (MB_FAILURE);
       }

       // create distribution factor array
       NcVar *fact_ns;
       INS_ID(wname, "dist_fact_ns%d", index);
       if ((fact_ns = ncFile->add_var(wname, ncDouble, num_nod_ns)) == NULL)
       {
         mWriteIface->report_error("WriteNCDF: failed to create node set %d distribution factor list", 
                                   node_set.id);
         return (MB_FAILURE);
       }

     }

   }

/* side set id array: */

  long non_empty_ss = 0;
    // need to go through nodesets to compute # nodesets, some might be empty
  std::vector<NeumannSetData>::iterator ss_it;
  for( ss_it = sideset_data.begin();
       ss_it != sideset_data.end(); ss_it++) {
    if (0 != (*ss_it).number_elements) non_empty_ss++;
  }

   if (non_empty_ss > 0) {
     NcDim *num_ss;
     if ((num_ss = ncFile->add_dim("num_side_sets", non_empty_ss)) == NULL)
     {
       mWriteIface->report_error("WriteNCDF: failed to define number of side sets");
       return (MB_FAILURE);
     }

     /* side set id status array: */
     NcVar *ss_idstat, *ss_idarr;
     if ((ss_idstat = ncFile->add_var("ss_status", ncLong, num_ss)) == NULL)
     {
       mWriteIface->report_error("WriteNCDF: failed to define side set status");
       return (MB_FAILURE);
     }

     /* side set id array: */
     if ((ss_idarr = ncFile->add_var("ss_prop1", ncLong, num_ss)) == NULL)
     {
       mWriteIface->report_error( "WriteNCDF: failed to define side set property");
       return (MB_FAILURE);
     }

/*   store property name as attribute of property array variable */
     if ((!ss_idarr->add_att("name", "ID")))
     {
       mWriteIface->report_error("WriteNCDF: failed to store side set property name ID");
       return (MB_FAILURE);
     }

     // now, define the arrays needed for each side set

     int index = 0;
     for(unsigned int i = 0; i < sideset_data.size(); i++)
     {
       NeumannSetData side_set = sideset_data[i];

       // dont define an empty set
       if (side_set.number_elements == 0)
         continue;

       index++;
       
       NcDim *num_side_ss;
       NcVar *elem_ss, *side_ss;
       INS_ID(wname, "num_side_ss%d", index);
       if ((num_side_ss = ncFile->add_dim(wname, (long)side_set.number_elements)) == NULL)
       {
         mWriteIface->report_error("WriteNCDF: failed to define number of sides in side set %d",
           side_set.id);
         return(MB_FAILURE);
       }

       INS_ID(wname, "elem_ss%d", index);
       if ((elem_ss = ncFile->add_var(wname, ncLong, num_side_ss)) == NULL)
       {
         mWriteIface->report_error("WriteNCDF: failed to create element list for side set %d",
                                   side_set.id);
         return(MB_FAILURE);            /* exit define mode and return */
       }
       INS_ID(wname, "side_ss%d", index);
       if ((side_ss = ncFile->add_var(wname, ncLong, num_side_ss)) == NULL)
       {
         mWriteIface->report_error("WriteNCDF: failed to create side list for side set %d",
                                   side_set.id);
         return(MB_FAILURE);         /* exit define mode and return */

       }
           
       //  sideset distribution factors
       NcDim *num_df_ss;
       INS_ID(wname, "num_df_ss%d", index);
       if ((num_df_ss = ncFile->add_dim(wname, (long)side_set.ss_dist_factors.size())) == NULL)
       {
         mWriteIface->report_error("WriteNCDF: failed to define number of dist factors in side set %d",
             side_set.id);
         return(MB_FAILURE);          /* exit define mode and return */
       }

/* create variable array in which to store the side set distribution factors */

       NcVar *fact_ss;
       INS_ID(wname, "dist_fact_ss%d", index);
       if ((fact_ss = ncFile->add_var(wname, ncDouble, num_df_ss)) == NULL)
       {
         mWriteIface->report_error("WriteNCDF: failed to create dist factors list for side set %d",
                                   side_set.id);
         return(MB_FAILURE);            /* exit define mode and return */
       }
     }
   }

/* node coordinate arrays: */

   NcVar *coord, *name_coord;
   if ((coord = ncFile->add_var("coord", ncDouble, num_dim, num_nodes)) == NULL)
   {
     mWriteIface->report_error("WriteNCDF: failed to define node coordinate array");
     return (MB_FAILURE);
   }
/* coordinate names array */

   if ((name_coord = ncFile->add_var("coor_names", ncChar, num_dim, dim_str)) == NULL)
   {
     mWriteIface->report_error("WriteNCDF: failed to define coordinate name array");
     return (MB_FAILURE);
   }

   // define genesis maps if required

   if (write_maps)
   {
     // element map
     NcVar *elem_map, *elem_map2, *node_map;
     if ((elem_map = ncFile->add_var("elem_map", ncLong, num_elem)) == NULL)
     {
       mWriteIface->report_error("WriteNCDF: failed to create element map array");
       return (MB_FAILURE);         /* exit define mode and return */
     }

     // create the element number map
     if ((elem_map2 = ncFile->add_var("elem_num_map", ncLong, num_elem)) == NULL)
     {
       mWriteIface->report_error("WriteNCDF: failed to create element numbering map");
     }
       // create node number map
     if ((node_map = ncFile->add_var("node_num_map", ncLong, num_nodes)) == NULL)
     {
       mWriteIface->report_error("WriteNCDF: failed to create node numbering map array");
       return(MB_FAILURE);         /* exit define mode and return */
     }
   }

   // define qa records to be used

   int num_qa_rec = mesh_info.qaRecords.size()/4;
   NcDim *num_qa;

   if ((num_qa = ncFile->add_dim("num_qa_rec",(long)num_qa_rec)) == NULL)
   {
     mWriteIface->report_error("WriteNCDF: failed to define qa record array size");
     return (MB_FAILURE);
   }

   // define qa array
   NcVar *qa_title;
   if ((qa_title = ncFile->add_var("qa_records", ncChar, num_qa, dim_four, dim_str)) == NULL)
   {
     mWriteIface->report_error("WriteNCDF: failed to define qa record array");
     return (MB_FAILURE);
   }

   return MB_SUCCESS;
}


ErrorCode WriteNCDF::open_file(const char* filename)
{
   // not a valid filname
   if(strlen((const char*)filename) == 0)
   {
     mWriteIface->report_error("Output Exodus filename not specified");
      return MB_FAILURE;
   }

   ncFile = new NcFile(filename, NcFile::Replace);

   // file couldn't be opened
   if(ncFile == NULL)
   {
     mWriteIface->report_error("Cannot open %s", filename);
     return MB_FAILURE;
   }
   return MB_SUCCESS;
}

ErrorCode WriteNCDF::get_sideset_elems(EntityHandle sideset, int current_sense,
                                           Range &forward_elems, Range &reverse_elems) 
{
  Range ss_elems, ss_meshsets;

    // get the sense tag; don't need to check return, might be an error if the tag
    // hasn't been created yet
  Tag sense_tag = 0;
  mdbImpl->tag_get_handle("SENSE", sense_tag);

    // get the entities in this set
  ErrorCode result = mdbImpl->get_entities_by_handle(sideset, ss_elems, true);
  if (MB_FAILURE == result) return result;
  
    // now remove the meshsets into the ss_meshsets; first find the first meshset,
  Range::iterator range_iter = ss_elems.begin();
  while (TYPE_FROM_HANDLE(*range_iter) != MBENTITYSET && range_iter != ss_elems.end())
    range_iter++;
  
    // then, if there are some, copy them into ss_meshsets and erase from ss_elems
  if (range_iter != ss_elems.end()) {
    std::copy(range_iter, ss_elems.end(), range_inserter(ss_meshsets));
    ss_elems.erase(range_iter, ss_elems.end());
  }
  

    // ok, for the elements, check the sense of this set and copy into the right range
    // (if the sense is 0, copy into both ranges)

    // need to step forward on list until we reach the right dimension
  Range::iterator dum_it = ss_elems.end();
  dum_it--;
  int target_dim = MBCN::Dimension(TYPE_FROM_HANDLE(*dum_it));
  dum_it = ss_elems.begin();
  while (target_dim != MBCN::Dimension(TYPE_FROM_HANDLE(*dum_it)) &&
         dum_it != ss_elems.end()) 
    dum_it++;

  if (current_sense == 1 || current_sense == 0)
    std::copy(dum_it, ss_elems.end(), range_inserter(forward_elems));
  if (current_sense == -1 || current_sense == 0)
    std::copy(dum_it, ss_elems.end(), range_inserter(reverse_elems));
  
    // now loop over the contained meshsets, getting the sense of those and calling this
    // function recursively
  for (range_iter = ss_meshsets.begin(); range_iter != ss_meshsets.end(); range_iter++) {

      // first get the sense; if it's not there, by convention it's forward
    int this_sense;
    if (0 == sense_tag ||
        MB_FAILURE == mdbImpl->tag_get_data(sense_tag, &(*range_iter), 1, &this_sense))
      this_sense = 1;
      
      // now get all the entities on this meshset, with the proper (possibly reversed) sense
    get_sideset_elems(*range_iter, this_sense*current_sense,
                      forward_elems, reverse_elems);
  }
  
  return result;
}


} // namespace moab

  
