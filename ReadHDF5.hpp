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
 * \class  ReadHDF5
 * \brief  Read mesh from MOAB HDF5 (.h5m) file.
 * \author Jason Kraftcheck
 * \date   18 April 2004
 */

#ifndef READ_HDF5_HPP
#define READ_HDF5_HPP

#include <stdlib.h>
#include <list>
#ifdef USE_MPI
#  include "moab_mpi.h"
#endif
#include "mhdf.h"
#include "moab/Forward.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Range.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/RangeMap.hpp"

namespace moab {

class ReadHDF5 : public ReaderIface
{
public:

  static ReaderIface* factory( Interface* );

  ReadHDF5( Interface* iface );
  
  virtual ~ReadHDF5();
  
  /** Export specified meshsets to file
   * \param filename     The filename to export.  Must end in <em>.mhdf</em>
   * \param export_sets  Array of handles to sets to export, or NULL to export all.
   * \param export_set_count Length of <code>export_sets</code> array.
   */
  ErrorCode load_file( const char* filename,
                         const EntityHandle* file_set,
                         const FileOptions& opts,
                         const ReaderIface::IDTag* subset_list = 0,
                         int subset_list_length = 0,
                         const Tag* file_id_tag = 0 );

  ErrorCode read_tag_values( const char* file_name,
                               const char* tag_name,
                               const FileOptions& opts,
                               std::vector<int>& tag_values_out,
                               const IDTag* subset_list = 0,
                               int subset_list_length = 0 );

protected:

  ErrorCode load_file_impl( const FileOptions& opts );

  ErrorCode load_file_partial( const ReaderIface::IDTag* subset_list,
                                 int subset_list_length,
                                 const FileOptions& opts );

  ErrorCode read_tag_values_all( int tag_index, std::vector<int>& results );
  ErrorCode read_tag_values_partial( int tag_index, const Range& file_ids,
                                       std::vector<int>& results );

private:
  ErrorCode init();
  
  inline int is_error( mhdf_Status& status ) {
    int i;
    if ((i = mhdf_isError(&status))) 
      readUtil->report_error( "%s", mhdf_message(&status) );
    return i;
  }
  
  //! The size of the data buffer (<code>dataBuffer</code>).
  int bufferSize;
  //! A memory buffer to use for all I/O operations.
  char* dataBuffer;

  //! Interface pointer passed to constructor
  Interface* iFace;
  
  //! The file handle from the mhdf library
  mhdf_FileHandle filePtr;
  
  //! File summary
  mhdf_FileDesc* fileInfo;
  
  //! Map from File ID to MOAB handle
  typedef RangeMap< long, EntityHandle > IDMap;
  IDMap idMap;
  
  //! Cache pointer to read util
  ReadUtilIface* readUtil;
  
  //! The type of an EntityHandle
  hid_t handleType;
  
  //! read/write property handle
  //! indepIO -> idependent IO during true parallel read
  //! collIO  -> collective IO during true parallel read
  //! Both are H5P_DEFAULT for serial IO and collective
  //! when reading the entire file on all processors.
  hid_t indepIO, collIO;
  
  ErrorCode set_up_read( const char* file_name, const FileOptions& opts );
  ErrorCode clean_up_read( const FileOptions& opts );
  
  
  //! Given a list of tags and values, get the file ids for the
  //! corresponding entities in the file.
  ErrorCode get_subset_ids( const ReaderIface::IDTag* subset_list,
                              int subset_list_length,
                              Range& file_ids_out );
  
  ErrorCode read_nodes( const Range& node_file_ids );
  
    // Read elements in fileInfo->elems[index]
  ErrorCode read_elems( int index );
  
    // Read subset of elements in fileInfo->elems[index]
  ErrorCode read_elems( int index, const Range& file_ids );
  
  //! Read element connectivity.
  ErrorCode read_elems( const mhdf_ElemDesc& elems, const Range& file_ids );
  
  
    // Read connectivity data for a list of element file ids.
    // passing back the file IDs for the element connectivity 
    // w/out actually creating any elements in MOAB.
  ErrorCode read_elems( int index, const Range& element_file_ids, Range& node_file_ids );

    // Scan all elements in group.  For each element for which all vertices
    // are contained in idMap (class member), read the element.  All new
    // elements are added to idMap.
    //
    // NOTE: Collective IO calls in parallel.
  ErrorCode read_node_adj_elems( const mhdf_ElemDesc& group );
  
    // Scan all elements in specified file table.  For each element for 
    // which all vertices are contained in idMap (class member), read the 
    // element.  All new elements are added to idMap.
    //
    // NOTE: Collective IO calls in parallel.
  ErrorCode read_node_adj_elems( const mhdf_ElemDesc& group,
                                   hid_t connectivity_handle );

  //! Read poly(gons|hedra)
  ErrorCode read_poly( const mhdf_ElemDesc& elems, const Range& file_ids );
  
  //! Read specified elements and any adjacent elements of lower dimension.
  //! Assumes vertices are already read in.
  ErrorCode read_elements_and_sides( const Range& file_ids );
  
  //! Read sets
  ErrorCode read_sets( const Range& set_file_ids );
  
  //! Read set contents
  ErrorCode read_set_contents( hid_t set_description_handle,
                                 hid_t set_contents_handle,
                                 const unsigned long data_len );
  
  //! Read set parents/children
  ErrorCode read_parents_children( bool parents, 
                                     hid_t set_description_handle,
                                     hid_t set_contents_handle,
                                     const unsigned long data_len );
  
  ErrorCode read_adjacencies( hid_t adjacency_table,
                                long table_length );
                                
  
  //! Create tag and read all data.
  ErrorCode read_tag( int index );
  
  //! Create new tag or varify type matches existing tag
  ErrorCode create_tag( const mhdf_TagDesc& info, Tag& handle, hid_t& type );
  
  
  //! Read dense tag for all entities 
  ErrorCode read_dense_tag( Tag tag_handle,
                              hid_t hdf_read_type,
                              hid_t data_table,
                              long start_id,
                              long count );

  
  //! Read sparse tag for all entities.
  ErrorCode read_sparse_tag( Tag tag_handle,
                               hid_t hdf_read_type,
                               hid_t ent_table,
                               hid_t val_table,
                               long num_entities );
  
  //! Read variable-length tag for all entities.
  ErrorCode read_var_len_tag( Tag tag_handle,
                                hid_t hdf_read_type,
                                hid_t ent_table,
                                hid_t val_table,
                                hid_t off_table,
                                long num_entities,
                                long num_values );
                               
  ErrorCode read_qa( EntityHandle file_set );
                               
  ErrorCode convert_id_to_handle( EntityHandle* in_out_array,
                                    size_t array_length );
                                    
  ErrorCode convert_range_to_handle( const EntityHandle* ranges,
                                       size_t num_ranges,
                                       Range& merge );
                                    
  static
  void convert_id_to_handle( EntityHandle* in_out_array,
                             size_t array_length,
                             const RangeMap<long,EntityHandle>& id_map );
                                    
  static
  void convert_id_to_handle( EntityHandle* in_out_array,
                             size_t array_length,
                             size_t& array_length_out,
                             const RangeMap<long,EntityHandle>& id_map );

  static
  void convert_range_to_handle( const EntityHandle* ranges,
                                size_t num_ranges,
                                const RangeMap<long,EntityHandle>& id_map,
                                Range& merge );
  
  /**\brief Search for entities with specified tag values 
   * 
   *\NOTE For parallel reads, this function does collective IO.
   *
   *\param tag_index  Index into info->tags specifying which tag to search.
   *\param sorted_values  List of tag values to check for, in ascending sorted
   *                  order.
   *\param file_ids_out  File IDs for entities with specified tag values.
   */
  ErrorCode search_tag_values( int tag_index,
                                 const std::vector<int>& sorted_values,
                                 Range& file_ids_out );
  
  /**\brief Search for entities with specified tag 
   * 
   *\NOTE For parallel reads, this function does collective IO.
   *
   *\param tag_index  Index into info->tags specifying which tag to search.
   *\param file_ids_out  File IDs for entities with specified tag values.
   */
  ErrorCode get_tagged_entities( int tag_index, Range& file_ids_out );
                                 
  /**\brief Search a table of tag data for a specified set of values.
   *
   * Search a table of tag values, returning the indices into the table
   * at which matches were found.
   *\NOTE For parallel reads, this function does collective IO.
   *
   *\param info       Summary of data contained in file.
   *\param tag_table     HDF5/mhdf handle for tag values
   *\param table_size    Number of values in table
   *\param sorted_values Sorted list of values to search for.
   *\param value_indices Output: Offsets into the table of data at which 
   *                       matching values were found.
   */
  ErrorCode search_tag_values( hid_t tag_table, 
                                 unsigned long table_size,
                                 const std::vector<int>& sorted_values,
                                 std::vector<EntityHandle>& value_indices );
  
  /**\brief Get the file IDs for nodes and elements contained in sets.
   *
   * Read the contents for the specified sets and return the file IDs
   * of all nodes and elements contained within those sets.
   *\param sets       Container of file IDs designating entity sets.
   *\param file_ids   Output: File IDs of entities contained in sets.
   */
  ErrorCode get_set_contents( const Range& sets, Range& file_ids );
 
  /** Given a list of file IDs for entity sets, find all contained
   *  or child sets (at any depth) and append them to the Range 
   *  of file IDs.
   */
  ErrorCode read_set_ids_recursive( Range& sets_in_out,
                                      bool containted_sets,
                                      bool child_sets );
  
  /** Find all sets containing one or more entities read from the file
   *  and added to idMap 
   */
  ErrorCode find_sets_containing( Range& sets_out );
 
  /**\brief Read sets from file into MOAB for partial read of file.
   *
   * Given the file IDs for entity sets (sets_in) and elements and 
   * nodes (id_map), read in all sets containing any of the elements
   * or nodes and all sets that are (recursively) children of any 
   * other set to be read (those in sets_in or those containging any
   * already-read element or node.)
   *\param sets_in    File IDs for sets to read (unconditionally)
   */
  ErrorCode read_sets_partial( const Range& sets_in );

  /** Find file IDs of sets containing any entities in the passed id_map */
  ErrorCode find_sets_containing( hid_t meta_handle,
                                    hid_t content_handle, 
                                    long content_len,
                                    Range& file_ids );  
 
  /** Given a list of file IDs for entity sets, read the list of 
   *  file IDs for all child entity sets.
   */
  ErrorCode read_child_ids( const Range& set_file_ids,
                              hid_t meta_handle,
                              hid_t child_handle,
                              Range& child_file_ids );
 
  /** Given a list of file IDs for entity sets, read the list of 
   *  file IDs for all contained entity sets.
   */
  ErrorCode read_contained_set_ids( const Range& set_file_ids,
                                      hid_t meta_handle,
                                      hid_t contents_handle,
                                      Range& containd_set_file_ids );
    
    /**\brief Create sets 
     *
     * For the list of entity sets designated by the file IDs contained
     * in file_ids, instantiate the sets in MOAB with consecutive handles.
     *
     *\param info      Summary of data contained in the file.
     *\param file_ids  List of file IDs designating which sets will be read.
     *\param set_meta_handle HDF5/mhdf handle for set metadata table.
     *\parma id_map    Map from file IDs to EntityHandle for entities read
     *                 from file.  Sets created by this function will be appended.
     *\param ranged_file_ids_out This function will add to this container
     *                 the file IDs for any sets for which the contents are
     *                 stored using ranged compression (any set for which
     *                 mhdf_SET_RANGE_BIT is set in the flags.)
     */
  ErrorCode read_sets( const Range& file_ids,
                         hid_t set_meta_handle, 
                         Range& ranged_file_ids_out,
                         EntityHandle& start_handle,
                         bool create = true );
   
   ErrorCode read_contents( const Range& set_file_ids,
                              EntityHandle start_handle,
                              hid_t set_meta_data_table,
                              hid_t set_contents_table,
                              long set_contents_length,
                              const Range& ranged_set_file_ids );

   ErrorCode read_parents( const Range& set_file_ids,
                             EntityHandle start_handle,
                             hid_t set_meta_data_table,
                             hid_t set_parents_table,
                             long set_parents_length );

   ErrorCode read_children( const Range& set_file_ids,
                              EntityHandle start_handle,
                              hid_t set_meta_data_table,
                              hid_t set_children_table,
                              long set_children_length );
   

  class ContentReader {
    public:
      virtual void 
      read_indices( long, long, long*, mhdf_Status& ) = 0;
      
      virtual void 
      read_contents( long, long, EntityHandle*, mhdf_Status& ) = 0;
      
      virtual ErrorCode 
      store_data( EntityHandle, long file_id, EntityHandle*, long len, bool ranged = false ) = 0;
   };
   
   /**\brief Read end-indexed variable length handle lists such as
     *       set parents, set children, set contents, and old-format
     *       poly connectivity.
     *
     * Read the contents of each entity set from the file, and add to
     * the corresponding entity set in MOAB any contained entities that
     * have been read from the file.
     *\param file_ids        File IDs of entity sets to read.
     *\param start_id        ID of first entity in table.
     *\param start_handle    EntityHandle for first entity set in file_ids
     *                       (assumes all sets in file_ids have sequential
     *                        handles.)
     *\param set_content_len Length of set contents table.
     *\param ranged_file_ids Subset of file_ids for which set contents
     *                       are stored in ranged format.
     */
  ErrorCode read_contents( ContentReader& tool,
                             const Range& file_ids,
                             const long start_id,
                             const EntityHandle start_handle,
                             const long num_sets,
                             const long set_content_len,
                             const Range& ranged_file_ids );
  
    /**\brief Store file IDS in tag values
     *
     * Copy fild ID from IDMap for each entity read from file
     * into a tag value on the entity.
     */
  ErrorCode store_file_ids( Tag tag );
  
    /**\brief Find index in mhdf_FileDesc* fileInfo for specified tag name
     *
     * Given a tag name, find its index in fileInfo and verify that
     * each tag value is a single integer.
     */
  ErrorCode find_int_tag( const char* name, int& index_out );
};

} // namespace moab

#endif
