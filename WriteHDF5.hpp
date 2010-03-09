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
 * \class  WriteHDF5
 * \brief  Write mesh database to TSTT HDF5 file.
 * \author Jason Kraftcheck
 * \date   01 April 2004
 */

#ifndef WRITE_HDF5_HPP
#define WRITE_HDF5_HPP

#include <list>
#ifdef USE_MPI // include this before HDF5 headers to avoid conflicts
#  include "MBmpi.h"
#endif
#include "mhdf.h"
#include "MBForward.hpp"
#include "MBRange.hpp"
#include "MBWriterIface.hpp"
#include "RangeMap.hpp"

class MBWriteUtilIface;

/* If this define is not set, node->entity adjacencies will not be written */
#undef MB_H5M_WRITE_NODE_ADJACENCIES

class MB_DLL_EXPORT WriteHDF5 : public MBWriterIface
{
public:

  static MBWriterIface* factory( MBInterface* );

  /** The type to use for entity IDs w/in the file.
   * 
   * NOTE:  If this is changed, the value of id_type 
   *        MUST be changed accordingly.
   */
  typedef MBEntityHandle id_t;
  
  /** HDF5 type corresponding to type of id_t */
  static const hid_t id_type;

  WriteHDF5( MBInterface* iface );
  
  virtual ~WriteHDF5();
  
  /** Export specified meshsets to file
   * \param filename     The filename to export. 
   * \param export_sets  Array of handles to sets to export, or NULL to export all.
   * \param export_set_count Length of <code>export_sets</code> array.
   */
  MBErrorCode write_file( const char* filename,
                          const bool overwrite,
                          const FileOptions& opts,
                          const MBEntityHandle* export_sets,
                          const int export_set_count,
                          const std::vector<std::string>& qa_records,
                          const MBTag* tag_list = 0,
                          int num_tags = 0,
                          int user_dimension = 3 );

  /** Create attributes holding the HDF5 type handle for the 
   *  type of a bunch of the default tags.
   */
  //static MBErrorCode register_known_tag_types( MBInterface* );

protected:
  
  MBErrorCode serial_create_file( const char* filename,
                                  bool overwrite,
                                  const std::vector<std::string>& qa_records,
                                  const MBTag* tag_list,
                                  int num_tags,
                                  int dimension = 3 );

  /** Function to create the file.  Virtual to allow override
   *  for parallel version.
   */
  virtual MBErrorCode parallel_create_file( const char* filename,
                                            bool overwrite,
                                            const std::vector<std::string>& qa_records,
                                            const MBTag* tag_list,
                                            int num_tags,
                                            int dimension = 3,
                                            int pcomm_no = 0 );


  /** Functions that the parallel version overrides*/
  virtual MBErrorCode write_shared_set_descriptions( hid_t ) 
    { return MB_SUCCESS;}
  virtual MBErrorCode write_shared_set_contents( hid_t )
    { return MB_SUCCESS;}
  virtual MBErrorCode write_shared_set_children( hid_t )
    { return MB_SUCCESS;}
  virtual MBErrorCode write_shared_set_parents( hid_t )
    { return MB_SUCCESS;}
  virtual MBErrorCode write_finished();
  virtual void tprint( const char* fmt, ... )
#ifdef __GNUC__
__attribute__((format(printf,2,3)))
#endif
  ;

 
  //! Gather tags
  MBErrorCode gather_tags( const MBTag* user_tag_list, int user_tag_list_length );

  /** Helper function for create-file
   *
   * Calculate the sum of the number of non-set adjacencies
   * of all entities in the passed range.
   */
  MBErrorCode count_adjacencies( const MBRange& elements, id_t& result );
  
  /** Helper function for create-file
   *
   * Create zero-ed tables where element connectivity and 
   * adjacency data will be stored.
   */
  MBErrorCode create_elem_tables( MBEntityType mb_type,
                                  int nodes_per_element,
                                  id_t number_elements,
                                  long& first_id_out );
  
  /** Helper function for create-file
   *
   * Calculate total length of set contents and child tables.
   */
  MBErrorCode count_set_size( const MBRange& sets,
                              long& contents_length_out,
                              long& children_length_out,
                              long& parents_length_out );
  
  /** Helper function for create-file
   *
   * Create zero-ed table where set descriptions will be written
   */
  MBErrorCode create_set_meta( id_t number_sets, long& first_id_out );

  /** Helper function for create-file
   *
   * Create zero-ed tables where set data will be written.
   */
  MBErrorCode create_set_tables( long contents_length, 
                                 long children_length,
                                 long parents_length );

  //! Write exodus-type QA info
  MBErrorCode write_qa( const std::vector<std::string>& list );


  //! Range of entities, grouped by type, to export 
  struct ExportSet 
  {
    //! The range of entities.
    MBRange range;
    //! The type of the entities in the range
    MBEntityType type;
    //! The number of nodes per entity - not used for nodes and sets
    int num_nodes;
    //! The first Id allocated by the mhdf library.  Entities in range have sequential IDs.
    id_t first_id;
    //! The offset at which to begin writting this processor's data.
    //! Always zero except for parallel IO.
    id_t offset;
    //! Offset for adjacency data.  Always zero except for parallel IO
    MBEntityID adj_offset;
    //! If doing parallel IO, largest number of entities to write
    //! for any processor (needed to do collective IO).  Zero if unused.
    long max_num_ents, max_num_adjs;
    
    bool operator<( const ExportSet& other ) const
      { return type < other.type || 
               (type == other.type && 
                type != MBPOLYGON &&
                type != MBPOLYHEDRON &&
                num_nodes < other.num_nodes); }
                
    const char* name() const;
  };
  
public:
  //! Tag to write to file.
  struct SparseTag
  {
    //! The tag handle
    MBTag tag_id;
    //! The list of entities with this tag
    MBRange range;
    //! The offset at which to begin writting this processor's data.
    //! Always zero except for parallel IO. 
    id_t offset;
    //! For variable-length tags, a second offset for the tag data table,
    //! separate from the offset used for the ID and Index tables.
    //! Always zero except for parallel IO. 
    id_t varDataOffset;
    //! Write tag data (for serial, is always equal to !range.empty())
    bool write;
    //! If doing parallel IO, largest number of tag values to write
    //! for any processor (needed to do collective IO).  Zero if unused.
    long max_num_ents;
    
    
    bool operator<(const SparseTag&) const;
  };
protected:
  
  //! The size of the data buffer (<code>dataBuffer</code>).
  size_t bufferSize;
  //! A memory buffer to use for all I/O operations.
  char* dataBuffer;

  //! MBInterface pointer passed to constructor
  MBInterface* iFace;
  //! Cached pointer to writeUtil interface.
  MBWriteUtilIface* writeUtil;
  
  //! The file handle from the mhdf library
  mhdf_FileHandle filePtr;
  
  //! Map from entity handles to file IDs
  RangeMap<MBEntityHandle,id_t> idMap;
  
  //! The list elements to export.
  std::list<ExportSet> exportList;
  //! The list of nodes to export
  ExportSet nodeSet;
  //! The list of sets to export
  ExportSet setSet;
  
  //! Offset into set contents table (zero except for parallel)
  unsigned long setContentsOffset;
  //! Offset into set children table (zero except for parallel)
  unsigned long setChildrenOffset, setParentsOffset;
  //! If doing parallel IO, largest number of values to write
  //! for any processor (needed to do collective IO).  Zero if unused.
  long maxNumSetContent, maxNumSetChildren, maxMumSetParents;
  //! Flags idicating if set data should be written.
  //! For the normal (non-parallel) case, these values
  //! will depend only on whether or not there is any
  //! data to be written.  For parallel-meshes, opening
  //! the data table is collective so the values must
  //! depend on whether or not any processor has meshsets
  //! to be written.
  bool writeSets, writeSetContents, writeSetChildren, writeSetParents;
  
  //! The list of tags to export
  std::list<SparseTag> tagList;

  //! True if doing parallel write
  bool parallelWrite;
  
  //! Property set to pass to H5Dwrite calls. 
  //! For serial, should be H5P_DEFAULTS.
  //! For parallel, may request collective IO.
  hid_t writeProp;

  void print_id_map() const;
  void print_id_map( std::ostream& str, const char* prefix = "" ) const;
  
  
private:

  //! Do the actual work of write_file.  Separated from write_file
  //! for easier resource cleanup.
  MBErrorCode write_file_impl( const char* filename,
                               const bool overwrite,
                               const FileOptions& opts,
                               const MBEntityHandle* export_sets,
                               const int export_set_count,
                               const std::vector<std::string>& qa_records,
                               const MBTag* tag_list,
                               int num_tags,
                               int user_dimension = 3 );

  MBErrorCode init();

  //! Get information about a meshset
  MBErrorCode get_set_info( MBEntityHandle set,
                            long& num_entities,
                            long& num_children,
                            long& num_parents,
                            unsigned long& flags );

protected:

  /** Helper function for create-file
   *
   * Write tag meta-info and create zero-ed table where
   * tag values will be written.
   *\param num_entities  Number of entities for which to write tag data.
   *\param var_len_total For variable-length tags, the total number of values
   *                     in the data table.
   */
  MBErrorCode create_tag( MBTag tag_id, 
                          unsigned long num_entities,
                          unsigned long var_len_total );
  
  /**\brief add entities to idMap */
  MBErrorCode assign_ids( const MBRange& entities, id_t first_id );
  
  /** Get possibly compacted list of IDs for passed entities
   *
   * For the passed range of entities, determine if IDs
   * can be compacted and write IDs to passed list.
   *
   * If the IDs are not compacted, the output list will contain
   * a simple ordered list of IDs.
   *
   * If IDs are compacted, the output list will contain 
   * {start,count} pairs.
   *
   * If the ID list is compacted, ranged_list will be 'true'.
   * Otherwise it will be 'false'.
   */
  MBErrorCode range_to_blocked_list( const MBRange& input_range,
                                     std::vector<id_t>& output_id_list , 
                                     bool& ranged_list );
  

  MBErrorCode range_to_id_list( const MBRange& input_range,
                                id_t* array );
  //! Get IDs for entities 
  MBErrorCode vector_to_id_list( const std::vector<MBEntityHandle>& input,
                                 std::vector<id_t>& output, 
                                 bool remove_non_written = false );

  /** When writing tags containing MBEntityHandles to file, need to convert tag
   *  data from MBEntityHandles to file IDs.  This function does that. 
   *
   * If the handle is not valid or does not correspond to an entity that will
   * be written to the file, the file ID is set to zero.
   *\param data  The data buffer.  As input, an array of MBEntityHandles.  As
   *             output an array of file IDS, where the size of each integral
   *             file ID is the same as the size of MBEntityHandle.
   *\param count The number of handles in the buffer.
   *\return true if at least one of the handles is valid and will be written to
   *             the file or at least one of the handles is NULL (zero). false
   *             otherwise
   */
  bool convert_handle_tag( MBEntityHandle* data, size_t count ) const;
  bool convert_handle_tag( const MBEntityHandle* source,
                           MBEntityHandle* dest, 
                           size_t count ) const;

  /** Get IDs of adjacent entities.
   * 
   * For all entities adjacent to the passed entity, if the
   * adjacent entity is to be exported (ID is not zero), append
   * the ID to the passed list.
   */
  MBErrorCode get_adjacencies( MBEntityHandle entity, std::vector<id_t>& adj );
                                
  //! get sum of lengths of tag values (as number of type) for 
  //! variable length tag data.
  MBErrorCode get_tag_data_length( const SparseTag& tag_info,
                                   unsigned long& result );
  
private:
  
  /** Get all mesh to export from given list of sets.
   *
   * Populate exportSets, nodeSet and setSet with lists of
   * entities to write.
   *
   * \param export_sets  The list of meshsets to export
   */
  MBErrorCode gather_mesh_info( const std::vector<MBEntityHandle>& export_sets );
  
  //! Same as gather_mesh_info, except for entire mesh
  MBErrorCode gather_all_mesh( );
  
  //! Initialize internal data structures from gathered mesh
  MBErrorCode initialize_mesh( const MBRange entities_by_dim[5] );
 
  /** Write out the nodes.
   *
   * Note: Assigns IDs to nodes.
   */
  MBErrorCode write_nodes( );
  
  /** Write out element connectivity.
   *
   * Write connectivity for passed set of elements.
   *
   * Note: Assigns element IDs.
   * Note: Must do write_nodes first so node IDs get assigned.
   */
  MBErrorCode write_elems( ExportSet& elemset );
  
  /** Write out meshsets
   * 
   * Write passed set of meshsets, including parent/child relations.
   *
   * Note: Must have written nodes and element connectivity
   *       so entities have assigned IDs.
   */
  MBErrorCode write_sets( );

  MBErrorCode write_parents_children( bool children );
  
  /** Write adjacency info for passed set of elements
   *
   * Note: Must have written element connectivity so elements
   *       have IDs assigned.
   */
  MBErrorCode write_adjacencies( const ExportSet& export_set );
  
  /** Write tag information and data.
   * 
   * Note: Must have already written nodes, elem connectivity and
   *       sets so that entities have IDs assigned.
   */
/*
  MBErrorCode write_tag( MBTag tag_handle );
  
  //! Write dense tag for all entities 
  MBErrorCode write_dense_tag( MBTag tag_handle,
                               hid_t hdf_write_type );

  //! Write dense tag for specified entity set
  MBErrorCode write_dense_tag( ExportSet& set,
                               MBTag tag_handle,
                               hid_t hdf_write_type );
*/  
  //! Write sparse tag for all entities.
  MBErrorCode write_sparse_tag( const SparseTag& tag_data );
                            
  //! Get element connectivity
  MBErrorCode get_connectivity( MBRange::const_iterator begin,
                                MBRange::const_iterator end,
                                int nodes_per_element,
                                id_t* id_data_out );
                                   
  //! Get size data for tag
  //!\param tag       MOAB tag ID
  //!\param moab_type Output: MBDataType for tag
  //!\param num_bytes Output: MOAB tag size (bits for bit tags).
  //!                         MB_VARIABLE_LENGTH for variable-length tags.
  //!\param elem_size Output: Size of type values per entity (e.g.
  //!                         sizeof(double) for MB_TYPE_DOUBLE data)
  //!                         One for bit and opaque tags.
  //!\param file_size Output: num_bytes/elem_size
  //!\param file_type Output: mhdf type enumeration
  //!\param hdf_type  Output: zero or handle for user-defined custom type
  //!                         (user-defined type available only for opaque
  //!                         data.)
  MBErrorCode get_tag_size( MBTag tag,
                            MBDataType& moab_type,
                            int& num_bytes,
                            int& elem_size,
                            int& file_size,
                            mhdf_TagDataType& file_type,
                            hid_t& hdf_type );
                            
  //! Write ID table for sparse tag
  MBErrorCode write_sparse_ids( const SparseTag& tag_data, hid_t table_handle );
  
  //! Write varialbe-length tag data
  MBErrorCode write_var_len_tag( const SparseTag& tag_info );
};

#endif
