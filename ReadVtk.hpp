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

#ifndef READ_VTK_HPP
#define READ_VTK_HPP

#include "MBForward.hpp"
#include "MBReaderIface.hpp"

#include <string>

class MBReadUtilIface;
class FileTokenizer;

class ReadVtk : public MBReaderIface
{
   
public:

  static MBReaderIface* factory( MBInterface* );

    //! load a file
  MBErrorCode load_file( const char *file_name,
                         const MBEntityHandle* file_set,
                         const FileOptions&,
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
  ReadVtk(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadVtk();

protected:

  MBErrorCode allocate_vertices( long num_vtx,
                                 MBEntityHandle& start_handle_out,
                                 double*& x_coord_array_out,
                                 double*& y_coord_array_out,
                                 double*& z_coord_array_out );

  MBErrorCode read_vertices( FileTokenizer& tokens,
                             long num_verts, 
                             MBEntityHandle& start_handle_out );

  MBErrorCode allocate_elements( long num_elements,
                                 int vert_per_element,
                                 MBEntityType type,
                                 MBEntityHandle& start_handle_out,
                                 MBEntityHandle*& conn_array_out,
                                 std::vector<MBRange>& append_to_this );

  MBErrorCode vtk_read_dataset( FileTokenizer& tokens,
                                MBRange& vertex_list,
                                std::vector<MBRange>& element_list );

  MBErrorCode vtk_read_structured_points( FileTokenizer& tokens, 
                                          MBRange& vertex_list,
                                          std::vector<MBRange>& elem_list );

  MBErrorCode vtk_read_structured_grid( FileTokenizer& tokens, 
                                        MBRange& vertex_list,
                                        std::vector<MBRange>& elem_list );

  MBErrorCode vtk_read_rectilinear_grid( FileTokenizer& tokens, 
                                         MBRange& vertex_list,
                                         std::vector<MBRange>& elem_list );
                                         
  MBErrorCode vtk_read_polydata( FileTokenizer& tokens, 
                                 MBRange& vertex_list,
                                 std::vector<MBRange>& elem_list );
                                 
  MBErrorCode vtk_read_polygons( FileTokenizer& tokens,
                                 MBEntityHandle first_vtx, 
                                 std::vector<MBRange>& elem_list );

  MBErrorCode vtk_read_unstructured_grid( FileTokenizer& tokens,
                                          MBRange& vertex_list,
                                          std::vector<MBRange>& elem_list  );

  MBErrorCode vtk_create_structured_elems( const long* dims, 
                                           MBEntityHandle first_vtx,
                                           std::vector<MBRange>& elem_list );

  MBErrorCode vtk_read_field( FileTokenizer& tokens );

  MBErrorCode vtk_read_attrib_data( FileTokenizer& tokens, 
                                    std::vector<MBRange>& entities );

  MBErrorCode vtk_read_tag_data( FileTokenizer& tokens, 
                                 int type, 
                                 size_t per_elem, 
                                 std::vector<MBRange>& entities,
                                 const char* name );

  MBErrorCode vtk_read_scalar_attrib( FileTokenizer& tokens, 
                                      std::vector<MBRange>& entities,
                                      const char* name);

  MBErrorCode vtk_read_color_attrib( FileTokenizer& tokens, 
                                     std::vector<MBRange>& entities,
                                     const char* name);

  MBErrorCode vtk_read_vector_attrib( FileTokenizer& tokens, 
                                      std::vector<MBRange>& entities,
                                      const char* name);

  MBErrorCode vtk_read_texture_attrib( FileTokenizer& tokens,
                                       std::vector<MBRange>& entities,
                                       const char* name);

  MBErrorCode vtk_read_tensor_attrib( FileTokenizer& tokens,
                                      std::vector<MBRange>& entities,
                                      const char* name );

  MBErrorCode vtk_read_field_attrib( FileTokenizer& tokens, 
                                     std::vector<MBRange>& entities,
                                     const char* name);

  MBErrorCode store_file_ids( MBTag tag,
                              const MBRange& vertices,
                              const std::vector<MBRange>& elements );

private:

  MBReadUtilIface* readMeshIface;

  //------------member variables ------------//

    //! interface instance
  MBInterface* mdbImpl;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;

    //! A field which, if present and having a single integer for storage, should be used to partition the mesh by range. Defaults to MATERIAL_SET_TAG_NAME
  std::string mPartitionTagName;
};

#endif




