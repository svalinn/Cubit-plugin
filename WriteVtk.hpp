/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */


#ifndef WRITE_VTK_HPP
#define WRITE_VTK_HPP

#include <iosfwd>

#include "MBForward.hpp"
#include "MBWriterIface.hpp"

class MBWriteUtilIface;

//class MB_DLL_EXPORT WriteVtk : public MBWriterIface
class WriteVtk : public MBWriterIface
{
 
public:

   //! Constructor
   WriteVtk(MBInterface *impl);

   //! Destructor
  virtual ~WriteVtk();
  
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

private:

    //! Get entities to write, given set list passed to \ref write_file
  MBErrorCode gather_mesh( const MBEntityHandle* set_list,
                           int num_sets, 
                           MBRange& nodes,
                           MBRange& elems );
    
    //! Write 4-line VTK file header
  MBErrorCode write_header( std::ostream& stream );
  
    //! Write node coordinates
  MBErrorCode write_nodes( std::ostream& stream, const MBRange& nodes );
  
    //! Write element connectivity
  MBErrorCode write_elems( std::ostream& stream, const MBRange& nodes, const MBRange& elems );
  
    //! Write all tags on either the list of nodes or the list of elements
  MBErrorCode write_tags( std::ostream& stream, bool nodes, const MBRange& entities,
                          const MBTag* tag_list, int num_tags );
  
    //! Write the tad description for the passed tag and call the template
    //! \ref write_tag function to write the tag data.
  MBErrorCode write_tag( std::ostream& stream, MBTag tag, const MBRange& entities, const MBRange& tagged_entities );
  
    //! Write tag data
  template <typename T> 
  MBErrorCode write_tag( std::ostream& stream, MBTag tag, const MBRange& entities, const MBRange& tagged_entities,
                         const int);

  MBErrorCode write_bit_tag( std::ostream& stream, MBTag tag, const MBRange& entities, const MBRange& tagged_entities );
    //! Write a list of values
  template <typename T>
  void write_data( std::ostream& stream, const std::vector<T>& data, unsigned vals_per_tag );

  MBInterface* mbImpl;
  MBWriteUtilIface* writeTool;
 
  bool mStrict; // If true, do not write data that cannot fit in strict VTK file format.
  
};

#endif
