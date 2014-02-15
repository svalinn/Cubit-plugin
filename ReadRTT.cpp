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
 * \class ReadRTT
 * \brief ReadRTT based on ReadNASTRAN
 *
 * See: 
 *
 * \author Andrew Davis
 */



#include "ReadRTT.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <assert.h>
#include <cmath>

#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "Internals.hpp" // for MB_START_ID
#include "moab/Range.hpp"
#include "moab/FileOptions.hpp"
#include "FileTokenizer.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"

namespace moab {

ReaderIface* ReadRTT::factory( Interface* iface ) { 
  return new ReadRTT( iface );
}

// constructor
ReadRTT::ReadRTT(Interface* impl)
  : MBI(impl) {
    assert(NULL != impl);
    MBI->query_interface(readMeshIface);
    assert(NULL != readMeshIface);
}

// destructor
ReadRTT::~ReadRTT() {
  if (readMeshIface) {
    MBI->release_interface(readMeshIface);
    readMeshIface = 0;
  }
}

ErrorCode ReadRTT::read_tag_values( const char*        /*file_name*/,
                                    const char*        /*tag_name*/,
                                    const FileOptions& /*opts*/,
                                    std::vector<int>&  /*tag_values_out*/,
                                    const SubsetList*  /*subset_list*/ )
{
  return MB_NOT_IMPLEMENTED;
}

// load the file as called by the Interface function
ErrorCode ReadRTT::load_file(const char                      *filename, 
                             const EntityHandle            *, 
                             const FileOptions             &,
                             const ReaderIface::SubsetList *subset_list,
                             const Tag*                     file_id_tag) {
  std::cout << "Reading an RTT file" << std::endl;
  // at this time there is no support for reading a subset of the file
  if (subset_list) {
    readMeshIface->report_error( "Reading subset of files not supported for RTT." );
    return MB_UNSUPPORTED_OPERATION;
  }
  std::cout << "Reading an RTT file" << std::endl;
  std::vector<node> node_data;
  std::cout << "Reading node data..." << std::endl;
  ErrorCode rval = ReadRTT::read_nodes(filename,node_data);
  std::cout << node_data.size() << std::endl;

  std::cout << "Reading facet data..." << std::endl;
  std::vector<facet> facet_data;
  rval = ReadRTT::read_facets(filename,facet_data);
  std::cout << facet_data.size() << std::endl;

  std::cout << "Reading tet data..." << std::endl;
  std::vector<tet> tet_data;
  rval = ReadRTT::read_tets(filename,tet_data);
  std::cout << tet_data.size() << std::endl;

  rval = ReadRTT::build_moab(node_data,facet_data,tet_data);

  return MB_SUCCESS;
}

/*
 * builds the moab representation of the mesh
 */
ErrorCode ReadRTT::build_moab(std::vector<node> node_data,
			      std::vector<facet> facet_data,
			      std::vector<tet> tet_data ){
  ErrorCode rval;
  EntityHandle file_set;
  rval = MBI->create_meshset( MESHSET_SET, file_set );
  if (MB_SUCCESS != rval)
    return rval;

  // create the vertices
  EntityHandle handle;
  std::vector<node>::iterator it;
  Range mb_coords;
  for ( it = node_data.begin() ; it != node_data.end() ; ++it) {
    node tmp = *it;
    double coords[3]={tmp.x,tmp.y,tmp.z};
    rval = MBI->create_vertex(coords,handle);
    mb_coords.insert(handle);
  }
  // add verts to set
  rval = MBI->add_entities(file_set,mb_coords);

  // create the facets
  EntityHandle triangle;
  std::vector<facet>::iterator it_f;
  Range mb_tris;
  for ( it_f = facet_data.begin() ; it_f != facet_data.end() ; ++it_f) {
    facet tmp = *it_f;
    EntityHandle tri_nodes[3]={mb_coords[tmp.connectivity[0]-1],
			       mb_coords[tmp.connectivity[1]-1],
			       mb_coords[tmp.connectivity[2]-1]};
    rval = MBI->create_element(MBTRI,tri_nodes,3,triangle);
    mb_tris.insert(triangle);
  }
  // add tris to set
  rval = MBI->add_entities(file_set,mb_tris);

  // create material number tag
  Tag mat_num_tag;
  //  int zero = 0;
  rval = MBI->tag_get_handle( "MATERIAL_NUMBER", 1, MB_TYPE_INTEGER,
			      mat_num_tag, MB_TAG_SPARSE|MB_TAG_CREAT); 

  // create the tets
  EntityHandle tetra;
  std::vector<tet>::iterator it_t;
  Range mb_tets;
  for ( it_t = tet_data.begin() ; it_t != tet_data.end() ; ++it_t) {
    tet tmp = *it_t;
    EntityHandle tet_nodes[4]={mb_coords[tmp.connectivity[0]-1],
			       mb_coords[tmp.connectivity[1]-1],
			       mb_coords[tmp.connectivity[2]-1],
			       mb_coords[tmp.connectivity[3]-1]};
    rval = MBI->create_element(MBTET,tet_nodes,4,tetra);
    int mat_number = tmp.material_number;
    rval = MBI->tag_set_data(mat_num_tag,&tetra,1,&mat_number);
    // set the tag data
    mb_tets.insert(tetra);
  }
  // add tris to set
  rval = MBI->add_entities(file_set,mb_tets);
  
  return MB_SUCCESS;
}

/*
 * Reads the node data fromt the filename pointed to
 */ 
ErrorCode ReadRTT::read_nodes(const char* filename, std::vector<node> &node_data ){
  std::string line; // the current line being read
  std::ifstream input_file (filename); // filestream for rttfile
  // file ok?
  if ( !input_file.good() )
    {
      std::cout << "Problems reading file = " << filename << std::endl;
    }
  // if it works
  if (input_file.is_open())
    {
      while ( std::getline (input_file,line) )
	{
	  if(line.compare("nodes\0") == 0)
	    {
	      // read lines until find end nodes
	      while( std::getline( input_file,line)) 
		{   
		  if(line.compare("end_nodes\0") == 0)
		    break;
		  node data = ReadRTT::get_node_data(line);
		  node_data.push_back(data);
		  //		  std::cout << data.id << " " << data.x << " " << " " << data.y << " " << data.z << std::endl;
		  //		  std::cout << line << std::endl;
		}
	    }
	}
      input_file.close();
    }
  return MB_SUCCESS;
}

/*
 * Reads the facet data fromt the filename pointed to
 */ 
ErrorCode ReadRTT::read_facets(const char* filename, std::vector<facet> &facet_data ){
  std::string line; // the current line being read
  std::ifstream input_file (filename); // filestream for rttfile
  // file ok?
  if ( !input_file.good() )
    {
      std::cout << "Problems reading file = " << filename << std::endl;
    }
  // if it works
  if (input_file.is_open())
    {
      while ( std::getline (input_file,line) )
	{
	  if(line.compare("sides\0") == 0)
	    {
	      // read lines until find end nodes
	      while( std::getline( input_file,line)) 
		{   
		  if(line.compare("end_sides\0") == 0)
		    break;
		  facet data = ReadRTT::get_facet_data(line);
		  facet_data.push_back(data);
		  //		  std::cout << data.id << " " << data.x << " " << " " << data.y << " " << data.z << std::endl;
		  //		  std::cout << line << std::endl;
		}
	    }
	}
      input_file.close();
    }
  return MB_SUCCESS;
}

/*
 * Reads the facet data fromt the filename pointed to
 */ 
ErrorCode ReadRTT::read_tets(const char* filename, std::vector<tet> &tet_data ){
  std::string line; // the current line being read
  std::ifstream input_file (filename); // filestream for rttfile
  // file ok?
  if ( !input_file.good() )
    {
      std::cout << "Problems reading file = " << filename << std::endl;
    }
  // if it works
  if (input_file.is_open())
    {
      while ( std::getline (input_file,line) )
	{
	  if(line.compare("cells\0") == 0)
	    {
	      // read lines until find end nodes
	      while( std::getline( input_file,line)) 
		{   
		  if(line.compare("end_cells\0") == 0)
		    break;
		  tet data = ReadRTT::get_tet_data(line);
		  tet_data.push_back(data);
		  //		  std::cout << data.id << " " << data.x << " " << " " << data.y << " " << data.z << std::endl;
		  //		  std::cout << line << std::endl;
		}
	    }
	}
      input_file.close();
    }
  return MB_SUCCESS;
}

/*
 * given the string nodedata, get the id number and coordinates of the node
 */
node ReadRTT::get_node_data(std::string nodedata) {
  node new_node;
  std::vector<std::string> tokens;
  tokens = ReadRTT::split_string(nodedata);

  std::vector<std::string>::iterator it;
  new_node.id=std::atoi(tokens[1].c_str());
  new_node.x=std::atof(tokens[2].c_str());
  new_node.y=std::atof(tokens[3].c_str());
  new_node.z=std::atof(tokens[4].c_str());
  return new_node;
}

/*
 * given the string nodedata, get the id number, connectivity and sense data
 */
facet ReadRTT::get_facet_data(std::string facetdata) {
  facet new_facet;
  std::vector<std::string> tokens;
  tokens = ReadRTT::split_string(facetdata);
  
  /*
  std::vector<std::string>::iterator it;
  for ( it = tokens.begin() ; it != tokens.end() ; ++it )
    std::cout << "a"<<*it<<"a"<<std::endl;
  exit(1);
  */

  new_facet.id = std::atoi(tokens[1].c_str());
  new_facet.connectivity[0] = std::atoi(tokens[3].c_str());
  new_facet.connectivity[1] = std::atoi(tokens[4].c_str());
  new_facet.connectivity[2] = std::atoi(tokens[5].c_str());
  new_facet.from = std::atoi(tokens[6].c_str());
  new_facet.to = std::atoi(tokens[7].c_str());

  //  new_node.id=std::atoi(tokens[1].c_str());
  //  new_node.x=std::atof(tokens[2].c_str());
  //  new_node.y=std::atof(tokens[3].c_str());
  //  new_node.z=std::atof(tokens[4].c_str());
  return new_facet;
}

/*
 * given the string tetdata, get the id number, connectivity and mat num of the tet
 */
tet ReadRTT::get_tet_data(std::string tetdata) {
  tet new_tet;
  std::vector<std::string> tokens;
  tokens = ReadRTT::split_string(tetdata);
  

  /*
  std::vector<std::string>::iterator it;
  for ( it = tokens.begin() ; it != tokens.end() ; ++it )
    std::cout << "a"<<*it<<"a"<<std::endl;
  exit(1);
  */

  new_tet.id = std::atoi(tokens[1].c_str());
  new_tet.connectivity[0] = std::atoi(tokens[3].c_str());
  new_tet.connectivity[1] = std::atoi(tokens[4].c_str());
  new_tet.connectivity[2] = std::atoi(tokens[5].c_str());
  new_tet.connectivity[3] = std::atoi(tokens[6].c_str());
  new_tet.material_number = std::atoi(tokens[7].c_str());
  
  //  new_node.id=std::atoi(tokens[1].c_str());
  //  new_node.x=std::atof(tokens[2].c_str());
  //  new_node.y=std::atof(tokens[3].c_str());
  //  new_node.z=std::atof(tokens[4].c_str());
  return new_tet;
}

/*
 * splits a string in a vector of strings split by spaces
 */
std::vector<std::string> ReadRTT::split_string(std::string string_to_split) {
  std::istringstream ss( string_to_split );
  std::vector<std::string> tokens;
  while (!ss.eof())         // See the WARNING above for WHY we're doing this!
    {
      std::string x;               // here's a nice, empty string
      std::getline( ss, x, ' ');  // try to read the next field into it
      //      std::cout << x << std::endl;      // print it out, EVEN IF WE ALREADY HIT EOF
      tokens.push_back(x);
    }

  // remove empty tokens
  std::vector<std::string>::iterator it;
  for ( it = tokens.begin() ; it != tokens.end() ; ++it ) {
    std::string string = *it;
    if(string.compare("\0") == 0 )
    //    if(string.find("\0") != std::string::npos )
      tokens.erase(it);
  }
  return tokens;
}

} // namespace moab
