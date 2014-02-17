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
#include <map>
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

  ReaderIface* ReadRTT::factory( Interface* iface ){
  return new ReadRTT( iface );
}

// constructor
ReadRTT::ReadRTT(Interface* impl)
  : MBI(impl),geom_tag(0), id_tag(0), name_tag(0), category_tag(0), faceting_tol_tag(0) {
    assert(NULL != impl);
    MBI->query_interface(readMeshIface);
    assert(NULL != readMeshIface);
    
    // this section copied from ReadCGM initalisation
    int negone = -1, zero = 0;
    ErrorCode rval;
    rval = MBI->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER,
				geom_tag, MB_TAG_SPARSE|MB_TAG_CREAT, &negone);
    assert(!rval);
    rval = MBI->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
				id_tag, MB_TAG_DENSE|MB_TAG_CREAT, &zero);
    assert(!rval);
    rval = MBI->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE,
				name_tag, MB_TAG_SPARSE|MB_TAG_CREAT );
    assert(!rval);
    rval = MBI->tag_get_handle( CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, MB_TYPE_OPAQUE,
				category_tag, MB_TAG_SPARSE|MB_TAG_CREAT );
    assert(!rval);
    rval = MBI->tag_get_handle("FACETING_TOL", 1, MB_TYPE_DOUBLE, faceting_tol_tag,
			     MB_TAG_SPARSE|MB_TAG_CREAT, &zero );
    assert(!rval);
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
  ErrorCode rval;

  // at this time there is no support for reading a subset of the file
  if (subset_list) {
    readMeshIface->report_error( "Reading subset of files not supported for RTT." );
    return MB_UNSUPPORTED_OPERATION;
  }
  std::cout << "Reading an RTT file" << std::endl;
  
  // read the side_flag data
  std::vector<side> side_data;
  std::cout << "Reading the side data..." << std::endl;
  rval = ReadRTT::read_sides(filename,side_data);
  std::cout << side_data.size() << std::endl;

  // read the cell data
  std::vector<cell> cell_data;
  std::cout << "Reading the cell data..." << std::endl;
  rval = ReadRTT::read_cells(filename,cell_data);
  std::cout << cell_data.size() << std::endl;

  // read the node data
  std::vector<node> node_data;
  std::cout << "Reading node data..." << std::endl;
  rval = ReadRTT::read_nodes(filename,node_data);
  std::cout << node_data.size() << std::endl;

  std::cout << "Reading facet data..." << std::endl;
  std::vector<facet> facet_data;
  rval = ReadRTT::read_facets(filename,facet_data);
  std::cout << facet_data.size() << std::endl;

  std::cout << "Reading tet data..." << std::endl;
  std::vector<tet> tet_data;
  rval = ReadRTT::read_tets(filename,tet_data);
  std::cout << tet_data.size() << std::endl;

  std::cout << "Generate topology ..." << std::endl;
  rval = ReadRTT::generate_topology(side_data,cell_data);
  
  std::cout << "Generate MOAB Data ..." << std::endl;
  rval = ReadRTT::build_moab(node_data,facet_data,tet_data);

  return MB_SUCCESS;
}

/*
 * builds the topology of the problem
 */
ErrorCode ReadRTT::generate_topology(std::vector<side> side_data,
				     std::vector<cell> cell_data){
  std::vector<EntityHandle> entmap[4]; // one for each dimension
  ErrorCode rval;

  const char geom_categories[][CATEGORY_TAG_SIZE] =
    {"Vertex\0", "Curve\0", "Surface\0", "Volume\0", "Group\0"};
  const char* const names[] = { "Vertex", "Curve", "Surface", "Volume"};

  //  rval = setup_basic_tags();

  // loop over surfaces
  int dim = 2;
  for ( unsigned int i = 0 ; i != side_data.size() ; i++ ) {
    EntityHandle handle;
    rval = MBI->create_meshset( dim == 1 ? MESHSET_ORDERED : MESHSET_SET, handle );
    if (rval != MB_SUCCESS )
      return rval;
    // collect the entity handles
    entmap[dim].push_back(handle);

    rval = MBI->tag_set_data( geom_tag, &handle, 1, &dim );
    if (MB_SUCCESS != rval)     
      return rval;              
    int id = side_data[i].id;         
    rval = MBI->tag_set_data( id_tag, &handle, 1, &id );
    if (MB_SUCCESS != rval)     
      return rval;              
    rval = MBI->tag_set_data( category_tag, &handle, 1, &geom_categories[dim] );
    if (MB_SUCCESS != rval)
      return rval;
  }

  // loop over volumes
  dim = 3;
  for ( unsigned int i = 0 ; i != cell_data.size() ; i++ ) {
    EntityHandle handle;
    rval = MBI->create_meshset( dim == 1 ? MESHSET_ORDERED : MESHSET_SET, handle );
    if (rval != MB_SUCCESS )
      return rval;
    // collect the entity handles
    entmap[dim].push_back(handle);

    rval = MBI->tag_set_data( geom_tag, &handle, 1, &dim );
    if (MB_SUCCESS != rval)     
      return rval;              
    int id = side_data[i].id;         
    rval = MBI->tag_set_data( id_tag, &handle, 1, &id );
    if (MB_SUCCESS != rval)     
      return rval;              
    rval = MBI->tag_set_data( category_tag, &handle, 1, &geom_categories[dim] );
    if (MB_SUCCESS != rval)
      return rval;
  }

  // generate parent child links

  return MB_SUCCESS;
}

/*
ErrorCode setup_basic_tags(){
  int negone = -1, zero = 0;
  ErrorCode rval;
  rval = MBI->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER,
                                  geom_tag, MB_TAG_SPARSE|MB_TAG_CREAT, &negone);
  assert(!rval);
  rval = MBI->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
                                  id_tag, MB_TAG_DENSE|MB_TAG_CREAT, &zero);
  assert(!rval);
  rval = MBI->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE,
                                  name_tag, MB_TAG_SPARSE|MB_TAG_CREAT );
  assert(!rval);
  rval = MBI->tag_get_handle( CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, MB_TYPE_OPAQUE,
                                  category_tag, MB_TAG_SPARSE|MB_TAG_CREAT );
  assert(!rval);
  rval = MBI->tag_get_handle("FACETING_TOL", 1, MB_TYPE_DOUBLE, faceting_tol_tag,
			     MB_TAG_SPARSE|MB_TAG_CREAT, &zero );
  assert(!rval);

  return rval;
}
*/

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


  // create sense tag
  Tag side_id_tag, surface_number_tag;
  //  int zero = 0;
  rval = MBI->tag_get_handle( "SIDEID_TAG", 1, MB_TYPE_INTEGER,
			      side_id_tag, MB_TAG_SPARSE|MB_TAG_CREAT); 
  rval = MBI->tag_get_handle( "SURFACE_NUMBER", 1, MB_TYPE_INTEGER,
			      surface_number_tag, MB_TAG_SPARSE|MB_TAG_CREAT); 

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
    // tag in sense
    rval = MBI->tag_set_data(side_id_tag,&triangle,1,&tmp.from);
    // tag out sense
    rval = MBI->tag_set_data(surface_number_tag,&triangle,1,&tmp.to);

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
 * reads the side data from the filename pointed to
 */
ErrorCode ReadRTT::read_sides(const char* filename, std::vector<side> &side_data ){
  std::string line; // the current line being read
  std::ifstream input_file (filename); // filestream for rttfile
  // file ok?
  if ( !input_file.good() )
    {
      std::cout << "Problems reading file = " << filename << std::endl;
      return MB_FAILURE;
    }
  // if it works
  if (input_file.is_open())
    {
      while ( std::getline (input_file,line) )
	{
	  if(line.compare("  2 FACES\0") == 0)
	    {
	      // read lines until find end nodes
	      while( std::getline( input_file,line)) 
		{   
		  if(line.compare("end_side_flags\0") == 0)
		    break;
		  side data = ReadRTT::get_side_data(line);
		  side_data.push_back(data);
		}
	    }
	}
      input_file.close();
    }
  return MB_SUCCESS;
}

/*
 * reads the cell data from the filename pointed to
 */
  ErrorCode ReadRTT::read_cells(const char* filename, std::vector<cell> &cell_data ){
  std::string line; // the current line being read
  std::ifstream input_file (filename); // filestream for rttfile
  // file ok?
  if ( !input_file.good() )
    {
      std::cout << "Problems reading file = " << filename << std::endl;
      return MB_FAILURE;
    }
  // if it works
  if (input_file.is_open())
    {
      while ( std::getline (input_file,line) )
	{
	  if(line.compare("  1 REGIONS\0") == 0)
	    {
	      // read lines until find end nodes
	      while( std::getline( input_file,line)) 
		{   
		  if(line.compare("end_cell_flags\0") == 0)
		    break;
		  cell data = ReadRTT::get_cell_data(line);
		  cell_data.push_back(data);
		}
	    }
	}
      input_file.close();
    }
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
 * given the string sidedata, get the id number, senses and names of the sides
 */
side ReadRTT::get_side_data(std::string sidedata) {
  side new_side;
  std::vector<std::string> tokens;
  tokens = ReadRTT::split_string(sidedata,' ');

  std::vector<std::string>::iterator it;
  //  for ( it = tokens.begin() ; it != tokens.end() ; ++it )
  //    std::cout << *it << std::endl;

  new_side.id = std::atoi(tokens[2].c_str());

  std::vector<std::string> cell_names = ReadRTT::split_string(tokens[3],'/');
  boundary new_bnd = ReadRTT::split_name(cell_names[0]);
  new_side.senses[0]=new_bnd.sense;
  new_side.names[0]=new_bnd.name;
  if (cell_names.size() > 1 )
    {
      boundary new_bnd = ReadRTT::split_name(cell_names[1]);
      new_side.senses[1]=new_bnd.sense;
      new_side.names[1]=new_bnd.name;
    }
  else
    {
      new_side.senses[1]=0;
      new_side.names[1]="\0";
    }

  return new_side;
}

/*
 * given the string celldata, get the id number and name of each cell
 */
cell ReadRTT::get_cell_data(std::string celldata) {
  cell new_cell;
  std::vector<std::string> tokens;
  tokens = ReadRTT::split_string(celldata,' ');

  std::vector<std::string>::iterator it;

  new_cell.id = std::atoi(tokens[2].c_str());
  new_cell.name = tokens[3];

  return new_cell;
}

/*
 * given the string nodedata, get the id number and coordinates of the node
 */
node ReadRTT::get_node_data(std::string nodedata) {
  node new_node;
  std::vector<std::string> tokens;
  tokens = ReadRTT::split_string(nodedata,' ');

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
  tokens = ReadRTT::split_string(facetdata,' ');
  
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
  tokens = ReadRTT::split_string(tetdata,' ');
  

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
 * splits string into sense and name, to later facilitate the building
 * of sense data, strips off the tailing @ if it exists
 */
boundary ReadRTT::split_name(std::string atilla_cellname) {
  boundary new_boundary;
  // default initialisation
  new_boundary.sense = 0;
  new_boundary.name = "\0";

  if( atilla_cellname.find("+") != std::string::npos )
    {			
      new_boundary.sense = 1;
      // look for the @# we do not want it
      std::size_t found = atilla_cellname.find("@");
      if(found != std::string::npos)
	new_boundary.name=atilla_cellname.substr(3,found);
      else
	new_boundary.name=atilla_cellname.substr(3,atilla_cellname.length());
    }
  else if ( atilla_cellname.find("-") != std::string::npos )
    {			
      new_boundary.sense = -1;
      new_boundary.name=atilla_cellname.substr(3,atilla_cellname.length());
    }
  return new_boundary;
}

/*
 * splits a string in a vector of strings split by spaces
 */
  std::vector<std::string> ReadRTT::split_string(std::string string_to_split, char split_char) {
  std::istringstream ss( string_to_split );
  std::vector<std::string> tokens;
  while (!ss.eof())         // See the WARNING above for WHY we're doing this!
    {
      std::string x;               // here's a nice, empty string
      std::getline( ss, x, split_char);  // try to read the next field into it
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
