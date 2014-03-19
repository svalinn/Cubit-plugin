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
// Filename      : ReadRTT.hpp                        
//                                
// Purpose       : RTT file reader
//                                             
// Creator       : Andrew Davis
//                                   
// Date          : 02/2014                
//                                                  
//-------------------------------------------------------------------------     
                                    
#ifndef READRTT_HPP                     
#define READRTT_HPP              
                                     
#ifndef IS_BUILDING_MB                   
  #error "ReadRTT.hpp isn't supposed to be included into an application"
#endif   

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

#include "moab/Interface.hpp"
#include "moab/ReaderIface.hpp"
#include "FileTokenizer.hpp"
#include "moab/RangeMap.hpp"

// structure to hold sense & vol data
struct boundary {
  int sense;
  std::string name;
};

// structure to hold side data
struct side {
  int id;
  int senses[2];
  std::string names[2];
};

// structure to hold cell data
struct cell {
  int id;
  std::string name;
};

// structure to hold node data
struct node {
  int id;
  double x,y,z;
};

// strucutre to hold facet data
struct facet {
  int id;
  int connectivity[3];
  int from;
  int to;
};

// strucutre to hold tet data
struct tet {
  int id;
  int connectivity[4];
  int material_number;
};

namespace moab {

class ReadUtilIface;
class GeomTopoTool;

class ReadRTT : public ReaderIface
{

public:
  // factory method
  static ReaderIface* factory( Interface* );
  
  ErrorCode load_file( const char* file_name,
                       const EntityHandle* file_set,
                       const FileOptions& opts,
                       const SubsetList* subset_list = 0,
                       const Tag* file_id_tag = 0 );
  // constructor
  ReadRTT(Interface* impl = NULL);

  // destructor
  virtual ~ReadRTT();

  ErrorCode read_tag_values( const char* file_name,
                             const char* tag_name,
                             const FileOptions& opts,
                             std::vector<int>& tag_values_out,
                             const SubsetList* subset_list = 0 );

protected:
  
private:  
  // read mesh interface
  ReadUtilIface* readMeshIface;
  // Moab Interface
  Interface* MBI;
  // geom tool instance 
  GeomTopoTool* myGeomTool;

  ErrorCode setup_basic_tags();

  ErrorCode generate_topology(std::vector<side> side_data,
			      std::vector<cell> cell_data,
			      std::map <int,EntityHandle> &surface_map);

  void generate_parent_child_links(int num_ents[4],std::vector<EntityHandle> entity_map[4],
				   std::vector<side> side_data, std::vector<cell> cell_data);
  void set_surface_senses(int num_ents[4], std::vector<EntityHandle> entity_map[4],
			  std::vector<side> side_data, std::vector<cell> cell_data);

  /**
   * creates the group data requried for dagmc, reflecting planes, material assignments etc
   * @param num_ents, vector of the number of entities in each dimension
   * @param entity_map, vector of vector of entitiy handles for each dimension
   * @param side_data, vector of side data
   * @param cell_data, vector of cell data
   */
  ErrorCode setup_group_data(int num_ents[4], std::vector<EntityHandle> entity_map[4],
			std::vector<side> side_data, std::vector<cell> cell_data);


  /**
   * create a group of a given name, mustkeep track of id
   * @param group_name, name of the group
   * @param id, integer id number
   * returns the entity handle of the group
   */
  EntityHandle create_group(std::string group_name, int id);


  ErrorCode build_moab(std::vector<node> node_data,
		       std::vector<facet> facet_data,
		       std::vector<tet> tet_data,
		       std::map <int,EntityHandle> surface_map);

  ErrorCode read_sides(const char* filename, std::vector<side> &side_data);
  ErrorCode read_cells(const char* filename, std::vector<cell> &cell_data);
  ErrorCode read_nodes(const char* filename, std::vector<node> &node_data);
  ErrorCode read_facets(const char* filename, std::vector<facet> &facet_data);
  ErrorCode read_tets(const char* filename, std::vector<tet> &tet_data);

  cell get_cell_data(std::string celldata);
  side get_side_data(std::string sidedata);
  node get_node_data(std::string nodedata);
  facet get_facet_data(std::string facetdata);
  tet get_tet_data(std::string tetdata);

  std::vector<std::string> split_string(std::string string_to_split, char split_char);
  boundary split_name(std::string atilla_cellname);

  /**
   * Count the number of unique surface numbers in the dataset, also get list of surface numbers
   * @param side_data, collection of all the side data in the mesh
   * @param surface_numbers, collection of surface numbers
   * returns the number of surface numbers
   */
  int count_sides(std::vector<side> side_data,
		  std::vector<int> &surface_numbers);


  Tag geom_tag,id_tag,name_tag,category_tag,faceting_tol_tag;
};

} // namespace moab

#endif
