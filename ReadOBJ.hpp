#ifndef READ_OBJ_HPP
#define READ_OBJ_HPP     
                                     
#ifndef IS_BUILDING_MB                   
  #error "ReadOBJ.hpp isn't supposed to be included into an application"
#endif   

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "moab/Interface.hpp"
#include "moab/ReaderIface.hpp"
#include "FileTokenizer.hpp"
#include "moab/RangeMap.hpp"
#include "MBTagConventions.hpp"

/* struct vertex is a structure that stores coordinates
   of vertices. This is a convenience structure making 
   i/o easier.
*/
struct vertex {
  int vertex_id;
  double coord[3];
};

/* struct face is a structure that stores connectivity.
   This is a convenience structure makin i/o easier.
*/
struct face {
  int face_id;
  moab::EntityHandle conn[3];
};
namespace moab {

class ReadUtilIface;
class GeomTopoTool;

class ReadOBJ : public ReaderIface
{
   
public:

  //! factory method 
  static ReaderIface* factory( Interface* );

  ErrorCode load_file( const char* file_name,
                       const EntityHandle* file_set,
                       const FileOptions& opts,
                       const SubsetList* subset_list = 0,
                       const Tag* file_id_tag = 0 );

  ErrorCode read_tag_values( const char* file_name,
                             const char* tag_name,
                             const FileOptions& opts,
                             std::vector<int>& tag_values_out,
                             const SubsetList* subset_list = 0 );
 

 
  //! Constructor
  ReadOBJ(Interface* impl = NULL);

   //! Destructor
  virtual ~ReadOBJ();

private:
  ReadUtilIface* readMeshIface;

  //! interface instance
  Interface* MBI;

  GeomTopoTool* myGeomTool;

  const char *fileName;
  Tag geom_tag,id_tag,name_tag,category_tag,faceting_tol_tag, geometry_resabs_tag, obj_name_tag, 
    sense_tag;
  

  /* The tokenize function takes a string as input and splits it into
   * a vector of strings based on the delimiter
   */
  static const char* delimiters;

  void tokenize( const std::string& str, std::vector<std::string>& tokens,
                 const char* delimiters );
  
  /*
   * The create_object funtion will create a new meshset for
   * each object that contains all vertices and faces 
   */
  ErrorCode create_new_object(std::string object_name, 
                              int object_id,
                              EntityHandle &curr_obj_meshset);
  
  /* create_new_vertex converts tokenized string input to 
     vertex structure
   */
  ErrorCode create_new_vertex (std::vector<std::string> v_tokens,
                                      EntityHandle &vertex_eh); 

  /* create_new_face converts tokenized string input to 
   * face structure
   */
  ErrorCode create_new_face (std::vector<std::string> f_tokens,
                                       const std::vector<EntityHandle>&vertex_list,
                                       EntityHandle &face_eh); 
 

  /*
   * The split_quad function creates 1 new vertex and 4 new tri faces
   * from a quad face.  
   */

  ErrorCode split_quad(std::vector<std::string> f_tokens,
                                       std::vector<EntityHandle>&vertex_list,
                                       EntityHandle &new_vertex_eh,
                                       Range &face_eh);

  ErrorCode create_center_vertex( Range quad_vert_eh, 
                                EntityHandle &new_vertex_eh);

  ErrorCode create_tri_faces( Range quad_vert_eh, 
                                     EntityHandle center_vertex_eh,
                                     Range &face_eh );

};

} // namespace moab



#endif
