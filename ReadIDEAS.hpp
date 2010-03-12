#include <iostream>
#include <fstream>
#include <vector>

#include "moab/ReaderIface.hpp"
#include "moab/Interface.hpp"

#define VERTEX_LIST       2411 
#define MAKE_TETRAHEDRA   2412

#define MAT_PROP_TABLE_TAG "mat_prop_table"
#define PHYS_PROP_TABLE_TAG "phys_prop_table"

namespace moab {

class ReadUtilIface;

class ReadIDEAS : public ReaderIface
{

public:

  static ReaderIface* factory( Interface* );

  ErrorCode load_file( const char* fname, 
			 const EntityHandle* meshset, 
			 const FileOptions&,
                         const ReaderIface::IDTag* subset_list = 0,
                         int subset_list_length = 0,
                         const Tag* file_id_tag = 0 );

  ErrorCode read_tag_values( const char* file_name,
                               const char* tag_name,
                               const FileOptions& opts,
                               std::vector<int>& tag_values_out,
                               const IDTag* subset_list = 0,
                               int subset_list_length = 0 );
  
  //! Constructor
  ReadIDEAS(Interface* impl = NULL);
  
  //! Destructor
  virtual ~ReadIDEAS() {}

protected:
  
  ErrorCode skip_header();
  ErrorCode create_vertices(EntityHandle& first_vertex, const Tag* file_id_tag);
  ErrorCode create_tetrahedral_elements(EntityHandle first_vertex, const Tag* file_id_tag);
  
private:
  
  std::ifstream file;
  
  // Read mesh interface
  ReadUtilIface* readMeshIface;
  
  // MOAB Interface
  Interface* MBI;

};

} // namespace moab
