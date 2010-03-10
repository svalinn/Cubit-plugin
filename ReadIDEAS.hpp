#include <iostream>
#include <fstream>
#include <vector>

#include "MBReaderIface.hpp"
#include "MBInterface.hpp"

#define VERTEX_LIST       2411 
#define MAKE_TETRAHEDRA   2412

#define MAT_PROP_TABLE_TAG "mat_prop_table"
#define PHYS_PROP_TABLE_TAG "phys_prop_table"

class MBReadUtilIface;

class ReadIDEAS : public MBReaderIface
{

public:

  static MBReaderIface* factory( MBInterface* );

  MBErrorCode load_file( const char* fname, 
			 const MBEntityHandle* meshset, 
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
  ReadIDEAS(MBInterface* impl = NULL);
  
  //! Destructor
  virtual ~ReadIDEAS() {}

protected:
  
  MBErrorCode skip_header();
  MBErrorCode create_vertices(MBEntityHandle& first_vertex, const MBTag* file_id_tag);
  MBErrorCode create_tetrahedral_elements(MBEntityHandle first_vertex, const MBTag* file_id_tag);
  
private:
  
  std::ifstream file;
  
  // Read mesh interface
  MBReadUtilIface* readMeshIface;
  
  // MOAB Interface
  MBInterface* MBI;

};
