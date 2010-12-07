//-------------------------------------------------------------------------
// Filename      : ReadNC.hpp
//
// Purpose       : Climate NC file reader
//
// Creator       : Tim Tautges
//-------------------------------------------------------------------------

#ifndef READNC_HPP
#define READNC_HPP

#ifndef IS_BUILDING_MB
//#error "ReadNC.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <map>
#include <string>
#include <netcdf.h>

#include "moab/Forward.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/Range.hpp"
#include "DebugOutput.hpp"

namespace moab {

class ReadUtilIface;


//! Output Exodus File for VERDE
class ReadNC : public ReaderIface
{
   
public:
  
  static ReaderIface* factory( Interface* );
  
    //! load an NC file
  ErrorCode load_file( const char* file_name,
                       const EntityHandle* file_set,
                       const FileOptions& opts,
                       const SubsetList* subset_list = 0,
                       const Tag* file_id_tag = 0 );

   //! Constructor
   ReadNC(Interface* impl = NULL);

   //! Destructor
  virtual ~ReadNC();

    /**
     *\brief Read tag values from a file.
     *
     * Read the list if all integer tag values from the file for
     * a tag that is a single integer value per entity.
     *
     *\param file_name      The file to read.
     *\param tag_name       The tag for which to read values
     *\param tag_values_out Output: The list of tag values.
     *\param subset_list    An array of tag name and value sets specifying
     *                      the subset of the file to read.  If multiple
     *                      tags are specified, the sets that match all
     *                      tags (intersection) should be read.
     *\param subset_list_length The length of the 'subset_list' array.
     */
  virtual ErrorCode read_tag_values( const char* file_name,
                                     const char* tag_name,
                                     const FileOptions& opts,
                                     std::vector<int>& tag_values_out,
                                     const SubsetList* subset_list = 0 );

private:

  struct AttData 
  {
    int attId;
    size_t attLen;
    int attVarId;
    nc_type attDataType;
    std::string attName;
  };

  struct VarData 
  {
    int varId;
    int numAtts;
    nc_type varDataType;
    std::vector<int> varDims;
    std::map<std::string,AttData> varAtts;
    std::string varName;
  };

  ReadUtilIface* readMeshIface;

  bool dimension_exists(const char *attrib_name);
  
  void reset();

    //! read the header information
  ErrorCode read_header();

    //! get all global attributes in the file
  ErrorCode get_attributes(int var_id, int num_atts, std::map<std::string,AttData> &atts,
                           const char *prefix="");
  
    //! get all dimensions in the file
  ErrorCode get_dimensions();

    //! get the variable names and other info defined for this file
  ErrorCode get_variables();
  
    //! parse min/max i/j/k in options, if any
  ErrorCode init_ijk_vals(const FileOptions &opts);

  ErrorCode read_coordinate(const char *var_name, int lmin, int lmax,
                            std::vector<double> &cvals);
  
    //! number of dimensions in this nc file
  int number_dimensions();

    //! create vertices for the file
  ErrorCode create_verts();

  //------------member variables ------------//

    //! interface instance
  Interface* mbImpl;
  
  int CPU_WORD_SIZE;
  int IO_WORD_SIZE;

    //! file name
  std::string fileName;

    //! file number assigned by netcdf
  int fileId;
  
    //! number of dimensions
  int numDims;

    //! dimensions
  std::map<std::string,int> dimVals;

    //! number of global attributes
  int numGAtts;

    //! global attribs
  std::map<std::string,AttData> globalAtts;
  
    //! number of variables
  int numVars;

    //! variable info
  std::map<std::string,VarData> varInfo;
  
    //! dimensions of grid in file
  int iMin, iMax, jMin, jMax, kMin, kMax;
  
    //! dimensions of my part of grid
  int ilMin, ilMax, jlMin, jlMax, klMin, klMax;

    //! the actual dimension and variable names used for the i/j/k dimensions
  std::string iName, jName, kName;

    //! values for i/j/k
  std::vector<double> ilVals, jlVals, klVals;
  
    //! number of the variable of unlimited dimension, if any
  int numUnLim;

    //! Meshset Handle for the mesh that is currently being read
  EntityHandle mCurrentMeshHandle;

    //! starting vertex handle for this read
  EntityHandle startVertex;

  //! Cached tags for reading.  Note that all these tags are defined when the
  //! core is initialized.
  Tag mGlobalIdTag;

  int max_line_length, max_str_length;

    //! range of entities in initial mesh, before this read
  Range initRange;

    //! offset of first vertex id
  int vertexOffset;

    //! debug stuff
  DebugOutput dbgOut;
};

// inline functions
inline int ReadNC::number_dimensions() 
{
   return numDims;
}

} // namespace moab

#endif




