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

#include "moab/Forward.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/Range.hpp"
#include "DebugOutput.hpp"

#ifdef USE_MPI
#  include "moab_mpi.h"
#  include "moab/ParallelComm.hpp"
#endif 

#ifdef PNETCDF_FILE
#  include "pnetcdf.h"
#  define NCFUNC(func) ncmpi_ ## func
#  define NCFUNCA(func) ncmpi_ ## func ## _all
//#  define NCASYNCH
#  ifdef NCASYNCH
#    define NCREQ , &requests[j]
#    define NCFUNCAG(func) ncmpi_iget ## func 
#    define NCWAIT
#  else
#    define NCREQ 
#    define NCFUNCAG(func) ncmpi_get ## func ## _all
#  endif
#  define NCDF_SIZE MPI_Offset
#  define NCDF_DIFF MPI_Offset
#else
#  include "netcdf.h"
#define NCREQ
#define NCGET get
#  define NCFUNC(func) nc_ ## func
#  define NCFUNCA(func) nc_ ## func
#  define NCFUNCAG(func) nc_get ## func
#  define NCDF_SIZE size_t
#  define NCDF_DIFF ptrdiff_t
#endif

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

  virtual ErrorCode read_tag_values( const char* file_name,
                                     const char* tag_name,
                                     const FileOptions& opts,
                                     std::vector<int>& tag_values_out,
                                     const SubsetList* subset_list = 0 );
private:

  class AttData 
  {
    public:
    AttData() : attId(-1), attLen(0), attVarId(-2) {}
    int attId;
    NCDF_SIZE attLen;
    int attVarId;
    nc_type attDataType;
    std::string attName;
  };

  class VarData 
  {
    public:
    VarData() : varId(-1), numAtts(-1), read(false) {}
    int varId;
    int numAtts;
    nc_type varDataType;
    std::vector<int> varDims;
    std::map<std::string,AttData> varAtts;
    std::string varName;
    bool read;
    std::vector<Tag> varTags;
    std::vector<void*> varDatas;
    std::vector<std::vector<NCDF_SIZE> > readDims, readCounts;
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
  ErrorCode init_ijkt_vals(const FileOptions &opts);

  int compute_partition_alljorkori(int np, int nr,
                                   int &ilMin, int &ilMax, int &jlMin, int &jlMax, 
                                   int &klMin, int &klMax);
  
  int compute_partition_alljkbal(int np, int nr,
                                 int &ilMin, int &ilMax, int &jlMin, int &jlMax, 
                                 int &klMin, int &klMax);
  
  int compute_partition_sqij(int np, int nr,
                             int &ilMin, int &ilMax, int &jlMin, int &jlMax, 
                             int &klMin, int &klMax);

  int compute_partition_sqjk(int np, int nr,
                             int &ilMin, int &ilMax, int &jlMin, int &jlMax, 
                             int &klMin, int &klMax);
  
  ErrorCode read_coordinate(const char *var_name, int lmin, int lmax,
                            std::vector<double> &cvals);
  
    //! number of dimensions in this nc file
  unsigned int number_dimensions();

    //! create vertices for the file
  ErrorCode create_verts_hexes(EntityHandle file_set, Range &hexes);

    //! check number of vertices and elements against what's already in file_set
  ErrorCode check_verts_hexes(EntityHandle file_set);
  
  ErrorCode parse_options(const FileOptions &opts,
                          std::vector<std::string> &var_names, 
                          std::vector<int> &tstep_nums,
                          std::vector<double> &tstep_vals,
                          bool &nomesh,
                          bool &novars,
                          std::string &partition_tag_name);
  
  ErrorCode read_variables(EntityHandle file_set, std::vector<std::string> &var_names,
                           std::vector<int> &tstep_nums);
  
  ErrorCode read_variable_allocate(std::vector<VarData> &vdatas,
                                   std::vector<int> &tstep_nums, 
                                   Range &verts);
  
  ErrorCode read_variable_setup(std::vector<std::string> &var_names,
                                std::vector<int> &tstep_nums, 
                                std::vector<VarData> &vdatas);
  
  ErrorCode convert_variable(EntityHandle file_set, VarData &var_data, int tstep_num);
    
  ErrorCode get_tag(VarData &var_data, int tstep_num, Tag &tagh);
  
    //! create nc conventional tags
  ErrorCode create_tags(EntityHandle file_set, const std::vector<int> &tstep_nums);

//------------member variables ------------//

    //! interface instance
  Interface* mbImpl;
  
  int CPU_WORD_SIZE;
  int IO_WORD_SIZE;

    //! file name
  std::string fileName;

    //! file number assigned by netcdf
  int fileId;
  
    //! dimensions
  std::vector<std::string> dimNames;
  std::vector<int> dimVals;
  std::string iName, jName, kName, tName;

    //! global attribs
  std::map<std::string,AttData> globalAtts;
  
    //! variable info
  std::map<std::string,VarData> varInfo;
  
    //! dimensions of grid in file
  int iMin, iMax, jMin, jMax, kMin, kMax, tMin, tMax;
  
    //! dimensions of my part of grid
  int ilMin, ilMax, jlMin, jlMax, klMin, klMax;

    //! values for i/j/k
  std::vector<double> ilVals, jlVals, klVals, tVals;

    //! dimension numbers for i, j, k, t
  int iDim, jDim, kDim, tDim;
  
    //! number of the dimension of unlimited dimension, if any
  int numUnLim;

    //! Meshset Handle for the mesh that is currently being read
  EntityHandle mCurrentMeshHandle;

    //! starting vertex and element handles for this read
  EntityHandle startVertex, startElem;

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

    //! are we reading in parallel?
  bool isParallel;

    //! partitioning method
  int partMethod;
  
#ifdef USE_MPI
  ParallelComm *myPcomm;
#endif
  
};

// inline functions
inline unsigned int ReadNC::number_dimensions() 
{
  return dimVals.size();
}

} // namespace moab

#endif




