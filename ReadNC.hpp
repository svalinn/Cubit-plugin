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
#include <set>
#include <string>

#include "moab/Forward.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/Range.hpp"
#include "moab/ScdInterface.hpp"
#include "DebugOutput.hpp"

#ifdef USE_MPI
#  include "moab_mpi.h"
#  include "moab/ParallelComm.hpp"
#endif 

#ifdef PNETCDF_FILE
#  include "pnetcdf.h"
#  define NCFUNC(func) ncmpi_ ## func
#  define NCFUNCA(func) ncmpi_ ## func ## _all
// keep it this way , introduce another macro, used so far only for ucd mesh
//#  define NCASYNCH
#  ifdef NCASYNCH
#    define NCREQ , &requests[j]
#    define NCFUNCAG(func) ncmpi_iget ## func 
#    define NCWAIT
#  else
#    define NCREQ2 , &requests[idxReq++]
#    define NCFUNCAG2(func) ncmpi_iget ## func
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
class ScdInterface;
class NCHelper;

//! Output Exodus File for VERDE
class ReadNC : public ReaderIface
{
  friend class NCHelper;
  friend class NCHelperEuler;
  friend class NCHelperFV;
  friend class NCHelperHOMME;
  friend class NCHelperMPAS;

public:

  static ReaderIface* factory(Interface*);

    //! load an NC file
  ErrorCode load_file(const char* file_name,
                       const EntityHandle* file_set,
                       const FileOptions& opts,
                       const SubsetList* subset_list = 0,
                       const Tag* file_id_tag = 0);

   //! Constructor
   ReadNC(Interface* impl = NULL);

   //! Destructor
  virtual ~ReadNC();

  virtual ErrorCode read_tag_values(const char* file_name,
                                    const char* tag_name,
                                    const FileOptions& opts,
                                    std::vector<int>& tag_values_out,
                                    const SubsetList* subset_list = 0);

  // ENTLOCNSEDGE for north/south edge
  // ENTLOCWEEDGE for west/east edge
  enum EntityLocation {ENTLOCVERT = 0, ENTLOCNSEDGE, ENTLOCEWEDGE, ENTLOCFACE, ENTLOCSET, ENTLOCEDGE, ENTLOCREGION};

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
    VarData() : varId(-1), numAtts(-1), read(false), entLoc(ENTLOCSET), numLev(1), sz(0), has_t(false) {}
    int varId;
    int numAtts;
    nc_type varDataType;
    std::vector<int> varDims; // the dimension indices making up this multi-dimensional variable
    std::map<std::string,AttData> varAtts;
    std::string varName;
    bool read;
    std::vector<Tag> varTags; // one tag for each timestep, varTags[t]
    std::vector<void*> varDatas;
    std::vector<std::vector<NCDF_SIZE> > readDims; // start value for this [t][dim]
    std::vector<std::vector<NCDF_SIZE> > readCounts; // number of data values for this [t][dim]
    int entLoc;
    int numLev;
    int sz;
    bool has_t;
  };

  ReadUtilIface* readMeshIface;

  bool dimension_exists(const char *attrib_name);

  void reset();

    //! read the header information
  ErrorCode read_header();

    //! get all global attributes in the file
  ErrorCode get_attributes(int var_id, int num_atts, std::map<std::string, AttData> &atts,
                           const char *prefix="");

    //! get all dimensions in the file
  ErrorCode get_dimensions(int file_id, std::vector<std::string> &dim_names, std::vector<int> &dim_vals);

    //! get the variable names and other info defined for this file
  ErrorCode get_variables();

  ErrorCode read_coordinate(const char *var_name, int lmin, int lmax,
                            std::vector<double> &cvals);

    //! number of dimensions in this nc file
  unsigned int number_dimensions();

    //! create vertices and quads for scd mesh
  ErrorCode create_scd_verts_quads(ScdInterface *scdi, EntityHandle file_set, Range &quads);

    //! make sure that localGid is properly initialized for ucd mesh
  ErrorCode check_ucd_localGid(EntityHandle file_set);

    //! check number of vertices and faces against what's already in file_set
  ErrorCode check_verts_faces(EntityHandle file_set);

  ErrorCode parse_options(const FileOptions &opts,
                          std::vector<std::string> &var_names, 
                          std::vector<int> &tstep_nums,
                          std::vector<double> &tstep_vals);

  ErrorCode read_variable_to_set_allocate(std::vector<VarData> &vdatas,
                                          std::vector<int> &tstep_nums);

  ErrorCode read_variable_to_set(EntityHandle file_set, std::vector<VarData> &vdatas,
				 std::vector<int> &tstep_nums, bool scd_mesh);

  ErrorCode read_variable_to_nonset(EntityHandle file_set, std::vector<VarData> &vdatas,
				    std::vector<int> &tstep_nums, bool scd_mesh);

#ifdef PNETCDF_FILE
  ErrorCode read_variable_to_nonset_async(EntityHandle file_set, std::vector<VarData> &vdatas,
              std::vector<int> &tstep_nums);
#endif

  ErrorCode read_variables(EntityHandle file_set, std::vector<std::string> &var_names,
                           std::vector<int> &tstep_nums, bool scd_mesh);

  ErrorCode read_variable_allocate(EntityHandle file_set, std::vector<VarData> &vdatas,
                                   std::vector<int> &tstep_nums, bool scd_mesh);

  ErrorCode read_variable_setup(std::vector<std::string> &var_names,
                                std::vector<int> &tstep_nums, 
                                std::vector<VarData> &vdatas,
                                std::vector<VarData> &vsetdatas,
                                bool scd_mesh);

  ErrorCode convert_variable(VarData &var_data, int tstep_num, bool scd_mesh);

  ErrorCode get_tag_to_set(VarData &var_data, int tstep_num, Tag &tagh);

  ErrorCode get_tag(VarData &var_data, int tstep_num, Tag &tagh, int num_lev);

    //! create nc conventional tags
  ErrorCode create_tags(ScdInterface *scdi, EntityHandle file_set, 
                        const std::vector<int> &tstep_nums);

    //! create a character string attString of attMap.  with '\0'
    //! terminating each attribute name, ';' separating the data type
    //! and value, and ';' separating one name/data type/value from
    //! the next'.  attLen stores the end postion for each name/data
    //! type/ value.
  ErrorCode create_attrib_string(const std::map<std::string, AttData>& attMap, 
				 std::string& attString,
				 std::vector<int>& attLen);

  //! create COORDS tag for quads coordinate
  ErrorCode create_quad_coordinate_tag(EntityHandle file_set);

  //! Init info for dimensions that don't have corresponding 
  //! coordinate variables - this info is used for creating tags
  void init_dims_with_no_cvars_info();

  ErrorCode load_BIL(std::string dir_name,
                     const EntityHandle* file_set,
                     const FileOptions& opts,
                     const Tag* file_id_tag);

  ErrorCode get_BIL_dir();

  bool BIL_mode_enabled(const char * file_name);

  template <typename T> ErrorCode kji_to_jik(size_t ni, size_t nj, size_t nk, void *dest, T *source) 
      {
        size_t nik = ni * nk, nij = ni * nj;
        T *tmp_data = reinterpret_cast<T*>(dest);
        for (std::size_t j = 0; j != nj; ++j)
          for (std::size_t i = 0; i != ni; ++i)
            for (std::size_t k = 0; k != nk; ++k) 
              tmp_data[j*nik + i*nk + k] = source[k*nij + j*ni + i];
        return MB_SUCCESS;
      }

  // this version takes as input the moab range, from which we actually need just the
  // size of each sequence, for a proper transpose of the data
  // we read one time step, one variable at a time, usually, so we will
  template <typename T> ErrorCode kji_to_jik_stride(size_t , size_t nj, size_t nk, void *dest, T *source)
      {
        std::size_t idxInSource = 0;// position of the start of the stride
        // for each subrange, we will transpose a matrix of size subrange*nj*nk (subrange takes
        //                                                                       the role of ni)
        T *tmp_data = reinterpret_cast<T*>(dest);
        for (
          Range::pair_iterator pair_iter = localGid.pair_begin();
          pair_iter != localGid.pair_end();
          pair_iter++)
        {
          std::size_t size_range = pair_iter->second - pair_iter->first + 1;
          std::size_t nik = size_range * nk, nij = size_range * nj;
          for (std::size_t j = 0; j != nj; ++j)
            for (std::size_t i = 0; i != size_range; ++i)
              for (std::size_t k = 0; k != nk; ++k)
                tmp_data[idxInSource + j*nik + i*nk + k] = source[idxInSource + k*nij + j*size_range + i];
          idxInSource += (size_range * nj * nk);
        }
        return MB_SUCCESS;
      }
//------------member variables ------------//

    //! interface instance
  Interface* mbImpl;

  int CPU_WORD_SIZE;
  int IO_WORD_SIZE;

    //! file name
  std::string fileName;

    //! file numbers assigned by netcdf
  int fileId, connectId;

    //! dimensions
  std::vector<std::string> dimNames;
  // these should be taken out when we fix the dummy var info things
  std::set<std::string> dummyVarNames;
  std::vector<int> dimVals;
  std::string iName, jName, kName, tName;
  std::string iCName, jCName;

    //! global attribs
  std::map<std::string,AttData> globalAtts;

    //! variable info
  std::map<std::string,VarData> varInfo;

    //! dimensions of grid in file
  int gDims[6], tMin, tMax;

    //! dimensions of my part of grid
  int lDims[6];

  //! center dimensions of grid in file
  int gCDims[6];

    //! center dimensions of my part of grid
  int lCDims[6];

    //! values for i/j/k
  std::vector<double> ilVals, jlVals, klVals, tVals;

  //! center values for i/j
  std::vector<double> ilCVals, jlCVals;

    //! dimension numbers for i, j, k, t
  int iDim, jDim, kDim, tDim;

    //! center dimension numbers for i, j
  int iCDim, jCDim;

    //! number of the dimension of unlimited dimension, if any
  int numUnLim;

    //! Meshset Handle for the mesh that is currently being read
  EntityHandle mCurrentMeshHandle;

    //! starting vertex and element handles for this read
  EntityHandle startVertex, startElem;

  //! Cached tags for reading.  Note that all these tags are defined when the
  //! core is initialized.
  Tag mGlobalIdTag;

  // this is a pointer to the file id tag that is passed from ReadParallel
  // it gets deleted at the end of resolve sharing, but it will have same data
  // as the global id tag
  // global id tag is preserved, and is needed later on.
  const Tag * mpFileIdTag;

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

  Range localGid; // used only by ucd mesh, e.g. HOMME grid

    //! whether mesh is locally periodic in i or j
  int locallyPeriodic[2];

    //! whether mesh is globally periodic in i or j
  int globallyPeriodic[2];

    //! parallel data object, to be cached with ScdBox
  ScdParData parData;

    //! directory where data is stored for BIL reader
  std::string BIL_dir;

#ifdef USE_MPI
  ParallelComm *myPcomm;
#endif

    // read option
  bool noMesh;

    // read option
  bool noVars;

    // read option
  bool spectralMesh;

    // read option
  std::string partitionTagName;

    //! Helper class instance
  NCHelper* myHelper;
};

// inline functions
inline unsigned int ReadNC::number_dimensions() 
{
  return dimVals.size();
}

} // namespace moab

#endif
