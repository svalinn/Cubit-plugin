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
//! keep it this way , introduce another macro, used so far only for ucd mesh
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
  friend class ScdNCHelper;
  friend class UcdNCHelper;
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

  //! ENTLOCNSEDGE for north/south edge
  //! ENTLOCWEEDGE for west/east edge
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
    VarData() : varId(-1), numAtts(-1), entLoc(ENTLOCSET), numLev(1), sz(0), has_t(false) {}
    int varId;
    int numAtts;
    nc_type varDataType;
    std::vector<int> varDims; // The dimension indices making up this multi-dimensional variable
    std::map<std::string, AttData> varAtts;
    std::string varName;
    std::vector<Tag> varTags; // One tag for each time step, varTags[t]
    std::vector<void*> varDatas;
    std::vector<std::vector<NCDF_SIZE> > readStarts; // Start value for this [t][dim]
    std::vector<std::vector<NCDF_SIZE> > readCounts; // Number of data values for this [t][dim]
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
  ErrorCode get_attributes(int var_id, int num_atts, std::map<std::string, AttData>& atts,
                           const char *prefix = "");

  //! get all dimensions in the file
  ErrorCode get_dimensions(int file_id, std::vector<std::string>& dim_names, std::vector<int>& dim_vals);

  //! get the variable names and other info defined for this file
  ErrorCode get_variables();

  ErrorCode read_coordinate(const char* var_name, int lmin, int lmax,
                            std::vector<double>& cvals);

  //! make sure that localGid is properly initialized for ucd mesh
  //ErrorCode check_ucd_localGid(EntityHandle file_set);

  ErrorCode parse_options(const FileOptions& opts,
                          std::vector<std::string>& var_names,
                          std::vector<int>& tstep_nums,
                          std::vector<double>& tstep_vals);

  ErrorCode get_tag_to_set(VarData& var_data, int tstep_num, Tag& tagh);

  ErrorCode get_tag_to_nonset(VarData& var_data, int tstep_num, Tag& tagh, int num_lev);

  //! create nc conventional tags
  ErrorCode create_conventional_tags(ScdInterface* scdi, EntityHandle file_set,
                                     const std::vector<int>& tstep_nums);

  //! create a character string attString of attMap.  with '\0'
  //! terminating each attribute name, ';' separating the data type
  //! and value, and ';' separating one name/data type/value from
  //! the next'.  attLen stores the end position for each name/data
  //! type/ value.
  ErrorCode create_attrib_string(const std::map<std::string, AttData>& attMap, 
				 std::string& attString,
				 std::vector<int>& attLen);

  //! create COORDS tag for quads coordinate
  ErrorCode create_quad_coordinate_tag(EntityHandle file_set);

  //! Init info for dimensions that don't have corresponding 
  //! coordinate variables - this info is used for creating tags
  void init_dims_with_no_cvars_info();

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
  //! these should be taken out when we fix the dummy var info things
  std::set<std::string> dummyVarNames;
  std::vector<int> dimVals;
  std::string iName, jName, kName, tName;
  std::string iCName, jCName;

  //! global attribs
  std::map<std::string, AttData> globalAtts;

  //! variable info
  std::map<std::string, VarData> varInfo;

  //! dimensions of grid in file
  int gDims[6], tMin, tMax;

  //! dimensions of my part of grid
  int lDims[6];

  //! center dimensions of grid in file
  int gCDims[6];

  //! center dimensions of my part of grid
  int lCDims[6];

  //! values for i/j/k/t
  std::vector<double> ilVals, jlVals, klVals, tVals;

  //! center values for i/j
  std::vector<double> ilCVals, jlCVals;

  //! dimension numbers for i, j, k, t
  int iDim, jDim, kDim, tDim;

  //! center dimension numbers for i, j
  int iCDim, jCDim;

  //! number of the dimension of unlimited dimension, if any
  int numUnLim;

  //! Cached tags for reading. Note that all these tags are defined when the
  //! core is initialized.
  Tag mGlobalIdTag;

  //! this is a pointer to the file id tag that is passed from ReadParallel
  //! it gets deleted at the end of resolve sharing, but it will have same data
  //! as the global id tag
  //! global id tag is preserved, and is needed later on.
  const Tag* mpFileIdTag;

  //! debug stuff
  DebugOutput dbgOut;

  //! are we reading in parallel?
  bool isParallel;

  //! partitioning method
  int partMethod;

  //! whether mesh is locally periodic in i or j
  int locallyPeriodic[2];

  //! whether mesh is globally periodic in i or j
  int globallyPeriodic[2];

  //! parallel data object, to be cached with ScdBox
  ScdParData parData;

#ifdef USE_MPI
  ParallelComm* myPcomm;
#endif

  //! Read options
  bool noMesh;
  bool noVars;
  bool spectralMesh;
  std::string partitionTagName;
  int gatherSetRank;

  //! Helper class instance
  NCHelper* myHelper;
};

} // namespace moab

#endif
