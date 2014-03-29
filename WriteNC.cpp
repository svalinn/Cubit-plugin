#include "WriteNC.hpp"
#include "moab/CN.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/FileOptions.hpp"
#include "NCWriteHelper.hpp"

#include <fstream>
#include <map>
#include <set>

#include <iostream>

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

WriterIface *WriteNC::factory( Interface* iface )
  { return new WriteNC( iface ); }

WriteNC::WriteNC(Interface *impl)
    : mbImpl(impl), dbgOut(stderr), partMethod(ScdParData::ALLJORKORI), scdi(NULL),

#ifdef USE_MPI
  myPcomm(NULL),
#endif

  noMesh(false), noVars(false), /*spectralMesh(false), noMixedElements(false), noEdges(false),*/
  gatherSetRank(-1), mGlobalIdTag(0), isParallel(false),
   myHelper(NULL)
{
  assert(impl != NULL);
  impl->query_interface(mWriteIface);
}

WriteNC::~WriteNC()
{
  mbImpl->release_interface(mWriteIface);
  if (myHelper != NULL)
    delete myHelper;
}

//! writes out a file
ErrorCode WriteNC::write_file(const char *file_name,
                                  const bool overwrite,
                                  const FileOptions& options,
                                  const EntityHandle *file_set,
                                  const int num_set,
                                  const std::vector<std::string>&,
                                  const Tag*,
                                  int,
                                  int )
{

  ErrorCode rval;
  // See if opts has variable(s) specified
  std::vector<std::string> var_names;
  std::vector<int> tstep_nums;
  std::vector<double> tstep_vals;

  // Get and cache predefined tag handles
  int dum_val = 0;
  rval = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, mGlobalIdTag, MB_TAG_DENSE , &dum_val);
  if (MB_SUCCESS != rval)
    return rval;

  // num set has to be 1, we will write only one set, the original file set used to load
  if (num_set!=1)
    ERRORR(MB_FAILURE, "we should write only one set, the file set used to read data into");

  rval = parse_options(options, var_names, tstep_nums, tstep_vals);
  ERRORR(rval, "Trouble parsing option string.");

  // important to create some data that will be used to write the file; dimensions, variables, etc
  // new variables still need to have some way of defining their dimensions
  // maybe it will be passed as write options
  rval = process_conventional_tags(*file_set);
  ERRORR(rval, "Trouble getting conventional tags.");

  // Create the file ; assume we will overwrite always, for the time being
  dbgOut.tprintf(1, "creating file %s\n", file_name);
  fileName = file_name;
  int success;

#ifdef PNETCDF_FILE
  int cmode= overwrite ? NC_CLOBBER : NC_NOCLOBBER;
  if (isParallel)
    success = NCFUNC(create)(myPcomm->proc_config().proc_comm(), file_name, cmode, MPI_INFO_NULL, &fileId);
  else
    success = NCFUNC(create)(MPI_COMM_SELF, file_name, cmode, MPI_INFO_NULL, &fileId);
#else
    // this is regular netcdf file
  success = NCFUNC(create)(file_name,  overwrite ? NC_CLOBBER : NC_NOCLOBBER, &fileId);
#endif
  ERRORS(success, "failed to create file");

  /* int nc_def_dim (int ncid, const char *name, size_t len, int *dimidp);
   * example:  status = nc_def_dim(fileId, "lat", 18L, &latid);
   * */

  /*
   * int nc_def_var (int ncid, const char *name, nc_type xtype,
                     int ndims, const int dimids[], int *varidp);
     example: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c/nc_005fdef_005fvar.html#nc_005fdef_005fvar
   */

    /*
     * Write an Entire Variable: nc_put_var_ type (double, int)
     * int nc_put_var_double(int ncid, int varid, const double *dp);
     */
  /*
   * Write an Array of Values: nc_put_vara_ type
   * int nc_put_vara_double(int ncid, int varid, const size_t start[],
                            const size_t count[], const double *dp);
   */


  return MB_SUCCESS;
}

ErrorCode WriteNC::parse_options(const FileOptions& opts, std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
                                std::vector<double>& tstep_vals)
{
  int tmpval;
  if (MB_SUCCESS == opts.get_int_option("DEBUG_IO", 1, tmpval)) {
    dbgOut.set_verbosity(tmpval);
    dbgOut.set_prefix("NCWrite");
  }

  ErrorCode rval = opts.get_strs_option("VARIABLE", var_names);
  if (MB_TYPE_OUT_OF_RANGE == rval)
    noVars = true;
  else
    noVars = false;
  opts.get_ints_option("TIMESTEP", tstep_nums);
  opts.get_reals_option("TIMEVAL", tstep_vals);
  rval = opts.get_null_option("NOMESH");
  if (MB_SUCCESS == rval)
    noMesh = true;

 /* these are not used yet, maybe later
  rval = opts.get_null_option("SPECTRAL_MESH");
  if (MB_SUCCESS == rval)
    spectralMesh = true;

  rval = opts.get_null_option("NO_MIXED_ELEMENTS");
  if (MB_SUCCESS == rval)
    noMixedElements = true;

  rval = opts.get_null_option("NO_EDGES");
  if (MB_SUCCESS == rval)
    noEdges = true;*/

  if (2 <= dbgOut.get_verbosity()) {
    if (!var_names.empty()) {
      std::cerr << "Variables requested: ";
      for (unsigned int i = 0; i < var_names.size(); i++)
        std::cerr << var_names[i];
      std::cerr << std::endl;
    }
    if (!tstep_nums.empty()) {
      std::cerr << "Timesteps requested: ";
      for (unsigned int i = 0; i < tstep_nums.size(); i++)
        std::cerr << tstep_nums[i];
      std::cerr << std::endl;
    }
    if (!tstep_vals.empty()) {
      std::cerr << "Time vals requested: ";
      for (unsigned int i = 0; i < tstep_vals.size(); i++)
        std::cerr << tstep_vals[i];
      std::cerr << std::endl;
    }
  }

  // the gather set will be important in parallel, this will be the rank that will accumulate the data
  // to be written in serial
  // improvement will be for writing in true parallel
  rval = opts.get_int_option("GATHER_SET", 0, gatherSetRank);
  if (MB_TYPE_OUT_OF_RANGE == rval) {
    mWriteIface->report_error("Invalid value for GATHER_SET option");
    return rval;
  }
// FIXME: copied from readnc, may need revise
#ifdef USE_MPI
  isParallel = (opts.match_option("PARALLEL","READ_PART") != MB_ENTITY_NOT_FOUND);

  if (!isParallel)
  // Return success here, since rval still has _NOT_FOUND from not finding option
  // in this case, myPcomm will be NULL, so it can never be used; always check for isParallel
  // before any use for myPcomm
    return MB_SUCCESS;

  int pcomm_no = 0;
  rval = opts.get_int_option("PARALLEL_COMM", pcomm_no);
  if (rval == MB_TYPE_OUT_OF_RANGE) {
    mWriteIface->report_error("Invalid value for PARALLEL_COMM option");
    return rval;
  }
  myPcomm = ParallelComm::get_pcomm(mbImpl, pcomm_no);
  if (0 == myPcomm) {
    myPcomm = new ParallelComm(mbImpl, MPI_COMM_WORLD);
  }
  const int rank = myPcomm->proc_config().proc_rank();
  dbgOut.set_rank(rank);

  int dum;
  rval = opts.match_option("PARTITION_METHOD", ScdParData::PartitionMethodNames, dum);
  if (rval == MB_FAILURE) {
    mWriteIface->report_error("Unknown partition method specified.");
    partMethod = ScdParData::ALLJORKORI;
  }
  else if (rval == MB_ENTITY_NOT_FOUND)
    partMethod = ScdParData::ALLJORKORI;
  else
    partMethod = dum;
#endif

  return MB_SUCCESS;
}

// this is the inverse process to create conventional tags
// will look at <pargal_source>/src/core/fileinfo.cpp, init dim, vars, atts
ErrorCode WriteNC::process_conventional_tags(EntityHandle fileSet)
{
  ErrorCode rval;
  // start copy
  Tag tag = 0;
  std::string tag_name = "__DIM_NAMES";
  const void * data = NULL;
  int num_dims = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, tag, MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __DIM_NAMES.");
  rval = mbImpl->tag_get_by_ptr(tag, &fileSet, 1, &data, &num_dims);
  ERRORR(rval, "Failed getting values for conventional tag __DIM_NAMES.");

  dbgOut.tprintf(1, "dim names size: %d\n", num_dims);

  std::size_t start = 0;
  const char * p = static_cast<const char *>(data);
  Tag dimValsTag = 0;
  tag_name = "__DIM_LENS";
  const void * valdata = NULL;
  int num_vals = 0;
  rval =  mbImpl->tag_get_handle(tag_name.c_str(), 0,  MB_TYPE_INTEGER,
      dimValsTag,  MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __DIM_LENS.");

  rval = mbImpl->tag_get_by_ptr(dimValsTag, &fileSet, 1, &valdata, &num_vals);
  ERRORR(rval, "Failed getting values for conventional tag __DIM_LENS.");

  dbgOut.tprintf(1, "num vals in dim lens tag %d\n", num_vals);
  const int * intp = static_cast<const int*>(valdata);
  int idx = -1;

  for (std::size_t i = 0; i != static_cast<std::size_t>(num_dims); ++i)
  {
    if (p[i] == '\0')
    {
      std::string dim_name(&p[start], i - start);
      ++idx;
      int sz = intp[idx];
      dimNames.push_back(dim_name);
      dimLens.push_back(sz);
      dbgOut.tprintf(2, "dimension %s has length %d\n", dim_name.c_str(), sz);
      // fixme: need info from moab to set unlimited dimension
      // currently assume each file has the same number of time dimensions
      /*if ((dim_name == "time") || (dim_name == "Time"))
        insert(dim_name,
            *(new pcdim(dim_name, sz * m_file_names.size(), true)));
      else
        insert(dim_name, *(new pcdim(dim_name, sz)));*/
      start = i + 1;
    }
  }

  Tag tagMeshType = 0;
  tag_name = "__MESH_TYPE";
  data = NULL;
  int sz = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, tagMeshType, MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __MESH_TYPE.");
  rval = mbImpl->tag_get_by_ptr(tagMeshType, &fileSet, 1, &data, &sz);
  ERRORR(rval, "Failed getting values for conventional tag __MESH_TYPE.");

  p = static_cast<const char *>(data);
  grid_type = std::string(&p[0], sz);
  dbgOut.tprintf(2, "mesh type: %s \n", grid_type.c_str());
  // read <__VAR_NAMES_LOCATIONS> tag
  Tag varLocTag = 0;
  tag_name = "__VAR_NAMES_LOCATIONS";
  const void * loc_data = NULL;
  int loc_sz = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0,
      MB_TYPE_INTEGER, varLocTag, MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __VAR_NAMES_LOCATIONS.");
  rval = mbImpl->tag_get_by_ptr(varLocTag, &fileSet, 1, &loc_data, &loc_sz);
  ERRORR(rval, "Failed getting values for conventional tag __VAR_NAMES_LOCATIONS.");
  const int * loc_p = static_cast<const int *>(loc_data);
  std::vector<int> varLoc(loc_sz);
  std::copy(loc_p, loc_p + loc_sz, varLoc.begin());

 /*
  std::map<int, std::string> locmap;
  locmap[0] = "VERTEX";
  locmap[1] = "NSEDGE";
  locmap[2] = "EWEDGE";
  if (grid_type == "MPAS")
  {
    locmap[3] = "POLYGON";
  }
  else
  {
    locmap[3] = "QUAD";
  }
  locmap[4] = "SET";
  locmap[5] = "EDGE";
  locmap[6] = "REGION";
  */
  int nthVar = 0;

  Tag tagVarNames = 0;
  tag_name = "__VAR_NAMES";
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, tagVarNames,
      MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __VAR_NAMES.");
  data = NULL;
  int num_vars = 0;
  rval = mbImpl->tag_get_by_ptr(tagVarNames, &fileSet, 1, &data, &num_vars);
  ERRORR(rval, "Failed getting values for conventional tag __VAR_NAMES.");
  p = static_cast<const char *>(data);
  start = 0;
  for (std::size_t i = 0; i != static_cast<std::size_t>(num_vars); ++i)
  {
    if (p[i] == '\0')
    {
      std::string var_name(&p[start], i - start);

      sz = 0;
      Tag dims_tag = 0;
      std::string dim_names = "__" + var_name + "_DIMS";
      rval = mbImpl->tag_get_handle(dim_names.c_str(), 0, MB_TYPE_OPAQUE,
          dims_tag, MB_TAG_ANY);
      ERRORR(rval, "FAILED to get tag for a variable dimensions");
      //fixme: doesn't handle variables have 0 dimension
      if (rval != moab::MB_SUCCESS)
      {
        start = i + 1;
        ++nthVar;
        continue;
      }
      rval = mbImpl->tag_get_length(dims_tag, sz);
      ERRORR(rval, " size of dimensions for variable");
      dbgOut.tprintf(2, "var name: %s has %d dimensions \n", var_name.c_str(), sz);
      //std::vector<const pcdim*> dims(sz, NULL);
      const void* ptr = NULL;
      rval = mbImpl->tag_get_by_ptr(dims_tag, &fileSet, 1, &ptr);
      //
      const Tag * ptags = static_cast<const moab::Tag*>(ptr);
      for (std::size_t j = 0; j != static_cast<std::size_t>(sz); ++j)
      {
        std::string dim_name;
        rval = mbImpl->tag_get_name(ptags[j], dim_name);
        ERRORR(rval, "name of tag for dimension");
        dbgOut.tprintf(2, "var name: %s has %d dimensions \n", var_name.c_str(), sz);
        //dims[j] = &(get_dim(dim_name));
      }
      /*insert(var_name,
          *(new fvar(var_name, moab::MB_TYPE_DOUBLE, dims,
              locmap[varLoc[nthVar++]])));*/
      start = i + 1;
    }
  }

  // attributes
#if 0
  std::vector<std::string> nameVec(num_vars+1);
        nameVec[0] = "GLOBAL";
        std::map<std::string, const fvar*>::iterator it = m_vars.begin();
        for (std::size_t i=1; it != m_vars.end(); ++it, ++i)
          nameVec[i] = it->first;

        for (std::size_t vec_counter = 0; vec_counter != nameVec.size(); ++vec_counter)
          {
      // read __<var_name>_ATTRIBS tag
      moab::Tag tag = 0;
      std::string tag_name = "__"+nameVec[vec_counter]+"_ATTRIBS";
      const void * data = NULL;
      int sz = 0;
      moab::ErrorCode rval = m_mb.tag_get_handle(tag_name.c_str(), 0, moab::MB_TYPE_OPAQUE, tag, moab::MB_TAG_ANY);
      if (rval != moab::MB_SUCCESS)
        throw pargal_except("Error: " + m_mb.get_error_string(rval),
                __FILE__, __LINE__, __PRETTY_FUNCTION__);

      rval = m_mb.tag_get_by_ptr(tag, &m_file_set, 1, &data, &sz);
      const char * p = static_cast<const char *>(data);
      std::string att_val(&p[0], sz);
      if (vec_counter == 0) nameVec[0]="MOAB_GLOBAL";
      const std::string& var_name = nameVec[vec_counter];

      // read __<var_name>_ATTRIBS_LEN tag
      moab::Tag attLenTag = 0;
      tag_name = tag_name + "_LEN";
      const void * len_data = NULL;
      int len_sz = 0;
      rval = m_mb.tag_get_handle(tag_name.c_str(), 0, moab::MB_TYPE_INTEGER, attLenTag, moab::MB_TAG_ANY);
      rval = m_mb.tag_get_by_ptr(attLenTag, &m_file_set, 1, &len_data, &len_sz);
      const int * len_p = static_cast<const int *>(len_data);
      std::vector<int> attLen(len_sz);
      std::copy(len_p, len_p+len_sz, attLen.begin());

      // create attribute
      insert(var_name, *(new pcatt(var_name, att_val, attLen)));
#endif
  return MB_SUCCESS;
}
#if 0


    void
    fileinfo::init_atts()
    {
      std::vector<std::string> nameVec(m_vars.size()+1);
      nameVec[0] = "GLOBAL";
      std::map<std::string, const fvar*>::iterator it = m_vars.begin();
      for (std::size_t i=1; it != m_vars.end(); ++it, ++i)
        nameVec[i] = it->first;

      for (std::size_t vec_counter = 0; vec_counter != nameVec.size(); ++vec_counter)
        {
    // read __<var_name>_ATTRIBS tag
    moab::Tag tag = 0;
    std::string tag_name = "__"+nameVec[vec_counter]+"_ATTRIBS";
    const void * data = NULL;
    int sz = 0;
    moab::ErrorCode rval = m_mb.tag_get_handle(tag_name.c_str(), 0, moab::MB_TYPE_OPAQUE, tag, moab::MB_TAG_ANY);
    if (rval != moab::MB_SUCCESS)
      throw pargal_except("Error: " + m_mb.get_error_string(rval),
              __FILE__, __LINE__, __PRETTY_FUNCTION__);

    rval = m_mb.tag_get_by_ptr(tag, &m_file_set, 1, &data, &sz);
    const char * p = static_cast<const char *>(data);
    std::string att_val(&p[0], sz);
    if (vec_counter == 0) nameVec[0]="MOAB_GLOBAL";
    const std::string& var_name = nameVec[vec_counter];

    // read __<var_name>_ATTRIBS_LEN tag
    moab::Tag attLenTag = 0;
    tag_name = tag_name + "_LEN";
    const void * len_data = NULL;
    int len_sz = 0;
    rval = m_mb.tag_get_handle(tag_name.c_str(), 0, moab::MB_TYPE_INTEGER, attLenTag, moab::MB_TAG_ANY);
    rval = m_mb.tag_get_by_ptr(attLenTag, &m_file_set, 1, &len_data, &len_sz);
    const int * len_p = static_cast<const int *>(len_data);
    std::vector<int> attLen(len_sz);
    std::copy(len_p, len_p+len_sz, attLen.begin());

    // create attribute
    insert(var_name, *(new pcatt(var_name, att_val, attLen)));
        }
    }

    void
    fileinfo::init_grid_type()
    {
      moab::Tag tag = 0;
      std::string tag_name = "__MESH_TYPE";
      const void * data = NULL;
      int sz = 0;
      moab::ErrorCode rval;
      rval = m_mb.tag_get_handle(tag_name.c_str(), 0, moab::MB_TYPE_OPAQUE, tag, moab::MB_TAG_ANY);
      rval = m_mb.tag_get_by_ptr(tag, &m_file_set, 1, &data, &sz);
      if (rval != moab::MB_SUCCESS)
        throw pargal_except("Error: " + m_mb.get_error_string(rval),
          __FILE__, __LINE__, __PRETTY_FUNCTION__);
      const char * p = static_cast<const char *>(data);
      m_grid_type = std::string(&p[0], sz);
    }

  return MB_SUCCESS;
}
#endif

} // end moab namespace
