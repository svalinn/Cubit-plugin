#include "ReadNC.hpp"
#include "NCHelper.hpp"

#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <dirent.h>

#include "moab/Core.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "FileOptions.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/SpectralMeshTool.hpp"

//#include "bil.h"

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
    if (err) {readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

ReaderIface* ReadNC::factory(Interface* iface) {
  return new ReadNC(iface);
}

ReadNC::ReadNC(Interface* impl) :
  mbImpl(impl), CPU_WORD_SIZE(-1), IO_WORD_SIZE(-1), fileId(-1), tMin(-1), tMax(-1), iDim(-1), jDim(-1), tDim(-1), iCDim(-1),
  jCDim(-1), numUnLim(-1), mCurrentMeshHandle(0), startVertex(0), startElem(0), mGlobalIdTag(0), mpFileIdTag(NULL), max_line_length(-1),
  max_str_length(-1), vertexOffset(0), dbgOut(stderr), isParallel(false), partMethod(-1),
#ifdef USE_MPI
  myPcomm(NULL), 
#endif
  noMesh(false), noVars(false), spectralMesh(false), myHelper(NULL)
{
  assert(impl != NULL);

  for (unsigned int i = 0; i < 6; i++) {
    gDims[i] = -1;
    lDims[i] = -1;
    gCDims[i] = -1;
    lCDims[i] = -1;
  }

  locallyPeriodic[0] = locallyPeriodic[1] = 0;
  globallyPeriodic[0] = globallyPeriodic[1] = 0;

  impl->query_interface(readMeshIface);
}

void ReadNC::reset() {
  CPU_WORD_SIZE = -1;
  IO_WORD_SIZE = -1;
  fileId = -1;
  tMin = tMax = -1;
  for (unsigned int i = 0; i < 6; i++) {
    gDims[i] = -1;
    lDims[i] = -1;
    gCDims[i] = -1;
    lCDims[i] = -1;
  }

  iDim = jDim = tDim = -1;
  iCDim = jCDim = -1;
  numUnLim = -1;
  mCurrentMeshHandle = 0;
  startVertex = startElem = 0;
  mGlobalIdTag = 0;
  max_line_length = -1;
  max_str_length = -1;
  vertexOffset = 0;
  dbgOut = stderr;
  mCurrentMeshHandle = 0;
  vertexOffset = 0;

#ifdef USE_MPI
  myPcomm = NULL;
#endif
}

ReadNC::~ReadNC() {
  mbImpl->release_interface(readMeshIface);
  if (myHelper != NULL)
    delete myHelper;
}

ErrorCode ReadNC::load_file(const char* file_name, const EntityHandle* file_set, const FileOptions& opts,
                            const ReaderIface::SubsetList* /*subset_list*/, const Tag* file_id_tag) {

  ErrorCode rval = MB_SUCCESS;

  //See if opts has variable(s) specified
  std::vector<std::string> var_names;
  std::vector<int> tstep_nums;
  std::vector<double> tstep_vals;
  /*
  if (file_id_tag)
    mGlobalIdTag = *file_id_tag;
  else {
    */
    //! get and cache predefined tag handles
  int dum_val = 0;
  rval = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, mGlobalIdTag, MB_TAG_DENSE | MB_TAG_CREAT, &dum_val);
  if (MB_SUCCESS != rval)
    return rval;
  mpFileIdTag = file_id_tag; // store the pointer to the tag ; if not null, set when global id tag
  // is set too, with the same data , duplicated

  std::string partition_tag_name;
  rval = parse_options(opts, var_names, tstep_nums, tstep_vals);
  ERRORR(rval, "Trouble parsing option string.");

  // Open the file
  dbgOut.tprintf(1, "Opening file %s\n", file_name);
  fileName = std::string(file_name);
  int success;

#ifdef PNETCDF_FILE
  if (isParallel)
    success = NCFUNC(open)(myPcomm->proc_config().proc_comm(), file_name, 0, MPI_INFO_NULL, &fileId);
  else
    success = NCFUNC(open)(MPI_COMM_SELF, file_name, 0, MPI_INFO_NULL, &fileId);
#else
  success = NCFUNC(open)(file_name, 0, &fileId);
#endif

  ERRORS(success, "Trouble opening file.");

  // BIL data

  if (BIL_mode_enabled(file_name)) {

    rval = get_BIL_dir();
    ERRORS(rval, "Failed to find directory with BIL data.");

    dbgOut.tprintf(1, "Reading BIL data from directory: %s\n", BIL_dir.c_str());

    rval = load_BIL(BIL_dir, file_set, opts, file_id_tag);
    ERRORR(rval, "Trouble reading BIL data.");

    return rval;
  }

  // end of BIL

  // Read the header (num dimensions, dimensions, num variables, global attribs)
  rval = read_header();
  ERRORR(rval, " ");

  // make sure there's a file set to put things in
  EntityHandle tmp_set;
  if (noMesh && !file_set) {
    ERRORR(MB_FAILURE, "NOMESH option requires non-NULL file set on input.\n");
  }
  else if (!file_set || (file_set && *file_set == 0)) {
    rval = mbImpl->create_meshset(MESHSET_SET, tmp_set);
    ERRORR(rval, "Trouble creating file set.");
  }
  else
    tmp_set = *file_set;

  // get the scd interface
  ScdInterface *scdi = NULL;
  rval = mbImpl->query_interface(scdi);
  if (!scdi)
    return MB_FAILURE;

  if (myHelper != NULL)
    delete myHelper;

  myHelper = NCHelper::get_nc_helper(this, fileId, opts);
  if (myHelper == NULL) {
    ERRORR(MB_FAILURE, "Failed to get NCHelper class instance.");
  }

  rval = myHelper->init_mesh_vals(opts, tmp_set);
  ERRORR(rval, "Trouble initializing mesh values.");

  // Create mesh vertex/edge/face sequences
  Range faces;
  if (noMesh && !noVars) {
    rval = check_verts_faces(tmp_set);
    ERRORR(rval, "Mesh characteristics didn't match from last read.\n");
  }
  else if (!noMesh) {
    rval = myHelper->create_mesh(scdi, opts, tmp_set, faces);
    ERRORR(rval, "Trouble creating mesh.");
  }

  // Read variables onto grid
  if (!noVars) {
    rval = myHelper->read_variables(tmp_set, var_names, tstep_nums);
    if (MB_FAILURE == rval)
      return rval;
  }
  else {
    // read dimension variable by default, the ones that are also variables
    std::vector<std::string> filteredDimNames;
    for (unsigned int i = 0; i < dimNames.size(); i++) {
      std::map<std::string, VarData>::iterator mit = varInfo.find(dimNames[i]);
      if (mit != varInfo.end())
        filteredDimNames.push_back(dimNames[i]);
    }
    rval = myHelper->read_variables(tmp_set, filteredDimNames, tstep_nums);
    if (MB_FAILURE == rval)
      return rval;
  }

#ifdef USE_MPI
  // create partition set, and populate with elements
  if (isParallel) {
    EntityHandle partn_set;
    rval = mbImpl->create_meshset(MESHSET_SET, partn_set);
    ERRORR(rval, "Trouble creating partition set.");

    rval = mbImpl->add_entities(partn_set, faces);
    ERRORR(rval, "Couldn't add new faces to partition set.");

    Range verts;
    rval = mbImpl->get_connectivity(faces, verts);
    ERRORR(rval, "Couldn't get verts of faces");

    rval = mbImpl->add_entities(partn_set, verts);
    ERRORR(rval, "Couldn't add new verts to partition set.");

    myPcomm->partition_sets().insert(partn_set);

    //write partition tag name on partition set
    Tag part_tag;
    rval = mbImpl->tag_get_handle(partitionTagName.c_str(), 1, MB_TYPE_INTEGER, part_tag);
    if (MB_SUCCESS != rval) {
      // fall back to the partition tag
      part_tag = myPcomm->partition_tag();
    }

    int dum_rank = myPcomm->proc_config().proc_rank();
    rval = mbImpl->tag_set_data(part_tag, &partn_set, 1, &dum_rank);
    if (MB_SUCCESS != rval)
      return rval;
  }
#endif

  // create nc conventional tags when loading header info only
  if (noMesh && noVars) {
    rval = create_tags(scdi, tmp_set, tstep_nums);
    ERRORR(rval, "Trouble creating nc conventional tags.");
  }

  mbImpl->release_interface(scdi);

  // close the file
  success = NCFUNC(close)(fileId);
  ERRORS(success, "Trouble closing file.");

  return MB_SUCCESS;
}
    
ErrorCode ReadNC::load_BIL(std::string, const EntityHandle*, const FileOptions&, const Tag*) {
  /*
   BIL_Init( MPI_COMM_WORLD );

   void ** buffer;

   DIR * dir;
   struct dirent * ent;
   dir = opendir(dir_name.c_str());
   if (dir != NULL) {
   while ((ent = readdir(dir)) != NULL) {
   if (strlen(ent->d_name) > 3) { //filter out . and ..

   dbgOut.tprintf(1,"reading block from %s\n",ent->d_name);

   int num_dims = 3;
   int time_d = 1;
   int lev_d  = 26;
   int ncol_d = 3458;
   int block_start[3] = {0,0,0};
   int block_size[3]  = {time_d, lev_d, ncol_d};
   const char * file_name = ent->d_name;
   const char * var_name = "T";

   BIL_Add_block_nc(num_dims, block_start, block_size,
   file_name, var_name, buffer);
   }
   }
   closedir (dir);
   }

   BIL_Read();

   BIL_Finalize();
   */
  return MB_SUCCESS;
}

ErrorCode ReadNC::get_BIL_dir() {
  std::map<std::string, AttData> dirAtt;
  ErrorCode result = get_attributes(NC_GLOBAL, 1, dirAtt);
  ERRORR(result, "Failed to get BIL_DIR attribute");

  std::string attname;
  std::map<std::string, AttData>::iterator attIt = dirAtt.find("BIL_DIR");

  unsigned int sz = attIt->second.attLen;
  char *att_data = (char *) malloc(sz + 1);
  att_data[sz] = '\000';
  int success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (char*) att_data);
  ERRORS(success, "Trouble getting BIL data directory.");
  BIL_dir = std::string(att_data);

  return MB_SUCCESS;
}

bool ReadNC::BIL_mode_enabled(const char* file_name) {
  std::string file_path = std::string(file_name);
  int idx = file_path.find_last_of("/");
  std::string file = file_path.substr(idx + 1);

  return (file == "BIL_DIR.nc");
}

ErrorCode ReadNC::parse_options(const FileOptions& opts, std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
                                std::vector<double>& tstep_vals) {
  int tmpval;
  if (MB_SUCCESS == opts.get_int_option("DEBUG_IO", 1, tmpval)) {
    dbgOut.set_verbosity(tmpval);
    dbgOut.set_prefix("NC ");
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

  rval = opts.get_null_option("SPECTRAL_MESH");
  if (MB_SUCCESS == rval)
    spectralMesh = true;

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

#ifdef USE_MPI
  //TODO handle options better
  //rval = opts.get_option("PARTITION", partitionTagName);
  /*
   part = (rval != MB_ENTITY_NOT_FOUND);
   rval = opts.match_option("PARALLEL", "READ_PART");
   parread = (rval != MB_ENTITY_NOT_FOUND);
   rval = opts.get_null_option("TRIVIAL_PARTITION");
   partriv = (rval != MB_ENTITY_NOT_FOUND);
   */
  isParallel = (opts.match_option("PARALLEL","READ_PART") != MB_ENTITY_NOT_FOUND);

  if (!isParallel)
  // return success here, since rval still has _NOT_FOUND from not finding option
  // in this case, myPcomm will be NULL, so it can never be used; always check for isParallel 
  // before any use for myPcomm
    return MB_SUCCESS;
  
  int pcomm_no = 0;
  rval = opts.get_int_option("PARALLEL_COMM", pcomm_no);
  if (rval == MB_TYPE_OUT_OF_RANGE) {
    readMeshIface->report_error("Invalid value for PARALLEL_COMM option");
    return rval;
  }
  myPcomm = ParallelComm::get_pcomm(mbImpl, pcomm_no);
  if (0 == myPcomm) {
    myPcomm = new ParallelComm(mbImpl, MPI_COMM_WORLD);
  }
  const int rank = myPcomm->proc_config().proc_rank();
  dbgOut.set_rank(rank);

  const char *part_options[] = {"alljorkori", "alljkbal", "sqij", "sqjk",
      "TRIVIAL_PARTITION" };
  int dum;
  rval = opts.match_option("PARTITION_METHOD", part_options, dum);
  if (rval == MB_FAILURE) {
    readMeshIface->report_error("Unknown partition method specified.");
    partMethod = ScdParData::ALLJORKORI;
  }
  else if (rval == MB_ENTITY_NOT_FOUND)
    partMethod = ScdParData::ALLJORKORI;
  else
    partMethod = dum;
#endif

  return MB_SUCCESS;
}

// In a script, the ReadNC class instance can get out of scope (and deleted). In that
// case, the localGid (initialized properly when the mesh was created) will be lost,
// so it has to be properly refilled with the Global Ids of the local vertices
ErrorCode ReadNC::check_ucd_localGid(EntityHandle tmp_set) {
  if (noMesh && localGid.empty()) {
    // we need to populate localGid range with the gids of vertices from the tmp_set
    // localGid is important in reading the variable data into the nodes
    // also, for our purposes, localGid is truly the GLOBAL_ID tag data, not other
    // file_id tags that could get passed around in other scenarios for parallel reading
    // for nodal_partition, this local gid is easier, should be initialized with only
    // the owned nodes

    // we need to get all vertices from tmp_set (it is the input set in no_mesh scenario)
    Range local_verts;
    ErrorCode rval = mbImpl->get_entities_by_dimension(tmp_set, 0, local_verts);
    if (MB_FAILURE == rval)
      return rval;

    if (!local_verts.empty()) {
      std::vector<int> gids(local_verts.size());

      // !IMPORTANT : this has to be the GLOBAL_ID tag
      rval = mbImpl->tag_get_data(mGlobalIdTag, local_verts, &gids[0]);
      if (MB_FAILURE == rval)
        return rval;

      // this will do a smart copy
      std::copy(gids.begin(), gids.end(), range_inserter(localGid));
    }
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::check_verts_faces(EntityHandle file_set) {
  // check parameters on this read against what was on the mesh from last read
  // get the number of vertices
  int num_verts;
  ErrorCode rval = mbImpl->get_number_entities_by_dimension(file_set, 0, num_verts);
  ERRORR(rval, "Trouble getting number of vertices.");

  // check against parameters
  //int expected_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2] + 1);
  //if (num_verts != expected_verts)
  //ERRORR(MB_FAILURE, "Number of vertices doesn't match.");

  // check the number of elements too
  int num_elems;
  rval = mbImpl->get_number_entities_by_dimension(file_set, (-1 == lDims[2] ? 2 : 3), num_elems);
  ERRORR(rval, "Trouble getting number of elements.");

  // check against parameters
  //int expected_elems = (lDims[3] - lDims[0]) * (lDims[4] - lDims[1]) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2]);
  //if (num_elems != expected_elems)
  //ERRORR(MB_FAILURE, "Number of elements doesn't match.");

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_tag_to_set(VarData& var_data, int tstep_num, Tag& tagh) {
  std::ostringstream tag_name;
  if ((!var_data.has_t) || (var_data.varDims.size() <= 1))
    tag_name << var_data.varName;
  else if (!tstep_num) {
    std::string tmp_name = var_data.varName + "0";
    tag_name << tmp_name.c_str();
  }
  else
    tag_name << var_data.varName << tstep_num;
  ErrorCode rval = MB_SUCCESS;
  tagh = 0;
  switch (var_data.varDataType) {
    case NC_BYTE:
    case NC_CHAR:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 0, MB_TYPE_OPAQUE, tagh, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
      break;
    case NC_DOUBLE:
    case NC_FLOAT:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 0, MB_TYPE_DOUBLE, tagh, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
      break;
    case NC_INT:
    case NC_SHORT:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 0, MB_TYPE_INTEGER, tagh, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
      break;
    default:
      std::cerr << "Unrecognized data type for tag " << tag_name << std::endl;
      rval = MB_FAILURE;
  }

  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.str().c_str());

  return rval;
}

ErrorCode ReadNC::get_tag(VarData& var_data, int tstep_num, Tag& tagh, int num_lev) {
  std::ostringstream tag_name;
  if (!tstep_num) {
    std::string tmp_name = var_data.varName + "0";
    tag_name << tmp_name.c_str();
  }
  else
    tag_name << var_data.varName << tstep_num;
  ErrorCode rval = MB_SUCCESS;
  tagh = 0;
  switch (var_data.varDataType) {
    case NC_BYTE:
    case NC_CHAR:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), num_lev, MB_TYPE_OPAQUE, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
      break;
    case NC_DOUBLE:
    case NC_FLOAT:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), num_lev, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
      break;
    case NC_INT:
    case NC_SHORT:
      rval = mbImpl->tag_get_handle(tag_name.str().c_str(), num_lev, MB_TYPE_INTEGER, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
      break;
    default:
      std::cerr << "Unrecognized data type for tag " << tag_name.str() << std::endl;
      rval = MB_FAILURE;
  }

  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.str().c_str());

  return rval;
}

void ReadNC::init_dims_with_no_cvars_info() {
  // hack: look at all dimensions, and see if we have one that does not appear in the list of varInfo names
  // right now, candidates are ncol and nbnd
  // for them, create dummy tags
  for (unsigned int i = 0; i < dimNames.size(); i++)
  {
    // if there is a var with this name, skip, we are fine; if not, create a varInfo...
    if (varInfo.find(dimNames[i]) != varInfo.end())
      continue; // we already have a variable with this dimension name

    int sizeTotalVar = varInfo.size();
    std::string var_name(dimNames[i]);
    VarData &data = varInfo[var_name];
    data.varName = std::string(var_name);
    data.varId =sizeTotalVar;
    data.varTags.resize(1, 0);
    data.varDataType = NC_DOUBLE; // could be int, actually, but we do not really need the type
    data.varDims.resize(1);
    data.varDims[0]= (int)i;
    data.numAtts=0;
    data.entLoc = ENTLOCSET;
    dbgOut.tprintf(2, "Dummy varInfo created for dimension %s\n", dimNames[i].c_str());
    dummyVarNames.insert(dimNames[i]);
  }
}

ErrorCode ReadNC::read_coordinate(const char* var_name, int lmin, int lmax, std::vector<double>& cvals) {
  std::map<std::string, VarData>::iterator vmit = varInfo.find(var_name);
  if (varInfo.end() == vmit)
    return MB_FAILURE;

  // check to make sure it's a float or double
  int fail;
  NCDF_SIZE tmin = lmin, tcount = lmax - lmin + 1;
  NCDF_DIFF dum_stride = 1;
  if (NC_DOUBLE == (*vmit).second.varDataType) {
    cvals.resize(tcount);
    fail = NCFUNCA(get_vars_double)(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &cvals[0]);
    if (fail)
      ERRORS(MB_FAILURE, "Failed to get coordinate values.");
  }
  else if (NC_FLOAT == (*vmit).second.varDataType) {
    std::vector<float> tcvals(tcount);
    fail = NCFUNCA(get_vars_float)(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &tcvals[0]);
    if (fail)
      ERRORS(MB_FAILURE, "Failed to get coordinate values.");
    std::copy(tcvals.begin(), tcvals.end(), cvals.begin());
  }
  else ERRORR(MB_FAILURE, "Wrong data type for coordinate variable.");

  return MB_SUCCESS;
}

ErrorCode ReadNC::read_header() {
  CPU_WORD_SIZE = sizeof(double);
  IO_WORD_SIZE = sizeof(double);

  dbgOut.tprint(1, "Reading header...\n");

  // get the global attributes
  int numgatts;
  int success;
  success = NCFUNC(inq_natts )(fileId, &numgatts);
  ERRORS(success, "Couldn't get number of global attributes.");

  // read attributes into globalAtts
  ErrorCode result = get_attributes(NC_GLOBAL, numgatts, globalAtts);
  ERRORR(result, "Getting attributes.");
  dbgOut.tprintf(1, "Read %u attributes\n", (unsigned int) globalAtts.size());

  // read in dimensions into dimVals
  result = get_dimensions(fileId, dimNames, dimVals);
  ERRORR(result, "Getting dimensions.");
  dbgOut.tprintf(1, "Read %u dimensions\n", (unsigned int) dimVals.size());

  // read in variables into varInfo
  result = get_variables();
  ERRORR(result, "Getting variables.");
  dbgOut.tprintf(1, "Read %u variables\n", (unsigned int) varInfo.size());

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_attributes(int var_id, int num_atts, std::map<std::string, AttData>& atts, const char* prefix) {
  char dum_name[120];

  for (int i = 0; i < num_atts; i++) {
    // get the name
    int success = NCFUNC(inq_attname)(fileId, var_id, i, dum_name);
    ERRORS(success, "Trouble getting attribute name.");

    AttData &data = atts[std::string(dum_name)];
    data.attName = std::string(dum_name);
    success = NCFUNC(inq_att)(fileId, var_id, dum_name, &data.attDataType, &data.attLen);
    ERRORS(success, "Trouble getting attribute info.");
    data.attVarId = var_id;

    dbgOut.tprintf(2, "%sAttribute %s: length=%u, varId=%d, type=%d\n", (prefix ? prefix : ""), data.attName.c_str(),
        (unsigned int) data.attLen, data.attVarId, data.attDataType);
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_dimensions(int file_id, std::vector<std::string>& dim_names, std::vector<int>& dim_vals) {
  // get the number of dimensions
  int num_dims;
  int success = NCFUNC(inq_ndims)(file_id, &num_dims);
  ERRORS(success, "Trouble getting number of dimensions.");

  if (num_dims > NC_MAX_DIMS) {
    readMeshIface->report_error("ReadNC: File contains %d dims but NetCDF library supports only %d\n", num_dims, (int) NC_MAX_DIMS);
    return MB_FAILURE;
  }

  char dim_name[NC_MAX_NAME + 1];
  NCDF_SIZE dum_len;
  dim_names.resize(num_dims);
  dim_vals.resize(num_dims);

  for (int i = 0; i < num_dims; i++) {
    success = NCFUNC(inq_dim)(file_id, i, dim_name, &dum_len);
    ERRORS(success, "Trouble getting dimension info.");

    dim_vals[i] = dum_len;
    dim_names[i] = std::string(dim_name);

    dbgOut.tprintf(2, "Dimension %s, length=%u\n", dim_name, (unsigned int) dum_len);
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_variables() {
  // first cache the number of time steps
  std::vector<std::string>::iterator vit = std::find(dimNames.begin(), dimNames.end(), "time");
  if (vit == dimNames.end())
    vit = std::find(dimNames.begin(), dimNames.end(), "t");

  int ntimes = 0;
  if (vit != dimNames.end())
    ntimes = dimVals[vit - dimNames.begin()];
  if (!ntimes)
    ntimes = 1;

  // get the number of variables
  int num_vars;
  int success = NCFUNC(inq_nvars)(fileId, &num_vars);
  ERRORS(success, "Trouble getting number of variables.");

  if (num_vars > NC_MAX_VARS) {
    readMeshIface->report_error("ReadNC: File contains %d vars but NetCDF library supports only %d\n", num_vars, (int) NC_MAX_VARS);
    return MB_FAILURE;
  }

  char var_name[NC_MAX_NAME + 1];
  int var_ndims;

  for (int i = 0; i < num_vars; i++) {
    // get the name first, so we can allocate a map iterate for this var
    success = NCFUNC(inq_varname )(fileId, i, var_name);
    ERRORS(success, "Trouble getting var name.");
    VarData &data = varInfo[std::string(var_name)];
    data.varName = std::string(var_name);
    data.varId = i;
    data.varTags.resize(ntimes, 0);

    // get the data type
    success = NCFUNC(inq_vartype)(fileId, i, &data.varDataType);
    ERRORS(success, "Trouble getting variable data type.");

    // get the number of dimensions, then the dimensions
    success = NCFUNC(inq_varndims)(fileId, i, &var_ndims);
    ERRORS(success, "Trouble getting number of dims of a variable.");
    data.varDims.resize(var_ndims);

    success = NCFUNC(inq_vardimid)(fileId, i, &data.varDims[0]);
    ERRORS(success, "Trouble getting variable dimensions.");

    // finally, get the number of attributes, then the attributes
    success = NCFUNC(inq_varnatts)(fileId, i, &data.numAtts);
    ERRORS(success, "Trouble getting number of dims of a variable.");

    // print debug info here so attribute info comes afterwards
    dbgOut.tprintf(2, "Variable %s: Id=%d, numAtts=%d, datatype=%d, num_dims=%u\n", data.varName.c_str(), data.varId, data.numAtts,
        data.varDataType, (unsigned int) data.varDims.size());

    ErrorCode rval = get_attributes(i, data.numAtts, data.varAtts, "   ");
    ERRORR(rval, "Trouble getting attributes for a variable.");
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::read_tag_values(const char*, const char*, const FileOptions&, std::vector<int>&, const SubsetList*) {
  return MB_FAILURE;
}

ErrorCode ReadNC::create_tags(ScdInterface* scdi, EntityHandle file_set, const std::vector<int>& tstep_nums) {
  ErrorCode rval;
  std::string tag_name;

  // <__NUM_DIMS>
  Tag numDimsTag = 0;
  tag_name = "__NUM_DIMS";
  int numDims = dimNames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 1, MB_TYPE_INTEGER, numDimsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_DIMS tag.");
  rval = mbImpl->tag_set_data(numDimsTag, &file_set, 1, &numDims);
  ERRORR(rval, "Trouble setting data for __NUM_DIMS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__NUM_VARS>
  Tag numVarsTag = 0;
  tag_name = "__NUM_VARS";
  int numVars = varInfo.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 1, MB_TYPE_INTEGER, numVarsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_VARS tag.");
  rval = mbImpl->tag_set_data(numVarsTag, &file_set, 1, &numVars);
  ERRORR(rval, "Trouble setting data for __NUM_VARS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__DIM_NAMES>
  Tag dimNamesTag = 0;
  tag_name = "__DIM_NAMES";
  std::string dimnames;
  unsigned int dimNamesSz = dimNames.size();
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    dimnames.append(dimNames[i]);
    dimnames.push_back('\0');
  }
  int dimnamesSz = dimnames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, dimNamesTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __DIM_NAMES tag.");
  const void* ptr = dimnames.c_str();
  rval = mbImpl->tag_set_by_ptr(dimNamesTag, &file_set, 1, &ptr, &dimnamesSz);
  ERRORR(rval, "Trouble setting data for __DIM_NAMES tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__DIM_VALUES>
  Tag dimValsTag = 0;
  tag_name = "__DIM_VALUES";
  int dimValsSz = (int)dimVals.size();

  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_INTEGER, dimValsTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __DIM_VALUES tag.");
  ptr = &(dimVals[0]);
  rval = mbImpl->tag_set_by_ptr(dimValsTag, &file_set, 1, &ptr, &dimValsSz);
  ERRORR(rval, "Trouble setting data for __DIM_VALUES tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__VAR_NAMES>
  Tag varNamesTag = 0;
  tag_name = "__VAR_NAMES";
  std::string varnames;
  std::map<std::string, VarData>::iterator mapIter;
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    varnames.append(mapIter->first);
    varnames.push_back('\0');
  }
  int varnamesSz = varnames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, varNamesTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __VAR_NAMES tag.");
  ptr = varnames.c_str();
  rval = mbImpl->tag_set_by_ptr(varNamesTag, &file_set, 1, &ptr, &varnamesSz);
  ERRORR(rval, "Trouble setting data for __VAR_NAMES tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // __<dim_name>_LOC_MINMAX
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    if (dimNames[i] == "time") {
      std::stringstream ss_tag_name;
      ss_tag_name << "__" << dimNames[i] << "_LOC_MINMAX";
      tag_name = ss_tag_name.str();
      Tag tagh = 0;
      std::vector<int> val(2, 0);
      val[0] = tMin;
      val[1] = tMax;
      rval = mbImpl->tag_get_handle(tag_name.c_str(), 2, MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT);
      ERRORR(rval, "Trouble creating __<dim_name>_LOC_MINMAX tag.");
      rval = mbImpl->tag_set_data(tagh, &file_set, 1, &val[0]);
      ERRORR(rval, "Trouble setting data for __<dim_name>_LOC_MINMAX tag.");
      if (MB_SUCCESS == rval)
        dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
    }
  }

  // __<dim_name>_LOC_VALS
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    if (dimNames[i] != "time")
      continue;
    std::vector<int> val;
    if (!tstep_nums.empty())
      val = tstep_nums;
    else {
      val.resize(tVals.size());
      for (unsigned int j = 0; j != tVals.size(); ++j)
        val[j] = j;
    }
    Tag tagh = 0;
    std::stringstream ss_tag_name;
    ss_tag_name << "__" << dimNames[i] << "_LOC_VALS";
    tag_name = ss_tag_name.str();
    rval = mbImpl->tag_get_handle(tag_name.c_str(), val.size(), MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<dim_name>_LOC_VALS tag.");
    rval = mbImpl->tag_set_data(tagh, &file_set, 1, &val[0]);
    ERRORR(rval, "Trouble setting data for __<dim_name>_LOC_VALS tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // __<var_name>_DIMS
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    Tag varNamesDimsTag = 0;
    std::stringstream ss_tag_name;
    ss_tag_name << "__" << mapIter->first << "_DIMS";
    tag_name = ss_tag_name.str();
    unsigned int varDimSz = varInfo[mapIter->first].varDims.size();
    if (varDimSz == 0)
      continue;
    varInfo[mapIter->first].varTags.resize(varDimSz, 0);
    for (unsigned int i = 0; i != varDimSz; ++i) {
      Tag tmptag = 0;
      std::string tmptagname = dimNames[varInfo[mapIter->first].varDims[i]];
      mbImpl->tag_get_handle(tmptagname.c_str(), 0, MB_TYPE_OPAQUE, tmptag, MB_TAG_ANY);
      varInfo[mapIter->first].varTags[i] = tmptag;
    }
    rval = mbImpl->tag_get_handle(tag_name.c_str(), varDimSz, MB_TYPE_HANDLE, varNamesDimsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<var_name>_DIMS tag.");
    rval = mbImpl->tag_set_data(varNamesDimsTag, &file_set, 1, &(varInfo[mapIter->first].varTags[0]));
    ERRORR(rval, "Trouble setting data for __<var_name>_DIMS tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // <PARTITION_METHOD>
  Tag part_tag = scdi->part_method_tag();
  if (!part_tag)
    ERRORR(MB_FAILURE, "Trouble getting partition method tag.");
  rval = mbImpl->tag_set_data(part_tag, &file_set, 1, &partMethod);
  ERRORR(rval, "Trouble setting data for PARTITION_METHOD tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__GLOBAL_ATTRIBS>
  tag_name = "__GLOBAL_ATTRIBS";
  Tag globalAttTag = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, globalAttTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __GLOBAL_ATTRIBS tag.");
  std::string gattVal;
  std::vector<int> gattLen;
  rval = create_attrib_string(globalAtts, gattVal, gattLen);
  ERRORR(rval, "Trouble creating attribute strings.");
  const void* gattptr = gattVal.c_str();
  int globalAttSz = gattVal.size();
  rval = mbImpl->tag_set_by_ptr(globalAttTag, &file_set, 1, &gattptr, &globalAttSz);
  ERRORR(rval, "Trouble setting data for __GLOBAL_ATTRIBS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__GLOBAL_ATTRIBS_LEN>
  tag_name = "__GLOBAL_ATTRIBS_LEN";
  Tag globalAttLenTag = 0;
  if (gattLen.size() == 0)
    gattLen.push_back(0);
  rval = mbImpl->tag_get_handle(tag_name.c_str(), gattLen.size(), MB_TYPE_INTEGER, globalAttLenTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __GLOBAL_ATTRIBS_LEN tag.");
  rval = mbImpl->tag_set_data(globalAttLenTag, &file_set, 1, &gattLen[0]);
  ERRORR(rval, "Trouble setting data for __GLOBAL_ATTRIBS_LEN tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // __<var_name>_ATTRIBS and __<var_name>_ATTRIBS_LEN
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    std::stringstream ssTagName;
    ssTagName << "__" << mapIter->first << "_ATTRIBS";
    tag_name = ssTagName.str();
    Tag varAttTag = 0;
    rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, varAttTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
    ERRORR(rval, "Trouble creating __<var_name>_ATTRIBS tag.");
    std::string varAttVal;
    std::vector<int> varAttLen;
    rval = create_attrib_string(mapIter->second.varAtts, varAttVal, varAttLen);
    ERRORR(rval, "Trouble creating attribute strings.");
    const void* varAttPtr = varAttVal.c_str();
    int varAttSz = varAttVal.size();
    rval = mbImpl->tag_set_by_ptr(varAttTag, &file_set, 1, &varAttPtr, &varAttSz);
    ERRORR(rval, "Trouble setting data for __<var_name>_ATTRIBS tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
    if (varAttLen.size() == 0)
      varAttLen.push_back(0);
    ssTagName << "_LEN";
    tag_name = ssTagName.str();
    Tag varAttLenTag = 0;
    rval = mbImpl->tag_get_handle(tag_name.c_str(), varAttLen.size(), MB_TYPE_INTEGER, varAttLenTag, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<var_name>_ATTRIBS_LEN tag.");
    rval = mbImpl->tag_set_data(varAttLenTag, &file_set, 1, &varAttLen[0]);
    ERRORR(rval, "Trouble setting data for __<var_name>_ATTRIBS_LEN tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // <__VAR_NAMES_LOCATIONS>
  tag_name = "__VAR_NAMES_LOCATIONS";
  Tag varNamesLocsTag = 0;
  std::vector<int> varNamesLocs(varInfo.size());
  rval = mbImpl->tag_get_handle(tag_name.c_str(), varNamesLocs.size(), MB_TYPE_INTEGER, varNamesLocsTag, MB_TAG_CREAT
      | MB_TAG_SPARSE);
  ERRORR(rval, "Trouble creating __VAR_NAMES_LOCATIONS tag.");
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    varNamesLocs[std::distance(varInfo.begin(), mapIter)] = mapIter->second.entLoc;
  }
  rval = mbImpl->tag_set_data(varNamesLocsTag, &file_set, 1, &varNamesLocs[0]);
  ERRORR(rval, "Trouble setting data for __VAR_NAMES_LOCATIONS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__MESH_TYPE>
  Tag meshTypeTag = 0;
  tag_name = "__MESH_TYPE";
  std::string meshTypeName = myHelper->get_mesh_type_name();

  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, meshTypeTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __MESH_TYPE tag.");
  ptr = meshTypeName.c_str();
  int leng= meshTypeName.size();
  rval = mbImpl->tag_set_by_ptr(meshTypeTag, &file_set, 1, &ptr, &leng);
  ERRORR(rval, "Trouble setting data for __MESH_TYPE tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  return MB_SUCCESS;
}

ErrorCode ReadNC::create_attrib_string(const std::map<std::string, AttData>& attMap, std::string& attVal, std::vector<int>& attLen) {
  int success;
  std::stringstream ssAtt;
  unsigned int sz = 0;
  std::map<std::string, AttData>::const_iterator attIt = attMap.begin();
  for (; attIt != attMap.end(); ++attIt) {
    ssAtt << attIt->second.attName;
    ssAtt << '\0';
    void* attData = NULL;
    switch (attIt->second.attDataType) {
      case NC_BYTE:
      case NC_CHAR:
        sz = attIt->second.attLen;
        attData = (char *) malloc(sz);
        success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (char*) attData);
        ERRORS(success, "Failed to read attribute char data.");
        ssAtt << "char;";
        break;
      case NC_DOUBLE:
        sz = attIt->second.attLen * sizeof(double);
        attData = (double *) malloc(sz);
        success = NCFUNC(get_att_double)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (double*) attData);
        ERRORS(success, "Failed to read attribute double data.");
        ssAtt << "double;";
        break;
      case NC_FLOAT:
        sz = attIt->second.attLen * sizeof(float);
        attData = (float *) malloc(sz);
        success = NCFUNC(get_att_float)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (float*) attData);
        ERRORS(success, "Failed to read attribute float data.");
        ssAtt << "float;";
        break;
      case NC_INT:
        sz = attIt->second.attLen * sizeof(int);
        attData = (int *) malloc(sz);
        success = NCFUNC(get_att_int)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (int*) attData);
        ERRORS(success, "Failed to read attribute int data.");
        ssAtt << "int;";
        break;
      case NC_SHORT:
        sz = attIt->second.attLen * sizeof(short);
        attData = (short *) malloc(sz);
        success = NCFUNC(get_att_short)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (short*) attData);
        ERRORS(success, "Failed to read attribute short data.");
        ssAtt << "short;";
        break;
      default:
        success = 1;
    }
    char* tmpc = (char *) attData;
    for (unsigned int counter = 0; counter != sz; ++counter)
      ssAtt << tmpc[counter];
    free(attData);
    ssAtt << ';';
    attLen.push_back(ssAtt.str().size() - 1);
  }
  attVal = ssAtt.str();

  return MB_SUCCESS;
}

ErrorCode ReadNC::create_quad_coordinate_tag(EntityHandle file_set) {
  Range ents;
  ErrorCode rval = mbImpl->get_entities_by_type(file_set, moab::MBQUAD, ents);
  ERRORR(rval, "Trouble getting QUAD entity.");

  std::size_t numOwnedEnts = 0;
#ifdef USE_MPI
  Range ents_owned;
  if (isParallel){
    rval = myPcomm->filter_pstatus(ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &ents_owned);
    ERRORR(rval, "Trouble getting owned QUAD entity.");
    numOwnedEnts = ents_owned.size();
  }
  else
  {
    numOwnedEnts = ents.size();
    ents_owned = ents;
  }
#else
  numOwnedEnts = ents.size();
#endif

  if (numOwnedEnts == 0)
    return MB_SUCCESS;

  assert(numOwnedEnts == ilCVals.size() * jlCVals.size());
  std::vector<double> coords(numOwnedEnts * 3);
  std::size_t pos = 0;
  for (std::size_t j = 0; j != jlCVals.size(); ++j) {
    for (std::size_t i = 0; i != ilCVals.size(); ++i) {
      pos = j * ilCVals.size() * 3 + i * 3;
      coords[pos] = ilCVals[i];
      coords[pos + 1] = jlCVals[j];
      coords[pos + 2] = 0.0;
    }
  }
  std::string tag_name = "COORDS";
  Tag tagh = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating COORDS tag.");

  void *data;
  int count;
#ifdef USE_MPI
  rval = mbImpl->tag_iterate(tagh, ents_owned.begin(), ents_owned.end(), count, data);
#else
  rval = mbImpl->tag_iterate(tagh, ents.begin(), ents.end(), count, data);
#endif
  ERRORR(rval, "Failed to get COORDS tag iterator.");
  assert(count == (int)numOwnedEnts);
  double* quad_data = (double*) data;
  std::copy(coords.begin(), coords.end(), quad_data);
  return MB_SUCCESS;
}

} // namespace moab
