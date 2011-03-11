#include "ReadNC.hpp"

#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <map>

#include "moab/Core.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "FileOptions.hpp"
#include "moab/ScdInterface.hpp"

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {readMeshIface->report_error(str); return rval;}
    
#define ERRORS(err, str) \
    if (err) {readMeshIface->report_error(str); return MB_FAILURE;}
    
namespace moab {

static inline bool strempty( const char* s ) { return !*s; }

ReaderIface* ReadNC::factory( Interface* iface )
  { return new ReadNC( iface ); }

ReadNC::ReadNC(Interface* impl)
        : mbImpl(impl), CPU_WORD_SIZE(-1), IO_WORD_SIZE(-1), fileId(-1), 
          iMin(-1), iMax(-1), jMin(-1), jMax(-1), kMin(-1), kMax(-1), tMin(-1), tMax(-1),
          ilMin(-1), ilMax(-1), jlMin(-1), jlMax(-1), klMin(-1), klMax(-1), 
          iDim(-1), jDim(-1), kDim(-1), tDim(-1), numUnLim(-1), mCurrentMeshHandle(0),
          startVertex(0), startElem(0), mGlobalIdTag(0), 
          max_line_length(-1), max_str_length(-1), vertexOffset(0), dbgOut(stderr),
          isParallel(false)
#ifdef USE_MPI
        , myPcomm(NULL)
#endif
{
  assert(impl != NULL);
  reset();
  
  impl->query_interface(readMeshIface);
}

void ReadNC::reset()
{
  CPU_WORD_SIZE = -1;
  IO_WORD_SIZE = -1;
  fileId = -1;
  iMin = iMax = jMin = jMax = kMin = kMax = tMin = tMax = -1;
  ilMin = ilMax = jlMin = jlMax = klMin = klMax = -1;
  iDim = jDim = kDim = tDim = -1;
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


ReadNC::~ReadNC() 
{
  mbImpl->release_interface(readMeshIface);
}
  
ErrorCode ReadNC::load_file(const char *file_name,
                            const EntityHandle* file_set,
                            const FileOptions& opts,
                            const ReaderIface::SubsetList* subset_list,
                            const Tag* file_id_tag)
{
  ErrorCode rval = MB_SUCCESS;

  //See if opts has variable(s) specified
  std::vector<std::string> var_names;
  std::vector<int> tstep_nums;
  std::vector<double> tstep_vals;

  if (file_id_tag)
    mGlobalIdTag = *file_id_tag;
  else {
      //! get and cache predefined tag handles
    int dum_val = 0;
    rval = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, mGlobalIdTag);
    if (MB_TAG_NOT_FOUND == rval)
      rval = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
                                MB_TYPE_INTEGER, mGlobalIdTag, &dum_val);
  }
  
  int tmpval;
  if (MB_SUCCESS == opts.get_int_option("DEBUG_IO", 1, tmpval)) {
    dbgOut.set_verbosity(tmpval);
    dbgOut.set_prefix("NC ");
  }
  
  bool nomesh = false;
  rval = parse_options(opts, var_names, tstep_nums, tstep_vals, nomesh);
  ERRORR(rval, "Trouble parsing option string.");

  // Open the file
  dbgOut.tprintf(1, "Opening file %s\n", file_name);
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
  
    // Read the header (num dimensions, dimensions, num variables, global attribs)
  rval = read_header();
  ERRORR(rval, " ");

    // make sure there's a file set to put things in
  EntityHandle tmp_set;
  if (!file_set || (file_set && *file_set == 0)) {
    rval = mbImpl->create_meshset(MESHSET_SET, tmp_set);
    ERRORR(rval, "Trouble creating file set.");
  }
  else tmp_set = *file_set;
  
    // Get bounds on ijk space
  rval = init_ijkt_vals(opts);
  ERRORR(rval, "Trouble initializing ijk values.");

    // Create structured mesh vertex/hex sequences
  Range hexes;
  rval = create_verts_hexes(tmp_set, hexes);
  ERRORR(rval, "Trouble creating vertices.");

    // Read variables onto grid
  rval = read_variables(tmp_set, var_names, tstep_nums, nomesh);
  if (MB_FAILURE == rval) return rval;

    // close the file
  success = NCFUNC(close)(fileId);
  ERRORS(success, "Trouble closing file.");

    // create partition set, and populate with elements
#ifdef USE_MPI  
  EntityHandle partn_set;
  rval = mbImpl->create_meshset(MESHSET_SET, partn_set);
  ERRORR(rval, "Trouble creating partition set.");
  myPcomm->partition_sets().insert(partn_set);
  rval = mbImpl->add_entities(partn_set, hexes);
  ERRORR(rval, "Couldn't add new hexes to partition set.");
#endif  

  return MB_SUCCESS;
}

ErrorCode ReadNC::parse_options(const FileOptions &opts,
                                std::vector<std::string> &var_names, 
                                std::vector<int> &tstep_nums,
                                std::vector<double> &tstep_vals,
                                bool &nomesh) 
{
  opts.get_strs_option("VARIABLE", var_names ); 
  opts.get_ints_option("TIMESTEP", tstep_nums); 
  opts.get_reals_option("TIMEVAL", tstep_vals);
  ErrorCode rval = opts.get_null_option("NOMESH");
  if (MB_SUCCESS == rval) nomesh = true;

  if (2 <= dbgOut.get_verbosity()) {
    if (!var_names.empty()) {
      std::cerr << "Variables requested: ";
      for (unsigned int i = 0; i < var_names.size(); i++) std::cerr << var_names[i];
      std::cerr << std::endl;
    }
    if (!tstep_nums.empty()) {
      std::cerr << "Timesteps requested: ";
      for (unsigned int i = 0; i < tstep_nums.size(); i++) std::cerr << tstep_nums[i];
      std::cerr << std::endl;
    }
    if (!tstep_vals.empty()) {
      std::cerr << "Time vals requested: ";
      for (unsigned int i = 0; i < tstep_vals.size(); i++) std::cerr << tstep_vals[i];
      std::cerr << std::endl;
    }
  }
  
#ifdef USE_MPI
  rval = opts.match_option("PARALLEL", "READ_PART");
  isParallel = (rval != MB_ENTITY_NOT_FOUND);

  if (!isParallel) return rval;
  
  int pcomm_no = 0;
  rval = opts.get_int_option("PARALLEL_COMM", pcomm_no);
  if (rval == MB_TYPE_OUT_OF_RANGE) {
    readMeshIface->report_error("Invalid value for PARALLEL_COMM option");
    return rval;
  }
  myPcomm = ParallelComm::get_pcomm(mbImpl, pcomm_no);
  if (0 == myPcomm) {
    myPcomm = new ParallelComm(mbImpl);
  }
  const int rank = myPcomm->proc_config().proc_rank();
  dbgOut.set_rank(rank);
#endif

  return MB_SUCCESS;
}
    
ErrorCode ReadNC::create_verts_hexes(EntityHandle tmp_set, Range &hexes) 
{
    // get the scd interface
  ScdInterface *scdi = NULL;
  ErrorCode rval = mbImpl->query_interface(scdi);
  if (!scdi) return MB_FAILURE;

  Range tmp_range;
  ScdBox *scd_box;
  rval = scdi->create_scd_sequence(HomCoord(ilMin, jlMin, (-1 != klMin ? klMin : 0), 1),
                                   HomCoord(ilMax, jlMax, (-1 != klMax ? klMax : 0), 1),
                                   MBVERTEX, 0, scd_box);
  if (MB_SUCCESS != rval) mbImpl->release_interface(scdi);
  ERRORR(rval, "Trouble creating scd vertex sequence.");

    // add the new vertices to the file set
  tmp_range.insert(scd_box->start_vertex(), scd_box->start_vertex()+scd_box->num_vertices()-1);
  rval = mbImpl->add_entities(tmp_set, tmp_range);
  if (MB_SUCCESS != rval) mbImpl->release_interface(scdi);
  ERRORR(rval, "Couldn't add new vertices to file set.");
  
    // get a ptr to global id memory
  Range::iterator viter = tmp_range.begin();
  void *data;
  rval = mbImpl->tag_iterate(mGlobalIdTag, viter, tmp_range.end(), data);
  if (MB_SUCCESS != rval) mbImpl->release_interface(scdi);
  ERRORR(rval, "Failed to get tag iterator.");
  int *gid_data = (int*)data;

    // set the vertex coordinates
  double *xc, *yc, *zc;
  rval = scd_box->get_coordinate_arrays(xc, yc, zc);
  if (MB_SUCCESS != rval) mbImpl->release_interface(scdi);
  ERRORR(rval, "Couldn't get vertex coordinate arrays.");

  int i, j, k, il, jl, kl;
  int dil = ilMax - ilMin + 1;
  int djl = jlMax - jlMin + 1;
  int di = iMax - iMin + 1;
  int dj = jMax - jMin + 1;
  assert(dil == (int)ilVals.size() && djl == (int)jlVals.size() && 
         (-1 == klMin || klMax-klMin+1 == (int)klVals.size()));
  for (kl = klMin; kl <= klMax; kl++) {
    k = kl - klMin;
    for (jl = jlMin; jl <= jlMax; jl++) {
      j = jl - jlMin;
      for (il = ilMin; il <= ilMax; il++) {
        i = il - ilMin;
        unsigned int pos = i + j*dil + k*dil*djl;
        xc[pos] = ilVals[i];
        yc[pos] = jlVals[j];
        zc[pos] = (-1 == klMin ? 0.0 : klVals[k]);
        *gid_data = (-1 != kl ? kl*di*dj : 0) + jl*di + il + 1;
        gid_data++;
      }
    }
  }

#ifndef NDEBUG
  int num_verts = (ilMax - ilMin + 1) * (jlMax - jlMin + 1) *
    (-1 == klMin ? 1 : klMax-klMin+1);
  std::vector<int> gids(num_verts);
  rval = mbImpl->tag_get_data(mGlobalIdTag, tmp_range, &gids[0]);
  if (MB_SUCCESS != rval) mbImpl->release_interface(scdi);
  ERRORR(rval, "Trouble getting gid values.");
  int vmin = *(std::min_element(gids.begin(), gids.end())),
      vmax = *(std::max_element(gids.begin(), gids.end()));
  dbgOut.tprintf(1, "Vertex gids %d-%d\n", vmin, vmax);
#endif  
    
    // create element sequence
  ScdBox *elem_box;
  rval = scdi->create_scd_sequence(HomCoord(ilMin, jlMin, (-1 != klMin ? klMin : 0), 1),
                                   HomCoord(ilMax, jlMax, (-1 != klMax ? klMax : 0), 1),
                                   (-1 != klMin ? MBHEX : MBQUAD), 0, elem_box);
  mbImpl->release_interface(scdi);
  ERRORR(rval, "Trouble creating scd element sequence.");
  
    // add vertex seq to element seq, forward orientation, unity transform
  rval = elem_box->add_vbox(scd_box,
                              // p1: imin,jmin
                            HomCoord(ilMin, jlMin, (-1 == klMin ? 0 : klMin)),
                            HomCoord(ilMin, jlMin, (-1 == klMin ? 0 : klMin)),
                              // p2: imax,jmin
                            HomCoord(ilMax, jlMin, (-1 == klMin ? 0 : klMin)),
                            HomCoord(ilMax, jlMin, (-1 == klMin ? 0 : klMin)),
                              // p3: imin,jmax
                            HomCoord(ilMin, jlMax, (-1 == klMin ? 0 : klMin)),
                            HomCoord(ilMin, jlMax, (-1 == klMin ? 0 : klMin)));
  ERRORR(rval, "Error constructing structured element sequence.");

    // add the new hexes to the file set
  tmp_range.insert(elem_box->start_element(), elem_box->start_element() + elem_box->num_elements()-1);
  rval = mbImpl->add_entities(tmp_set, tmp_range);
  ERRORR(rval, "Couldn't add new vertices to file set.");
  
  if (2 <= dbgOut.get_verbosity()) {
    assert(elem_box->boundary_complete());
    EntityHandle dum_ent = elem_box->start_element();
    rval = mbImpl->list_entities(&dum_ent, 1);
    ERRORR(rval, "Trouble listing first hex.");
  
    std::vector<EntityHandle> connect;
    rval = mbImpl->get_connectivity(&dum_ent, 1, connect);
    ERRORR(rval, "Trouble getting connectivity.");
  
    rval = mbImpl->list_entities(&connect[0], connect.size());
    ERRORR(rval, "Trouble listing element connectivity.");
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::read_variables(EntityHandle file_set, std::vector<std::string> &var_names,
                                 std::vector<int> &tstep_nums, bool nomesh) 
{
  std::vector<VarData> vdatas;
  std::map<std::string,VarData>::iterator mit;
  ErrorCode rval = MB_SUCCESS;
  
  if (var_names.empty()) {
    for (mit = varInfo.begin(); mit != varInfo.end(); mit++) {
      VarData vd = (*mit).second;
      std::vector<int> tmp_v;
      if (-1 != tMin && 
          std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) 
        tmp_v.push_back(tDim);
      if (-1 != klMin && 
          std::find(vd.varDims.begin(), vd.varDims.end(), kDim) != vd.varDims.end()) 
        tmp_v.push_back(kDim);
      tmp_v.push_back(jDim);
      tmp_v.push_back(iDim);
      if (tmp_v == vd.varDims)
        vdatas.push_back(vd);
    }
  }
  else {
    for (unsigned int i = 0; i < var_names.size(); i++) {
      mit = varInfo.find(var_names[i]);
      if (mit != varInfo.end()) vdatas.push_back((*mit).second);
      else ERRORR(MB_FAILURE, "Couldn't find variable.");
    }
  }
  
  if (tstep_nums.empty() && -1 != tMin)
      // no timesteps input, get them all
    for (int i = tMin; i <= tMax; i++) tstep_nums.push_back(i);
    
  
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = read_variable(file_set, vdatas[i], tstep_nums[t]);
      if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
    }
  }

    // debug output, if requested
  if (1 == dbgOut.get_verbosity()) {
    dbgOut.printf(1, "Read variables: %s", vdatas.begin()->varName.c_str());
    for (unsigned int i = 1; i < vdatas.size(); i++)
      dbgOut.printf(1, ", %s ", vdatas[i].varName.c_str());
    dbgOut.tprintf(1, "\n");
  }
  
  return rval;
}

ErrorCode ReadNC::read_variable(EntityHandle file_set,
                                VarData &var_data, int tstep_num) 
{
    // get the tag to read into
  ErrorCode rval;
  if (!var_data.varTags[tstep_num]) {
    rval = get_tag(var_data, tstep_num, var_data.varTags[tstep_num]);
    ERRORR(rval, "Trouble getting tag.");
  }
  
    // assume point-based values for now?
  if (-1 != tstep_num) {
    if (-1 == tDim || dimVals[tDim] <= tstep_num || tstep_num < 0) {
      ERRORR(MB_INDEX_OUT_OF_RANGE, "Wrong value for timestep number.");
    }
    else if (var_data.varDims[0] != tDim) {
      ERRORR(MB_INDEX_OUT_OF_RANGE, "Non-default timestep number given for time-independent variable.");
    }
  }
  
    // set up the dimensions and counts
  std::vector<NCDF_SIZE> dims, counts;
    // first time
  if (-1 != tstep_num) {
    dims.push_back(tstep_num);
    counts.push_back(1);
  }
    // then z/y/x
  bool have_k = false;
  if (std::find(var_data.varDims.begin(), var_data.varDims.end(), kDim) != var_data.varDims.end()) {
    dims.push_back(klMin); 
    counts.push_back(klMax - klMin + 1);
    have_k = true;
  }
    
  bool have_ij = true;
#ifdef PNETCDF_FILE
    // whether we actually read anything depends on parallel and whether there's a k
  if (!have_k && -1 != kDim && 0 != myPcomm->proc_config().proc_rank())
    have_ij = false;
#endif    

  dims.push_back(jlMin);
  dims.push_back(ilMin);
  counts.push_back(have_ij ? jlMax-jlMin+1 : 0);
  counts.push_back(have_ij ? ilMax-ilMin+1 : 0);

  assert(dims.size() == var_data.varDims.size());
  
    // get vertices in set
  Range verts;
  rval = mbImpl->get_entities_by_dimension(file_set, 0, verts);
  ERRORR(rval, "Trouble getting vertices in set.");
  assert("Should only have a single vertex subrange, since they were read in one shot" &&
         verts.psize() == 1);
  
    // get ptr to tag space
  Range::iterator viter = verts.begin();
  void *data;
  rval = mbImpl->tag_iterate(var_data.varTags[tstep_num], viter, verts.end(), data);
  ERRORR(rval, "Failed to get tag iterator.");
  assert(viter == verts.end());
  
    // finally, read into that space
  int success, *idata;
  double *ddata;
  float *fdata;
  short *sdata;
  
  switch (var_data.varDataType) {
    case NC_BYTE:
    case NC_CHAR:
        success = NCFUNCA(get_vara_text)(fileId, var_data.varId, &dims[0], &counts[0], (char*)data);
        ERRORS(success, "Failed to read char data.");
        break;
    case NC_DOUBLE:
        success = NCFUNCA(get_vara_double)(fileId, var_data.varId, &dims[0], &counts[0], (double*)data);
        ERRORS(success, "Failed to read double data.");
        break;
    case NC_FLOAT:
        success = NCFUNCA(get_vara_float)(fileId, var_data.varId, &dims[0], &counts[0], (float*)data);
        ERRORS(success, "Failed to read float data.");
        ddata = (double*)data;
        fdata = (float*)data;
          // convert in-place
        for (int i = verts.size()-1; i >= 0; i--) 
          ddata[i] = fdata[i];
        break;
    case NC_INT:
        success = NCFUNCA(get_vara_int)(fileId, var_data.varId, &dims[0], &counts[0], (int*)data);
        ERRORS(success, "Failed to read int data.");
        break;
    case NC_SHORT:
        success = NCFUNCA(get_vara_short)(fileId, var_data.varId, &dims[0], &counts[0], (short*)data);
        ERRORS(success, "Failed to read short data.");
        idata = (int*)data;
        sdata = (short*)data;
          // convert in-place
        for (int i = verts.size()-1; i >= 0; i--) 
          idata[i] = sdata[i];
        break;
    default:
        success = 1;
  }
  
  if (2 <= dbgOut.get_verbosity() && !success) {
    double dmin, dmax;
    int imin, imax;
    switch (var_data.varDataType) {
      case NC_DOUBLE:
      case NC_FLOAT:
          ddata = (double*)data;
          if (!verts.empty())
            dmin = ddata[0], dmax = ddata[0];
          for (unsigned int i = 1; i < verts.size(); i++) {
            if (ddata[i] < dmin) dmin = ddata[i];
            if (ddata[i] > dmax) dmax = ddata[i];
          }
          dbgOut.tprintf(2, "Variable %s (double): min = %f, max = %f\n", var_data.varName.c_str(), dmin, dmax);
          break;
      case NC_INT:
      case NC_SHORT:
          idata = (int*)data;
          if (!verts.empty())
            imin = idata[0], imax = idata[0];
          for (unsigned int i = 1; i < verts.size(); i++) {
            if (idata[i] < imin) imin = idata[i];
            if (idata[i] > imax) imax = idata[i];
          }
          dbgOut.tprintf(2, "Variable %s (int): min = %d, max = %d\n", var_data.varName.c_str(), imin, imax);
          break;
      case NC_NAT:
      case NC_BYTE:
      case NC_CHAR:
          break;
    }
  }

  if (success) ERRORR(MB_FAILURE, "Trouble reading variable.");

  return rval;
}
    
ErrorCode ReadNC::get_tag(VarData &var_data, int tstep_num, Tag &tagh) 
{
  std::ostringstream tag_name;
  tag_name << var_data.varName;
  if (0 != tstep_num) tag_name << tstep_num;
  ErrorCode rval = MB_SUCCESS;
  tagh = 0;
  switch (var_data.varDataType) {
    case NC_BYTE:
    case NC_CHAR:
        rval = mbImpl->tag_create(tag_name.str().c_str(), 1, MB_TAG_DENSE, MB_TYPE_OPAQUE, tagh, NULL, true);
        break;
    case NC_DOUBLE:
    case NC_FLOAT:
        rval = mbImpl->tag_create(tag_name.str().c_str(), sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, tagh, NULL, true);
        break;
    case NC_INT:
    case NC_SHORT:
        rval = mbImpl->tag_create(tag_name.str().c_str(), sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tagh, NULL, true);
        break;
    default:
        std::cerr << "Unrecognized data type for tag " << tag_name << std::endl;
        rval = MB_FAILURE;
  }
  
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.str().c_str());

  return rval;
}

ErrorCode ReadNC::init_ijkt_vals(const FileOptions &opts) 
{
    // look for names of i/j/k dimensions
  iMin = iMax = -1;
  std::vector<std::string>::iterator vit;
  unsigned int idx;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lon")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "x1")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "x")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find i variable.");
  iDim = idx;
  iMax = dimVals[idx]-1;
  iMin = 0;
  iName = dimNames[idx];

  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lat")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "y1")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "y")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find j variable.");
  jDim = idx;
  jMax = dimVals[idx]-1;
  jMin = 0;
  jName = dimNames[idx];
  
  kMin = kMax = -1;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lev")) != dimNames.end()) {
    idx = vit-dimNames.begin();
    kMax = dimVals[idx]-1, kMin = 0, kName = std::string("lev");
    kDim = idx;
  }
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "z")) != dimNames.end()) {
    idx = vit-dimNames.begin();
    kMax = dimVals[idx]-1, kMin = 0, kName = std::string("z");
    kDim = idx;
  }

  if ((vit = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "t")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find time variable.");
  tDim = idx;
  tMax = dimVals[idx]-1;
  tMin = 0;
  tName = dimNames[idx];

    // initialize parameter bounds
  ilMin = iMin; ilMax = iMax;
  jlMin = jMin; jlMax = jMax;
  klMin = kMin; klMax = kMax;
  
    // parse options to get subset
#ifdef USE_MPI
  if (isParallel) {
      // partition *the elements* over the parametric space; 1d partition for now, in the k parameter 
      // (or j if there isn't one)
    if (-1 != klMin) {
      int dk = (kMax - kMin) / myPcomm->proc_config().proc_size();
      unsigned int extra = (kMax - kMin) % myPcomm->proc_config().proc_size();
      klMin = kMin + myPcomm->proc_config().proc_rank()*dk + 
          std::min(myPcomm->proc_config().proc_rank(), extra);
      klMax = klMin + dk + (myPcomm->proc_config().proc_rank() < extra ? 1 : 0);

      ilMin = iMin; ilMax = iMax;
      jlMin = jMin; jlMax = jMax;
    }
    else {
      int dj = (jMax - jMin) / myPcomm->proc_config().proc_size();
      unsigned int extra = (jMax - jMin) % myPcomm->proc_config().proc_size();
      jlMin = jMin + myPcomm->proc_config().proc_rank()*dj + 
          std::min(myPcomm->proc_config().proc_rank(), extra);
      jlMax = jlMin + dj + (myPcomm->proc_config().proc_rank() < extra ? 1 : 0);

      ilMin = iMin; ilMax = iMax;
    }
  }
#endif
    
  opts.get_int_option("IMIN", ilMin);
  opts.get_int_option("IMAX", ilMax);
  opts.get_int_option("JMIN", jlMin);
  opts.get_int_option("JMAX", jlMax);

  if (-1 != kMin) {
    opts.get_int_option("KMIN", klMin);
    opts.get_int_option("KMAX", klMax);
  }
  
    // now get actual coordinate values for these dimensions
    // first allocate space...
  if (-1 != ilMin) ilVals.resize(ilMax - ilMin + 1);
  if (-1 != jlMin) jlVals.resize(jlMax - jlMin + 1);
  if (-1 != klMin) klVals.resize(klMax - klMin + 1);
  if (-1 != tMin) tVals.resize(tMax - tMin + 1);

    // ... then read actual values
  std::map<std::string,VarData>::iterator vmit;
  ErrorCode rval;
  if (ilMin != -1) {
    if ((vmit = varInfo.find(iName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(iName.c_str(), ilMin, ilMax, ilVals);
      ERRORR(rval, "Trouble reading x variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find x coordinate.");
    }
  }
  
  if (jlMin != -1) {
    if ((vmit = varInfo.find(jName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(jName.c_str(), jlMin, jlMax, jlVals);
      ERRORR(rval, "Trouble reading y variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find y coordinate.");
    }
  }
  
  if (klMin != -1) {
    if ((vmit = varInfo.find(kName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(kName.c_str(), klMin, klMax, klVals);
      ERRORR(rval, "Trouble reading z variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find z coordinate.");
    }
  }
  
  if (tMin != -1) {
    if ((vmit = varInfo.find(tName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(tName.c_str(), tMin, tMax, tVals);
      ERRORR(rval, "Trouble reading time variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find time coordinate.");
    }
  }

  dbgOut.tprintf(1, "I=%d-%d, J=%d-%d, K=%d-%d\n", ilMin, ilMax, jlMin, jlMax, klMin, klMax);
  dbgOut.tprintf(1, "%d elements, %d vertices\n", (ilMax-ilMin)*(jlMax-jlMin)*(klMax-klMin),
                 (ilMax-ilMin+1)*(jlMax-jlMin+1)*(klMax-klMin+1));
  
  return MB_SUCCESS;
}
    
ErrorCode ReadNC::read_coordinate(const char *var_name, int lmin, int lmax,
                                  std::vector<double> &cvals) 
{
  std::map<std::string,VarData>::iterator vmit = varInfo.find(var_name);
  if (varInfo.end() == vmit) return MB_FAILURE;
  
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
  else
    ERRORR(MB_FAILURE, "Wrong data type for coordinate variable.");

  return MB_SUCCESS;
}
      
ErrorCode ReadNC::read_header()
{
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
  dbgOut.tprintf(1, "Read %u attributes\n", (unsigned int)globalAtts.size());

    // read in dimensions into dimVals
  result = get_dimensions();
  ERRORR(result, "Getting dimensions.");
  dbgOut.tprintf(1, "Read %u dimensions\n", (unsigned int)dimVals.size());

    // read in variables into varInfo
  result = get_variables();
  ERRORR(result, "Getting variables.");
  dbgOut.tprintf(1, "Read %u variables\n", (unsigned int)varInfo.size());

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_attributes(int var_id, int num_atts, std::map<std::string,AttData> &atts,
                                 const char *prefix) 
{

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

    dbgOut.tprintf(2, "%sAttribute %s: length=%u, varId=%d, type=%d\n",
                   (prefix ? prefix : ""), data.attName.c_str(), (unsigned int)data.attLen, 
                   data.attVarId, data.attDataType);
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::get_dimensions()
{
    // get the number of dimensions
  int num_dims;
  int success = NCFUNC(inq_ndims)(fileId, &num_dims);
  ERRORS(success, "Trouble getting number of dimensions.");

  if (num_dims > NC_MAX_DIMS) {
    readMeshIface->report_error("ReadNC: File contains %d dims but NetCDF library supports only %d\n",
                                num_dims, (int)NC_MAX_DIMS);
    return MB_FAILURE;
  }

  char dim_name[NC_MAX_NAME+1];
  NCDF_SIZE dum_len;
  dimNames.resize(num_dims);
  dimVals.resize(num_dims);
  
  for (int i = 0; i < num_dims; i++) {
    success = NCFUNC(inq_dim)(fileId, i, dim_name, &dum_len);
    ERRORS(success, "Trouble getting dimension info.");
    
    dimVals[i] = dum_len;
    dimNames[i] = std::string(dim_name);

    dbgOut.tprintf(2, "Dimension %s, length=%u\n",
                   dim_name, (unsigned int)dum_len);
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::get_variables() 
{
    // first cache the number of time steps
  std::vector<std::string>::iterator vit = std::find(dimNames.begin(), dimNames.end(), "time");
  if (vit == dimNames.end())
    vit = std::find(dimNames.begin(), dimNames.end(), "t");

  int ntimes = 0;
  if (vit != dimNames.end()) dimVals[vit-dimNames.begin()];
  if (!ntimes) ntimes = 1;
  
    // get the number of variables
  int num_vars;
  int success = NCFUNC(inq_nvars)(fileId, &num_vars);
  ERRORS(success, "Trouble getting number of variables.");

  if (num_vars > NC_MAX_VARS) {
    readMeshIface->report_error("ReadNC: File contains %d vars but NetCDF library supports only %d\n",
                                num_vars, (int)NC_MAX_VARS);
    return MB_FAILURE;
  }
  
  char var_name[NC_MAX_NAME+1];
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
    dbgOut.tprintf(2, "Variable %s: Id=%d, numAtts=%d, datatype=%d, num_dims=%d\n",
                   data.varName.c_str(), data.varId, data.numAtts, data.varDataType, 
                   (unsigned int)data.varDims.size());

    ErrorCode rval = get_attributes(i, data.numAtts, data.varAtts, "   ");
    ERRORR(rval, "Trouble getting attributes for a variable.");

  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::read_tag_values( const char* file_name,
                                   const char* tag_name,
                                   const FileOptions& opts,
                                   std::vector<int>& tag_values_out,
                                   const SubsetList* subset_list) 
{
  return MB_FAILURE;
}

} // namespace moab
