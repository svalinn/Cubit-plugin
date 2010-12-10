#include "ReadNC.hpp"
#include "netcdf.h"

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
#include "SequenceManager.hpp"
#include "StructuredElementSeq.hpp"
#include "VertexSequence.hpp"
#include "FileOptions.hpp"

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
          max_line_length(-1), max_str_length(-1), vertexOffset(0), dbgOut(stderr)
{
  assert(impl != NULL);
  reset();
  
  void* ptr = 0;
  impl->query_interface( "ReadUtilIface", &ptr );
  readMeshIface = reinterpret_cast<ReadUtilIface*>(ptr);

  //! get and cache predefined tag handles
  int dum_val = 0;
  ErrorCode result = impl->tag_get_handle(GLOBAL_ID_TAG_NAME, mGlobalIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER, mGlobalIdTag, &dum_val);
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
}


ReadNC::~ReadNC() 
{
  std::string iface_name = "ReadUtilIface";
  mbImpl->release_interface(iface_name, readMeshIface);
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
  int success = nc_open(file_name, 0, &fileId);
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
  rval = create_verts_hexes(tmp_set);
  ERRORR(rval, "Trouble creating vertices.");

    // Read variables onto grid
  rval = read_variables(tmp_set, var_names, tstep_nums, nomesh);
  if (MB_FAILURE == rval) return rval;

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
  
  return MB_SUCCESS;
}
    
ErrorCode ReadNC::create_verts_hexes(EntityHandle tmp_set) 
{
  Core *tmpImpl = dynamic_cast<Core*>(mbImpl);
  VertexSequence *vert_seq;
  EntitySequence *dum_seq;

    // get the seq manager from gMB
  SequenceManager *seq_mgr = tmpImpl->sequence_manager();

  Range tmp_range;
  ErrorCode rval = seq_mgr->create_scd_sequence(ilMin, jlMin, (-1 != klMin ? klMin : 0),
                                                ilMax, jlMax, (-1 != klMax ? klMax : 0),
                                                MBVERTEX, (moab::EntityID)0, startVertex, dum_seq);
  ERRORR(rval, "Trouble creating scd vertex sequence.");

    // add the new vertices to the file set
  tmp_range.insert(startVertex, startVertex+(ilMax-ilMin+1)*(jlMax-jlMin+1)*
                   (-1 == klMin ? 1 : klMax-klMin+1));
  rval = mbImpl->add_entities(tmp_set, tmp_range);
  ERRORR(rval, "Couldn't add new vertices to file set.");
  
  vert_seq = dynamic_cast<VertexSequence*>(dum_seq);
  
    // set the vertex coordinates
  double *xc, *yc, *zc;
  rval = vert_seq->get_coordinate_arrays(xc, yc, zc);
  ERRORR(rval, "Couldn't get vertex coordinate arrays.");

  int i, j, k, il, jl, kl;
  int di = ilMax - ilMin + 1;
  int dj = jlMax - jlMin + 1;

  assert(di == (int)ilVals.size() && dj == (int)jlVals.size() && 
         (-1 == klMin || klMax-klMin+1 == (int)klVals.size()));
  for (kl = klMin; kl <= klMax; kl++) {
    k = kl - klMin;
    for (jl = jlMin; jl <= jlMax; jl++) {
      j = jl - jlMin;
      for (il = ilMin; il <= ilMax; il++) {
        i = il - ilMin;
        unsigned int pos = i + j*di + k*di*dj;
        xc[pos] = ilVals[i];
        yc[pos] = jlVals[j];
        zc[pos] = (-1 == klMin ? 0.0 : klVals[k]);
      }
    }
  }
    
    // create element sequence
  rval = seq_mgr->create_scd_sequence(ilMin, jlMin, (-1 != klMin ? klMin : 0),
                                      ilMax, jlMax, (-1 != klMax ? klMax : 0),
                                      (-1 != klMin ? MBHEX : MBQUAD), (moab::EntityID)0, startElem, dum_seq);
  ERRORR(rval, "Trouble creating scd element sequence.");
  
  StructuredElementSeq *elem_seq = dynamic_cast<StructuredElementSeq*>(dum_seq);
  assert (MB_FAILURE != rval && dum_seq != NULL && elem_seq != NULL);
  
    // add vertex seq to element seq, forward orientation, unity transform
  ScdVertexData *dum_data = dynamic_cast<ScdVertexData*>(vert_seq->data());
  rval = elem_seq->sdata()->add_vsequence(dum_data,
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
  tmp_range.insert(startElem, startElem+(ilMax-ilMin)*(jlMax-jlMin)*
                   (-1 == klMin ? 1 : klMax-klMin));
  rval = mbImpl->add_entities(tmp_set, tmp_range);
  ERRORR(rval, "Couldn't add new vertices to file set.");
  
  if (2 <= dbgOut.get_verbosity()) {
    assert(elem_seq->boundary_complete());
    rval = mbImpl->list_entities(&startElem, 1);
    ERRORR(rval, "Trouble listing first hex.");
  
    std::vector<EntityHandle> connect;
    rval = mbImpl->get_connectivity(&startElem, 1, connect);
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
      tmp_v.push_back(tDim);
      if (-1 != klMin) tmp_v.push_back(kDim);
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
  
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(1, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = read_variable(file_set, vdatas[i], tstep_nums[t]);
      if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
    }
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
  std::vector<size_t> dims, counts;
    // first time
  if (-1 != tstep_num) {
    dims.push_back(tstep_num);
    counts.push_back(1);
  }
    // then z/y/x
  if (-1 != klMin) {
    dims.push_back(klMin); counts.push_back(klMax-klMin+1);
  }
  dims.push_back(jlMin); counts.push_back(jlMax-jlMin+1);
  dims.push_back(ilMin); counts.push_back(ilMax-ilMin+1);

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
    //assert(viter == verts.end());
  
    // finally, read into that space
  int success, *idata;
  double *ddata;
  float *fdata;
  short *sdata;
  
  switch (var_data.varDataType) {
    case NC_BYTE:
    case NC_CHAR:
        success = nc_get_vara_text(fileId, var_data.varId, &dims[0], &counts[0], (char*)data);
        break;
    case NC_DOUBLE:
        success = nc_get_vara_double(fileId, var_data.varId, &dims[0], &counts[0], (double*)data);
        break;
    case NC_FLOAT:
        success = nc_get_vara_float(fileId, var_data.varId, &dims[0], &counts[0], (float*)data);
        ddata = (double*)data;
        fdata = (float*)data;
          // convert in-place
        for (unsigned int i = verts.size(); i > 0; i--) 
          ddata[i] = fdata[i];
        break;
    case NC_INT:
        success = nc_get_vara_int(fileId, var_data.varId, &dims[0], &counts[0], (int*)data);
        break;
    case NC_SHORT:
        success = nc_get_vara_short(fileId, var_data.varId, &dims[0], &counts[0], (short*)data);
        idata = (int*)data;
        sdata = (short*)data;
          // convert in-place
        for (unsigned int i = verts.size(); i > 0; i--) 
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
          for (unsigned int i = 0; i < verts.size(); i++) {
            if (!i) dmin = ddata[i], dmax = ddata[i];
            else {
              if (ddata[i] < dmin) dmin = ddata[i];
              if (ddata[i] > dmax) dmax = ddata[i];
            }
          }
          dbgOut.tprintf(2, "Variable %s (double): min = %f, max = %f\n", var_data.varName.c_str(), dmin, dmax);
          break;
      case NC_INT:
      case NC_SHORT:
          idata = (int*)data;
          for (unsigned int i = 0; i < verts.size(); i++) {
            if (!i) imin = idata[i], imax = idata[i];
            else {
              if (idata[i] < imin) imin = idata[i];
              if (idata[i] > imax) imax = idata[i];
            }
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
  
  if (MB_SUCCESS == rval) dbgOut.tprintf(1, "Tag created for variable %s\n", tag_name.str().c_str());

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
  else ERRORR(MB_FAILURE, "Couldn't find i variable.");
  iDim = idx;
  iMax = dimVals[idx]-1;
  iMin = 0;
  iName = dimNames[idx];

  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lat")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "y1")) != dimNames.end()) 
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

  if ((vit = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find time variable.");
  tDim = idx;
  tMax = dimVals[idx]-1;
  tMin = 0;
  
    // parse options to get subset
    // FOR NOW: just take whole thing
  ilMin = iMin, ilMax = iMax;
  jlMin = jMin, jlMax = jMax;
  klMin = kMin, klMax = kMax;
  
    // now get actual coordinate values for these dimensions
    // first allocate space...
  if (-1 != ilMin) ilVals.resize(ilMax - ilMin + 1);
  if (-1 != jlMin) jlVals.resize(jlMax - jlMin + 1);
  if (-1 != klMin) klVals.resize(klMax - klMin + 1);
  if (-1 != tMin) tVals.resize(tMax - tMin + 1);

    // ... then read actual values
  ErrorCode rval;
  std::map<std::string,VarData>::iterator vmit;
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
    if ((vmit = varInfo.find("time")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate("time", tMin, tMax, tVals);
      ERRORR(rval, "Trouble reading time variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find time coordinate.");
    }
  }
  
  return MB_SUCCESS;
}
    
ErrorCode ReadNC::read_coordinate(const char *var_name, int lmin, int lmax,
                                  std::vector<double> &cvals) 
{
  std::map<std::string,VarData>::iterator vmit = varInfo.find(var_name);
  if (varInfo.end() == vmit) return MB_FAILURE;
  
    // check to make sure it's a float or double
  int success;
  size_t tmin = lmin, tcount = lmax - lmin + 1;
  ptrdiff_t dum_stride = 1;
  if (NC_DOUBLE == (*vmit).second.varDataType) {
    cvals.resize(tcount);
    success = nc_get_vars_double(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &cvals[0]);
    ERRORS(success, "Failed to get coordinate values.");
  }
  else if (NC_FLOAT == (*vmit).second.varDataType) {
    std::vector<float> tcvals(tcount);
    success = nc_get_vars_float(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &tcvals[0]);
    ERRORS(success, "Failed to get coordinate values.");
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
  int success = nc_inq_natts (fileId, &numgatts);
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
    int success = nc_inq_attname(fileId, var_id, i, dum_name);
    ERRORS(success, "Trouble getting attribute name.");
    
    AttData &data = atts[std::string(dum_name)];
    data.attName = std::string(dum_name);
    success = nc_inq_att(fileId, var_id, dum_name, &data.attDataType, &data.attLen);
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
  int success = nc_inq_ndims(fileId, &num_dims);
  ERRORS(success, "Trouble getting number of dimensions.");

  if (num_dims > NC_MAX_DIMS) {
    readMeshIface->report_error("ReadNC: File contains %d dims but NetCDF library supports only %d\n",
                                num_dims, (int)NC_MAX_DIMS);
    return MB_FAILURE;
  }

  char dim_name[NC_MAX_NAME+1];
  size_t dum_len;
  dimNames.resize(num_dims);
  dimVals.resize(num_dims);
  
  for (int i = 0; i < num_dims; i++) {
    success = nc_inq_dim(fileId, i, dim_name, &dum_len);
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
  assert("Should always be a time variable" && vit != dimNames.end());
  int ntimes = dimVals[vit-dimNames.begin()];
  if (!ntimes) ntimes = 1;
  
    // get the number of variables
  int num_vars;
  int success = nc_inq_nvars(fileId, &num_vars);
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
    success = nc_inq_varname (fileId, i, var_name);
    ERRORS(success, "Trouble getting var name.");
    VarData &data = varInfo[std::string(var_name)];
    data.varName = std::string(var_name);
    data.varId = i;
    data.varTags.resize(ntimes, 0);

      // get the data type
    success = nc_inq_vartype(fileId, i, &data.varDataType);
    ERRORS(success, "Trouble getting variable data type.");

      // get the number of dimensions, then the dimensions
    success = nc_inq_varndims(fileId, i, &var_ndims);
    ERRORS(success, "Trouble getting number of dims of a variable.");
    data.varDims.resize(var_ndims);

    success = nc_inq_vardimid(fileId, i, &data.varDims[0]);
    ERRORS(success, "Trouble getting variable dimensions.");

      // finally, get the number of attributes, then the attributes
    success = nc_inq_varnatts(fileId, i, &data.numAtts);
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
