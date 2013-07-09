#include "NCHelper.hpp"
#include "NCHelperEuler.hpp"
#include "NCHelperFV.hpp"
#include "NCHelperHOMME.hpp"
#include "NCHelperMPAS.hpp"

#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {_readNC->readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
    if (err) {_readNC->readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

NCHelper* NCHelper::get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts)
{
  // Check if CF convention is being followed
  bool is_CF = false;

  std::map<std::string, ReadNC::AttData>& globalAtts = readNC->globalAtts;
  std::map<std::string, ReadNC::AttData>::iterator attIt = globalAtts.find("conventions");
  if (attIt == globalAtts.end())
    attIt = globalAtts.find("Conventions");

  if (attIt != globalAtts.end()) {
    unsigned int sz = attIt->second.attLen;
    std::string att_data;
    att_data.resize(sz + 1);
    att_data[sz] = '\000';
    int success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0]);
    if (0 == success && att_data.find("CF") != std::string::npos)
      is_CF = true;
  }

  if (is_CF) {
    if (NCHelperEuler::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperEuler(readNC, fileId);
    else if (NCHelperFV::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperFV(readNC, fileId);
    else if (NCHelperHOMME::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperHOMME(readNC, fileId, opts);
  }
  else {
    if (NCHelperMPAS::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperMPAS(readNC, fileId, opts);
  }

  // Unknown NetCDF grid (will fill this in later for POP, CICE and CLM)
  return NULL;
}

ErrorCode NCHelper::read_variable_to_set_allocate(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  std::vector<int>& dimVals = _readNC->dimVals;
  int tDim = _readNC->tDim;
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = MB_SUCCESS;

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    if ((std::find(vdatas[i].varDims.begin(), vdatas[i].varDims.end(), tDim) != vdatas[i].varDims.end()))
      vdatas[i].has_t = true;

    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);

      // get the tag to read into
      if (!vdatas[i].varTags[t]) {
        rval = _readNC->get_tag_to_set(vdatas[i], tstep_nums[t], vdatas[i].varTags[t]);
        ERRORR(rval, "Trouble getting tag.");
      }

      // assume point-based values for now?
      if (-1 == tDim || dimVals[tDim] <= (int) t)
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Wrong value for timestep number.");

      // set up the dimensions and counts
      // first variable dimension is time, if it exists
      if (vdatas[i].has_t)
      {
        if (vdatas[i].varDims.size() != 1)
        {
          vdatas[i].readDims[t].push_back(tstep_nums[t]);
          vdatas[i].readCounts[t].push_back(1);
        }
        else
        {
          vdatas[i].readDims[t].push_back(0);
          vdatas[i].readCounts[t].push_back(tstep_nums.size());
        }
      }

      // Set up other dimensions and counts
      if (vdatas[i].varDims.empty()) {
        // Scalar variable
        vdatas[i].readDims[t].push_back(0);
        vdatas[i].readCounts[t].push_back(1);
      }
      else {
        for (unsigned int idx = 0; idx != vdatas[i].varDims.size(); idx++){
          if (tDim != vdatas[i].varDims[idx]){
            // Push other variable dimensions, except time, which was already pushed
            vdatas[i].readDims[t].push_back(0);
            vdatas[i].readCounts[t].push_back(dimVals[vdatas[i].varDims[idx]]);
          }
        }
      }
      std::size_t sz = 1;
      for (std::size_t idx = 0; idx != vdatas[i].readCounts[t].size(); idx++)
        sz *= vdatas[i].readCounts[t][idx];
      vdatas[i].sz = sz;
      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR:
          vdatas[i].varDatas[t] = new char[sz];
          break;
        case NC_DOUBLE:
        case NC_FLOAT:
          vdatas[i].varDatas[t] = new double[sz];
          break;
        case NC_INT:
        case NC_SHORT:
          vdatas[i].varDatas[t] = new int[sz];
          break;
        default:
          std::cerr << "Unrecognized data type for tag " << std::endl;
          rval = MB_FAILURE;
      }
      if (vdatas[i].varDims.size() <= 1)
        break;
    }
  }

  return rval;
}

ErrorCode NCHelper::read_variable_to_set(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  std::set<std::string>& dummyVarNames = _readNC->dummyVarNames;;
  Interface*& mbImpl = _readNC->mbImpl;
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_variable_to_set_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables to set.");

  // finally, read into that space
  int success;
  std::vector<int> requests(vdatas.size() * tstep_nums.size()), statuss(vdatas.size() * tstep_nums.size());
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    if (dummyVarNames.find(vdatas[i].varName) != dummyVarNames.end() )
       continue;// this is a dummy one, we don't have it; we created it for the dummy tag
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      void* data = vdatas[i].varDatas[t];

      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR:
          success = NCFUNCAG(_vara_text)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (char*) data NCREQ);
          ERRORS(success, "Failed to read char data.");
          break;
        case NC_DOUBLE:
          success = NCFUNCAG(_vara_double)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (double*) data NCREQ);
          ERRORS(success, "Failed to read double data.");
          break;
        case NC_FLOAT: {
          success = NCFUNCAG(_vara_float)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (float*) data NCREQ);
          ERRORS(success, "Failed to read float data.");
          break;
        }
        case NC_INT:
          success = NCFUNCAG(_vara_int)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (int*) data NCREQ);
          ERRORS(success, "Failed to read int data.");
          break;
        case NC_SHORT:
          success = NCFUNCAG(_vara_short)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              (short*) data NCREQ);
          ERRORS(success, "Failed to read short data.");
          break;
        default:
          success = 1;
      }

      if (success)
        ERRORR(MB_FAILURE, "Trouble reading variable.");
      if (vdatas[i].varDims.size() <= 1)
        break;
    }
  }

#ifdef NCWAIT
  int success = ncmpi_wait_all(fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Converting variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = convert_variable(vdatas[i], t);
      if (MB_SUCCESS != tmp_rval)
        rval = tmp_rval;
      if (vdatas[i].varDims.size() <= 1)
        break;
    }
  }
  // debug output, if requested
  if (1 == dbgOut.get_verbosity()) {
    dbgOut.printf(1, "Read variables: %s", vdatas.begin()->varName.c_str());
    for (unsigned int i = 1; i < vdatas.size(); i++)
      dbgOut.printf(1, ", %s ", vdatas[i].varName.c_str());
    dbgOut.tprintf(1, "\n");
  }

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Setting data for variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = mbImpl->tag_set_by_ptr(vdatas[i].varTags[t], &file_set, 1, &(vdatas[i].varDatas[t]), &vdatas[i].sz);
      if (MB_SUCCESS != tmp_rval)
        rval = tmp_rval;
      if (vdatas[i].varDims.size() <= 1)
        break;
    }
  }

  return rval;
}

ErrorCode NCHelper::convert_variable(ReadNC::VarData& var_data, int tstep_num)
{
  DebugOutput& dbgOut = _readNC->dbgOut;

  // Get ptr to tag space
  void* data = var_data.varDatas[tstep_num];

  std::size_t sz = 1;
  for (std::size_t idx = 0; idx != var_data.readCounts[tstep_num].size(); idx++)
    sz *= var_data.readCounts[tstep_num][idx];

  // Finally, read into that space
  int success = 0;
  int* idata;
  double* ddata;
  float* fdata;
  short* sdata;

  switch (var_data.varDataType) {
    case NC_FLOAT:
      ddata = (double*) var_data.varDatas[tstep_num];
      fdata = (float*) var_data.varDatas[tstep_num];
      // Convert in-place
      for (int i = sz - 1; i >= 0; i--)
        ddata[i] = fdata[i];
      break;
    case NC_SHORT:
      idata = (int*) var_data.varDatas[tstep_num];
      sdata = (short*) var_data.varDatas[tstep_num];
      // Convert in-place
      for (int i = sz - 1; i >= 0; i--)
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
        ddata = (double*) data;
        if (sz == 0)
          break;

        dmin = dmax = ddata[0];
        for (unsigned int i = 1; i < sz; i++) {
          if (ddata[i] < dmin)
            dmin = ddata[i];
          if (ddata[i] > dmax)
            dmax = ddata[i];
        }
        dbgOut.tprintf(2, "Variable %s (double): min = %f, max = %f\n", var_data.varName.c_str(), dmin, dmax);
        break;
      case NC_INT:
      case NC_SHORT:
        idata = (int*) data;
        if (sz == 0)
          break;

        imin = imax = idata[0];
        for (unsigned int i = 1; i < sz; i++) {
          if (idata[i] < imin)
            imin = idata[i];
          if (idata[i] > imax)
            imax = idata[i];
        }
        dbgOut.tprintf(2, "Variable %s (int): min = %d, max = %d\n", var_data.varName.c_str(), imin, imax);
        break;
      case NC_NAT:
      case NC_BYTE:
      case NC_CHAR:
        break;
      default: // Default case added to remove compiler warnings
        success = 1;
    }
  }

  return MB_SUCCESS;
}

ErrorCode ScdNCHelper::check_existing_mesh(EntityHandle file_set) {
  Interface*& mbImpl = _readNC->mbImpl;
  int (&lDims)[6] = _readNC->lDims;

  // Get the number of vertices
  int num_verts;
  ErrorCode rval = mbImpl->get_number_entities_by_dimension(file_set, 0, num_verts);
  ERRORR(rval, "Trouble getting number of vertices.");

  /*
  // Check against parameters
  // When ghosting is used, this check might fail (to be updated later)
  if (num_verts > 0)
  {
    int expected_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2] + 1);
    if (num_verts != expected_verts) {
      ERRORR(MB_FAILURE, "Number of vertices doesn't match.");
    }
  }
  */

  // Check the number of elements too
  int num_elems;
  rval = mbImpl->get_number_entities_by_dimension(file_set, (-1 == lDims[2] ? 2 : 3), num_elems);
  ERRORR(rval, "Trouble getting number of elements.");

  /*
  // Check against parameters
  // The expected number of elements calculated below is incorrect (to be updated later)
  if (num_elems > 0)
  {
    int expected_elems = (lDims[3] - lDims[0]) * (lDims[4] - lDims[1]) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2]);
    if (num_elems != expected_elems) {
      ERRORR(MB_FAILURE, "Number of elements doesn't match.");
    }
  }
  */

  return MB_SUCCESS;
}

ErrorCode ScdNCHelper::create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& faces)
{
  Interface*& mbImpl = _readNC->mbImpl;
  int (&gDims)[6] = _readNC->gDims;
  int (&lDims)[6] = _readNC->lDims;
  std::vector<double>& ilVals = _readNC->ilVals;
  std::vector<double>& jlVals = _readNC->jlVals;
  std::vector<double>& klVals = _readNC->klVals;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  DebugOutput& dbgOut = _readNC->dbgOut;
  int (&locallyPeriodic)[2] = _readNC->locallyPeriodic;
  int (&globallyPeriodic)[2] = _readNC->globallyPeriodic;
  ScdParData& parData = _readNC->parData;

  Range tmp_range;
  ScdBox *scd_box;

  ErrorCode rval = scdi->construct_box(HomCoord(lDims[0], lDims[1], lDims[2], 1), HomCoord(lDims[3], lDims[4], lDims[5], 1), NULL,
      0, scd_box, locallyPeriodic, &parData);
  ERRORR(rval, "Trouble creating scd vertex sequence.");

  // Add box set and new vertices, elements to the file set
  tmp_range.insert(scd_box->start_vertex(), scd_box->start_vertex() + scd_box->num_vertices() - 1);
  tmp_range.insert(scd_box->start_element(), scd_box->start_element() + scd_box->num_elements() - 1);
  tmp_range.insert(scd_box->box_set());
  rval = mbImpl->add_entities(file_set, tmp_range);
  ERRORR(rval, "Couldn't add new vertices to file set.");

  dbgOut.tprintf(1, "scdbox %d quads, %d vertices\n", scd_box->num_elements(), scd_box->num_vertices());

  // Get a ptr to global id memory
  void* data;
  int count;
  const Range::iterator topv = tmp_range.upper_bound(tmp_range.begin(), tmp_range.end(), scd_box->start_vertex()
      + scd_box->num_vertices());
  rval = mbImpl->tag_iterate(mGlobalIdTag, tmp_range.begin(), topv, count, data);
  ERRORR(rval, "Failed to get tag iterator.");
  assert(count == scd_box->num_vertices());
  int* gid_data = (int*) data;

  // Set the vertex coordinates
  double *xc, *yc, *zc;
  rval = scd_box->get_coordinate_arrays(xc, yc, zc);
  ERRORR(rval, "Couldn't get vertex coordinate arrays.");

  int i, j, k, il, jl, kl, itmp;
  int dil = lDims[3] - lDims[0] + 1;
  int djl = lDims[4] - lDims[1] + 1;
  int di = gDims[3] - gDims[0] + 1;
  int dj = gDims[4] - gDims[1] + 1;
  assert(dil == (int)ilVals.size() && djl == (int)jlVals.size() &&
      (-1 == lDims[2] || lDims[5]-lDims[2]+1 == (int)klVals.size()));
#define INDEX(i, j, k) ()
  for (kl = lDims[2]; kl <= lDims[5]; kl++) {
    k = kl - lDims[2];
    for (jl = lDims[1]; jl <= lDims[4]; jl++) {
      j = jl - lDims[1];
      for (il = lDims[0]; il <= lDims[3]; il++) {
        i = il - lDims[0];
        unsigned int pos = i + j * dil + k * dil * djl;
        xc[pos] = ilVals[i];
        yc[pos] = jlVals[j];
        zc[pos] = (-1 == lDims[2] ? 0.0 : klVals[k]);
        itmp = (!locallyPeriodic[0] && globallyPeriodic[0] && il == gDims[3] ? gDims[0] : il);
        *gid_data = (-1 != kl ? kl * di * dj : 0) + jl * di + itmp + 1;
        gid_data++;
      }
    }
  }
#undef INDEX

#ifndef NDEBUG
  int num_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2] + 1);
  std::vector<int> gids(num_verts);
  Range verts(scd_box->start_vertex(), scd_box->start_vertex() + scd_box->num_vertices() - 1);
  rval = mbImpl->tag_get_data(mGlobalIdTag, verts, &gids[0]);
  ERRORR(rval, "Trouble getting gid values.");
  int vmin = *(std::min_element(gids.begin(), gids.end())), vmax = *(std::max_element(gids.begin(), gids.end()));
  dbgOut.tprintf(1, "Vertex gids %d-%d\n", vmin, vmax);
#endif

  // add elements to the range passed in
  faces.insert(scd_box->start_element(), scd_box->start_element() + scd_box->num_elements() - 1);

  if (2 <= dbgOut.get_verbosity()) {
    assert(scd_box->boundary_complete());
    EntityHandle dum_ent = scd_box->start_element();
    rval = mbImpl->list_entities(&dum_ent, 1);
    ERRORR(rval, "Trouble listing first hex.");

    std::vector<EntityHandle> connect;
    rval = mbImpl->get_connectivity(&dum_ent, 1, connect);
    ERRORR(rval, "Trouble getting connectivity.");

    rval = mbImpl->list_entities(&connect[0], connect.size());
    ERRORR(rval, "Trouble listing element connectivity.");
  }

  Range edges;
  mbImpl->get_adjacencies(faces, 1, true, edges, Interface::UNION);

  return MB_SUCCESS;
}

ErrorCode ScdNCHelper::read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
{
  std::vector<ReadNC::VarData> vdatas;
  std::vector<ReadNC::VarData> vsetdatas;

  ErrorCode rval = read_scd_variable_setup(var_names, tstep_nums, vdatas, vsetdatas);
  ERRORR(rval, "Trouble setting up read variable.");

  // create COORDS tag for quads
  rval = _readNC->create_quad_coordinate_tag(file_set);
  ERRORR(rval, "Trouble creating coordinate tags to entities quads");

  if (!vsetdatas.empty()) {
    rval = read_variable_to_set(file_set, vsetdatas, tstep_nums);
    ERRORR(rval, "Trouble read variables to set.");
  }

  if (!vdatas.empty()) {
    rval = read_scd_variable_to_nonset(file_set, vdatas, tstep_nums);
    ERRORR(rval, "Trouble read variables to entities verts/edges/faces.");
  }

  return MB_SUCCESS;
}

ErrorCode ScdNCHelper::read_scd_variable_setup(std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
                                               std::vector<ReadNC::VarData>& vdatas, std::vector<ReadNC::VarData>& vsetdatas)
{
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  int& tMin = _readNC->tMin;
  int& tMax = _readNC->tMax;
  int& iDim = _readNC->iDim;
  int& jDim = _readNC->jDim;
  int& tDim = _readNC->tDim;
  int& iCDim = _readNC->iCDim;
  int& jCDim = _readNC->jCDim;

  std::map<std::string, ReadNC::VarData>::iterator mit;

  // If empty read them all
  if (var_names.empty()) {
    for (mit = varInfo.begin(); mit != varInfo.end(); ++mit) {
      ReadNC::VarData vd = (*mit).second;
      if ((std::find(vd.varDims.begin(), vd.varDims.end(), iCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
          vd.varDims.end(), jCDim) != vd.varDims.end()))
        vdatas.push_back(vd);
      else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
         vd.varDims.end(), iCDim) != vd.varDims.end()))
        vdatas.push_back(vd);
      else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
          vd.varDims.end(), iDim) != vd.varDims.end()))
        vdatas.push_back(vd);
      else
        vsetdatas.push_back(vd);
    }
  }
  else {
    for (unsigned int i = 0; i < var_names.size(); i++) {
      mit = varInfo.find(var_names[i]);
      if (mit != varInfo.end()) {
        ReadNC::VarData vd = (*mit).second;
        if ((std::find(vd.varDims.begin(), vd.varDims.end(), iCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
           vd.varDims.end(), jCDim) != vd.varDims.end()))
          vdatas.push_back(vd);
        else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), iCDim) != vd.varDims.end()))
          vdatas.push_back(vd);
        else if ((std::find(vd.varDims.begin(), vd.varDims.end(), jCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), iDim) != vd.varDims.end()))
          vdatas.push_back(vd);
        else
          vsetdatas.push_back(vd);
      }
      else ERRORR(MB_FAILURE, "Couldn't find variable.");
    }
  }

  if (tstep_nums.empty() && -1 != tMin) {
    // no timesteps input, get them all
    for (int i = tMin; i <= tMax; i++)
      tstep_nums.push_back(i);
  }
  if (!tstep_nums.empty()) {
    for (unsigned int i = 0; i < vdatas.size(); i++) {
      vdatas[i].varTags.resize(tstep_nums.size(), 0);
      vdatas[i].varDatas.resize(tstep_nums.size());
      vdatas[i].readDims.resize(tstep_nums.size());
      vdatas[i].readCounts.resize(tstep_nums.size());
    }
    for (unsigned int i = 0; i < vsetdatas.size(); i++) {
      if ((std::find(vsetdatas[i].varDims.begin(), vsetdatas[i].varDims.end(), tDim) != vsetdatas[i].varDims.end())
          && (vsetdatas[i].varDims.size() != 1)) {
        vsetdatas[i].varTags.resize(tstep_nums.size(), 0);
        vsetdatas[i].varDatas.resize(tstep_nums.size());
        vsetdatas[i].readDims.resize(tstep_nums.size());
        vsetdatas[i].readCounts.resize(tstep_nums.size());
      }
      else {
        vsetdatas[i].varTags.resize(1, 0);
        vsetdatas[i].varDatas.resize(1);
        vsetdatas[i].readDims.resize(1);
        vsetdatas[i].readCounts.resize(1);
      }
    }
  }

  return MB_SUCCESS;
}

ErrorCode ScdNCHelper::read_scd_variable_to_nonset_allocate(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::vector<std::string>& dimNames = _readNC->dimNames;
  std::vector<int>& dimVals = _readNC->dimVals;
  int (&lDims)[6] = _readNC->lDims;
  int (&lCDims)[6] = _readNC->lCDims;
  int& tDim = _readNC->tDim;
  DebugOutput& dbgOut = _readNC->dbgOut;
  bool& isParallel = _readNC->isParallel;
 #ifdef USE_MPI
  ParallelComm*& myPcomm = _readNC->myPcomm;
#endif

  ErrorCode rval = MB_SUCCESS;

  Range* range = NULL;

  // get vertices in set
  Range verts;
  rval = mbImpl->get_entities_by_dimension(file_set, 0, verts);
  ERRORR(rval, "Trouble getting vertices in set.");
  assert("Should only have a single vertex subrange, since they were read in one shot" &&
      verts.psize() == 1);

  Range edges;
  rval = mbImpl->get_entities_by_dimension(file_set, 1, edges);
  ERRORR(rval, "Trouble getting edges in set.");

  // Get faces in set
  Range faces;
  rval = mbImpl->get_entities_by_dimension(file_set, 2, faces);
  ERRORR(rval, "Trouble getting faces in set.");
  assert("Should only have a single face subrange, since they were read in one shot" &&
      faces.psize() == 1);

#ifdef USE_MPI
  moab::Range faces_owned;
  if (isParallel)
  {
    rval = myPcomm->filter_pstatus(faces, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &faces_owned);
    ERRORR(rval, "Trouble getting owned faces in set.");
  }
  else
    faces_owned = faces; // not running in parallel, but still with MPI
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);

      std::vector<std::string>::iterator vit;
      int idx_lev = -1;
      int idx_ilev = -1;
      if ((vit = std::find(dimNames.begin(), dimNames.end(), "lev")) != dimNames.end())
        idx_lev = vit - dimNames.begin();
      if ((vit = std::find(dimNames.begin(), dimNames.end(), "ilev")) != dimNames.end())
        idx_ilev = vit - dimNames.begin();
      if (std::find(vdatas[i].varDims.begin(), vdatas[i].varDims.end(), idx_lev) != vdatas[i].varDims.end())
        vdatas[i].numLev = dimVals[idx_lev];
      else if (std::find(vdatas[i].varDims.begin(), vdatas[i].varDims.end(), idx_ilev) != vdatas[i].varDims.end())
        vdatas[i].numLev = dimVals[idx_ilev];

      // Get the tag to read into
      if (!vdatas[i].varTags[t]) {
        rval = _readNC->get_tag_to_nonset(vdatas[i], tstep_nums[t], vdatas[i].varTags[t], vdatas[i].numLev);
        ERRORR(rval, "Trouble getting tag.");
      }

      // Assume point-based values for now?
      if (-1 == tDim || dimVals[tDim] <= (int) t) {
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Wrong value for timestep number.");
      }
      else if (vdatas[i].varDims[0] != tDim) {
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Non-default timestep number given for time-independent variable.");
      }

      // Set up the dimensions and counts
      // First time
      vdatas[i].readDims[t].push_back(tstep_nums[t]);
      vdatas[i].readCounts[t].push_back(1);

      // then z/y/x
      if (vdatas[i].numLev != 1) {
        vdatas[i].readDims[t].push_back(0);
        vdatas[i].readCounts[t].push_back(vdatas[i].numLev);
      }

      switch (vdatas[i].entLoc) {
        case ReadNC::ENTLOCVERT:
          // vertices
          // only structured mesh has j parameter that multiplies i to get total # vertices
          vdatas[i].readDims[t].push_back(lDims[1]);
          vdatas[i].readCounts[t].push_back(lDims[4] - lDims[1] + 1);
          vdatas[i].readDims[t].push_back(lDims[0]);
          vdatas[i].readCounts[t].push_back(lDims[3] - lDims[0] + 1);
          assert(vdatas[i].readDims[t].size() == vdatas[i].varDims.size());
          range = &verts;
          break;
        case ReadNC::ENTLOCNSEDGE:
          ERRORR(MB_FAILURE, "Reading edge data not implemented yet.");
          break;
        case ReadNC::ENTLOCEWEDGE:
          ERRORR(MB_FAILURE, "Reading edge data not implemented yet.");
          break;
        case ReadNC::ENTLOCFACE:
          // faces
          vdatas[i].readDims[t].push_back(lCDims[1]);
          vdatas[i].readDims[t].push_back(lCDims[0]);
          vdatas[i].readCounts[t].push_back(lCDims[4] - lCDims[1] + 1);
          vdatas[i].readCounts[t].push_back(lCDims[3] - lCDims[0] + 1);
          assert(vdatas[i].readDims[t].size() == vdatas[i].varDims.size());
#ifdef USE_MPI
          range = &faces_owned;
#else
          range = &faces;
#endif
          break;
        case ReadNC::ENTLOCSET:
          // set
          break;
        default:
          ERRORR(MB_FAILURE, "Unrecognized entity location type.");
          break;
      }

      // Get ptr to tag space
      void* data;
      int count;
      rval = mbImpl->tag_iterate(vdatas[i].varTags[t], range->begin(), range->end(), count, data);
      ERRORR(rval, "Failed to get tag iterator.");
      assert((unsigned)count == range->size());
      vdatas[i].varDatas[t] = data;
    }
  }

  return rval;
}

ErrorCode ScdNCHelper::read_scd_variable_to_nonset(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_scd_variable_to_nonset_allocate(file_set, vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables.");

  // Finally, read into that space
  int success;
  std::vector<int> requests(vdatas.size() * tstep_nums.size()), statuss(vdatas.size() * tstep_nums.size());
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      std::size_t sz = 1;
      for (std::size_t idx = 0; idx != vdatas[i].readCounts[t].size(); idx++)
        sz *= vdatas[i].readCounts[t][idx];
      void* data = vdatas[i].varDatas[t];
      size_t ni = vdatas[i].readCounts[t][2];
      size_t nj = vdatas[i].readCounts[t][3];
      size_t nk = vdatas[i].readCounts[t][1];

      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR: {
          std::vector<char> tmpchardata(sz);
          success = NCFUNCAG(_vara_text)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpchardata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = _readNC->kji_to_jik(ni, nj, nk, data, &tmpchardata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpchardata.size(); idx++)
              ((char*) data)[idx] = tmpchardata[idx];
          }
          ERRORS(success, "Failed to read char data.");
          break;
        }
        case NC_DOUBLE: {
          std::vector<double> tmpdoubledata(sz);
          success = NCFUNCAG(_vara_double)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpdoubledata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // switch from k varying slowest to k varying fastest
            success = _readNC->kji_to_jik(ni, nj, nk, data, &tmpdoubledata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpdoubledata.size(); idx++)
              ((double*) data)[idx] = tmpdoubledata[idx];
          }
          ERRORS(success, "Failed to read double data.");
          break;
        }
        case NC_FLOAT: {
          std::vector<float> tmpfloatdata(sz);
          success = NCFUNCAG(_vara_float)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpfloatdata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = _readNC->kji_to_jik(ni, nj, nk, data, &tmpfloatdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpfloatdata.size(); idx++)
              ((float*) data)[idx] = tmpfloatdata[idx];
          }
          ERRORS(success, "Failed to read float data.");
          break;
        }
        case NC_INT: {
          std::vector<int> tmpintdata(sz);
          success = NCFUNCAG(_vara_int)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpintdata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = _readNC->kji_to_jik(ni, nj, nk, data, &tmpintdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpintdata.size(); idx++)
              ((int*) data)[idx] = tmpintdata[idx];
          }
          ERRORS(success, "Failed to read int data.");
          break;
        }
        case NC_SHORT: {
          std::vector<short> tmpshortdata(sz);
          success = NCFUNCAG(_vara_short)(_fileId, vdatas[i].varId, &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0],
              &tmpshortdata[0] NCREQ);
          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = _readNC->kji_to_jik(ni, nj, nk, data, &tmpshortdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpshortdata.size(); idx++)
              ((short*) data)[idx] = tmpshortdata[idx];
          }
          ERRORS(success, "Failed to read short data.");
          break;
        }
        default:
          success = 1;
      }

      if (success)
        ERRORR(MB_FAILURE, "Trouble reading variable.");
    }
  }

#ifdef NCWAIT
  int success = ncmpi_wait_all(fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Converting variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = convert_variable(vdatas[i], t);
      if (MB_SUCCESS != tmp_rval)
        rval = tmp_rval;
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

ErrorCode UcdNCHelper::read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
{
  std::vector<ReadNC::VarData> vdatas;
  std::vector<ReadNC::VarData> vsetdatas;

  ErrorCode rval = read_ucd_variable_setup(var_names, tstep_nums, vdatas, vsetdatas);
  ERRORR(rval, "Trouble setting up read variable.");

  if (!vsetdatas.empty()) {
    rval = read_variable_to_set(file_set, vsetdatas, tstep_nums);
    ERRORR(rval, "Trouble read variables to set.");
  }

  if (!vdatas.empty()) {
#ifdef PNETCDF_FILE
    // in serial, we will use the old read, everything is contiguous
    // in parallel, we will use async read in pnetcdf
    // the other mechanism is not working, forget about it
    rval = read_ucd_variable_to_nonset_async(file_set, vdatas, tstep_nums);
#else
    rval = read_ucd_variable_to_nonset(file_set, vdatas, tstep_nums);
#endif
    ERRORR(rval, "Trouble read variables to entities verts/edges/faces.");
  }

  return MB_SUCCESS;
}

} // namespace moab
