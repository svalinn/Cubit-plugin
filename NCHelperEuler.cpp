#include "NCHelperEuler.hpp"
#include "moab/ReadUtilIface.hpp"
#include "FileOptions.hpp"

#include <cmath>
#include <sstream>

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {_readNC->readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
    if (err) {_readNC->readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

bool NCHelperEuler::can_read_file(ReadNC* readNC, int fileId)
{
  std::vector<std::string>& dimNames = readNC->dimNames;

  // If dimension names "lon" AND "lat' exist then it could be either the Eulerian Spectral grid or the FV grid
  if ((std::find(dimNames.begin(), dimNames.end(), std::string("lon")) != dimNames.end()) && (std::find(dimNames.begin(),
      dimNames.end(), std::string("lat")) != dimNames.end())) {
    // If dimension names "lon" AND "lat" AND "slon" AND "slat" exist then it should be the FV grid
    if ((std::find(dimNames.begin(), dimNames.end(), std::string("slon")) != dimNames.end()) && (std::find(dimNames.begin(),
        dimNames.end(), std::string("slat")) != dimNames.end()))
      return false;

    // Make sure it is CAM grid
    std::map<std::string, ReadNC::AttData>::iterator attIt = readNC->globalAtts.find("source");
    if (attIt == readNC->globalAtts.end()) {
      readNC->readMeshIface->report_error("%s", "File does not have source global attribute.");
      return false;
    }
    unsigned int sz = attIt->second.attLen;
    std::string att_data;
    att_data.resize(sz + 1);
    att_data[sz] = '\000';
    int success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0]);
    if (success != 0) {
      readNC->readMeshIface->report_error("%s", "Failed to read source global attribute char data.");
      return false;
    }
    if (att_data.find("CAM") == std::string::npos)
      return false;

    return true;
  }

  return false;
}

ErrorCode NCHelperEuler::init_mesh_vals(const FileOptions& opts, EntityHandle file_set)
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::vector<std::string>& dimNames = _readNC->dimNames;
  std::vector<int>& dimVals = _readNC->dimVals;
  std::string& iName = _readNC->iName;
  std::string& jName = _readNC->jName;
  std::string& tName = _readNC->tName;
  std::string& iCName = _readNC->iCName;
  std::string& jCName = _readNC->jCName;
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  int& tMin = _readNC->tMin;
  int& tMax = _readNC->tMax;
  int (&gDims)[6] = _readNC->gDims;
  int (&lDims)[6] = _readNC->lDims;
  int (&lCDims)[6] = _readNC->lCDims;
  int (&gCDims)[6] = _readNC->gCDims;
  std::vector<double>& ilVals = _readNC->ilVals;
  std::vector<double>& jlVals = _readNC->jlVals;
  std::vector<double>& tVals = _readNC->tVals;
  std::vector<double>& ilCVals = _readNC->ilCVals;
  std::vector<double>& jlCVals = _readNC->jlCVals;
  int& tDim = _readNC->tDim;
  int& iCDim = _readNC->iCDim;
  int& jCDim = _readNC->jCDim;
  DebugOutput& dbgOut = _readNC->dbgOut;
  bool& isParallel = _readNC->isParallel;
  int& partMethod = _readNC->partMethod;
  int (&locallyPeriodic)[2] = _readNC->locallyPeriodic;
  int (&globallyPeriodic)[2] = _readNC->globallyPeriodic;
  ScdParData& parData = _readNC->parData;
#ifdef USE_MPI
  ParallelComm*& myPcomm = _readNC->myPcomm;
#endif

  // look for names of center i/j dimensions
  std::vector<std::string>::iterator vit;
  unsigned int idx;
  iCName = std::string("lon");
  iName = std::string("slon");
  if ((vit = std::find(dimNames.begin(), dimNames.end(), iCName.c_str())) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find center i variable.");
  }
  iCDim = idx;

  // decide on i periodicity using math for now
  std::vector<double> tilVals(dimVals[idx]);
  ErrorCode rval = _readNC->read_coordinate(iCName.c_str(), 0, dimVals[idx] - 1, tilVals);
  ERRORR(rval, "Trouble reading lon variable.");
  if (std::fabs(2 * (*(tilVals.rbegin())) - *(tilVals.rbegin() + 1) - 360) < 0.001)
    globallyPeriodic[0] = 1;

  // now we can set gCDims and gDims for i
  gCDims[0] = 0;
  gDims[0] = 0;
  gCDims[3] = dimVals[idx] - 1; // these are stored directly in file
  gDims[3] = gCDims[3] + (globallyPeriodic[0] ? 0 : 1); // only if not periodic is vertex param max > elem param max

  // now j
  jCName = std::string("lat");
  jName = std::string("slat");
  if ((vit = std::find(dimNames.begin(), dimNames.end(), jCName.c_str())) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find center j variable.");
  }
  jCDim = idx;

  // for Eul models, will always be non-periodic in j
  gCDims[1] = 0;
  gDims[1] = 0;
  gCDims[4] = dimVals[idx] - 1;
  gDims[4] = gCDims[4] + 1;

  // try a truly 2d mesh
  gDims[2] = -1;
  gDims[5] = -1;

  // look for time dimensions
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end())
    idx = vit - dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "t")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find time variable.");
  }
  tDim = idx;
  tMax = dimVals[idx] - 1;
  tMin = 0;
  tName = dimNames[idx];

  // parse options to get subset
  if (isParallel) {
#ifdef USE_MPI
    for (int i = 0; i < 6; i++)
      parData.gDims[i] = gDims[i];
    for (int i = 0; i < 3; i++)
      parData.gPeriodic[i] = globallyPeriodic[i];
    parData.partMethod = partMethod;
    int pdims[3];

    rval = ScdInterface::compute_partition(myPcomm->proc_config().proc_size(),
        myPcomm->proc_config().proc_rank(),
        parData, lDims, locallyPeriodic, pdims);
    if (MB_SUCCESS != rval)
      return rval;
    for (int i = 0; i < 3; i++)
      parData.pDims[i] = pdims[i];

    dbgOut.tprintf(1, "Partition: %dx%d (out of %dx%d)\n",
        lDims[3] - lDims[0] + 1, lDims[4] - lDims[1] + 1,
        gDims[3] - gDims[0] + 1, gDims[4] - gDims[1] + 1);
    if (myPcomm->proc_config().proc_rank() == 0)
      dbgOut.tprintf(1, "Contiguous chunks of size %d bytes.\n", 8 * (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1));
#endif
  }
  else {
    for (int i = 0; i < 6; i++)
      lDims[i] = gDims[i];
    locallyPeriodic[0] = globallyPeriodic[0];
  }

  opts.get_int_option("IMIN", lDims[0]);
  opts.get_int_option("IMAX", lDims[3]);
  opts.get_int_option("JMIN", lDims[1]);
  opts.get_int_option("JMAX", lDims[4]);

  // now get actual coordinate values for vertices and cell centers; first resize
  if (locallyPeriodic[0]) {
    // if locally periodic, doesn't matter what global periodicity is, # vertex coords = # elem coords
    ilVals.resize(lDims[3] - lDims[0] + 1);
    ilCVals.resize(lDims[3] - lDims[0] + 1);
    lCDims[3] = lDims[3];
  }
  else {
    if (!locallyPeriodic[0] && globallyPeriodic[0] && lDims[3] > gDims[3]) {
      // globally periodic and I'm the last proc, get fewer vertex coords than vertices in i
      ilVals.resize(lDims[3] - lDims[0] + 1);
      ilCVals.resize(lDims[3] - lDims[0]);
      lCDims[3] = lDims[3] - 1;
    }
    else {
      ilVals.resize(lDims[3] - lDims[0] + 1);
      ilCVals.resize(lDims[3] - lDims[0]);
      lCDims[3] = lDims[3] - 1;
    }
  }

  lCDims[0] = lDims[0];
  lCDims[4] = lDims[4] - 1;
  lCDims[1] = lDims[1];

  if (-1 != lDims[1]) {
    jlVals.resize(lDims[4] - lDims[1] + 1);
    jlCVals.resize(lCDims[4] - lCDims[1] + 1);
  }

  if (-1 != tMin)
    tVals.resize(tMax - tMin + 1);

  // now read coord values
  std::map<std::string, ReadNC::VarData>::iterator vmit;
  if (!ilCVals.empty()) {
    if ((vmit = varInfo.find(iCName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = _readNC->read_coordinate(iCName.c_str(), lDims[0], lDims[0] + ilCVals.size() - 1, ilCVals);
      ERRORR(rval, "Trouble reading lon variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find lon coordinate.");
    }
  }

  if (!jlCVals.empty()) {
    if ((vmit = varInfo.find(jCName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = _readNC->read_coordinate(jCName.c_str(), lDims[1], lDims[1] + jlCVals.size() - 1, jlCVals);
      ERRORR(rval, "Trouble reading lat variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find lat coordinate.");
    }
  }

  if (lDims[0] != -1) {
    if ((vmit = varInfo.find(iCName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      double dif = (ilCVals[1] - ilCVals[0]) / 2;
      std::size_t i;
      for (i = 0; i != ilCVals.size(); i++)
        ilVals[i] = ilCVals[i] - dif;
      // the last one is needed only if not periodic
      if (!locallyPeriodic[0])
        ilVals[i] = ilCVals[i - 1] + dif;
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find x coordinate.");
    }
  }

  if (lDims[1] != -1) {
    if ((vmit = varInfo.find(jCName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      if (!isParallel || ((gDims[4] - gDims[1]) == (lDims[4] - lDims[1]))) {
        std::string gwName("gw");
        std::vector<double> gwVals(lDims[4] - lDims[1] - 1);
        rval = _readNC->read_coordinate(gwName.c_str(), lDims[1], lDims[4] - 2, gwVals);
        ERRORR(rval, "Trouble reading gw variable.");
        // copy the correct piece
        jlVals[0] = -(M_PI / 2) * 180 / M_PI;
        unsigned int i = 0;
        double gwSum = -1;
        for (i = 1; i != gwVals.size() + 1; i++) {
          gwSum = gwSum + gwVals[i - 1];
          jlVals[i] = std::asin(gwSum) * 180 / M_PI;
        }
        jlVals[i] = 90.0; // using value of i after loop exits.
      }
      else {
        std::string gwName("gw");
        double gwSum = 0;

        // If this is the first row
        if (lDims[1] == gDims[1]) {
          std::vector<double> gwVals(lDims[4]);
          rval = _readNC->read_coordinate(gwName.c_str(), 0, lDims[4] - 1, gwVals);
          ERRORR(rval, "Trouble reading gw variable.");
          // copy the correct piece
          jlVals[0] = -(M_PI / 2) * 180 / M_PI;
          gwSum = -1;
          for (std::size_t i = 1; i != jlVals.size(); i++) {
            gwSum = gwSum + gwVals[i - 1];
            jlVals[i] = std::asin(gwSum) * 180 / M_PI;
          }
        }
        // or if it's the last row
        else if (lDims[4] == gDims[4]) {
          std::vector<double> gwVals(lDims[4] - 1);
          rval = _readNC->read_coordinate(gwName.c_str(), 0, lDims[4] - 2, gwVals);
          ERRORR(rval, "Trouble reading gw variable.");
          // copy the correct piece
          gwSum = -1;
          for (int j = 0; j != lDims[1] - 1; j++) {
            gwSum = gwSum + gwVals[j];
          }
          std::size_t i = 0;
          for (; i != jlVals.size() - 1; i++) {
            gwSum = gwSum + gwVals[lDims[1] - 1 + i];
            jlVals[i] = std::asin(gwSum) * 180 / M_PI;
          }
          jlVals[i] = 90.0; // using value of i after loop exits.
        }
        // it's in the middle
        else {
          int start = lDims[1] - 1;
          int end = lDims[4] - 1;
          std::vector<double> gwVals(end);
          rval = _readNC->read_coordinate(gwName.c_str(), 0, end - 1, gwVals);
          ERRORR(rval, "Trouble reading gw variable.");
          gwSum = -1;
          for (int j = 0; j != start - 1; j++) {
            gwSum = gwSum + gwVals[j];
          }
          std::size_t i = 0;
          for (; i != jlVals.size(); i++) {
            gwSum = gwSum + gwVals[start - 1 + i];
            jlVals[i] = std::asin(gwSum) * 180 / M_PI;
          }
        }
      }
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find y coordinate.");
    }
  }

  if (tMin != -1) {
    if ((vmit = varInfo.find(tName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = _readNC->read_coordinate(tName.c_str(), tMin, tMax, tVals);
      ERRORR(rval, "Trouble reading time variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find time coordinate.");
    }
  }

  dbgOut.tprintf(1, "I=%d-%d, J=%d-%d\n", lDims[0], lDims[3], lDims[1], lDims[4]);
  dbgOut.tprintf(1, "%d elements, %d vertices\n", (lDims[3] - lDims[0]) * (lDims[4] - lDims[1]), (lDims[3] - lDims[0] + 1)
      * (lDims[4] - lDims[1] + 1));

  // determine the entity location type of a variable
  std::map<std::string, ReadNC::VarData>::iterator mit;
  for (mit = varInfo.begin(); mit != varInfo.end(); ++mit) {
    ReadNC::VarData& vd = (*mit).second;
    if ((std::find(vd.varDims.begin(), vd.varDims.end(), iCDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
        vd.varDims.end(), jCDim) != vd.varDims.end()))
      vd.entLoc = ReadNC::ENTLOCFACE;
  }

  // <coordinate_dim_name>
  std::vector<std::string> ijdimNames(4);
  ijdimNames[0] = "__slon";
  ijdimNames[1] = "__slat";
  ijdimNames[2] = "__lon";
  ijdimNames[3] = "__lat";

  std::string tag_name;
  int val_len = 0;
  for (unsigned int i = 0; i != ijdimNames.size(); i++) {
    tag_name = ijdimNames[i];
    void * val = NULL;
    if (tag_name == "__slon") {
      val = &ilVals[0];
      val_len = ilVals.size();
    }
    else if (tag_name == "__slat") {
      val = &jlVals[0];
      val_len = jlVals.size();
    }
    else if (tag_name == "__lon") {
      val = &ilCVals[0];
      val_len = ilCVals.size();
    }
    else if (tag_name == "__lat") {
      val = &jlCVals[0];
      val_len = jlCVals.size();
    }
    Tag tagh = 0;
    DataType data_type;

    // assume all has same data type as lon
    switch (varInfo["lon"].varDataType) {
      case NC_BYTE:
      case NC_CHAR:
      case NC_DOUBLE:
        data_type = MB_TYPE_DOUBLE;
        break;
      case NC_FLOAT:
        data_type = MB_TYPE_DOUBLE;
        break;
      case NC_INT:
        data_type = MB_TYPE_INTEGER;
        break;
      case NC_SHORT:
      default:
        std::cerr << "Unrecognized data type for tag " << tag_name << std::endl;
        ERRORR(MB_FAILURE, "Unrecognized data type");
        break;
    }
    rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, data_type, tagh, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
    ERRORR(rval, "Trouble creating <coordinate_dim_name> tag.");
    rval = mbImpl->tag_set_by_ptr(tagh, &file_set, 1, &val, &val_len);
    ERRORR(rval, "Trouble setting data for <coordinate_dim_name> tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // __<coordinate_dim_name>_LOC_MINMAX
  for (unsigned int i = 0; i != ijdimNames.size(); i++) {
    std::stringstream ss_tag_name;
    ss_tag_name << ijdimNames[i] << "_LOC_MINMAX";
    tag_name = ss_tag_name.str();
    Tag tagh = 0;
    std::vector<int> val(2, 0);
    if (ijdimNames[i] == "__slon") {
      val[0] = lDims[0];
      val[1] = lDims[3];
    }
    else if (ijdimNames[i] == "__slat") {
      val[0] = lDims[1];
      val[1] = lDims[4];
    }
    else if (ijdimNames[i] == "__lon") {
      val[0] = lCDims[0];
      val[1] = lCDims[3];
    }
    else if (ijdimNames[i] == "__lat") {
      val[0] = lCDims[1];
      val[1] = lCDims[4];
    }
    rval = mbImpl->tag_get_handle(tag_name.c_str(), 2, MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<coordinate_dim_name>_LOC_MINMAX tag.");
    rval = mbImpl->tag_set_data(tagh, &file_set, 1, &val[0]);
    ERRORR(rval, "Trouble setting data for __<coordinate_dim_name>_LOC_MINMAX tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // __<coordinate_dim_name>_GLOBAL_MINMAX
  for (unsigned int i = 0; i != ijdimNames.size(); i++) {
    std::stringstream ss_tag_name;
    ss_tag_name << ijdimNames[i] << "_GLOBAL_MINMAX";
    tag_name = ss_tag_name.str();
    Tag tagh = 0;
    std::vector<int> val(2, 0);
    if (ijdimNames[i] == "__slon") {
      val[0] = gDims[0];
      val[1] = gDims[3];
    }
    else if (ijdimNames[i] == "__slat") {
      val[0] = gDims[1];
      val[1] = gDims[4];
    }
    else if (ijdimNames[i] == "__lon") {
      val[0] = gCDims[0];
      val[1] = gCDims[3];
    }
    else if (ijdimNames[i] == "__lat") {
      val[0] = gCDims[1];
      val[1] = gCDims[4];
    }
    rval = mbImpl->tag_get_handle(tag_name.c_str(), 2, MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<coordinate_dim_name>_GLOBAL_MINMAX tag.");
    rval = mbImpl->tag_set_data(tagh, &file_set, 1, &val[0]);
    ERRORR(rval, "Trouble setting data for __<coordinate_dim_name>_GLOBAL_MINMAX tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // hack: create dummy tags, if needed, for variables like nbnd
  // with no corresponding variables
  _readNC->init_dims_with_no_cvars_info();

  return MB_SUCCESS;
}

ErrorCode NCHelperEuler::create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads)
{
  return _readNC->create_scd_verts_quads(scdi, file_set, quads);
}

} // namespace moab
