#include "NCHelper.hpp"
#include "NCHelperEuler.hpp"
#include "NCHelperFV.hpp"
#include "NCHelperHOMME.hpp"
#include "NCHelperMPAS.hpp"

#include <sstream>

#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) {_readNC->readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
  if (err) {_readNC->readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

NCHelper* NCHelper::get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
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
      return new (std::nothrow) NCHelperEuler(readNC, fileId, opts, fileSet);
    else if (NCHelperFV::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperFV(readNC, fileId, opts, fileSet);
    else if (NCHelperHOMME::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperHOMME(readNC, fileId, opts, fileSet);
  }
  else {
    if (NCHelperMPAS::can_read_file(readNC))
      return new (std::nothrow) NCHelperMPAS(readNC, fileId, opts, fileSet);
  }

  // Unknown NetCDF grid (will fill this in later for POP, CICE and CLM)
  return NULL;
}

ErrorCode NCHelper::create_conventional_tags(const std::vector<int>& tstep_nums) {
  Interface*& mbImpl = _readNC->mbImpl;
  std::vector<std::string>& dimNames = _readNC->dimNames;
  std::vector<int>& dimLens = _readNC->dimLens;
  std::map<std::string, ReadNC::AttData>& globalAtts = _readNC->globalAtts;
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  DebugOutput& dbgOut = _readNC->dbgOut;
  int& partMethod = _readNC->partMethod;
  ScdInterface* scdi = _readNC->scdi;

  ErrorCode rval;
  std::string tag_name;

  // <__NUM_DIMS>
  Tag numDimsTag = 0;
  tag_name = "__NUM_DIMS";
  int numDims = dimNames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 1, MB_TYPE_INTEGER, numDimsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_DIMS tag.");
  rval = mbImpl->tag_set_data(numDimsTag, &_fileSet, 1, &numDims);
  ERRORR(rval, "Trouble setting data for __NUM_DIMS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__NUM_VARS>
  Tag numVarsTag = 0;
  tag_name = "__NUM_VARS";
  int numVars = varInfo.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 1, MB_TYPE_INTEGER, numVarsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_VARS tag.");
  rval = mbImpl->tag_set_data(numVarsTag, &_fileSet, 1, &numVars);
  ERRORR(rval, "Trouble setting data for __NUM_VARS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__DIM_NAMES>
  Tag dimNamesTag = 0;
  tag_name = "__DIM_NAMES";
  std::string dimnames;
  unsigned int dimNamesSz = dimNames.size();
  for (unsigned int i = 0; i != dimNamesSz; i++) {
    dimnames.append(dimNames[i]);
    dimnames.push_back('\0');
  }
  int dimnamesSz = dimnames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, dimNamesTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __DIM_NAMES tag.");
  const void* ptr = dimnames.c_str();
  rval = mbImpl->tag_set_by_ptr(dimNamesTag, &_fileSet, 1, &ptr, &dimnamesSz);
  ERRORR(rval, "Trouble setting data for __DIM_NAMES tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__DIM_LENS>
  Tag dimLensTag = 0;
  tag_name = "__DIM_LENS";
  int dimLensSz = dimLens.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_INTEGER, dimLensTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __DIM_LENS tag.");
  ptr = &(dimLens[0]);
  rval = mbImpl->tag_set_by_ptr(dimLensTag, &_fileSet, 1, &ptr, &dimLensSz);
  ERRORR(rval, "Trouble setting data for __DIM_LENS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__VAR_NAMES>
  Tag varNamesTag = 0;
  tag_name = "__VAR_NAMES";
  std::string varnames;
  std::map<std::string, ReadNC::VarData>::iterator mapIter;
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    varnames.append(mapIter->first);
    varnames.push_back('\0');
  }
  int varnamesSz = varnames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, varNamesTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __VAR_NAMES tag.");
  ptr = varnames.c_str();
  rval = mbImpl->tag_set_by_ptr(varNamesTag, &_fileSet, 1, &ptr, &varnamesSz);
  ERRORR(rval, "Trouble setting data for __VAR_NAMES tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // __<dim_name>_LOC_MINMAX (for time)
  for (unsigned int i = 0; i != dimNamesSz; i++) {
    if (dimNames[i] == "time" || dimNames[i] == "Time" || dimNames[i] == "t") {
      std::stringstream ss_tag_name;
      ss_tag_name << "__" << dimNames[i] << "_LOC_MINMAX";
      tag_name = ss_tag_name.str();
      Tag tagh = 0;
      std::vector<int> val(2, 0);
      val[0] = 0;
      val[1] = nTimeSteps - 1;
      rval = mbImpl->tag_get_handle(tag_name.c_str(), 2, MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT);
      ERRORR(rval, "Trouble creating __<dim_name>_LOC_MINMAX tag.");
      rval = mbImpl->tag_set_data(tagh, &_fileSet, 1, &val[0]);
      ERRORR(rval, "Trouble setting data for __<dim_name>_LOC_MINMAX tag.");
      if (MB_SUCCESS == rval)
        dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
    }
  }

  // __<dim_name>_LOC_VALS (for time)
  for (unsigned int i = 0; i != dimNamesSz; i++) {
    if (dimNames[i] == "time" || dimNames[i] == "Time" || dimNames[i] == "t") {
      std::vector<int> val;
      if (!tstep_nums.empty())
        val = tstep_nums;
      else {
        val.resize(tVals.size());
        for (unsigned int j = 0; j != tVals.size(); j++)
          val[j] = j;
      }
      Tag tagh = 0;
      std::stringstream ss_tag_name;
      ss_tag_name << "__" << dimNames[i] << "_LOC_VALS";
      tag_name = ss_tag_name.str();
      rval = mbImpl->tag_get_handle(tag_name.c_str(), val.size(), MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT);
      ERRORR(rval, "Trouble creating __<dim_name>_LOC_VALS tag.");
      rval = mbImpl->tag_set_data(tagh, &_fileSet, 1, &val[0]);
      ERRORR(rval, "Trouble setting data for __<dim_name>_LOC_VALS tag.");
      if (MB_SUCCESS == rval)
        dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
    }
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
    for (unsigned int i = 0; i != varDimSz; i++) {
      Tag tmptag = 0;
      std::string tmptagname = dimNames[varInfo[mapIter->first].varDims[i]];
      mbImpl->tag_get_handle(tmptagname.c_str(), 0, MB_TYPE_OPAQUE, tmptag, MB_TAG_ANY);
      varInfo[mapIter->first].varTags[i] = tmptag;
    }
    rval = mbImpl->tag_get_handle(tag_name.c_str(), varDimSz, MB_TYPE_HANDLE, varNamesDimsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<var_name>_DIMS tag.");
    rval = mbImpl->tag_set_data(varNamesDimsTag, &_fileSet, 1, &(varInfo[mapIter->first].varTags[0]));
    ERRORR(rval, "Trouble setting data for __<var_name>_DIMS tag.");
    if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());
  }

  // <PARTITION_METHOD>
  Tag part_tag = scdi->part_method_tag();
  if (!part_tag)
    ERRORR(MB_FAILURE, "Trouble getting partition method tag.");
  rval = mbImpl->tag_set_data(part_tag, &_fileSet, 1, &partMethod);
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
  rval = mbImpl->tag_set_by_ptr(globalAttTag, &_fileSet, 1, &gattptr, &globalAttSz);
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
  rval = mbImpl->tag_set_data(globalAttLenTag, &_fileSet, 1, &gattLen[0]);
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
    rval = mbImpl->tag_set_by_ptr(varAttTag, &_fileSet, 1, &varAttPtr, &varAttSz);
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
    rval = mbImpl->tag_set_data(varAttLenTag, &_fileSet, 1, &varAttLen[0]);
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
  rval = mbImpl->tag_set_data(varNamesLocsTag, &_fileSet, 1, &varNamesLocs[0]);
  ERRORR(rval, "Trouble setting data for __VAR_NAMES_LOCATIONS tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  // <__MESH_TYPE>
  Tag meshTypeTag = 0;
  tag_name = "__MESH_TYPE";
  std::string meshTypeName = get_mesh_type_name();

  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, meshTypeTag, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __MESH_TYPE tag.");
  ptr = meshTypeName.c_str();
  int leng = meshTypeName.size();
  rval = mbImpl->tag_set_by_ptr(meshTypeTag, &_fileSet, 1, &ptr, &leng);
  ERRORR(rval, "Trouble setting data for __MESH_TYPE tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());

  return MB_SUCCESS;
}

ErrorCode NCHelper::read_variable_setup(std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
                                        std::vector<ReadNC::VarData>& vdatas, std::vector<ReadNC::VarData>& vsetdatas)
{
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  std::map<std::string, ReadNC::VarData>::iterator mit;

  // If empty read them all (except ignored variables)
  if (var_names.empty()) {
    for (mit = varInfo.begin(); mit != varInfo.end(); ++mit) {
      ReadNC::VarData vd = (*mit).second;

      // If read all variables at once, skip ignored ones
      // Upon creation of dummy variables, tag values have already been set
      if (ignoredVarNames.find(vd.varName) != ignoredVarNames.end() ||
          dummyVarNames.find(vd.varName) != dummyVarNames.end())
         continue;

      if (vd.entLoc == ReadNC::ENTLOCSET)
        vsetdatas.push_back(vd);
      else
        vdatas.push_back(vd);
    }
  }
  else {
    // Read specified variables (might include ignored ones)
    for (unsigned int i = 0; i < var_names.size(); i++) {
      mit = varInfo.find(var_names[i]);
      if (mit != varInfo.end()) {
        ReadNC::VarData vd = (*mit).second;

        // Upon creation of dummy variables, tag values have already been set
        if (dummyVarNames.find(vd.varName) != dummyVarNames.end())
           continue;

        if (vd.entLoc == ReadNC::ENTLOCSET)
          vsetdatas.push_back(vd);
        else
          vdatas.push_back(vd);
      }
      else {
        ERRORR(MB_FAILURE, "Couldn't find variable.");
      }
    }
  }

  if (tstep_nums.empty() && nTimeSteps > 0) {
    // No timesteps input, get them all
    for (int i = 0; i < nTimeSteps; i++)
      tstep_nums.push_back(i);
  }

  if (!tstep_nums.empty()) {
    for (unsigned int i = 0; i < vdatas.size(); i++) {
      vdatas[i].varTags.resize(tstep_nums.size(), 0);
      vdatas[i].varDatas.resize(tstep_nums.size());
      vdatas[i].readStarts.resize(tstep_nums.size());
      vdatas[i].readCounts.resize(tstep_nums.size());
    }

    for (unsigned int i = 0; i < vsetdatas.size(); i++) {
      if ((std::find(vsetdatas[i].varDims.begin(), vsetdatas[i].varDims.end(), tDim) != vsetdatas[i].varDims.end())
          && (vsetdatas[i].varDims.size() != 1)) {
        vsetdatas[i].varTags.resize(tstep_nums.size(), 0);
        vsetdatas[i].varDatas.resize(tstep_nums.size());
        vsetdatas[i].readStarts.resize(tstep_nums.size());
        vsetdatas[i].readCounts.resize(tstep_nums.size());
      }
      else {
        vsetdatas[i].varTags.resize(1, 0);
        vsetdatas[i].varDatas.resize(1);
        vsetdatas[i].readStarts.resize(1);
        vsetdatas[i].readCounts.resize(1);
      }
    }
  }

  return MB_SUCCESS;
}

ErrorCode NCHelper::read_variable_to_set(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  Interface*& mbImpl = _readNC->mbImpl;
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_variable_to_set_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables to set.");

  // Finally, read into that space
  int success;
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      void* data = vdatas[i].varDatas[t];

      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR:
          success = NCFUNCAG(_vara_text)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              (char*) data);
          ERRORS(success, "Failed to read char data.");
          break;
        case NC_DOUBLE:
          success = NCFUNCAG(_vara_double)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              (double*) data);
          ERRORS(success, "Failed to read double data.");
          break;
        case NC_FLOAT:
          success = NCFUNCAG(_vara_float)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              (float*) data);
          ERRORS(success, "Failed to read float data.");
          break;
        case NC_INT:
          success = NCFUNCAG(_vara_int)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              (int*) data);
          ERRORS(success, "Failed to read int data.");
          break;
        case NC_SHORT:
          success = NCFUNCAG(_vara_short)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              (short*) data);
          ERRORS(success, "Failed to read short data.");
          break;
        default:
          success = 1;
      }

      if (success)
        ERRORR(MB_FAILURE, "Trouble reading variable.");

      dbgOut.tprintf(2, "Converting variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      rval = convert_variable(vdatas[i], t);
      ERRORR(rval, "Failed to convert variable.");

      dbgOut.tprintf(2, "Setting data for variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      rval = mbImpl->tag_set_by_ptr(vdatas[i].varTags[t], &_fileSet, 1, &data, &vdatas[i].sz);
      ERRORR(rval, "Failed to set data for variable.");

      // Memory pointed by pointer data can be deleted, as tag_set_by_ptr() has already copied the tag values
      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR:
          delete[] (char*) data;
          break;
        case NC_DOUBLE:
        case NC_FLOAT:
          delete[] (double*) data;
          break;
        case NC_INT:
        case NC_SHORT:
          delete[] (int*) data;
          break;
        default:
          break;
      }
      vdatas[i].varDatas[t] = NULL;

      if (vdatas[i].varDims.size() <= 1 || !vdatas[i].has_t)
        break;
    }
  }

  // Debug output, if requested
  if (1 == dbgOut.get_verbosity()) {
    dbgOut.printf(1, "Read variables: %s", vdatas.begin()->varName.c_str());
    for (unsigned int i = 1; i < vdatas.size(); i++)
      dbgOut.printf(1, ", %s ", vdatas[i].varName.c_str());
    dbgOut.tprintf(1, "\n");
  }

  return rval;
}

ErrorCode NCHelper::convert_variable(ReadNC::VarData& var_data, int tstep_num)
{
  DebugOutput& dbgOut = _readNC->dbgOut;

  // Get ptr to tag space
  void* data = var_data.varDatas[tstep_num];

  // Get variable size
  std::size_t sz = var_data.sz;
  assert(sz > 0);

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

ErrorCode NCHelper::read_coordinate(const char* var_name, int lmin, int lmax, std::vector<double>& cvals)
{
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  std::map<std::string, ReadNC::VarData>::iterator vmit = varInfo.find(var_name);
  if (varInfo.end() == vmit)
    return MB_FAILURE;

  assert(lmin >= 0 && lmax >= lmin);
  NCDF_SIZE tstart = lmin;
  NCDF_SIZE tcount = lmax - lmin + 1;
  NCDF_DIFF dum_stride = 1;
  int fail;

  // Check size
  if ((std::size_t)tcount != cvals.size())
    cvals.resize(tcount);

  // Check to make sure it's a float or double
  if (NC_DOUBLE == (*vmit).second.varDataType) {
    fail = NCFUNCAG(_vars_double)(_fileId, (*vmit).second.varId, &tstart, &tcount, &dum_stride, &cvals[0]);
    if (fail)
      ERRORS(MB_FAILURE, "Failed to get coordinate values.");
  }
  else if (NC_FLOAT == (*vmit).second.varDataType) {
    std::vector<float> tcvals(tcount);
    fail = NCFUNCAG(_vars_float)(_fileId, (*vmit).second.varId, &tstart, &tcount, &dum_stride, &tcvals[0]);
    if (fail)
      ERRORS(MB_FAILURE, "Failed to get coordinate values.");
    std::copy(tcvals.begin(), tcvals.end(), cvals.begin());
  }
  else {
    ERRORR(MB_FAILURE, "Wrong data type for coordinate variable.");
  }

  return MB_SUCCESS;
}

ErrorCode NCHelper::get_tag_to_set(ReadNC::VarData& var_data, int tstep_num, Tag& tagh)
{
  Interface*& mbImpl = _readNC->mbImpl;
  DebugOutput& dbgOut = _readNC->dbgOut;

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

ErrorCode NCHelper::get_tag_to_nonset(ReadNC::VarData& var_data, int tstep_num, Tag& tagh, int num_lev)
{
  Interface*& mbImpl = _readNC->mbImpl;
  DebugOutput& dbgOut = _readNC->dbgOut;

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

ErrorCode NCHelper::create_attrib_string(const std::map<std::string, ReadNC::AttData>& attMap, std::string& attVal, std::vector<int>& attLen)
{
  int success;
  std::stringstream ssAtt;
  unsigned int sz = 0;
  std::map<std::string, ReadNC::AttData>::const_iterator attIt = attMap.begin();
  for (; attIt != attMap.end(); ++attIt) {
    ssAtt << attIt->second.attName;
    ssAtt << '\0';
    void* attData = NULL;
    switch (attIt->second.attDataType) {
      case NC_BYTE:
      case NC_CHAR:
        sz = attIt->second.attLen;
        attData = (char *) malloc(sz);
        success = NCFUNC(get_att_text)(_fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (char*) attData);
        ERRORS(success, "Failed to read attribute char data.");
        ssAtt << "char;";
        break;
      case NC_DOUBLE:
        sz = attIt->second.attLen * sizeof(double);
        attData = (double *) malloc(sz);
        success = NCFUNC(get_att_double)(_fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (double*) attData);
        ERRORS(success, "Failed to read attribute double data.");
        ssAtt << "double;";
        break;
      case NC_FLOAT:
        sz = attIt->second.attLen * sizeof(float);
        attData = (float *) malloc(sz);
        success = NCFUNC(get_att_float)(_fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (float*) attData);
        ERRORS(success, "Failed to read attribute float data.");
        ssAtt << "float;";
        break;
      case NC_INT:
        sz = attIt->second.attLen * sizeof(int);
        attData = (int *) malloc(sz);
        success = NCFUNC(get_att_int)(_fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (int*) attData);
        ERRORS(success, "Failed to read attribute int data.");
        ssAtt << "int;";
        break;
      case NC_SHORT:
        sz = attIt->second.attLen * sizeof(short);
        attData = (short *) malloc(sz);
        success = NCFUNC(get_att_short)(_fileId, attIt->second.attVarId, attIt->second.attName.c_str(), (short*) attData);
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

ErrorCode NCHelper::create_dummy_variables()
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::vector<std::string>& dimNames = _readNC->dimNames;
  std::vector<int>& dimLens = _readNC->dimLens;
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  DebugOutput& dbgOut = _readNC->dbgOut;

  // Hack: look at all dimensions, and see if we have one that does not appear in the list of varInfo names
  // Right now, candidates are from unstructured meshes, such as ncol (HOMME) and nCells (MPAS)
  // For each of them, create a dummy variable with a sparse tag to store the dimension length
  for (unsigned int i = 0; i < dimNames.size(); i++) {
    // If there is a variable with this dimension name, skip
    if (varInfo.find(dimNames[i]) != varInfo.end())
      continue;

    // Create a dummy variable
    int sizeTotalVar = varInfo.size();
    std::string var_name(dimNames[i]);
    ReadNC::VarData& data = varInfo[var_name];
    data.varName = std::string(var_name);
    data.varId = sizeTotalVar;
    data.varTags.resize(1, 0);
    data.varDataType = NC_INT;
    data.varDims.resize(1);
    data.varDims[0] = (int)i;
    data.numAtts = 0;
    data.entLoc = ReadNC::ENTLOCSET;
    dummyVarNames.insert(dimNames[i]);
    dbgOut.tprintf(2, "Dummy variable created for dimension %s\n", dimNames[i].c_str());

    // Create a corresponding sparse tag
    Tag tagh;
    ErrorCode rval = mbImpl->tag_get_handle(dimNames[i].c_str(), 0, MB_TYPE_INTEGER, tagh,
                                            MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN);
    ERRORR(rval, "Failed to create tag for a dummy dimension variable.");

    // Tag value is the dimension length
    const void* ptr = &dimLens[i];
    // Tag size is 1
    int size = 1;
    rval = mbImpl->tag_set_by_ptr(tagh, &_fileSet, 1, &ptr, &size);
    ERRORR(rval, "Failed to set data for dimension tag.");

    dbgOut.tprintf(2, "Sparse tag created for dimension %s\n", dimNames[i].c_str());
  }

  return MB_SUCCESS;
}

ErrorCode NCHelper::read_variable_to_set_allocate(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  std::vector<int>& dimLens = _readNC->dimLens;
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = MB_SUCCESS;

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    if ((std::find(vdatas[i].varDims.begin(), vdatas[i].varDims.end(), tDim) != vdatas[i].varDims.end()))
      vdatas[i].has_t = true;

    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);

      // Get the tag to read into
      if (!vdatas[i].varTags[t]) {
        rval = get_tag_to_set(vdatas[i], tstep_nums[t], vdatas[i].varTags[t]);
        ERRORR(rval, "Trouble getting tag.");
      }

      // Assume point-based values for now?
      if (-1 == tDim || dimLens[tDim] <= (int) t)
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Wrong value for timestep number.");

      // Set up the dimensions and counts
      // First variable dimension is time, if it exists
      if (vdatas[i].has_t) {
        if (vdatas[i].varDims.size() != 1) {
          vdatas[i].readStarts[t].push_back(tstep_nums[t]);
          vdatas[i].readCounts[t].push_back(1);
        }
        else {
          vdatas[i].readStarts[t].push_back(0);
          vdatas[i].readCounts[t].push_back(tstep_nums.size());
        }
      }

      // Set up other dimensions and counts
      if (vdatas[i].varDims.empty()) {
        // Scalar variable
        vdatas[i].readStarts[t].push_back(0);
        vdatas[i].readCounts[t].push_back(1);
      }
      else {
        for (unsigned int idx = 0; idx != vdatas[i].varDims.size(); idx++){
          if (tDim != vdatas[i].varDims[idx]){
            // Push other variable dimensions, except time, which was already pushed
            vdatas[i].readStarts[t].push_back(0);
            vdatas[i].readCounts[t].push_back(dimLens[vdatas[i].varDims[idx]]);
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
          std::cerr << "Unrecognized data type for set variable tag values" << std::endl;
          rval = MB_FAILURE;
      }

      if (vdatas[i].varDims.size() <= 1 || !vdatas[i].has_t)
        break;
    }
  }

  return rval;
}

ErrorCode ScdNCHelper::check_existing_mesh() {
  Interface*& mbImpl = _readNC->mbImpl;

  // Get the number of vertices
  int num_verts;
  ErrorCode rval = mbImpl->get_number_entities_by_dimension(_fileSet, 0, num_verts);
  ERRORR(rval, "Trouble getting number of vertices.");

  /*
  // Check against parameters
  // When ghosting is used, this check might fail (to be updated later)
  if (num_verts > 0) {
    int expected_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2] + 1);
    if (num_verts != expected_verts) {
      ERRORR(MB_FAILURE, "Number of vertices doesn't match.");
    }
  }
  */

  // Check the number of elements too
  int num_elems;
  rval = mbImpl->get_number_entities_by_dimension(_fileSet, (-1 == lCDims[2] ? 2 : 3), num_elems);
  ERRORR(rval, "Trouble getting number of elements.");

  /*
  // Check against parameters
  // When ghosting is used, this check might fail (to be updated later)
  if (num_elems > 0) {
    int expected_elems = (lCDims[3] - lCDims[0] + 1) * (lCDims[4] - lCDims[1] + 1) * (-1 == lCDims[2] ? 1 : (lCDims[5] - lCDims[2] + 1));
    if (num_elems != expected_elems) {
      ERRORR(MB_FAILURE, "Number of elements doesn't match.");
    }
  }
  */

  return MB_SUCCESS;
}

ErrorCode ScdNCHelper::create_mesh(Range& faces)
{
  Interface*& mbImpl = _readNC->mbImpl;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
  DebugOutput& dbgOut = _readNC->dbgOut;
  ScdInterface* scdi = _readNC->scdi;
  ScdParData& parData = _readNC->parData;

  Range tmp_range;
  ScdBox* scd_box;

  ErrorCode rval = scdi->construct_box(HomCoord(lDims[0], lDims[1], lDims[2], 1), HomCoord(lDims[3], lDims[4], lDims[5], 1), 
                                       NULL, 0, scd_box, locallyPeriodic, &parData, true);
  ERRORR(rval, "Trouble creating scd vertex sequence.");

  // Add verts to tmp_range first, so we can duplicate global ids in vertex ids
  tmp_range.insert(scd_box->start_vertex(), scd_box->start_vertex() + scd_box->num_vertices() - 1);

  if (mpFileIdTag) {
    int count;
    void* data;
    rval = mbImpl->tag_iterate(*mpFileIdTag, tmp_range.begin(), tmp_range.end(), count, data);
    ERRORR(rval, "Failed to get tag iterator on file id tag.");
    assert(count == scd_box->num_vertices());
    int* fid_data = (int*) data;
    rval = mbImpl->tag_iterate(mGlobalIdTag, tmp_range.begin(), tmp_range.end(), count, data);
    ERRORR(rval, "Failed to get tag iterator on global id tag.");
    assert(count == scd_box->num_vertices());
    int* gid_data = (int*) data;
    for (int i = 0; i < count; i++)
      fid_data[i] = gid_data[i];
  }

  // Then add box set and elements to the range, then to the file set
  tmp_range.insert(scd_box->start_element(), scd_box->start_element() + scd_box->num_elements() - 1);
  tmp_range.insert(scd_box->box_set());
  rval = mbImpl->add_entities(_fileSet, tmp_range);
  ERRORR(rval, "Couldn't add new vertices to file set.");

  dbgOut.tprintf(1, "scdbox %d quads, %d vertices\n", scd_box->num_elements(), scd_box->num_vertices());

  // Set the vertex coordinates
  double *xc, *yc, *zc;
  rval = scd_box->get_coordinate_arrays(xc, yc, zc);
  ERRORR(rval, "Couldn't get vertex coordinate arrays.");

  int i, j, k, il, jl, kl;
  int dil = lDims[3] - lDims[0] + 1;
  int djl = lDims[4] - lDims[1] + 1;
  assert(dil == (int)ilVals.size() && djl == (int)jlVals.size() &&
      (-1 == lDims[2] || lDims[5] - lDims[2] + 1 == (int)levVals.size()));
//#define INDEX(i, j, k) ()
  for (kl = lDims[2]; kl <= lDims[5]; kl++) {
    k = kl - lDims[2];
    for (jl = lDims[1]; jl <= lDims[4]; jl++) {
      j = jl - lDims[1];
      for (il = lDims[0]; il <= lDims[3]; il++) {
        i = il - lDims[0];
        unsigned int pos = i + j * dil + k * dil * djl;
        xc[pos] = ilVals[i];
        yc[pos] = jlVals[j];
        zc[pos] = (-1 == lDims[2] ? 0.0 : levVals[k]);
      }
    }
  }
//#undef INDEX

#ifndef NDEBUG
  int num_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2] + 1);
  std::vector<int> gids(num_verts);
  Range verts(scd_box->start_vertex(), scd_box->start_vertex() + scd_box->num_vertices() - 1);
  rval = mbImpl->tag_get_data(mGlobalIdTag, verts, &gids[0]);
  ERRORR(rval, "Trouble getting gid values.");
  int vmin = *(std::min_element(gids.begin(), gids.end())), vmax = *(std::max_element(gids.begin(), gids.end()));
  dbgOut.tprintf(1, "Vertex gids %d-%d\n", vmin, vmax);
#endif

  // Add elements to the range passed in
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

ErrorCode ScdNCHelper::read_variables(std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
{
  std::vector<ReadNC::VarData> vdatas;
  std::vector<ReadNC::VarData> vsetdatas;

  ErrorCode rval = read_variable_setup(var_names, tstep_nums, vdatas, vsetdatas);
  ERRORR(rval, "Trouble setting up read variable.");

  // Create COORDS tag for quads
  rval = create_quad_coordinate_tag();
  ERRORR(rval, "Trouble creating coordinate tags to entities quads");

  if (!vsetdatas.empty()) {
    rval = read_variable_to_set(vsetdatas, tstep_nums);
    ERRORR(rval, "Trouble read variables to set.");
  }

  if (!vdatas.empty()) {
    rval = read_scd_variable_to_nonset(vdatas, tstep_nums);
    ERRORR(rval, "Trouble read variables to entities verts/edges/faces.");
  }

  return MB_SUCCESS;
}

ErrorCode ScdNCHelper::read_scd_variable_to_nonset_allocate(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::vector<int>& dimLens = _readNC->dimLens;
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = MB_SUCCESS;

  Range* range = NULL;

  // Get vertices in set
  Range verts;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 0, verts);
  ERRORR(rval, "Trouble getting vertices in set.");
  assert("Should only have a single vertex subrange, since they were read in one shot" &&
      verts.psize() == 1);

  Range edges;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 1, edges);
  ERRORR(rval, "Trouble getting edges in set.");

  // Get faces in set
  Range faces;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 2, faces);
  ERRORR(rval, "Trouble getting faces in set.");
  assert("Should only have a single face subrange, since they were read in one shot" &&
      faces.psize() == 1);

#ifdef USE_MPI
  moab::Range faces_owned;
  bool& isParallel = _readNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _readNC->myPcomm;
    rval = myPcomm->filter_pstatus(faces, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &faces_owned);
    ERRORR(rval, "Trouble getting owned faces in set.");
  }
  else
    faces_owned = faces; // Not running in parallel, but still with MPI
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);

      // Get the tag to read into
      if (!vdatas[i].varTags[t]) {
        rval = get_tag_to_nonset(vdatas[i], tstep_nums[t], vdatas[i].varTags[t], vdatas[i].numLev);
        ERRORR(rval, "Trouble getting tag.");
      }

      // Assume point-based values for now?
      if (-1 == tDim || dimLens[tDim] <= (int) t) {
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Wrong value for timestep number.");
      }
      else if (vdatas[i].varDims[0] != tDim) {
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Non-default timestep number given for time-independent variable.");
      }

      // Set up the dimensions and counts
      // First: time
      vdatas[i].readStarts[t].push_back(tstep_nums[t]);
      vdatas[i].readCounts[t].push_back(1);

      // Next: numLev, even if it is 1
      vdatas[i].readStarts[t].push_back(0);
      vdatas[i].readCounts[t].push_back(vdatas[i].numLev);

      // Finally: y and x
      switch (vdatas[i].entLoc) {
        case ReadNC::ENTLOCVERT:
          // Vertices
          vdatas[i].readStarts[t].push_back(lDims[1]);
          vdatas[i].readCounts[t].push_back(lDims[4] - lDims[1] + 1);
          vdatas[i].readStarts[t].push_back(lDims[0]);
          vdatas[i].readCounts[t].push_back(lDims[3] - lDims[0] + 1);
          assert(vdatas[i].readStarts[t].size() == vdatas[i].varDims.size());
          range = &verts;
          break;
        case ReadNC::ENTLOCNSEDGE:
          ERRORR(MB_FAILURE, "Reading edge data not implemented yet.");
          break;
        case ReadNC::ENTLOCEWEDGE:
          ERRORR(MB_FAILURE, "Reading edge data not implemented yet.");
          break;
        case ReadNC::ENTLOCFACE:
          // Faces
          vdatas[i].readStarts[t].push_back(lCDims[1]);
          vdatas[i].readCounts[t].push_back(lCDims[4] - lCDims[1] + 1);
          vdatas[i].readStarts[t].push_back(lCDims[0]);
          vdatas[i].readCounts[t].push_back(lCDims[3] - lCDims[0] + 1);
          assert(vdatas[i].readStarts[t].size() == vdatas[i].varDims.size());
#ifdef USE_MPI
          range = &faces_owned;
#else
          range = &faces;
#endif
          break;
        case ReadNC::ENTLOCSET:
          // Set
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

    // Calculate variable size
    std::size_t sz = 1;
    for (std::size_t idx = 0; idx != vdatas[i].readCounts[0].size(); idx++)
      sz *= vdatas[i].readCounts[0][idx];
    vdatas[i].sz = sz;
  }

  return rval;
}

ErrorCode ScdNCHelper::read_scd_variable_to_nonset(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_scd_variable_to_nonset_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables.");

  // Finally, read into that space
  int success;
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    std::size_t sz = vdatas[i].sz;

    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      void* data = vdatas[i].varDatas[t];
      size_t ni = vdatas[i].readCounts[t][2];
      size_t nj = vdatas[i].readCounts[t][3];
      size_t nk = vdatas[i].readCounts[t][1];

      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR: {
          std::vector<char> tmpchardata(sz);
          success = NCFUNCAG(_vara_text)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              &tmpchardata[0]);
          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpchardata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpchardata.size(); idx++)
              ((char*) data)[idx] = tmpchardata[idx];
          }
          ERRORS(success, "Failed to read char data.");
          break;
        }
        case NC_DOUBLE: {
          std::vector<double> tmpdoubledata(sz);
          success = NCFUNCAG(_vara_double)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              &tmpdoubledata[0]);
          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpdoubledata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpdoubledata.size(); idx++)
              ((double*) data)[idx] = tmpdoubledata[idx];
          }
          ERRORS(success, "Failed to read double data.");
          break;
        }
        case NC_FLOAT: {
          std::vector<float> tmpfloatdata(sz);
          success = NCFUNCAG(_vara_float)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              &tmpfloatdata[0]);
          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpfloatdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpfloatdata.size(); idx++)
              ((float*) data)[idx] = tmpfloatdata[idx];
          }
          ERRORS(success, "Failed to read float data.");
          break;
        }
        case NC_INT: {
          std::vector<int> tmpintdata(sz);
          success = NCFUNCAG(_vara_int)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              &tmpintdata[0]);
          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpintdata[0]);
          else {
            for (std::size_t idx = 0; idx != tmpintdata.size(); idx++)
              ((int*) data)[idx] = tmpintdata[idx];
          }
          ERRORS(success, "Failed to read int data.");
          break;
        }
        case NC_SHORT: {
          std::vector<short> tmpshortdata(sz);
          success = NCFUNCAG(_vara_short)(_fileId, vdatas[i].varId, &vdatas[i].readStarts[t][0], &vdatas[i].readCounts[t][0],
              &tmpshortdata[0]);
          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik(ni, nj, nk, data, &tmpshortdata[0]);
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

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Converting variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = convert_variable(vdatas[i], t);
      if (MB_SUCCESS != tmp_rval)
        rval = tmp_rval;
    }
  }

  // Debug output, if requested
  if (1 == dbgOut.get_verbosity()) {
    dbgOut.printf(1, "Read variables: %s", vdatas.begin()->varName.c_str());
    for (unsigned int i = 1; i < vdatas.size(); i++)
      dbgOut.printf(1, ", %s ", vdatas[i].varName.c_str());
    dbgOut.tprintf(1, "\n");
  }

  return rval;
}

ErrorCode ScdNCHelper::create_quad_coordinate_tag() {
  Interface*& mbImpl = _readNC->mbImpl;

  Range ents;
  ErrorCode rval = mbImpl->get_entities_by_type(_fileSet, moab::MBQUAD, ents);
  ERRORR(rval, "Trouble getting QUAD entity.");

  std::size_t numOwnedEnts = 0;
#ifdef USE_MPI
  Range ents_owned;
  bool& isParallel = _readNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _readNC->myPcomm;
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

ErrorCode UcdNCHelper::read_variables(std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
{
  std::vector<ReadNC::VarData> vdatas;
  std::vector<ReadNC::VarData> vsetdatas;

  ErrorCode rval = read_variable_setup(var_names, tstep_nums, vdatas, vsetdatas);
  ERRORR(rval, "Trouble setting up read variable.");

  if (!vsetdatas.empty()) {
    rval = read_variable_to_set(vsetdatas, tstep_nums);
    ERRORR(rval, "Trouble read variables to set.");
  }

  if (!vdatas.empty()) {
#ifdef PNETCDF_FILE
    // With pnetcdf support, we will use async read
    rval = read_ucd_variable_to_nonset_async(vdatas, tstep_nums);
#else
    // Without pnetcdf support, we will use old read
    rval = read_ucd_variable_to_nonset(vdatas, tstep_nums);
#endif
    ERRORR(rval, "Trouble read variables to entities verts/edges/faces.");
  }

  return MB_SUCCESS;
}

} // namespace moab
