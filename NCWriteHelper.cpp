/*
 * NCWriteHelper.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: iulian
 */

#include "NCWriteHelper.hpp"
#include "NCWriteEuler.hpp"
#include "NCWriteFV.hpp"
#include "NCWriteHOMME.hpp"
#include "NCWriteMPAS.hpp"

#include "moab/WriteUtilIface.hpp"

#include <sstream>

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { _writeNC->mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { _writeNC->mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

//! Get appropriate helper instance for WriteNC class; based on some info in the file set
NCWriteHelper* NCWriteHelper::get_nc_helper(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
{
  std::string& grid_type = writeNC->grid_type;
  if (grid_type == "CAM_EUL")
    return new (std::nothrow) NCWriteEuler(writeNC, fileId, opts, fileSet);
  else if (grid_type == "CAM_FV")
    return new (std::nothrow) NCWriteFV(writeNC, fileId, opts, fileSet);
  else if (grid_type == "CAM_SE")
    return new (std::nothrow) NCWriteHOMME(writeNC, fileId, opts, fileSet);
  else if (grid_type == "MPAS")
    return new (std::nothrow) NCWriteMPAS(writeNC, fileId, opts, fileSet);

  // Unknown NetCDF grid
  return NULL;
}

ErrorCode NCWriteHelper::collect_variable_data(std::vector<std::string>& var_names)
{
  Interface*& mbImpl = _writeNC->mbImpl;
  std::vector<std::string>& dimNames = _writeNC->dimNames;
  std::vector<int>& dimLens = _writeNC->dimLens;
  std::set<std::string>& usedCoordinates = _writeNC->usedCoordinates;
  std::set<std::string>& dummyVarNames = _writeNC->dummyVarNames;
  std::map<std::string, WriteNC::VarData>& varInfo = _writeNC->varInfo;
  DebugOutput& dbgOut = _writeNC->dbgOut;

  ErrorCode rval;

  usedCoordinates.clear();

  for (size_t i = 0; i < var_names.size(); i++) {
    std::string varname = var_names[i];
    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(varname);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find one variable.");

    WriteNC::VarData& currentVarData = vit->second;

    currentVarData.has_tsteps = false;
    if ((std::find(currentVarData.varDims.begin(), currentVarData.varDims.end(), tDim) != currentVarData.varDims.end())
        && (currentVarData.varDims.size() > 1)) // So it is not time itself
      currentVarData.has_tsteps = true;

    currentVarData.numLev = 1;
    if ((std::find(currentVarData.varDims.begin(), currentVarData.varDims.end(), tDim) != currentVarData.varDims.end()))
      currentVarData.numLev = nLevels;

    dbgOut.tprintf(2, "    for variable %s varDims.size %d \n", varname.c_str(), (int)currentVarData.varDims.size());
    for (size_t j = 0; j < currentVarData.varDims.size(); j++) {
      std::string dimName = dimNames[currentVarData.varDims[j]];
      vit = varInfo.find(dimName);
      if (vit == varInfo.end())
        ERRORR(MB_FAILURE, "Can't find one coordinate variable.");

      usedCoordinates.insert(dimName); // Collect those used, we will need to write them to the file
      dbgOut.tprintf(2, "    for variable %s need dimension %s with length %d\n", varname.c_str(), dimName.c_str(), dimLens[currentVarData.varDims[j]]);
    }

    if (currentVarData.has_tsteps) {
      int index = 0;
      // FIXME: Should use tstep_nums (from writing options) later
      while (true) {
        Tag indexedTag = 0;
        std::stringstream ssTagNameWithIndex;
        ssTagNameWithIndex << varname << index;
        rval = mbImpl->tag_get_handle(ssTagNameWithIndex.str().c_str(), indexedTag);
        if (MB_SUCCESS != rval)
          break;
        dbgOut.tprintf(2, "    found indexed tag %d with name %s\n", index, ssTagNameWithIndex.str().c_str());
        currentVarData.varTags.push_back(indexedTag);
        index++; // We should get out of the loop at some point

        // The type of the tag is fixed though
        DataType type;
        rval = mbImpl->tag_get_data_type(indexedTag, type);
        ERRORR(rval, "Can't get tag type.");

        currentVarData.varDataType = NC_DOUBLE;
        if (MB_TYPE_INTEGER == type)
          currentVarData.varDataType = NC_INT;
      }
    }
    else {
      // Get the tag with varname
      Tag tag = 0;
      rval = mbImpl->tag_get_handle(varname.c_str(), tag);
      ERRORR(rval, "Can't find one tag.");
      currentVarData.varTags.push_back(tag); // Really, only one for these
      const void* data;
      int size;
      rval = mbImpl->tag_get_by_ptr(tag, &_fileSet, 1, &data, &size);
      ERRORR(rval, "Can't get tag values.");

      // Find the type of tag, and use it
      DataType type;
      rval = mbImpl->tag_get_data_type(tag, type);
      ERRORR(rval, "Can't get tag type.");

      currentVarData.varDataType = NC_DOUBLE;
      if (MB_TYPE_INTEGER == type)
        currentVarData.varDataType = NC_INT;

      assert(currentVarData.memoryHogs.size() == 0); // Nothing so far
      currentVarData.memoryHogs.push_back((void*)data);

      if (currentVarData.varDims.empty()) {
        // Scalar variable
        currentVarData.writeStarts.push_back(0);
        currentVarData.writeCounts.push_back(1);
      }
      else {
        for (unsigned int idx = 0; idx != currentVarData.varDims.size(); idx++){
          currentVarData.writeStarts.push_back(0);
          currentVarData.writeCounts.push_back(dimLens[currentVarData.varDims[idx]]);
        }
      }
    }
  }

  // Check that for used coordinates we have found the tags
  for (std::set<std::string>::iterator setIt = usedCoordinates.begin();
      setIt != usedCoordinates.end(); ++setIt) {
    const std::string& coordName = *setIt;

    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(coordName);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find one coordinate variable.");

    WriteNC::VarData& varCoordData = vit->second;
    Tag coordTag = 0;
    rval = mbImpl->tag_get_handle(coordName.c_str(), coordTag);
    ERRORR(rval, "Can't find one tag.");
    varCoordData.varTags.push_back(coordTag); // Really, only one for these

    const void* data;
    int sizeCoordinate;
    rval = mbImpl->tag_get_by_ptr(coordTag, &_fileSet, 1, &data, &sizeCoordinate);
    ERRORR(rval, "Can't get coordinate values.");
    dbgOut.tprintf(2, "    found coordinate tag with name %s and length %d\n", coordName.c_str(),
        sizeCoordinate);

    // Get dimension length (the only dimension of this coordinate variable, with the same name)
    assert(1 == varCoordData.varDims.size());
    int coordDimLen = dimLens[varCoordData.varDims[0]];

    if (dummyVarNames.find(coordName) != dummyVarNames.end()) {
      // For a dummy coordinate variable, the tag size is always 1
      // The number of coordinates should be set to dimension length, instead of 1
      assert(1 == sizeCoordinate);
      sizeCoordinate = coordDimLen;
    }
    else {
      // The number of coordinates should be exactly the same as dimension length
      assert(sizeCoordinate == coordDimLen);
    }

    // This is the length
    varCoordData.sz = sizeCoordinate;
    varCoordData.writeStarts.resize(1);
    varCoordData.writeStarts[0] = 0;
    varCoordData.writeCounts.resize(1);
    varCoordData.writeCounts[0] = sizeCoordinate;

    // Find the type of tag, and use it
    DataType type;
    rval = mbImpl->tag_get_data_type(coordTag, type);
    ERRORR(rval, "Can't get tag type.");

    varCoordData.varDataType = NC_DOUBLE;
    if (MB_TYPE_INTEGER == type)
      varCoordData.varDataType = NC_INT;

    assert(0 == varCoordData.memoryHogs.size()); // Nothing so far
    varCoordData.memoryHogs.push_back((void*)data);
  }

  return MB_SUCCESS;
}

ErrorCode NCWriteHelper::init_file(std::vector<std::string>& var_names)
{
  std::vector<std::string>& dimNames = _writeNC->dimNames;
  std::set<std::string>& usedCoordinates = _writeNC->usedCoordinates;
  std::set<std::string>& dummyVarNames = _writeNC->dummyVarNames;
  std::map<std::string, WriteNC::VarData>& varInfo = _writeNC->varInfo;
  std::map<std::string, WriteNC::AttData>& globalAtts = _writeNC->globalAtts;
  DebugOutput& dbgOut = _writeNC->dbgOut;

  // First initialize all coordinates, then fill VarData for actual variables (and dimensions)
  // Check that for used coordinates we have found the tags
  for (std::set<std::string>::iterator setIt = usedCoordinates.begin();
      setIt != usedCoordinates.end(); ++setIt) {
    const std::string& coordName = *setIt;

    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(coordName);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find one coordinate variable.");

    WriteNC::VarData& varCoordData = vit->second;
    varCoordData.varDims.resize(1);

    /* int nc_def_dim (int ncid, const char *name, size_t len, int *dimidp);
       * example:  status = nc_def_dim(fileId, "lat", 18L, &latid);
    */

    // Actually define a dimension
    if (NCFUNC(def_dim)(_fileId, coordName.c_str(), (size_t)varCoordData.sz,
        &varCoordData.varDims[0]) != NC_NOERR)
     ERRORR(MB_FAILURE, "Failed to generate dimension.");

    dbgOut.tprintf(2, "    for coordName %s dim id is %d \n", coordName.c_str(), (int)varCoordData.varDims[0]);

    // Create a variable with the same name, and its only dimension the one we just defined
    /*
     * int nc_def_var (int ncid, const char *name, nc_type xtype,
                       int ndims, const int dimids[], int *varidp);
       example: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c/nc_005fdef_005fvar.html#nc_005fdef_005fvar
     */

    // Skip dummy coordinate variables (e.g. ncol)
    if (dummyVarNames.find(coordName) != dummyVarNames.end())
      continue;

    // Define a coordinate variable
    if (NCFUNC(def_var)(_fileId, coordName.c_str(), varCoordData.varDataType,
        1, &(varCoordData.varDims[0]), &varCoordData.varId) != NC_NOERR)
      ERRORR(MB_FAILURE, "Failed to create coordinate variable.");

    dbgOut.tprintf(2, "    for coordName %s variable id is %d \n", coordName.c_str(), varCoordData.varId);
  }

  // Now look at requested variables, and update from the index in dimNames to the actual dimension id
  for (size_t i = 0; i < var_names.size(); i++) {
    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(var_names[i]);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find variable requested.");

    WriteNC::VarData& variableData = vit->second;
    int numDims = (int)variableData.varDims.size();
    // The index is for dimNames; we need to find out the actual dimension id (from above)
    for (int j = 0; j < numDims; j++) {
      std::string dimName = dimNames[variableData.varDims[j]];
      std::map<std::string, WriteNC::VarData>::iterator vit2 = varInfo.find(dimName);
      if (vit2 == varInfo.end())
        ERRORR(MB_FAILURE, "Can't find coordinate variable requested.");

      WriteNC::VarData& coordData = vit2->second;
      // Index in dimNames to actual dimension id
      variableData.varDims[j] = coordData.varDims[0]; // This one, being a coordinate, is the only one
      dbgOut.tprintf(2, "          dimension with index %d name %s has ID %d \n",
          j, dimName.c_str(), variableData.varDims[j]);
    }

    // Define the variable now:
    if (NCFUNC(def_var)(_fileId, var_names[i].c_str(), variableData.varDataType,
        (int)variableData.varDims.size(), &(variableData.varDims[0]),
        &variableData.varId) != NC_NOERR)
      ERRORR(MB_FAILURE, "Failed to create coordinate variable.");

    dbgOut.tprintf(2, "    for variable %s variable id is %d \n", var_names[i].c_str(), variableData.varId);
    // Now define the variable, with all dimensions
  }

  // Define global attributes (exactly copied from the original file for the time being)
  // Should we modify some of them (e.g. revision_Id) later?
  std::map<std::string, WriteNC::AttData>::iterator attIt;
  for (attIt = globalAtts.begin(); attIt != globalAtts.end(); ++attIt) {
    const std::string& attName = attIt->first;
    WriteNC::AttData& attData = attIt->second;
    NCDF_SIZE& attLen = attData.attLen;
    nc_type& attDataType = attData.attDataType;
    const std::string& attValue = attData.attValue;

    switch (attDataType) {
      case NC_BYTE:
      case NC_CHAR:
        if (NC_NOERR != NCFUNC(put_att_text)(_fileId, NC_GLOBAL, attName.c_str(), attLen, attValue.c_str()))
          ERRORR(MB_FAILURE, "Failed to define text type attribute.");
        break;
      case NC_DOUBLE:
        if (NC_NOERR != NCFUNC(put_att_double)(_fileId, NC_GLOBAL, attName.c_str(), NC_DOUBLE, 1, (double*)attValue.c_str()))
          ERRORR(MB_FAILURE, "Failed to define double type attribute.");
        break;
      case NC_FLOAT:
        if (NC_NOERR != NCFUNC(put_att_float)(_fileId, NC_GLOBAL, attName.c_str(), NC_FLOAT, 1, (float*)attValue.c_str()))
          ERRORR(MB_FAILURE, "Failed to define float type attribute.");
        break;
      case NC_INT:
        if (NC_NOERR != NCFUNC(put_att_int)(_fileId, NC_GLOBAL, attName.c_str(), NC_INT, 1, (int*)attValue.c_str()))
          ERRORR(MB_FAILURE, "Failed to define int type attribute.");
        break;
      case NC_SHORT:
        if (NC_NOERR != NCFUNC(put_att_short)(_fileId, NC_GLOBAL, attName.c_str(), NC_SHORT, 1, (short*)attValue.c_str()))
          ERRORR(MB_FAILURE, "Failed to define short type attribute.");
        break;
      default:
        ERRORR(MB_FAILURE, "Unknown attribute data type.");
    }
  }

  // Take it out of define mode
  if (NC_NOERR != NCFUNC(enddef)(_fileId))
    ERRORR(MB_FAILURE, "Failed to close define mode.");

  return MB_SUCCESS;
}

ErrorCode ScdNCWriteHelper::collect_mesh_info()
{
  Interface*& mbImpl = _writeNC->mbImpl;
  std::vector<std::string>& dimNames = _writeNC->dimNames;
  std::vector<int>& dimLens = _writeNC->dimLens;

  ErrorCode rval;

  // Look for time dimension
  std::vector<std::string>::iterator vecIt;
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end())
    tDim = vecIt - dimNames.begin();
  else if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "t")) != dimNames.end())
    tDim = vecIt - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'time' or 't' dimension.");
  }
  nTimeSteps = dimLens[tDim];

  // Get number of levels
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "lev")) != dimNames.end())
    levDim = vecIt - dimNames.begin();
  else if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "ilev")) != dimNames.end())
    levDim = vecIt - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'lev' or 'ilev' dimension.");
  }
  nLevels = dimLens[levDim];

  // __<dim_name>_LOC_MINMAX (for slon, slat, lon and lat)
  Tag convTag = 0;
  rval = mbImpl->tag_get_handle("__slon_LOC_MINMAX", 0, MB_TYPE_INTEGER, convTag, MB_TAG_ANY);
  ERRORR(rval, "Trouble getting conventional tag __slon_LOC_MINMAX.");
  int val[2];
  rval = mbImpl->tag_get_data(convTag, &_fileSet, 1, val);
  ERRORR(rval, "Trouble getting values for conventional tag __slon_LOC_MINMAX.");
  lDims[0] = val[0];
  lDims[3] = val[1];

  rval = mbImpl->tag_get_handle("__slat_LOC_MINMAX", 0, MB_TYPE_INTEGER, convTag, MB_TAG_ANY);
  ERRORR(rval, "Trouble getting conventional tag __slat_LOC_MINMAX.");
  rval = mbImpl->tag_get_data(convTag, &_fileSet, 1, val);
  ERRORR(rval, "Trouble getting values for conventional tag __slat_LOC_MINMAX.");
  lDims[1] = val[0];
  lDims[4] = val[1];

  rval = mbImpl->tag_get_handle("__lon_LOC_MINMAX", 0, MB_TYPE_INTEGER, convTag, MB_TAG_ANY);
  ERRORR(rval, "Trouble getting conventional tag __lon_LOC_MINMAX.");
  rval = mbImpl->tag_get_data(convTag, &_fileSet, 1, val);
  ERRORR(rval, "Trouble getting values for conventional tag __lon_LOC_MINMAX.");
  lCDims[0] = val[0];
  lCDims[3] = val[1];

  rval = mbImpl->tag_get_handle("__lat_LOC_MINMAX", 0, MB_TYPE_INTEGER, convTag, MB_TAG_ANY);
  ERRORR(rval, "Trouble getting conventional tag __lat_LOC_MINMAX.");
  rval = mbImpl->tag_get_data(convTag, &_fileSet, 1, val);
  ERRORR(rval, "Trouble getting values for conventional tag __lat_LOC_MINMAX.");
  lCDims[1] = val[0];
  lCDims[4] = val[1];

  return MB_SUCCESS;
}

ErrorCode ScdNCWriteHelper::collect_variable_data(std::vector<std::string>& var_names)
{
  NCWriteHelper::collect_variable_data(var_names);

  std::map<std::string, WriteNC::VarData>& varInfo = _writeNC->varInfo;

  for (size_t i = 0; i < var_names.size(); i++) {
    std::string varname = var_names[i];
    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(varname);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find one variable.");

    WriteNC::VarData& currentVarData = vit->second;
    if (currentVarData.has_tsteps) {
      // Support non-set variables with 4 dimensions like (time, lev, lat, lon)
      assert(4 == currentVarData.varDims.size());

      // Time should be the first dimension
      assert(tDim == currentVarData.varDims[0]);

      // Set up writeStarts and writeCounts
      currentVarData.writeStarts.resize(4);
      currentVarData.writeCounts.resize(4);

      // First: time
      currentVarData.writeStarts[0] = 0; // This value is timestep dependent, will be set later
      currentVarData.writeCounts[0] = 1;

      // Next: lev
      currentVarData.writeStarts[1] = 0;
      currentVarData.writeCounts[1] = currentVarData.numLev;

      // Finally: lat (or slat) and lon (or slon)
      switch (currentVarData.entLoc) {
        case WriteNC::ENTLOCFACE:
          // Faces
          currentVarData.writeStarts[2] = lCDims[1];
          currentVarData.writeCounts[2] = lCDims[4] - lCDims[1] + 1;
          currentVarData.writeStarts[3] = lCDims[0];
          currentVarData.writeCounts[3] = lCDims[3] - lCDims[0] + 1;
          break;
        default:
          ERRORR(MB_FAILURE, "Not implemented yet.");
      }
    }

    // Get variable size
    currentVarData.sz = 1;
    for (std::size_t idx = 0; idx != currentVarData.writeCounts.size(); idx++)
      currentVarData.sz *= currentVarData.writeCounts[idx];
  }

  return MB_SUCCESS;
}

ErrorCode ScdNCWriteHelper::write_values(std::vector<std::string>& var_names)
{
  Interface*& mbImpl = _writeNC->mbImpl;
  std::set<std::string>& usedCoordinates = _writeNC->usedCoordinates;
  std::set<std::string>& dummyVarNames = _writeNC->dummyVarNames;
  std::map<std::string, WriteNC::VarData>& varInfo = _writeNC->varInfo;

  ErrorCode rval;

  // Start with coordinates
  for (std::set<std::string>::iterator setIt = usedCoordinates.begin();
      setIt != usedCoordinates.end(); ++setIt) {
    const std::string& coordName = *setIt;

    // Skip dummy coordinate variables (if any)
    if (dummyVarNames.find(coordName) != dummyVarNames.end())
      continue;

    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(coordName);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find one coordinate variable.");

    WriteNC::VarData& varCoordData = vit->second;

    int success = 0;
    switch (varCoordData.varDataType) {
      case NC_DOUBLE:
        success = NCFUNCAP(_vara_double)(_fileId, varCoordData.varId, &varCoordData.writeStarts[0],
                  &varCoordData.writeCounts[0], (double*)(varCoordData.memoryHogs[0]));
        ERRORS(success, "Failed to write double data.");
        break;
      case NC_INT:
        success = NCFUNCAP(_vara_int)(_fileId, varCoordData.varId, &varCoordData.writeStarts[0],
                  &varCoordData.writeCounts[0], (int*)(varCoordData.memoryHogs[0]));
        ERRORS(success, "Failed to write int data.");
        break;
      default:
        success = 1;
        break;
    }
 }

  // Now look at requested var_names; if they have time, we will have a list, and write one at a time
  // if not, just write regularly
  for (size_t i = 0; i < var_names.size(); i++) {
    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(var_names[i]);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find variable requested.");

    WriteNC::VarData& variableData = vit->second;
    int numTimeSteps = (int)variableData.varTags.size();
    if (variableData.has_tsteps) {
      // Time should be the first dimension
      assert(tDim == variableData.varDims[0]);

      // Get entities of this variable
      Range ents;
      switch (variableData.entLoc) {
        case WriteNC::ENTLOCFACE:
          // Faces
          rval = mbImpl->get_entities_by_dimension(_fileSet, 2, ents);
          ERRORR(rval, "Can't get entities for faces.");
          break;
        default:
          ERRORR(MB_FAILURE, "Not implemented yet.");
      }

      // A typical variable has 4 dimensions as (time, lev, lat, lon)
      // At each timestep, we need to transpose tag format (lat, lon, lev) back
      // to NC format (lev, lat, lon) for writing
      size_t ni = variableData.writeCounts[3]; // lon
      size_t nj = variableData.writeCounts[2]; // lat
      size_t nk = variableData.writeCounts[1]; // lev

      variableData.writeCounts[0] = 1; // We will write one time step
      // FIXME: Should use tstep_nums (from writing options) later
      for (int j = 0; j < numTimeSteps; j++) {
        // We will write one time step, and count will be one; start will be different
        // We will write values directly from tag_iterate, but we should also transpose for level
        // so that means deep copy for transpose
        variableData.writeStarts[0] = j; // This is time, again
        int count;
        void* dataptr;
        rval = mbImpl->tag_iterate(variableData.varTags[j], ents.begin(), ents.end(), count, dataptr);
        assert(count == (int)ents.size());

        // Now write from memory directly
        int success = 0;
        switch (variableData.varDataType) {
          case NC_DOUBLE: {
            std::vector<double> tmpdoubledata(ni*nj*nk);
            // Transpose (lat, lon, lev) back to (lev, lat, lon)
            jik_to_kji(ni, nj, nk, &tmpdoubledata[0], (double*)(dataptr));
            success = NCFUNCAP(_vara_double)(_fileId, variableData.varId,
                      &variableData.writeStarts[0], &variableData.writeCounts[0],
                      &tmpdoubledata[0]);
            ERRORS(success, "Failed to write double data.");
            break;
          }
          default:
            ERRORR(MB_FAILURE, "Not implemented yet.");
        }
      }
    }
    else {
      int success = 0;
      switch (variableData.varDataType) {
        case NC_DOUBLE:
          success = NCFUNCAP(_vara_double)(_fileId, variableData.varId, &variableData.writeStarts[0],
                    &variableData.writeCounts[0], (double*)(variableData.memoryHogs[0]));
          ERRORS(success, "Failed to write double data.");
          break;
        default:
          ERRORR(MB_FAILURE, "Not implemented yet.");
      }
    }
  }

  return MB_SUCCESS;
}

} /* namespace moab */
