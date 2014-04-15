/*
 * NCWriteHOMME.cpp
 *
 *  Created on: April 9, 2014
 */

#include "NCWriteHOMME.hpp"
#include "moab/WriteUtilIface.hpp"

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { _writeNC->mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { _writeNC->mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

NCWriteHOMME::~NCWriteHOMME()
{
  // TODO Auto-generated destructor stub
}

ErrorCode NCWriteHOMME::write_values(std::vector<std::string>& var_names, EntityHandle fileSet)
{
  Interface*& mbImpl = _writeNC->mbImpl;
  std::set<std::string>& usedCoordinates = _writeNC->usedCoordinates;
  std::set<std::string>& dummyVarNames = _writeNC->dummyVarNames;
  std::map<std::string, WriteNC::VarData>& varInfo = _writeNC->varInfo;

  // Start with coordinates
  for (std::set<std::string>::iterator setIt = usedCoordinates.begin();
      setIt != usedCoordinates.end(); ++setIt) {
    std::string coordName = *setIt; // Deep copy

    // Skip dummy coordinate variables (e.g. ncol)
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
  // Need to transpose from lev dimension

  // Write data, to be implemented ...

  return MB_SUCCESS;
}

} /* namespace moab */
