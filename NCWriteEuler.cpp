/*
 * NCWriteEuler.cpp
 *
 *  Created on: Mar 28, 2014
 */

#include "NCWriteEuler.hpp"
#include "moab/WriteUtilIface.hpp"

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { _writeNC->mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { _writeNC->mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

NCWriteEuler::~NCWriteEuler()
{
  // TODO Auto-generated destructor stub
}

ErrorCode NCWriteEuler::write_values(std::vector<std::string>& var_names, EntityHandle fileSet)
{
  Interface*& mbImpl = _writeNC->mbImpl;
  std::set<std::string>& usedCoordinates = _writeNC->usedCoordinates;
  std::set<std::string>& dummyVarNames = _writeNC->dummyVarNames;
  std::map<std::string, WriteNC::VarData>& varInfo = _writeNC->varInfo;

  // Start with coordinates
  for (std::set<std::string>::iterator setIt = usedCoordinates.begin();
      setIt != usedCoordinates.end(); ++setIt) {
    std::string coordName = *setIt; // Deep copy

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
  // Need to transpose from lev dimension
  ErrorCode rval;

  // For each variable tag in the indexed lists, write a time step data
  // Assume the first dimension is time (need to check); if not, just write regularly
  for (size_t i = 0; i < var_names.size(); i++) {
    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(var_names[i]);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find variable requested.");

    WriteNC::VarData& variableData = vit->second;
    int numTimeSteps = (int)variableData.varTags.size();
    if (variableData.has_tsteps) {
      // Get entities of this variable
      Range ents;
      switch (variableData.entLoc) {
        case WriteNC::ENTLOCFACE:
          // Faces
          rval = mbImpl->get_entities_by_dimension(fileSet, 2, ents);
          ERRORR(rval, "Can't get entities for faces.");
          break;
        default:
          ERRORR(MB_FAILURE, "Unexpected entity location type for CAM-EUL non-set variable.");
      }

      // A typical variable has 4 dimensions as (time, lev, lat, lon)
      // At each timestep, we need to transpose tag format (lat, lon, lev) back
      // to NC format (lev, lat, lon) for writing
      size_t ni = variableData.writeCounts[3]; // lon
      size_t nj = variableData.writeCounts[2]; // lat
      size_t nk = variableData.writeCounts[1]; // lev

      variableData.writeCounts[0] = 1; // We will write one time step
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
          case NC_INT: {
            std::vector<int> tmpintdata(ni*nj*nk);
            // Transpose (lat, lon, lev) back to (lev, lat, lon)
            jik_to_kji(ni, nj, nk, &tmpintdata[0], (int*)(dataptr));
            success = NCFUNCAP(_vara_int)(_fileId, variableData.varId,
                      &variableData.writeStarts[0], &variableData.writeCounts[0],
                      &tmpintdata[0]);
            ERRORS(success, "Failed to write int data.");
            break;
          }
          default:
            success = 1;
            break;
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
        case NC_INT:
          success = NCFUNCAP(_vara_int)(_fileId, variableData.varId, &variableData.writeStarts[0],
                    &variableData.writeCounts[0], (int*)(variableData.memoryHogs[0]));
          ERRORS(success, "Failed to write int data.");
          break;
        default:
          success = 1;
          break;
      }
    }
  }

  return MB_SUCCESS;
}

} /* namespace moab */
