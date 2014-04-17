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

ErrorCode NCWriteHOMME::collect_mesh_info()
{
  Interface*& mbImpl = _writeNC->mbImpl;
  std::vector<std::string>& dimNames = _writeNC->dimNames;
  std::vector<int>& dimLens = _writeNC->dimLens;
  Tag& mGlobalIdTag = _writeNC->mGlobalIdTag;

  ErrorCode rval;

  // Look for time dimension
  std::vector<std::string>::iterator vecIt;
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end())
    tDim = vecIt - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'time' dimension.");
  }
  nTimeSteps = dimLens[tDim];

  // Get number of levels
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "lev")) != dimNames.end())
    levDim = vecIt - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'lev' dimension.");
  }
  nLevels = dimLens[levDim];

  Range local_verts;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 0, local_verts);
  ERRORR(rval, "Trouble getting local vertices in current file set.");
  assert(!local_verts.empty());

  std::vector<int> gids(local_verts.size());
  rval = mbImpl->tag_get_data(mGlobalIdTag, local_verts, &gids[0]);
  ERRORR(rval, "Trouble getting global IDs on local vertices.");

  // Restore localGidVerts
  std::copy(gids.rbegin(), gids.rend(), range_inserter(localGidVerts));
  nLocalVertices = localGidVerts.size();

  return MB_SUCCESS;
}

ErrorCode NCWriteHOMME::collect_variable_data(std::vector<std::string>& var_names)
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
      // Support non-set variables with 3 dimensions like (time, lev, ncol)
      assert(3 == currentVarData.varDims.size());

      // Time should be the first dimension
      assert(tDim == currentVarData.varDims[0]);

      // Set up writeStarts and writeCounts
      currentVarData.writeStarts.resize(3);
      currentVarData.writeCounts.resize(3);

      // First: time
      currentVarData.writeStarts[0] = 0; // This value is timestep dependent, will be set later
      currentVarData.writeCounts[0] = 1;

      // Next: lev
      currentVarData.writeStarts[1] = 0;
      currentVarData.writeCounts[1] = currentVarData.numLev;

      // Finally: ncol
      switch (currentVarData.entLoc) {
        case WriteNC::ENTLOCVERT:
          // Vertices
          // Start from the first localGidVerts
          // Actually, this will be reset later for writing
          currentVarData.writeStarts[2] = localGidVerts[0] - 1;
          currentVarData.writeCounts[2] = nLocalVertices;
          break;
        default:
          ERRORR(MB_FAILURE, "Unexpected entity location type for HOMME non-set variable.");
      }
    }

    // Get variable size
    currentVarData.sz = 1;
    for (std::size_t idx = 0; idx != currentVarData.writeCounts.size(); idx++)
      currentVarData.sz *= currentVarData.writeCounts[idx];
  }

  return MB_SUCCESS;
}

ErrorCode NCWriteHOMME::write_values(std::vector<std::string>& var_names)
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
        case WriteNC::ENTLOCVERT:
          // Vertices
          rval = mbImpl->get_entities_by_dimension(_fileSet, 0, ents);
          ERRORR(rval, "Can't get entities for vertices.");
          break;
        default:
          ERRORR(MB_FAILURE, "Unexpected entity location type for HOMME non-set variable.");
      }

      // A typical variable has 3 dimensions as (time, lev, ncol)
      // At each timestep, we need to transpose tag format (ncol, lev) back
      // to NC format (lev, ncol) for writing
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

#ifdef PNETCDF_FILE
        size_t nb_writes = localGidVerts.psize();
        std::vector<int> requests(nb_writes), statuss(nb_writes);
        size_t idxReq = 0;
#endif

        // Now write from memory directly
        int success = 0;
        switch (variableData.varDataType) {
          case NC_DOUBLE: {
            std::vector<double> tmpdoubledata(nLocalVertices * variableData.numLev);
            // Transpose (ncol, lev) back to (lev, ncol)
            jik_to_kji(nLocalVertices, 1, variableData.numLev, &tmpdoubledata[0], (double*)(dataptr));

            size_t indexInDoubleArray = 0;
            size_t ic = 0;
            for (Range::pair_iterator pair_iter = localGidVerts.pair_begin();
                pair_iter != localGidVerts.pair_end(); ++pair_iter, ic++) {
              EntityHandle starth = pair_iter->first;
              EntityHandle endh = pair_iter->second;
              variableData.writeStarts[2] = (NCDF_SIZE)(starth - 1);
              variableData.writeCounts[2] = (NCDF_SIZE)(endh - starth + 1);

              // Do a partial write, in each subrange
#ifdef PNETCDF_FILE
              // Wait outside this loop
              success = NCFUNCREQP(_vara_double)(_fileId, variableData.varId,
                  &(variableData.writeStarts[0]), &(variableData.writeCounts[0]),
                             &(tmpdoubledata[indexInDoubleArray]), &requests[idxReq++]);
#else
              success = NCFUNCAP(_vara_double)(_fileId, variableData.varId,
                  &(variableData.writeStarts[0]), &(variableData.writeCounts[0]),
                             &(tmpdoubledata[indexInDoubleArray]));
#endif
              ERRORS(success, "Failed to read double data in loop");
              // We need to increment the index in double array for the
              // next subrange
              indexInDoubleArray += (endh - starth + 1) * variableData.numLev;
            }
            assert(ic == localGidVerts.psize());
#ifdef PNETCDF_FILE
            success = ncmpi_wait_all(_fileId, requests.size(), &requests[0], &statuss[0]);
            ERRORS(success, "Failed on wait_all.");
#endif
            break;
          }
          default:
            ERRORR(MB_FAILURE, "Not implemented yet.");
        }
      }
    } // if (variableData.has_tsteps)
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
