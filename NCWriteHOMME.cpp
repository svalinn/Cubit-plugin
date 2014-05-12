/*
 * NCWriteHOMME.cpp
 *
 *  Created on: April 9, 2014
 */

#include "NCWriteHOMME.hpp"
#include "moab/WriteUtilIface.hpp"
#include "MBTagConventions.hpp"

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

  // Get local vertices
  rval = mbImpl->get_entities_by_dimension(_fileSet, 0, localVertsOwned);
  ERRORR(rval, "Trouble getting local vertices in current file set.");
  assert(!localVertsOwned.empty());

#ifdef USE_MPI
  bool& isParallel = _writeNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _writeNC->myPcomm;
    int rank = myPcomm->proc_config().proc_rank();
    int procs = myPcomm->proc_config().proc_size();
    if (procs > 1) {
      unsigned int num_local_verts = localVertsOwned.size();
      rval = myPcomm->filter_pstatus(localVertsOwned, PSTATUS_NOT_OWNED, PSTATUS_NOT);
      ERRORR(rval, "Trouble getting owned vertices in set.");

      // Assume that PARALLEL_RESOLVE_SHARED_ENTS option is set
      // Verify that not all local vertices are owned by the last processor
      if (procs - 1 == rank)
        assert("PARALLEL_RESOLVE_SHARED_ENTS option is set" && localVertsOwned.size() < num_local_verts);
    }
  }
#endif

  std::vector<int> gids(localVertsOwned.size());
  rval = mbImpl->tag_get_data(mGlobalIdTag, localVertsOwned, &gids[0]);
  ERRORR(rval, "Trouble getting global IDs on local vertices.");

  // Get localGidVertsOwned
  std::copy(gids.rbegin(), gids.rend(), range_inserter(localGidVertsOwned));

  return MB_SUCCESS;
}

ErrorCode NCWriteHOMME::collect_variable_data(std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
{
  NCWriteHelper::collect_variable_data(var_names, tstep_nums);

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
          currentVarData.writeStarts[2] = localGidVertsOwned[0] - 1;
          currentVarData.writeCounts[2] = localGidVertsOwned.size();
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

ErrorCode NCWriteHOMME::write_nonset_variables(std::vector<WriteNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  Interface*& mbImpl = _writeNC->mbImpl;

  int success;
  int num_local_verts_owned = localVertsOwned.size();

  // For each indexed variable tag, write a time step data
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    WriteNC::VarData& variableData = vdatas[i];

    // Time should be the first dimension
    assert(tDim == variableData.varDims[0]);

    // Assume this variable is on vertices for the time being
    switch (variableData.entLoc) {
      case WriteNC::ENTLOCVERT:
        // Vertices
        break;
      default:
        ERRORR(MB_FAILURE, "Unexpected entity location type for HOMME non-set variable.");
    }

    // A typical HOMME non-set variable has 3 dimensions as (time, lev, ncol)
    // At each timestep, we need to transpose tag format (ncol, lev) back
    // to NC format (lev, ncol) for writing
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      // We will write one time step, and count will be one; start will be different
      // Use tag_get_data instead of tag_iterate to copy tag data, as localVertsOwned
      // might not be contiguous. We should also transpose for level so that means
      // deep copy for transpose
      variableData.writeStarts[0] = t; // This is start for time
      std::vector<double> tag_data(num_local_verts_owned * variableData.numLev);
      ErrorCode rval = mbImpl->tag_get_data(variableData.varTags[t], localVertsOwned, &tag_data[0]);
      ERRORR(rval, "Trouble getting tag data on owned vertices.");

#ifdef PNETCDF_FILE
      size_t nb_writes = localGidVertsOwned.psize();
      std::vector<int> requests(nb_writes), statuss(nb_writes);
      size_t idxReq = 0;
#endif

      // Now transpose and write copied tag data
      // Use nonblocking put (request aggregation)
      switch (variableData.varDataType) {
        case NC_DOUBLE: {
          std::vector<double> tmpdoubledata(num_local_verts_owned * variableData.numLev);
          // Transpose (ncol, lev) back to (lev, ncol)
          jik_to_kji(num_local_verts_owned, 1, variableData.numLev, &tmpdoubledata[0], &tag_data[0]);

          size_t indexInDoubleArray = 0;
          size_t ic = 0;
          for (Range::pair_iterator pair_iter = localGidVertsOwned.pair_begin();
              pair_iter != localGidVertsOwned.pair_end(); ++pair_iter, ic++) {
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
          assert(ic == localGidVertsOwned.psize());
#ifdef PNETCDF_FILE
          success = ncmpi_wait_all(_fileId, requests.size(), &requests[0], &statuss[0]);
          ERRORS(success, "Failed on wait_all.");
#endif
          break;
        }
        default:
          ERRORR(MB_NOT_IMPLEMENTED, "Writing with current data type not implemented yet.");
      }
    }
  }

  return MB_SUCCESS;
}

} /* namespace moab */
