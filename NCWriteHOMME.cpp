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
    std::vector<int>& varDims = currentVarData.varDims;

    // Skip set variables, which were already processed in NCWriteHelper::collect_variable_data()
    if (WriteNC::ENTLOCSET == currentVarData.entLoc)
      continue;

    // Set up writeStarts and writeCounts (maximum number of dimensions is 3)
    currentVarData.writeStarts.resize(3);
    currentVarData.writeCounts.resize(3);
    unsigned int dim_idx = 0;

    // First: time
    if (currentVarData.has_tsteps) {
      // Non-set variables with timesteps
      // 3 dimensions like (time, lev, ncol)
      // 2 dimensions like (time, ncol)
      assert(3 == varDims.size() || 2 == varDims.size());

      // Time should be the first dimension
      assert(tDim == varDims[0]);

      currentVarData.writeStarts[dim_idx] = 0; // This value is timestep dependent, will be set later
      currentVarData.writeCounts[dim_idx] = 1;
      dim_idx++;
    }
    else {
      // Non-set variables without timesteps
      // 2 dimensions like (lev, ncol)
      // 1 dimension like (ncol)
      assert(2 == varDims.size() || 1 == varDims.size());
    }

    // Next: lev
    if (currentVarData.numLev > 0) {
      // Non-set variables with levels
      // 3 dimensions like (time, lev, ncol)
      // 2 dimensions like (lev, ncol)
      assert(3 == varDims.size() || 2 == varDims.size());

      currentVarData.writeStarts[dim_idx] = 0;
      currentVarData.writeCounts[dim_idx] = currentVarData.numLev;
      dim_idx++;
    }
    else {
      // Non-set variables without levels
      // 2 dimensions like (time, ncol)
      // 1 dimension like (ncol)
      assert(2 == varDims.size() || 1 == varDims.size());
    }

    // Finally: ncol
    switch (currentVarData.entLoc) {
      case WriteNC::ENTLOCVERT:
        // Vertices
        // Start from the first localGidVerts
        // Actually, this will be reset later for writing
        currentVarData.writeStarts[dim_idx] = localGidVertsOwned[0] - 1;
        currentVarData.writeCounts[dim_idx] = localGidVertsOwned.size();
        break;
      default:
        ERRORR(MB_FAILURE, "Unexpected entity location type for HOMME non-set variable.");
    }
    dim_idx++;

    // Get variable size
    currentVarData.sz = 1;
    for (std::size_t idx = 0; idx < dim_idx; idx++)
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

    // Assume this variable is on vertices for the time being
    switch (variableData.entLoc) {
      case WriteNC::ENTLOCVERT:
        // Vertices
        break;
      default:
        ERRORR(MB_FAILURE, "Unexpected entity location type for HOMME non-set variable.");
    }

    unsigned int num_timesteps;
    unsigned int ncol_idx = 0;
    if (variableData.has_tsteps) {
      // Non-set variables with timesteps
      // 3 dimensions like (time, lev, ncol)
      // 2 dimensions like (time, ncol)
      num_timesteps = tstep_nums.size();
      ncol_idx++;
    }
    else {
      // Non-set variables without timesteps
      // 2 dimensions like (lev, ncol)
      // 1 dimension like (ncol)
      num_timesteps = 1;
    }

    unsigned int num_lev;
    if (variableData.numLev > 0) {
      // Non-set variables with levels
      // 3 dimensions like (time, lev, ncol)
      // 2 dimensions like (lev, ncol)
      num_lev = variableData.numLev;
      ncol_idx++;
    }
    else {
      // Non-set variables without levels
      // 2 dimensions like (time, ncol)
      // 1 dimension like (ncol)
      num_lev = 1;
    }

    // At each timestep, we need to transpose tag format (ncol, lev) back
    // to NC format (lev, ncol) for writing
    for (unsigned int t = 0; t < num_timesteps; t++) {
      // We will write one time step, and count will be one; start will be different
      // Use tag_get_data instead of tag_iterate to copy tag data, as localVertsOwned
      // might not be contiguous. We should also transpose for level so that means
      // deep copy for transpose
      if (tDim == variableData.varDims[0])
        variableData.writeStarts[0] = t; // This is start for time
      std::vector<double> tag_data(num_local_verts_owned * num_lev);
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
          std::vector<double> tmpdoubledata(num_local_verts_owned * num_lev);
          if (num_lev > 1)
            // Transpose (ncol, lev) back to (lev, ncol)
            jik_to_kji(num_local_verts_owned, 1, num_lev, &tmpdoubledata[0], &tag_data[0]);

          size_t indexInDoubleArray = 0;
          size_t ic = 0;
          for (Range::pair_iterator pair_iter = localGidVertsOwned.pair_begin();
              pair_iter != localGidVertsOwned.pair_end(); ++pair_iter, ic++) {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second;
            variableData.writeStarts[ncol_idx] = (NCDF_SIZE)(starth - 1);
            variableData.writeCounts[ncol_idx] = (NCDF_SIZE)(endh - starth + 1);

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
            indexInDoubleArray += (endh - starth + 1) * num_lev;
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
