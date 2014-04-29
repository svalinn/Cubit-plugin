/*
 * NCWriteMPAS.cpp
 *
 *  Created on: April 9, 2014
 */

#include "NCWriteMPAS.hpp"
#include "moab/WriteUtilIface.hpp"
#include "MBTagConventions.hpp"

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { _writeNC->mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { _writeNC->mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

NCWriteMPAS::~NCWriteMPAS()
{
  // TODO Auto-generated destructor stub
}

ErrorCode NCWriteMPAS::collect_mesh_info()
{
  Interface*& mbImpl = _writeNC->mbImpl;
  std::vector<std::string>& dimNames = _writeNC->dimNames;
  std::vector<int>& dimLens = _writeNC->dimLens;
  Tag& mGlobalIdTag = _writeNC->mGlobalIdTag;

  ErrorCode rval;

  // Look for time dimension
  std::vector<std::string>::iterator vecIt;
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "Time")) != dimNames.end())
    tDim = vecIt - dimNames.begin();
  else if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end())
    tDim = vecIt - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'Time' or 'time' dimension.");
  }
  nTimeSteps = dimLens[tDim];

  // Get number of levels
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "nVertLevels")) != dimNames.end())
    levDim = vecIt - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nVertLevels' dimension.");
  }
  nLevels = dimLens[levDim];

  Range local_verts;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 0, local_verts);
  ERRORR(rval, "Trouble getting local vertices in current file set.");
  assert(!local_verts.empty());

  // Depends on whether NO_EDGES read option is set or not
  Range local_edges;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 1, local_edges);
  ERRORR(rval, "Trouble getting local edges in current file set.");
  noEdges = local_edges.empty();

  Range local_cells;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 2, local_cells);
  ERRORR(rval, "Trouble getting local cells in current file set.");
  assert(!local_cells.empty());

#ifdef USE_MPI
  bool& isParallel = _writeNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _writeNC->myPcomm;
    int rank = myPcomm->proc_config().proc_rank();
    int procs = myPcomm->proc_config().proc_size();
    if (procs > 1) {
      rval = myPcomm->filter_pstatus(local_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &localVertsOwned);
      ERRORR(rval, "Trouble getting owned vertices in set.");
      // Assume that PARALLEL_RESOLVE_SHARED_ENTS option is set
      // We should avoid writing in parallel with overlapped data
      if (procs - 1 == rank)
        assert("PARALLEL_RESOLVE_SHARED_ENTS option is set" && localVertsOwned.size() < local_verts.size());

      if (!noEdges) {
        rval = myPcomm->filter_pstatus(local_edges, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &localEdgesOwned);
        ERRORR(rval, "Trouble getting owned edges in set.");
      }

      rval = myPcomm->filter_pstatus(local_cells, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &localCellsOwned);
      ERRORR(rval, "Trouble getting owned cells in set.");
    }
    else {
      localVertsOwned = local_verts;
      localEdgesOwned = local_edges;
      localCellsOwned = local_cells;
    }
  }
  else {
    // Not running in parallel, but still with MPI
    localVertsOwned = local_verts;
    localEdgesOwned = local_edges;
    localCellsOwned = local_cells;
  }
#else
  localVertsOwned = local_verts;
  localEdgesOwned = local_edges;
  localCellsOwned = local_cells;
#endif

  std::vector<int> gids(localVertsOwned.size());
  rval = mbImpl->tag_get_data(mGlobalIdTag, localVertsOwned, &gids[0]);
  ERRORR(rval, "Trouble getting global IDs on local vertices.");

  // Restore localGidVerts
  std::copy(gids.rbegin(), gids.rend(), range_inserter(localGidVertsOwned));
  nLocalVerticesOwned = localGidVertsOwned.size();

  gids.resize(localEdgesOwned.size());
  rval = mbImpl->tag_get_data(mGlobalIdTag, localEdgesOwned, &gids[0]);
  ERRORR(rval, "Trouble getting global IDs on local edges.");

  // Restore localGidEdges
  std::copy(gids.rbegin(), gids.rend(), range_inserter(localGidEdgesOwned));
  nLocalEdgesOwned = localGidEdgesOwned.size();

  gids.resize(localCellsOwned.size());
  rval = mbImpl->tag_get_data(mGlobalIdTag, localCellsOwned, &gids[0]);
  ERRORR(rval, "Trouble getting global IDs on local cells.");

  // Restore localGidCells
  std::copy(gids.rbegin(), gids.rend(), range_inserter(localGidCellsOwned));
  nLocalCellsOwned = localGidCellsOwned.size();

  return MB_SUCCESS;
}

ErrorCode NCWriteMPAS::collect_variable_data(std::vector<std::string>& var_names)
{
  NCWriteHelper::collect_variable_data(var_names);

  std::vector<std::string>& dimNames = _writeNC->dimNames;
  std::vector<int>& dimLens = _writeNC->dimLens;

  // Dimension numbers for other optional levels
  std::vector<unsigned int> opt_lev_dims;

  unsigned int lev_idx;
  std::vector<std::string>::iterator vecIt;

  // Get number of vertex levels P1
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "nVertLevelsP1")) != dimNames.end()) {
    lev_idx = vecIt - dimNames.begin();
    opt_lev_dims.push_back(lev_idx);
  }

  // Get number of vertex levels P2
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "nVertLevelsP2")) != dimNames.end()) {
    lev_idx = vecIt - dimNames.begin();
    opt_lev_dims.push_back(lev_idx);
  }

  // Get number of soil levels
  if ((vecIt = std::find(dimNames.begin(), dimNames.end(), "nSoilLevels")) != dimNames.end()) {
    lev_idx = vecIt - dimNames.begin();
    opt_lev_dims.push_back(lev_idx);
  }

  std::map<std::string, WriteNC::VarData>& varInfo = _writeNC->varInfo;

  for (size_t i = 0; i < var_names.size(); i++) {
    std::string varname = var_names[i];
    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(varname);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find one variable.");

    WriteNC::VarData& currentVarData = vit->second;
    std::vector<int>& varDims = currentVarData.varDims;

    // If nVertLevels dimension is not found, try other optional levels such as nVertLevelsP1
    if (std::find(varDims.begin(), varDims.end(), levDim) == varDims.end()) {
      for (unsigned int j = 0; j < opt_lev_dims.size(); j++) {
        if (std::find(varDims.begin(), varDims.end(), opt_lev_dims[j]) != varDims.end()) {
          currentVarData.numLev = dimLens[opt_lev_dims[j]];
          break;
        }
      }
    }

    if (currentVarData.has_tsteps) {
      // Support non-set variables with 3 dimensions like (Time, nCells, nVertLevels), or
      // 2 dimensions like (Time, nCells)
      assert(3 == currentVarData.varDims.size() || 2 == currentVarData.varDims.size());

      // Time should be the first dimension
      assert(tDim == currentVarData.varDims[0]);

      // Set up writeStarts and writeCounts
      currentVarData.writeStarts.resize(3);
      currentVarData.writeCounts.resize(3);

      // First: Time
      currentVarData.writeStarts[0] = 0; // This value is timestep dependent, will be set later
      currentVarData.writeCounts[0] = 1;

      // Next: nVertices / nCells / nEdges
      switch (currentVarData.entLoc) {
        case WriteNC::ENTLOCVERT:
          // Vertices
          // Start from the first localGidVerts
          // Actually, this will be reset later for writing
          currentVarData.writeStarts[1] = localGidVertsOwned[0] - 1;
          currentVarData.writeCounts[1] = nLocalVerticesOwned;
          break;
        case WriteNC::ENTLOCFACE:
          // Faces
          // Start from the first localGidCells
          // Actually, this will be reset later for writing
          currentVarData.writeStarts[1] = localGidCellsOwned[0] - 1;
          currentVarData.writeCounts[1] = nLocalCellsOwned;
          break;
        case WriteNC::ENTLOCEDGE:
          // Edges
          // Start from the first localGidEdges
          // Actually, this will be reset later for writing
          currentVarData.writeStarts[1] = localGidEdgesOwned[0] - 1;
          currentVarData.writeCounts[1] = nLocalEdgesOwned;
          break;
        default:
          ERRORR(MB_FAILURE, "Unexpected entity location type for MPAS non-set variable.");
      }

      // Finally: nVertLevels or other optional levels, it is possible that there
      // is no level dimension for this non-set variable
      currentVarData.writeStarts[2] = 0;
      currentVarData.writeCounts[2] = currentVarData.numLev;
    }

    // Get variable size
    currentVarData.sz = 1;
    for (std::size_t idx = 0; idx != currentVarData.writeCounts.size(); idx++)
      currentVarData.sz *= currentVarData.writeCounts[idx];
  }

  return MB_SUCCESS;
}

ErrorCode NCWriteMPAS::write_values(std::vector<std::string>& var_names)
{
  Interface*& mbImpl = _writeNC->mbImpl;
  std::set<std::string>& usedCoordinates = _writeNC->usedCoordinates;
  std::set<std::string>& dummyVarNames = _writeNC->dummyVarNames;
  std::map<std::string, WriteNC::VarData>& varInfo = _writeNC->varInfo;

  int success;
  Range* pLocalEntsOwned = NULL;
  Range* pLocalGidEntsOwned = NULL;

  // Now look at requested var_names; if they have time, we will have a list, and write one at a time
  // For each variable tag in the indexed lists, write a time step data
  // Assume the first dimension is time (need to check); if not, just write regularly
  for (size_t i = 0; i < var_names.size(); i++) {
    std::map<std::string, WriteNC::VarData>::iterator vit = varInfo.find(var_names[i]);
    if (vit == varInfo.end())
      ERRORR(MB_FAILURE, "Can't find variable requested.");

    WriteNC::VarData& variableData = vit->second;
    int numTimeSteps = (int)variableData.varTags.size();
    if (variableData.has_tsteps) {
      // Time should be the first dimension
      assert(tDim == variableData.varDims[0]);

      // Assume this variable is on vertices for the time being
      switch (variableData.entLoc) {
        case WriteNC::ENTLOCVERT:
          // Vertices
          pLocalEntsOwned = &localVertsOwned;
          pLocalGidEntsOwned = &localGidVertsOwned;
          break;
        case WriteNC::ENTLOCEDGE:
          // Edges
          pLocalEntsOwned = &localEdgesOwned;
          pLocalGidEntsOwned = &localGidEdgesOwned;
          break;
        case WriteNC::ENTLOCFACE:
          // Cells
          pLocalEntsOwned = &localCellsOwned;
          pLocalGidEntsOwned = &localGidCellsOwned;
          break;
        default:
          ERRORR(MB_FAILURE, "Unexpected entity location type for MPAS non-set variable.");
      }

      // A typical variable has 3 dimensions as (Time, nCells, nVertLevels)
      // FIXME: Should use tstep_nums (from writing options) later
      for (int j = 0; j < numTimeSteps; j++) {
        // We will write one time step, and count will be one; start will be different
        // Use tag_get_data instead of tag_iterate to get values, as localEntsOwned
        // might not be contiguous.
        variableData.writeStarts[0] = j; // This is time, again
        std::vector<double> tag_data(pLocalEntsOwned->size() * variableData.numLev);
        ErrorCode rval = mbImpl->tag_get_data(variableData.varTags[j], *pLocalEntsOwned, &tag_data[0]);
        ERRORR(rval, "Trouble getting tag data on owned vertices.");

#ifdef PNETCDF_FILE
        size_t nb_writes = pLocalGidEntsOwned->psize();
        std::vector<int> requests(nb_writes), statuss(nb_writes);
        size_t idxReq = 0;
#endif

        // Now write from memory directly
        switch (variableData.varDataType) {
          case NC_DOUBLE: {
            size_t indexInDoubleArray = 0;
            size_t ic = 0;
            for (Range::pair_iterator pair_iter = pLocalGidEntsOwned->pair_begin();
                pair_iter != pLocalGidEntsOwned->pair_end(); ++pair_iter, ic++) {
              EntityHandle starth = pair_iter->first;
              EntityHandle endh = pair_iter->second;
              variableData.writeStarts[1] = (NCDF_SIZE)(starth - 1);
              variableData.writeCounts[1] = (NCDF_SIZE)(endh - starth + 1);

              // Do a partial write, in each subrange
#ifdef PNETCDF_FILE
              // Wait outside this loop
              success = NCFUNCREQP(_vara_double)(_fileId, variableData.varId,
                  &(variableData.writeStarts[0]), &(variableData.writeCounts[0]),
                             &(tag_data[indexInDoubleArray]), &requests[idxReq++]);
#else
              success = NCFUNCAP(_vara_double)(_fileId, variableData.varId,
                  &(variableData.writeStarts[0]), &(variableData.writeCounts[0]),
                             &(tag_data[indexInDoubleArray]));
#endif
              ERRORS(success, "Failed to read double data in loop");
              // We need to increment the index in double array for the
              // next subrange
              indexInDoubleArray += (endh - starth + 1) * variableData.numLev;
            }
            assert(ic == pLocalGidEntsOwned->psize());
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

  // Write coordinates used by requested var_names
  // Use independent I/O mode put, since this write is only for the root processor
  // CAUTION: if the NetCDF ID is from a previous call to ncmpi_create rather than ncmpi_open,
  // all processors need to call ncmpi_begin_indep_data(). If only the root processor does so,
  // ncmpi_begin_indep_data() call will be blocked forever :(
#ifdef PNETCDF_FILE
  // Enter independent I/O mode
  success = NCFUNC(begin_indep_data)(_fileId);
  ERRORS(success, "Failed to begin independent I/O mode.");
#endif

  int rank = 0;
#ifdef USE_MPI
  bool& isParallel = _writeNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _writeNC->myPcomm;
    rank = myPcomm->proc_config().proc_rank();
  }
#endif
  if (0 == rank) {
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

      switch (varCoordData.varDataType) {
        case NC_DOUBLE:
          // Independent I/O mode put
          success = NCFUNCP(_vara_double)(_fileId, varCoordData.varId, &varCoordData.writeStarts[0],
                    &varCoordData.writeCounts[0], (double*)(varCoordData.memoryHogs[0]));
          ERRORS(success, "Failed to write double data.");
          break;
        case NC_INT:
          // Independent I/O mode put
          success = NCFUNCP(_vara_int)(_fileId, varCoordData.varId, &varCoordData.writeStarts[0],
                    &varCoordData.writeCounts[0], (int*)(varCoordData.memoryHogs[0]));
          ERRORS(success, "Failed to write int data.");
          break;
        default:
          ERRORR(MB_FAILURE, "Not implemented yet.");
      }
    }
  }

#ifdef PNETCDF_FILE
  // End independent I/O mode
  success = NCFUNC(end_indep_data)(_fileId);
  ERRORS(success, "Failed to end independent I/O mode.");
#endif

  return MB_SUCCESS;
}

} /* namespace moab */
