#include "NCHelperMPAS.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/FileOptions.hpp"
#include "moab/SpectralMeshTool.hpp"
#include "MBTagConventions.hpp"

#if HAVE_ZOLTAN
#include "MBZoltan.hpp"
#endif

#include <cmath>

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) {_readNC->readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
  if (err) {_readNC->readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

const int DEFAULT_MAX_EDGES_PER_CELL = 10;

// This is used for Zoltan only related stuff
#if HAVE_ZOLTAN
#define GET_VAR(name, id, dims) \
{ \
  id = -1; \
  int gvfail = NCFUNC(inq_varid)(_fileId, name, &id); \
  if (NC_NOERR == gvfail) { \
    int ndims; \
    gvfail = NCFUNC(inq_varndims)(_fileId, id, &ndims); \
    if (NC_NOERR == gvfail) { \
      dims.resize(ndims); \
      gvfail = NCFUNC(inq_vardimid)(_fileId, id, &dims[0]); \
    } \
  } \
}

#define GET_1D_DBL_VAR_RANGE(name, id, startIndex, endIndex, vals) \
{ \
  std::vector<int> dum_dims; \
  GET_VAR(name, id, dum_dims); \
  if (-1 != id) { \
    if (0 == rank) \
      std::cout << "var: " << name << " id: " << id << "\n"; \
    NCDF_SIZE ntmp; \
    int dvfail = NCFUNC(inq_dimlen)(_fileId, dum_dims[0], &ntmp); \
    if (0 == rank) \
      std::cout << "Prepare reading var " << name << " start: " << startIndex << \
        " end: " << endIndex << " Actual size: " << ntmp << "\n" ; \
    if (startIndex > (int)ntmp || endIndex > (int)ntmp || startIndex >= endIndex) { \
      std::cout << "Bad input\n"; \
      return MB_FAILURE; \
    } \
    vals.resize(endIndex - startIndex); \
    NCDF_SIZE ntmp1 = startIndex; \
    NCDF_SIZE count = endIndex - startIndex; \
    dvfail = NCFUNC(get_vara_double_all)(_fileId, id, &ntmp1, &count, &vals[0]); \
    if (NC_NOERR != dvfail) { \
      std::cout << "ReadNCDF: Problem getting variable " << name << "\n"; \
      return MB_FAILURE; \
    } \
  } \
}
#endif

NCHelperMPAS::NCHelperMPAS(ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
: UcdNCHelper(readNC, fileId, opts, fileSet)
, maxEdgesPerCell(DEFAULT_MAX_EDGES_PER_CELL)
, numCellGroups(0)
{
}

bool NCHelperMPAS::can_read_file(ReadNC* readNC)
{
  std::vector<std::string>& dimNames = readNC->dimNames;

  // If dimension name "vertexDegree" exists then it should be the MPAS grid
  if (std::find(dimNames.begin(), dimNames.end(), std::string("vertexDegree")) != dimNames.end())
    return true;

  return false;
}

ErrorCode NCHelperMPAS::init_mesh_vals()
{
  std::vector<std::string>& dimNames = _readNC->dimNames;
  std::vector<int>& dimLens = _readNC->dimLens;
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;

  ErrorCode rval;
  unsigned int idx;
  std::vector<std::string>::iterator vit;

  // Get max edges per cell reported in the MPAS file header
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "maxEdges")) != dimNames.end()) {
    idx = vit - dimNames.begin();
    maxEdgesPerCell = dimLens[idx];
    if (maxEdgesPerCell > DEFAULT_MAX_EDGES_PER_CELL) {
      ERRORR(MB_FAILURE, "maxEdgesPerCell read from the MPAS file header has exceeded DEFAULT_MAX_EDGES_PER_CELL.");
    }
  }

  // Look for time dimension
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "Time")) != dimNames.end())
    idx = vit - dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'Time' or 'time' dimension.");
  }
  tDim = idx;
  nTimeSteps = dimLens[idx];

  // Get number of cells
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nCells")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nCells' dimension.");
  }
  cDim = idx;
  nCells = dimLens[idx];

  // Get number of edges
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nEdges")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nEdges' dimension.");
  }
  eDim = idx;
  nEdges = dimLens[idx];

  // Get number of vertices
  vDim = -1;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nVertices")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nVertices' dimension.");
  }
  vDim = idx;
  nVertices = dimLens[idx];

  // Get number of levels
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nVertLevels")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nVertLevels' dimension.");
  }
  levDim = idx;
  nLevels = dimLens[idx];

  std::map<std::string, ReadNC::VarData>::iterator vmit;

  // Store time coordinate values in tVals
  if (nTimeSteps > 0) {
    if ((vmit = varInfo.find("xtime")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate("xtime", 0, nTimeSteps - 1, tVals);
      ERRORR(rval, "Trouble reading 'xtime' variable.");
    }
    else {
      // If expected time variable is not available, set dummy time coordinate values to tVals
      for (int t = 0; t < nTimeSteps; t++)
        tVals.push_back((double)t);
    }
  }

  // Determine the entity location type of a variable
  for (vmit = varInfo.begin(); vmit != varInfo.end(); ++vmit) {
    ReadNC::VarData& vd = (*vmit).second;
    if ((std::find(vd.varDims.begin(), vd.varDims.end(), vDim) != vd.varDims.end()) &&
      (std::find(vd.varDims.begin(), vd.varDims.end(), levDim) != vd.varDims.end()))
      vd.entLoc = ReadNC::ENTLOCVERT;
    else if ((std::find(vd.varDims.begin(), vd.varDims.end(), eDim) != vd.varDims.end()) &&
      (std::find(vd.varDims.begin(), vd.varDims.end(), levDim) != vd.varDims.end()))
      vd.entLoc = ReadNC::ENTLOCEDGE;
    else if ((std::find(vd.varDims.begin(), vd.varDims.end(), cDim) != vd.varDims.end()) &&
      (std::find(vd.varDims.begin(), vd.varDims.end(), levDim) != vd.varDims.end()))
      vd.entLoc = ReadNC::ENTLOCFACE;
  }

  return MB_SUCCESS;
}

// When noMesh option is used on this read, the old ReadNC class instance for last read can get out
// of scope (and deleted). The old instance initialized some variables properly when the mesh was
// created, but they are now lost. The new instance (will not create the mesh with noMesh option)
// has to restore them based on the existing mesh from last read
ErrorCode NCHelperMPAS::check_existing_mesh()
{
  Interface*& mbImpl = _readNC->mbImpl;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  bool& noMesh = _readNC->noMesh;

  if (noMesh) {
    ErrorCode rval;

    // Restore numCellGroups
    if (0 == numCellGroups) {
      Tag numCellGroupsTag;
      rval = mbImpl->tag_get_handle("__NUM_CELL_GROUPS", 1, MB_TYPE_INTEGER, numCellGroupsTag);
      if (MB_SUCCESS == rval)
        rval = mbImpl->tag_get_data(numCellGroupsTag, &_fileSet, 1, &numCellGroups);
    }

    if (localGidVerts.empty()) {
      // Get all vertices from tmp_set (it is the input set in no_mesh scenario)
      Range local_verts;
      rval = mbImpl->get_entities_by_dimension(_fileSet, 0, local_verts);
      if (MB_FAILURE == rval)
        return rval;

      if (!local_verts.empty()) {
        std::vector<int> gids(local_verts.size());

        // !IMPORTANT : this has to be the GLOBAL_ID tag
        rval = mbImpl->tag_get_data(mGlobalIdTag, local_verts, &gids[0]);
        if (MB_FAILURE == rval)
          return rval;

        // Restore localGidVerts
        std::copy(gids.rbegin(), gids.rend(), range_inserter(localGidVerts));
        nLocalVertices = localGidVerts.size();
      }
    }

    if (localGidEdges.empty()) {
      // Get all edges from _fileSet (it is the input set in no_mesh scenario)
      Range local_edges;
      rval = mbImpl->get_entities_by_dimension(_fileSet, 1, local_edges);
      if (MB_FAILURE == rval)
        return rval;

      if (!local_edges.empty()) {
        std::vector<int> gids(local_edges.size());

        // !IMPORTANT : this has to be the GLOBAL_ID tag
        rval = mbImpl->tag_get_data(mGlobalIdTag, local_edges, &gids[0]);
        if (MB_FAILURE == rval)
          return rval;

        // Restore localGidEdges
        std::copy(gids.rbegin(), gids.rend(), range_inserter(localGidEdges));
        nLocalEdges = localGidEdges.size();
      }
    }

    if (localGidCells.empty()) {
      // Get all cells from tmp_set (it is the input set in no_mesh scenario)
      Range local_cells;
      rval = mbImpl->get_entities_by_dimension(_fileSet, 2, local_cells);
      if (MB_FAILURE == rval)
        return rval;

      if (!local_cells.empty()) {
        std::vector<int> gids(local_cells.size());

        // !IMPORTANT : this has to be the GLOBAL_ID tag
        rval = mbImpl->tag_get_data(mGlobalIdTag, local_cells, &gids[0]);
        if (MB_FAILURE == rval)
          return rval;

        // Restore localGidCells
        std::copy(gids.rbegin(), gids.rend(), range_inserter(localGidCells));
        nLocalCells = localGidCells.size();

        if (numCellGroups > 1) {
          // Restore cellHandleToGlobalID
          Range::const_iterator rit;
          int i;
          for (rit = local_cells.begin(), i = 0; rit != local_cells.end(); ++rit, i++)
            cellHandleToGlobalID[*rit] = gids[i];
        }
      }
    }
  }

  return MB_SUCCESS;
}

ErrorCode NCHelperMPAS::create_mesh(Range& faces)
{
  Interface*& mbImpl = _readNC->mbImpl;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
  int& gatherSetRank = _readNC->gatherSetRank;
  bool& noMixedElements = _readNC->noMixedElements;
  bool& noEdges = _readNC->noEdges;
  DebugOutput& dbgOut = _readNC->dbgOut;

  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  bool& isParallel = _readNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _readNC->myPcomm;
    rank = myPcomm->proc_config().proc_rank();
    procs = myPcomm->proc_config().proc_size();
  }
#endif

  // Need to know whether we'll be creating gather mesh
  bool create_gathers = false;
  if (rank == gatherSetRank)
    create_gathers = true;

  // Trivial partition
  // Compute the number of local cells on this proc
  nLocalCells = int(std::floor(1.0 * nCells / procs));
  // start_cell_idx is the starting global cell index in the MPAS file for this proc
  int start_cell_idx = rank * nLocalCells;
  // iextra = # cells extra after equal split over procs
  int iextra = nCells % procs;
  if (rank < iextra)
    nLocalCells++;
  start_cell_idx += std::min(rank, iextra);

#ifdef USE_MPI
#ifdef HAVE_ZOLTAN
  int& partMethod = _readNC->partMethod;
  if (partMethod == ScdParData::RCBZOLTAN && procs >= 2) // it does not make sense to partition
    // if the number of processors is less than 2; trivial partition is good enough
  {
    // Zoltan partition using RCB; maybe more studies would be good, as to which partition
    // is better
    int temp_dim;
    MBZoltan* mbZTool = new MBZoltan(mbImpl, false, 0, NULL);

    std::vector<double> x, y, z;
    int end_cell_index = start_cell_idx + nLocalCells;
    GET_1D_DBL_VAR_RANGE("xCell", temp_dim, start_cell_idx, end_cell_index, x);
    GET_1D_DBL_VAR_RANGE("yCell", temp_dim, start_cell_idx, end_cell_index, y);
    GET_1D_DBL_VAR_RANGE("zCell", temp_dim, start_cell_idx, end_cell_index, z);

    ErrorCode rval = mbZTool->repartition(x, y, z, start_cell_idx + 1, "RCB", localGidCells);
    //delete mbZTool;
    if (rval != MB_SUCCESS) {
      std::cout << "Error in partitioning\n";
      return MB_FAILURE;
    }

    dbgOut.tprintf(1, "After partitioning, localGidCells.psize() = %d\n", (int)localGidCells.psize());
    dbgOut.tprintf(1, "                    localGidCells.size() = %d\n", (int)localGidCells.size());

    // This is important: local cells are now redistributed, so nLocalCells is different!
    nLocalCells = localGidCells.size();
  }
  else {
#endif /* this is end for HAVE_ZOLTAN */
#endif  /* if use mpi */
    /* without zoltan, only trivial partition is possible */
    start_cell_idx++; // 0 based -> 1 based
    localGidCells.insert(start_cell_idx, start_cell_idx + nLocalCells - 1);
#ifdef USE_MPI
#ifdef HAVE_ZOLTAN
  }
#endif /* end for HAVE_ZOLTAN */
#endif /* end for USE_MPI */

  // Read number of edges on each local cell, to calculate actual maxEdgesPerCell
  int nEdgesOnCellVarId;
  int success = NCFUNC(inq_varid)(_fileId, "nEdgesOnCell", &nEdgesOnCellVarId);
  ERRORS(success, "Failed to get variable id of nEdgesOnCell.");
  std::vector<int> num_edges_on_local_cells(nLocalCells);
#ifdef PNETCDF_FILE
  size_t nb_reads = localGidCells.psize();
  std::vector<int> requests(nb_reads);
  std::vector<int> statuss(nb_reads);
  size_t idxReq = 0;
#endif
  size_t indexInArray = 0;
  for (Range::pair_iterator pair_iter = localGidCells.pair_begin();
       pair_iter != localGidCells.pair_end();
       pair_iter++) {
    EntityHandle starth = pair_iter->first;
    EntityHandle endh = pair_iter->second;
    NCDF_SIZE tmp_start = (NCDF_SIZE) (starth - 1);
    NCDF_SIZE tmp_count = (NCDF_SIZE) (endh - starth + 1);

    // Do a partial read in each subrange
#ifdef PNETCDF_FILE
    success = NCFUNCREQG(_vara_int)(_fileId, nEdgesOnCellVarId, &tmp_start, &tmp_count,
                                      &(num_edges_on_local_cells[indexInArray]), &requests[idxReq++]);
#else
    success = NCFUNCAG(_vara_int)(_fileId, nEdgesOnCellVarId, &tmp_start, &tmp_count,
                                      &(num_edges_on_local_cells[indexInArray]));
#endif
    ERRORS(success, "Failed to read nEdgesOnCell data in a loop");

    // Increment the index for next subrange
    indexInArray += (endh - starth + 1);
  }

#ifdef PNETCDF_FILE
  // Wait outside the loop
  success = NCFUNC(wait_all)(_fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  // Get local maxEdgesPerCell on this proc
  int local_max_edges_per_cell = *(std::max_element(num_edges_on_local_cells.begin(), num_edges_on_local_cells.end()));
  maxEdgesPerCell = local_max_edges_per_cell;

  // If parallel, do a MPI_Allreduce to get a common global maxEdgesPerCell used across all procs
#ifdef USE_MPI
  if (procs > 1) {
    int global_max_edges_per_cell;
    ParallelComm*& myPcomm = _readNC->myPcomm;
    MPI_Allreduce(&local_max_edges_per_cell, &global_max_edges_per_cell, 1, MPI_INTEGER, MPI_MAX, myPcomm->proc_config().proc_comm());
    assert(local_max_edges_per_cell <= global_max_edges_per_cell);
    maxEdgesPerCell = global_max_edges_per_cell;
    if (0==rank)
      dbgOut.tprintf(1, "  global_max_edges_per_cell = %d\n", global_max_edges_per_cell);
  }
#endif

  // Read edges on each local cell, to get localGidEdges later
  int edgesOnCellVarId;
  success = NCFUNC(inq_varid)(_fileId, "edgesOnCell", &edgesOnCellVarId);
  ERRORS(success, "Failed to get variable id of edgesOnCell.");
  std::vector<int> edges_on_local_cells(nLocalCells * maxEdgesPerCell);
  dbgOut.tprintf(1, "   edges_on_local_cells.size() = %d\n", (int)edges_on_local_cells.size());
#ifdef PNETCDF_FILE
  idxReq = 0;
#endif
  indexInArray = 0;
  for (Range::pair_iterator pair_iter = localGidCells.pair_begin();
       pair_iter != localGidCells.pair_end();
       pair_iter++) {
    EntityHandle starth = pair_iter->first;
    EntityHandle endh = pair_iter->second;
    NCDF_SIZE tmp_starts[2] = {static_cast<NCDF_SIZE>(starth - 1), 0};
    NCDF_SIZE tmp_counts[2] = {static_cast<NCDF_SIZE>(endh - starth + 1), static_cast<NCDF_SIZE>(maxEdgesPerCell)};

    // Do a partial read in each subrange
#ifdef PNETCDF_FILE
    success = NCFUNCREQG(_vara_int)(_fileId, edgesOnCellVarId, tmp_starts, tmp_counts,
                                      &(edges_on_local_cells[indexInArray]), &requests[idxReq++]);
#else
    success = NCFUNCAG(_vara_int)(_fileId, edgesOnCellVarId, tmp_starts, tmp_counts,
                                      &(edges_on_local_cells[indexInArray]));
#endif
    ERRORS(success, "Failed to read edgesOnCell data in a loop");

    // Increment the index for next subrange
    indexInArray += (endh - starth + 1) * maxEdgesPerCell;
  }
#ifdef PNETCDF_FILE
  // Wait outside the loop
  success = NCFUNC(wait_all)(_fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  // Collect local edges
  std::sort(edges_on_local_cells.begin(), edges_on_local_cells.end());
  std::copy(edges_on_local_cells.rbegin(), edges_on_local_cells.rend(), range_inserter(localGidEdges));
  nLocalEdges = localGidEdges.size();

  dbgOut.tprintf(1, "   localGidEdges.psize() = %d\n", (int)localGidEdges.psize());
  dbgOut.tprintf(1, "   localGidEdges.size() = %d\n", (int)localGidEdges.size());

  // Read vertices on each local cell, to get localGidVerts and cell connectivity later
  int verticesOnCellVarId;
  success = NCFUNC(inq_varid)(_fileId, "verticesOnCell", &verticesOnCellVarId);
  ERRORS(success, "Failed to get variable id of verticesOnCell.");
  std::vector<int> vertices_on_local_cells(nLocalCells * maxEdgesPerCell);
#ifdef PNETCDF_FILE
  idxReq = 0;
#endif
  indexInArray = 0;
  for (Range::pair_iterator pair_iter = localGidCells.pair_begin();
       pair_iter != localGidCells.pair_end();
       pair_iter++) {
    EntityHandle starth = pair_iter->first;
    EntityHandle endh = pair_iter->second;
    NCDF_SIZE tmp_starts[2] = {static_cast<NCDF_SIZE>(starth - 1), 0};
    NCDF_SIZE tmp_counts[2] = {static_cast<NCDF_SIZE>(endh - starth + 1), maxEdgesPerCell};

    // Do a partial read in each subrange
#ifdef PNETCDF_FILE
    success = NCFUNCREQG(_vara_int)(_fileId, verticesOnCellVarId, tmp_starts, tmp_counts,
                                      &(vertices_on_local_cells[indexInArray]), &requests[idxReq++]);
#else
    success = NCFUNCAG(_vara_int)(_fileId, verticesOnCellVarId, tmp_starts, tmp_counts,
                                      &(vertices_on_local_cells[indexInArray]));
#endif
    ERRORS(success, "Failed to read verticesOnCell data in a loop");

    // Increment the index for next subrange
    indexInArray += (endh - starth + 1) * maxEdgesPerCell;
  }

#ifdef PNETCDF_FILE
  // Wait outside the loop
  success = NCFUNC(wait_all)(_fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  // To keep unsorted vertices_on_local_cells for cell connectivity, we reuse edges_on_local_cells for sorting
  std::vector<int>& vertices_on_local_cells_sorted = edges_on_local_cells;
  std::copy(vertices_on_local_cells.begin(), vertices_on_local_cells.end(), vertices_on_local_cells_sorted.begin());

  // Collect local vertices
  std::sort(vertices_on_local_cells_sorted.begin(), vertices_on_local_cells_sorted.end());
  std::copy(vertices_on_local_cells_sorted.rbegin(), vertices_on_local_cells_sorted.rend(), range_inserter(localGidVerts));
  nLocalVertices = localGidVerts.size();

  dbgOut.tprintf(1, "   localGidVerts.psize() = %d\n", (int)localGidVerts.psize());
  dbgOut.tprintf(1, "   localGidVerts.size() = %d\n", (int)localGidVerts.size());

  // Create local vertices
  // We can temporarily use the memory storage allocated before it will be populated
  EntityHandle start_vertex;
  std::vector<double*> arrays;
  Range tmp_range;
  ErrorCode rval = _readNC->readMeshIface->get_node_coords(3, nLocalVertices, 0, start_vertex, arrays,
                                                          // Might have to create gather mesh later
                                                          (create_gathers ? nLocalVertices + nVertices : nLocalVertices));
  ERRORR(rval, "Couldn't create local vertices in MPAS mesh.");
  tmp_range.insert(start_vertex, start_vertex + nLocalVertices - 1);

  // Get ptr to GID memory for local vertices
  Range local_verts_range(start_vertex, start_vertex + nLocalVertices - 1);
  int count = 0;
  void* data = NULL;
  rval = mbImpl->tag_iterate(mGlobalIdTag, local_verts_range.begin(), local_verts_range.end(), count, data);
  ERRORR(rval, "Failed to get global id tag iterator on local vertices.");
  assert(count == nLocalVertices);
  int* gid_data = (int*) data;
  std::copy(localGidVerts.begin(), localGidVerts.end(), gid_data);

  // Duplicate GID data, which will be used to resolve sharing
  if (mpFileIdTag) {
    rval = mbImpl->tag_iterate(*mpFileIdTag, local_verts_range.begin(), local_verts_range.end(), count, data);
    ERRORR(rval, "Failed to get file id tag iterator on local vertices.");
    assert(count == nLocalVertices);
    gid_data = (int*) data;
    std::copy(localGidVerts.begin(), localGidVerts.end(), gid_data);
  }

  // Create local edges
  // We can temporarily use the memory storage allocated before it will be populated
  EntityHandle start_edge;
  EntityHandle* conn_arr_edges = NULL;
  if (!noEdges) {
    rval = _readNC->readMeshIface->get_element_connect(nLocalEdges, 2, MBEDGE, 0, start_edge, conn_arr_edges,
                                                    // Might have to create gather mesh later
                                                    (create_gathers ? nLocalEdges + nEdges : nLocalEdges));
    ERRORR(rval, "Couldn't create local edges in MPAS mesh.");
    Range local_edges_range(start_edge, start_edge + nLocalEdges - 1);
    tmp_range.insert(start_edge, start_edge + nLocalEdges - 1);

    // Get ptr to GID memory for edges
    rval = mbImpl->tag_iterate(mGlobalIdTag, local_edges_range.begin(), local_edges_range.end(), count, data);
    ERRORR(rval, "Failed to get global id tag iterator on local edges.");
    assert(count == (int) nLocalEdges);
    gid_data = (int*) data;
    std::copy(localGidEdges.begin(), localGidEdges.end(), gid_data);
  }

  EntityHandle start_element = 0;

  // Used when NO_MIXED_ELEMENTS option is NOT set
  EntityHandle* conn_arr_local_cells_with_n_edges[DEFAULT_MAX_EDGES_PER_CELL + 1];
  // put in subranges; each subrange (for numEdges) will contain, in order, the
  // gids for cells with that many edges (numEdges)
  Range local_cells_with_n_edges[DEFAULT_MAX_EDGES_PER_CELL + 1];
  std::vector<int> num_edges_on_cell_groups;
  std::vector<EntityHandle> start_element_on_cell_groups;

  // Used when NO_MIXED_ELEMENTS option is set
  EntityHandle* conn_arr_local_cells = NULL;

  // Create local cells
  // We can temporarily use the memory storage allocated before it will be populated
  if (noMixedElements) {
    // Only one group of cells (each cell is represented by a polygon with maxEdgesPerCell edges)
    numCellGroups = 1;

    // Create cells for this cell group
    rval = _readNC->readMeshIface->get_element_connect(nLocalCells, maxEdgesPerCell, MBPOLYGON, 0, start_element, conn_arr_local_cells,
                                                      // Might have to create gather mesh later
                                                      (create_gathers ? nLocalCells + nCells : nLocalCells));
    ERRORR(rval, "Couldn't create local cells in MPAS mesh.");
    Range local_cell_range(start_element, start_element + nLocalCells - 1);
    tmp_range.insert(start_element, start_element + nLocalCells - 1);
    faces.insert(start_element, start_element + nLocalCells - 1);

    // Get ptr to GID memory for local cells
    rval = mbImpl->tag_iterate(mGlobalIdTag, local_cell_range.begin(), local_cell_range.end(), count, data);
    ERRORR(rval, "Failed to get global id tag iterator on local cells.");
    assert(count == nLocalCells);
    gid_data = (int*) data;
    std::copy(localGidCells.begin(), localGidCells.end(), gid_data);
  }
  else {
    // Divide local cells into groups based on the number of edges
    for (int i = 0; i < nLocalCells; i++) {
      int num_edges = num_edges_on_local_cells[i];
      local_cells_with_n_edges[num_edges].insert(localGidCells[i]); // Global cell index
      // by construction, this subrange will contain cells in order
    }

    for (int i = 3; i <= maxEdgesPerCell; i++) {
      if (local_cells_with_n_edges[i].size() > 0)
        num_edges_on_cell_groups.push_back(i);
    }
    numCellGroups = num_edges_on_cell_groups.size();

    // Create cells for each non-empty cell group
    for (int i = 0; i < numCellGroups; i++) {
      int num_edges_per_cell = num_edges_on_cell_groups[i];
      int num_group_cells = local_cells_with_n_edges[num_edges_per_cell].size();

      rval = _readNC->readMeshIface->get_element_connect(num_group_cells, num_edges_per_cell, MBPOLYGON, 0, start_element,
                                                         conn_arr_local_cells_with_n_edges[num_edges_per_cell], num_group_cells);
      ERRORR(rval, "Couldn't create local cells in MPAS mesh.");
      start_element_on_cell_groups.push_back(start_element);
      Range local_cell_range(start_element, start_element + num_group_cells - 1);
      tmp_range.insert(start_element, start_element + num_group_cells - 1);
      faces.insert(start_element, start_element + num_group_cells - 1);

      // Get ptr to gid memory for local cells
      rval = mbImpl->tag_iterate(mGlobalIdTag, local_cell_range.begin(), local_cell_range.end(), count, data);
      ERRORR(rval, "Failed to get global id tag iterator on local cells.");
      assert(count == num_group_cells);
      gid_data = (int*) data;
      std::copy(local_cells_with_n_edges[num_edges_per_cell].begin(), local_cells_with_n_edges[num_edges_per_cell].end(), gid_data);
    }
  }

  // Create a temporary map from vertex GIDs to local handles
  // Utilize the memory storage pointed by arrays[0]
  // EntityHandle* vert_gid_to_local_handle_tmp_map = (EntityHandle*) arrays[0];
  // you don't need it:
  // so vertex handle start_vertex+idx has the global ID localGidVerts[idx]
  //  if I need the vertex handle for gid, use Range::index();
  // so, from gid to eh: int index = localGidVerts.index(gid)
  // eh = start_vertex+idx;
  /*Range::const_iterator rit;
  int vert_idx;
  for (rit = localGidVerts.begin(), vert_idx = 0; rit != localGidVerts.end(); ++rit, vert_idx++)
    vert_gid_to_local_handle_tmp_map[*rit - 1] = start_vertex + vert_idx;*/

  int verticesOnEdgeVarId;
  if (!noEdges) {
    // Read vertices on each local edge, to get edge connectivity
    // Utilize the memory storage pointed by conn_arr_edges
    success = NCFUNC(inq_varid)(_fileId, "verticesOnEdge", &verticesOnEdgeVarId);
    ERRORS(success, "Failed to get variable id of verticesOnEdge.");
    int* vertices_on_local_edges = (int*) conn_arr_edges;
#ifdef PNETCDF_FILE
    nb_reads = localGidEdges.psize();
    requests.resize(nb_reads);
    statuss.resize(nb_reads);
    idxReq = 0;
#endif
    indexInArray = 0;
    for (Range::pair_iterator pair_iter = localGidEdges.pair_begin();
        pair_iter != localGidEdges.pair_end();
        pair_iter++) {
      EntityHandle starth = pair_iter->first;
      EntityHandle endh = pair_iter->second;
      NCDF_SIZE tmp_starts[2] = {static_cast<NCDF_SIZE>(starth - 1), 0};
      NCDF_SIZE tmp_counts[2] = {static_cast<NCDF_SIZE>(endh - starth + 1), 2};

      // Do a partial read in each subrange
#ifdef PNETCDF_FILE
      success = NCFUNCREQG(_vara_int)(_fileId, verticesOnEdgeVarId, tmp_starts, tmp_counts,
                                      &(vertices_on_local_edges[indexInArray]), &requests[idxReq++]);
#else
      success = NCFUNCAG(_vara_int)(_fileId, verticesOnEdgeVarId, tmp_starts, tmp_counts,
                                      &(vertices_on_local_edges[indexInArray]));
#endif
      ERRORS(success, "Failed to read verticesOnEdge data in a loop");

      // Increment the index for next subrange
      indexInArray += (endh - starth + 1) * 2;
    }

#ifdef PNETCDF_FILE
    // Wait outside the loop
    success = NCFUNC(wait_all)(_fileId, requests.size(), &requests[0], &statuss[0]);
    ERRORS(success, "Failed on wait_all.");
#endif

    // Populate connectivity for local edges
    // Convert in-place from int to EntityHandle type (backward)
    for (int edge_vert = nLocalEdges * 2 - 1; edge_vert >= 0; edge_vert--) {
      EntityHandle global_vert_id = vertices_on_local_edges[edge_vert];
      int idx_vertex = localGidVerts.index(global_vert_id);
      assert(idx_vertex != -1);
      conn_arr_edges[edge_vert] = start_vertex + idx_vertex;
    }
  }

  // Populate connectivity for local cells
  int cell_idx;
  if (noMixedElements) {
    // Set connectivity array with proper local vertices handles
    for (cell_idx = 0; cell_idx < nLocalCells; cell_idx++) {
      int num_edges = num_edges_on_local_cells[cell_idx];
      for (int i = 0; i < num_edges; i++) {
        EntityHandle global_vert_id = vertices_on_local_cells[cell_idx * maxEdgesPerCell + i];
        int idx_vertex=localGidVerts.index(global_vert_id);
        assert(idx_vertex!=-1);
        conn_arr_local_cells[cell_idx * maxEdgesPerCell + i] =  start_vertex+idx_vertex;
      }

      // Padding: fill connectivity array with last vertex handle
      if (num_edges < maxEdgesPerCell) {
        EntityHandle last_vert_id = conn_arr_local_cells[cell_idx * maxEdgesPerCell + num_edges - 1];
        for (int i = num_edges; i < maxEdgesPerCell; i++)
          conn_arr_local_cells[cell_idx * maxEdgesPerCell + i] = last_vert_id;
      }
    }
  } // if (noMixedElements)
  else {
    // Create a temporary map from cell GID to local index
    // Utilize the memory storage pointed by arrays[1]
    // again, not correct because of overflow
    // here we have to use a different map, maybe RangeMap
    /*EntityHandle* cell_gid_to_local_index_tmp_map = (EntityHandle*) arrays[1];
    for (rit = localGidCells.begin(), cell_idx = 0; rit != localGidCells.end(); ++rit, cell_idx++)
      cell_gid_to_local_index_tmp_map[*rit - 1] = cell_idx;*/

    // For each non-empty cell group, set connectivity array with proper local vertices handles
    for (int i = 0; i < numCellGroups; i++) {
      int num_edges_per_cell = num_edges_on_cell_groups[i];
      int num_cells = (int)local_cells_with_n_edges[num_edges_per_cell].size();
      start_element = start_element_on_cell_groups[i];

      for (int j = 0; j < num_cells; j++) {
        cell_idx = local_cells_with_n_edges[num_edges_per_cell][j]; // Global cell index

        if (numCellGroups > 1)
          cellHandleToGlobalID[start_element + j] = cell_idx;
        // find the index in the range: local_cells_with_n_edges[num_edges_per_cell];

        // overflow averted :(
        // cell_idx = cell_gid_to_local_index_tmp_map[cell_idx - 1]; // Local cell index
        cell_idx = localGidCells.index(cell_idx);
        assert(cell_idx!=-1);
        for (int k = 0; k < num_edges_per_cell; k++) {
          EntityHandle global_vert_id = vertices_on_local_cells[cell_idx * maxEdgesPerCell + k];
          int idx_vertex=localGidVerts.index(global_vert_id);
          assert(idx_vertex!=-1);
          conn_arr_local_cells_with_n_edges[num_edges_per_cell][j * num_edges_per_cell + k] =
              start_vertex+idx_vertex;
        }
      }
    }
  }

#ifdef PNETCDF_FILE
  nb_reads = localGidVerts.psize();
  requests.resize(nb_reads);
  statuss.resize(nb_reads);
#endif

  // Read x coordinates for local vertices
  // Note: vert_gid_to_local_handle_tmp_map is no longer used
  double* xptr = arrays[0];
  int xVertexVarId;
  success = NCFUNC(inq_varid)(_fileId, "xVertex", &xVertexVarId);
  ERRORS(success, "Failed to get variable id of xVertex.");
#ifdef PNETCDF_FILE
  idxReq = 0;
#endif
  indexInArray = 0;
  for (Range::pair_iterator pair_iter = localGidVerts.pair_begin();
       pair_iter != localGidVerts.pair_end();
       pair_iter++) {
    EntityHandle starth = pair_iter->first;
    EntityHandle endh = pair_iter->second;
    NCDF_SIZE tmp_start = (NCDF_SIZE) (starth - 1);
    NCDF_SIZE tmp_count = (NCDF_SIZE) (endh - starth + 1);

    // Do a partial read in each subrange
#ifdef PNETCDF_FILE
    success = NCFUNCREQG(_vara_double)(_fileId, xVertexVarId, &tmp_start, &tmp_count,
                                      &(xptr[indexInArray]), &requests[idxReq++]);
#else
    success = NCFUNCAG(_vara_double)(_fileId, xVertexVarId, &tmp_start, &tmp_count,
                                      &(xptr[indexInArray]));
#endif
    ERRORS(success, "Failed to read xVertex data in a loop");

    // Increment the index for next subrange
    indexInArray += (endh - starth + 1);
  }

#ifdef PNETCDF_FILE
  // Wait outside the loop
  success = NCFUNC(wait_all)(_fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  // Read y coordinates for local vertices
  // Note: cell_gid_to_local_index_tmp_ma is no longer used
  double* yptr = arrays[1];
  int yVertexVarId;
  success = NCFUNC(inq_varid)(_fileId, "yVertex", &yVertexVarId);
  ERRORS(success, "Failed to get variable id of yVertex.");
#ifdef PNETCDF_FILE
  idxReq = 0;
#endif
  indexInArray = 0;
  for (Range::pair_iterator pair_iter = localGidVerts.pair_begin();
       pair_iter != localGidVerts.pair_end();
       pair_iter++) {
    EntityHandle starth = pair_iter->first;
    EntityHandle endh = pair_iter->second;
    NCDF_SIZE tmp_start = (NCDF_SIZE) (starth - 1);
    NCDF_SIZE tmp_count = (NCDF_SIZE) (endh - starth + 1);

    // Do a partial read in each subrange
#ifdef PNETCDF_FILE
    success = NCFUNCREQG(_vara_double)(_fileId, yVertexVarId, &tmp_start, &tmp_count,
                                      &(yptr[indexInArray]), &requests[idxReq++]);
#else
    success = NCFUNCAG(_vara_double)(_fileId, yVertexVarId, &tmp_start, &tmp_count,
                                      &(yptr[indexInArray]));
#endif
    ERRORS(success, "Failed to read yVertex data in a loop");

    // Increment the index for next subrange
    indexInArray += (endh - starth + 1);
  }

#ifdef PNETCDF_FILE
  // Wait outside the loop
  success = NCFUNC(wait_all)(_fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  // Read z coordinates for local vertices
  double* zptr = arrays[2];
  int zVertexVarId;
  success = NCFUNC(inq_varid)(_fileId, "zVertex", &zVertexVarId);
  ERRORS(success, "Failed to get variable id of zVertex.");
#ifdef PNETCDF_FILE
  idxReq = 0;
#endif
  indexInArray = 0;
  for (Range::pair_iterator pair_iter = localGidVerts.pair_begin();
       pair_iter != localGidVerts.pair_end();
       pair_iter++) {
    EntityHandle starth = pair_iter->first;
    EntityHandle endh = pair_iter->second;
    NCDF_SIZE tmp_start = (NCDF_SIZE) (starth - 1);
    NCDF_SIZE tmp_count = (NCDF_SIZE) (endh - starth + 1);

    // Do a partial read in each subrange
#ifdef PNETCDF_FILE
    success = NCFUNCREQG(_vara_double)(_fileId, zVertexVarId, &tmp_start, &tmp_count,
                                      &(zptr[indexInArray]), &requests[idxReq++]);
#else
    success = NCFUNCAG(_vara_double)(_fileId, zVertexVarId, &tmp_start, &tmp_count,
                                      &(zptr[indexInArray]));
#endif
    ERRORS(success, "Failed to read zVertex data in a loop");

    // Increment the index for next subrange
    indexInArray += (endh - starth + 1);
  }

#ifdef PNETCDF_FILE
  // Wait outside the loop
  success = NCFUNC(wait_all)(_fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  // Set tag for numCellGroups
  Tag numCellGroupsTag = 0;
  rval = mbImpl->tag_get_handle("__NUM_CELL_GROUPS", 1, MB_TYPE_INTEGER, numCellGroupsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_CELL_GROUPS tag.");
  rval = mbImpl->tag_set_data(numCellGroupsTag, &_fileSet, 1, &numCellGroups);
  ERRORR(rval, "Trouble setting data for __NUM_CELL_GROUPS tag.");

  // Add new vertices, elements and edges to the file set
  rval = _readNC->mbImpl->add_entities(_fileSet, tmp_range);
  ERRORR(rval, "Couldn't add new vertices/faces/edges to file set.");

  dbgOut.tprintf(1, "   local entities in set: tmp_range.size() = %d\n",  (int)tmp_range.size());
  if (create_gathers) {
    EntityHandle gather_set;
    rval = _readNC->readMeshIface->create_gather_set(gather_set);
    ERRORR(rval, "Trouble creating gather set.");

    // Create gather vertices
    arrays.clear();
    // Don't need to specify allocation number here, because we know enough vertices were created before
    rval = _readNC->readMeshIface->get_node_coords(3, nVertices, 0, start_vertex, arrays);
    ERRORR(rval, "Couldn't create vertices in MPAS mesh for gather set.");
    Range gather_verts_range(start_vertex, start_vertex + nVertices - 1);

    xptr = arrays[0];
    yptr = arrays[1];
    zptr = arrays[2];
    NCDF_SIZE tmp_start = 0;
    NCDF_SIZE tmp_count = static_cast<NCDF_SIZE>(nVertices);

    // Read x coordinates for gather vertices
#ifdef PNETCDF_FILE
    // Enter independent I/O mode, since this read is only for the gather processor
    success = NCFUNC(begin_indep_data)(_fileId);
    ERRORS(success, "Failed to begin independent I/O mode.");
    success = NCFUNCG(_vara_double)(_fileId, xVertexVarId, &tmp_start, &tmp_count, xptr);
    ERRORS(success, "Failed to read xVertex data.");
    success = NCFUNC(end_indep_data)(_fileId);
    ERRORS(success, "Failed to end independent I/O mode.");
#else
    success = NCFUNCG(_vara_double)(_fileId, xVertexVarId, &tmp_start, &tmp_count, xptr);
    ERRORS(success, "Failed to read xVertex data.");
#endif

    // Read y coordinates for gather vertices
#ifdef PNETCDF_FILE
    // Enter independent I/O mode, since this read is only for the gather processor
    success = NCFUNC(begin_indep_data)(_fileId);
    ERRORS(success, "Failed to begin independent I/O mode.");
    success = NCFUNCG(_vara_double)(_fileId, yVertexVarId, &tmp_start, &tmp_count, yptr);
    ERRORS(success, "Failed to read yVertex data.");
    success = NCFUNC(end_indep_data)(_fileId);
    ERRORS(success, "Failed to end independent I/O mode.");
#else
    success = NCFUNCG(_vara_double)(_fileId, yVertexVarId, &tmp_start, &tmp_count, yptr);
    ERRORS(success, "Failed to read yVertex data.");
#endif

    // Read z coordinates for gather vertices
#ifdef PNETCDF_FILE
    // Enter independent I/O mode, since this read is only for the gather processor
    success = NCFUNC(begin_indep_data)(_fileId);
    ERRORS(success, "Failed to begin independent I/O mode.");
    success = NCFUNCG(_vara_double)(_fileId, zVertexVarId, &tmp_start, &tmp_count, zptr);
    ERRORS(success, "Failed to read zVertex data.");
    success = NCFUNC(end_indep_data)(_fileId);
    ERRORS(success, "Failed to end independent I/O mode.");
#else
    success = NCFUNCG(_vara_double)(_fileId, zVertexVarId, &tmp_start, &tmp_count, zptr);
    ERRORS(success, "Failed to read zVertex data.");
#endif

    // Get ptr to GID memory for gather vertices
    rval = mbImpl->tag_iterate(mGlobalIdTag, gather_verts_range.begin(), gather_verts_range.end(), count, data);
    ERRORR(rval, "Failed to get global id tag iterator on gather vertices.");
    assert(count == nVertices);
    gid_data = (int*) data;
    for (int j = 1; j <= nVertices; j++)
      gid_data[j - 1] = j;

    // Set the file id tag too, it should be bigger something not interfering with global id
    if (mpFileIdTag) {
      rval = mbImpl->tag_iterate(*mpFileIdTag, gather_verts_range.begin(), gather_verts_range.end(), count, data);
      ERRORR(rval, "Failed to get file id tag iterator on gather vertices.");
      assert(count == nVertices);
      gid_data = (int*) data;
      for (int j = 1; j <= nVertices; j++)
        gid_data[j - 1] = nVertices + j; // Bigger than global id tag
    }

    rval = mbImpl->add_entities(gather_set, gather_verts_range);
    ERRORR(rval, "Couldn't add vertices to gather set.");

    // Read number of edges on each gather cell
    std::vector<int> num_edges_on_gather_cells(nCells);
    tmp_start = 0;
    tmp_count = static_cast<NCDF_SIZE>(nCells);
#ifdef PNETCDF_FILE
    // Enter independent I/O mode, since this read is only for the gather processor
    success = NCFUNC(begin_indep_data)(_fileId);
    ERRORS(success, "Failed to begin independent I/O mode.");
    success = NCFUNCG(_vara_int)(_fileId, nEdgesOnCellVarId, &tmp_start, &tmp_count, &num_edges_on_gather_cells[0]);
    ERRORS(success, "Failed to read nEdgesOnCell data.");
    success = NCFUNC(end_indep_data)(_fileId);
    ERRORS(success, "Failed to end independent I/O mode.");
#else
    success = NCFUNCG(_vara_int)(_fileId, nEdgesOnCellVarId, &tmp_start, &tmp_count, &num_edges_on_gather_cells[0]);
    ERRORS(success, "Failed to read nEdgesOnCell data.");
#endif

    // Create gather edges
    EntityHandle* conn_arr_gather_edges;
    if (!noEdges) {
      // Don't need to specify allocation number here, because we know enough edges were created before
      rval = _readNC->readMeshIface->get_element_connect(nEdges, 2, MBEDGE, 0, start_edge, conn_arr_gather_edges);
      ERRORR(rval, "Couldn't create edges in MPAS mesh for gather set.");
      Range gather_edges_range(start_edge, start_edge + nEdges - 1);
      rval = mbImpl->add_entities(gather_set, gather_edges_range);
      ERRORR(rval, "Couldn't add edges to gather set.");

      // Read vertices on each edge
      std::vector<int> vertices_on_gather_edges(nEdges * 2);
      NCDF_SIZE tmp_starts[2] = {0, 0};
      NCDF_SIZE tmp_counts[2] = {static_cast<NCDF_SIZE>(nEdges), 2};
#ifdef PNETCDF_FILE
      // Enter independent I/O mode, since this read is only for the gather processor
      success = NCFUNC(begin_indep_data)(_fileId);
      ERRORS(success, "Failed to begin independent I/O mode.");
      success = NCFUNCG(_vara_int)(_fileId, verticesOnEdgeVarId, tmp_starts, tmp_counts, &vertices_on_gather_edges[0]);
      ERRORS(success, "Failed to read verticesOnEdge data.");
      success = NCFUNC(end_indep_data)(_fileId);
      ERRORS(success, "Failed to end independent I/O mode.");
#else
      success = NCFUNCG(_vara_int)(_fileId, verticesOnEdgeVarId, tmp_starts, tmp_counts, &vertices_on_gather_edges[0]);
      ERRORS(success, "Failed to read verticesOnEdge data.");
#endif

      std::copy(vertices_on_gather_edges.begin(), vertices_on_gather_edges.end(), conn_arr_gather_edges);
      for (int i = 0; i < 2 * nEdges; i++)
        // Connectivity array is shifted by where the gather verts start
        conn_arr_gather_edges[i] += start_vertex - 1;
    }

    // Read vertices on each gather cell (connectivity)
    std::vector<int> vertices_on_gather_cells(nCells * maxEdgesPerCell);
    NCDF_SIZE tmp_starts[2] = {0, 0};
    NCDF_SIZE tmp_counts[2] = {static_cast<NCDF_SIZE>(nCells), maxEdgesPerCell};
#ifdef PNETCDF_FILE
    // Enter independent I/O mode, since this read is only for the gather processor
    success = NCFUNC(begin_indep_data)(_fileId);
    ERRORS(success, "Failed to begin independent I/O mode.");
    success = NCFUNCG(_vara_int)(_fileId, verticesOnCellVarId, tmp_starts, tmp_counts, &vertices_on_gather_cells[0]);
    ERRORS(success, "Failed to read verticesOnCell data.");
    success = NCFUNC(end_indep_data)(_fileId);
    ERRORS(success, "Failed to end independent I/O mode.");
#else
    success = NCFUNCG(_vara_int)(_fileId, verticesOnCellVarId, tmp_starts, tmp_counts, &vertices_on_gather_cells[0]);
    ERRORS(success, "Failed to read verticesOnCell data.");
#endif

    Range gather_cells_range;
    if (noMixedElements) {
      // Create gather cells
      EntityHandle* conn_arr_gather_cells;
      // Don't need to specify allocation number here, because we know enough cells were created before
      rval = _readNC->readMeshIface->get_element_connect(nCells, maxEdgesPerCell, MBPOLYGON, 0, start_element, conn_arr_gather_cells);
      ERRORR(rval, "Couldn't create cells in MPAS mesh for gather set.");
      gather_cells_range.insert(start_element, start_element + nCells - 1);

      for (cell_idx = 0; cell_idx < nCells; cell_idx++) {
        int num_edges = num_edges_on_gather_cells[cell_idx];
        for (int i = 0; i < num_edges; i++)
          // Connectivity array is shifted by where the gather verts start
          conn_arr_gather_cells[cell_idx * maxEdgesPerCell + i] = (start_vertex - 1) + vertices_on_gather_cells[cell_idx * maxEdgesPerCell + i];

        // Padding: fill connectivity array with last vertex handle
        EntityHandle last_vert_id = conn_arr_gather_cells[cell_idx * maxEdgesPerCell + num_edges - 1];
        for (int i = num_edges; i < maxEdgesPerCell; i++)
          conn_arr_gather_cells[cell_idx * maxEdgesPerCell + i] = last_vert_id;
      }
    }
    else {
      // Divide gather cells into groups based on the number of edges
      std::vector<int> gather_cells_with_n_edges[DEFAULT_MAX_EDGES_PER_CELL + 1];
      for (int i = 0; i < nCells; i++) {
        int num_edges = num_edges_on_gather_cells[i];
        gather_cells_with_n_edges[num_edges].push_back(i + 1); // 0 based -> 1 based
      }

      // Create gather cells
      EntityHandle* conn_arr_gather_cells_with_n_edges[DEFAULT_MAX_EDGES_PER_CELL + 1];
      for (int num_edges_per_cell = 3; num_edges_per_cell <= maxEdgesPerCell; num_edges_per_cell++) {
        int num_cells = gather_cells_with_n_edges[num_edges_per_cell].size();
        if (num_cells > 0) {
          rval = _readNC->readMeshIface->get_element_connect(num_cells, num_edges_per_cell, MBPOLYGON, 0, start_element,
                                                             conn_arr_gather_cells_with_n_edges[num_edges_per_cell], num_cells);
          ERRORR(rval, "Couldn't create cells in MPAS mesh for gather set.");
          gather_cells_range.insert(start_element, start_element + num_cells - 1);
          for (int j = 0; j < num_cells; j++) {
            cell_idx = gather_cells_with_n_edges[num_edges_per_cell][j];
            cell_idx--; // 1 based -> 0 based
            for (int k = 0; k < num_edges_per_cell; k++)
              // Connectivity array is shifted by where the gather verts start
              conn_arr_gather_cells_with_n_edges[num_edges_per_cell][j * num_edges_per_cell + k] =
                (start_vertex - 1) + vertices_on_gather_cells[cell_idx * maxEdgesPerCell + k];
          }
        }
      }
    }

    rval = mbImpl->add_entities(gather_set, gather_cells_range);
    ERRORR(rval, "Couldn't add cells to gather set.");
  }

  return MB_SUCCESS;
}

ErrorCode NCHelperMPAS::read_ucd_variable_setup(std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
                                                 std::vector<ReadNC::VarData>& vdatas, std::vector<ReadNC::VarData>& vsetdatas)
{
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  std::map<std::string, ReadNC::VarData>::iterator mit;

  // If empty read them all
  if (var_names.empty()) {
    for (mit = varInfo.begin(); mit != varInfo.end(); ++mit) {
      ReadNC::VarData vd = (*mit).second;
      if (3 == vd.varDims.size()) {
        if ((std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), cDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(), vd.varDims.end(), levDim)
            != vd.varDims.end()))
          vdatas.push_back(vd); // 3D data (Time, nCells, nVertLevels) read here
        else if ((std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), eDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(), vd.varDims.end(), levDim)
            != vd.varDims.end()))
          vdatas.push_back(vd); // 3D data (Time, nEdges, nVertLevels) read here
        else if ((std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
            vd.varDims.end(), vDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(), vd.varDims.end(), levDim)
            != vd.varDims.end()))
          vdatas.push_back(vd); // 3D data (Time, nVertices, nVertLevels) read here
      }
      else if (1 == vd.varDims.size())
        vsetdatas.push_back(vd);
    }
  }
  else {
    for (unsigned int i = 0; i < var_names.size(); i++) {
      mit = varInfo.find(var_names[i]);
      if (mit != varInfo.end()) {
        ReadNC::VarData vd = (*mit).second;
        if (3 == vd.varDims.size()) {
          if ((std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
              vd.varDims.end(), cDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(), vd.varDims.end(), levDim)
              != vd.varDims.end()))
            vdatas.push_back(vd); // 3D data (Time, nCells, nVertLevels) read here
          else if ((std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
              vd.varDims.end(), eDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(), vd.varDims.end(), levDim)
              != vd.varDims.end()))
            vdatas.push_back(vd); // 3D data (Time, nEdges, nVertLevels) read here
          else if ((std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
              vd.varDims.end(), vDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(), vd.varDims.end(), levDim)
              != vd.varDims.end()))
            vdatas.push_back(vd); // 3D data (Time, nVertices, nVertLevels) read here
        }
        else if (1 == vd.varDims.size())
          vsetdatas.push_back(vd);
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

ErrorCode NCHelperMPAS::read_ucd_variable_to_nonset_allocate(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::vector<int>& dimLens = _readNC->dimLens;
  bool& noEdges = _readNC->noEdges;
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = MB_SUCCESS;

  Range* range = NULL;

  // Get vertices in set
  Range verts;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 0, verts);
  ERRORR(rval, "Trouble getting vertices in set.");
  assert("Should only have a single vertex subrange, since they were read in one shot" &&
      verts.psize() == 1);

  // Get edges in set
  Range edges;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 1, edges);
  ERRORR(rval, "Trouble getting edges in set.");

  // Get faces in set
  Range faces;
  rval = mbImpl->get_entities_by_dimension(_fileSet, 2, faces);
  ERRORR(rval, "Trouble getting faces in set.");
  // Note, for MPAS faces.psize() can be more than 1

#ifdef USE_MPI
  bool& isParallel = _readNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _readNC->myPcomm;
    rval = myPcomm->filter_pstatus(faces, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &facesOwned);
    ERRORR(rval, "Trouble getting owned faces in set.");
  }
  else
    facesOwned = faces; // not running in parallel, but still with MPI
#else
  facesOwned = faces;
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    if (noEdges && vdatas[i].entLoc == ReadNC::ENTLOCEDGE)
      continue;

    vdatas[i].numLev = nLevels;

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
      // First: Time
      vdatas[i].readStarts[t].push_back(tstep_nums[t]);
      vdatas[i].readCounts[t].push_back(1);

      // Next: nCells or nEdges or nVertices
      switch (vdatas[i].entLoc) {
        case ReadNC::ENTLOCVERT:
          // Vertices
          vdatas[i].readStarts[t].push_back(localGidVerts[0] - 1);
          vdatas[i].readCounts[t].push_back(nLocalVertices);
          range = &verts;
          break;
        case ReadNC::ENTLOCFACE:
          // Faces
          vdatas[i].readStarts[t].push_back(localGidCells[0] - 1);
          vdatas[i].readCounts[t].push_back(nLocalCells);
          range = &facesOwned;
          break;
        case ReadNC::ENTLOCEDGE:
          // Edges
          vdatas[i].readStarts[t].push_back(localGidEdges[0] - 1);
          vdatas[i].readCounts[t].push_back(nLocalEdges);
          range = &edges;
          break;
        default:
          ERRORR(MB_FAILURE, "Unexpected entity location type for MPAS non-set variable.");
          break;
      }

      // Last, numLev, even if it is 1
      vdatas[i].readStarts[t].push_back(0);
      vdatas[i].readCounts[t].push_back(vdatas[i].numLev);
      assert(vdatas[i].readStarts[t].size() == vdatas[i].varDims.size());

      // Get ptr to tag space
      if (vdatas[i].entLoc == ReadNC::ENTLOCFACE && numCellGroups > 1) {
        // For a cell variable that is NOT on one contiguous chunk of faces, defer its tag space allocation
        vdatas[i].varDatas[t] = NULL;
      }
      else {
        assert(1 == range->psize());
        void* data;
        int count;
        rval = mbImpl->tag_iterate(vdatas[i].varTags[t], range->begin(), range->end(), count, data);
        ERRORR(rval, "Failed to get tag iterator.");
        assert((unsigned)count == range->size());
        vdatas[i].varDatas[t] = data;
      }
    }

    // Calculate variable size
    std::size_t sz = 1;
    for (std::size_t idx = 0; idx != vdatas[i].readCounts[0].size(); idx++)
      sz *= vdatas[i].readCounts[0][idx];
    vdatas[i].sz = sz;
  }

  return rval;
}

#ifdef PNETCDF_FILE
ErrorCode NCHelperMPAS::read_ucd_variable_to_nonset_async(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  Interface*& mbImpl = _readNC->mbImpl;
  bool& noEdges = _readNC->noEdges;
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_ucd_variable_to_nonset_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables.");

  // Finally, read into that space
  int success;
  Range* pLocalGid = NULL;

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    if (noEdges && vdatas[i].entLoc == ReadNC::ENTLOCEDGE)
      continue;

    switch (vdatas[i].entLoc) {
      case ReadNC::ENTLOCVERT:
        pLocalGid = &localGidVerts;
        break;
      case ReadNC::ENTLOCFACE:
        pLocalGid = &localGidCells;
        break;
      case ReadNC::ENTLOCEDGE:
        pLocalGid = &localGidEdges;
        break;
      default:
        ERRORR(MB_FAILURE, "Unrecognized entity location type.");
        break;
    }

    std::size_t sz = vdatas[i].sz;

    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      // We will synchronize all these reads with the other processors,
      // so the wait will be inside this double loop; is it too much?
      size_t nb_reads = pLocalGid->psize();
      std::vector<int> requests(nb_reads), statuss(nb_reads);
      size_t idxReq = 0;

      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        case NC_DOUBLE: {
          std::vector<double> tmpdoubledata(sz);

          // In the case of ucd mesh, and on multiple proc,
          // we need to read as many times as subranges we have in the
          // localGid range;
          // basically, we have to give a different point
          // for data to start, for every subrange :(
          size_t nbDims = vdatas[i].readStarts[t].size();

          // Assume that the last dimension is for the nVertLevels
          size_t indexInDoubleArray = 0;
          size_t ic = 0;
          for (Range::pair_iterator pair_iter = pLocalGid->pair_begin();
              pair_iter != pLocalGid->pair_end();
              pair_iter++, ic++) {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second; // inclusive
            vdatas[i].readStarts[t][nbDims - 2] = (NCDF_SIZE) (starth - 1);
            vdatas[i].readCounts[t][nbDims - 2] = (NCDF_SIZE) (endh - starth + 1);

            // Do a partial read, in each subrange
            // wait outside this loop
            success = NCFUNCREQG(_vara_double)(_fileId, vdatas[i].varId,
                &(vdatas[i].readStarts[t][0]), &(vdatas[i].readCounts[t][0]),
                            &(tmpdoubledata[indexInDoubleArray]), &requests[idxReq++]);
            ERRORS(success, "Failed to read double data in loop");
            // We need to increment the index in double array for the
            // next subrange
            indexInDoubleArray += (endh - starth + 1) * 1 * vdatas[i].numLev;
          }
          assert(ic == pLocalGid->psize());

          success = NCFUNC(wait_all)(_fileId, requests.size(), &requests[0], &statuss[0]);
          ERRORS(success, "Failed on wait_all.");

          if (vdatas[i].entLoc == ReadNC::ENTLOCFACE && numCellGroups > 1) {
            // For a cell variable that is NOT on one contiguous chunk of faces, allocate tag space for
            // each cell group, and utilize cellHandleToGlobalID map to read tag data
            Range::iterator iter = facesOwned.begin();
            while (iter != facesOwned.end()) {
              int count;
              void* ptr;
              rval = mbImpl->tag_iterate(vdatas[i].varTags[t], iter, facesOwned.end(), count, ptr);
              ERRORR(rval, "Failed to get tag iterator.");

              for (int j = 0; j < count; j++) {
                int cell_idx = cellHandleToGlobalID[*(iter + j)]; // Global cell index
                cell_idx -= localGidCells[0]; // Local cell index
                for (int level = 0; level < vdatas[i].numLev; level++)
                  ((double*) ptr)[j * vdatas[i].numLev + level] = tmpdoubledata[cell_idx * vdatas[i].numLev + level];
              }

              iter += count;
            }
          }
          else {
            void* data = vdatas[i].varDatas[t];
            for (std::size_t idx = 0; idx != tmpdoubledata.size(); idx++)
              ((double*) data)[idx] = tmpdoubledata[idx];
          }

          break;
        }
        case NC_FLOAT: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        case NC_INT: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        case NC_SHORT: {
          ERRORR(MB_FAILURE, "not implemented");
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
    if (noEdges && vdatas[i].entLoc == ReadNC::ENTLOCEDGE)
      continue;

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
#else
ErrorCode NCHelperMPAS::read_ucd_variable_to_nonset(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  Interface*& mbImpl = _readNC->mbImpl;
  bool& noEdges = _readNC->noEdges;
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_ucd_variable_to_nonset_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables.");

  // Finally, read into that space
  int success;
  Range* pLocalGid = NULL;

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    if (noEdges && vdatas[i].entLoc == ReadNC::ENTLOCEDGE)
      continue;

    switch (vdatas[i].entLoc) {
      case ReadNC::ENTLOCVERT:
        pLocalGid = &localGidVerts;
        break;
      case ReadNC::ENTLOCFACE:
        pLocalGid = &localGidCells;
        break;
      case ReadNC::ENTLOCEDGE:
        pLocalGid = &localGidEdges;
        break;
      default:
        ERRORR(MB_FAILURE, "Unrecognized entity location type.");
        break;
    }

    std::size_t sz = vdatas[i].sz;

    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        case NC_DOUBLE: {
          // Copy from float case
          std::vector<double> tmpdoubledata(sz);

          // In the case of ucd mesh, and on multiple proc,
          // we need to read as many times as subranges we have in the
          // localGid range;
          // basically, we have to give a different point
          // for data to start, for every subrange :(
          size_t nbDims = vdatas[i].readStarts[t].size();

          // Assume that the last dimension is for the nVertLevels
          size_t indexInDoubleArray = 0;
          size_t ic = 0;
          for (Range::pair_iterator pair_iter = pLocalGid->pair_begin();
              pair_iter != pLocalGid->pair_end();
              pair_iter++, ic++) {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second; // Inclusive
            vdatas[i].readStarts[t][nbDims - 2] = (NCDF_SIZE) (starth - 1);
            vdatas[i].readCounts[t][nbDims - 2] = (NCDF_SIZE) (endh - starth + 1);

            success = NCFUNCAG(_vara_double)(_fileId, vdatas[i].varId,
                &(vdatas[i].readStarts[t][0]), &(vdatas[i].readCounts[t][0]),
                            &(tmpdoubledata[indexInDoubleArray]));
            ERRORS(success, "Failed to read double data in loop");
            // We need to increment the index in double array for the
            // next subrange
            indexInDoubleArray += (endh - starth + 1) * 1 * vdatas[i].numLev;
          }
          assert(ic == pLocalGid->psize());

          if (vdatas[i].entLoc == ReadNC::ENTLOCFACE && numCellGroups > 1) {
            // For a cell variable that is NOT on one contiguous chunk of faces, allocate tag space for
            // each cell group, and utilize cellHandleToGlobalID map to read tag data
            Range::iterator iter = facesOwned.begin();
            while (iter != facesOwned.end()) {
              int count;
              void* ptr;
              rval = mbImpl->tag_iterate(vdatas[i].varTags[t], iter, facesOwned.end(), count, ptr);
              ERRORR(rval, "Failed to get tag iterator.");

              for (int j = 0; j < count; j++) {
                int cell_idx = cellHandleToGlobalID[*(iter + j)]; // Global cell index
                cell_idx -= localGidCells[0]; // Local cell index
                for (int level = 0; level < vdatas[i].numLev; level++)
                  ((double*) ptr)[j * vdatas[i].numLev + level] = tmpdoubledata[cell_idx * vdatas[i].numLev + level];
              }

              iter += count;
            }
          }
          else {
            void* data = vdatas[i].varDatas[t];
            for (std::size_t idx = 0; idx != tmpdoubledata.size(); idx++)
              ((double*) data)[idx] = tmpdoubledata[idx];
          }

          break;
        }
        case NC_FLOAT: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        case NC_INT: {
          ERRORR(MB_FAILURE, "not implemented");
          break;
        }
        case NC_SHORT: {
          ERRORR(MB_FAILURE, "not implemented");
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
    if (noEdges && vdatas[i].entLoc == ReadNC::ENTLOCEDGE)
      continue;

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
#endif

} // namespace moab
