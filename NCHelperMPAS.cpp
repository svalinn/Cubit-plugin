#include "NCHelperMPAS.hpp"
#include "moab/ReadUtilIface.hpp"
#include "FileOptions.hpp"
#include "moab/SpectralMeshTool.hpp"
#include "MBTagConventions.hpp"

#include <cmath>

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) {_readNC->readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
  if (err) {_readNC->readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

const int MAX_EDGES_PER_CELL = 10;

NCHelperMPAS::NCHelperMPAS(ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
: UcdNCHelper(readNC, fileId, opts, fileSet)
, maxCellEdges(MAX_EDGES_PER_CELL)
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
  std::vector<int>& dimVals = _readNC->dimVals;
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;

  ErrorCode rval;
  unsigned int idx;
  std::vector<std::string>::iterator vit;

  // Get max edges per cell
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "maxEdges")) != dimNames.end()) {
    idx = vit - dimNames.begin();
    maxCellEdges = dimVals[idx];
    if (maxCellEdges > MAX_EDGES_PER_CELL) {
      ERRORR(MB_FAILURE, "maxEdges read from MPAS file has exceeded the limit");
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
  nTimeSteps = dimVals[idx];

  // Get number of cells
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nCells")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nCells' dimension.");
  }
  cDim = idx;
  nCells = dimVals[idx];

  // Get number of edges
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nEdges")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nEdges' dimension.");
  }
  eDim = idx;
  nEdges = dimVals[idx];

  // Get number of vertices
  vDim = -1;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nVertices")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nVertices' dimension.");
  }
  vDim = idx;
  nVertices = dimVals[idx];

  // Get number of levels
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "nVertLevels")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'nVertLevels' dimension.");
  }
  levDim = idx;
  nLevels = dimVals[idx];

  // Store xVertex values in xVertVals
  std::map<std::string, ReadNC::VarData>::iterator vmit;
  if ((vmit = varInfo.find("xVertex")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate("xVertex", 0, nVertices - 1, xVertVals);
    ERRORR(rval, "Trouble reading 'xVertex' variable.");
  }
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'xVertex' variable.");
  }

  // Store yVertex values in yVertVals
  if ((vmit = varInfo.find("yVertex")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate("yVertex", 0, nVertices - 1, yVertVals);
    ERRORR(rval, "Trouble reading 'yVertex' variable.");
  }
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'yVertex' variable.");
  }

  // Store zVertex values in zVertVals
  if ((vmit = varInfo.find("zVertex")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate("zVertex", 0, nVertices - 1, zVertVals);
    ERRORR(rval, "Trouble reading 'zVertex' variable.");
  }
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'zVertex' variable.");
  }

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

  // Read vertices on each edge
  int verticesOnEdgeVarId;
  int success = NCFUNC(inq_varid)(_fileId, "verticesOnEdge", &verticesOnEdgeVarId);
  ERRORS(success, "Failed to get variable id of verticesOnEdge.");
  NCDF_SIZE tmp_starts[2] = {0, 0};
  NCDF_SIZE tmp_counts[2] = {static_cast<size_t>(nEdges), 2};
  verticesOnEdge.resize(nEdges * 2);
  success = NCFUNCAG(_vara_int)(_fileId, verticesOnEdgeVarId, tmp_starts, tmp_counts, &verticesOnEdge[0] NCREQ);
  ERRORS(success, "Failed to read variable values of verticesOnEdge.");

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

        // Restore cellHandleToGlobalID
        Range::const_iterator rit;
        int i;
        for (rit = local_cells.begin(), i = 0; rit != local_cells.end(); ++rit, i++)
          cellHandleToGlobalID[*rit] = gids[i];
      }
    }

    // Restore numCellGroups
    if (0 == numCellGroups) {
      Tag numCellGroupsTag;
      rval = mbImpl->tag_get_handle("__NUM_CELL_GROUPS", 1, MB_TYPE_INTEGER, numCellGroupsTag);
      rval = mbImpl->tag_get_data(numCellGroupsTag, &_fileSet, 1, &numCellGroups);
    }
  }

  return MB_SUCCESS;
}

ErrorCode NCHelperMPAS::create_mesh(Range& faces)
{
  Interface*& mbImpl = _readNC->mbImpl;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  const Tag*& mpFileIdTag = _readNC->mpFileIdTag;

  int rank = 0, procs = 1;
#ifdef USE_MPI
  bool& isParallel = _readNC->isParallel;
  if (isParallel) {
    ParallelComm*& myPcomm = _readNC->myPcomm;
    rank = myPcomm->proc_config().proc_rank();
    procs = myPcomm->proc_config().proc_size();
  }
#endif

  ErrorCode rval;

  // Compute the number of local cells on this proc
  nLocalCells = int(std::floor(1.0 * nCells / procs));
  // start_cell_idx is the starting cell index in the MPAS file for this proc
  int start_cell_idx = rank * nLocalCells;
  // iextra = # cells extra after equal split over procs
  int iextra = nCells % procs;
  if (rank < iextra)
    nLocalCells++;
  start_cell_idx += std::min(rank, iextra);
  start_cell_idx++; // 0 based -> 1 based

  localGidCells.insert(start_cell_idx, start_cell_idx + nLocalCells - 1);

  // Read number of edges on each cell
  int nEdgesOnCellVarId;
  int success = NCFUNC(inq_varid)(_fileId, "nEdgesOnCell", &nEdgesOnCellVarId);
  ERRORS(success, "Failed to get variable id of nEdgesOnCell.");
  NCDF_SIZE tmp_starts_1[1] = {static_cast<size_t>(start_cell_idx - 1)};
  NCDF_SIZE tmp_counts_1[1] = {static_cast<size_t>(nLocalCells)};
  std::vector<int> num_edges_on_cell(nLocalCells);
  success = NCFUNCAG(_vara_int)(_fileId, nEdgesOnCellVarId, tmp_starts_1, tmp_counts_1, &num_edges_on_cell[0] NCREQ);
  ERRORS(success, "Failed to read variable values of nEdgesOnCell.");

  // Read vertices on each cell (connectivity)
  int verticesOnCellVarId;
  success = NCFUNC(inq_varid)(_fileId, "verticesOnCell", &verticesOnCellVarId);
  ERRORS(success, "Failed to get variable id of verticesOnCell.");
  NCDF_SIZE tmp_starts_2[2] = {static_cast<size_t>(start_cell_idx - 1), 0};
  NCDF_SIZE tmp_counts_2[2] = {static_cast<size_t>(nLocalCells), maxCellEdges};
  std::vector<int> vertices_on_cell(nLocalCells * maxCellEdges);
  success = NCFUNCAG(_vara_int)(_fileId, verticesOnCellVarId, tmp_starts_2, tmp_counts_2, &vertices_on_cell[0] NCREQ);
  ERRORS(success, "Failed to read variable values of verticesOnCell.");

  // Read edges on each cell
  int edgesOnCellVarId;
  success = NCFUNC(inq_varid)(_fileId, "edgesOnCell", &edgesOnCellVarId);
  ERRORS(success, "Failed to get variable id of edgesOnCell.");
  NCDF_SIZE tmp_starts_3[2] = {static_cast<size_t>(start_cell_idx - 1), 0};
  NCDF_SIZE tmp_counts_3[2] = {static_cast<size_t>(nLocalCells), maxCellEdges};
  std::vector<int> edges_on_cell(nLocalCells * maxCellEdges);
  success = NCFUNCAG(_vara_int)(_fileId, edgesOnCellVarId, tmp_starts_3, tmp_counts_3, &edges_on_cell[0] NCREQ);
  ERRORS(success, "Failed to read variable values of edgesOnCell.");

  // Divide cells into groups based on the number of edges
  std::vector<int> cells_with_n_edges[MAX_EDGES_PER_CELL + 1];
  for (int i = 0; i < nLocalCells; i++) {
    int cell_index = start_cell_idx + i; // Global cell index
    int num_edges = num_edges_on_cell[i];
    cells_with_n_edges[num_edges].push_back(cell_index);
  }

  // For each non-empty cell group, create cells and set connectivity array initially based on file ids
  EntityHandle start_element;
  EntityHandle* conn_arr_cells_with_n_edges[MAX_EDGES_PER_CELL + 1];
  Range tmp_range;
  void* data;
  int count;
  int* gid_data;
  std::set<int> local_vertices_set;
  std::set<int> local_edges_set;
  numCellGroups = 0;
  for (int i = 3; i <= maxCellEdges; i++) {
    int num_cells = cells_with_n_edges[i].size();
    if (num_cells > 0) {
      numCellGroups++;

      rval = _readNC->readMeshIface->get_element_connect(num_cells, i, MBPOLYGON, 0, start_element, conn_arr_cells_with_n_edges[i], num_cells);
      tmp_range.insert(start_element, start_element + num_cells - 1);
      faces.insert(start_element, start_element + num_cells - 1);

      // Get ptr to gid memory for cells
      Range cell_range(start_element, start_element + num_cells - 1);
      rval = mbImpl->tag_iterate(mGlobalIdTag, cell_range.begin(), cell_range.end(), count, data);
      ERRORR(rval, "Failed to get tag iterator on global id tag.");
      assert(count == (int) num_cells);
      gid_data = (int*) data;
      std::copy(cells_with_n_edges[i].begin(), cells_with_n_edges[i].end(), gid_data);

      for (int j = 0; j < num_cells; j++) {
        int cell_idx = cells_with_n_edges[i][j]; // Global cell index
        cellHandleToGlobalID[start_element + j] = cell_idx;
        cell_idx -= start_cell_idx; // Local cell index
        for (int k = 0; k < i; k++) {
          conn_arr_cells_with_n_edges[i][i * j + k] = vertices_on_cell[cell_idx * maxCellEdges + k];
          local_vertices_set.insert(vertices_on_cell[cell_idx * maxCellEdges + k]);
          local_edges_set.insert(edges_on_cell[cell_idx * maxCellEdges + k]);
        }
      }
    }
  }

  // Set tag for numCellGroups
  Tag numCellGroupsTag = 0;
  rval = mbImpl->tag_get_handle("__NUM_CELL_GROUPS", 1, MB_TYPE_INTEGER, numCellGroupsTag, MB_TAG_SPARSE | MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_CELL_GROUPS tag.");
  rval = mbImpl->tag_set_data(numCellGroupsTag, &_fileSet, 1, &numCellGroups);
  ERRORR(rval, "Trouble setting data for __NUM_CELL_GROUPS tag.");

  // Collect localGid for vertices
  std::copy(local_vertices_set.rbegin(), local_vertices_set.rend(), range_inserter(localGidVerts));
  nLocalVertices = localGidVerts.size();

  // Create local vertices
  EntityHandle start_vertex;
  std::vector<double*> arrays;
  rval = _readNC->readMeshIface->get_node_coords(3, nLocalVertices, 0, start_vertex, arrays);
  tmp_range.insert(start_vertex, start_vertex + nLocalVertices - 1);

  // Set coordinates for local vertices
  double* xptr = arrays[0];
  double* yptr = arrays[1];
  double* zptr = arrays[2];
  Range::const_iterator rit;
  int vert_idx;
  for (vert_idx = 0, rit = localGidVerts.begin(); vert_idx < nLocalVertices; vert_idx++, ++rit) {
    assert(*rit < xVertVals.size() + 1);
    xptr[vert_idx] = xVertVals[(*rit) - 1];
    yptr[vert_idx] = yVertVals[(*rit) - 1];
    zptr[vert_idx] = zVertVals[(*rit) - 1];
  }

  // Get ptr to gid memory for vertices
  Range vert_range(start_vertex, start_vertex + nLocalVertices - 1);
  rval = mbImpl->tag_iterate(mGlobalIdTag, vert_range.begin(), vert_range.end(), count, data);
  ERRORR(rval, "Failed to get tag iterator on global id tag.");
  assert(count == (int) nLocalVertices);
  gid_data = (int*) data;
  std::copy(localGidVerts.begin(), localGidVerts.end(), gid_data);

  // Duplicate global id data, which will be used to resolve sharing
  if (mpFileIdTag) {
    rval = mbImpl->tag_iterate(*mpFileIdTag, vert_range.begin(), vert_range.end(), count, data);
    ERRORR(rval, "Failed to get tag iterator on file id tag.");
    assert(count == (int) nLocalVertices);
    gid_data = (int*) data;
    std::copy(localGidVerts.begin(), localGidVerts.end(), gid_data);
  }

  // Create map from file ids to vertex handles
  std::map<EntityHandle, EntityHandle> vert_handles;
  for (rit = localGidVerts.begin(), vert_idx = 0; rit != localGidVerts.end(); ++rit, vert_idx++)
    vert_handles[*rit] = start_vertex + vert_idx;

  // For each non-empty cell group, set connectivity array with proper local vertices handles
  for (int i = 3; i <= maxCellEdges; i++) {
    int num_cells = cells_with_n_edges[i].size();
    if (num_cells > 0) {
      for (int j = 0; j < num_cells; j++) {
        for (int k = 0; k < i; k++) {
          EntityHandle global_vert_id = conn_arr_cells_with_n_edges[i][i * j + k];
          conn_arr_cells_with_n_edges[i][i * j + k] = vert_handles[global_vert_id];
        }
      }
    }
  }

  // Collect localGid for edges
  std::copy(local_edges_set.rbegin(), local_edges_set.rend(), range_inserter(localGidEdges));
  nLocalEdges = localGidEdges.size();

  // Create local edges
  EntityHandle start_edge;
  EntityHandle* conn_arr_edges;
  rval = _readNC->readMeshIface->get_element_connect(nLocalEdges, 2, MBEDGE, 0, start_edge, conn_arr_edges);
  tmp_range.insert(start_edge, start_edge + nLocalEdges - 1);

  // Set vertices for local edges
  int edge_idx;
  for (rit = localGidEdges.begin(), edge_idx = 0; rit != localGidEdges.end(); ++rit, edge_idx += 2) {
    EntityHandle gloabl_edge_id = *rit;
    EntityHandle gloabl_vert_id_1 = verticesOnEdge[(gloabl_edge_id - 1) * 2 + 0];
    EntityHandle gloabl_vert_id_2 = verticesOnEdge[(gloabl_edge_id - 1) * 2 + 1];
    conn_arr_edges[edge_idx] = vert_handles[gloabl_vert_id_1];
    conn_arr_edges[edge_idx + 1] = vert_handles[gloabl_vert_id_2];
  }

  // Get ptr to gid memory for edges
  Range edge_range(start_edge, start_edge + nLocalEdges - 1);
  rval = mbImpl->tag_iterate(mGlobalIdTag, edge_range.begin(), edge_range.end(), count, data);
  ERRORR(rval, "Failed to get tag iterator on global id tag.");
  assert(count == (int) nLocalEdges);
  gid_data = (int*) data;
  std::copy(localGidEdges.begin(), localGidEdges.end(), gid_data);

  // Add new vertices, elements and edges to the file set
  rval = _readNC->mbImpl->add_entities(_fileSet, tmp_range);
  ERRORR(rval, "Couldn't add new vertices/faces/edges to file set.");

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
      if (3 == vd.varDims.size())
      {
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
  std::vector<int>& dimVals = _readNC->dimVals;
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
    vdatas[i].numLev = nLevels;

    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      // Get the tag to read into
      if (!vdatas[i].varTags[t]) {
        rval = get_tag_to_nonset(vdatas[i], tstep_nums[t], vdatas[i].varTags[t], vdatas[i].numLev);
        ERRORR(rval, "Trouble getting tag.");
      }

      // Assume point-based values for now?
      if (-1 == tDim || dimVals[tDim] <= (int) t) {
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
      if (vdatas[i].entLoc == ReadNC::ENTLOCFACE && numCellGroups > 1)
        vdatas[i].varDatas[t] = NULL;
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
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_ucd_variable_to_nonset_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables.");

  // Finally, read into that space
  int success;
  Range* pLocalGid = NULL;

  for (unsigned int i = 0; i < vdatas.size(); i++) {
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
            success = NCFUNCAG2(_vara_double)(_fileId, vdatas[i].varId,
                &(vdatas[i].readStarts[t][0]), &(vdatas[i].readCounts[t][0]),
                            &(tmpdoubledata[indexInDoubleArray]) NCREQ2);
            ERRORS(success, "Failed to read double data in loop");
            // We need to increment the index in double array for the
            // next subrange
            indexInDoubleArray += (endh - starth + 1) * 1 * vdatas[i].numLev;
          }
          assert(ic == pLocalGid->psize());

          success = ncmpi_wait_all(_fileId, requests.size(), &requests[0], &statuss[0]);
          ERRORS(success, "Failed on wait_all.");

          // Local cells are grouped by the number of edges on each cell, e.g. 5, 6 or 7
          // Tags created for cell variables may have to read data from different groups
          if (vdatas[i].entLoc == ReadNC::ENTLOCFACE && numCellGroups > 1) {
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
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_ucd_variable_to_nonset_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables.");

  // Finally, read into that space
  int success;
  Range* pLocalGid = NULL;

  std::vector<int> requests(vdatas.size() * tstep_nums.size()), statuss(vdatas.size() * tstep_nums.size());
  for (unsigned int i = 0; i < vdatas.size(); i++) {
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
                            &(tmpdoubledata[indexInDoubleArray]) NCREQ);
            ERRORS(success, "Failed to read double data in loop");
            // We need to increment the index in double array for the
            // next subrange
            indexInDoubleArray += (endh - starth + 1) * 1 * vdatas[i].numLev;
          }
          assert(ic == pLocalGid->psize());

          // Local cells are grouped by the number of edges on each cell, e.g. 5, 6 or 7
          // Tags created for cell variables may have to read data from different groups
          if (vdatas[i].entLoc == ReadNC::ENTLOCFACE && numCellGroups > 1) {
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

#ifdef NCWAIT
  int success = ncmpi_wait_all(fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

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
#endif

} // namespace moab
