#include "NCHelperHOMME.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/FileOptions.hpp"
#include "moab/SpectralMeshTool.hpp"

#include <cmath>

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) {_readNC->readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
  if (err) {_readNC->readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

NCHelperHOMME::NCHelperHOMME(ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
: UcdNCHelper(readNC, fileId, opts, fileSet),
_spectralOrder(-1), connectId(-1)
{
  // Calculate spectral order
  std::map<std::string, ReadNC::AttData>::iterator attIt = readNC->globalAtts.find("np");
  if (attIt != readNC->globalAtts.end()) {
    int success = NCFUNC(get_att_int)(readNC->fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &_spectralOrder);
    if (success != 0)
      readNC->readMeshIface->report_error("%s", "Failed to read np global attribute int data.");
    else
      _spectralOrder--; // Spectral order is one less than np
  }
}

bool NCHelperHOMME::can_read_file(ReadNC* readNC, int fileId)
{
  // If global attribute "np" exists then it should be the HOMME grid
  if (readNC->globalAtts.find("np") != readNC->globalAtts.end()) {
    // Make sure it is CAM grid
    std::map<std::string, ReadNC::AttData>::iterator attIt = readNC->globalAtts.find("source");
    if (attIt == readNC->globalAtts.end()) {
      readNC->readMeshIface->report_error("%s", "File does not have source global attribute.");
      return false;
    }
    unsigned int sz = attIt->second.attLen;
    std::string att_data;
    att_data.resize(sz + 1);
    att_data[sz] = '\000';
    int success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0]);
    if (success != 0) {
      readNC->readMeshIface->report_error("%s", "Failed to read source global attribute char data.");
      return false;
    }
    if (att_data.find("CAM") == std::string::npos)
      return false;

    return true;
  }

  return false;
}

ErrorCode NCHelperHOMME::init_mesh_vals()
{
  std::vector<std::string>& dimNames = _readNC->dimNames;
  std::vector<int>& dimLens = _readNC->dimLens;
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;

  ErrorCode rval;
  unsigned int idx;
  std::vector<std::string>::iterator vit;

  // Look for time dimension
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end())
    idx = vit - dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "t")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'time' or 't' dimension.");
  }
  tDim = idx;
  nTimeSteps = dimLens[idx];

  // Get number of vertices (labeled as number of columns)
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "ncol")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'ncol' dimension.");
  }
  vDim = idx;
  nVertices = dimLens[idx];

  // Set number of cells
  nCells = nVertices - 2;

  // Get number of levels
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lev")) != dimNames.end())
    idx = vit - dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "ilev")) != dimNames.end())
    idx = vit - dimNames.begin();
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'lev' or 'ilev' dimension.");
  }
  levDim = idx;
  nLevels = dimLens[idx];

  // Store lon values in xVertVals
  std::map<std::string, ReadNC::VarData>::iterator vmit;
  if ((vmit = varInfo.find("lon")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate("lon", 0, nVertices - 1, xVertVals);
    ERRORR(rval, "Trouble reading 'lon' variable.");
  }
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'lon' variable.");
  }

  // Store lat values in yVertVals
  if ((vmit = varInfo.find("lat")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate("lat", 0, nVertices - 1, yVertVals);
    ERRORR(rval, "Trouble reading 'lat' variable.");
  }
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'lat' variable.");
  }

  // Store lev values in levVals
  if ((vmit = varInfo.find("lev")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate("lev", 0, nLevels - 1, levVals);
    ERRORR(rval, "Trouble reading 'lev' variable.");

    // Decide whether down is positive
    char posval[10];
    int success = NCFUNC(get_att_text)(_fileId, (*vmit).second.varId, "positive", posval);
    if (0 == success && !strcmp(posval, "down")) {
      for (std::vector<double>::iterator dvit = levVals.begin(); dvit != levVals.end(); ++dvit)
        (*dvit) *= -1.0;
    }
  }
  else {
    ERRORR(MB_FAILURE, "Couldn't find 'lev' variable.");
  }

  // Store time coordinate values in tVals
  if ((vmit = varInfo.find("time")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate("time", 0, nTimeSteps - 1, tVals);
    ERRORR(rval, "Trouble reading 'time' variable.");
  }
  else if ((vmit = varInfo.find("t")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate("t", 0, nTimeSteps - 1, tVals);
    ERRORR(rval, "Trouble reading 't' variable.");
  }
  else {
    // If expected time variable is not available, set dummy time coordinate values to tVals
    for (int t = 0; t < nTimeSteps; t++)
      tVals.push_back((double)t);
  }

  // For each variable, determine the entity location type and number of levels
  std::map<std::string, ReadNC::VarData>::iterator mit;
  for (mit = varInfo.begin(); mit != varInfo.end(); ++mit) {
    ReadNC::VarData& vd = (*mit).second;

    vd.entLoc = ReadNC::ENTLOCSET;
    if (std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) {
      if (std::find(vd.varDims.begin(), vd.varDims.end(), vDim) != vd.varDims.end())
        vd.entLoc = ReadNC::ENTLOCVERT;
    }

    vd.numLev = 1;
    if (std::find(vd.varDims.begin(), vd.varDims.end(), levDim) != vd.varDims.end())
      vd.numLev = nLevels;
  }

  // Hack: create dummy tags for dimensions (like ncol) with no corresponding coordinate variables
  init_dims_with_no_coord_vars_info();

  return MB_SUCCESS;
}

// When noMesh option is used on this read, the old ReadNC class instance for last read can get out
// of scope (and deleted). The old instance initialized localGidVerts properly when the mesh was
// created, but it is now lost. The new instance (will not create the mesh with noMesh option) has
// to restore it based on the existing mesh from last read
ErrorCode NCHelperHOMME::check_existing_mesh()
{
  Interface*& mbImpl = _readNC->mbImpl;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  bool& noMesh = _readNC->noMesh;

  if (noMesh && localGidVerts.empty()) {
    // We need to populate localGidVerts range with the gids of vertices from the tmp_set
    // localGidVerts is important in reading the variable data into the nodes
    // also, for our purposes, localGidVerts is truly the GLOBAL_ID tag data, not other
    // file_id tags that could get passed around in other scenarios for parallel reading
    // for nodal_partition, this local gid is easier, should be initialized with only
    // the owned nodes

    // We need to get all vertices from tmp_set (it is the input set in no_mesh scenario)
    Range local_verts;
    ErrorCode rval = mbImpl->get_entities_by_dimension(_fileSet, 0, local_verts);
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

  return MB_SUCCESS;
}

ErrorCode NCHelperHOMME::create_mesh(Range& faces)
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::string& fileName = _readNC->fileName;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
  DebugOutput& dbgOut = _readNC->dbgOut;
  bool& spectralMesh = _readNC->spectralMesh;
  int& gatherSetRank = _readNC->gatherSetRank;

  // Need to get/read connectivity data before creating elements
  std::string conn_fname;

  // Try to open the connectivity file through CONN option, if used
  ErrorCode rval = _opts.get_str_option("CONN", conn_fname);
  if (MB_SUCCESS != rval) {
    // Default convention for reading HOMME is a file HommeMapping.nc in same dir as data file
    conn_fname = std::string(fileName);
    size_t idx = conn_fname.find_last_of("/");
    if (idx != std::string::npos)
      conn_fname = conn_fname.substr(0, idx).append("/HommeMapping.nc");
    else
      conn_fname = "HommeMapping.nc";
  }

  int success = 0;

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

#ifdef PNETCDF_FILE
#ifdef USE_MPI
  if (isParallel) {
    ParallelComm*& myPcomm = _readNC->myPcomm;
    success = NCFUNC(open)(myPcomm->proc_config().proc_comm(), conn_fname.c_str(), 0, MPI_INFO_NULL, &connectId);
  }
  else
    success = NCFUNC(open)(MPI_COMM_SELF, conn_fname.c_str(), 0, MPI_INFO_NULL, &connectId);
#endif
#else
  success = NCFUNC(open)(conn_fname.c_str(), 0, &connectId);
#endif
  ERRORS(success, "Failed on open.");

  std::vector<std::string> conn_names;
  std::vector<int> conn_vals;
  rval = _readNC->get_dimensions(connectId, conn_names, conn_vals);
  ERRORR(rval, "Failed to get dimensions for connectivity.");

  // Read connectivity into temporary variable
  int num_fine_quads = 0;
  int num_coarse_quads = 0;
  int start_idx = 0;
  std::vector<std::string>::iterator vit;
  int idx = 0;
  if ((vit = std::find(conn_names.begin(), conn_names.end(), "ncells")) != conn_names.end())
    idx = vit - conn_names.begin();
  else if ((vit = std::find(conn_names.begin(), conn_names.end(), "ncenters")) != conn_names.end())
    idx = vit - conn_names.begin();
  else {
    ERRORR(MB_FAILURE, "Failed to get number of quads.");
  }
  int num_quads = conn_vals[idx];
  if (num_quads != nCells) {
    dbgOut.tprintf(1, "Warning: number of quads from %s and cells from %s are inconsistent; num_quads = %d, nCells = %d.\n",
        conn_fname.c_str(), fileName.c_str(), num_quads, nCells);
  }

  // Get the connectivity into tmp_conn2 and permute into tmp_conn
  int cornerVarId;
  success = NCFUNC(inq_varid)(connectId, "element_corners", &cornerVarId);
  ERRORS(success, "Failed to get variable id.");
  NCDF_SIZE tmp_starts[2] = {0, 0};
  NCDF_SIZE tmp_counts[2] = {4, static_cast<NCDF_SIZE>(num_quads)};
  std::vector<int> tmp_conn(4 * num_quads), tmp_conn2(4 * num_quads);
  success = NCFUNCAG(_vara_int)(connectId, cornerVarId, tmp_starts, tmp_counts, &tmp_conn2[0]);
  ERRORS(success, "Failed to get temporary connectivity.");
  success = NCFUNC(close)(connectId);
  ERRORS(success, "Failed on close.");
  // Permute the connectivity
  for (int i = 0; i < num_quads; i++) {
    tmp_conn[4 * i] = tmp_conn2[i];
    tmp_conn[4 * i + 1] = tmp_conn2[i + 1 * num_quads];
    tmp_conn[4 * i + 2] = tmp_conn2[i + 2 * num_quads];
    tmp_conn[4 * i + 3] = tmp_conn2[i + 3 * num_quads];
  }

  // Need to know whether we'll be creating gather mesh later, to make sure we allocate enough space
  // in one shot
  bool create_gathers = false;
  if (rank == gatherSetRank)
    create_gathers = true;

  // Compute the number of local quads, accounting for coarse or fine representation
  // spectral_unit is the # fine quads per coarse quad, or spectralOrder^2
  int spectral_unit = (spectralMesh ? _spectralOrder * _spectralOrder : 1);
  // num_coarse_quads is the number of quads instantiated in MOAB; if !spectralMesh, num_coarse_quads = num_fine_quads
  num_coarse_quads = int(std::floor(1.0 * num_quads / (spectral_unit * procs)));
  // start_idx is the starting index in the HommeMapping connectivity list for this proc, before converting to coarse quad representation
  start_idx = 4 * rank * num_coarse_quads * spectral_unit;
  // iextra = # coarse quads extra after equal split over procs
  int iextra = num_quads % (procs * spectral_unit);
  if (rank < iextra)
    num_coarse_quads++;
  start_idx += 4 * spectral_unit * std::min(rank, iextra);
  // num_fine_quads is the number of quads in the connectivity list in HommeMapping file assigned to this proc
  num_fine_quads = spectral_unit * num_coarse_quads;

  // Now create num_coarse_quads
  EntityHandle* conn_arr;
  EntityHandle start_vertex;
  Range tmp_range;

  // Read connectivity into that space
  EntityHandle* sv_ptr = NULL;
  EntityHandle start_quad;
  SpectralMeshTool smt(mbImpl, _spectralOrder);
  if (!spectralMesh) {
    rval = _readNC->readMeshIface->get_element_connect(num_coarse_quads, 4,
                                              MBQUAD, 0, start_quad, conn_arr,
                                              // Might have to create gather mesh later
                                              (create_gathers ? num_coarse_quads + num_quads : num_coarse_quads));
    ERRORR(rval, "Failed to create quads.");
    tmp_range.insert(start_quad, start_quad + num_coarse_quads - 1);
    std::copy(&tmp_conn[start_idx], &tmp_conn[start_idx + 4 * num_fine_quads], conn_arr);
    std::copy(conn_arr, conn_arr + 4 * num_fine_quads, range_inserter(localGidVerts));
  }
  else {
    rval = smt.create_spectral_elems(&tmp_conn[0], num_fine_quads, 2, tmp_range, start_idx, &localGidVerts);
    ERRORR(rval, "Failed to create spectral elements.");
    int count, v_per_e;
    rval = mbImpl->connect_iterate(tmp_range.begin(), tmp_range.end(), conn_arr, v_per_e, count);
    ERRORR(rval, "Failed to get connectivity of spectral elements.");
    rval = mbImpl->tag_iterate(smt.spectral_vertices_tag(true), tmp_range.begin(), tmp_range.end(),
                               count, (void*&)sv_ptr);
    ERRORR(rval, "Failed to get fine connectivity of spectral elements.");
  }

  // Create vertices
  nLocalVertices = localGidVerts.size();
  std::vector<double*> arrays;
  rval = _readNC->readMeshIface->get_node_coords(3, nLocalVertices, 0, start_vertex, arrays,
                                        // Might have to create gather mesh later
                                        (create_gathers ? nLocalVertices + nVertices : nLocalVertices));
  ERRORR(rval, "Couldn't create vertices in HOMME mesh.");

  // Set vertex coordinates
  Range::iterator rit;
  double* xptr = arrays[0];
  double* yptr = arrays[1];
  double* zptr = arrays[2];
  int i;
  for (i = 0, rit = localGidVerts.begin(); i < nLocalVertices; i++, ++rit) {
    assert(*rit < xVertVals.size() + 1);
    xptr[i] = xVertVals[(*rit) - 1]; // lon
    yptr[i] = yVertVals[(*rit) - 1]; // lat
  }

  // Convert lon/lat/rad to x/y/z
  const double pideg = acos(-1.0) / 180.0;
  for (i = 0; i < nLocalVertices; i++) {
    double cosphi = cos(pideg * yptr[i]);
    double zmult = sin(pideg * yptr[i]);
    double xmult = cosphi * cos(xptr[i] * pideg);
    double ymult = cosphi * sin(xptr[i] * pideg);
    double rad = 8000.0 + levVals[0];
    xptr[i] = rad * xmult;
    yptr[i] = rad * ymult;
    zptr[i] = rad * zmult;
  }

  // Get ptr to gid memory for vertices
  Range vert_range(start_vertex, start_vertex + nLocalVertices - 1);
  void* data;
  int count;
  rval = mbImpl->tag_iterate(mGlobalIdTag, vert_range.begin(), vert_range.end(), count, data);
  ERRORR(rval, "Failed to get tag iterator.");
  assert(count == nLocalVertices);
  int* gid_data = (int*) data;
  std::copy(localGidVerts.begin(), localGidVerts.end(), gid_data);

  // Duplicate global id data, which will be used to resolve sharing
  if (mpFileIdTag) {
    rval = mbImpl->tag_iterate(*mpFileIdTag, vert_range.begin(), vert_range.end(), count, data);
    ERRORR(rval, "Failed to get tag iterator on file id tag.");
    assert(count == nLocalVertices);
    gid_data = (int*) data;
    std::copy(localGidVerts.begin(), localGidVerts.end(), gid_data);
  }

  // Create map from file ids to vertex handles, used later to set connectivity
  std::map<EntityHandle, EntityHandle> vert_handles;
  for (rit = localGidVerts.begin(), i = 0; rit != localGidVerts.end(); ++rit, i++)
    vert_handles[*rit] = start_vertex + i;

  // Compute proper handles in connectivity using offset
  for (int q = 0; q < 4 * num_coarse_quads; q++) {
    conn_arr[q] = vert_handles[conn_arr[q]];
    assert(conn_arr[q]);
  }
  if (spectralMesh) {
    int verts_per_quad = (_spectralOrder + 1) * (_spectralOrder + 1);
    for (int q = 0; q < verts_per_quad * num_coarse_quads; q++) {
      sv_ptr[q] = vert_handles[sv_ptr[q]];
      assert(sv_ptr[q]);
    }
  }

  // Add new vertices and elements to the set
  faces.merge(tmp_range);
  tmp_range.insert(start_vertex, start_vertex + nLocalVertices - 1);
  rval = mbImpl->add_entities(_fileSet, tmp_range);
  ERRORR(rval, "Couldn't add new vertices and quads to file set.");

  // Mark the set with the spectral order
  Tag sporder;
  rval = mbImpl->tag_get_handle("SPECTRAL_ORDER", 1, MB_TYPE_INTEGER, sporder, MB_TAG_CREAT | MB_TAG_SPARSE);
  ERRORR(rval, "Couldn't create spectral order tag.");
  rval = mbImpl->tag_set_data(sporder, &_fileSet, 1, &_spectralOrder);
  ERRORR(rval, "Couldn't set value for spectral order tag.");

  if (create_gathers) {
    EntityHandle gather_set;
    rval = _readNC->readMeshIface->create_gather_set(gather_set);
    ERRORR(rval, "Trouble creating gather set.");

    // Create vertices
    arrays.clear();
    // Don't need to specify allocation number here, because we know enough verts were created before
    rval = _readNC->readMeshIface->get_node_coords(3, nVertices, 0, start_vertex, arrays);
    ERRORR(rval, "Couldn't create vertices in HOMME mesh for gather set.");

    xptr = arrays[0];
    yptr = arrays[1];
    zptr = arrays[2];
    for (i = 0; i < nVertices; i++) {
      double cosphi = cos(pideg * yVertVals[i]);
      double zmult = sin(pideg * yVertVals[i]);
      double xmult = cosphi * cos(xVertVals[i] * pideg);
      double ymult = cosphi * sin(xVertVals[i] * pideg);
      double rad = 8000.0 + levVals[0];
      xptr[i] = rad * xmult;
      yptr[i] = rad * ymult;
      zptr[i] = rad * zmult;
    }

    // Get ptr to gid memory for vertices
    Range gather_verts(start_vertex, start_vertex + nVertices - 1);
    rval = mbImpl->tag_iterate(mGlobalIdTag, gather_verts.begin(), gather_verts.end(), count, data);
    ERRORR(rval, "Failed to get tag iterator.");
    assert(count == nVertices);
    gid_data = (int*) data;
    for (int j = 1; j <= nVertices; j++)
      gid_data[j - 1] = j;
    // Set the file id tag too, it should be bigger something not interfering with global id
    if (mpFileIdTag) {
      rval = mbImpl->tag_iterate(*mpFileIdTag, gather_verts.begin(), gather_verts.end(), count, data);
      ERRORR(rval, "Failed to get tag iterator in file id tag.");
      assert(count == nVertices);
      gid_data = (int*) data;
      for (int j = 1; j <= nVertices; j++)
        gid_data[j - 1] = nVertices + j; // bigger than global id tag
    }

    rval = mbImpl->add_entities(gather_set, gather_verts);
    ERRORR(rval, "Couldn't add vertices to gather set.");

    // Create quads
    Range gather_quads;
    // Don't need to specify allocation number here, because we know enough quads were created before
    rval = _readNC->readMeshIface->get_element_connect(num_quads, 4,
                                              MBQUAD, 0, start_quad, conn_arr);
    ERRORR(rval, "Failed to create quads.");
    gather_quads.insert(start_quad, start_quad + num_quads - 1);
    std::copy(&tmp_conn[0], &tmp_conn[4 * num_quads], conn_arr);
    for (i = 0; i != 4 * num_quads; i++)
      conn_arr[i] += start_vertex - 1; // Connectivity array is shifted by where the gather verts start
    rval = mbImpl->add_entities(gather_set, gather_quads);
    ERRORR(rval, "Couldn't add quads to gather set.");
  }

  return MB_SUCCESS;
}

ErrorCode NCHelperHOMME::read_ucd_variable_to_nonset_allocate(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
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

      // Finally: nVertices
      switch (vdatas[i].entLoc) {
        case ReadNC::ENTLOCVERT:
          // Vertices
          // We will start from the first localGidVerts, actually; we will reset that
          // later on, anyway, in a loop
          vdatas[i].readStarts[t].push_back(localGidVerts[0] - 1);
          vdatas[i].readCounts[t].push_back(nLocalVertices);
          assert(vdatas[i].readStarts[t].size() == vdatas[i].varDims.size());
          range = &verts;
          break;
        default:
          ERRORR(MB_FAILURE, "Unexpected entity location type for HOMME non-set variable.");
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

#ifdef PNETCDF_FILE
ErrorCode NCHelperHOMME::read_ucd_variable_to_nonset_async(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_ucd_variable_to_nonset_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables.");

  // Finally, read into that space
  int success;

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    std::size_t sz = vdatas[i].sz;

    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      // We will synchronize all these reads with the other processors,
      // so the wait will be inside this double loop; is it too much?
      size_t nb_reads = localGidVerts.psize();
      std::vector<int> requests(nb_reads), statuss(nb_reads);
      size_t idxReq = 0;
      void* data = vdatas[i].varDatas[t];
      size_t ni = vdatas[i].readCounts[t][2];
      size_t nj = 1; // For HOMME, nj holds # quads, so here should set to 1
      size_t nk = vdatas[i].readCounts[t][1];

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
          // localGidVerts range;
          // basically, we have to give a different point
          // for data to start, for every subrange :(
          size_t nbDims = vdatas[i].readStarts[t].size();
          // assume that the last dimension is for the ncol,
          // node varying variable

          size_t indexInDoubleArray = 0;
          size_t ic = 0;
          for (Range::pair_iterator pair_iter = localGidVerts.pair_begin();
              pair_iter != localGidVerts.pair_end();
              pair_iter++, ic++) {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second; // inclusive
            vdatas[i].readStarts[t][nbDims - 1] = (NCDF_SIZE) (starth - 1);
            vdatas[i].readCounts[t][nbDims - 1] = (NCDF_SIZE) (endh - starth + 1);

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
          assert(ic == localGidVerts.psize());

          success = ncmpi_wait_all(_fileId, requests.size(), &requests[0], &statuss[0]);
          ERRORS(success, "Failed on wait_all.");

          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik_stride(ni, nj, nk, data, &tmpdoubledata[0], localGidVerts);
          else {
            for (std::size_t idx = 0; idx != tmpdoubledata.size(); idx++)
              ((double*) data)[idx] = tmpdoubledata[idx];
          }
          ERRORS(success, "Failed to read double data.");
          break;
        }
        case NC_FLOAT: {
          std::vector<float> tmpfloatdata(sz);

          // In the case of ucd mesh, and on multiple proc,
          // we need to read as many times as subranges we have in the
          // localGidVerts range;
          // basically, we have to give a different point
          // for data to start, for every subrange :(
          size_t nbDims = vdatas[i].readStarts[t].size();

          // Assume that the last dimension is for the ncol, number of vertices
          size_t indexInFloatArray = 0;
          size_t ic = 0;
          for (Range::pair_iterator pair_iter = localGidVerts.pair_begin();
              pair_iter != localGidVerts.pair_end();
              pair_iter++, ic++) {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second; // inclusive
            vdatas[i].readStarts[t][nbDims - 1] = (NCDF_SIZE) (starth - 1);
            vdatas[i].readCounts[t][nbDims - 1] = (NCDF_SIZE) (endh - starth + 1);

            // Do a partial read, in each subrange
            // wait outside this loop
            success = NCFUNCREQG(_vara_float)(_fileId, vdatas[i].varId,
                &(vdatas[i].readStarts[t][0]), &(vdatas[i].readCounts[t][0]),
                            &(tmpfloatdata[indexInFloatArray]), &requests[idxReq++]);
            ERRORS(success, "Failed to read float data in loop");
            // We need to increment the index in float array for the
            // next subrange
            indexInFloatArray += (endh - starth + 1) * 1 * vdatas[i].numLev;
          }
          assert(ic == localGidVerts.psize());

          success = ncmpi_wait_all(_fileId, requests.size(), &requests[0], &statuss[0]);
          ERRORS(success, "Failed on wait_all.");

          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik_stride(ni, nj, nk, data, &tmpfloatdata[0], localGidVerts);
          else {
            for (std::size_t idx = 0; idx != tmpfloatdata.size(); idx++)
              ((float*) data)[idx] = tmpfloatdata[idx];
          }
          ERRORS(success, "Failed to read float data.");
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
ErrorCode NCHelperHOMME::read_ucd_variable_to_nonset(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums)
{
  DebugOutput& dbgOut = _readNC->dbgOut;

  ErrorCode rval = read_ucd_variable_to_nonset_allocate(vdatas, tstep_nums);
  ERRORR(rval, "Trouble allocating read variables.");

  // Finally, read into that space
  int success;
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    std::size_t sz = vdatas[i].sz;

    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      void* data = vdatas[i].varDatas[t];
      size_t ni = vdatas[i].readCounts[t][2];
      size_t nj = 1; // For HOMME, nj holds # quads, so here should set to 1
      size_t nk = vdatas[i].readCounts[t][1];

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
          // localGidVerts range;
          // basically, we have to give a different point
          // for data to start, for every subrange :(
          size_t nbDims = vdatas[i].readStarts[t].size();

          // Assume that the last dimension is for the ncol, number of vertices
          size_t indexInDoubleArray = 0;
          size_t ic = 0;
          for (Range::pair_iterator pair_iter = localGidVerts.pair_begin();
              pair_iter != localGidVerts.pair_end();
              pair_iter++, ic++) {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second; // Inclusive
            vdatas[i].readStarts[t][nbDims - 1] = (NCDF_SIZE) (starth - 1);
            vdatas[i].readCounts[t][nbDims - 1] = (NCDF_SIZE) (endh - starth + 1);

            success = NCFUNCAG(_vara_double)(_fileId, vdatas[i].varId,
                &(vdatas[i].readStarts[t][0]), &(vdatas[i].readCounts[t][0]),
                            &(tmpdoubledata[indexInDoubleArray]));
            ERRORS(success, "Failed to read float data in loop");
            // We need to increment the index in double array for the
            // next subrange
            indexInDoubleArray += (endh - starth + 1) * 1 * vdatas[i].numLev;
          }
          assert(ic == localGidVerts.psize());

          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik_stride(ni, nj, nk, data, &tmpdoubledata[0], localGidVerts);
          else {
            for (std::size_t idx = 0; idx != tmpdoubledata.size(); idx++)
              ((double*) data)[idx] = tmpdoubledata[idx];
          }
          ERRORS(success, "Failed to read double data.");
          break;
        }
        case NC_FLOAT: {
          std::vector<float> tmpfloatdata(sz);

          // In the case of ucd mesh, and on multiple proc,
          // we need to read as many times as subranges we have in the
          // localGidVerts range;
          // basically, we have to give a different point
          // for data to start, for every subrange :(
          size_t nbDims = vdatas[i].readStarts[t].size();

          // Assume that the last dimension is for the ncol
          size_t indexInFloatArray = 0;
          size_t ic = 0;
          for (Range::pair_iterator pair_iter = localGidVerts.pair_begin();
              pair_iter != localGidVerts.pair_end();
              pair_iter++, ic++) {
            EntityHandle starth = pair_iter->first;
            EntityHandle endh = pair_iter->second; // Inclusive
            vdatas[i].readStarts[t][nbDims-1] = (NCDF_SIZE) (starth - 1);
            vdatas[i].readCounts[t][nbDims-1] = (NCDF_SIZE) (endh - starth + 1);

            success = NCFUNCAG(_vara_float)(_fileId, vdatas[i].varId,
                &(vdatas[i].readStarts[t][0]), &(vdatas[i].readCounts[t][0]),
                            &(tmpfloatdata[indexInFloatArray]));
            ERRORS(success, "Failed to read float data in loop");
            // We need to increment the index in float array for the
            // next subrange
            indexInFloatArray += (endh - starth + 1) * 1 * vdatas[i].numLev;
          }
          assert(ic == localGidVerts.psize());

          if (vdatas[i].numLev != 1)
            // Switch from k varying slowest to k varying fastest
            success = kji_to_jik_stride(ni, nj, nk, data, &tmpfloatdata[0], localGidVerts);
          else {
            for (std::size_t idx = 0; idx != tmpfloatdata.size(); idx++)
              ((float*) data)[idx] = tmpfloatdata[idx];
          }
          ERRORS(success, "Failed to read float data.");
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
#endif

} // namespace moab
