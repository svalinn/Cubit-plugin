#include "NCHelperHOMME.hpp"
#include "moab/ReadUtilIface.hpp"
#include "FileOptions.hpp"
#include "moab/SpectralMeshTool.hpp"

#include <cmath>

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {_readNC->readMeshIface->report_error("%s", str); return rval;}

#define ERRORS(err, str) \
    if (err) {_readNC->readMeshIface->report_error("%s", str); return MB_FAILURE;}

namespace moab {

bool NCHelperHOMME::can_read_file(ReadNC* readNC, const FileOptions& opts, int& spectralOrder)
{
  // If global attribute "np" exists then it's the HOMME grid
  std::map<std::string, ReadNC::AttData>::iterator attIt = readNC->globalAtts.find("np");
  if (attIt != readNC->globalAtts.end())
  {
    int success = NCFUNC(get_att_int)(readNC->fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &spectralOrder);
    if (success != 0) {
      readNC->readMeshIface->report_error("%s", "Failed to read np global attribute int data.");
      return false;
    }

    spectralOrder--; // Spectral order is one less than np

    if (MB_SUCCESS == opts.match_option("PARTITION_METHOD", "NODAL_PARTITION"))
      readNC->partMethod = -1;

    return true;
  }

  return false;
}

ErrorCode NCHelperHOMME::init_mesh_vals(const FileOptions& opts, EntityHandle file_set)
{
  ErrorCode rval = init_HOMMEucd_vals();
  if (MB_SUCCESS != rval)
    _readNC->readMeshIface->report_error("%s", "Trouble initializing HOMME grid.");

  return rval;
}

ErrorCode NCHelperHOMME::create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads)
{
  ErrorCode rval = create_HOMMEucd_verts_quads(opts, file_set, quads);
  if (MB_SUCCESS != rval)
    _readNC->readMeshIface->report_error("%s", "Trouble creating vertices and quads for HOMME grid.");

  return rval;
}

ErrorCode NCHelperHOMME::init_HOMMEucd_vals()
{
  std::vector<std::string>& dimNames = _readNC->dimNames;
  std::vector<int>& dimVals = _readNC->dimVals;
  std::string& kName = _readNC->kName;
  std::string& tName = _readNC->tName;
  std::map<std::string, ReadNC::VarData>& varInfo = _readNC->varInfo;
  int& tMin = _readNC->tMin;
  int& tMax = _readNC->tMax;
  int (&gDims)[6] = _readNC->gDims;
  int (&lDims)[6] = _readNC->lDims;
  int& iDim = _readNC->iDim;
  int& kDim = _readNC->kDim;
  int& tDim = _readNC->tDim;
  std::vector<double>& ilVals = _readNC->ilVals;
  std::vector<double>& jlVals = _readNC->jlVals;
  std::vector<double>& klVals = _readNC->klVals;
  std::vector<double>& tVals = _readNC->tVals;

  ErrorCode rval;
  unsigned int idx;
  std::vector<std::string>::iterator vit;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end())
    idx = vit - dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "t")) != dimNames.end())
    idx = vit - dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find time variable.");
  tDim = idx;
  tMax = dimVals[idx] - 1;
  tMin = 0;
  tName = dimNames[idx];

  // get number of vertices (labeled as number of columns) and levels
  gDims[0] = gDims[3] = -1;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "ncol")) != dimNames.end()) {
    idx = vit - dimNames.begin();
    gDims[3] = dimVals[idx] - 1;
    gDims[0] = 0;
    iDim = idx;
  }
  if (-1 == gDims[0])
    return MB_FAILURE;

  // set j coordinate to the number of quads
  gDims[1] = gDims[0];
  gDims[4] = gDims[3] - 2;

  gDims[2] = gDims[5] = -1;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lev")) != dimNames.end()) {
    idx = vit - dimNames.begin();
    gDims[5] = dimVals[idx] - 1, gDims[2] = 0, kName = std::string("lev");
    kDim = idx;
  }
  if (-1 == gDims[2])
    return MB_FAILURE;

  // read coordinate data
  std::map<std::string, ReadNC::VarData>::iterator vmit;
  if (gDims[0] != -1) {
    if ((vmit = varInfo.find("lon")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = _readNC->read_coordinate("lon", gDims[0], gDims[3], ilVals);
      ERRORR(rval, "Trouble reading x variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find x coordinate.");
    }
  }

  // store lat values in jlVals parameterized by j
  if (gDims[1] != -1) {
    if ((vmit = varInfo.find("lat")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = _readNC->read_coordinate("lat", gDims[0], gDims[3], jlVals);
      ERRORR(rval, "Trouble reading y variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find y coordinate.");
    }
  }

  if (gDims[2] != -1) {
    if ((vmit = varInfo.find("lev")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = _readNC->read_coordinate("lev", gDims[2], gDims[5], klVals);
      ERRORR(rval, "Trouble reading z variable.");

      // decide whether down is positive
      char posval[10];
      int success = NCFUNC(get_att_text)(_readNC->fileId, (*vmit).second.varId, "positive", posval);
      if (0 == success && !strcmp(posval, "down")) {
        for (std::vector<double>::iterator dvit = klVals.begin(); dvit != klVals.end(); ++dvit)
          (*dvit) *= -1.0;
      }
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find z coordinate.");
    }
  }

  if (tMin != -1) {
    if ((vmit = varInfo.find(tName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = _readNC->read_coordinate(tName.c_str(), tMin, tMax, tVals);
      ERRORR(rval, "Trouble reading time variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find time coordinate.");
    }
  }

  if ((vmit = varInfo.find(tName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = _readNC->read_coordinate(tName.c_str(), tMin, tMax, tVals);
    ERRORR(rval, "Trouble reading time variable.");
  }
  else {
    ERRORR(MB_FAILURE, "Couldn't find time coordinate.");
  }

  // determine the entity location type of a variable
  std::map<std::string, ReadNC::VarData>::iterator mit;
  for (mit = varInfo.begin(); mit != varInfo.end(); ++mit) {
    ReadNC::VarData& vd = (*mit).second;
    if ((std::find(vd.varDims.begin(), vd.varDims.end(), iDim) != vd.varDims.end()) && (std::find(vd.varDims.begin(),
        vd.varDims.end(), kDim) != vd.varDims.end()))
      vd.entLoc = ReadNC::ENTLOCNODE;
  }

  std::copy(gDims, gDims + 6, lDims);

  // don't read coordinates of columns until we actually create the mesh

  // hack: create dummy tags, if needed, for variables like ncol and nbnd
  // with no corresponding variables
  _readNC->init_dims_with_no_cvars_info();

  return MB_SUCCESS;
}

ErrorCode NCHelperHOMME::create_HOMMEucd_verts_quads(const FileOptions& opts, EntityHandle file_set, Range& quads)
{
  Interface*& mbImpl = _readNC->mbImpl;
  std::string& fileName = _readNC->fileName;
  int& connectId = _readNC->connectId;
  int (&gDims)[6] = _readNC->gDims;
  int (&lDims)[6] = _readNC->lDims;
  std::vector<double>& ilVals = _readNC->ilVals;
  std::vector<double>& jlVals = _readNC->jlVals;
  std::vector<double>& klVals = _readNC->klVals;
  Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
  const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
  DebugOutput& dbgOut = _readNC->dbgOut;
  bool& isParallel = _readNC->isParallel;
  Range& localGid = _readNC->localGid;
#ifdef USE_MPI
  ParallelComm*& myPcomm = _readNC->myPcomm;
#endif
  bool& spectralMesh = _readNC->spectralMesh;

  // need to get/read connectivity data before creating elements
  std::string conn_fname;

  // try to open the connectivity file through CONN option, if used
  ErrorCode rval = opts.get_str_option("CONN", conn_fname);
  if (MB_SUCCESS != rval) {
    // default convention for reading HOMME is a file HommeMapping.nc in same dir as data file
    conn_fname = std::string(fileName);
    size_t idx = conn_fname.find_last_of("/");
    if (idx != std::string::npos)
      conn_fname = conn_fname.substr(0, idx).append("/HommeMapping.nc");
    else
      conn_fname = "HommeMapping.nc";
  }

  int success;

  int rank, procs;
#ifdef PNETCDF_FILE
  if (isParallel) {
    success = NCFUNC(open)(myPcomm->proc_config().proc_comm(), conn_fname.c_str(), 0, MPI_INFO_NULL, &connectId);
    rank = myPcomm->proc_config().proc_rank();
    procs = myPcomm->proc_config().proc_size();
  }
  else {
    success = NCFUNC(open)(MPI_COMM_SELF, conn_fname.c_str(), 0, MPI_INFO_NULL, &connectId);
    rank = 0;
    procs = 1;
  }
#else
  success = NCFUNC(open)(conn_fname.c_str(), 0, &connectId);
  rank = 0;
  procs = 1;
#endif
  ERRORS(success, "Failed on open.");

  std::vector<std::string> conn_names;
  std::vector<int> conn_vals;
  rval = _readNC->get_dimensions(connectId, conn_names, conn_vals);
  ERRORR(rval, "Failed to get dimensions for connectivity.");

  if (conn_vals[0] != gDims[3] - gDims[0] + 1 - 2) {
    dbgOut.tprintf(1, "Warning: number of quads from %s and vertices from %s are inconsistent; nverts = %d, nquads = %d.\n",
        conn_fname.c_str(), fileName.c_str(), gDims[3] - gDims[0] + 1, conn_vals[0]);
  }

  // read connectivity into temporary variable
  int num_fine_quads, num_coarse_quads, start_idx;
  std::vector<std::string>::iterator vit;
  int idx;
  if ((vit = std::find(conn_names.begin(), conn_names.end(), "ncells")) != conn_names.end())
  {
    idx = vit - conn_names.begin();
  }
  else if ((vit = std::find(conn_names.begin(), conn_names.end(), "ncenters")) != conn_names.end())
  {
    idx = vit - conn_names.begin();
  }
  else
  {
    ERRORR(MB_FAILURE, "Failed to get number of quads.");
  }
  int num_quads = conn_vals[idx];

  // get the connectivity into tmp_conn2 and permute into tmp_conn
  int cornerVarId;
  success = NCFUNC(inq_varid)(connectId, "element_corners", &cornerVarId);
  ERRORS(success, "Failed to get variable id.");
  NCDF_SIZE tmp_dims[2] = { 0, 0 }, tmp_counts[2] = { 4, static_cast<size_t>(num_quads) };
  std::vector<int> tmp_conn(4 * num_quads), tmp_conn2(4 * num_quads);
  success = NCFUNCAG(_vara_int)(connectId, cornerVarId, tmp_dims, tmp_counts, &tmp_conn2[0] NCREQ);
  ERRORS(success, "Failed to get temporary connectivity.");
  success = NCFUNC(close)(connectId);
  ERRORS(success, "Failed on close.");
  // permute the connectivity
  for (int i = 0; i < num_quads; i++) {
    tmp_conn[4*i] = tmp_conn2[i];
    tmp_conn[4*i + 1] = tmp_conn2[i + 1 * num_quads];
    tmp_conn[4*i + 2] = tmp_conn2[i + 2 * num_quads];
    tmp_conn[4*i + 3] = tmp_conn2[i + 3 * num_quads];
  }

  // need to know whether we'll be creating gather mesh later, to make sure we allocate enough space
  // in one shot
  bool create_gathers = true;
#ifdef USE_MPI
  if (isParallel)
    if (myPcomm->proc_config().proc_rank() != 0)
      create_gathers = false;
#endif

  // compute the number of local quads, accounting for coarse or fine representation
  // spectral_unit is the # fine quads per coarse quad, or spectralOrder^2
  int spectral_unit = (spectralMesh ? _spectralOrder * _spectralOrder : 1);
  // num_coarse_quads is the number of quads instantiated in MOAB; if !spectralMesh, num_coarse_quads = num_fine_quads
  num_coarse_quads = int(std::floor(1.0 * num_quads / (spectral_unit * procs)));
  // start_idx is the starting index in the HommeMapping connectivity list for this proc, before converting to coarse quad representation
  start_idx = 4 * rank * num_coarse_quads * spectral_unit;
  // iextra = # coarse quads extra after equal split over procs
  int iextra = num_quads % (procs * spectral_unit);
  if (rank < iextra) num_coarse_quads++;
  start_idx += 4 * spectral_unit * std::min(rank, iextra);
  // num_fine_quads is the number of quads in the connectivity list in HommeMapping file assigned to this proc
  num_fine_quads = spectral_unit * num_coarse_quads;

  // now create num_coarse_quads
  EntityHandle *conn_arr;
  EntityHandle start_vertex;
  Range tmp_range;

  // read connectivity into that space
  EntityHandle *sv_ptr = NULL, start_quad;
  SpectralMeshTool smt(mbImpl, _spectralOrder);
  if (!spectralMesh) {
    rval = _readNC->readMeshIface->get_element_connect(num_coarse_quads, 4,
                                              MBQUAD, 0, start_quad, conn_arr,
                                                // might have to create gather mesh later
                                              (create_gathers ? num_coarse_quads + num_quads : num_coarse_quads));
    ERRORR(rval, "Failed to create quads.");
    tmp_range.insert(start_quad, start_quad + num_coarse_quads - 1);
    std::copy(&tmp_conn[start_idx], &tmp_conn[start_idx + 4 * num_fine_quads], conn_arr);
    std::copy(conn_arr, conn_arr + 4 * num_fine_quads, range_inserter(localGid));
  }
  else {
    rval = smt.create_spectral_elems(&tmp_conn[0], num_fine_quads, 2, tmp_range, start_idx, &localGid);
    ERRORR(rval, "Failed to create spectral elements.");
    int count, v_per_e;
    rval = mbImpl->connect_iterate(tmp_range.begin(), tmp_range.end(), conn_arr, v_per_e, count);
    ERRORR(rval, "Failed to get connectivity of spectral elements.");
    rval = mbImpl->tag_iterate(smt.spectral_vertices_tag(true), tmp_range.begin(), tmp_range.end(),
                               count, (void*&)sv_ptr);
    ERRORR(rval, "Failed to get fine connectivity of spectral elements.");
  }

  // on this proc, I get columns ldims[1]..ldims[4], inclusive; need to find which vertices those correpond to
  unsigned int num_local_verts = localGid.size();
  unsigned int num_total_verts = gDims[3] - gDims[0] + 1;

  // create vertices
  std::vector<double*> arrays;
  rval = _readNC->readMeshIface->get_node_coords(3, num_local_verts, 0, start_vertex, arrays,
                                          // might have to create gather mesh later
                                        (create_gathers ? num_local_verts+num_total_verts : num_local_verts));
  ERRORR(rval, "Couldn't create vertices in ucd mesh.");

  // set vertex coordinates
  Range::iterator rit;
  double *xptr = arrays[0], *yptr = arrays[1], *zptr = arrays[2];
  int i;
  for (i = 0, rit = localGid.begin(); i < (int)num_local_verts; i++, ++rit) {
    assert(*rit < ilVals.size()+1);
    xptr[i] = ilVals[(*rit) - 1];
    yptr[i] = jlVals[(*rit) - 1];
    zptr[i] = klVals[lDims[2]];
  }

  const double pideg = acos(-1.0) / 180.0;
  for (i = 0; i < (int)num_local_verts; i++) {
    double cosphi = cos(pideg * yptr[i]);
    double zmult = sin(pideg * yptr[i]), xmult = cosphi * cos(xptr[i] * pideg), ymult = cosphi * sin(xptr[i] * pideg);
    double rad = 8.0e3 + klVals[lDims[2]];
    xptr[i] = rad * xmult, yptr[i] = rad * ymult, zptr[i] = rad * zmult;
  }

  // get ptr to gid memory for vertices
  Range vert_range(start_vertex, start_vertex + num_local_verts - 1);
  void *data;
  int count;
  rval = mbImpl->tag_iterate(mGlobalIdTag, vert_range.begin(), vert_range.end(), count, data);
  ERRORR(rval, "Failed to get tag iterator.");
  assert(count == (int) num_local_verts);
  int *gid_data = (int*) data;
  std::copy(localGid.begin(), localGid.end(), gid_data);
  // duplicate global id data, which will be used to resolve sharing
  if (mpFileIdTag)
  {
    rval = mbImpl->tag_iterate(*mpFileIdTag, vert_range.begin(), vert_range.end(), count, data);
    ERRORR(rval, "Failed to get tag iterator on file id tag.");
    assert(count == (int) num_local_verts);
    gid_data = (int*) data;
    std::copy(localGid.begin(), localGid.end(), gid_data);
  }

  // create map from file ids to vertex handles, used later to set connectivity
  std::map<EntityHandle, EntityHandle> vert_handles;
  for (rit = localGid.begin(), i = 0; rit != localGid.end(); ++rit, i++) {
    vert_handles[*rit] = start_vertex + i;
  }

  // compute proper handles in connectivity using offset
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

  // add new vertices and elements to the set
  quads.merge(tmp_range);
  tmp_range.insert(start_vertex, start_vertex + num_local_verts - 1);
  rval = mbImpl->add_entities(file_set, tmp_range);
  ERRORR(rval, "Couldn't add new vertices and quads/hexes to file set.");

  // mark the set with the spectral order
  Tag sporder;
  rval = mbImpl->tag_get_handle("SPECTRAL_ORDER", 1, MB_TYPE_INTEGER, sporder, MB_TAG_CREAT | MB_TAG_SPARSE);
  ERRORR(rval, "Couldn't create spectral order tag.");
  rval = mbImpl->tag_set_data(sporder, &file_set, 1, &_spectralOrder);
  ERRORR(rval, "Couldn't set value for spectral order tag.");

#ifdef USE_MPI
  if (isParallel && myPcomm->proc_config().proc_rank() == 0) {
#endif
    EntityHandle gather_set;
    rval = mbImpl->create_meshset(MESHSET_SET, gather_set);
    ERRORR(rval, "Trouble creating gather set.");

    // create vertices
    arrays.clear();
    // don't need to specify allocation number here, because we know enough verts were created before
    rval = _readNC->readMeshIface->get_node_coords(3, num_total_verts, 0, start_vertex, arrays);
    ERRORR(rval, "Couldn't create vertices in ucd mesh for gather set.");

    xptr = arrays[0], yptr = arrays[1], zptr = arrays[2];
    for (i = 0; i < (int)num_total_verts; i++) {
      double cosphi = cos(pideg * jlVals[i]);
      double zmult = sin(pideg * jlVals[i]);
      double xmult = cosphi * cos(ilVals[i] * pideg);
      double ymult = cosphi * sin(ilVals[i] * pideg);
      double rad = 8.0e3 + klVals[lDims[2]];
      xptr[i] = rad * xmult;
      yptr[i] = rad * ymult;
      zptr[i] = rad * zmult;
    }

    // get ptr to gid memory for vertices
    Range gather_verts(start_vertex, start_vertex + num_total_verts - 1);
    rval = mbImpl->tag_iterate(mGlobalIdTag, gather_verts.begin(), gather_verts.end(), count, data);
    ERRORR(rval, "Failed to get tag iterator.");
    assert(count == (int) num_total_verts);
    gid_data = (int*) data;
    for (int j = 1; j <= (int) num_total_verts; j++)
      gid_data[j - 1] = j;
    // set the file id tag too, it should be bigger something not interfering with global id
    if (mpFileIdTag) {
      rval = mbImpl->tag_iterate(*mpFileIdTag, gather_verts.begin(), gather_verts.end(), count, data);
      ERRORR(rval, "Failed to get tag iterator in file id tag.");
      assert(count == (int) num_total_verts);
      gid_data = (int*) data;
      for (int j = 1; j <= (int) num_total_verts; j++)
        gid_data[j - 1] = num_total_verts + j; // bigger than global id tag
    }

    rval = mbImpl->add_entities(gather_set, gather_verts);
    ERRORR(rval, "Couldn't add vertices to gather set.");

    // create quads
    Range gather_quads;
    // don't need to specify allocation number here, because we know enough quads were created before
    rval = _readNC->readMeshIface->get_element_connect(num_quads, 4,
                                              MBQUAD, 0, start_quad, conn_arr);
    ERRORR(rval, "Failed to create quads.");
    gather_quads.insert(start_quad, start_quad + num_quads - 1);
    std::copy(&tmp_conn[0], &tmp_conn[4 * num_quads], conn_arr);
    for (i = 0; i != 4 * num_quads; i++)
      conn_arr[i] += start_vertex - 1; // connectivity array is shifted by where the gather verts start
    rval = mbImpl->add_entities(gather_set, gather_quads);
    ERRORR(rval, "Couldn't add quads to gather set.");

    Tag gathersettag;
    rval = mbImpl->tag_get_handle("GATHER_SET", 1, MB_TYPE_INTEGER, gathersettag,
				  MB_TAG_CREAT | MB_TAG_SPARSE);
    ERRORR(rval, "Couldn't create gather set tag.");
    int gatherval = 1;
    rval = mbImpl->tag_set_data(gathersettag, &gather_set, 1, &gatherval);
    ERRORR(rval, "Couldn't set value for gather set tag.");

#ifdef USE_MPI
  }
#endif

  return MB_SUCCESS;
}

} // namespace moab
