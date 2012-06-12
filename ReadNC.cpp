#include "ReadNC.hpp"

#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <map>

#include "moab/Core.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "FileOptions.hpp"
#include "moab/ScdInterface.hpp"

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {readMeshIface->report_error("%s", str); return rval;}
    
#define ERRORS(err, str) \
    if (err) {readMeshIface->report_error("%s", str); return MB_FAILURE;}
    
namespace moab {

ReaderIface* ReadNC::factory( Interface* iface )
  { return new ReadNC( iface ); }

ReadNC::ReadNC(Interface* impl)
        : mbImpl(impl), CPU_WORD_SIZE(-1), IO_WORD_SIZE(-1), fileId(-1), 
          tMin(-1), tMax(-1),
          iDim(-1), jDim(-1), kDim(-1), tDim(-1), numUnLim(-1), mCurrentMeshHandle(0),
          startVertex(0), startElem(0), mGlobalIdTag(0), 
          max_line_length(-1), max_str_length(-1), vertexOffset(0), dbgOut(stderr),
          isParallel(false), partMethod(-1), ucdMesh(false)
#ifdef USE_MPI
        , myPcomm(NULL)
#endif
{
  assert(impl != NULL);

  for (unsigned int i = 0; i < 6; i++) {
    gDims[i] = -1;
    lDims[i] = -1;
  }
    
  impl->query_interface(readMeshIface);
}

void ReadNC::reset()
{
  CPU_WORD_SIZE = -1;
  IO_WORD_SIZE = -1;
  fileId = -1;
  tMin = tMax = -1;
  for (unsigned int i = 0; i < 6; i++) {
    gDims[i] = -1;
    lDims[i] = -1;
  }
  
  iDim = jDim = kDim = tDim = -1;
  numUnLim = -1;
  mCurrentMeshHandle = 0;
  startVertex = startElem = 0;
  mGlobalIdTag = 0;
  max_line_length = -1;
  max_str_length = -1;
  vertexOffset = 0;
  dbgOut = stderr;
  mCurrentMeshHandle = 0;
  vertexOffset = 0; 
  ucdMesh = false;
  
#ifdef USE_MPI
  myPcomm = NULL;
#endif
}


ReadNC::~ReadNC() 
{
  mbImpl->release_interface(readMeshIface);
}
  
ErrorCode ReadNC::load_file(const char *file_name,
                            const EntityHandle* file_set,
                            const FileOptions& opts,
                            const ReaderIface::SubsetList* /*subset_list*/,
                            const Tag* file_id_tag)
{
  ErrorCode rval = MB_SUCCESS;

  //See if opts has variable(s) specified
  std::vector<std::string> var_names;
  std::vector<int> tstep_nums;
  std::vector<double> tstep_vals;

  if (file_id_tag)
    mGlobalIdTag = *file_id_tag;
  else {
      //! get and cache predefined tag handles
    int dum_val = 0;
    rval = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
                                  mGlobalIdTag, MB_TAG_DENSE|MB_TAG_CREAT, &dum_val);
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  bool nomesh = false, novars = false;
  std::string partition_tag_name;
  rval = parse_options(opts, var_names, tstep_nums, tstep_vals, nomesh, novars, partition_tag_name);
  ERRORR(rval, "Trouble parsing option string.");

  // Open the file
  dbgOut.tprintf(1, "Opening file %s\n", file_name);
  int success;
  
#ifdef PNETCDF_FILE
  if (isParallel)
    success = NCFUNC(open)(myPcomm->proc_config().proc_comm(), file_name, 0, MPI_INFO_NULL, &fileId);
  else
    success = NCFUNC(open)(MPI_COMM_SELF, file_name, 0, MPI_INFO_NULL, &fileId);
#else
  success = NCFUNC(open)(file_name, 0, &fileId);
#endif

  ERRORS(success, "Trouble opening file.");
  
    // Read the header (num dimensions, dimensions, num variables, global attribs)
  rval = read_header();
  ERRORR(rval, " ");

    // make sure there's a file set to put things in
  EntityHandle tmp_set;
  if (nomesh && !file_set) {
    ERRORR(MB_FAILURE, "NOMESH option requires non-NULL file set on input.\n");
  }
  else if (!file_set || (file_set && *file_set == 0)) {
    rval = mbImpl->create_meshset(MESHSET_SET, tmp_set);
    ERRORR(rval, "Trouble creating file set.");
  }
  else tmp_set = *file_set;
  
    // get the scd interface
  ScdInterface *scdi = NULL;
  rval = mbImpl->query_interface(scdi);
  if (!scdi) return MB_FAILURE;

    // Get bounds on ijk space
  rval = init_ijkt_vals(opts, scdi);

  if (MB_SUCCESS != rval) {
      // try reading ucd format (HOMME) data
    rval = init_ucd_mesh(opts);
    if (MB_SUCCESS == rval) {
      ucdMesh = true;
    }
  }
  ERRORR(rval, "Trouble initializing ijk values.");

    // Create mesh vertex/hex sequences
  Range hexes;
  if (nomesh && !novars) {
    rval = check_verts_hexes(tmp_set);
    ERRORR(rval, "Mesh characteristics didn't match from last read.\n");
  }
  else if (!nomesh) {
    if (ucdMesh) 
      rval = create_ucd_verts_hexes(opts, tmp_set, hexes);
    else
      rval = create_verts_hexes(scdi, tmp_set, hexes);
    ERRORR(rval, "Trouble creating vertices.");
  }

    // Read variables onto grid
  if (!novars) {
    rval = read_variables(tmp_set, var_names, tstep_nums);
    if (MB_FAILURE == rval) return rval;
  }
  
#ifdef USE_MPI
    // create partition set, and populate with elements
  if (isParallel) {
    EntityHandle partn_set;
    rval = mbImpl->create_meshset(MESHSET_SET, partn_set);
    ERRORR(rval, "Trouble creating partition set.");
    myPcomm->partition_sets().insert(partn_set);
    rval = mbImpl->add_entities(partn_set, hexes);
    ERRORR(rval, "Couldn't add new hexes to partition set.");

    Tag part_tag;
    rval = mbImpl->tag_get_handle( partition_tag_name.c_str(), 1, MB_TYPE_INTEGER, part_tag );
    if (MB_SUCCESS != rval) {
        // fall back to the partition tag
      part_tag = myPcomm->partition_tag();
    }

    int dum_rank = myPcomm->proc_config().proc_rank();
    rval = mbImpl->tag_set_data(part_tag, &partn_set, 1, &dum_rank);
    if (MB_SUCCESS != rval) return rval;
  }
#endif
  
  mbImpl->release_interface(scdi);
  ERRORR(rval, "Trouble creating scd element sequence.");
  
    // create nc conventional tags when loading header info only
  if (nomesh && novars) {
    rval = create_tags(scdi, tmp_set, tstep_nums);
    ERRORR(rval, "Trouble creating nc conventional tags.");
  }

    // close the file
  success = NCFUNC(close)(fileId);
  ERRORS(success, "Trouble closing file.");

  return MB_SUCCESS;
}

ErrorCode ReadNC::parse_options(const FileOptions &opts,
                                std::vector<std::string> &var_names, 
                                std::vector<int> &tstep_nums,
                                std::vector<double> &tstep_vals,
                                bool &nomesh,
                                bool &novars,
                                std::string &partition_tag_name) 
{
  int tmpval;
  if (MB_SUCCESS == opts.get_int_option("DEBUG_IO", 1, tmpval)) {
    dbgOut.set_verbosity(tmpval);
    dbgOut.set_prefix("NC ");
  }
  
  ErrorCode rval = opts.get_strs_option("VARIABLE", var_names ); 
  if (MB_TYPE_OUT_OF_RANGE == rval) novars = true;
  else novars = false;
  opts.get_ints_option("TIMESTEP", tstep_nums); 
  opts.get_reals_option("TIMEVAL", tstep_vals);
  rval = opts.get_null_option("NOMESH");
  if (MB_SUCCESS == rval) nomesh = true;

  if (2 <= dbgOut.get_verbosity()) {
    if (!var_names.empty()) {
      std::cerr << "Variables requested: ";
      for (unsigned int i = 0; i < var_names.size(); i++) std::cerr << var_names[i];
      std::cerr << std::endl;
    }
    if (!tstep_nums.empty()) {
      std::cerr << "Timesteps requested: ";
      for (unsigned int i = 0; i < tstep_nums.size(); i++) std::cerr << tstep_nums[i];
      std::cerr << std::endl;
    }
    if (!tstep_vals.empty()) {
      std::cerr << "Time vals requested: ";
      for (unsigned int i = 0; i < tstep_vals.size(); i++) std::cerr << tstep_vals[i];
      std::cerr << std::endl;
    }
  }
  
#ifdef USE_MPI
  rval = opts.get_option("PARTITION", partition_tag_name);
  isParallel = (rval != MB_ENTITY_NOT_FOUND);
  rval = opts.match_option("PARALLEL", "READ_PART");
  isParallel = (rval != MB_ENTITY_NOT_FOUND);
  

  if (!isParallel) 
      // return success here, since rval still has _NOT_FOUND from not finding option
    return MB_SUCCESS;
  
  int pcomm_no = 0;
  rval = opts.get_int_option("PARALLEL_COMM", pcomm_no);
  if (rval == MB_TYPE_OUT_OF_RANGE) {
    readMeshIface->report_error("Invalid value for PARALLEL_COMM option");
    return rval;
  }
  myPcomm = ParallelComm::get_pcomm(mbImpl, pcomm_no);
  if (0 == myPcomm) {
    myPcomm = new ParallelComm(mbImpl, MPI_COMM_WORLD);
  }
  const int rank = myPcomm->proc_config().proc_rank();
  dbgOut.set_rank(rank);

  const char *part_options[] = {"alljorkori", "alljkbal", "sqij", "sqjk"};
  int dum;
  rval = opts.match_option("PARTITION_METHOD", part_options, dum);
  if (rval == MB_FAILURE) {
    readMeshIface->report_error("Unknown partition method specified.");
    partMethod = ScdInterface::ALLJORKORI;
  }
  else if (rval == MB_ENTITY_NOT_FOUND)
    partMethod = ScdInterface::ALLJORKORI;
  else
    partMethod = dum;
#endif

  return MB_SUCCESS;
}
    
ErrorCode ReadNC::check_verts_hexes(EntityHandle file_set) 
{
    // check parameters on this read against what was on the mesh from last read
    // get the number of vertices 
  int num_verts;
  ErrorCode rval = mbImpl->get_number_entities_by_dimension(file_set, 0, num_verts);
  ERRORR(rval, "Trouble getting number of vertices.");
  
    // check against parameters
  int expected_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2] + 1);
  if (num_verts != expected_verts)
    ERRORR(MB_FAILURE, "Number of vertices doesn't match.");
  
    // check the number of elements too
  int num_elems;
  rval = mbImpl->get_number_entities_by_dimension(file_set, (-1 == lDims[2] ? 2 : 3), num_elems);
  ERRORR(rval, "Trouble getting number of elements.");
  
    // check against parameters
  //int expected_elems = (lDims[3] - lDims[0]) * (lDims[4] - lDims[1]) * (-1 == lDims[2] ? 1 : lDims[5] - lDims[2]);
  //if (num_elems != expected_elems)
  //ERRORR(MB_FAILURE, "Number of elements doesn't match.");
  
  return MB_SUCCESS;
}
  
ErrorCode ReadNC::create_verts_hexes(ScdInterface *scdi, EntityHandle tmp_set, Range &hexes) 
{
  Range tmp_range;
  ScdBox *scd_box;
  ErrorCode rval = scdi->construct_box(HomCoord(lDims[0], lDims[1], (-1 != lDims[2] ? lDims[2] : 0), 1),
                                       HomCoord(lDims[3], lDims[4], (-1 != lDims[5] ? lDims[5] : 0), 1),
                                       NULL, 0, scd_box, true);
  ERRORR(rval, "Trouble creating scd vertex sequence.");

    // set the global box parameters
  scd_box->set_global_box_dims(gDims);

    // set the partitioning method
  scd_box->part_method(partMethod);
  
    // add box set and new vertices, elements to the file set
  tmp_range.insert(scd_box->start_vertex(), scd_box->start_vertex()+scd_box->num_vertices()-1);
  tmp_range.insert(scd_box->start_element(), scd_box->start_element()+scd_box->num_elements()-1);
  tmp_range.insert(scd_box->box_set());
  rval = mbImpl->add_entities(tmp_set, tmp_range);
  ERRORR(rval, "Couldn't add new vertices to file set.");
  
    // get a ptr to global id memory
  void *data;
  int count;
  const Range::iterator topv = tmp_range.upper_bound(tmp_range.begin(), tmp_range.end(),
                                                     scd_box->start_vertex() + scd_box->num_vertices());
  rval = mbImpl->tag_iterate(mGlobalIdTag, tmp_range.begin(), topv, 
                             count, data);
  ERRORR(rval, "Failed to get tag iterator.");
  assert(count == scd_box->num_vertices());
  int *gid_data = (int*)data;

    // set the vertex coordinates
  double *xc, *yc, *zc;
  rval = scd_box->get_coordinate_arrays(xc, yc, zc);
  ERRORR(rval, "Couldn't get vertex coordinate arrays.");

  int i, j, k, il, jl, kl;
  int dil = lDims[3] - lDims[0] + 1;
  int djl = lDims[4] - lDims[1] + 1;
  int di = gDims[3] - gDims[0] + 1;
  int dj = gDims[4] - gDims[1] + 1;
  assert(dil == (int)ilVals.size() && djl == (int)jlVals.size() && 
         (-1 == lDims[2] || lDims[5]-lDims[2]+1 == (int)klVals.size()));
  for (kl = lDims[2]; kl <= lDims[5]; kl++) {
    k = kl - lDims[2];
    for (jl = lDims[1]; jl <= lDims[4]; jl++) {
      j = jl - lDims[1];
      for (il = lDims[0]; il <= lDims[3]; il++) {
        i = il - lDims[0];
        unsigned int pos = i + j*dil + k*dil*djl;
        xc[pos] = ilVals[i];
        yc[pos] = jlVals[j];
        zc[pos] = (-1 == lDims[2] ? 0.0 : klVals[k]);
        *gid_data = (-1 != kl ? kl*di*dj : 0) + jl*di + il + 1;
        gid_data++;
      }
    }
  }

#ifndef NDEBUG
  int num_verts = (lDims[3] - lDims[0] + 1) * (lDims[4] - lDims[1] + 1) *
    (-1 == lDims[2] ? 1 : lDims[5]-lDims[2]+1);
  std::vector<int> gids(num_verts);
  Range verts(scd_box->start_vertex(), scd_box->start_vertex()+scd_box->num_vertices()-1);
  rval = mbImpl->tag_get_data(mGlobalIdTag, verts, &gids[0]);
  ERRORR(rval, "Trouble getting gid values.");
  int vmin = *(std::min_element(gids.begin(), gids.end())),
      vmax = *(std::max_element(gids.begin(), gids.end()));
  dbgOut.tprintf(1, "Vertex gids %d-%d\n", vmin, vmax);
#endif  
    
    // add elements to the range passed in
  hexes.insert(scd_box->start_element(), scd_box->start_element() + scd_box->num_elements()-1);
  
  if (2 <= dbgOut.get_verbosity()) {
    assert(scd_box->boundary_complete());
    EntityHandle dum_ent = scd_box->start_element();
    rval = mbImpl->list_entities(&dum_ent, 1);
    ERRORR(rval, "Trouble listing first hex.");
  
    std::vector<EntityHandle> connect;
    rval = mbImpl->get_connectivity(&dum_ent, 1, connect);
    ERRORR(rval, "Trouble getting connectivity.");
  
    rval = mbImpl->list_entities(&connect[0], connect.size());
    ERRORR(rval, "Trouble listing element connectivity.");
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::create_ucd_verts_hexes(const FileOptions &opts,
                                         EntityHandle tmp_set, Range &hexes) 
{
    // need to get/read connectivity data before creating elements

    // open the connectivity file
  std::string conn_fname;
  ErrorCode rval = opts.get_str_option("CONN", conn_fname);
  if (MB_SUCCESS != rval) return rval;

  int success;
#ifdef PNETCDF_FILE
  if (isParallel)
    success = NCFUNC(open)(myPcomm->proc_config().proc_comm(), conn_fname.c_str(), 0, MPI_INFO_NULL, &connectId);
  else
    success = NCFUNC(open)(MPI_COMM_SELF, conn_fname.c_str(), 0, MPI_INFO_NULL, &connectId);
#else
  success = NCFUNC(open)(conn_fname.c_str(), 0, &connectId);
#endif
  ERRORS(success, "Failed on open.");

  std::vector<std::string> conn_names;
  std::vector<int> conn_vals;
  rval = get_dimensions(connectId, conn_names, conn_vals);
  ERRORR(rval, "Failed to get dimensions for connectivity.");
  
    // hardwire for ncenters, ncorners
  if (2 != conn_names.size() || conn_names[0] != std::string("ncenters") || conn_names[1] != std::string("ncorners"))
    ERRORR(MB_FAILURE, "Connectivity file didn't have correct dimension names.");

  if (conn_vals[0] != gDims[3]-gDims[0]+1+2) {
    dbgOut.tprintf(1, "Warning: number of quads from %s and vertices from %s are inconsistent; nverts = %d, nquads = %d.\n",
                   conn_fname.c_str(), fileName.c_str(), gDims[3]-gDims[0]+1, conn_vals[0]);
  }
  
  unsigned int num_hex_layers = lDims[5] - lDims[2];

    // read connectivity into temporary variable
  int num_quads = lDims[4]-lDims[1]+1;

  std::vector<int> tmp_conn(4*num_quads);
  NCDF_SIZE tmp_dims[2] = {0, lDims[0]}, tmp_counts[2] = {4, num_quads};
  success = NCFUNCAG(_vara_int)(connectId, 0, tmp_dims, tmp_counts, &tmp_conn[0] NCREQ);

    // on this proc, I get columns ldims[1]..ldims[4], inclusive; need to find which vertices those correpond to
  Range tmp_range;
  std::copy(tmp_conn.begin(), tmp_conn.end(), range_inserter(tmp_range));
  unsigned int num_verts = tmp_range.size();

   // create vertices
  ReadUtilIface *iface;
  rval = mbImpl->query_interface(iface);
  std::vector<double*> arrays;
  EntityHandle start_vertex, start_hex;
  rval = iface->get_node_coords(3, num_verts*(num_hex_layers+1), 0, start_vertex, arrays);
  ERRORR(rval, "Couldn't create vertices in ucd mesh.");

    // set vertex coordinates
  unsigned int i, k;
  Range::iterator rit;
    // start with first layer
  double *xptr = arrays[0], *yptr = arrays[1], *zptr = arrays[2];
  for (i = 0, rit = tmp_range.begin(); i < num_verts; i++, rit++) {
    xptr[i] = ilVals[(*rit)-1];
    yptr[i] = jlVals[(*rit)-1];
    zptr[i] = klVals[lDims[2]];
  }
  zptr += num_verts;
  for (i = lDims[2]+1; (int)i <= lDims[5]; i++) {
    std::copy(xptr, xptr+num_verts, xptr+num_verts);
    std::copy(yptr, yptr+num_verts, yptr+num_verts);
    std::fill(zptr, zptr+num_verts, klVals[i]);
    xptr += num_verts;
    yptr += num_verts;
    zptr += num_verts;
  }
  xptr = arrays[0], yptr = arrays[1], zptr = arrays[2];
  const double pideg = acos(-1.0)/180.0;
  for (i = 0; i < num_verts; i++) {
    double cosphi = cos(pideg*yptr[i]);
    double zmult = sin(pideg*yptr[i]), xmult = cosphi * cos(xptr[i]*pideg),
        ymult = cosphi * sin(xptr[i]*pideg);
    double rad = 8.0e3 + klVals[lDims[2]];
    xptr[i] = rad * xmult, yptr[i] = rad * ymult, zptr[i] = rad * zmult;
    for (k = 1; k < num_hex_layers+1; k++) {
      double thisrad = 8.0e3 + ((int)k+lDims[2] < gDims[5] ? klVals[lDims[2]+k] : klVals[gDims[5]]+1.0);
      xptr[i+k*num_verts] = thisrad * xmult;
      yptr[i+k*num_verts] = thisrad * ymult;
      zptr[i+k*num_verts] = thisrad * zmult;
    }
  }
    
    // get ptr to gid memory for vertices
  Range vert_range(start_vertex, start_vertex+num_verts-1);
  void *data;
  int count;
  rval = mbImpl->tag_iterate(mGlobalIdTag, vert_range.begin(), vert_range.end(), count, data);
  ERRORR(rval, "Failed to get tag iterator.");
  assert(count == (int) num_verts);
  int *gid_data = (int*)data;

    // copy first num_verts into gids
  std::copy(tmp_range.begin(), tmp_range.end(), gid_data);
  
    // generate subsequent layers by adding to first; use global # verts, not local number
  unsigned int global_num_verts = gDims[3] - gDims[0] + 1;
  for (i = lDims[2]+1; (int)i <= lDims[5]; i++) {
    std::transform(gid_data, gid_data+num_verts, gid_data+num_verts, std::bind2nd(std::plus<int>(),global_num_verts) );
    gid_data += num_verts;
  }
  
    // create map from file ids to vertex handles, used later to set connectivity
  std::map<unsigned int,EntityHandle> vert_handles;
  for (rit = tmp_range.begin(), i = 0; rit != tmp_range.end(); rit++, i++) {
    vert_handles[*rit] = start_vertex+i;
  }
  
    // now create hexes
  EntityHandle *conn_array;
  rval = readMeshIface->get_element_connect(num_quads*num_hex_layers, 8, MBHEX, 0, start_hex, conn_array);
  ERRORR(rval, "Failed to create hexes.");
  
    // set hex connectivity
  for (unsigned int l = 0; l < num_hex_layers; l++) {
    for (int q = 0; q < num_quads; q++) {
        // first layer 
      conn_array[0] = vert_handles[tmp_conn[q  ]]+l*num_verts;
      conn_array[1] = vert_handles[tmp_conn[q+num_quads]]+l*num_verts;
      conn_array[2] = vert_handles[tmp_conn[q+2*num_quads]]+l*num_verts;
      conn_array[3] = vert_handles[tmp_conn[q+3*num_quads]]+l*num_verts;
        // second layer
      conn_array[4] = conn_array[0] + num_verts;
      conn_array[5] = conn_array[1] + num_verts;
      conn_array[6] = conn_array[2] + num_verts;
      conn_array[7] = conn_array[3] + num_verts;
      conn_array += 8;
    }
  }

    // reuse tmp_range
  tmp_range.clear();
  tmp_range.insert(start_vertex, start_vertex+num_verts*(num_hex_layers+1)-1);
  tmp_range.insert(start_hex, start_hex+num_quads*num_hex_layers-1);
  rval = mbImpl->add_entities(tmp_set, tmp_range);
  ERRORR(rval, "Couldn't add new vertices and hexes to file set.");

    // add hexes to range passed in
  hexes.insert(start_hex, start_hex + num_quads*num_hex_layers-1);
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::read_variable_setup(std::vector<std::string> &var_names,
                                      std::vector<int> &tstep_nums, 
                                      std::vector<VarData> &vdatas) 
{
      
  std::map<std::string,VarData>::iterator mit;
  if (var_names.empty()) {
    for (mit = varInfo.begin(); mit != varInfo.end(); mit++) {
      VarData vd = (*mit).second;
      std::vector<int> tmp_v;
      if (-1 != tMin && 
          std::find(vd.varDims.begin(), vd.varDims.end(), tDim) != vd.varDims.end()) 
        tmp_v.push_back(tDim);
      if (-1 != lDims[2] && 
          std::find(vd.varDims.begin(), vd.varDims.end(), kDim) != vd.varDims.end()) 
        tmp_v.push_back(kDim);
      tmp_v.push_back(jDim);
      tmp_v.push_back(iDim);
      if (tmp_v == vd.varDims)
        vdatas.push_back(vd);
    }
  }
  else {
    for (unsigned int i = 0; i < var_names.size(); i++) {
      mit = varInfo.find(var_names[i]);
      if (mit != varInfo.end()) vdatas.push_back((*mit).second);
      else ERRORR(MB_FAILURE, "Couldn't find variable.");
    }
  }
  
  if (tstep_nums.empty() && -1 != tMin) {
    // no timesteps input, get them all
    for (int i = tMin; i <= tMax; i++) tstep_nums.push_back(i);
  }
  if (!tstep_nums.empty()) {
    for (unsigned int i = 0; i < vdatas.size(); i++) {
      vdatas[i].varTags.resize(tstep_nums.size(), 0);      
      vdatas[i].varDatas.resize(tstep_nums.size());
      vdatas[i].readDims.resize(tstep_nums.size());
      vdatas[i].readCounts.resize(tstep_nums.size());
    }
  }

  return MB_SUCCESS;
}

ErrorCode ReadNC::read_variable_allocate(std::vector<VarData> &vdatas,
                                         std::vector<int> &tstep_nums, 
                                         Range &verts) 
{
  ErrorCode rval = MB_SUCCESS;
  
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Reading variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);

        // get the tag to read into
      if (!vdatas[i].varTags[t]) {
        rval = get_tag(vdatas[i], tstep_nums[t], vdatas[i].varTags[t]);
        ERRORR(rval, "Trouble getting tag.");
      }
  
        // assume point-based values for now?
      if (-1 == tDim || dimVals[tDim] <= (int)t) {
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Wrong value for timestep number.");
      }
      else if (vdatas[i].varDims[0] != tDim) {
        ERRORR(MB_INDEX_OUT_OF_RANGE, "Non-default timestep number given for time-independent variable.");
      }
  
        // set up the dimensions and counts
        // first time
      vdatas[i].readDims[t].push_back(tstep_nums[t]);
      vdatas[i].readCounts[t].push_back(1);

        // then z/y/x
      bool have_k = false;
      if (std::find(vdatas[i].varDims.begin(), vdatas[i].varDims.end(), kDim) != vdatas[i].varDims.end()) {
        vdatas[i].readDims[t].push_back(lDims[2]); 
        vdatas[i].readCounts[t].push_back(lDims[5] - lDims[2] + 1);
        have_k = true;
      }
    
      bool have_ij = true;
#ifdef PNETCDF_FILE
        // whether we actually read anything depends on parallel and whether there's a k
      if (!have_k && -1 != kDim && 0 != myPcomm->proc_config().proc_rank())
        have_ij = false;
#else
      if (have_k); // to get rid of compiler warning
#endif    

      vdatas[i].readDims[t].push_back(lDims[1]);
      vdatas[i].readDims[t].push_back(lDims[0]);
      vdatas[i].readCounts[t].push_back(have_ij ? lDims[4]-lDims[1]+1 : 0);
      vdatas[i].readCounts[t].push_back(have_ij ? lDims[3]-lDims[0]+1 : 0);

      assert(vdatas[i].readDims[t].size() == vdatas[i].varDims.size());
  
        // get ptr to tag space
      void *data;
      int count;
      rval = mbImpl->tag_iterate(vdatas[i].varTags[t], verts.begin(), verts.end(), count, data);
      ERRORR(rval, "Failed to get tag iterator.");
      assert((unsigned)count == verts.size());
      vdatas[i].varDatas[t] = data;
    }
  }

  return rval;
}

ErrorCode ReadNC::read_variables(EntityHandle file_set, std::vector<std::string> &var_names,
                                 std::vector<int> &tstep_nums) 
{
  std::vector<VarData> vdatas;
  ErrorCode rval = read_variable_setup(var_names, tstep_nums, vdatas);
  ERRORR(rval, "Trouble setting up read variable.");
  
    // get vertices in set
  Range verts;
  rval = mbImpl->get_entities_by_dimension(file_set, 0, verts);
  ERRORR(rval, "Trouble getting vertices in set.");
  assert("Should only have a single vertex subrange, since they were read in one shot" &&
         verts.psize() == 1);

  rval = read_variable_allocate(vdatas, tstep_nums, verts);
  ERRORR(rval, "Trouble allocating read variables.");
  
  
    // finally, read into that space
  int success;
  std::vector<int> requests(vdatas.size()*tstep_nums.size()), statuss(vdatas.size()*tstep_nums.size());
  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      void *data = vdatas[i].varDatas[t];
      
      switch (vdatas[i].varDataType) {
        case NC_BYTE:
        case NC_CHAR:
            success = NCFUNCAG(_vara_text)(fileId, vdatas[i].varId, 
                                           &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0], 
                                           (char*)data NCREQ);
            ERRORS(success, "Failed to read char data.");
            break;
        case NC_DOUBLE:
            success = NCFUNCAG(_vara_double)(fileId, vdatas[i].varId, 
                                             &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0], 
                                             (double*)data NCREQ);
            ERRORS(success, "Failed to read double data.");
            break;
        case NC_FLOAT:
            success = NCFUNCAG(_vara_float)(fileId, vdatas[i].varId, 
                                            &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0], 
                                            (float*)data NCREQ);
            ERRORS(success, "Failed to read float data.");
            break;
        case NC_INT:
            success = NCFUNCAG(_vara_int)(fileId, vdatas[i].varId, 
                                          &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0], 
                                          (int*)data NCREQ);
            ERRORS(success, "Failed to read int data.");
            break;
        case NC_SHORT:
            success = NCFUNCAG(_vara_short)(fileId, vdatas[i].varId, 
                                            &vdatas[i].readDims[t][0], &vdatas[i].readCounts[t][0], 
                                            (short*)data NCREQ);
            ERRORS(success, "Failed to read short data.");
            break;
        default:
            success = 1;
      }

      if (success) ERRORR(MB_FAILURE, "Trouble reading variable.");
    }
  }

#ifdef NCWAIT
  int success = ncmpi_wait_all(fileId, requests.size(), &requests[0], &statuss[0]);
  ERRORS(success, "Failed on wait_all.");
#endif

  for (unsigned int i = 0; i < vdatas.size(); i++) {
    for (unsigned int t = 0; t < tstep_nums.size(); t++) {
      dbgOut.tprintf(2, "Converting variable %s, time step %d\n", vdatas[i].varName.c_str(), tstep_nums[t]);
      ErrorCode tmp_rval = convert_variable(file_set, vdatas[i], t);
      if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
    }
  }
    // debug output, if requested
  if (1 == dbgOut.get_verbosity()) {
    dbgOut.printf(1, "Read variables: %s", vdatas.begin()->varName.c_str());
    for (unsigned int i = 1; i < vdatas.size(); i++)
      dbgOut.printf(1, ", %s ", vdatas[i].varName.c_str());
    dbgOut.tprintf(1, "\n");
  }
  
  return rval;
}

ErrorCode ReadNC::convert_variable(EntityHandle file_set, VarData &var_data, int tstep_num) 
{
  ErrorCode rval;

    // get vertices in set
  Range verts;
  rval = mbImpl->get_entities_by_dimension(file_set, 0, verts);
  ERRORR(rval, "Trouble getting vertices in set.");
  assert("Should only have a single vertex subrange, since they were read in one shot" &&
         verts.psize() == 1);
  
    // get ptr to tag space
  void *data = var_data.varDatas[tstep_num];
  
    // finally, read into that space
  int success = 0, *idata;
  double *ddata;
  float *fdata;
  short *sdata;

  switch (var_data.varDataType) {
    case NC_FLOAT:
        ddata = (double*)var_data.varDatas[tstep_num];
        fdata = (float*)var_data.varDatas[tstep_num];
          // convert in-place
        for (int i = verts.size()-1; i >= 0; i--) 
          ddata[i] = fdata[i];
        break;
    case NC_SHORT:
        idata = (int*)var_data.varDatas[tstep_num];
        sdata = (short*)var_data.varDatas[tstep_num];
          // convert in-place
        for (int i = verts.size()-1; i >= 0; i--) 
          idata[i] = sdata[i];
        break;
    default:
        success = 1;
  }

  if (2 <= dbgOut.get_verbosity() && !success) {
    double dmin, dmax;
    int imin, imax;
    switch (var_data.varDataType) {
      case NC_DOUBLE:
      case NC_FLOAT:
          ddata = (double*)data;
          if (verts.empty()) break;

          dmin = dmax = ddata[0];
          for (unsigned int i = 1; i < verts.size(); i++) {
            if (ddata[i] < dmin) dmin = ddata[i];
            if (ddata[i] > dmax) dmax = ddata[i];
          }
          dbgOut.tprintf(2, "Variable %s (double): min = %f, max = %f\n", var_data.varName.c_str(), dmin, dmax);
          break;
      case NC_INT:
      case NC_SHORT:
          idata = (int*)data;
          if (verts.empty()) break;

          imin = imax = idata[0];
          for (unsigned int i = 1; i < verts.size(); i++) {
            if (idata[i] < imin) imin = idata[i];
            if (idata[i] > imax) imax = idata[i];
          }
          dbgOut.tprintf(2, "Variable %s (int): min = %d, max = %d\n", var_data.varName.c_str(), imin, imax);
          break;
      case NC_NAT:
      case NC_BYTE:
      case NC_CHAR:
          break;
    }
  }

  return rval;
}
    
ErrorCode ReadNC::get_tag(VarData &var_data, int tstep_num, Tag &tagh) 
{
  std::ostringstream tag_name;
  if (!tstep_num) {
    std::string tmp_name = var_data.varName + "0";
    tag_name << tmp_name.c_str();
  }
  else
    tag_name << var_data.varName << tstep_num;
  ErrorCode rval = MB_SUCCESS;
  tagh = 0;
  switch (var_data.varDataType) {
    case NC_BYTE:
    case NC_CHAR:
        rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 1, MB_TYPE_OPAQUE, tagh, MB_TAG_DENSE|MB_TAG_CREAT);
        break;
    case NC_DOUBLE:
    case NC_FLOAT:
        rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 1, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE|MB_TAG_CREAT);
        break;
    case NC_INT:
    case NC_SHORT:
        rval = mbImpl->tag_get_handle(tag_name.str().c_str(), 1, MB_TYPE_INTEGER, tagh, MB_TAG_DENSE|MB_TAG_CREAT);
        break;
    default:
        std::cerr << "Unrecognized data type for tag " << tag_name << std::endl;
        rval = MB_FAILURE;
  }
  
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.str().c_str());

  return rval;
}

ErrorCode ReadNC::init_ijkt_vals(const FileOptions &opts, ScdInterface *scdi) 
{
    // look for names of i/j/k dimensions
  gDims[0] = gDims[3] = -1;
  std::vector<std::string>::iterator vit;
  unsigned int idx;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lon")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "x1")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "x")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find i variable.");
  iDim = idx;
  gDims[3] = dimVals[idx]-1;
  gDims[0] = 0;
  iName = dimNames[idx];

  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lat")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "y1")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "y")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find j variable.");
  jDim = idx;
  gDims[4] = dimVals[idx]-1;
  gDims[1] = 0;
  jName = dimNames[idx];
  
  gDims[2] = gDims[5] = -1;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "lev")) != dimNames.end()) {
    idx = vit-dimNames.begin();
    gDims[5] = dimVals[idx]-1, gDims[2] = 0, kName = std::string("lev");
    kDim = idx;
  }
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "z")) != dimNames.end()) {
    idx = vit-dimNames.begin();
    gDims[5] = dimVals[idx]-1, gDims[2] = 0, kName = std::string("z");
    kDim = idx;
  }

  if ((vit = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "t")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find time variable.");
  tDim = idx;
  tMax = dimVals[idx]-1;
  tMin = 0;
  tName = dimNames[idx];

    // initialize parameter bounds
  std::copy(gDims, gDims+6, lDims);
  
    // parse options to get subset
  ErrorCode rval;
#ifdef USE_MPI
  if (isParallel) {
    rval = ScdInterface::compute_partition(partMethod, myPcomm->proc_config().proc_size(), 
                                           myPcomm->proc_config().proc_rank(), 
                                           gDims, lDims);
    if (MB_SUCCESS != rval) return rval;

    dbgOut.tprintf(1, "Partition: %dx%dx%d (out of %dx%dx%d)\n", 
                   lDims[3]-lDims[0]+1, lDims[4]-lDims[1]+1, lDims[5]-lDims[2]+1,
                   gDims[3]-gDims[0]+1, gDims[4]-gDims[1]+1, gDims[5]-gDims[2]+1);
    if (myPcomm->proc_config().proc_rank() == 0) 
      dbgOut.tprintf(1, "Contiguous chunks of size %d bytes.\n", 8*(lDims[3]-lDims[0]+1)*(lDims[4]-lDims[1]+1));
  }
#endif
    
  opts.get_int_option("IMIN", lDims[0]);
  opts.get_int_option("IMAX", lDims[3]);
  opts.get_int_option("JMIN", lDims[1]);
  opts.get_int_option("JMAX", lDims[4]);

  if (-1 != gDims[2]) {
    opts.get_int_option("KMIN", lDims[2]);
    opts.get_int_option("KMAX", lDims[5]);
  }
  
    // now get actual coordinate values for these dimensions
    // first allocate space...
  if (-1 != lDims[0]) ilVals.resize(lDims[3] - lDims[0] + 1);
  if (-1 != lDims[1]) jlVals.resize(lDims[4] - lDims[1] + 1);
  if (-1 != lDims[2]) klVals.resize(lDims[5] - lDims[2] + 1);
  if (-1 != tMin) tVals.resize(tMax - tMin + 1);

    // ... then read actual values
  std::map<std::string,VarData>::iterator vmit;
  if (lDims[0] != -1) {
    if ((vmit = varInfo.find(iName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(iName.c_str(), lDims[0], lDims[3], ilVals);
      ERRORR(rval, "Trouble reading x variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find x coordinate.");
    }
  }
  
  if (lDims[1] != -1) {
    if ((vmit = varInfo.find(jName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(jName.c_str(), lDims[1], lDims[4], jlVals);
      ERRORR(rval, "Trouble reading y variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find y coordinate.");
    }
  }
  
  if (lDims[2] != -1) {
    if ((vmit = varInfo.find(kName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(kName.c_str(), lDims[2], lDims[5], klVals);
      ERRORR(rval, "Trouble reading z variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find z coordinate.");
    }
  }
  
  if (tMin != -1) {
    if ((vmit = varInfo.find(tName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(tName.c_str(), tMin, tMax, tVals);
      ERRORR(rval, "Trouble reading time variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find time coordinate.");
    }
  }

  dbgOut.tprintf(1, "I=%d-%d, J=%d-%d, K=%d-%d\n", lDims[0], lDims[3], lDims[1], lDims[4], lDims[2], lDims[5]);
  dbgOut.tprintf(1, "%d elements, %d vertices\n", (lDims[3]-lDims[0])*(lDims[4]-lDims[1])*(lDims[5]-lDims[2]),
                 (lDims[3]-lDims[0]+1)*(lDims[4]-lDims[1]+1)*(lDims[5]-lDims[2]+1));
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::init_ucd_mesh(const FileOptions &opts) 
{
  ErrorCode rval;
  unsigned int idx;
  std::vector<std::string>::iterator vit;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "time")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else if ((vit = std::find(dimNames.begin(), dimNames.end(), "t")) != dimNames.end()) 
    idx = vit-dimNames.begin();
  else ERRORR(MB_FAILURE, "Couldn't find time variable.");
  tDim = idx;
  tMax = dimVals[idx]-1;
  tMin = 0;
  tName = dimNames[idx];

    // get number of vertices (labeled as number of columns) and levels
  gDims[0] = gDims[3] = -1;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "ncol")) != dimNames.end()) {
    idx = vit-dimNames.begin();
    gDims[3] = dimVals[idx] - 1;
    gDims[0] = 0;
    iDim = idx;
  }
  if (-1 == gDims[0]) return MB_FAILURE;

    // set j coordinate to the number of quads
  gDims[1] = gDims[0];
  gDims[4] = gDims[3] - 2;
  
  gDims[2] = gDims[5] = -1;
  if ((vit = std::find(dimNames.begin(), dimNames.end(), "ilev")) != dimNames.end()) {
    idx = vit-dimNames.begin();
    gDims[5] = dimVals[idx]-1, gDims[2] = 0, kName = std::string("ilev");
    kDim = idx;
  }
  if (-1 == gDims[2]) return MB_FAILURE;

    // read coordinate data
  std::map<std::string,VarData>::iterator vmit;
  if (gDims[0] != -1) {
    if ((vmit = varInfo.find("lon")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate("lon", gDims[0], gDims[3], ilVals);
      ERRORR(rval, "Trouble reading x variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find x coordinate.");
    }
  }
  
    // store lon values in jlVals parameterized by i
  if (gDims[1] != -1) {
    if ((vmit = varInfo.find("lat")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate("lat", gDims[0], gDims[3], jlVals);
      ERRORR(rval, "Trouble reading y variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find y coordinate.");
    }
  }
  
  if (gDims[2] != -1) {
    if ((vmit = varInfo.find("ilev")) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate("ilev", gDims[2], gDims[5], klVals);
      ERRORR(rval, "Trouble reading z variable.");

        // decide whether down is positive
      char posval[10];
      int success = NCFUNC(get_att_text)(fileId, (*vmit).second.varId, "positive", posval);
      if (0 == success && !strcmp(posval, "down")) {
        for (std::vector<double>::iterator vit = klVals.begin(); vit != klVals.end(); vit++)
          (*vit) *= -1.0;
      }
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find z coordinate.");
    }
  }
  
  if (tMin != -1) {
    if ((vmit = varInfo.find(tName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
      rval = read_coordinate(tName.c_str(), tMin, tMax, tVals);
      ERRORR(rval, "Trouble reading time variable.");
    }
    else {
      ERRORR(MB_FAILURE, "Couldn't find time coordinate.");
    }
  }

  if ((vmit = varInfo.find(tName)) != varInfo.end() && (*vmit).second.varDims.size() == 1) {
    rval = read_coordinate(tName.c_str(), tMin, tMax, tVals);
    ERRORR(rval, "Trouble reading time variable.");
  }
  else {
    ERRORR(MB_FAILURE, "Couldn't find time coordinate.");
  }

    // partition; must use the sqjk method here, others don't make sense
#ifdef USE_MPI
  if (isParallel) {
    rval = ScdInterface::compute_partition(ScdInterface::SQJK, 
                                           myPcomm->proc_config().proc_size(), 
                                           myPcomm->proc_config().proc_rank(), 
                                           gDims, lDims);
    if (MB_SUCCESS != rval) return rval;

    dbgOut.tprintf(1, "Partition: %dx%dx%d (out of %dx%dx%d)\n", 
                   lDims[3]-lDims[0]+1, lDims[4]-lDims[1]+1, lDims[5]-lDims[2]+1,
                   gDims[3]-gDims[0]+1, gDims[4]-gDims[1]+1, gDims[5]-gDims[2]+1);
  }
  else
    std::copy(gDims, gDims+6, lDims);
#else
    // local is same as global
  std::copy(gDims, gDims+6, lDims);
#endif

    // don't read coordinates of columns until we actually create the mesh
  return MB_SUCCESS;
}

ErrorCode ReadNC::read_coordinate(const char *var_name, int lmin, int lmax,
                                  std::vector<double> &cvals) 
{
  std::map<std::string,VarData>::iterator vmit = varInfo.find(var_name);
  if (varInfo.end() == vmit) return MB_FAILURE;
  
    // check to make sure it's a float or double
  int fail;
  NCDF_SIZE tmin = lmin, tcount = lmax - lmin + 1;
  NCDF_DIFF dum_stride = 1;
  if (NC_DOUBLE == (*vmit).second.varDataType) {
    cvals.resize(tcount);
    fail = NCFUNCA(get_vars_double)(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &cvals[0]);
    if (fail)
      ERRORS(MB_FAILURE, "Failed to get coordinate values.");
  }
  else if (NC_FLOAT == (*vmit).second.varDataType) {
    std::vector<float> tcvals(tcount);
    fail = NCFUNCA(get_vars_float)(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &tcvals[0]);
    if (fail)
      ERRORS(MB_FAILURE, "Failed to get coordinate values.");
    std::copy(tcvals.begin(), tcvals.end(), cvals.begin());
  }
  else
    ERRORR(MB_FAILURE, "Wrong data type for coordinate variable.");

  return MB_SUCCESS;
}
      
ErrorCode ReadNC::read_header()
{
  CPU_WORD_SIZE = sizeof(double);
  IO_WORD_SIZE = sizeof(double);
  
  dbgOut.tprint(1, "Reading header...\n");

    // get the global attributes
  int numgatts;
  int success;
  success = NCFUNC(inq_natts )(fileId, &numgatts);
  ERRORS(success, "Couldn't get number of global attributes.");

    // read attributes into globalAtts
  ErrorCode result = get_attributes(NC_GLOBAL, numgatts, globalAtts);
  ERRORR(result, "Getting attributes.");
  dbgOut.tprintf(1, "Read %u attributes\n", (unsigned int)globalAtts.size());

    // read in dimensions into dimVals
  result = get_dimensions(fileId, dimNames, dimVals);
  ERRORR(result, "Getting dimensions.");
  dbgOut.tprintf(1, "Read %u dimensions\n", (unsigned int)dimVals.size());

    // read in variables into varInfo
  result = get_variables();
  ERRORR(result, "Getting variables.");
  dbgOut.tprintf(1, "Read %u variables\n", (unsigned int)varInfo.size());

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_attributes(int var_id, int num_atts, std::map<std::string,AttData> &atts,
                                 const char *prefix) 
{

  char dum_name[120];

  for (int i = 0; i < num_atts; i++) {
      // get the name
    int success = NCFUNC(inq_attname)(fileId, var_id, i, dum_name);
    ERRORS(success, "Trouble getting attribute name.");
    
    AttData &data = atts[std::string(dum_name)];
    data.attName = std::string(dum_name);
    success = NCFUNC(inq_att)(fileId, var_id, dum_name, &data.attDataType, &data.attLen);
    ERRORS(success, "Trouble getting attribute info.");
    data.attVarId = var_id;

    dbgOut.tprintf(2, "%sAttribute %s: length=%u, varId=%d, type=%d\n",
                   (prefix ? prefix : ""), data.attName.c_str(), (unsigned int)data.attLen, 
                   data.attVarId, data.attDataType);
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::get_dimensions(int file_id, std::vector<std::string> &dim_names, std::vector<int> &dim_vals)
{
    // get the number of dimensions
  int num_dims;
  int success = NCFUNC(inq_ndims)(file_id, &num_dims);
  ERRORS(success, "Trouble getting number of dimensions.");

  if (num_dims > NC_MAX_DIMS) {
    readMeshIface->report_error("ReadNC: File contains %d dims but NetCDF library supports only %d\n",
                                num_dims, (int)NC_MAX_DIMS);
    return MB_FAILURE;
  }

  char dim_name[NC_MAX_NAME+1];
  NCDF_SIZE dum_len;
  dim_names.resize(num_dims);
  dim_vals.resize(num_dims);
  
  for (int i = 0; i < num_dims; i++) {
    success = NCFUNC(inq_dim)(file_id, i, dim_name, &dum_len);
    ERRORS(success, "Trouble getting dimension info.");
    
    dim_vals[i] = dum_len;
    dim_names[i] = std::string(dim_name);

    dbgOut.tprintf(2, "Dimension %s, length=%u\n",
                   dim_name, (unsigned int)dum_len);
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::get_variables() 
{
    // first cache the number of time steps
  std::vector<std::string>::iterator vit = std::find(dimNames.begin(), dimNames.end(), "time");
  if (vit == dimNames.end())
    vit = std::find(dimNames.begin(), dimNames.end(), "t");

  int ntimes = 0;
  if (vit != dimNames.end()) dimVals[vit-dimNames.begin()];
  if (!ntimes) ntimes = 1;
  
    // get the number of variables
  int num_vars;
  int success = NCFUNC(inq_nvars)(fileId, &num_vars);
  ERRORS(success, "Trouble getting number of variables.");

  if (num_vars > NC_MAX_VARS) {
    readMeshIface->report_error("ReadNC: File contains %d vars but NetCDF library supports only %d\n",
                                num_vars, (int)NC_MAX_VARS);
    return MB_FAILURE;
  }
  
  char var_name[NC_MAX_NAME+1];
  int var_ndims;
  
  for (int i = 0; i < num_vars; i++) {
      // get the name first, so we can allocate a map iterate for this var
    success = NCFUNC(inq_varname )(fileId, i, var_name);
    ERRORS(success, "Trouble getting var name.");
    VarData &data = varInfo[std::string(var_name)];
    data.varName = std::string(var_name);
    data.varId = i;
    data.varTags.resize(ntimes, 0);

      // get the data type
    success = NCFUNC(inq_vartype)(fileId, i, &data.varDataType);
    ERRORS(success, "Trouble getting variable data type.");

      // get the number of dimensions, then the dimensions
    success = NCFUNC(inq_varndims)(fileId, i, &var_ndims);
    ERRORS(success, "Trouble getting number of dims of a variable.");
    data.varDims.resize(var_ndims);

    success = NCFUNC(inq_vardimid)(fileId, i, &data.varDims[0]);
    ERRORS(success, "Trouble getting variable dimensions.");

      // finally, get the number of attributes, then the attributes
    success = NCFUNC(inq_varnatts)(fileId, i, &data.numAtts);
    ERRORS(success, "Trouble getting number of dims of a variable.");

      // print debug info here so attribute info comes afterwards
    dbgOut.tprintf(2, "Variable %s: Id=%d, numAtts=%d, datatype=%d, num_dims=%u\n",
                   data.varName.c_str(), data.varId, data.numAtts, data.varDataType, 
                   (unsigned int)data.varDims.size());

    ErrorCode rval = get_attributes(i, data.numAtts, data.varAtts, "   ");
    ERRORR(rval, "Trouble getting attributes for a variable.");

  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::read_tag_values( const char* ,
                                   const char* ,
                                   const FileOptions& ,
                                   std::vector<int>& ,
                                   const SubsetList* ) 
{
  return MB_FAILURE;
}

ErrorCode ReadNC::create_tags(ScdInterface *scdi, EntityHandle file_set, 
                              const std::vector<int>& tstep_nums)
{
  ErrorCode rval;
  std::string tag_name;
  
  // <__NUM_DIMS>
  Tag numDimsTag = 0;  
  tag_name = "__NUM_DIMS";
  int numDims = dimNames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 1, MB_TYPE_INTEGER, numDimsTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_DIMS tag.");
  rval = mbImpl->tag_set_data(numDimsTag, &file_set, 1, &numDims);
  ERRORR(rval, "Trouble setting data for __NUM_DIMS tag.");
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
  
  // <__NUM_VARS>
  Tag numVarsTag = 0;
  tag_name = "__NUM_VARS";
  int numVars = varInfo.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 1, MB_TYPE_INTEGER, numVarsTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __NUM_VARS tag.");
  rval = mbImpl->tag_set_data(numVarsTag, &file_set, 1, &numVars);
  ERRORR(rval, "Trouble setting data for __NUM_VARS tag.");
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
  
  // <__DIM_NAMES>
  Tag dimNamesTag = 0;
  tag_name = "__DIM_NAMES";
  std::string dimnames;
  unsigned int dimNamesSz = dimNames.size();
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    dimnames.append(dimNames[i]);
    dimnames.push_back('\0');
  }
  int dimnamesSz = dimnames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, dimNamesTag, MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN); 
  ERRORR(rval, "Trouble creating __DIM_NAMES tag.");
  const void* ptr = dimnames.c_str();
  rval = mbImpl->tag_set_by_ptr(dimNamesTag, &file_set, 1, &ptr, &dimnamesSz);
  ERRORR(rval, "Trouble setting data for __DIM_NAMES tag.");
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    

  // <__VAR_NAMES>
  Tag varNamesTag = 0;
  tag_name = "__VAR_NAMES";
  std::string varnames;
  std::map<std::string,VarData>::iterator mapIter;
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    varnames.append(mapIter->first);
    varnames.push_back('\0');
  }
  int varnamesSz = varnames.size();
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, varNamesTag, MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __VAR_NAMES tag.");
  ptr = varnames.c_str();
  rval = mbImpl->tag_set_by_ptr(varNamesTag, &file_set, 1, &ptr, &varnamesSz);
  ERRORR(rval, "Trouble setting data for __VAR_NAMES tag.");
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
  
  // <dim_name>
  dimNamesSz = dimNames.size();
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    tag_name = dimNames[i];
    void * val = NULL;
    if (tag_name == "lon") 
      val = &ilVals[0];
    else if (tag_name == "lat")
      val = &jlVals[0];
    else if (tag_name == "lev")
      val = &klVals[0];
    else if (tag_name == "time")
      val = &tVals[0];
    else
      continue;

    Tag tagh = 0; 
    DataType data_type;
    int val_len = dimVals[i];
    switch (varInfo[tag_name].varDataType) {
    case NC_BYTE:
    case NC_CHAR:
    case NC_DOUBLE:
    case NC_FLOAT:
      data_type = MB_TYPE_DOUBLE;
      break;
    case NC_INT:
      data_type = MB_TYPE_INTEGER;
      break;
    case NC_SHORT:
    default:
      std::cerr << "Unrecognized data type for tag " << tag_name << std::endl;
      ERRORR(MB_FAILURE, "Unrecognized data type");
    }
    rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, data_type, tagh, MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN);
    ERRORR(rval, "Trouble creating <dim_name> tag.");
    rval = mbImpl->tag_set_by_ptr(tagh, &file_set, 1, &val, &val_len);
    ERRORR(rval, "Trouble setting data for <dim_name> tag.");
    if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
  }

  // __<dim_name>_LOC_MINMAX
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    if (dimNames[i] == "lon" || dimNames[i] == "lat" || dimNames[i] == "lev") {
      std::stringstream ss_tag_name;
      ss_tag_name << "__" << dimNames[i] << "_LOC_MINMAX";
      tag_name = ss_tag_name.str();
      Tag tagh = 0; 
      std::vector<int> val(2, 0);
      if (dimNames[i] == "lon") {
	val[0] = lDims[0]; 
	val[1] = lDims[3]; 
      }
      else if (dimNames[i] == "lat") {
	val[0] = lDims[1]; 
	val[1] = lDims[4];
      }
      else if (dimNames[i] == "lev") {
	val[0] = lDims[2]; 
	val[1] = lDims[5];
      }
      else if (dimNames[i] == "time") {
	val[0] = tMin; 
	val[1] = tMax;
      }
      rval = mbImpl->tag_get_handle(tag_name.c_str(), 2, MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE|MB_TAG_CREAT);
      ERRORR(rval, "Trouble creating __<dim_name>_LOC_MINMAX tag.");
      rval = mbImpl->tag_set_data(tagh, &file_set, 1, &val[0]);
      ERRORR(rval, "Trouble setting data for __<dim_name>_LOC_MINMAX tag.");
      if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
    }
  }

  // __<dim_name>_LOC_VALS
  for (unsigned int i = 0; i != dimNamesSz; ++i) {
    if (dimNames[i] != "time")
      continue;
    std::vector<int> val;
    if (!tstep_nums.empty())
      val = tstep_nums;
    else {
      val.resize(tVals.size());
      for (unsigned int i = 0; i != tVals.size(); ++i)
	val[i] = i;
    }
    Tag tagh = 0; 
    std::stringstream ss_tag_name;
    ss_tag_name << "__" << dimNames[i] << "_LOC_VALS";
    tag_name = ss_tag_name.str();
    rval = mbImpl->tag_get_handle(tag_name.c_str(), val.size(), MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE|MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<dim_name>_LOC_VALS tag.");
    rval = mbImpl->tag_set_data(tagh, &file_set, 1, &val[0]);
    ERRORR(rval, "Trouble setting data for __<dim_name>_LOC_VALS tag.");
    if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
  }

  // __<var_name>_DIMS
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    Tag varNamesDimsTag = 0;
    std::stringstream ss_tag_name;
    ss_tag_name << "__" << mapIter->first << "_DIMS";
    tag_name = ss_tag_name.str();
    unsigned int varDimSz = varInfo[mapIter->first].varDims.size();
    if (varDimSz == 0)
      continue;
    varInfo[mapIter->first].varTags.resize(varDimSz, 0);
    for (unsigned int i = 0; i != varDimSz; ++i) {
      Tag tmptag = 0;
      std::string tmptagname = dimNames[varInfo[mapIter->first].varDims[i]];
      mbImpl->tag_get_handle(tmptagname.c_str(), 0, MB_TYPE_OPAQUE, tmptag, MB_TAG_ANY);
      varInfo[mapIter->first].varTags[i] = tmptag;
    }
    rval = mbImpl->tag_get_handle(tag_name.c_str(), varDimSz, MB_TYPE_HANDLE, varNamesDimsTag, MB_TAG_SPARSE|MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<var_name>_DIMS tag.");
    rval = mbImpl->tag_set_data(varNamesDimsTag, &file_set, 1, &(varInfo[mapIter->first].varTags[0]));
    ERRORR(rval, "Trouble setting data for __<var_name>_DIMS tag.");
    if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
  }
  
  // <PARTITION_METHOD>
  Tag part_tag = scdi->part_method_tag();
  if (!part_tag) 
    ERRORR(MB_FAILURE, "Trouble getting partition method tag.");
  rval = mbImpl->tag_set_data(part_tag, &file_set, 1, &partMethod);
  ERRORR(rval, "Trouble setting data for PARTITION_METHOD tag.");
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    

  // <__GLOBAL_ATTRIBS>
  tag_name = "__GLOBAL_ATTRIBS";
  Tag globalAttTag = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, globalAttTag, MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN);
  ERRORR(rval, "Trouble creating __GLOBAL_ATTRIBS tag.");
  std::string gattVal;
  std::vector<int> gattLen;
  rval = create_attrib_string(globalAtts, gattVal, gattLen);
  ERRORR(rval, "Trouble creating attribute strings.");
  const void* gattptr = gattVal.c_str();
  int globalAttSz = gattVal.size();
  rval = mbImpl->tag_set_by_ptr(globalAttTag, &file_set, 1, &gattptr, &globalAttSz);
  ERRORR(rval, "Trouble setting data for __GLOBAL_ATTRIBS tag.");
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
  
  // <__GLOBAL_ATTRIBS_LEN>
  tag_name ="__GLOBAL_ATTRIBS_LEN";
  Tag globalAttLenTag = 0;
  if (gattLen.size() == 0)
    gattLen.push_back(0);
  rval = mbImpl->tag_get_handle(tag_name.c_str(), gattLen.size(), MB_TYPE_INTEGER, 
				globalAttLenTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  ERRORR(rval, "Trouble creating __GLOBAL_ATTRIBS_LEN tag.");
  rval = mbImpl->tag_set_data(globalAttLenTag, &file_set, 1, &gattLen[0]);
  ERRORR(rval, "Trouble setting data for __GLOBAL_ATTRIBS_LEN tag.");
  if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    

  // __<var_name>_ATTRIBS and __<var_name>_ATTRIBS_LEN
  for (mapIter = varInfo.begin(); mapIter != varInfo.end(); ++mapIter) {
    std::stringstream ssTagName;
    ssTagName << "__" << mapIter->first << "_ATTRIBS";
    tag_name = ssTagName.str();
    Tag varAttTag = 0;
    rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, varAttTag, MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN);
    ERRORR(rval, "Trouble creating __<var_name>_ATTRIBS tag.");
    std::string varAttVal;
    std::vector<int> varAttLen;
    rval = create_attrib_string(mapIter->second.varAtts, varAttVal, varAttLen);
    ERRORR(rval, "Trouble creating attribute strings.");
    const void* varAttPtr = varAttVal.c_str();
    int varAttSz = varAttVal.size(); 
    rval = mbImpl->tag_set_by_ptr(varAttTag, &file_set, 1, &varAttPtr, &varAttSz);
    ERRORR(rval, "Trouble setting data for __<var_name>_ATTRIBS tag.");
    if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
    if (varAttLen.size() == 0)
      varAttLen.push_back(0);
    ssTagName << "_LEN";
    tag_name = ssTagName.str();
    Tag varAttLenTag = 0;
    rval = mbImpl->tag_get_handle(tag_name.c_str(), varAttLen.size(), MB_TYPE_INTEGER, 
				  varAttLenTag, MB_TAG_SPARSE|MB_TAG_CREAT);
    ERRORR(rval, "Trouble creating __<var_name>_ATTRIBS_LEN tag.");
    rval = mbImpl->tag_set_data(varAttLenTag, &file_set, 1, &varAttLen[0]);
    ERRORR(rval, "Trouble setting data for __<var_name>_ATTRIBS_LEN tag.");
    if (MB_SUCCESS == rval) dbgOut.tprintf(2, "Tag created for variable %s\n", tag_name.c_str());    
  } 
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::create_attrib_string(const std::map<std::string, AttData>& attMap, 
				       std::string& attVal,
				       std::vector<int>& attLen)
{
  int success;
  std::stringstream ssAtt;
  unsigned int sz = 0;
  std::map<std::string,AttData>::const_iterator attIt = attMap.begin();
  for (;attIt != attMap.end(); ++attIt) {
    ssAtt << attIt->second.attName;
    ssAtt << '\0';
    void* attData = NULL;
    switch (attIt->second.attDataType) {
    case NC_BYTE:
    case NC_CHAR:
      sz = attIt->second.attLen;
      attData = (char *) malloc(sz);
      success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, 
				     attIt->second.attName.c_str(), (char*)attData);
      ERRORS(success, "Failed to read attribute char data.");
      ssAtt << "char;";
      break;
    case NC_DOUBLE:
      sz = attIt->second.attLen * sizeof(double);
      attData = (double *) malloc(sz);
      success = NCFUNC(get_att_double)(fileId, attIt->second.attVarId,
				       attIt->second.attName.c_str(), (double*)attData);
      ERRORS(success, "Failed to read attribute double data.");
      ssAtt << "double;";
      break;
    case NC_FLOAT:
      sz = attIt->second.attLen * sizeof(float);
      attData = (float *) malloc(sz);
      success = NCFUNC(get_att_float)(fileId, attIt->second.attVarId,
				      attIt->second.attName.c_str(), (float*)attData);
      ERRORS(success, "Failed to read attribute float data.");
      ssAtt << "float;";
      break;
    case NC_INT:
      sz = attIt->second.attLen * sizeof(int);
      attData = (int *) malloc(sz);
      success = NCFUNC(get_att_int)(fileId, attIt->second.attVarId,
				    attIt->second.attName.c_str(), (int*)attData);
      ERRORS(success, "Failed to read attribute int data.");
      ssAtt << "int;";
      break;
    case NC_SHORT:
      sz = attIt->second.attLen * sizeof(short);
      attData = (short *) malloc(sz);
      success = NCFUNC(get_att_short)(fileId, attIt->second.attVarId,
				      attIt->second.attName.c_str(), (short*)attData);
      ERRORS(success, "Failed to read attribute short data.");
      ssAtt << "short;";
      break;
    default:
      success = 1;
    }    
    char* tmpc = (char *) attData;
    for (unsigned int counter = 0; counter != sz; ++counter)
      ssAtt << tmpc[counter];
    free(attData);
    ssAtt << ';';
    attLen.push_back(ssAtt.str().size()-1);
  }
  attVal = ssAtt.str();
  
  return MB_SUCCESS;
}


} // namespace moab
