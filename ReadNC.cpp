#include "ReadNC.hpp"
#include "netcdf.h"

#include <algorithm>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <sstream>
#include <map>

#include "moab/Core.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "SequenceManager.hpp"
#include "VertexSequence.hpp"

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {readMeshIface->report_error(str); return rval;}
    
#define ERRORS(err, str) \
    if (err) {readMeshIface->report_error(str); return MB_FAILURE;}
    


namespace moab {

ReaderIface* ReadNC::factory( Interface* iface )
  { return new ReadNC( iface ); }

ReadNC::ReadNC(Interface* impl)
    : mbImpl(impl), max_line_length(-1), max_str_length(-1)
{
  assert(impl != NULL);
  reset();
  
  void* ptr = 0;
  impl->query_interface( "ReadUtilIface", &ptr );
  readMeshIface = reinterpret_cast<ReadUtilIface*>(ptr);

  // initialize in case tag_get_handle fails below
  mGlobalIdTag     = 0;

  //! get and cache predefined tag handles
  int dum_val = 0;
  ErrorCode result = impl->tag_get_handle(GLOBAL_ID_TAG_NAME, mGlobalIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, 
                              MB_TYPE_INTEGER, mGlobalIdTag, &dum_val);
}

void ReadNC::reset()
{
  numDims = -1;
  mCurrentMeshHandle = 0;
  vertexOffset = 0; 
}


ReadNC::~ReadNC() 
{
  std::string iface_name = "ReadUtilIface";
  mbImpl->release_interface(iface_name, readMeshIface);
}
  
ErrorCode ReadNC::load_file(const char *file_name,
                            const EntityHandle* file_set,
                            const FileOptions& opts,
                            const ReaderIface::SubsetList* subset_list,
                            const Tag* file_id_tag)
{
  ErrorCode rval = MB_SUCCESS;

  //See if opts has variable(s) specified
  std::vector<std::string> var_names;
  std::vector<int> tstep_nums;
  std::vector<double> tstep_vals;
//  bool nomesh = false;
  
//  rval = parse_options(opts, var_names, tstep_nums, tstep_vals, nomesh);
  ERRORR(rval, "Trouble parsing option string.");

  // 0. Open the file.
  int success = nc_open(file_name, 0, &fileId);
  ERRORS(success, "Trouble opening file.");
  
    // 1. Read the header (num dimensions, dimensions, num variables, global attribs)
  rval = read_header();
  ERRORR(rval, " ");

  rval = init_ijk_vals(opts);
  ERRORR(rval, "Trouble initializing ijk values.");

  rval = create_verts();
  ERRORR(rval, "Trouble creating vertices.");
/*  
  status = mdbImpl->get_entities_by_handle(0, initRange);
  if (MB_FAILURE == status) return status;

    // 2. Read/create nodes and elements
  status = read_nodes_elements(file_id_tag);
  if (MB_FAILURE == status) return status;
 
    // 3. Read variables onto grid
  status = read_variables(file_id_tag);
  if (MB_FAILURE == status) return status;
*/
  return MB_SUCCESS;
}

/*
ErrorCode ReadNC::parse_options(FileOptions &opts,
                                std::vector<std::string> &var_names, 
                                std::vector<int> &tstep_nums,
                                std::vector<double> &tstep_vals,
                                bool &nomesh) 
{
  std::string s;
  ErrorCode rval = opts.get_str_option("variable", s ); 
  char *tmp_str;
  if (MB_SUCCESS == rval) {
    if (s == "MOAB_ALL_VARIABLES") {
      allVars = true;
    }
    else {
      tmp_str = strdup(s.c_str());
      for (char* i = strtok(tmp_str, ","); i; i = strtok(NULL, ",")) 
        if (!strempty(i)) // skip empty strings
          var_names.push_back(std::string(i));
    }
  }
  
  rval = opts.get_ints_option("timestep", tstep_nums); 
      
*/

ErrorCode ReadNC::create_verts() 
{
  Core *tmpImpl = dynamic_cast<Core*>(mbImpl);
  SequenceManager *seq_mgr = tmpImpl->sequence_manager();
  VertexSequence *vert_seq;
  EntitySequence *dum_seq;
  ErrorCode rval = seq_mgr->create_scd_sequence(ilMin, jlMin, (-1 != klMin ? klMin : 0),
                                                ilMax, jlMax, (-1 != klMax ? klMax : 0),
                                                MBVERTEX, (moab::EntityID)0, startVertex, dum_seq);
  ERRORR(rval, "Trouble creating scd vertex sequence.");
  
  vert_seq = dynamic_cast<VertexSequence*>(dum_seq);
  
    // set the vertex coordinates
  double *xc, *yc, *zc;
  rval = vert_seq->get_coordinate_arrays(xc, yc, zc);
  ERRORR(rval, "Couldn't get vertex coordinate arrays.");

  return MB_SUCCESS;
}

ErrorCode ReadNC::init_ijk_vals(const FileOptions &opts) 
{
    // look for names of i/j/k dimensions
  iMin = iMax = -1;
  std::map<std::string,int>::iterator mit;
  if ((mit = dimVals.find("lon")) != dimVals.end()) 
    iMax = (*mit).second-1, iMin = 0, iName = (*mit).first;
  else if ((mit = dimVals.find("x1")) != dimVals.end()) 
    iMax = (*mit).second-1, iMin = 0, iName = (*mit).first;
  else ERRORR(MB_FAILURE, "Couldn't find i variable.");

  jMin = jMax = -1;
  if ((mit = dimVals.find("lat")) != dimVals.end()) 
    jMax = (*mit).second-1, jMin = 0, jName = (*mit).first;
  else if ((mit = dimVals.find("y1")) != dimVals.end()) 
    jMax = (*mit).second-1, jMin = 0, jName = (*mit).first;
  else ERRORR(MB_FAILURE, "Couldn't find j variable.");
  
  kMin = kMax = -1;
  if ((mit = dimVals.find("lev")) != dimVals.end()) 
    kMax = (*mit).second-1, kMin = 0, kName = (*mit).first;

    // parse options to get subset
    // FOR NOW: just take whole thing
  ilMin = iMin, ilMax = iMax;
  jlMin = jMin, jlMax = jMax;
  klMin = kMin, klMax = kMax;
  
    // now get actual coordinate values for these dimensions
    // first allocate space...
  if (-1 != ilMin) ilVals.resize(ilMax - ilMin + 1);
  if (-1 != jlMin) jlVals.resize(jlMax - jlMin + 1);
  if (-1 != klMin) klVals.resize(klMax - klMin + 1);

    // ... then read actual values
  ErrorCode rval;
  if (ilMin != -1) {
    rval = read_coordinate("lon", ilMin, ilMax, ilVals);
    if (MB_SUCCESS != rval) 
      rval = read_coordinate("x1", ilMin, ilMax, ilVals);
    ERRORR(rval, "Couldn't find x coordinate.");
  }
  
  if (jlMin != -1) {
    rval = read_coordinate("lat", jlMin, jlMax, jlVals);
    if (MB_SUCCESS != rval) 
      rval = read_coordinate("y1", jlMin, jlMax, jlVals);
    ERRORR(rval, "Couldn't find y coordinate.");
  }

  if (klMin != -1) {
    rval = read_coordinate("lev", klMin, klMax, klVals);
  }
  
  return MB_SUCCESS;
}
    
ErrorCode ReadNC::read_coordinate(const char *var_name, int lmin, int lmax,
                                  std::vector<double> &cvals) 
{
  std::map<std::string,VarData>::iterator vmit = varInfo.find(var_name);
  if (varInfo.end() == vmit) return MB_FAILURE;
  
    // check to make sure it's a float or double
  int success;
  size_t tmin = lmin, tcount = lmax - lmin + 1;
  ptrdiff_t dum_stride = 1;
  if (NC_DOUBLE == (*vmit).second.varDataType) {
    cvals.resize(tcount);
    success = nc_get_vars_double(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &cvals[0]);
    ERRORS(success, "Failed to get coordinate values.");
  }
  else if (NC_FLOAT == (*vmit).second.varDataType) {
    std::vector<float> tcvals(tcount);
    success = nc_get_vars_float(fileId, (*vmit).second.varId, &tmin, &tcount, &dum_stride, &tcvals[0]);
    ERRORS(success, "Failed to get coordinate values.");
    std::copy(tcvals.begin(), tcvals.end(), cvals.begin());
  }
  else
    ERRORR(MB_FAILURE, "Wrong data type for coordinate variable.");

  return MB_SUCCESS;
}
      
ErrorCode ReadNC::read_header()
{
  CPU_WORD_SIZE = sizeof(double);  // With ExodusII version 2, all floats
  IO_WORD_SIZE = sizeof(double);   // should be changed to doubles
  
    // get the global attributes
  int numgatts;
  int success = nc_inq_natts (fileId, &numgatts);
  ERRORS(success, "Couldn't get number of global attributes.");
  
  ErrorCode result = get_attributes(numgatts, NC_GLOBAL, globalAtts);
  ERRORR(result, "Getting attributes.");

    // read in dimensions
  result = get_dimensions();
  ERRORR(result, "Getting dimensions.");

    // read in variables
  result = get_variables();
  ERRORR(result, "Getting variables.");

  return MB_SUCCESS;
}

ErrorCode ReadNC::get_attributes(int var_id, int num_atts, std::map<std::string,AttData> &atts) 
{

  char dum_name[120];

  for (int i = 0; i < num_atts; i++) {
      // get the name
    int success = nc_inq_attname(fileId, NC_GLOBAL, i, dum_name);
    ERRORS(success, "Trouble getting attribute name.");
    
    AttData &data = globalAtts[std::string(dum_name)];
    success = nc_inq_att(fileId, NC_GLOBAL, dum_name, &data.attDataType, &data.attLen);
    ERRORS(success, "Trouble getting attribute info.");
    data.attVarId = NC_GLOBAL;
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::get_dimensions() 
{
    // get the number of dimensions
  int success = nc_inq_ndims(fileId, &numDims);
  ERRORS(success, "Trouble getting number of dimensions.");

  if (numDims > NC_MAX_DIMS) {
    readMeshIface->report_error("ReadNC: File contains %d dims but NetCDF library supports only %d\n",
                                numDims, (int)NC_MAX_DIMS);
    return MB_FAILURE;
  }

  char dim_name[NC_MAX_NAME+1];
  size_t dum_len;
  
  for (int i = 0; i < numDims; i++) {
    success = nc_inq_dim(fileId, i, dim_name, &dum_len);
    ERRORS(success, "Trouble getting dimension info.");
    
    dimVals[std::string(dim_name)] = dum_len;
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::get_variables() 
{
    // get the number of variables
  int success = nc_inq_nvars(fileId, &numVars);
  ERRORS(success, "Trouble getting number of variables.");

  if (numVars > NC_MAX_VARS) {
    readMeshIface->report_error("ReadNC: File contains %d vars but NetCDF library supports only %d\n",
                                numVars, (int)NC_MAX_VARS);
    return MB_FAILURE;
  }
  
  char var_name[NC_MAX_NAME+1];
  int var_ndims, var_natts;
  
  for (int i = 0; i < numVars; i++) {
      // get the name first, so we can allocate a map iterate for this var
    success = nc_inq_varname (fileId, i, var_name);
    ERRORS(success, "Trouble getting var name.");
    VarData &data = varInfo[std::string(var_name)];

      // get the data type
    success = nc_inq_vartype(fileId, i, &data.varDataType);
    ERRORS(success, "Trouble getting variable data type.");

      // get the number of dimensions, then the dimensions
    success = nc_inq_varndims(fileId, i, &var_ndims);
    ERRORS(success, "Trouble getting number of dims of a variable.");
    data.varDims.resize(var_ndims);

    success = nc_inq_vardimid(fileId, i, &data.varDims[0]);
    ERRORS(success, "Trouble getting variable dimensions.");

      // finally, get the number of attributes, then the attributes
    success = nc_inq_varnatts(fileId, i, &var_natts);
    ERRORS(success, "Trouble getting number of dims of a variable.");

    ErrorCode rval = get_attributes(i, var_natts, data.varAtts);
    ERRORR(rval, "Trouble getting attributes for a variable.");
  }
  
  return MB_SUCCESS;
}

ErrorCode ReadNC::read_tag_values( const char* file_name,
                                   const char* tag_name,
                                   const FileOptions& opts,
                                   std::vector<int>& tag_values_out,
                                   const SubsetList* subset_list) 
{
  return MB_FAILURE;
}

} // namespace moab
