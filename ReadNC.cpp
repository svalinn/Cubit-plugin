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

#include "moab/Interface.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

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
  ErrorCode rval;

    /*
  int num_blocks = 0;
  const int* blocks_to_load = 0;
  if (subset_list) {
    if (subset_list->tag_list_length > 1 ||
        !strcmp( subset_list->tag_list[0].tag_name, MATERIAL_SET_TAG_NAME) ) {
      readMeshIface->report_error( "ExodusII reader supports subset read only by material ID." );
      return MB_UNSUPPORTED_OPERATION;
    }
    if (subset_list->num_parts) {
      readMeshIface->report_error( "ExodusII reader does not support mesh partitioning");
      return MB_UNSUPPORTED_OPERATION;
    }
    blocks_to_load = subset_list->tag_list[0].tag_values;
    num_blocks = subset_list->tag_list[0].num_tag_values;
  }
  
    // this function directs the reading of an exoii file, but doesn't do any of
    // the actual work
  
  //See if opts has tdata.
  ErrorCode rval;
  std::string s;
  rval = opts.get_str_option("tdata", s ); 
  if(MB_SUCCESS == rval && !s.empty())
    return update(exodus_file_name, opts, num_blocks, blocks_to_load, *file_set); 

  reset();
    */

  // 0. Open the file.
  int success = nc_open(file_name, 0, &fileId);
  ERRORS(success, "Trouble opening file.");
  
    // 1. Read the header (num dimensions, dimensions, num variables, global attribs)
  rval = read_header();
  ERRORR(rval, " ");

  rval = init_ijk_vals(opts);
  ERRORR(rval, "Trouble initializing ijk values.");
  
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
