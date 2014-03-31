#include "WriteNC.hpp"
#include "moab/CN.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/FileOptions.hpp"
#include "NCWriteHelper.hpp"

#include <fstream>
#include <map>
#include <set>

#include <iostream>
#include <sstream>

#define ERRORR(rval, str) \
  if (MB_SUCCESS != rval) { mWriteIface->report_error("%s", str); return rval; }

#define ERRORS(err, str) \
  if (err) { mWriteIface->report_error("%s", str); return MB_FAILURE; }

namespace moab {

WriterIface *WriteNC::factory( Interface* iface )
  { return new WriteNC( iface ); }

WriteNC::WriteNC(Interface *impl)
    : mbImpl(impl), dbgOut(stderr), partMethod(ScdParData::ALLJORKORI), scdi(NULL),

#ifdef USE_MPI
  myPcomm(NULL),
#endif

  noMesh(false), noVars(false), /*spectralMesh(false), noMixedElements(false), noEdges(false),*/
  gatherSetRank(-1), mGlobalIdTag(0), isParallel(false),
   myHelper(NULL)
{
  assert(impl != NULL);
  impl->query_interface(mWriteIface);
}

WriteNC::~WriteNC()
{
  mbImpl->release_interface(mWriteIface);
  if (myHelper != NULL)
    delete myHelper;
}

//! writes out a file
ErrorCode WriteNC::write_file(const char *file_name,
                                  const bool overwrite,
                                  const FileOptions& options,
                                  const EntityHandle *file_set,
                                  const int num_set,
                                  const std::vector<std::string>&,
                                  const Tag*,
                                  int,
                                  int )
{

  ErrorCode rval;
  // See if opts has variable(s) specified
  std::vector<std::string> var_names;
  std::vector<int> tstep_nums;
  std::vector<double> tstep_vals;

  // Get and cache predefined tag handles
  int dum_val = 0;
  rval = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, mGlobalIdTag, MB_TAG_DENSE , &dum_val);
  if (MB_SUCCESS != rval)
    return rval;

  // num set has to be 1, we will write only one set, the original file set used to load
  if (num_set!=1)
    ERRORR(MB_FAILURE, "we should write only one set, the file set used to read data into");

  rval = parse_options(options, var_names, tstep_nums, tstep_vals);
  ERRORR(rval, "Trouble parsing option string.");

  // important to create some data that will be used to write the file; dimensions, variables, etc
  // new variables still need to have some way of defining their dimensions
  // maybe it will be passed as write options
  rval = process_conventional_tags(*file_set);
  ERRORR(rval, "Trouble getting conventional tags.");

  // Create the file ; assume we will overwrite always, for the time being
  dbgOut.tprintf(1, "creating file %s\n", file_name);
  fileName = file_name;
  int success;

#ifdef PNETCDF_FILE
  int cmode= overwrite ? NC_CLOBBER : NC_NOCLOBBER;
  if (isParallel)
    success = NCFUNC(create)(myPcomm->proc_config().proc_comm(), file_name, cmode, MPI_INFO_NULL, &fileId);
  else
    success = NCFUNC(create)(MPI_COMM_SELF, file_name, cmode, MPI_INFO_NULL, &fileId);
#else
    // this is regular netcdf file
  success = NCFUNC(create)(file_name,  overwrite ? NC_CLOBBER : NC_NOCLOBBER, &fileId);
#endif
  ERRORS(success, "failed to create file");
  //
  rval = collect_variable_data(var_names, tstep_nums, tstep_vals, *file_set);
  ERRORR(rval, "Trouble collecting data.");

  rval = initialize_file(var_names);
  ERRORR(rval, "failed to initialize file.");

  rval = write_values( var_names, *file_set);
  ERRORR(rval, "failed to write values ");

  success =  NCFUNC(close)(fileId);
  ERRORS(success, "failed to close file");

  return MB_SUCCESS;
}

ErrorCode WriteNC::parse_options(const FileOptions& opts, std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
                                std::vector<double>& tstep_vals)
{
  int tmpval;
  if (MB_SUCCESS == opts.get_int_option("DEBUG_IO", 1, tmpval)) {
    dbgOut.set_verbosity(tmpval);
    dbgOut.set_prefix("NCWrite");
  }

  ErrorCode rval = opts.get_strs_option("VARIABLE", var_names);
  if (MB_TYPE_OUT_OF_RANGE == rval)
    noVars = true;
  else
    noVars = false;
  opts.get_ints_option("TIMESTEP", tstep_nums);
  opts.get_reals_option("TIMEVAL", tstep_vals);
  rval = opts.get_null_option("NOMESH");
  if (MB_SUCCESS == rval)
    noMesh = true;

 /* these are not used yet, maybe later
  rval = opts.get_null_option("SPECTRAL_MESH");
  if (MB_SUCCESS == rval)
    spectralMesh = true;

  rval = opts.get_null_option("NO_MIXED_ELEMENTS");
  if (MB_SUCCESS == rval)
    noMixedElements = true;

  rval = opts.get_null_option("NO_EDGES");
  if (MB_SUCCESS == rval)
    noEdges = true;*/

  if (2 <= dbgOut.get_verbosity()) {
    if (!var_names.empty()) {
      std::cerr << "Variables requested: ";
      for (unsigned int i = 0; i < var_names.size(); i++)
        std::cerr << var_names[i];
      std::cerr << std::endl;
    }
    if (!tstep_nums.empty()) {
      std::cerr << "Timesteps requested: ";
      for (unsigned int i = 0; i < tstep_nums.size(); i++)
        std::cerr << tstep_nums[i];
      std::cerr << std::endl;
    }
    if (!tstep_vals.empty()) {
      std::cerr << "Time vals requested: ";
      for (unsigned int i = 0; i < tstep_vals.size(); i++)
        std::cerr << tstep_vals[i];
      std::cerr << std::endl;
    }
  }

  // the gather set will be important in parallel, this will be the rank that will accumulate the data
  // to be written in serial
  // improvement will be for writing in true parallel
  rval = opts.get_int_option("GATHER_SET", 0, gatherSetRank);
  if (MB_TYPE_OUT_OF_RANGE == rval) {
    mWriteIface->report_error("Invalid value for GATHER_SET option");
    return rval;
  }
// FIXME: copied from readnc, may need revise
#ifdef USE_MPI
  isParallel = (opts.match_option("PARALLEL","READ_PART") != MB_ENTITY_NOT_FOUND);

  if (!isParallel)
  // Return success here, since rval still has _NOT_FOUND from not finding option
  // in this case, myPcomm will be NULL, so it can never be used; always check for isParallel
  // before any use for myPcomm
    return MB_SUCCESS;

  int pcomm_no = 0;
  rval = opts.get_int_option("PARALLEL_COMM", pcomm_no);
  if (rval == MB_TYPE_OUT_OF_RANGE) {
    mWriteIface->report_error("Invalid value for PARALLEL_COMM option");
    return rval;
  }
  myPcomm = ParallelComm::get_pcomm(mbImpl, pcomm_no);
  if (0 == myPcomm) {
    myPcomm = new ParallelComm(mbImpl, MPI_COMM_WORLD);
  }
  const int rank = myPcomm->proc_config().proc_rank();
  dbgOut.set_rank(rank);

  int dum;
  rval = opts.match_option("PARTITION_METHOD", ScdParData::PartitionMethodNames, dum);
  if (rval == MB_FAILURE) {
    mWriteIface->report_error("Unknown partition method specified.");
    partMethod = ScdParData::ALLJORKORI;
  }
  else if (rval == MB_ENTITY_NOT_FOUND)
    partMethod = ScdParData::ALLJORKORI;
  else
    partMethod = dum;
#endif

  return MB_SUCCESS;
}

// this is the inverse process to create conventional tags
// will look at <pargal_source>/src/core/fileinfo.cpp, init dim, vars, atts
ErrorCode WriteNC::process_conventional_tags(EntityHandle fileSet)
{
  ErrorCode rval;
  // start copy
  Tag tag = 0;
  std::string tag_name = "__DIM_NAMES";
  const void * data = NULL;
  int num_dims = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, tag, MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __DIM_NAMES.");
  rval = mbImpl->tag_get_by_ptr(tag, &fileSet, 1, &data, &num_dims);
  ERRORR(rval, "Failed getting values for conventional tag __DIM_NAMES.");

  dbgOut.tprintf(1, "dim names size: %d\n", num_dims);

  std::size_t start = 0;
  const char * p = static_cast<const char *>(data);
  Tag dimValsTag = 0;
  tag_name = "__DIM_LENS";
  const void * valdata = NULL;
  int num_vals = 0;
  rval =  mbImpl->tag_get_handle(tag_name.c_str(), 0,  MB_TYPE_INTEGER,
      dimValsTag,  MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __DIM_LENS.");

  rval = mbImpl->tag_get_by_ptr(dimValsTag, &fileSet, 1, &valdata, &num_vals);
  ERRORR(rval, "Failed getting values for conventional tag __DIM_LENS.");

  dbgOut.tprintf(1, "num vals in dim lens tag %d\n", num_vals);
  const int * intp = static_cast<const int*>(valdata);
  int idx = -1;

  for (std::size_t i = 0; i != static_cast<std::size_t>(num_dims); ++i)
  {
    if (p[i] == '\0')
    {
      std::string dim_name(&p[start], i - start);
      ++idx;
      int sz = intp[idx];
      dimNames.push_back(dim_name);
      dimLens.push_back(sz);
      dbgOut.tprintf(2, "dimension %s has length %d\n", dim_name.c_str(), sz);
      // fixme: need info from moab to set unlimited dimension
      // currently assume each file has the same number of time dimensions
      /*if ((dim_name == "time") || (dim_name == "Time"))
        insert(dim_name,
            *(new pcdim(dim_name, sz * m_file_names.size(), true)));
      else
        insert(dim_name, *(new pcdim(dim_name, sz)));*/
      start = i + 1;
    }
  }

  Tag tagMeshType = 0;
  tag_name = "__MESH_TYPE";
  data = NULL;
  int sz = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, tagMeshType, MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __MESH_TYPE.");
  rval = mbImpl->tag_get_by_ptr(tagMeshType, &fileSet, 1, &data, &sz);
  ERRORR(rval, "Failed getting values for conventional tag __MESH_TYPE.");

  p = static_cast<const char *>(data);
  grid_type = std::string(&p[0], sz);
  dbgOut.tprintf(2, "mesh type: %s \n", grid_type.c_str());
  // read <__VAR_NAMES_LOCATIONS> tag
  Tag varLocTag = 0;
  tag_name = "__VAR_NAMES_LOCATIONS";
  const void * loc_data = NULL;
  int loc_sz = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0,
      MB_TYPE_INTEGER, varLocTag, MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __VAR_NAMES_LOCATIONS.");
  rval = mbImpl->tag_get_by_ptr(varLocTag, &fileSet, 1, &loc_data, &loc_sz);
  ERRORR(rval, "Failed getting values for conventional tag __VAR_NAMES_LOCATIONS.");
  const int * loc_p = static_cast<const int *>(loc_data);
  std::vector<int> varLoc(loc_sz);
  std::copy(loc_p, loc_p + loc_sz, varLoc.begin());

 /*
  std::map<int, std::string> locmap;
  locmap[0] = "VERTEX";
  locmap[1] = "NSEDGE";
  locmap[2] = "EWEDGE";
  if (grid_type == "MPAS")
  {
    locmap[3] = "POLYGON";
  }
  else
  {
    locmap[3] = "QUAD";
  }
  locmap[4] = "SET";
  locmap[5] = "EDGE";
  locmap[6] = "REGION";
  */
  int nthVar = 0;

  Tag tagVarNames = 0;
  tag_name = "__VAR_NAMES";
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, tagVarNames,
      MB_TAG_ANY);
  ERRORR(rval, "Failed getting conventional tag __VAR_NAMES.");
  data = NULL;
  int num_vars = 0;
  rval = mbImpl->tag_get_by_ptr(tagVarNames, &fileSet, 1, &data, &num_vars);
  dbgOut.tprintf(2, "var names  has %d names \n", num_vars );
  ERRORR(rval, "Failed getting values for conventional tag __VAR_NAMES.");
  p = static_cast<const char *>(data);
  start = 0;
  int idxVar=0;
  for (std::size_t i = 0; i != static_cast<std::size_t>(num_vars); ++i)
  {
    if (p[i] == '\0')
    {
      std::string var_name(&p[start], i - start);

      dbgOut.tprintf(2, "var name: %s index %d \n", var_name.c_str(), idxVar);
      // process var name:
      // this will create/initiate map; we will populate variableDataStruct wit info about dims, tags, etc
      // reference & is important; otherwise variableDataStruct will go out of scope, and deleted :(
      VarData  & variableDataStruct = varInfo[var_name];

      dbgOut.tprintf(2, "at var name %s varInfo size %d \n", var_name.c_str(), varInfo.size() );

      sz = 0;
      Tag dims_tag = 0;
      std::string dim_names = "__" + var_name + "_DIMS";
      rval = mbImpl->tag_get_handle(dim_names.c_str(), 0, MB_TYPE_OPAQUE,
          dims_tag, MB_TAG_ANY);
      ERRORR(rval, "FAILED to get tag for a variable dimensions");
      //fixme: doesn't handle variables have 0 dimension
      if (rval != moab::MB_SUCCESS)
      {
        start = i + 1;
        ++nthVar;
        continue;
      }
      rval = mbImpl->tag_get_length(dims_tag, sz);
      ERRORR(rval, " size of dimensions for variable");
      dbgOut.tprintf(2, "var name: %s has %d dimensions \n", var_name.c_str(), sz);

      variableDataStruct.varDims.resize(sz);
      //std::vector<const pcdim*> dims(sz, NULL);
      const void* ptr = NULL;
      rval = mbImpl->tag_get_by_ptr(dims_tag, &fileSet, 1, &ptr);
      //

      const Tag * ptags = static_cast<const moab::Tag*>(ptr);
      for (std::size_t j = 0; j != static_cast<std::size_t>(sz); ++j)
      {
        std::string dim_name;
        rval = mbImpl->tag_get_name(ptags[j], dim_name);
        ERRORR(rval, "name of tag for dimension");
        dbgOut.tprintf(2, "var name: %s has %s as dimension \n", var_name.c_str(), dim_name.c_str() );
        std::vector<std::string>::iterator vit=std::find(dimNames.begin(), dimNames.end(), dim_name);
        if (vit==dimNames.end())
          ERRORR(MB_FAILURE, "dimension not found\n");
        variableDataStruct.varDims[j]= (int)(vit-dimNames.begin()) ; // will be used for writing
        // this will have to change to actual file dimension, for writing

        // do we have a variable for each dimension? I mean, a tag?
        //dims[j] = &(get_dim(dim_name));
      }

      // attributes for this variable
      std::stringstream ssTagName;
      ssTagName << "__" << var_name << "_ATTRIBS";
      tag_name = ssTagName.str();
      Tag varAttTag = 0;
      rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, varAttTag, MB_TAG_SPARSE | MB_TAG_VARLEN);
      ERRORR(rval, "Trouble getting __<var_name>_ATTRIBS tag.");
      std::string varAttVal;
      std::vector<int> varAttLen;
      const void* varAttPtr = 0;
      int varAttSz = 0;
      rval = mbImpl->tag_get_by_ptr(varAttTag, &fileSet, 1, &varAttPtr, &varAttSz);
      ERRORR(rval, "Trouble setting data for __<var_name>_ATTRIBS tag.");
      if (MB_SUCCESS == rval)
        dbgOut.tprintf(2, "Tag retrieved for variable %s\n", tag_name.c_str());

      ssTagName << "_LEN";
      tag_name = ssTagName.str();
      Tag varAttLenTag = 0;
      rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_INTEGER, varAttLenTag, MB_TAG_ANY);
      ERRORR(rval, "Trouble getting __<var_name>_ATTRIBS_LEN tag.");
      int varAttLenSz=0;
      rval = mbImpl->tag_get_length(varAttLenTag, varAttLenSz);
      ERRORR(rval, "Trouble getting __<var_name>_ATTRIBS_LEN length.");
      varAttLen.resize(varAttLenSz);

      rval = mbImpl->tag_get_data(varAttLenTag, &fileSet, 1, &varAttLen[0]);
      ERRORR(rval, "Trouble getting data for __<var_name>_ATTRIBS_LEN tag.");

      rval = process_concatenated_attribute(varAttPtr, varAttSz, varAttLen, variableDataStruct.varAtts);
      ERRORR(rval, " trouble processing global attributes ");

      if (MB_SUCCESS == rval)
        dbgOut.tprintf(2, "Tag metadata for variable %s\n", tag_name.c_str());
      // end attribute
      start = i + 1;
      idxVar++;
    }
  }

  // attributes
  tag_name = "__GLOBAL_ATTRIBS";
  Tag globalAttTag = 0;
  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_OPAQUE, globalAttTag, MB_TAG_SPARSE | MB_TAG_VARLEN);
  ERRORR(rval, "Trouble getting __GLOBAL_ATTRIBS tag.");
  std::string gattVal;
  std::vector<int> gattLen;

  const void* gattptr;
  int globalAttSz=0;
  rval = mbImpl->tag_get_by_ptr(globalAttTag, &fileSet, 1, &gattptr, &globalAttSz);
  ERRORR(rval, "Trouble getting data for __GLOBAL_ATTRIBS tag.");

  if (MB_SUCCESS == rval)
      dbgOut.tprintf(2, "Tag value retrieved for %s size %d\n", tag_name.c_str(), globalAttSz);

    // <__GLOBAL_ATTRIBS_LEN>
  tag_name = "__GLOBAL_ATTRIBS_LEN";
  Tag globalAttLenTag = 0;

  rval = mbImpl->tag_get_handle(tag_name.c_str(), 0, MB_TYPE_INTEGER, globalAttLenTag, MB_TAG_ANY);
  ERRORR(rval, "Trouble getting __GLOBAL_ATTRIBS_LEN tag.");
  int sizeGAtt =0;
  rval = mbImpl->tag_get_length(globalAttLenTag, sizeGAtt);
  ERRORR(rval, "Trouble getting length of __GLOBAL_ATTRIBS_LEN tag ");
  gattLen.resize(sizeGAtt);
  rval = mbImpl->tag_get_data(globalAttLenTag, &fileSet, 1, &gattLen[0]);
  ERRORR(rval, "Trouble setting data for __GLOBAL_ATTRIBS_LEN tag.");
  if (MB_SUCCESS == rval)
    dbgOut.tprintf(2, "Tag retrieved for variable %s\n", tag_name.c_str());

  rval = process_concatenated_attribute(gattptr, globalAttSz, gattLen, globalAtts);
  ERRORR(rval, " trouble processing global attributes ");

  return MB_SUCCESS;
}
// reverse process from create_attrib_string
ErrorCode WriteNC::process_concatenated_attribute(const void * gattptr, int globalAttSz, std::vector<int> & gattLen,
      std::map<std::string, AttData> & attributes)
{

  std::size_t start = 0;
  std::size_t att_counter = 0;
  std::string concatString( (char*)gattptr, (char*)gattptr+globalAttSz);

  for (std::size_t i = 0; i != (size_t)globalAttSz; ++i)
  {
    if (concatString[i] == '\0')
    {
      std::string att_name(&concatString[start], i - start);
      start = i + 1;
      while (concatString[i] != ';')
        ++i;
      std::string data_type(&concatString[start], i - start);
      ++i;
      start = i;
      i = gattLen[att_counter];
      if (concatString[i] != ';')
        ERRORR(MB_FAILURE, "Error parsing attributes ");

      std::string data_val(&concatString[start], i - start);
      start = i + 1;
      AttData attrib = attributes[att_name];
      attrib.attValue = data_val;
      ++att_counter;
      dbgOut.tprintf(2, "       Process attribute %s with value %s \n",att_name.c_str(), data_val.c_str() );
    }
  }

  return MB_SUCCESS;
}

ErrorCode WriteNC::collect_variable_data( std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
     std::vector<double>& tstep_vals, EntityHandle fileSet)
{
  // in general, in netcdf, variables that have the same name as their only dimension are called
  // coordinate variables
  // for the time being, check if all dimensions for variables are coordinate variables
  ErrorCode rval;

  usedCoordinates.clear();

  for (size_t i=0; i<var_names.size(); i++)
  {
    std::string varname=var_names[i];
    std::map<std::string, VarData>::iterator vit = varInfo.find(varname);
    if (vit==varInfo.end())
      ERRORR(MB_FAILURE, "can't find one variable");

    size_t sizeVar = 1;// get multiplied by dim lengths
    VarData  & currentVarData = vit->second;
    dbgOut.tprintf(2, "    for variable %s varDims.size %d \n", varname.c_str(), (int)currentVarData.varDims.size() );
    for (size_t j =0; j<currentVarData.varDims.size(); j++)
    {
      std::string dimName= dimNames[ currentVarData.varDims[j] ];
      vit = varInfo.find(dimName);
      if (vit==varInfo.end())
        ERRORR(MB_FAILURE, "can't find one coordinate variable");

      if( (dimName == "time" || dimName == "Time" || dimName == "t") &&
          currentVarData.varDims.size()>1 ) // so it is not time itself
        currentVarData.has_tsteps=true;

      // probably will have to look at tstep_vals to match them
      sizeVar *= dimLens[j];
      usedCoordinates.insert(dimName); // collect those used, we will need to write them to the file
      dbgOut.tprintf(2, "    for variable %s need dimension %s with length %d\n", varname.c_str(), dimName.c_str(), dimLens[ currentVarData.varDims[j] ] );
    }

    currentVarData.sz=sizeVar;

    if (currentVarData.has_tsteps)
    {

      int index=0;
      while(1)
      {
        Tag indexedTag;
        std::stringstream ssTagNameWithIndex;
        ssTagNameWithIndex << varname << index;
        rval = mbImpl->tag_get_handle(ssTagNameWithIndex.str().c_str(), indexedTag);
        if (rval!=MB_SUCCESS)
          break;
        dbgOut.tprintf(2, "    found indexed tag %d with name %s\n", index,  ssTagNameWithIndex.str().c_str());
        currentVarData.varTags.push_back(indexedTag);
        index++; // we should get out of the loop at some point
        // we will have to collect data for these tags; maybe even allocate memory again

        // the type of the tag is fixed though
        DataType type;
        rval = mbImpl->tag_get_data_type(indexedTag, type);
        ERRORR(rval, "can't get tag type");

        currentVarData.varDataType = NC_DOUBLE;
        if (type==MB_TYPE_INTEGER)
          currentVarData.varDataType = NC_INT;
      }
    }
    else
    {
      // get the tag with varname
      Tag coordtag=0;
      rval = mbImpl->tag_get_handle(varname.c_str(), coordtag);
      ERRORR(rval, "can't find one tag");
      currentVarData.varTags.push_back(coordtag); // really, only one for these
      const void * data;
      int sizeCoordinate;
      rval = mbImpl->tag_get_by_ptr(coordtag, &fileSet, 1, &data, &sizeCoordinate);
      ERRORR(rval, "can't get coordinate values");
      assert(currentVarData.memoryHogs.size()==0);// nothing so far
      currentVarData.memoryHogs.push_back((void*)data);
    }

  }

  // check that for used coordinates we have found the tags
  for (std::set<std::string>::iterator setIt = usedCoordinates.begin(); setIt!=usedCoordinates.end(); setIt++)
  {
    std::string coordName=*setIt; // deep copy

    std::map<std::string, VarData>::iterator vit = varInfo.find(coordName);
    if (vit==varInfo.end())
      ERRORR(MB_FAILURE, "can't find one coordinate variable");

    VarData & varCoordData = vit->second;
    Tag coordtag=0;
    rval = mbImpl->tag_get_handle(coordName.c_str(), coordtag);
    ERRORR(rval, "can't find one tag");
    varCoordData.varTags.push_back(coordtag); // really, only one for these

    const void * data;
    int sizeCoordinate;
    rval = mbImpl->tag_get_by_ptr(coordtag, &fileSet, 1, &data, &sizeCoordinate);
    ERRORR(rval, "can't get coordinate values");
    dbgOut.tprintf(2, "    found coordinate tag with name %s and length %d\n", coordName.c_str(),
        sizeCoordinate);
    // this is the length
    varCoordData.sz = sizeCoordinate;
    varCoordData.writeStarts.resize(1);
    varCoordData.writeStarts[0]=0;
    varCoordData.writeCounts.resize(1);
    varCoordData.writeCounts[0]=sizeCoordinate;
    // find the type of tag, and use it
    DataType type;
    rval = mbImpl->tag_get_data_type(coordtag, type);
    ERRORR(rval, "can't get tag type");

    varCoordData.varDataType = NC_DOUBLE;
    if (type==MB_TYPE_INTEGER)
      varCoordData.varDataType = NC_INT;

    assert(varCoordData.memoryHogs.size()==0);// nothing so far
    varCoordData.memoryHogs.push_back((void*)data);
  }
  return MB_SUCCESS;
}

ErrorCode WriteNC::initialize_file( std::vector<std::string> & var_names)
{
  // first initialize all coordinates, then fill VarData for actual variables (and dimensions)
  // check that for used coordinates we have found the tags
  for (std::set<std::string>::iterator setIt = usedCoordinates.begin(); setIt!=usedCoordinates.end(); setIt++)
  {
    std::string  coordName=*setIt; // deep copy

    std::map<std::string, VarData>::iterator vit = varInfo.find(coordName);
    if (vit==varInfo.end())
      ERRORR(MB_FAILURE, "can't find one coordinate variable");

    VarData & varCoordData = vit->second;

    varCoordData.varDims.resize(1);

    /* int nc_def_dim (int ncid, const char *name, size_t len, int *dimidp);
       * example:  status = nc_def_dim(fileId, "lat", 18L, &latid);
    */


    // actually define a dimension
    if (NCFUNC(def_dim)(fileId, coordName.c_str() , (size_t)varCoordData.sz ,
        &varCoordData.varDims[0]) != NC_NOERR)
     ERRORR(MB_FAILURE, "failed to generate dimension");

    dbgOut.tprintf(2, "    for coordName %s dim id is %d \n", coordName.c_str(), (int)varCoordData.varDims[0] );

    // create a variable with the same name, and its only dimension the one we just defined
    /*
     * int nc_def_var (int ncid, const char *name, nc_type xtype,
                       int ndims, const int dimids[], int *varidp);
       example: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c/nc_005fdef_005fvar.html#nc_005fdef_005fvar
     */

    if (NCFUNC(def_var)(fileId, coordName.c_str(), varCoordData.varDataType,
        1, &(varCoordData.varDims[0]), &varCoordData.varId) != NC_NOERR)
      ERRORR(MB_FAILURE, "failed to create coordinate variable");

    dbgOut.tprintf(2, "    for coordName %s variable id is %d \n", coordName.c_str(), varCoordData.varId );

  }
  // now look at requested variables, and update from the index in dimNames to the actual dimension id

  for (size_t i=0; i<var_names.size(); i++)
  {
    std::map<std::string, VarData>::iterator vit = varInfo.find(var_names[i]);
    if (vit==varInfo.end())
      ERRORR(MB_FAILURE, "can't find variable requested");

    VarData & variableData = vit->second;
    int numDims = (int) variableData.varDims.size();
    // the index is for dimNames; we need to find out the actual dimension id (from above)
    for (int j=0; j< numDims; j++)
    {
      std::string dimName = dimNames[variableData.varDims[j]];
      std::map<std::string, VarData>::iterator vit2 = varInfo.find(dimName);
      if (vit2==varInfo.end())
        ERRORR(MB_FAILURE, "can't find coordinate variable requested");

      VarData & coordData = vit2->second;
      variableData.varDims[j] = coordData.varDims[0]; // this one, being a coordinate, is the only one
      dbgOut.tprintf(2, "          dimension with index %d name %s has ID %d \n",
          j, dimName.c_str(), variableData.varDims[j]);

      variableData.writeStarts.push_back(0); // assume we will write all, so start at 0 for all dimensions
      variableData.writeCounts.push_back(coordData.sz); // again, write all; times will be one at a time
    }
    // define the variable now :
    if (NCFUNC(def_var)(fileId, var_names[i].c_str(), variableData.varDataType,
            (int)variableData.varDims.size() , &(variableData.varDims[0]),
            &variableData.varId) != NC_NOERR)
    ERRORR(MB_FAILURE, "failed to create coordinate variable");

    dbgOut.tprintf(2, "    for variable %s variable id is %d \n", var_names[i].c_str(), variableData.varId );
    // now define the variable, with all dimensions
  }

  // take it out of define mode
  if (NC_NOERR != NCFUNC(enddef)(fileId))
    ERRORR(MB_FAILURE, "failed to close define mode");

  return MB_SUCCESS;
}
ErrorCode WriteNC::write_values(std::vector<std::string> & var_names, EntityHandle fileSet)
{

  /*
  * Write an Array of Values: nc_put_vara_ type
  * int nc_put_vara_double(int ncid, int varid, const size_t start[],
                           const size_t count[], const double *dp);
  */
  // start with coordinates
  for (std::set<std::string>::iterator setIt = usedCoordinates.begin(); setIt!=usedCoordinates.end(); setIt++)
 {
   std::string  coordName=*setIt; // deep copy

   std::map<std::string, VarData>::iterator vit = varInfo.find(coordName);
   if (vit==varInfo.end())
     ERRORR(MB_FAILURE, "can't find one coordinate variable");

   VarData & varCoordData = vit->second;

   int success =0;
   switch (varCoordData.varDataType) {
      case NC_DOUBLE:
        success = NCFUNCAP(_vara_double)(fileId, varCoordData.varId,  &varCoordData.writeStarts[0],
            &varCoordData.writeCounts[0], (double*) (varCoordData.memoryHogs[0])  );
        ERRORS(success, "Failed to write double data.");
        break;
      case NC_INT:
        success = NCFUNCAP(_vara_int)(fileId, varCoordData.varId,  &varCoordData.writeStarts[0],
            &varCoordData.writeCounts[0], (int*) (varCoordData.memoryHogs[0]) );
        ERRORS(success, "Failed to write int data.");
        break;
      default:
        success = 1;
        break;
   }

 }
   // now look at requested var_names; if they have time, we will have a list, and write one at a time
  // we may also need to gather, and transpose stuff
  Range ents2d;
  // get all entities of dimension 2 from set; assume now location is on cells;
  // need to reorder stuff in the order from the file, also transpose from lev dimension
  ErrorCode rval = mbImpl->get_entities_by_dimension(fileSet, 2, ents2d);
  ERRORR(rval, "can't get entities for 2d");

  // for each variabletag in the indexed lists, write a time step data
  // assume the first dimension is time (need to check); if not, just write regularly
  //
  for (size_t i=0; i<var_names.size(); i++)
  {
    std::map<std::string, VarData>::iterator vit = varInfo.find(var_names[i]);
    if (vit==varInfo.end())
      ERRORR(MB_FAILURE, "can't find variable requested");

    VarData & variableData = vit->second;
    int numTimeSteps = (int)variableData.varTags.size();
    if (variableData.has_tsteps)
    {
      variableData.writeCounts[0] = 1; // we will write one time step
      for(int j=0; j<numTimeSteps; j++)
      {
        // we will write one time step, and count will be one; start will be different
        // we will write values directly from tag_iterate, but we should also transpose for level
        // so that means deep copy for transpose
        // !!!!!!!!!!!!!!
        //
        //  FIXME !!!!!!!!!!!
        variableData.writeStarts[0] = j; // this is time, again
        int count;
        void * dataptr;
        rval = mbImpl->tag_iterate(variableData.varTags[j], ents2d.begin(), ents2d.end(), count, dataptr );
        assert(count == ents2d.size());

        // now write from memory directly   FIXME: we need to transpose and gather for multiple processors
        int success = 0;
        switch (variableData.varDataType)
        {
        case NC_DOUBLE:
          success = NCFUNCAP(_vara_double)(fileId, variableData.varId,
              &variableData.writeStarts[0], &variableData.writeCounts[0],
              (double*) (dataptr));
          ERRORS(success, "Failed to write double data.");
          break;
        case NC_INT:
          success =NCFUNCAP(_vara_int)(fileId, variableData.varId,
              &variableData.writeStarts[0], &variableData.writeCounts[0],
              (int*) (dataptr));
          ERRORS(success, "Failed to write int data.");
          break;
        default:
          success = 1;
          break;
        }

      }
    }
    else
    {
      // FIXME
    }
  }

  return MB_SUCCESS;
}

} // end moab namespace
