/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

#ifndef WRITENC_HPP_
#define WRITENC_HPP_


#ifndef IS_BUILDING_MB
//#error "ReadNC.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <map>
#include <set>
#include <string>

#include "moab/WriterIface.hpp"
#include "moab/ScdInterface.hpp"
#include "DebugOutput.hpp"

#ifdef USE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

#ifdef PNETCDF_FILE
#include "pnetcdf.h"
#define NCFUNC(func) ncmpi_ ## func

//! Collective I/O mode get
#define NCFUNCAG(func) ncmpi_get ## func ## _all

//! Collective I/O mode put
#define NCFUNCAP(func) ncmpi_put ## func ## _all

//! Independent I/O mode get
#define NCFUNCG(func) ncmpi_get ## func

//! Nonblocking get (request aggregation), used so far only for ucd mesh
#define NCFUNCREQG(func) ncmpi_iget ## func

#define NCDF_SIZE MPI_Offset
#define NCDF_DIFF MPI_Offset
#else
#include "netcdf.h"
#define NCFUNC(func) nc_ ## func
#define NCFUNCAG(func) nc_get ## func
#define NCFUNCAP(func) nc_put ## func
#define NCFUNCG(func) nc_get ## func
#define NCDF_SIZE size_t
#define NCDF_DIFF ptrdiff_t
#endif

namespace moab {

class WriteUtilIface;
class ScdInterface;
class NCWriteHelper;

/**
 * \brief Export NC files.
 */

class WriteNC : public WriterIface
{
  friend class NCWriteHelper;
  friend class NCWriteEuler;

public:

    //! factory method
  static WriterIface* factory( Interface* );

   //! Constructor
  WriteNC(Interface *impl);

   //! Destructor
  virtual ~WriteNC();


    //! writes out a file
  ErrorCode write_file(const char *file_name,
                         const bool overwrite,
                         const FileOptions& opts,
                         const EntityHandle *output_list,
                         const int num_sets,
                         const std::vector<std::string>& qa_list,
                         const Tag* tag_list = NULL,
                         int num_tags = 0,
                         int export_dimension = 3);


private:

  enum EntityLocation {ENTLOCVERT = 0, ENTLOCNSEDGE, ENTLOCEWEDGE, ENTLOCFACE, ENTLOCSET, ENTLOCEDGE, ENTLOCREGION};

  class AttData
  {
    public:
    AttData() : attId(-1), attLen(0), attVarId(-2) {}
    int attId;
    NCDF_SIZE attLen;
    int attVarId;
    nc_type attDataType;
    std::string attValue;
  };

  class VarData
  {
    public:VarData() : varId(-1), numAtts(-1), entLoc(ENTLOCSET), numLev(1), sz(0), has_tsteps(false) {}
    int varId;
    int numAtts;
    nc_type varDataType;
    std::vector<int> varDims; // The dimension indices making up this multi-dimensional variable
    std::map<std::string, AttData> varAtts;
    std::string varName;
    std::vector<Tag> varTags; // Tags created for this variable, e.g. one tag per timestep
    std::vector<void*> memoryHogs; // these will point to the real data; fill before writing the data
    std::vector<NCDF_SIZE> writeStarts; // Starting index for reading data values along each dimension
    std::vector<NCDF_SIZE> writeCounts; // Number of data values to be read along each dimension
    int entLoc;
    int numLev;
    int sz;
    bool has_tsteps; // Indicate whether timestep numbers are appended to tag names
  };

  // this info will be reconstructed from metadata stored on conventional fileSet tags
  //! Dimension names
  std::vector<std::string> dimNames;

  //! Dimension lengths
  std::vector<int> dimLens;

  // will collect used dimensions (ccordinate variables)
  std::set<std::string> usedCoordinates;

  //! Global attribs
  std::map<std::string, AttData> globalAtts;

  //! Variable info
  std::map<std::string, VarData> varInfo;

  ErrorCode parse_options(const FileOptions& opts, std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
                                  std::vector<double>& tstep_vals);
  /*
   * map out the header, from tags on file set; it is the inverse process from
   * ErrorCode NCHelper::create_conventional_tags
   */
  ErrorCode process_conventional_tags(EntityHandle fileSet);

  ErrorCode process_concatenated_attribute(const void * gattptr, int globalAttSz, std::vector<int> & gattLen,
      std::map<std::string, AttData> & globalAtts);

  // will collect data; it should be only on gather processor, but for the time being, collect
  // for everybody
  ErrorCode collect_variable_data( std::vector<std::string>& var_names, std::vector<int>& tstep_nums,
      std::vector<double>& tstep_vals, EntityHandle fileSet);

  // initialize file: this is where all defines are done
  // the VarData dimension ids are filled up after define
  ErrorCode initialize_file( std::vector<std::string> & var_names); // these are from options

  // take the info from VarData and write first the coordinates, then the actual variables
  ErrorCode write_values(std::vector<std::string> & var_names, EntityHandle fileSet);
    // interface instance
  Interface *mbImpl;
  WriteUtilIface* mWriteIface;

  // File var
  const char *fileName;
  int IndexFile;
  //! File numbers assigned by (p)netcdf
  int fileId;


  //! Debug stuff
  DebugOutput dbgOut;

  //! Partitioning method
  int partMethod;

  //! Scd interface
  ScdInterface* scdi;

  //! Parallel data object, to be cached with ScdBox
  ScdParData parData;

#ifdef USE_MPI
  ParallelComm* myPcomm;
#endif
  //! write options
  //
  bool noMesh;
  bool noVars;
  /*
   *  not used yet, maybe later
  bool spectralMesh;
  bool noMixedElements;
  bool noEdges;*/
  int gatherSetRank;

  //! Cached tags for writing. this will be important for ordering the data, in parallel
  //
  Tag mGlobalIdTag;

  //! Are we writing in parallel? (probably in the future)
  bool isParallel;

  // CAM Euler, etc,
  std::string grid_type;
  //! Helper class instance
  NCWriteHelper* myHelper;

};

} // namespace moab

#endif
