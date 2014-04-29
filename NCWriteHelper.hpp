/*
 * NCWriteHelper.hpp
 *
 * Purpose       : Climate NC writer file helper; abstract, will be implemented for each type
 *
 *  Created on: Mar 28, 2014
 */

#ifndef NCWRITEHELPER_HPP_
#define NCWRITEHELPER_HPP_
#include "WriteNC.hpp"

namespace moab {

class NCWriteHelper
{
public:
  NCWriteHelper(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
: _writeNC(writeNC), _fileId(fileId), _opts(opts), _fileSet(fileSet),
  nTimeSteps(0), nLevels(1), tDim(-1), levDim(-1) {}
  virtual ~NCWriteHelper() {};

  //! Get appropriate helper instance for WriteNC class based on some info in the file set
  static NCWriteHelper* get_nc_helper(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet);

  //! Collect necessary info about local mesh
  virtual ErrorCode collect_mesh_info() = 0;

  //! Collect data for specified variables
  virtual ErrorCode collect_variable_data(std::vector<std::string>& var_names);

  //! Take the info from VarData and write first the coordinates, then the actual variables
  virtual ErrorCode write_values(std::vector<std::string>& var_names) = 0;

  //! Initialize file: this is where all defines are done
  //! The VarData dimension ids are filled up after define
  ErrorCode init_file(std::vector<std::string>& var_names, std::vector<std::string>& desired_names, bool _append);

protected:
  template <typename T> void jik_to_kji(size_t ni, size_t nj, size_t nk, T* dest, T* source)
  {
    size_t nik = ni * nk, nij = ni * nj;
    for (std::size_t k = 0; k != nk; k++)
      for (std::size_t j = 0; j != nj; j++)
        for (std::size_t i = 0; i != ni; i++)
          dest[k*nij + j*ni + i] = source[j*nik + i*nk + k];
  }

protected:
  //! Allow NCWriteHelper to directly access members of WriteNC
  WriteNC* _writeNC;

  //! Cache some information from WriteNC
  int _fileId;
  const FileOptions& _opts;
  EntityHandle _fileSet;

  //! Dimensions of time and level
  int nTimeSteps, nLevels;

  //! Dimension numbers for time and level
  int tDim, levDim;
};

//! Child helper class for scd mesh, e.g. CAM_EL or CAM_FV
class ScdNCWriteHelper : public NCWriteHelper
{
public:
  ScdNCWriteHelper(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
: NCWriteHelper(writeNC, fileId, opts, fileSet)
  {
    for (unsigned int i = 0; i < 6; i++) {
      lDims[i] = -1;
      lCDims[i] = -1;
    }
  }
  virtual ~ScdNCWriteHelper() {}

private:
  //! Implementation of NCWriteHelper::collect_mesh_info()
  virtual ErrorCode collect_mesh_info();

  //! Collect data for specified variables
  virtual ErrorCode collect_variable_data(std::vector<std::string>& var_names);

  //! Implementation of NCWriteHelper::write_values()
  virtual ErrorCode write_values(std::vector<std::string>& var_names);

protected:
  //! Dimensions of my local part of grid
  int lDims[6];

  //! Center dimensions of my local part of grid
  int lCDims[6];
};

//! Child helper class for ucd mesh, e.g. CAM_SE (HOMME) or MPAS
class UcdNCWriteHelper : public NCWriteHelper
{
public:
  UcdNCWriteHelper(WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet)
: NCWriteHelper(writeNC, fileId, opts, fileSet),
  nLocalCellsOwned(0), nLocalEdgesOwned(0), nLocalVerticesOwned(0),
  cDim(-1), eDim(-1), vDim(-1) {}
  virtual ~UcdNCWriteHelper() {}

protected:
  //! Dimensions of my local owned part of grid
  int nLocalCellsOwned;
  int nLocalEdgesOwned;
  int nLocalVerticesOwned;

  //! Dimension numbers for nCells, nEdges and nVertices
  int cDim, eDim, vDim;

  //! Local owned cells, edges and vertices
  Range localCellsOwned, localEdgesOwned, localVertsOwned;

  //! Local global ID for owned cells, edges and vertices
  Range localGidCellsOwned, localGidEdgesOwned, localGidVertsOwned;
};

} // namespace moab

#endif
