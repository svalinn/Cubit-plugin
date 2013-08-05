//-------------------------------------------------------------------------
// Filename      : NCHelper.hpp
//
// Purpose       : Climate NC file helper
//
// Creator       : Danqing Wu
//-------------------------------------------------------------------------

#ifndef NCHELPER_HPP
#define NCHELPER_HPP

#include "ReadNC.hpp"

namespace moab {

//! Helper class to isolate reading of several different nc file formats
class NCHelper
{
public:
  NCHelper(ReadNC* readNC, int fileId) : _readNC(readNC), _fileId(fileId) {}
  virtual ~NCHelper() {}

  //! Get appropriate helper instance for ReadNC class
  static NCHelper* get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts);

  //! Interfaces to be implemented in child classes
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set) = 0;
  virtual ErrorCode check_existing_mesh(EntityHandle file_set) = 0;
  virtual ErrorCode create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& faces) = 0;
  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums) = 0;
  virtual std::string get_mesh_type_name() = 0;

protected:
  //! Read set variables, common to scd mesh and ucd mesh
  ErrorCode read_variable_to_set(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums);

  //! Convert variables in place
  ErrorCode convert_variable(ReadNC::VarData& var_data, int tstep_num);

private:
  //! Used by read_variable_to_set()
  ErrorCode read_variable_to_set_allocate(std::vector<ReadNC::VarData>& vdatas, std::vector<int>& tstep_nums);

protected:
  ReadNC* _readNC;
  int _fileId;
};

//! Child helper class for scd mesh, e.g. CAM_EL or CAM_FV
class ScdNCHelper : public NCHelper
{
public:
  ScdNCHelper(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}
  virtual ~ScdNCHelper() {}

private:
  //! Implementation of NCHelper::check_existing_mesh()
  virtual ErrorCode check_existing_mesh(EntityHandle file_set);
  //! Implementation of NCHelper::create_mesh()
  virtual ErrorCode create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& faces);
  //! Implementation of NCHelper::read_variables()
  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums);

  //! Separate set and non-set variables for scd mesh
  ErrorCode read_scd_variable_setup(std::vector<std::string>& var_names,
                                    std::vector<int>& tstep_nums,
                                    std::vector<ReadNC::VarData>& vdatas,
                                    std::vector<ReadNC::VarData>& vsetdatas);

  //! Read non-set variables for scd mesh
  ErrorCode read_scd_variable_to_nonset_allocate(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                                 std::vector<int>& tstep_nums);
  ErrorCode read_scd_variable_to_nonset(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                        std::vector<int>& tstep_nums);

  template <typename T> ErrorCode kji_to_jik(size_t ni, size_t nj, size_t nk, void* dest, T* source)
  {
    size_t nik = ni * nk, nij = ni * nj;
    T* tmp_data = reinterpret_cast<T*>(dest);
    for (std::size_t j = 0; j != nj; j++)
      for (std::size_t i = 0; i != ni; i++)
        for (std::size_t k = 0; k != nk; k++)
          tmp_data[j*nik + i*nk + k] = source[k*nij + j*ni + i];
    return MB_SUCCESS;
  }
};

//! Child helper class for ucd mesh, e.g. CAM_SE (HOMME) or MPAS
class UcdNCHelper : public NCHelper
{
public:
  UcdNCHelper(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId),
  cDim(-1), eDim(-1), vDim(-1), levDim(-1),
  nCells(0), nEdges(0), nVertices(0),
  nLocalCells(0), nLocalEdges(0), nLocalVertices(0) {}
  virtual ~UcdNCHelper() {}

private:
  //! Implementation of NCHelper::read_variables()
  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names,
                                   std::vector<int> &tstep_nums);

  //! Separate set and non-set variables for ucd mesh (implemented differently in child classes)
  virtual ErrorCode read_ucd_variable_setup(std::vector<std::string>& var_names,
                                            std::vector<int>& tstep_nums,
                                            std::vector<ReadNC::VarData>& vdatas,
                                            std::vector<ReadNC::VarData>& vsetdatas) = 0;

  //! Read non-set variables for ucd mesh (implemented differently in child classes)
  virtual ErrorCode read_ucd_variable_to_nonset_allocate(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                                         std::vector<int>& tstep_nums) = 0;
#ifdef PNETCDF_FILE
  virtual ErrorCode read_ucd_variable_to_nonset_async(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                                      std::vector<int>& tstep_nums) = 0;
#else
  virtual ErrorCode read_ucd_variable_to_nonset(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                                std::vector<int>& tstep_nums) = 0;
#endif

protected:
  //! This version takes as input the moab range, from which we actually need just the
  //! size of each sequence, for a proper transpose of the data
  template <typename T> ErrorCode kji_to_jik_stride(size_t , size_t nj, size_t nk, void* dest, T* source, Range& localGid)
  {
    std::size_t idxInSource = 0; // Position of the start of the stride
    // For each subrange, we will transpose a matrix of size
    // subrange*nj*nk (subrange takes the role of ni)
    T* tmp_data = reinterpret_cast<T*>(dest);
    for (Range::pair_iterator pair_iter = localGid.pair_begin();
        pair_iter != localGid.pair_end(); ++pair_iter) {
      std::size_t size_range = pair_iter->second - pair_iter->first + 1;
      std::size_t nik = size_range * nk, nij = size_range * nj;
      for (std::size_t j = 0; j != nj; j++)
        for (std::size_t i = 0; i != size_range; i++)
          for (std::size_t k = 0; k != nk; k++)
            tmp_data[idxInSource + j*nik + i*nk + k] = source[idxInSource + k*nij + j*size_range + i];
      idxInSource += (size_range*nj*nk);
    }
    return MB_SUCCESS;
  }

  //! Dimension numbers for nCells, nEdges, nVertices, nLevels
  int cDim, eDim, vDim, levDim;

  //! Coordinate values for vertices
  std::vector<double> xVertVals, yVertVals, zVertVals;

  int nCells;
  int nEdges;
  int nVertices;

  int nLocalCells;
  int nLocalEdges;
  int nLocalVertices;

  Range localGidCells, localGidEdges, localGidVerts;
};

} // namespace moab

#endif
