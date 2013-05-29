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

  static NCHelper* get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts);

  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set) = 0;

  virtual ErrorCode create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads) = 0;

  virtual ErrorCode init_local_gid(EntityHandle file_set) = 0;

  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums) = 0;

  virtual std::string get_mesh_type_name() = 0;

protected:
  ReadNC* _readNC;

  int _fileId;
};

//! Child helper class for Eulerian Spectral grid (CAM_EUL)
class NCHEuler : public NCHelper
{
public:
  NCHEuler(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}

  static bool can_read_file(ReadNC* readNC);

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);

  virtual ErrorCode create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads)
  { return _readNC->create_scd_verts_quads(scdi, file_set, quads); }

  virtual ErrorCode init_local_gid(EntityHandle file_set)
  { return MB_SUCCESS; }

  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
  { return _readNC->read_variables(file_set, var_names, tstep_nums, true); }

  virtual std::string get_mesh_type_name()
  { return "CAM_EUL"; }
};

//! Child helper class for Finite Volume grid (CAM_FV)
class NCHFV : public NCHelper
{
public:
  NCHFV(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}

  static bool can_read_file(ReadNC* readNC);

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);

  virtual ErrorCode create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads)
  { return _readNC->create_scd_verts_quads(scdi, file_set, quads); }

  virtual ErrorCode init_local_gid(EntityHandle file_set)
  { return MB_SUCCESS; }

  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
  { return _readNC->read_variables(file_set, var_names, tstep_nums, true); }

  virtual std::string get_mesh_type_name()
  { return "CAM_FV"; }
};

//! Child helper class for HOMME grid (CAM_SE)
class NCHHomme : public NCHelper
{
public:
  NCHHomme(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}

  static bool can_read_file(ReadNC* readNC, const FileOptions& opts);

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);

  virtual ErrorCode create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads)
  { return _readNC->create_ucd_verts_quads(opts, file_set, quads); }

  virtual ErrorCode init_local_gid(EntityHandle file_set)
  { return _readNC->init_local_gid(file_set); }

  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
  { return _readNC->read_variables(file_set, var_names, tstep_nums, false); }

  virtual std::string get_mesh_type_name()
  { return "CAM_SE"; }
};

//! Child helper class for unknown CAM grid (CAM_UNKNOWN)
class NCHUnknownCam : public NCHelper
{
public:
  NCHUnknownCam(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);

  virtual ErrorCode create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads)
  { return MB_FAILURE; }

  virtual ErrorCode init_local_gid(EntityHandle file_set)
  { return MB_FAILURE; }

  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
  { return MB_FAILURE; }

  virtual std::string get_mesh_type_name()
  { return "CAM_UNKNOWN"; }
};

//! Child helper class for unknown grid (NOT_CAM)
class NCHNotCam : public NCHelper
{
public:
  NCHNotCam(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);

  virtual ErrorCode create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads)
  { return MB_FAILURE; }

  virtual ErrorCode init_local_gid(EntityHandle file_set)
  { return MB_FAILURE; }

  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums)
  { return MB_FAILURE; }

  virtual std::string get_mesh_type_name()
  { return "NOT_CAM"; }
};

} // namespace moab

#endif
