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
  virtual ErrorCode create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& faces) = 0;
  virtual std::string get_mesh_type_name() = 0;
  virtual bool is_scd_mesh() = 0;

protected:
  ReadNC* _readNC;
  int _fileId;
};

} // namespace moab

#endif
