//-------------------------------------------------------------------------
// Filename      : NCHelperFV.hpp
//
// Purpose       : Climate NC file helper for Finite Volume grid
//
// Creator       : Danqing Wu
//-------------------------------------------------------------------------

#ifndef NCHELPERFV_HPP
#define NCHELPERFV_HPP

#include "NCHelper.hpp"

namespace moab {

//! Child helper class for Finite Volume grid (CAM_FV)
class NCHelperFV : public NCHelper
{
public:
  NCHelperFV(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}
  static bool can_read_file(ReadNC* readNC, int fileId);

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);
  virtual ErrorCode create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads);
  virtual std::string get_mesh_type_name() { return "CAM_FV"; }
  virtual bool is_scd_mesh() { return true; }
};

} // namespace moab

#endif
