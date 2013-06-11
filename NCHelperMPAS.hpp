//-------------------------------------------------------------------------
// Filename      : NCHelperMPAS.hpp
//
// Purpose       : Climate NC file helper for MPAS grid
//
// Creator       : Danqing Wu
//-------------------------------------------------------------------------

#ifndef NCHELPERMPAS_HPP
#define NCHELPERMPAS_HPP

#include "NCHelper.hpp"

namespace moab {

//! Child helper class for MPAS grid
class NCHelperMPAS : public NCHelper
{
public:
  NCHelperMPAS(ReadNC* readNC, int fileId, const FileOptions& opts);
  static bool can_read_file(ReadNC* readNC, int fileId);

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);
  virtual ErrorCode create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& faces);
  virtual std::string get_mesh_type_name() { return "MPAS"; }
  virtual bool is_scd_mesh() { return false; }
};

} // namespace moab

#endif
