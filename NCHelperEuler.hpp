//-------------------------------------------------------------------------
// Filename      : NCHelperEuler.hpp
//
// Purpose       : Climate NC file helper for Eulerian Spectral grid
//
// Creator       : Danqing Wu
//-------------------------------------------------------------------------

#ifndef NCHELPEREULER_HPP
#define NCHELPEREULER_HPP

#include "NCHelper.hpp"

namespace moab {

//! Child helper class for Eulerian Spectral grid (CAM_EUL)
class NCHelperEuler : public NCHelper
{
public:
  NCHelperEuler(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}

  static bool can_read_file(ReadNC* readNC, int fileId);

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);

  virtual ErrorCode create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads);

  virtual std::string get_mesh_type_name() { return "CAM_EUL"; }

  virtual bool is_scd_mesh() { return true; }
};

} // namespace moab

#endif
