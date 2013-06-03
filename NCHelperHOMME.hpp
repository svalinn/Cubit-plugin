//-------------------------------------------------------------------------
// Filename      : NCHelperHOMME.hpp
//
// Purpose       : Climate NC file helper for HOMME grid
//
// Creator       : Danqing Wu
//-------------------------------------------------------------------------

#ifndef NCHELPERHOMME_HPP
#define NCHELPERHOMME_HPP

#include "NCHelper.hpp"

namespace moab {

//! Child helper class for HOMME grid (CAM_SE)
class NCHelperHOMME : public NCHelper
{
public:
  NCHelperHOMME(ReadNC* readNC, int fileId, int spectralOrder) : NCHelper(readNC, fileId), _spectralOrder(spectralOrder) {}
  static bool can_read_file(ReadNC* readNC, const FileOptions& opts, int& spectralOrder);

private:
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set);
  virtual ErrorCode create_verts_quads(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& quads);
  virtual std::string get_mesh_type_name() { return "CAM_SE"; }
  virtual bool is_scd_mesh() { return false; }

  ErrorCode init_HOMMEucd_vals();
  ErrorCode create_HOMMEucd_verts_quads(const FileOptions& opts, EntityHandle file_set, Range& quads);

private:
  int _spectralOrder; // read from variable 'np'
};

} // namespace moab

#endif
