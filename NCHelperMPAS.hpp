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
class NCHelperMPAS : public UcdNCHelper
{
public:
  NCHelperMPAS(ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet);
  static bool can_read_file(ReadNC* readNC);

private:
  //! Implementation of NCHelper::init_mesh_vals()
  virtual ErrorCode init_mesh_vals();
  //! Implementation of NCHelper::check_existing_mesh()
  virtual ErrorCode check_existing_mesh();
  //! Implementation of NCHelper::create_mesh()
  virtual ErrorCode create_mesh(Range& quads);
  //! Implementation of NCHelper::get_mesh_type_name()
  virtual std::string get_mesh_type_name() { return "MPAS"; }

  //! Implementation of UcdNCHelper::read_ucd_variable_to_nonset_allocate()
  virtual ErrorCode read_ucd_variable_to_nonset_allocate(std::vector<ReadNC::VarData>& vdatas,
                                                         std::vector<int>& tstep_nums);
  //! Implementation of UcdNCHelper::read_ucd_variable_setup()
  virtual ErrorCode read_ucd_variable_setup(std::vector<std::string>& var_names,
                                            std::vector<int>& tstep_nums,
                                            std::vector<ReadNC::VarData>& vdatas,
                                            std::vector<ReadNC::VarData>& vsetdatas);
#ifdef PNETCDF_FILE
  //! Implementation of UcdNCHelper::read_ucd_variable_to_nonset_async()
  virtual ErrorCode read_ucd_variable_to_nonset_async(std::vector<ReadNC::VarData>& vdatas,
                                                      std::vector<int>& tstep_nums);
#else
  //! Implementation of UcdNCHelper::read_ucd_variable_to_nonset()
  virtual ErrorCode read_ucd_variable_to_nonset(std::vector<ReadNC::VarData>& vdatas,
                                                std::vector<int>& tstep_nums);
#endif

private:
  int maxEdgesPerCell;
  int numCellGroups;
  std::map<EntityHandle, int> cellHandleToGlobalID;
  Range facesOwned;
};

} // namespace moab

#endif
