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

  static NCHelper* get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts);

  //! Interfaces to be implemented by child classes
  virtual ErrorCode init_mesh_vals(const FileOptions& opts, EntityHandle file_set) = 0;
  virtual ErrorCode create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& faces) = 0;
  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums) = 0;
  virtual std::string get_mesh_type_name() = 0;

protected:
  ReadNC* _readNC;
  int _fileId;
};

//! Child helper class for structured mesh, e.g. CAM_EL or CAM_FV
class ScdNCHelper : public NCHelper
{
public:
  ScdNCHelper(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}
  virtual ~ScdNCHelper() {}

private:
  //! Implementation of NCHelper::create_mesh()
  virtual ErrorCode create_mesh(ScdInterface* scdi, const FileOptions& opts, EntityHandle file_set, Range& faces);
  //! Implementation of NCHelper::read_variables()
  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names, std::vector<int>& tstep_nums);

  //! These functions are used by read_variables(), fully implemented
  ErrorCode read_scd_variable_setup(std::vector<std::string>& var_names,
                                    std::vector<int>& tstep_nums,
                                    std::vector<ReadNC::VarData>& vdatas,
                                    std::vector<ReadNC::VarData>& vsetdatas);
  ErrorCode read_scd_variable_to_set_allocate(std::vector<ReadNC::VarData>& vdatas,
                                              std::vector<int>& tstep_nums);
  ErrorCode read_scd_variable_to_set(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                     std::vector<int>& tstep_nums);
  ErrorCode read_scd_variable_to_nonset_allocate(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                                 std::vector<int>& tstep_nums);
  ErrorCode read_scd_variable_to_nonset(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                        std::vector<int>& tstep_nums);
  ErrorCode convert_scd_variable(ReadNC::VarData& var_data, int tstep_num);
};

//! Child helper class for unstructured mesh, e.g. CAM_SE (HOMME) or MPAS
class UcdNCHelper : public NCHelper
{
public:
  UcdNCHelper(ReadNC* readNC, int fileId) : NCHelper(readNC, fileId) {}
  virtual ~UcdNCHelper() {}

private:
  //! Implementation of NCHelper::read_variables()
  virtual ErrorCode read_variables(EntityHandle file_set, std::vector<std::string>& var_names,
                                   std::vector<int> &tstep_nums);

  //! These functions are used by read_variables(), partially implemented
  virtual ErrorCode read_ucd_variable_setup(std::vector<std::string>& var_names,
                                            std::vector<int>& tstep_nums,
                                            std::vector<ReadNC::VarData>& vdatas,
                                            std::vector<ReadNC::VarData>& vsetdatas) = 0;
  ErrorCode read_ucd_variable_to_set_allocate(std::vector<ReadNC::VarData>& vdatas,
                                              std::vector<int>& tstep_nums);
  ErrorCode read_ucd_variable_to_set(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                     std::vector<int>& tstep_nums);
  virtual ErrorCode read_ucd_variable_to_nonset_allocate(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                                         std::vector<int>& tstep_nums) = 0;
#ifdef PNETCDF_FILE
  virtual ErrorCode read_ucd_variable_to_nonset_async(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                                      std::vector<int>& tstep_nums) = 0;
#else
  virtual ErrorCode read_ucd_variable_to_nonset(EntityHandle file_set, std::vector<ReadNC::VarData>& vdatas,
                                                std::vector<int>& tstep_nums) = 0;
#endif
  virtual ErrorCode convert_ucd_variable(ReadNC::VarData& var_data, int tstep_num) = 0;
};

} // namespace moab

#endif
