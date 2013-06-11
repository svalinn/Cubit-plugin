#include "NCHelper.hpp"
#include "NCHelperEuler.hpp"
#include "NCHelperFV.hpp"
#include "NCHelperHOMME.hpp"
#include "NCHelperMPAS.hpp"

namespace moab {

NCHelper* NCHelper::get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts)
{
  // Check if CF convention is being followed
  bool is_CF = false;

  std::map<std::string, ReadNC::AttData>& globalAtts = readNC->globalAtts;
  std::map<std::string, ReadNC::AttData>::iterator attIt = globalAtts.find("conventions");
  if (attIt == globalAtts.end())
    attIt = globalAtts.find("Conventions");

  if (attIt != globalAtts.end()) {
    unsigned int sz = attIt->second.attLen;
    std::string att_data;
    att_data.resize(sz + 1);
    att_data[sz] = '\000';
    int success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0]);
    if (0 == success && att_data.find("CF") != std::string::npos)
      is_CF = true;
  }

  if (is_CF) {
    if (NCHelperEuler::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperEuler(readNC, fileId);
    else if (NCHelperFV::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperFV(readNC, fileId);
    else if (NCHelperHOMME::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperHOMME(readNC, fileId, opts);
  }
  else {
    if (NCHelperMPAS::can_read_file(readNC, fileId))
      return new (std::nothrow) NCHelperMPAS(readNC, fileId, opts);
  }

  // Unknown NetCDF grid (will fill this in later for POP, CICE and CLM)
  return NULL;
}

} // namespace moab
