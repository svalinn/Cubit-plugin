#include "NCHelper.hpp"
#include "NCHelperEuler.hpp"
#include "NCHelperFV.hpp"
#include "NCHelperHOMME.hpp"
#include "moab/ReadUtilIface.hpp"

namespace moab {

NCHelper* NCHelper::get_nc_helper(ReadNC* readNC, int fileId, const FileOptions& opts)
{
  std::map<std::string, ReadNC::AttData>::iterator attIt = readNC->globalAtts.find("source");
  if (attIt == readNC->globalAtts.end()) {
    readNC->readMeshIface->report_error("%s", "File does not have source global attribute.");
    return NULL;
  }

  unsigned int sz = attIt->second.attLen;
  std::string att_data;
  att_data.resize(sz + 1);
  att_data[sz] = '\000';
  int success = NCFUNC(get_att_text)(fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0]);
  if (success != 0) {
	  readNC->readMeshIface->report_error("%s", "Failed to read source global attribute char data.");
    return NULL;
  }

  // If a CAM file, which type?
  if (att_data.find("CAM") != std::string::npos) {
    int spectralOrder = -1;

    if (NCHelperEuler::can_read_file(readNC))
      return new (std::nothrow) NCHelperEuler(readNC, fileId);
    else if (NCHelperFV::can_read_file(readNC))
      return new (std::nothrow) NCHelperFV(readNC, fileId);
    else if (NCHelperHOMME::can_read_file(readNC, opts, spectralOrder))
      return new (std::nothrow) NCHelperHOMME(readNC, fileId, spectralOrder);
    else // Unknown CAM type
      return NULL;
  }

  // Not a CAM file, will fill this in later for POP, CICE and CLM
  return NULL;
}

} // namespace moab
