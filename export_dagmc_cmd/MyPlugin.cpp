#include "MyPlugin.hpp"
#include "DAGMCExportCommand.hpp"
#include "iGeom_test.hpp"
#include "MCNPImp.hpp"

MyPlugin::MyPlugin()
{}

MyPlugin::~MyPlugin()
{}

std::vector<std::string> MyPlugin::get_keys()
{
  std::vector<std::string> keys;
#ifdef BUILD_DAGMC_EXPORT
  keys.push_back("DAGMCExportCommand");
#endif
#ifdef BUILD_IGEOM_TESTS
  keys.push_back("iGeom_test");
#endif
#ifdef BUILD_MCNP_IMPORT
  keys.push_back("MCNPImp");
#endif

  return keys;
}

CubitCommand* MyPlugin::create_command(const std::string &key)
{
  // NOTE: The internals of Trelis will take owernship of the command,
  // and delete it when it is time to clean up.

#ifdef BUILD_DAGMC_EXPORT
  if(key == "DAGMCExportCommand")
    return new DAGMCExportCommand();
#endif

#ifdef BUILD_IGEOM_TESTS
  if(key == "iGeom_test")
    return new iGeom_test();
#endif
  
#ifdef BUILD_MCNP_IMPORT
  if(key == "MCNPImp")
    return new MCNPImp();
#endif

  return NULL;
}
