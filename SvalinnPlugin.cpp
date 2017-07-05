#include "SvalinnPlugin.hpp"
#include "CubitInterface.hpp"

//Not including this when you aren't building the DAGMC export means you don't need to have MakeWatertight.hpp
#ifdef BUILD_DAGMC_EXPORT
#include "export_dagmc_cmd/DAGMCExportCommand.hpp"
#endif

#ifdef BUILD_MCNP_IMPORT
#include "MCNPImp.hpp"
#endif

#ifdef BUILD_IGEOM_TEST
#include "iGeom_test.hpp"
#endif

SvalinnPlugin::SvalinnPlugin()
{

  CubitMessageHandler *console = CubitInterface::get_cubit_message_handler();
  if (console) {
    std::ostringstream load_message;
    load_message.str("");
    load_message << "Loaded Svalinn plugin." << std::endl;
    CubitInterface::get_cubit_message_handler()->print_error(load_message.str().c_str());
  }

}

SvalinnPlugin::~SvalinnPlugin()
{}

std::vector<std::string> SvalinnPlugin::get_keys()
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

CubitCommand* SvalinnPlugin::create_command(const std::string &key)
{
  // NOTE: The internals of Trelis will take owernship of the command,
  // and delete it when it is time to clean up.
  
  // Trelis will crash if included keys haven't been built.

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
