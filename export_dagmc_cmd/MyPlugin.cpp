#include "MyPlugin.hpp"
#include "DAGMCExportCommand.hpp"
#include "iGeom_test.hpp"
#include "mcnp2cad.hpp"

MyPlugin::MyPlugin()
{}

MyPlugin::~MyPlugin()
{}

std::vector<std::string> MyPlugin::get_keys()
{
  std::vector<std::string> keys;
  keys.push_back("DAGMCExportCommand");
  keys.push_back("iGeom_test");
  keys.push_back("MCNP2CAD");

  return keys;
}

CubitCommand* MyPlugin::create_command(const std::string &key)
{
  // NOTE: The internals of Trelis will take owernship of the command,
  // and delete it when it is time to clean up.

  if(key == "DAGMCExportCommand")
    return new DAGMCExportCommand();

  if(key == "iGeom_test")
    return new iGeom_test();
  
  if(key == "MCNP2CAD")
    return new MCNP2CAD();

  return NULL;
}
