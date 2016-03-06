#include "MyPlugin.hpp"
#include "DAGMCExportCommand.hpp"

MyPlugin::MyPlugin()
{}

MyPlugin::~MyPlugin()
{}

std::vector<std::string> MyPlugin::get_keys()
{
  std::vector<std::string> keys;
  keys.push_back("DAGMCExportCommand");

  return keys;
}

CubitCommand* MyPlugin::create_command(const std::string &key)
{
  // NOTE: The internals of Trelis will take owernship of the command,
  // and delete it when it is time to clean up.

  if(key == "DAGMCExportCommand")
    return new DAGMCExportCommand();

  return NULL;
}
