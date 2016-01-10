#include "MyPlugin.hpp"
#include "MyVersionCommand.hpp"
#include "MyExportCommand.hpp"

MyPlugin::MyPlugin()
{}

MyPlugin::~MyPlugin()
{}

std::vector<std::string> MyPlugin::get_keys()
{
  std::vector<std::string> keys;
  keys.push_back("MyVersionCommand");
  keys.push_back("MyExportCommand");

  return keys;
}

CubitCommand* MyPlugin::create_command(const std::string &key)
{
  // NOTE: The internals of Trelis will take owernship of the command,
  // and delete it when it is time to clean up.

  if(key == "MyVersionCommand")
    return new MyVersionCommand();

  else if(key == "MyExportCommand")
    return new MyExportCommand();

  return NULL;
}
