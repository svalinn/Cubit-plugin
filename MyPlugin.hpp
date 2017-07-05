#ifndef MYPLUGIN_HPP
#define MYPLUGIN_HPP

#include "CubitCommandInterface.hpp"
#include "CubitPluginExport.hpp"

/*!
 * \brief The MyPlugin class derives from CubitCommandInterface
 * which gives access to the parser and allows you to implement
 * your own commands.
 */
class MyPlugin : public CubitCommandInterface
{
public:
  MyPlugin();
  ~MyPlugin();

  //! Used to get a list of unique keys that identify the commands this plugin
  //! will create.
  std::vector<std::string> get_keys();

  //! Returns the command corresponding to 'key'.
  CubitCommand* create_command(const std::string &key);
};

//! This macro is required to identify this as a valid Cubit plugin. The plugin
//! will NOT be loaded if this macro is not present.
CUBIT_PLUGIN(MyPlugin)

#endif // MYPLUGIN_HPP
