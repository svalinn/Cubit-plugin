#ifndef MYVERSIONCOMMAND_HPP
#define MYVERSIONCOMMAND_HPP

#include "CubitCommandInterface.hpp"

/*!
 * \brief The MyVersionCommand class reimplements the "version" command to
 * display information about the version of MyPlugin in addition to the
 * version information displayed by Trelis.
 */
class MyVersionCommand : public CubitCommand
{
public:
  MyVersionCommand();
  ~MyVersionCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);
};

#endif // MYVERSIONCOMMAND_HPP
