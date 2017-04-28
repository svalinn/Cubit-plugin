#ifndef MCNPIMP_HPP
#define MCNPIMP_HPP
#include <iostream>
#include <fstream>

#include "CubitCommandInterface.hpp"
#include "CubitMessageHandler.hpp"

// CGM includes
#include "moab/Interface.hpp"
#include "RefEntity.hpp"

class MCNPImp: public CubitCommand
{
public:
  MCNPImp();
  ~MCNPImp();

  // These three are needed by the plugin for the syntax of the command to work
  // get_syntax_help() and get_help() return empty vectors.
  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();

  // Analogous to the main function.
  bool execute(CubitCommandData &data);

  // Takes the data from Trelis and sets global variables for mcnp2cad options.
  void parse_options(CubitCommandData &data);

};

#endif //MCNPIMP_HPP
