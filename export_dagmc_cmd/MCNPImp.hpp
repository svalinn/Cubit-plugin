#ifndef MCNPIMP_HPP
#define MCNPIMP_HPP

#include "CubitCommandInterface.hpp"
#include "CubitMessageHandler.hpp"

// CGM includes
#include "moab/Interface.hpp"
#include "RefEntity.hpp"
#include <iostream>
#include <fstream>

class MCNPImp: public CubitCommand
{
public:
  MCNPImp();
  ~MCNPImp();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute( CubitCommandData &data );
  bool parse_options( CubitCommandData &data );
};

#endif //MCNPIMP_HPP
