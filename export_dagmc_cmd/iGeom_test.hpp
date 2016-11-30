#ifndef IGEOM_TEST_HPP
#define IGEOM_TEST_HPP

#include "CubitCommandInterface.hpp"
#include "CubitMessageHandler.hpp"

// CGM includes
#include "RefEntity.hpp"

class iGeom_test: public CubitCommand
{
public:
  iGeom_test();
  ~iGeom_test();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

private:
  double radius;

};

#endif // MYEXPORTCOMMAND_HPP
