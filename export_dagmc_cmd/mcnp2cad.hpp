#ifndef MCNP2CAD_HPP
#define MCNP2CAD_HPP

//MCNP2CAD includes

#include "CubitCommandInterface.hpp"
#include "CubitMessageHandler.hpp"

// CGM includes
#include "moab/Interface.hpp"
#include "RefEntity.hpp"
#include <iostream>
#include <fstream>

extern std::ofstream record;

class MCNP2CAD: public CubitCommand
{
public:
  MCNP2CAD();
  ~MCNP2CAD();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute( CubitCommandData &data );
  bool parse_options( CubitCommandData &data/*, moab::EntityHandle* file_set*/ );
  class GeometryContext;
/*
protected:

  void teardown();

*/
//private:

 // CubitMessageHandler* console;

  std::ostringstream message;

  bool din, dout;

};


#endif // MCNP2CAD_HPP
