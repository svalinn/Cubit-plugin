#ifndef MCNP2CAD_HPP
#define MCNP2CAD_HPP

//MCNP2CAD includes

#include "CubitCommandInterface.hpp"
#include "CubitMessageHandler.hpp"

// CGM includes
#include "RefEntity.hpp"

class MCNP2CAD: public CubitCommand
{
public:
  MCNP2CAD();
  ~MCNP2CAD();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute( CubitCommandData &data );
  void parse_options( CubitCommandData &data );
/*
protected:

  void teardown();

*/
private:

 // CubitMessageHandler* console;

  std::ostringstream message;

/*
  int norm_tol;
  double faceting_tol;
  double len_tol;
  bool verbose_warnings;
  bool fatal_on_curves;
  */

};

#endif // MCNP2CAD_HPP
