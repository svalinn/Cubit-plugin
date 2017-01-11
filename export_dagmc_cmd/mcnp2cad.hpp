#ifndef MCNP2CAD_HPP
#define MCNP2CAD_HPP

//MCNP2CAD includes

#include "CubitCommandInterface.hpp"
#include "CubitMessageHandler.hpp"

// CGM includes
#include "moab/Interface.hpp"
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
  moab::ErrorCode parse_options( CubitCommandData &data, moab::EntityHandle* file_segt );
  class GeometryContext;
/*
protected:

  void teardown();

*/
//private:

 // CubitMessageHandler* console;

  std::ostringstream message;

  static double specific_tol;
  bool verbose, debug, din, dout, extraEff, skipMats, skipMerge,
       skipImps, skipNums, skipGrave, skipImprint, uwuwNames;

};

#endif // MCNP2CAD_HPP
