#include "mcnp2cad/options.hpp"
#include "MCNPImp.hpp"
#include "mcnp2cad.hpp"

MCNPImp::MCNPImp()
{
  // set default values
  Gopt.verbose = false;
  Gopt.debug = false;
  Gopt.din = false;
  Gopt.dout = false;
  Gopt.infinite_lattice_extra_effort = false;
  Gopt.tag_materials = true;
  Gopt.merge_geom = true;
  Gopt.tag_importances = true;
  Gopt.tag_cell_IDs = true;
  Gopt.make_graveyard = true;
  Gopt.imprint_geom = true;
  Gopt.uwuw_names = false;
  Gopt.override_tolerance = false;
  Gopt.input_file = "";

}

MCNPImp::~MCNPImp()
{}

std::vector<std::string> MCNPImp::get_syntax()
{
  // Define the syntax for the command. Note the syntax is a modified BNF
  // format. Full documentation on the command specification syntax can be
  // found in the documentation.
  std::string syntax =
      "import MCNP "
      "<string:label='filename',help='<filename>'> "
      "[verbose] [debug] [debug_output] [debug_input] "
      "[extra_effort] [skip_mats] [skip_merge] "
      "[skip_imps] [skip_nums] [skip_graveyard] "
      "[skip_imprint] [uwuw_names] "
      "[tol <value:label='specific_tolerance',help='<specific tolerance>'>] ";

  std::vector<std::string> syntax_list;
  syntax_list.push_back(syntax);

  return syntax_list;
}

std::vector<std::string> MCNPImp::get_syntax_help()
{
  std::vector<std::string> help;
  return help;
}

std::vector<std::string> MCNPImp::get_help()
{
  std::vector<std::string> help;
  return help;
}

bool MCNPImp::execute(CubitCommandData &data)
{
  std::string filename;
  data.get_string("filename",filename);

  bool rval;

  rval = parse_options(data);
//  TODO: find a way to check whether the parsing worked as expected
//  CHK_MB_ERR_RET("Error parsing options: ",rval);

  return import_mcnp(filename);
}


bool MCNPImp::parse_options(CubitCommandData &data)
{

  // read parsed command for tolerance
  data.get_value("specific_tolerance",Gopt.specific_tolerance);
  record << "Setting specific tolerance to " << Gopt.specific_tolerance << std::endl;
  //If tolerance was specified
  if( data.find_keyword("tol") ){
    Gopt.override_tolerance = true;
    if( Gopt.specific_tolerance <= 0.0 || Gopt.specific_tolerance > .1 ){
      std::cerr << "Warning: you seem to have specified an unusual tolerance (" 
                << Gopt.specific_tolerance << ")." << std::endl;
      record << "Warning: you seem to have specified an unusual tolerance (" 
             << Gopt.specific_tolerance << ")." << std::endl;
    }
  }
  
  // read parsed boolean commands
  Gopt.verbose = data.find_keyword("verbose");
  Gopt.debug = data.find_keyword("debug");
  Gopt.infinite_lattice_extra_effort = data.find_keyword("extra_effort");
  Gopt.tag_materials = !data.find_keyword("skip_mats");
  Gopt.tag_importances = !data.find_keyword("skip_imps");
  Gopt.tag_cell_IDs = !data.find_keyword("skip_nums");
  Gopt.make_graveyard = !data.find_keyword("skip_graveyard");
  Gopt.merge_geom = !data.find_keyword("skip_merge");
  Gopt.imprint_geom = !data.find_keyword("skip_imprint");
  Gopt.uwuw_names = data.find_keyword("uwuw_names");
  Gopt.dout = data.find_keyword("debug_output");
  Gopt.din  = data.find_keyword("debug_input");

  if( Gopt.merge_geom && !Gopt.imprint_geom ) {
    record << "Warning: cannot merge geometry without imprinting, will skip merge too." << std::endl;
    std::cerr << "Warning: cannot merge geometry without imprinting, will skip merge too." << std::endl;
  }
  
  record << "Reading input file..." << std::endl;

  // if debug_input and not debug, set debugging to be true for InputDeck::build() call only
  
  if( Gopt.din && !OPT_DEBUG ){
    Gopt.debug = true;
  }
  else{ Gopt.din = false; }


  return true;
}
