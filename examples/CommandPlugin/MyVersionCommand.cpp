#include "MyVersionCommand.hpp"
#include "CubitInterface.hpp"
#include "CubitMessageHandler.hpp"

MyVersionCommand::MyVersionCommand()
{}

MyVersionCommand::~MyVersionCommand()
{}

std::vector<std::string> MyVersionCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back("version");

  return syntax_list;
}

std::vector<std::string> MyVersionCommand::get_syntax_help()
{
  std::vector<std::string> help;
  return help;
}

std::vector<std::string> MyVersionCommand::get_help()
{
  std::vector<std::string> help;
  return help;
}

bool MyVersionCommand::execute(CubitCommandData &data)
{
  std::string cubit_version = CubitInterface::get_version();
  std::string build_number = CubitInterface::get_build_number();
  std::string revision_date = CubitInterface::get_revision_date();
  std::string geom_version = CubitInterface::get_acis_version();
  std::string exodus_version = CubitInterface::get_exodus_version();
  std::string graphics_version = CubitInterface::get_graphics_version();

  CubitMessageHandler* console = CubitInterface::get_cubit_message_handler();

  std::string output;

  output = "\tMyPlugin Version 0.1\n";
  console->print_message(output.c_str());

  output = "\tCubit Version " + cubit_version + " Build " + build_number + "\n";
  console->print_message(output.c_str());

  output = "\t" + graphics_version + "\n";
  console->print_message(output.c_str());

  output = "\t" + geom_version + "\n";
  console->print_message(output.c_str());

  output = "\tExodus Version " + exodus_version + "\n";
  console->print_message(output.c_str());

  output = "\tRevised " + revision_date + "\n\n",
  console->print_message(output.c_str());

  return true;
}
