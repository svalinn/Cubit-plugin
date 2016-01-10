#ifndef MYEXPORTCOMMAND_HPP
#define MYEXPORTCOMMAND_HPP

#include "CubitCommandInterface.hpp"

class MeshExportInterface;

/*!
 * \brief The MyExportCommand class implements a sample export command
 * that also contains examples of how to use the mesh export interface.
 */
class MyExportCommand: public CubitCommand
{
public:
  MyExportCommand();
  ~MyExportCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

protected:
  //! \brief open a file
  //! \param[in] file the name of the file to opened
  //! \param[out] output_file the file stream that is used for writing
  //! \return true if the file open succeeded
  bool open_file(const std::string &file, std::ofstream& output_file);

  //! \brief write data to the file
  //! \param output_file the file stream that is used for writing
  //! \param[in] iface the reference to the MeshExportInterface
  //! \return true if file write succeeded
  bool write_file(std::ofstream& output_file, MeshExportInterface *iface);

  //! \brief close the file
  //! \param output_file the file stream that is used for writing
  //! \returns always returns true
  bool close_file(std::ofstream& output_file);

  //! \brief method used internally to write_file to write the mesh connectivity
  //! \param output_file the file stream that is used for writing
  //! \param[in] iface the reference to the MeshExportInterface
  //! \return true if file write succeeded
  bool write_connectivity(std::ofstream& output_file,MeshExportInterface *iface);
};

#endif // MYEXPORTCOMMAND_HPP
