#include "MyExportCommand.hpp"
#include "MeshExportInterface.hpp"
#include "CubitInterface.hpp"
#include "CubitMessageHandler.hpp"

#include <fstream>

MyExportCommand::MyExportCommand()
{}

MyExportCommand::~MyExportCommand()
{}

std::vector<std::string> MyExportCommand::get_syntax()
{
  // Define the syntax for the command. Note the syntax is a modified BNF
  // format. Full documentation on the command specification syntax can be
  // found in the documentation.
  std::string syntax =
      "export sample <string:label='filename',help='<filename>'> "
      "[overwrite]";

  std::vector<std::string> syntax_list;
  syntax_list.push_back(syntax);

  return syntax_list;
}

std::vector<std::string> MyExportCommand::get_syntax_help()
{
  std::vector<std::string> help;
  return help;
}

std::vector<std::string> MyExportCommand::get_help()
{
  std::vector<std::string> help;
  return help;
}

bool MyExportCommand::execute(CubitCommandData &data)
{
  std::ofstream output_file;
  MeshExportInterface *iface;

  // The command syntax specified that we should have a string labeled "filename" in the
  // the command. We would not have gotten to this point without having that required
  // piece of data. Retrieve it here
  std::string filename;
  data.get_string("filename", filename);

  // The syntax also indicated that we could have an optional piece of data asking if we
  // wanted to overwrite an existing file. This will return true if the overwrite option
  // was specified in the command
  bool overwrite = data.find_keyword("overwrite");

  // TODO: check if the given file exists. If it does and 'overwrite' was not specified,
  // print an error and return from this function.

  // Get the MeshExport interface from CubitInterface
  iface = dynamic_cast<MeshExportInterface*>(CubitInterface::get_interface("MeshExport"));
  if(!iface)
  {
    CubitMessageHandler* console = CubitInterface::get_cubit_message_handler();
    console->print_error("Unable to get mesh export interface.\n");
    return false;
  }

  // Now we can use the interface as before
  iface->set_use_sequential_ids(false);

  // Open the file and write to it
  bool result =
      open_file(filename, output_file) && write_file(output_file, iface);

  // Close the file
  close_file(output_file);

  // Make sure that you release the interface after accessing it
  CubitInterface::release_interface(iface);

  return result;
}

bool MyExportCommand::open_file(const std::string& file, std::ofstream& output_file)
{
  output_file.open(file.c_str());

  return true;
}

bool MyExportCommand::write_file(std::ofstream& output_file, MeshExportInterface *iface)
{
  // Initialize the exporter
  iface->initialize_export();

  // Fetch and output the coordinates
  int number_nodes = iface->get_num_nodes();

  if(!number_nodes)
  {
    CubitMessageHandler* console = CubitInterface::get_cubit_message_handler();
    console->print_message("WARNING: No nodes in model...\n");
    return false;
  }

  output_file << "Test Output for MeshExportInterface\n";
  output_file << "Number of Nodes:  " << number_nodes << "\n";
  output_file << "List of Nodes:\n";

  int buf_size = 100;
  std::vector<double> xcoords(buf_size), ycoords(buf_size), zcoords(buf_size);
  std::vector<int> node_ids(buf_size);

  int start = 0, number_found = 0;
  while ((number_found = iface->get_coords(start, buf_size, xcoords, ycoords, zcoords, node_ids)))
  {
    // Write out the coordinates
    for(int i = 0; i < number_found; i++)
    {
      output_file << node_ids[i] << "  " <<
                     xcoords[i] << "  " <<
                     ycoords[i] << "  " <<
                     zcoords[i] << "\n";
    }
    start += number_found;
  }

  // Write the connectivity
  bool result = write_connectivity(output_file, iface);

  return result;
}

bool MyExportCommand::close_file(std::ofstream& output_file)
{
  output_file.close();
  return true;
}

bool MyExportCommand::write_connectivity(std::ofstream& output_file,MeshExportInterface *iface)
{
  output_file << " List of Elements" << std::endl;
  output_file << " Element Type, ID, Block, Connectivity" << std::endl;

  // Get the list of blocks
  std::vector<BlockHandle> blocks;
  iface->get_block_list(blocks);

  // Test for get_block_size
  int num_elements = 0;
  for (size_t i = 0; i < blocks.size(); i++)
  {
    int nelems;
    bool status = iface->get_block_size(blocks[i], nelems);
    if (status)
      num_elements += nelems;
  }

  // Get a batch of elements in an initialized block
  int buf_size = 100;
  std::vector<ElementType>   element_type(buf_size);
  std::vector<ElementHandle> handles(buf_size);

  // Elements in a buffer set will be of the same element type and in the same block
  for (size_t i = 0; i < blocks.size(); i++)
  {
    int num_elems;
    int start_index = 0;
    BlockHandle block = blocks[i];
    while( (num_elems = iface->get_block_elements(start_index, buf_size, block, element_type, handles)) > 0)
    {
      // Get ids for the element handles
      std::vector<int> ids(num_elems);
      iface->get_element_ids(num_elems, handles, ids);

      int block_id = iface->id_from_handle(block);

      // Write out the connectivity
      for (int i = 0; i < num_elems; i++)
      {
        std::vector<int> conn(27);
        int num_nodes = iface->get_connectivity(handles[i], conn);

        output_file << (int) element_type[i] << " " << ids[i] << " " << block_id << " ";
        for (int j = 0; j < num_nodes; j++)
          output_file << conn[j] << " ";
        output_file << std::endl;
      }
      start_index += num_elems;
    }
  }

  return true;
}
