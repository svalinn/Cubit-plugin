/*!
 * \file   MyObserver.hpp
 * \brief
 *   Provides an example of how to listen to events from Cubit.
 * \author Michael Plooster
 * \date   4 Aug 2014
 */

#ifndef MYOBSERVER_HPP
#define MYOBSERVER_HPP

#include "CubitInterface.hpp"

// The MyObserver class re-implements methods that we are interested in.
// These methods are called by Cubit in response to events. For example,
// when a user saves their model to a .cub file, a save event is generated
// and the method "notify_cub_file_saved" is executed.
class MyObserver : public CubitInterface::ObserverBase
{
public:
  MyObserver();
  ~MyObserver();

  // Provides an example of reading our custom data from a .cub file.
  // This function is executed after Cubit is done reading the file.
  void notify_cub_file_read(std::string filename);

  // Provides an example of writing our custom data to a .cub file.
  // This function is executed after Cubit is done writing its data
  // to the file.
  void notify_cub_file_saved(std::string filename);
};

#endif // MYOBSERVER_HPP
