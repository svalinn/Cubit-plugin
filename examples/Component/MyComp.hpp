/*!
 *  \file MyComp.hpp
 *  \brief
 *    Provides an example of creating a custom component.
 *  \author Michael Plooster
 *  \date 11 Feb 2014
 */

#ifndef MYCOMP_HPP
#define MYCOMP_HPP

#include "Component.hpp"

class MenuManager;
class ToolbarManager;
class PanelManager;
class ExportManager;

class MyObserver;

class MyComp : public Component
{
public:
  MyComp();
  ~MyComp();

  // Methods inherited from Component
  void start_up(int withGUI);
  void clean_up();

private:
  void setup_menus();
  void cleanup_menus();

  void setup_toolbars();
  void cleanup_toolbars();

  void setup_command_panels();
  void cleanup_command_panels();

  //! Add export types to the file menu
  void add_exports();
  void cleanup_exports();

  void setup_observers(int withGUI);
  void cleanup_observers();

  MenuManager*    myMenus;
  ToolbarManager* myToolbars;
  PanelManager*   myPanels;
  ExportManager*  myExportManager;

  MyObserver* mListener;
};

#endif // MYCOMP_HPP
