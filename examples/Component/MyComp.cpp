/*!
 *  \file MyComp.cpp
 *  \brief
 *    Provides an example of creating a custom component.
 *  \author Michael Plooster
 *  \date 11 Feb 2014
 */

#include "MyComp.hpp"
#include "MenuManager.hpp"
#include "ToolbarManager.hpp"
#include "PanelManager.hpp"
#include "ExportManager.hpp"
#include "MyObserver.hpp"

// Default constructor. Remember to include the component name (should match
// the module name in mycomp.i).
MyComp::MyComp() :
  Component("mycomp"),
  myMenus(NULL),
  myToolbars(NULL),
  myPanels(NULL),
  myExportManager(NULL),
  mListener(NULL)
{}

MyComp::~MyComp()
{
  if(myMenus)
    delete myMenus;

  if(myToolbars)
    delete myToolbars;

  if(myPanels)
    delete myPanels;

  if(myExportManager)
    delete myExportManager;

  if(mListener)
    delete mListener;
}

void MyComp::start_up(int withGUI)
{
  if(withGUI)
  {
    setup_menus();
    setup_toolbars();
    setup_command_panels();
    add_exports();
  }

  setup_observers(withGUI);
}

void MyComp::clean_up()
{
  cleanup_menus();
  cleanup_toolbars();
  cleanup_command_panels();
  cleanup_exports();
  cleanup_observers();

  // Let the framework know you are done.
  clean_up_complete();
}

void MyComp::setup_menus()
{
  if(!myMenus)
    myMenus = new MenuManager;

  myMenus->add_to_existing_menu();
  myMenus->setup_custom_menu();
}

void MyComp::cleanup_menus()
{
  if(myMenus)
    myMenus->remove_menu_items();
}

void MyComp::setup_toolbars()
{
  if(!myToolbars)
    myToolbars = new ToolbarManager;

  myToolbars->add_to_existing_toolbar();
  myToolbars->setup_custom_toolbar();
}

void MyComp::cleanup_toolbars()
{
  if(myToolbars)
    myToolbars->remove_toolbar_items();
}

void MyComp::setup_command_panels()
{
  if(!myPanels)
    myPanels = new PanelManager;

  myPanels->initialize();
}

void MyComp::cleanup_command_panels()
{
  if(myPanels)
    myPanels->clear();
}

void MyComp::add_exports()
{
  if(!myExportManager)
  {
    myExportManager = new ExportManager();
    myExportManager->add_export_types();
  }
}

void MyComp::cleanup_exports()
{
  if(myExportManager)
    myExportManager->remove_export_types();
}

void MyComp::setup_observers(int withGUI)
{
  if(!mListener)
    mListener = new MyObserver();
}

void MyComp::cleanup_observers()
{
  if(mListener)
  {
    delete mListener;
    mListener = NULL;
  }
}
