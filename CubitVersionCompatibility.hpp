#ifndef SVALINN_COMPATIBILITY_HPP
#define SVALINN_COMPATIBILITY_HPP

#include "CubitVersion.h"


#if CUBIT_VERSION_MAJOR <= 17
    //#include "CubitInterface.hpp"
    #define CUBIT_INTERFACE_HEADER "CubitInterface.hpp"
    #define MSG_HANDLER CubitInterface::get_cubit_message_handler()
#else
    //#include "CubitCoreformInterface.hpp"
    #define CUBIT_INTERFACE_HEADER "CubitCoreformInterface.hpp"
    #define MSG_HANDLER CubitMessage::get_message_handler()
#endif

#endif  // SVALINN_COMPATIBILITY_HPP
