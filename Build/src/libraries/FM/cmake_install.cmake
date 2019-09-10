# Install script for directory: /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/..")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xFlameMaster_Librariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/lib" TYPE STATIC_LIBRARY FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/FM/libFM.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/lib/libFM.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/lib/libFM.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/lib/libFM.a")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xFlameMaster_Librariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/include" TYPE FILE FILES
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM/Interrupt.h"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM/ListTool.h"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM/Spline.h"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM/MapMan.h"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM/TofZ.h"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM/BetaPDF.h"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM/SmallNewton.h"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM/paramtr.h"
    )
endif()

