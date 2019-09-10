# Install script for directory: /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/Tools/ListTool/ListTool

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xListToolx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/bin" TYPE EXECUTABLE FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/ListTool/ListTool")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/bin/ListTool" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/bin/ListTool")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/../Bin/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/bin/ListTool")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/bin/ListTool")
    endif()
  endif()
endif()

