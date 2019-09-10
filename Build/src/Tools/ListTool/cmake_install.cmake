# Install script for directory: /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/Tools/ListTool

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/cmake/ListToolTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/cmake/ListToolTargets.cmake"
         "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/CMakeFiles/Export/Bin/cmake/ListToolTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/cmake/ListToolTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/cmake/ListToolTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/cmake" TYPE FILE FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/CMakeFiles/Export/Bin/cmake/ListToolTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/cmake" TYPE FILE FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/CMakeFiles/Export/Bin/cmake/ListToolTargets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xListToolx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/cmake" TYPE FILE FILES
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/install_config/ListTool/ListToolConfig.cmake"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/lib/cmake/ListTool/ListToolConfigVersion.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/List/cmake_install.cmake")
  include("/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/FLReader/cmake_install.cmake")
  include("/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/ListTool/cmake_install.cmake")
  include("/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/WSS/cmake_install.cmake")

endif()

