# Install script for directory: /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xSourceScriptsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/bin" TYPE FILE PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ FILES
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/install_config/FlameMaster/Source.bash"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/install_config/FlameMaster/Source.zsh"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/install_config/FlameMaster/Source.csh"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/install_config/FlameMaster/Source.BAT"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xFlameMaster_Examplesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Run" TYPE DIRECTORY FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/examples/" USE_SOURCE_PERMISSIONS)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xFlameMaster_Documentationx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/lib/Doc" TYPE DIRECTORY FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/doc/")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/cmake_install.cmake")
  include("/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/cmake_install.cmake")
  include("/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/ScanMan/cmake_install.cmake")
  include("/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
