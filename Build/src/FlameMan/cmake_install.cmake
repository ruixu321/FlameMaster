# Install script for directory: /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/FlameMan

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xFlameManx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/bin" TYPE EXECUTABLE FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/FlameMan")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/bin/FlameMan" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/bin/FlameMan")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/../Bin/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/bin/FlameMan")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/bin/FlameMan")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xFlameManx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/cmake/FlameManTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/cmake/FlameManTargets.cmake"
         "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/CMakeFiles/Export/Bin/cmake/FlameManTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/cmake/FlameManTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/Bin/cmake/FlameManTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/cmake" TYPE FILE FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/CMakeFiles/Export/Bin/cmake/FlameManTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/cmake" TYPE FILE FILES "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/CMakeFiles/Export/Bin/cmake/FlameManTargets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xFlameManx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/Bin/cmake" TYPE FILE FILES
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/install_config/FlameMan/FlameManConfig.cmake"
    "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/lib/cmake/FlameMan/FlameManConfigVersion.cmake"
    )
endif()

