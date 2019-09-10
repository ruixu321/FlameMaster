cmake_minimum_required(VERSION 2.8.12)
project(FlameMaster)
include(${CMAKE_SOURCE_DIR}/cmake/Macros.cmake)

# The version number.
set (FlameMaster_VERSION_MAJOR 4)
set (FlameMaster_VERSION_MINOR 0)
set (FlameMaster_VERSION_PATCH 0)
set(FlameMaster_VERSION "${FlameMaster_VERSION_MAJOR}.${FlameMaster_VERSION_MINOR}.${FlameMaster_VERSION_PATCH}")

#add_compile_options(-std=c++11)
########################################################################
#             GLOBAL OPTIONS FOR THE FLAMEMASTER PACKAGE               #
########################################################################

option(COMPILE_FORTRAN_SRC "Compile optional fortran code
like dassl, the adiabatic flame temperature module, and a mechanism converter tool" OFF)

option(COMPILE_SCANMANNEW "Compile the new ScanMan" OFF)

option(WRITE_LOCAL_LEWIS "Write Lewis numbers at every grid point" OFF)

option(PREM_UPWIND_CONVEC "Use upwind scheme for the convective term. OFF means central differences are used instead. This option is only used in premixed flames." ON)

option(BILIN_OMEGA "Compile with bilinear interpolation for collision integrals (omega)" OFF)

option(NEWTON_PERFORMANCE "Compile with time measurement of each Newton step" OFF)


########################################################################
#        DIRECTORIES/LOCATIONS USED IN THE FLAMEMASTER PACKAGE         #
########################################################################

# FlameMaster has always been installed locally
# The following lines ensure that this is still the case
SET(FLAMEMASTER_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/.." CACHE PATH "Package prefix prepended to
install directories")
SET(CMAKE_INSTALL_PREFIX "${FLAMEMASTER_INSTALL_PREFIX}" CACHE INTERNAL "Invisible prefix
prepended to install directories" FORCE)

if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

set(FlameMaster_DATA_PATH ${CMAKE_BINARY_DIR}/../Data
  CACHE PATH "Absolute path to directory used by ScanMan and FlameMan to find input")

if(IS_ABSOLUTE FlameMaster_DATA_PATH)
    message(FATAL_ERROR "The path FlameMaster_DATA_PATH has to be an absolute path!\n"
                        "Value: ${FlameMaster_DATA_PATH}")
endif(IS_ABSOLUTE FlameMaster_DATA_PATH)

if(NOT IS_DIRECTORY "${FlameMaster_DATA_PATH}")
    message(STATUS "Creating data directory: '${FlameMaster_DATA_PATH}'")
    file(MAKE_DIRECTORY ${FlameMaster_DATA_PATH})
endif(NOT IS_DIRECTORY "${FlameMaster_DATA_PATH}")

# follow conventions for cmake installation directories?
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR cmake)
else()
  set(DEF_INSTALL_CMAKE_DIR cmake) # normally this is 'lib/cmake', but
                                   # this makes it harder use find_package
                                   # with a simple HINT like
                                   # <prefix>/<to>/<some>/<FlameMaster>/<install>/<dir>
endif()
# Remove 'lib' prefix for shared libraries on Windows
if (WIN32)
  #set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBRARIES OFF)
if(MINGW)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif(MINGW)
  set(CMAKE_SHARED_LIBRARY_PREFIX "")
endif ()

# Used to modify the installation directories
option(USE_FLAMEMASTER_PREFIX "This is an option to create a nicer DragNDrop Installer for MacOS" OFF)
if(USE_FLAMEMASTER_PREFIX)
	set(extra_prefix "FlameMaster/")
else()
	set(extra_prefix "")
endif(USE_FLAMEMASTER_PREFIX)

set(INSTALL_LIB_DIR ${extra_prefix}Bin/lib CACHE PATH "Installation directory for libraries" FORCE)
set(INSTALL_BIN_DIR ${extra_prefix}Bin/bin CACHE PATH "Installation directory for executables" FORCE)
set(INSTALL_INCLUDE_DIR ${extra_prefix}Bin/include CACHE PATH "Installation directory for header files" FORCE)
set(INSTALL_EXAMPLES_DIR ${extra_prefix}Run CACHE PATH "Installation directory for CreateBinFile, ScanMan, and FlameMan examples" FORCE)
set(INSTALL_DOC_DIR ${extra_prefix}${INSTALL_LIB_DIR}/Doc CACHE PATH "Installation directory for (some) documentation" FORCE)
set(INSTALL_CMAKE_DIR ${extra_prefix}Bin/${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files" FORCE)

set(PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE EXAMPLES DOC CMAKE)
#foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(INSTALL_${p}_DIR_REL ${${var}})
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)

########################################################################
#                   COMPILER AND SYSTEM DETECTION                      #
########################################################################

unset(SYSTEM_LINUX)
unset(SYSTEM_DARWIN)
unset(SYSTEM_AIX)
unset(SYSTEM_WINDOWS)
unset(SYSTEM_CYGWIN)

message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYSTEM_LINUX)
    set(SYSTEM_LINUX 1 CACHE INTERNAL "Is system Linux?")
endif()
if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_definitions(-DSYSTEM_DARWIN)
    set(SYSTEM_DARWIN 1 CACHE INTERNAL "Is system Mac?")
endif()
if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    add_definitions(-DSYSTEM_AIX)
    set(SYSTEM_AIX 1 CACHE INTERNAL "Is system AIX?")
endif()
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions(-DSYSTEM_WINDOWS)
    set(SYSTEM_WINDOWS 1 CACHE INTERNAL "Is system WINDOWS?")
endif()
if(CYGWIN)
    add_definitions(-DSYSTEM_CYGWIN)
    set(SYSTEM_CYGWIN 1 CACHE INTERNAL "Is cygwin used?")
endif()

# set variables that provide information about the system
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")
# Check Compiler
unset(COMP_CLANG)
unset(COMP_GNU)
unset(COMP_INTEL)
unset(COMP_MSVC)
unset(COMP_UNSPEC)

if (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    add_definitions(-DCLANG)
    set(COMP_CLANG 1 CACHE INTERNAL "Is comppiler Clang?")
elseif (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    add_definitions(-DGNU)
    set(COMP_GNU 1 CACHE INTERNAL "Is comppiler GNU?")
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    add_definitions(-DINTEL)
    set(COMP_INTEL 1 CACHE INTERNAL "Is comppiler Intel?")
elseif (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    add_definitions(-DMSVC)
    set(COMP_MSVC 1 CACHE INTERNAL "Is comppiler MSVC?")
else()
    add_definitions(-DUNSPECCOMP)
    set(COMP_UNSPEC 1 CACHE INTERNAL "Is comppiler something else?")
endif()

unset(COMP_MINGW)
if(MINGW)
    add_definitions(-DMINGW)
    set(COMP_MINGW 1 CACHE INTERNAL "Is comppiler also MINGW? (GNU flag is always true in this case!)")
endif(MINGW)

########################################################################
#                            COMPILER FLAGS                            #
########################################################################

if (CMAKE_CXX_COMPILER_ID MATCHES GNU OR
    CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(warnings "-pedantic") #-Wall -Wextra
elseif (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    #set(warnings "/W4 /WX /EHsc")
	set(warnings /wd4005)
	set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON)
else()
    message(STATUS "Using default flags for compiler")
endif()

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${warnings}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings}")

########################################################################
#             WINDOWS AND SUNDIALS RELATED OPTIONS/SETTINGS            #
########################################################################

set(CMAKE_LEGACY_CYGWIN_WIN32 0)

option(SUNDIALS_USE_STATIC_LIBRARIES "Use static linking for sundials libraries" ON)
if(WIN32)
  if(MINGW)
    option(CMAKE_STATIC_LINK "Use static linking (recommended for Windows)" ON)
  endif(MINGW)
  if(CMAKE_STATIC_LINK AND MINGW)
    set(CMAKE_EXE_LINKER_FLAGS "-static" CACHE STRING "Compiler options for linking executables" FORCE)
  endif(CMAKE_STATIC_LINK AND MINGW)
endif(WIN32)


########################################################################
#               FIND STANDARD LIBRARY HEADERS AND ALIKE                #
########################################################################

include(CheckIncludeFiles)
CHECK_INCLUDE_FILES (unistd.h HAVE_UNISTD_H)
if(HAVE_UNISTD_H)
    message(STATUS "unistd.h header found")
else()
    message(STATUS "unistd.h header not found")
endif(HAVE_UNISTD_H)

CHECK_INCLUDE_FILES (direct.h HAVE_DIRECT_H)
if(HAVE_DIRECT_H)
    message(STATUS "direct.h header found")
else()
    message(STATUS "direct.h header not found")
endif(HAVE_DIRECT_H)

CHECK_INCLUDE_FILES (alloca.h HAVE_ALLOCA_H)

if(HAVE_ALLOCA_H)
    message(STATUS "alloca.h header found; compiling alloca.c is not required")
else()
    message(STATUS "alloca.h header not found; compiling alloca.c is required")
endif(HAVE_ALLOCA_H)

CHECK_INCLUDE_FILES (string.h HAVE_STRING_H)
if(HAVE_STRING_H)
    message(STATUS "string.h header found")
else()
    message(STATUS "string.h header not found")
endif(HAVE_STRING_H)

CHECK_INCLUDE_FILES (getopt.h HAVE_GETOPT_H)
if(HAVE_GETOPT_H)
  message(STATUS "getopt.h header found")
else()
if(NOT CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  message(FATAL_ERROR "The function getopt() as defined in the POSIX standard is not avaiable.")
endif(NOT CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  message(STATUS "getopt.h header not found; using our own implementation for MSVC")
endif(HAVE_GETOPT_H)

CHECK_INCLUDE_FILES (stdlib.h HAVE_STDLIB_H)
if(HAVE_STDLIB_H)
  message(STATUS "stdlib.h header found")
else()
  message(STATUS "stdlib.h header not found; using implementation in the flamemaster package")
endif(HAVE_STDLIB_H)

########################################################################
#           AUTOMATIC INSTALLTION OF THIRD PARTY PROGRAMS              #
########################################################################


# Path to ThirdParty directory
set(INSTALL_THIRDPARTY_DIR ${PROJECT_BINARY_DIR}/../ThirdParty CACHE PATH
  "Installation directory for third party libraries")

IF(COMPILE_SCANMANNEW)
	option(INSTALL_CEREAL "Install cereal library for serialization." ON)
ENDIF(COMPILE_SCANMANNEW)
option(INSTALL_SUNDIALS "Install the sundials package which provides solvers used in FlameMaster. Do not switch this ON if you already installed sundials." OFF)
option(SUNDIALS_NO_IDAS_SEARCH
"Do not search for IDAS (part of sundials).
Currently, we do not use IDAS, but this might cause problems in the future.
Sundials from Zypper (on openSUSE) does not provide IDAS. You can use this option. "
  OFF)

if(WIN32 AND NOT CYGWIN)
  option(INSTALL_WINFLEXBISON "Install winflexbison executables for parsers in ScanMan and FlameMan." ON)
endif(WIN32 AND NOT CYGWIN)

if(WIN32 AND NOT CYGWIN AND INSTALL_WINFLEXBISON)
  #download and unpack winflexbison at configure time
  set(flex_bison_path "${INSTALL_THIRDPARTY_DIR}/winflexbison/")
  set(myflex "win_flex.exe")
  set(mybison "win_bison.exe")
  if(NOT EXISTS "${flex_bison_path}${myflex}" OR NOT EXISTS "${flex_bison_path}${mybison}")
    configure_file("${PROJECT_SOURCE_DIR}/../cmake/thirdParty/WinFlexBison.cmake.in"
                 "${INSTALL_THIRDPARTY_DIR}/winflexbison/CMakeLists.txt")
	execute_process(COMMAND ${CMAKE_COMMAND} .
		WORKING_DIRECTORY ${flex_bison_path} )
  else()
	message(STATUS "WinFlexBison found")
  endif()
  find_file(FLEX_EXECUTABLE "${myflex}" "${flex_bison_path}")
  find_file(BISON_EXECUTABLE "${mybison}" "${flex_bison_path}")
endif(WIN32 AND NOT CYGWIN AND INSTALL_WINFLEXBISON)

if(INSTALL_CEREAL)
  #download and unpack cereal at configure time
  configure_file("${PROJECT_SOURCE_DIR}/../cmake/thirdParty/Cereal.cmake.in"
                 "${INSTALL_THIRDPARTY_DIR}/cereal/CMakeLists.txt")
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY ${INSTALL_THIRDPARTY_DIR}/cereal )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}
    WORKING_DIRECTORY ${INSTALL_THIRDPARTY_DIR}/cereal )
endif(INSTALL_CEREAL)

if(INSTALL_SUNDIALS)
  #download and unpack cereal at sundials time
  configure_file("${PROJECT_SOURCE_DIR}/../cmake/thirdParty/Sundials.cmake.in"
                 "${INSTALL_THIRDPARTY_DIR}/sundials/CMakeLists.txt")
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY ${INSTALL_THIRDPARTY_DIR}/sundials )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}
    WORKING_DIRECTORY ${INSTALL_THIRDPARTY_DIR}/sundials )
    message(STATUS "FlameMaster installed sundials for you in: ${INSTALL_THIRDPARTY_DIR}/sundials\n
        This installation is automatically used if you:\n
        1) do not set the PATH_SUNDIALS_ROOT environment variable to a valid sundials installation path
        2) do not set the PATH_SUNDIALS_ROOT CMakeCache variable to a valid sundials installation path
        3) did not find a different sundials installation in a previous configuration step (\"cmake ../Repository ...\")
        4) set the INSTALL_SUNDIALS option to on (e.g. \"cmake ../Repository -DINSTALL_SUNDIALS=ON\")
        ")
    set(SUNDIALS_USER_PATHS ${SUNDIALS_USER_PATHS}
        ${INSTALL_THIRDPARTY_DIR}/sundials
        )
else()
  message(STATUS "Sundials will not be installed by FlameMaster")
  PRINT_SUNDIALS_HINT_INSTALL()
endif(INSTALL_SUNDIALS)

########################################################################
#                   FILES THAT REQUIRE CONFIGURATION                   #
########################################################################

file(READ "${CMAKE_SOURCE_DIR}/config/DescriptionTemplate1.in" DescriptionTemplate1 )
file(READ "${CMAKE_SOURCE_DIR}/config/DescriptionTemplate2.in" DescriptionTemplate2 )

set(ScanMan_TEST_MECH_INPUT_PATH_CH4 ${CMAKE_SOURCE_DIR}/examples/ScanMan/CH4/CH4.72.mech)
set(ScanMan_TEST_MECH_INPUT_PATH ${CMAKE_BINARY_DIR}/src/ScanManNew/sandiego20141004_mechCK.txt)
set(ScanMan_TEST_THERM_INPUT_PATH ${CMAKE_BINARY_DIR}/src/ScanManNew/sandiego20120907_therm.txt)
set(ScanMan_TEST_TRANS_INPUT_PATH ${CMAKE_BINARY_DIR}/src/ScanManNew/sandiego20120907_trans.txt)

# configure a header file to pass CMake settings to the source code
configure_file (
  "${CMAKE_SOURCE_DIR}/config/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h"
  )

set(BASH_SOURCE_FM Source.bash)
set(ZSH_SOURCE_FM Source.zsh)
set(CSH_SOURCE_FM Source.csh)
set(BAT_SOURCE_FM Source.BAT)
set(SourceBashTemplate "${CMAKE_SOURCE_DIR}/config/SourceTemplate.bash.in") # input path
set(SourceZshTemplate "${CMAKE_SOURCE_DIR}/config/SourceTemplate.zsh.in")
set(SourceCshTemplate "${CMAKE_SOURCE_DIR}/config/SourceTemplate.csh.in")
set(SourceBATTemplate "${CMAKE_SOURCE_DIR}/config/SourceTemplate.BAT.in")
set(SourceBashBuild "${CMAKE_CURRENT_BINARY_DIR}/install_config/${PROJECT_NAME}/${BASH_SOURCE_FM}")
set(SourceZshBuild "${CMAKE_CURRENT_BINARY_DIR}/install_config/${PROJECT_NAME}/${ZSH_SOURCE_FM}")
set(SourceCshBuild "${CMAKE_CURRENT_BINARY_DIR}/install_config/${PROJECT_NAME}/${CSH_SOURCE_FM}")
set(SourceBATBuild "${CMAKE_CURRENT_BINARY_DIR}/install_config/${PROJECT_NAME}/${BAT_SOURCE_FM}")

# This has to be appropriate for the value of ${INSTALL_BIN_DIR_REL}  which
# is location of the 'source' scripts that can be used to set up the
# environment. If you use Bin/bin/ REL_SOURCE should be ../..
# TO-DO: automate this!
set(REL_SOURCE "../..")

configure_file(
  "${CMAKE_SOURCE_DIR}/config/configure_source_scripts.cmake.in"
  "${PROJECT_BINARY_DIR}/configure_source_scripts.cmake" @ONLY)

# configure and install a script that can be used to source/run the installation
# if(WIN32)
#	message(STATUS "We should provide a batch script to run FlameMaster on Windows")
#else()
add_custom_command(OUTPUT ${SourceBashBuild} ${SourceZshBuild} ${SourceCshBuild} ${SourceBATBuild}
  COMMAND ${CMAKE_COMMAND} -P
    "${PROJECT_BINARY_DIR}/configure_source_scripts.cmake"
  DEPENDS
    "${PROJECT_BINARY_DIR}/configure_source_scripts.cmake"
    ${SourceBashTemplate}
    ${SourceZshTemplate}
    ${SourceCshTemplate}
    ${SourceBATTemplate}
  COMMENT "Configuring scripts that can be used to conveniently set up the environment for the FlameMaster package"
  VERBATIM)

add_custom_target(configure_source_sripts ALL
  DEPENDS ${SourceBashBuild} ${SourceZshBuild} ${SourceCshBuild} ${SourceBATBuild})

# If you change the installation destination of the source scripts
# make sure that you adjust REL_SOURCE appropriately
install(FILES
	${SourceBashBuild} ${SourceZshBuild} ${SourceCshBuild} ${SourceBATBuild}
  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
	DESTINATION ${INSTALL_BIN_DIR_REL} COMPONENT "SourceScripts")
#endif()

########################################################################
#                             INSTALLATION                             #
########################################################################

# add the binary tree to the search path for include files so that we will find Config.h
include_directories("${PROJECT_BINARY_DIR}")

# RPATH settings# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/Bin/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/Bin/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/Bin/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

install(DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ DESTINATION ${INSTALL_EXAMPLES_DIR_REL} USE_SOURCE_PERMISSIONS COMPONENT "FlameMaster_Examples")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/doc/ DESTINATION ${INSTALL_DOC_DIR_REL} COMPONENT "FlameMaster_Documentation")

########################################################################
#                             UNIT TESTING                             #
########################################################################

option(TESTS "Install google-test framework for unit tests and enable units tests." OFF)

if(TESTS)
  #download and unpack gtest at configure time
  cmake_minimum_required(VERSION 3.0.0)
  cmake_policy(SET CMP0057 NEW)
  # Prevent GoogleTest from overriding our compiler/linker options
  # when building with Visual Studio
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  configure_file("${PROJECT_SOURCE_DIR}/../cmake/thirdParty/GTest.cmake.in"
                 "${INSTALL_THIRDPARTY_DIR}/googletest/src/CMakeLists.txt")
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY ${INSTALL_THIRDPARTY_DIR}/googletest/src )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}
    WORKING_DIRECTORY ${INSTALL_THIRDPARTY_DIR}/googletest/src )

endif(TESTS)

# Section for Googletest
IF(TESTS)
  set(TEST_UTILITY_DIR ${CMAKE_BINARY_DIR}/test_utility CACHE INTERNAL "Path for test utilities shared by FlameMan and ScanMan" FORCE)
  set(PREFILE_01 "TEST_01.pre" CACHE INTERNAL "Data contained in the prefile is used to test pure species properties \
    such as thermodynamic and transport properties" FORCE)

  if (TESTS OR EXISTS ${INSTALL_THIRDPARTY_DIR}/googletest/)
     # find and include GoogleTest to run tests
     set(CMAKE_PREFIX_PATH ${INSTALL_THIRDPARTY_DIR}/googletest/)
  endif(TESTS OR EXISTS ${INSTALL_THIRDPARTY_DIR}/googletest/)

  find_package(GTest REQUIRED)

  enable_testing() # statement has to be before 'add_subdirectory' to define tests in each sublevel CMakeLists

ENDIF(TESTS)

# NOTE: The order matters! The most independent ones go first.
add_subdirectory (${PROJECT_SOURCE_DIR}/libraries)
add_subdirectory (${PROJECT_SOURCE_DIR}/FlameMan)
add_subdirectory (${PROJECT_SOURCE_DIR}/ScanMan)
IF(COMPILE_SCANMANNEW)
	add_subdirectory (${PROJECT_SOURCE_DIR}/ScanManNew)
ENDIF(COMPILE_SCANMANNEW)
add_subdirectory (${PROJECT_SOURCE_DIR}/Tools)

########################################################################
#                         CPACK IMPLEMENTATION                         #
########################################################################

include(InstallRequiredSystemLibraries)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "FlameMaster is a C++ Program Package for 0D Combustion and 1D Laminar Flame Calculations")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${FlameMaster_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${FlameMaster_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${FlameMaster_VERSION_PATCH}")
set(CPACK_PACKAGE_CONTACT "Raymond Langer")


SET(CPACK_PACKAGE_INSTALL_DIRECTORY "FlameMaster")
set(CPACK_PACKAGE_VENDOR "Institute for Combustion Technology")
IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backslashes.
  #SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/icons\\\\FlameMasterL.ico")
  #SET(CPACK_NSIS_INSTALLED_ICON_NAME "${CMAKE_SOURCE_DIR}/icons\\\\FlameMasterL.ico")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "Bin\\\\bin\\\\FlameMan.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "FlameMaster")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.itv.rwth-aachen.de")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.itv.rwth-aachen.de\\\\downloads\\\\flamemaster")
  SET(CPACK_NSIS_CONTACT "r.langer@itv.rwth-aachen.de")
  SET(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
  SET(CPACK_NSIS_MODIFY_PATH ON)
  #SET(CPACK_NSIS_EXECUTABLES_DIRECTORY "bin")
  SET(CPACK_NSIS_INSTALL_DIRECTORY)
ELSE(WIN32 AND NOT UNIX)
	SET(CPACK_BUNDLE_NAME "FlameMaster")
	#SET(CPACK_BUNDLE_PLIST "icons/Info.plist")
  #SET(CPACK_BUNDLE_ICON "${CMAKE_SOURCE_DIR}/icons/FlameMasterL.ico")
  SET(CPACK_SET_DESTDIR "FlameMaster")
  SET(CPACK_STRIP_FILES "Bin/bin/FlameMan")
  SET(CPACK_SOURCE_STRIP_FILES "")
ENDIF(WIN32 AND NOT UNIX)
SET(CPACK_PACKAGE_EXECUTABLES "FlameMan" "ScanMan" "CreateBinFile" "Converter")
include(CPack)
