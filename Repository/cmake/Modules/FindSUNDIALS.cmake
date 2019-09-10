# - Try to find SUNDIALS
#

# Find SUNDIALS, the SUite of Nonlinear and DIfferential/ALgebraic equation Solvers.
#
# The module will optionally accept the COMPONENTS argument.  If no COMPONENTS
# are specified, then the find module will default to find all the SUNDIALS
# libraries.  If one or more COMPONENTS are specified, the module will attempt to
# find the specified components.
#
# On UNIX systems, this module will read the variable SUNDIALS_USE_STATIC_LIBRARIES
# to determine whether or not to prefer a static link to a dynamic link for SUNDIALS
# and all of it's dependencies.  To use this feature, make sure that the
# SUNDIALS_USE_STATIC_LIBRARIES variable is set before the call to find_package.
#
# To provide the module with a hint about where to find your SUNDIALS installation,
# you can set the environment variable PATH_SUNDIALS_ROOT. The FindSUNDIALS module will
# then look in this path when searching for SUNDIALS paths and libraries.
#
# This module will define the following variables:
#  SUNDIALS_INCLUDE_DIRS    - Location of the SUNDIALS includes
#  SUNDIALS_FOUND           - true if SUNDIALS was found on the system
#  SUNDIALS_LIBRARIES       - Required libraries for all requested components
#  SUNDIALS_VERSION         - Contains the sundials version string (e.g. "2.7.0" without ")
#  SUNDIALS_VERSION_MAJOR   - Contains the sundials version string (e.g. "2" without ")
#  SUNDIALS_VERSION_MINOR   - Contains the sundials version string (e.g. "7" without ")
#  SUNDIALS_VERSION_PATCH   - Contains the sundials version string (e.g. "0" without ")

option(SUNDI_DIAGNOS_SEARCH_LIB_FILES "Search SUNDIALS libraries in user specified directories. Only for troubleshooting!" OFF)
option(SUNDI_DIAGNOS_SEARCH_HEADER_FILES "Search SUNDIALS headers in user specified directories. Only for troubleshooting!" OFF)

include(${CMAKE_SOURCE_DIR}/cmake/Macros.cmake)
include(FindPackageHandleStandardArgs)

# Option that allows users to specify custom SUNDIALS path
if (NOT "$ENV{PATH_SUNDIALS_ROOT}" STREQUAL "")
    list (APPEND SUNDIALS_USER_PATHS "$ENV{PATH_SUNDIALS_ROOT}")
    list (APPEND SUNDIALS_USER_PATHS "$ENV{PATH_SUNDIALS_ROOT}")
    message(STATUS "environment variable PATH_SUNDIALS_ROOT set to: $ENV{PATH_SUNDIALS_ROOT}")
else ()
    message(STATUS "Environment variable PATH_SUNDIALS_ROOT is not used")
    PRINT_SUNDIALS_HINT_ENV()
endif ()

if (NOT PATH_SUNDIALS_ROOT)
   set (PATH_SUNDIALS_ROOT "" CACHE STRING "Optional path to the SUNDIALS lib and include direrctory" FORCE)
else ()
   list (APPEND SUNDIALS_USER_PATHS "${PATH_SUNDIALS_ROOT}")
   list (APPEND SUNDIALS_USER_PATHS "${PATH_SUNDIALS_ROOT}")
endif ()

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  set(SUNDIALS_USER_PATHS ${SUNDIALS_USER_PATHS} /home/itv/lib/sundials/2.4.0_intel12)
elseif(CMAKE_CXX_COMPILER_ID MATCHES GNU OR CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(SUNDIALS_USER_PATHS ${SUNDIALS_USER_PATHS} /home/itv/lib/sundials/2.4.0_gcc44/)
endif()

set(SUNDIALS_USER_PATHS ${SUNDIALS_USER_PATHS}
#    "C:/Users/itv/Downloads/sundials-2.6.2/sundials-2.6.2-install"
#    $ENV{HOME}/.local
#    $ENV{HOME}
#    $ENV{HOME}/sundials
#    $ENV{HOME}/sundials/installation
#    $ENV{HOME}/sundials/install
#    /usr
)

if(NOT DEFINED SUNDIALS_INCLUDE_DIR 
    OR SUNDIALS_INCLUDE_DIR STREQUAL "SUNDIALS_INCLUDE_DIR-NOTFOUND")
    message(STATUS "Searching sundials in: ${SUNDIALS_USER_PATHS}")
else()
    message(STATUS "Sundials location already found and stored in the cache during a previous configuration step (Remove files in the Build directory to search again or modify the cache yourself using the ccmake gui)")
endif()

set(CMAKE_PREFIX_PATH ${SUNDIALS_USER_PATHS} ${CMAKE_PREFIX_PATH})

# find package manager version of sundials
if (APPLE)
	list(APPEND CMAKE_PREFIX_PATH "/usr/local")
	list(APPEND CMAKE_PREFIX_PATH "/opt")
	list(APPEND CMAKE_PREFIX_PATH "/opt/local")
    #list (APPEND SUNDIALS_USER_PATHS "/opt/local/sundials")
    #list (APPEND SUNDIALS_USER_PATHS "/opt/sundials")
endif()
if (WIN32)
	# Windows...
endif()
if (UNIX)
	list(APPEND CMAKE_PREFIX_PATH "/usr/local")
	list(APPEND CMAKE_PREFIX_PATH "/opt")
	list(APPEND CAMKE_PREFIX_PATH "/usr")
endif()

# List of the valid SUNDIALS components
set( SUNDIALS_VALID_COMPONENTS
    sundials_cvode
    sundials_cvodes
    sundials_ida
    sundials_kinsol
    sundials_nvecserial
)

if(NOT SUNDIALS_NO_IDAS_SEARCH)
    list(APPEND 
    SUNDIALS_VALID_COMPONENTS
    sundials_idas)
endif(NOT SUNDIALS_NO_IDAS_SEARCH)

if( NOT SUNDIALS_FIND_COMPONENTS )
    set( SUNDIALS_WANT_COMPONENTS ${SUNDIALS_VALID_COMPONENTS} )
else()
    # add the extra specified components, ensuring that they are valid.
    foreach( _COMPONENT ${SUNDIALS_FIND_COMPONENTS} )
        string (TOLOWER ${_COMPONENT} _COMPONENT_LOWER)
        list( FIND SUNDIALS_VALID_COMPONENTS ${_COMPONENT_LOWER} COMPONENT_LOCATION )
        if( ${COMPONENT_LOCATION} EQUAL -1 )
            message( FATAL_ERROR
                "\"${_COMPONENT_LOWER}\" is not a valid SUNDIALS component." )
        else()
            list( APPEND SUNDIALS_WANT_COMPONENTS ${_COMPONENT_LOWER} )
        endif()
    endforeach()
endif()

set(header_search sundials_types.h)
# find the SUNDIALS include directories
find_path( SUNDIALS_INCLUDE_DIR ${header_search}
    ENV
        PATH_SUNDIALS_ROOT
    PATHS
        ${SUNDIALS_USER_PATHS}
    PATH_SUFFIXES
        include
        include/sundials
)

set( SUNDIALS_INCLUDE_DIRS
     "${SUNDIALS_INCLUDE_DIR}/.."
     "${SUNDIALS_INCLUDE_DIR}/../cvode"
     "${SUNDIALS_INCLUDE_DIR}/../cvodes"
     "${SUNDIALS_INCLUDE_DIR}/../ida"
     "${SUNDIALS_INCLUDE_DIR}/../idas"
     "${SUNDIALS_INCLUDE_DIR}/../kinsol"
     "${SUNDIALS_INCLUDE_DIR}/../nvector"
    "${SUNDIALS_INCLUDE_DIR}"
)

if(SUNDI_DIAGNOS_SEARCH_LIB_FILES OR SUNDI_DIAGNOS_SEARCH_HEADER_FILES)
    message(STATUS "Run diagnostic search... this my take some minutes... Later USE UP AND DOWN ARROWS TO to see all hints!")
endif()

# diagnostics for header files search
if(SUNDIALS_INCLUDE_DIR STREQUAL "SUNDIALS_INCLUDE_DIR-NOTFOUND" AND SUNDI_DIAGNOS_SEARCH_HEADER_FILES)
	set(found_files "")
	message("Header not found, searched for: ${header_search}")
	message("Note: The following search is much more tolerant than the search for the right header file search...")
	message("Trying to find promising files in (user specified) search paths. This might give you a hint what's wrong...")
	foreach(DIR ${SUNDIALS_USER_PATHS})
        message("Checking directory (recursively): ${DIR}")
		SEARCH_FILES(files "${DIR}/sundials*.h")
        list(APPEND found_files ${files})
    endforeach()
    message("Found the following files:") 
	foreach(result_path ${found_files})  
		message("\t${result_path}")
	endforeach()
    if(NOT found_files STREQUAL "")
    	message(FATAL_ERROR "Maybe the header_search ('${header_search}') is not up-to-date!")
	else()
		message(FATAL_ERROR "Did not find any promising files. Did you miss the right directory?")
    endif()
endif()

# find the SUNDIALS libraries
foreach( LIB ${SUNDIALS_WANT_COMPONENTS} )
    if( UNIX AND SUNDIALS_USE_STATIC_LIBRARIES )
        # According to bug 1643 on the CMake bug tracker, this is the
        # preferred method for searching for a static library.
        # See http://www.cmake.org/Bug/view.php?id=1643.  We search
        # first for the full static library name, but fall back to a
        # generic search on the name if the static search fails.
        set( THIS_LIBRARY_SEARCH lib${LIB}.a ${LIB} )
    elseif(SUNDIALS_USE_STATIC_LIBRARIES AND MINGW)
        set( THIS_LIBRARY_SEARCH lib${LIB}.a )
    else()
        set( THIS_LIBRARY_SEARCH ${LIB} )
    endif()
    find_library( SUNDIALS_${LIB}_LIBRARY
        NAMES ${THIS_LIBRARY_SEARCH} ${THIS_LIBRARY_SEARCH}.so
        ENV
            PATH_SUNDIALS_ROOT
        PATHS
            ${SUNDIALS_USER_PATHS}
        PATH_SUFFIXES
            lib
            Lib
    )
    
    if(SUNDIALS_${LIB}_LIBRARY STREQUAL "SUNDIALS_${LIB}_LIBRARY-NOTFOUND" AND SUNDI_DIAGNOS_SEARCH_LIB_FILES)
    	set(found_files "")
    	message("Library not found: ${LIB}")
    	message("Note: The following search is much more tolerant than the search for libraries...")
		message("Trying to find promising files in (user specified) search paths. This might give you a hint what's wrong...")
		foreach(DIR ${SUNDIALS_USER_PATHS})
            message("Checking directory (recursively): ${DIR}")
			SEARCH_FILES(files "${DIR}/lib${THIS_LIBRARY_SEARCH}*")
            list(APPEND found_files ${files})
        endforeach()
	    message("Found the following files:") 
		foreach(result_path ${found_files})  
			message("\t${result_path}")
		endforeach()
        if(NOT found_files STREQUAL "")
        	message(FATAL_ERROR "Maybe there are no symbolic links. If you used a package manager to install sundials"
                    " you might try to install dev packages as well!")
		else()
			message(FATAL_ERROR "Did not find any promising files. Did you miss the right directory?")
        endif()
    endif()
    list( APPEND SUNDIALS_LIBRARIES ${SUNDIALS_${LIB}_LIBRARY} )
    mark_as_advanced(SUNDIALS_${LIB}_LIBRARY)
endforeach()

file(STRINGS "${SUNDIALS_INCLUDE_DIR}/sundials_config.h" _sundials_VERSION_HPP_CONTENTS REGEX "#define SUNDIALS_(PACKAGE_)?VERSION ")
set(_VERSION_NUMBER "([0-9]+)")
string(REGEX MATCH "${_VERSION_NUMBER}\\.${_VERSION_NUMBER}\\.${_VERSION_NUMBER}" SUNDIALS_VERSION "${_sundials_VERSION_HPP_CONTENTS}")
string(REPLACE "." ";" VERSION_LIST ${SUNDIALS_VERSION})
list(GET VERSION_LIST 0 SUNDIALS_VERSION_MAJOR)
list(GET VERSION_LIST 1 SUNDIALS_VERSION_MINOR)
list(GET VERSION_LIST 2 SUNDIALS_VERSION_PATCH)

unset(_sundials_VERSION_HPP_CONTENTS)

find_package_handle_standard_args( SUNDIALS REQUIRED_VARS SUNDIALS_LIBRARIES SUNDIALS_INCLUDE_DIRS VERSION_VAR SUNDIALS_VERSION)

if(NOT SUNDIALS_FOUND)
    PRINT_SUNDIALS_HINT_ENV()
    PRINT_SUNDIALS_HINT_INSTALL()
endif()

mark_as_advanced(
    SUNDIALS_LIBRARIES
    SUNDIALS_INCLUDE_DIR
    SUNDIALS_INCLUDE_DIRS
    SUNDIALS_VERSION
    SUNDIALS_VERSION_MAJOR
    SUNDIALS_VERSION_MINOR
    SUNDIALS_VERSION_PATCH
)
