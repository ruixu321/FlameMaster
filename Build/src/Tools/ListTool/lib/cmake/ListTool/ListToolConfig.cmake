# - Config file for the ListTool package
# It defines the following variables:

#  ListTool_INCLUDE_DIRS               - include directories for Tools public headers
#  ListTool_LIBS_DIRS                  - location of the ListTool libraries (files if static libraries are used)
#  ListTool_BIN_DIR                    - location of the ListTool (files if static libraries are used)
#  ListTool_VERSION                    - version of the library package
#  ListTool_ConfigPackageLocation      - location of the ListToolConfig.cmake and ListToolConfigVersion.cmake files

# Imported Targets:
#  ListTool_LIBRARIES           - (all) libraries one might want to link against (consider specific 
#                                 libraries below)
#  ListTool_BINARY              - ListTool binary for reformating FlameMaster output


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was ListToolConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../../../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

find_package(FlameMasterLibs REQUIRED HINTS ${PACKAGE_PREFIX_DIR}) # try FlameMasterLibs from the same installation

set(cfg_hint "The problem is most likely that the generation of ListToolConfig.cmake in Tools/ListTool/CMakeLists.txt is not up-to-date.")
set(ListTool_VERSION "0.1.0")
message("-- Found ListTool, Version ${ListTool_VERSION}")
if(0)
	set_and_check(ListTool_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/Bin/include")
	set_and_check(ListTool_LIBS_DIRS "${PACKAGE_PREFIX_DIR}/Bin/lib")
	set_and_check(ListTool_BIN_DIR "")
	set_and_check(ListTool_ConfigPackageLocation "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/lib/cmake/ListTool")
	if (NOT PACKAGE_PREFIX_DIR STREQUAL "")
		message("-- Using 'ListTool' from '${PACKAGE_PREFIX_DIR}'")
	endif()
else()
	# check include directory (should already exist, since build config will not be installed)
	FOREACH(directory ${ListTool_src_header_direcotries})
		if(NOT EXISTS "${directory}")
	    	message(FATAL_ERROR "Directory from '${ListTool_src_header_direcotries}' not found: ${direcotry}")
		endif()
	ENDFOREACH()
	set(ListTool_INCLUDE_DIRS "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/Tools/ListTool/FLReader;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/Tools/ListTool/List;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/Tools/ListTool/ListTool;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/Tools/ListTool/WSS")
	# check again when find_package is used
	if(ListTool_INCLUDE_DIRS STREQUAL "")
		message(FATAL_ERROR "ListTool_INCLUDE_DIRS should contain source header directories, but it is empty."
							" ${cfg_hint}")
	endif()
	FOREACH(directory ${ListTool_INCLUDE_DIRS})
		if(NOT EXISTS "${directory}")
	    	message(FATAL_ERROR "Directory from 'ListTool_INCLUDE_DIRS' not found: ${direcotry}"
					"It seems like this build directory cannot find the source directory of the ListTool anymore!"
					" Did you change the location of the build or the source tree after you executed cmake?"
					" ${cfg_hint}")
		endif()
	ENDFOREACH()
	# library directory (cannot be checked at configuration time)
	set(ListTool_LIBS_DIRS "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/List;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/FLReader;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/WSS")
	# check when find_package is used
	if(ListTool_LIBS_DIRS STREQUAL "")
		message(FATAL_ERROR "ListTool_LIBS_DIRS should contain library directories, but it is empty."
							" ${cfg_hint}")
	endif()
	FOREACH(file ${ListTool_LIBS_DIRS})
		if(NOT EXISTS "${file}")
	    	message(FATAL_ERROR "Directory from 'ListTool_LIBS_DIRS' not found: ${file}"
					"It seems like this build directory does not contain expected directories or files!"
					" ${cfg_hint}")
		endif()
	ENDFOREACH()
	# set binary/binaries
	set(ListTool_BIN_DIR "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/ListTool")
	FOREACH(file ${ListTool_BIN_DIR})
		if(NOT EXISTS "${file}")
	    	message(FATAL_ERROR "Directory from 'ListTool_BIN_DIR' not found: ${file}"
					"It seems like this build directory does not contain the expected executable!"
					" ${cfg_hint}")
		endif()
	ENDFOREACH()
	set_and_check(ListTool_ConfigPackageLocation "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/ListTool/lib/cmake/ListTool")
	message("-- Using 'ListTool' from a build tree")
	if (NOT ListTool_LIBS_DIRS STREQUAL "")
		message("\t libraries:      '${ListTool_LIBS_DIRS}'")
	endif()
	if (NOT ListTool_INCLUDE_DIRS STREQUAL "")
		message("\t public headers: '${ListTool_INCLUDE_DIRS}'")
	endif()
	if (NOT ListTool_BIN_DIR STREQUAL "")
		message("\t binary dir.:    '${ListTool_BIN_DIR}'")
	endif()
	if (NOT ListTool_ConfigPackageLocation STREQUAL "")
		message("\t target import:  '${ListTool_ConfigPackageLocation}/ListToolTargets.cmake'")
	endif()
endif()

include(${ListTool_ConfigPackageLocation}/ListToolTargets.cmake)
set(ListTool_LIBRARIES ListTool::FLReader ListTool::WSS ListTool::List)
set(ListTool_BINARY ListTool::ListTool)
foreach(t ${ListTool_LIBRARIES})
	if(NOT TARGET ${t})
		message(FATAL_ERROR "Expected target does not exist: '${t}'\n"
			" Note that targets are case sensitive."
			" ${cfg_hint}")
	endif()
endforeach()
foreach(t ${ListTool_BINARY})
	if(NOT TARGET ${t})
		message(FATAL_ERROR "Expected target does not exist: '${t}'\n"
			" Note that targets are case sensitive."
			" ${cfg_hint}")
	endif()
endforeach()
