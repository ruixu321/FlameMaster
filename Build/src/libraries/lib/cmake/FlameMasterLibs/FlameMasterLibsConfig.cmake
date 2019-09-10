# - Config file for the FlameMasterLibs package
# It defines the following variables:

#  FlameMasterLibs_INCLUDE_DIRS - include directories for FlameMasterLibs public headers
#  FlameMasterLibs_LIBS_DIRS    - location of the FlameMasterLibs (files if static libraries are used)
#  FlameMasterLibs_VERSION      - version of the library package

# Imported Targets:
#  FlameMasterLibs_LIBRARIES    - (all) libraries one might want to link against (consider specific libraries below)
#  FMLibs::Alligator            - library for memory allocation
#  FMLibs::ArrayMan             - library for (basic) linear algebra, vectors, matrices, and linear solver
#  FMLibs::FM                   - FlameMaster libraries (some mathematical operations, output postprocessing...)
#  FMLibs::Newton               - library for a newton solver
#  FMLibs::dassl                - library for (deprecated) Differential- Algebraic System Solver (DASSL)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was FlameMasterLibsConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../../../../" ABSOLUTE)

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
set(cfg_hint "The problem is most likely that the generation of FlameMasterLibsLibsConfig.cmake in lib/CMakeLists.txt is not up-to-date.")
set(FlameMasterLibs_VERSION "0.1.0")
message("-- Found FlameMaster Libraries, Version ${FlameMasterLibs_VERSION}")
if(0)
	set_and_check(FlameMasterLibs_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/Bin/include")
	set_and_check(FlameMasterLibs_LIBS_DIRS "${PACKAGE_PREFIX_DIR}/Bin/lib")
	set_and_check(FlameMasterLibs_ConfigPackageLocation "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/lib/cmake/FlameMasterLibs")
	if (NOT PACKAGE_PREFIX_DIR STREQUAL "")
		message("-- Using 'FlameManLibs' from '${PACKAGE_PREFIX_DIR}'")
	endif()
else()
	# check include directory (should already exist, since build config will not be installed)
	FOREACH(directory ${FMLibs_src_header_direcotries})
		if(NOT EXISTS "${directory}")
	    	message(FATAL_ERROR "Directory from 'FMLibs_src_header_direcotries' not found: ${direcotry}")
		endif()
	ENDFOREACH()
	set(FlameMasterLibs_INCLUDE_DIRS "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Alligator;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/ArrayMan;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Dassl;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/FM;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton")
	# check again when find_package is used
	if(FlameMasterLibs_INCLUDE_DIRS STREQUAL "")
		message(FATAL_ERROR "FlameMasterLibs_INCLUDE_DIRS should contain source header directories, but it is empty."
							" ${cfg_hint}")
	endif()
	FOREACH(directory ${FlameMasterLibs_INCLUDE_DIRS})
		if(NOT EXISTS "${directory}")
	    	message(FATAL_ERROR "Directory from 'FlameMasterLibs_INCLUDE_DIRS' not found: ${direcotry}"
					"It seems like this build directory cannot find the source directory of the FlameMasterLibs anymore!"
					" Did you change the location of the build or the source tree after you executed cmake?"
					" ${cfg_hint}")
		endif()
	ENDFOREACH()
	# library directory (cannot be checked at configuration time)
	set(FlameMasterLibs_LIBS_DIRS "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Alligator;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/ArrayMan;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/FM;/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton")
	# check when find_package is used
	if(FlameMasterLibs_LIBS_DIRS STREQUAL "")
		message(FATAL_ERROR "FlameMasterLibs_LIBS_DIRS should contain library directories, but it is empty."
							" ${cfg_hint}")
	endif()
	FOREACH(file ${FlameMasterLibs_LIBS_DIRS})
		if(NOT EXISTS "${file}")
	    	message(FATAL_ERROR "Directory from 'FlameMasterLibs_LIBS_DIRS' not found: ${file}"
					"It seems like this build directory does not contain expected directories or files!"
					" ${cfg_hint}")
		endif()
	ENDFOREACH()
	set_and_check(FlameMasterLibs_ConfigPackageLocation "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/lib/cmake/FlameMasterLibs")
	message("-- Using 'FlameManLibs' from a build tree")
	if (NOT FlameMasterLibs_LIBS_DIRS STREQUAL "")
		message("\t libraries:      '${FlameMasterLibs_LIBS_DIRS}'")
	endif()
	if (NOT FlameMasterLibs_INCLUDE_DIRS STREQUAL "")
		message("\t public headers: '${FlameMasterLibs_INCLUDE_DIRS}'")
	endif()
	if (NOT FlameMasterLibs_ConfigPackageLocation STREQUAL "")
		message("\t target import:  '${FlameMasterLibs_ConfigPackageLocation}/FlameMasterLibsTargets.cmake'")
	endif()
endif()

include(${FlameMasterLibs_ConfigPackageLocation}/FlameMasterLibsTargets.cmake)
set(FlameMasterLibs_LIBRARIES FMLibs::Alligator FMLibs::ArrayMan FMLibs::FM FMLibs::Newton)
if(OFF)
	set(FlameMasterLibs_LIBRARIES ${FlameMasterLibs_LIBRARIES} FMLibs::Dassl)
endif(OFF)

foreach(t ${FlameMasterLibs_LIBRARIES})
	if(NOT TARGET ${t})
		message(FATAL_ERROR "Expected target does not exist: '${t}'\n"
			" Note that targets are case sensitive."
			" ${cfg_hint}")
	endif()
endforeach()
