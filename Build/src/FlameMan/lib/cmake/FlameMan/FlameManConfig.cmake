# - Config file for the FlameMan package
# It defines the following variables:

#  FlameMan_BIN_DIR                      - location of FlameMan
#  FlameMan_VERSION                      - version of the binary
#  FlameMan_ConfigPackageLocation        - location of the configuration files (= FlameMan_DIR)
# Imported Targets:
#  FlameMan_BINS                         - all executables
#  FlameMan::FlameMan                    - executable target


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was FlameManConfig.cmake.in                            ########

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
set(cfg_hint "The problem is most likely that the generation of FlameManConfig.cmake in FlameMan/CMakeLists.txt is not up-to-date.")
set(FlameMan_VERSION "4.0.0")
message("-- Found FlameMan, Version ${FlameMan_VERSION}")
if(0)
	set_and_check(FlameMan_BIN_DIR "")
	set_and_check(FlameMan_ConfigPackageLocation "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/lib/cmake/FlameMan")
	if (NOT PACKAGE_PREFIX_DIR STREQUAL "")
		message("-- Using 'FlameMan' from '${PACKAGE_PREFIX_DIR}'")
	endif()
else()
	set(FlameMan_BIN_DIR "")
	# check when find_package is used whether directory exists
	if(NOT EXISTS "${FlameMan_BIN_DIR}")
		message(FATAL_ERROR "FlameMan_BIN_DIR not found: ${FlameMan_BIN_DIR}"
							" ${cfg_hint}")
	endif()
	set_and_check(FlameMan_ConfigPackageLocation "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/FlameMan/lib/cmake/FlameMan")
	message("-- Using 'FlameMan' from a build tree")
	if (NOT FlameMan_BIN_DIR STREQUAL "")
		message("\t binary dir.:        '${FlameMan_BIN_DIR}'")
	endif()
	if (NOT FlameMasterLibs_ConfigPackageLocation STREQUAL "")
		message("\t target import:      '${FlameMan_ConfigPackageLocation}/FlameManTargets.cmake'")
	endif()
	if (NOT FlameMan_ConfigPackageLocation STREQUAL "")
		message("\t package location:   '${FlameMan_ConfigPackageLocation}'")
	endif()
endif()

include(${FlameMan_ConfigPackageLocation}/FlameManTargets.cmake)
# check whether configuration file is still up to date
set(FlameMan_BINS FlameMan::FlameMan )
foreach(t ${FlameMan_BINS})
	if(NOT TARGET ${t})
		message(FATAL_ERROR "Expected target does not exist: '${t}'\n"
			" Note that targets are case sensitive."
			" ${cfg_hint}")
	endif()
endforeach()
