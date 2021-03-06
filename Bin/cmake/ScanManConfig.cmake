# - Config file for the ScanMan package
# It defines the following variables:

#  ScanMan_BIN_DIR                      - location of ScanMan
#  ScanMan_VERSION                      - version of the binary
#  ScanMan_ConfigPackageLocation        - location of the configuration files (= ScanMan_DIR)
# Imported Targets:
#  ScanMan_BINS                         - all executables
#  ScanMan::ScanMan                     - executable target


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was ScanManConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../" ABSOLUTE)

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
set(cfg_hint "The problem is most likely that the generation of ScanManConfig.cmake in ScanMan/CMakeLists.txt is not up-to-date.")
set(ScanMan_VERSION "0.1.0")
message("-- Found ScanMan, Version ${ScanMan_VERSION}")
if(1)
	set_and_check(ScanMan_BIN_DIR "${PACKAGE_PREFIX_DIR}/Bin/bin")
	set_and_check(ScanMan_ConfigPackageLocation "${PACKAGE_PREFIX_DIR}/Bin/cmake")
	if (NOT PACKAGE_PREFIX_DIR STREQUAL "")
		message("-- Using 'ScanMan' from '${PACKAGE_PREFIX_DIR}'")
	endif()
else()
	set(ScanMan_BIN_DIR "")
	# check when find_package is used whether directory exists
	if(NOT EXISTS "${ScanMan_BIN_DIR}")
		message(FATAL_ERROR "ScanMan_BIN_DIR not found: ${ScanMan_BIN_DIR}"
							" ${cfg_hint}")
	endif()
	set_and_check(ScanMan_ConfigPackageLocation "${PACKAGE_PREFIX_DIR}/Bin/cmake")
	message("-- Using 'ScanMan' from a build tree")
	if (NOT ScanMan_BIN_DIR STREQUAL "")
		message("\t binary dir.:         '${ScanMan_BIN_DIR}'")
	endif()
	if (NOT ScanMasterLibs_ConfigPackageLocation STREQUAL "")
		message("\t target import:       '${ScanMan_ConfigPackageLocation}/ScanManTargets.cmake'")
	endif()
	if (NOT ScanMan_ConfigPackageLocation STREQUAL "")
		message("\t package location:    '${ScanMan_ConfigPackageLocation}'")
	endif()
endif()

include(${ScanMan_ConfigPackageLocation}/ScanManTargets.cmake)
# check whether configuration file is still up to date
set(ScanMan_BINS ScanMan::ScanMan )
foreach(t ${ScanMan_BINS})
	if(NOT TARGET ${t})
		message(FATAL_ERROR "Expected target does not exist: '${t}'\n"
			" Note that targets are case sensitive."
			" ${cfg_hint}")
	endif()
endforeach()
