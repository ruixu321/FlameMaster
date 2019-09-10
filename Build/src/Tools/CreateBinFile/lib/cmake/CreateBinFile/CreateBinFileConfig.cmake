# - Config file for the CreateBinFile package
# It defines the following variables:

#  CreateBinFile_BIN_DIR                      - location of CreateBinFile
#  CreateBinFile_VERSION                      - version of the binary
#  CreateBinFile_ConfigPackageLocation        - location of the configuration files ( = CreateFileBin_DIR)
# Imported Targets:
#  CreateBinFile_BINS                         - (all) executable targets
#  CreateBinFile::CreateBinFile               - executable target


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was CreateBinFileConfig.cmake.in                            ########

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
set(cfg_hint "The problem is most likely that the generation of CreateBinFileConfig.cmake in CreateBinFile/CMakeLists.txt is not up-to-date.")
set(CreateBinFile_VERSION "0.1.0")
message("-- Found CreateBinFile, Version ${CreateBinFile_VERSION}")
if(0)
	set_and_check(CreateBinFile_BIN_DIR "${PACKAGE_PREFIX_DIR}/Bin/bin")
	set_and_check(CreateBinFile_ConfigPackageLocation "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/CreateBinFile/lib/cmake/CreateBinFile")
	if (NOT PACKAGE_PREFIX_DIR STREQUAL "")
		message("-- Using 'CreateBinFile' from '${PACKAGE_PREFIX_DIR}'")
	endif()
else()
	set(CreateBinFile_BIN_DIR "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/CreateBinFile")
	# check when find_package is used whether directory exists
	if(NOT EXISTS "${CreateBinFile_BIN_DIR}")
		message(FATAL_ERROR "CreateBinFile_BIN_DIR not found: ${CreateBinFile_BIN_DIR}"
							" ${cfg_hint}")
	endif()
	set_and_check(CreateBinFile_ConfigPackageLocation "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/Tools/CreateBinFile/lib/cmake/CreateBinFile")
	message("-- Using 'CreateBinFile' from a build tree")
	if (NOT CreateBinFile_BIN_DIR STREQUAL "")
		message("\t binary dir.:       '${CreateBinFile_BIN_DIR}'")
	endif()
	if (NOT CreateBinFile_ConfigPackageLocation STREQUAL "")
		message("\t target import:     '${CreateBinFile_ConfigPackageLocation}/CreateBinFileTargets.cmake'")
	endif()
	if (NOT CreateBinFile_ConfigPackageLocation STREQUAL "")
		message("\t package location:  '${CreateBinFile_ConfigPackageLocation}'")
	endif()
endif()

include(${CreateBinFile_ConfigPackageLocation}/CreateBinFileTargets.cmake)
# check whether configuration file is still up to date
set(CreateBinFile_BINS CreateBinFile::CreateBinFile)
foreach(t ${CreateBinFile_BINS})
	if(NOT TARGET ${t})
		message(FATAL_ERROR "Expected target does not exist: '${t}'\n"
			" Note that targets are case sensitive."
			" ${cfg_hint}")
	endif()
endforeach()
