#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "CreateBinFile::CreateBinFile" for configuration "Debug"
set_property(TARGET CreateBinFile::CreateBinFile APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(CreateBinFile::CreateBinFile PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/bin/CreateBinFile"
  )

list(APPEND _IMPORT_CHECK_TARGETS CreateBinFile::CreateBinFile )
list(APPEND _IMPORT_CHECK_FILES_FOR_CreateBinFile::CreateBinFile "${_IMPORT_PREFIX}/Bin/bin/CreateBinFile" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
