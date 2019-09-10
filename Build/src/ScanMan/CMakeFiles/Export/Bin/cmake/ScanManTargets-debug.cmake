#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ScanMan::ScanMan" for configuration "Debug"
set_property(TARGET ScanMan::ScanMan APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(ScanMan::ScanMan PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/bin/ScanMan"
  )

list(APPEND _IMPORT_CHECK_TARGETS ScanMan::ScanMan )
list(APPEND _IMPORT_CHECK_FILES_FOR_ScanMan::ScanMan "${_IMPORT_PREFIX}/Bin/bin/ScanMan" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
