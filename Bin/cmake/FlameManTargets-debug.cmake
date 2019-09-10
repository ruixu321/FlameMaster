#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "FlameMan::FlameMan" for configuration "Debug"
set_property(TARGET FlameMan::FlameMan APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(FlameMan::FlameMan PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/bin/FlameMan"
  )

list(APPEND _IMPORT_CHECK_TARGETS FlameMan::FlameMan )
list(APPEND _IMPORT_CHECK_FILES_FOR_FlameMan::FlameMan "${_IMPORT_PREFIX}/Bin/bin/FlameMan" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
