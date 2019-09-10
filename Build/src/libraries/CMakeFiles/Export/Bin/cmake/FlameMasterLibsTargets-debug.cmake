#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "FMLibs::Alligator" for configuration "Debug"
set_property(TARGET FMLibs::Alligator APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(FMLibs::Alligator PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/lib/libAlligator.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS FMLibs::Alligator )
list(APPEND _IMPORT_CHECK_FILES_FOR_FMLibs::Alligator "${_IMPORT_PREFIX}/Bin/lib/libAlligator.a" )

# Import target "FMLibs::ArrayMan" for configuration "Debug"
set_property(TARGET FMLibs::ArrayMan APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(FMLibs::ArrayMan PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/lib/libArrayMan.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS FMLibs::ArrayMan )
list(APPEND _IMPORT_CHECK_FILES_FOR_FMLibs::ArrayMan "${_IMPORT_PREFIX}/Bin/lib/libArrayMan.a" )

# Import target "FMLibs::FM" for configuration "Debug"
set_property(TARGET FMLibs::FM APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(FMLibs::FM PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C;CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/lib/libFM.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS FMLibs::FM )
list(APPEND _IMPORT_CHECK_FILES_FOR_FMLibs::FM "${_IMPORT_PREFIX}/Bin/lib/libFM.a" )

# Import target "FMLibs::Newton" for configuration "Debug"
set_property(TARGET FMLibs::Newton APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(FMLibs::Newton PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/lib/libNewton.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS FMLibs::Newton )
list(APPEND _IMPORT_CHECK_FILES_FOR_FMLibs::Newton "${_IMPORT_PREFIX}/Bin/lib/libNewton.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
