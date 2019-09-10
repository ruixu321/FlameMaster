#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ListTool::List" for configuration "Debug"
set_property(TARGET ListTool::List APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(ListTool::List PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/lib/libList.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS ListTool::List )
list(APPEND _IMPORT_CHECK_FILES_FOR_ListTool::List "${_IMPORT_PREFIX}/Bin/lib/libList.a" )

# Import target "ListTool::FLReader" for configuration "Debug"
set_property(TARGET ListTool::FLReader APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(ListTool::FLReader PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/lib/libFLReader.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS ListTool::FLReader )
list(APPEND _IMPORT_CHECK_FILES_FOR_ListTool::FLReader "${_IMPORT_PREFIX}/Bin/lib/libFLReader.a" )

# Import target "ListTool::ListTool" for configuration "Debug"
set_property(TARGET ListTool::ListTool APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(ListTool::ListTool PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/bin/ListTool"
  )

list(APPEND _IMPORT_CHECK_TARGETS ListTool::ListTool )
list(APPEND _IMPORT_CHECK_FILES_FOR_ListTool::ListTool "${_IMPORT_PREFIX}/Bin/bin/ListTool" )

# Import target "ListTool::WSS" for configuration "Debug"
set_property(TARGET ListTool::WSS APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(ListTool::WSS PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/Bin/lib/libWSS.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS ListTool::WSS )
list(APPEND _IMPORT_CHECK_FILES_FOR_ListTool::WSS "${_IMPORT_PREFIX}/Bin/lib/libWSS.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
