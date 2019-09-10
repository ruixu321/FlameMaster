file(REMOVE_RECURSE
  "FlameMan.pdb"
  "FlameMan"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/FlameMan.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
