# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.11.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.11.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build

# Include any dependencies generated for this target.
include src/libraries/Newton/CMakeFiles/Newton.dir/depend.make

# Include the progress variables for this target.
include src/libraries/Newton/CMakeFiles/Newton.dir/progress.make

# Include the compile flags for this target's objects.
include src/libraries/Newton/CMakeFiles/Newton.dir/flags.make

src/libraries/Newton/CMakeFiles/Newton.dir/AdaptiveGrid.C.o: src/libraries/Newton/CMakeFiles/Newton.dir/flags.make
src/libraries/Newton/CMakeFiles/Newton.dir/AdaptiveGrid.C.o: /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/AdaptiveGrid.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/libraries/Newton/CMakeFiles/Newton.dir/AdaptiveGrid.C.o"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Newton.dir/AdaptiveGrid.C.o -c /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/AdaptiveGrid.C

src/libraries/Newton/CMakeFiles/Newton.dir/AdaptiveGrid.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Newton.dir/AdaptiveGrid.C.i"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/AdaptiveGrid.C > CMakeFiles/Newton.dir/AdaptiveGrid.C.i

src/libraries/Newton/CMakeFiles/Newton.dir/AdaptiveGrid.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Newton.dir/AdaptiveGrid.C.s"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/AdaptiveGrid.C -o CMakeFiles/Newton.dir/AdaptiveGrid.C.s

src/libraries/Newton/CMakeFiles/Newton.dir/Newton.C.o: src/libraries/Newton/CMakeFiles/Newton.dir/flags.make
src/libraries/Newton/CMakeFiles/Newton.dir/Newton.C.o: /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/Newton.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/libraries/Newton/CMakeFiles/Newton.dir/Newton.C.o"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Newton.dir/Newton.C.o -c /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/Newton.C

src/libraries/Newton/CMakeFiles/Newton.dir/Newton.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Newton.dir/Newton.C.i"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/Newton.C > CMakeFiles/Newton.dir/Newton.C.i

src/libraries/Newton/CMakeFiles/Newton.dir/Newton.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Newton.dir/Newton.C.s"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/Newton.C -o CMakeFiles/Newton.dir/Newton.C.s

src/libraries/Newton/CMakeFiles/Newton.dir/NewtonUt.C.o: src/libraries/Newton/CMakeFiles/Newton.dir/flags.make
src/libraries/Newton/CMakeFiles/Newton.dir/NewtonUt.C.o: /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/NewtonUt.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/libraries/Newton/CMakeFiles/Newton.dir/NewtonUt.C.o"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Newton.dir/NewtonUt.C.o -c /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/NewtonUt.C

src/libraries/Newton/CMakeFiles/Newton.dir/NewtonUt.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Newton.dir/NewtonUt.C.i"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/NewtonUt.C > CMakeFiles/Newton.dir/NewtonUt.C.i

src/libraries/Newton/CMakeFiles/Newton.dir/NewtonUt.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Newton.dir/NewtonUt.C.s"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton/NewtonUt.C -o CMakeFiles/Newton.dir/NewtonUt.C.s

# Object files for target Newton
Newton_OBJECTS = \
"CMakeFiles/Newton.dir/AdaptiveGrid.C.o" \
"CMakeFiles/Newton.dir/Newton.C.o" \
"CMakeFiles/Newton.dir/NewtonUt.C.o"

# External object files for target Newton
Newton_EXTERNAL_OBJECTS =

src/libraries/Newton/libNewton.a: src/libraries/Newton/CMakeFiles/Newton.dir/AdaptiveGrid.C.o
src/libraries/Newton/libNewton.a: src/libraries/Newton/CMakeFiles/Newton.dir/Newton.C.o
src/libraries/Newton/libNewton.a: src/libraries/Newton/CMakeFiles/Newton.dir/NewtonUt.C.o
src/libraries/Newton/libNewton.a: src/libraries/Newton/CMakeFiles/Newton.dir/build.make
src/libraries/Newton/libNewton.a: src/libraries/Newton/CMakeFiles/Newton.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libNewton.a"
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && $(CMAKE_COMMAND) -P CMakeFiles/Newton.dir/cmake_clean_target.cmake
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Newton.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/libraries/Newton/CMakeFiles/Newton.dir/build: src/libraries/Newton/libNewton.a

.PHONY : src/libraries/Newton/CMakeFiles/Newton.dir/build

src/libraries/Newton/CMakeFiles/Newton.dir/clean:
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton && $(CMAKE_COMMAND) -P CMakeFiles/Newton.dir/cmake_clean.cmake
.PHONY : src/libraries/Newton/CMakeFiles/Newton.dir/clean

src/libraries/Newton/CMakeFiles/Newton.dir/depend:
	cd /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/src/libraries/Newton /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton /Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/libraries/Newton/CMakeFiles/Newton.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/libraries/Newton/CMakeFiles/Newton.dir/depend

