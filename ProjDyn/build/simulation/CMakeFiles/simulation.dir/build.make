# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.15.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.15.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/huajian/Geo3D/project/ProjDyn

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/huajian/Geo3D/project/ProjDyn/build

# Include any dependencies generated for this target.
include simulation/CMakeFiles/simulation.dir/depend.make

# Include the progress variables for this target.
include simulation/CMakeFiles/simulation.dir/progress.make

# Include the compile flags for this target's objects.
include simulation/CMakeFiles/simulation.dir/flags.make

simulation/CMakeFiles/simulation.dir/main.cpp.o: simulation/CMakeFiles/simulation.dir/flags.make
simulation/CMakeFiles/simulation.dir/main.cpp.o: ../simulation/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huajian/Geo3D/project/ProjDyn/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object simulation/CMakeFiles/simulation.dir/main.cpp.o"
	cd /Users/huajian/Geo3D/project/ProjDyn/build/simulation && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simulation.dir/main.cpp.o -c /Users/huajian/Geo3D/project/ProjDyn/simulation/main.cpp

simulation/CMakeFiles/simulation.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulation.dir/main.cpp.i"
	cd /Users/huajian/Geo3D/project/ProjDyn/build/simulation && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huajian/Geo3D/project/ProjDyn/simulation/main.cpp > CMakeFiles/simulation.dir/main.cpp.i

simulation/CMakeFiles/simulation.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulation.dir/main.cpp.s"
	cd /Users/huajian/Geo3D/project/ProjDyn/build/simulation && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huajian/Geo3D/project/ProjDyn/simulation/main.cpp -o CMakeFiles/simulation.dir/main.cpp.s

# Object files for target simulation
simulation_OBJECTS = \
"CMakeFiles/simulation.dir/main.cpp.o"

# External object files for target simulation
simulation_EXTERNAL_OBJECTS =

simulation/simulation: simulation/CMakeFiles/simulation.dir/main.cpp.o
simulation/simulation: simulation/CMakeFiles/simulation.dir/build.make
simulation/simulation: externals/tetgen/libtetgen.1.0.dylib
simulation/simulation: nanogui/libnanogui.a
simulation/simulation: externals/surface_mesh/libsurface_mesh.1.0.dylib
simulation/simulation: simulation/CMakeFiles/simulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/huajian/Geo3D/project/ProjDyn/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable simulation"
	cd /Users/huajian/Geo3D/project/ProjDyn/build/simulation && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
simulation/CMakeFiles/simulation.dir/build: simulation/simulation

.PHONY : simulation/CMakeFiles/simulation.dir/build

simulation/CMakeFiles/simulation.dir/clean:
	cd /Users/huajian/Geo3D/project/ProjDyn/build/simulation && $(CMAKE_COMMAND) -P CMakeFiles/simulation.dir/cmake_clean.cmake
.PHONY : simulation/CMakeFiles/simulation.dir/clean

simulation/CMakeFiles/simulation.dir/depend:
	cd /Users/huajian/Geo3D/project/ProjDyn/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/huajian/Geo3D/project/ProjDyn /Users/huajian/Geo3D/project/ProjDyn/simulation /Users/huajian/Geo3D/project/ProjDyn/build /Users/huajian/Geo3D/project/ProjDyn/build/simulation /Users/huajian/Geo3D/project/ProjDyn/build/simulation/CMakeFiles/simulation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : simulation/CMakeFiles/simulation.dir/depend
