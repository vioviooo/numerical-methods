# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/vioviooo/Desktop/numerical-methods

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/vioviooo/Desktop/numerical-methods/build

# Include any dependencies generated for this target.
include 3.4/CMakeFiles/task_3_4.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include 3.4/CMakeFiles/task_3_4.dir/compiler_depend.make

# Include the progress variables for this target.
include 3.4/CMakeFiles/task_3_4.dir/progress.make

# Include the compile flags for this target's objects.
include 3.4/CMakeFiles/task_3_4.dir/flags.make

3.4/CMakeFiles/task_3_4.dir/main.cpp.o: 3.4/CMakeFiles/task_3_4.dir/flags.make
3.4/CMakeFiles/task_3_4.dir/main.cpp.o: /Users/vioviooo/Desktop/numerical-methods/3.4/main.cpp
3.4/CMakeFiles/task_3_4.dir/main.cpp.o: 3.4/CMakeFiles/task_3_4.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/vioviooo/Desktop/numerical-methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object 3.4/CMakeFiles/task_3_4.dir/main.cpp.o"
	cd /Users/vioviooo/Desktop/numerical-methods/build/3.4 && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT 3.4/CMakeFiles/task_3_4.dir/main.cpp.o -MF CMakeFiles/task_3_4.dir/main.cpp.o.d -o CMakeFiles/task_3_4.dir/main.cpp.o -c /Users/vioviooo/Desktop/numerical-methods/3.4/main.cpp

3.4/CMakeFiles/task_3_4.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/task_3_4.dir/main.cpp.i"
	cd /Users/vioviooo/Desktop/numerical-methods/build/3.4 && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/vioviooo/Desktop/numerical-methods/3.4/main.cpp > CMakeFiles/task_3_4.dir/main.cpp.i

3.4/CMakeFiles/task_3_4.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/task_3_4.dir/main.cpp.s"
	cd /Users/vioviooo/Desktop/numerical-methods/build/3.4 && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/vioviooo/Desktop/numerical-methods/3.4/main.cpp -o CMakeFiles/task_3_4.dir/main.cpp.s

# Object files for target task_3_4
task_3_4_OBJECTS = \
"CMakeFiles/task_3_4.dir/main.cpp.o"

# External object files for target task_3_4
task_3_4_EXTERNAL_OBJECTS =

3.4/task_3_4: 3.4/CMakeFiles/task_3_4.dir/main.cpp.o
3.4/task_3_4: 3.4/CMakeFiles/task_3_4.dir/build.make
3.4/task_3_4: 3.4/CMakeFiles/task_3_4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/vioviooo/Desktop/numerical-methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable task_3_4"
	cd /Users/vioviooo/Desktop/numerical-methods/build/3.4 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/task_3_4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
3.4/CMakeFiles/task_3_4.dir/build: 3.4/task_3_4
.PHONY : 3.4/CMakeFiles/task_3_4.dir/build

3.4/CMakeFiles/task_3_4.dir/clean:
	cd /Users/vioviooo/Desktop/numerical-methods/build/3.4 && $(CMAKE_COMMAND) -P CMakeFiles/task_3_4.dir/cmake_clean.cmake
.PHONY : 3.4/CMakeFiles/task_3_4.dir/clean

3.4/CMakeFiles/task_3_4.dir/depend:
	cd /Users/vioviooo/Desktop/numerical-methods/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/vioviooo/Desktop/numerical-methods /Users/vioviooo/Desktop/numerical-methods/3.4 /Users/vioviooo/Desktop/numerical-methods/build /Users/vioviooo/Desktop/numerical-methods/build/3.4 /Users/vioviooo/Desktop/numerical-methods/build/3.4/CMakeFiles/task_3_4.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : 3.4/CMakeFiles/task_3_4.dir/depend

