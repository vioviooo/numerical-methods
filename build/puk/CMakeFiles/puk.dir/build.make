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
include puk/CMakeFiles/puk.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include puk/CMakeFiles/puk.dir/compiler_depend.make

# Include the progress variables for this target.
include puk/CMakeFiles/puk.dir/progress.make

# Include the compile flags for this target's objects.
include puk/CMakeFiles/puk.dir/flags.make

puk/CMakeFiles/puk.dir/main.cpp.o: puk/CMakeFiles/puk.dir/flags.make
puk/CMakeFiles/puk.dir/main.cpp.o: /Users/vioviooo/Desktop/numerical-methods/puk/main.cpp
puk/CMakeFiles/puk.dir/main.cpp.o: puk/CMakeFiles/puk.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/vioviooo/Desktop/numerical-methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object puk/CMakeFiles/puk.dir/main.cpp.o"
	cd /Users/vioviooo/Desktop/numerical-methods/build/puk && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT puk/CMakeFiles/puk.dir/main.cpp.o -MF CMakeFiles/puk.dir/main.cpp.o.d -o CMakeFiles/puk.dir/main.cpp.o -c /Users/vioviooo/Desktop/numerical-methods/puk/main.cpp

puk/CMakeFiles/puk.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/puk.dir/main.cpp.i"
	cd /Users/vioviooo/Desktop/numerical-methods/build/puk && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/vioviooo/Desktop/numerical-methods/puk/main.cpp > CMakeFiles/puk.dir/main.cpp.i

puk/CMakeFiles/puk.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/puk.dir/main.cpp.s"
	cd /Users/vioviooo/Desktop/numerical-methods/build/puk && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/vioviooo/Desktop/numerical-methods/puk/main.cpp -o CMakeFiles/puk.dir/main.cpp.s

# Object files for target puk
puk_OBJECTS = \
"CMakeFiles/puk.dir/main.cpp.o"

# External object files for target puk
puk_EXTERNAL_OBJECTS =

puk/puk: puk/CMakeFiles/puk.dir/main.cpp.o
puk/puk: puk/CMakeFiles/puk.dir/build.make
puk/puk: puk/CMakeFiles/puk.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/vioviooo/Desktop/numerical-methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable puk"
	cd /Users/vioviooo/Desktop/numerical-methods/build/puk && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/puk.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
puk/CMakeFiles/puk.dir/build: puk/puk
.PHONY : puk/CMakeFiles/puk.dir/build

puk/CMakeFiles/puk.dir/clean:
	cd /Users/vioviooo/Desktop/numerical-methods/build/puk && $(CMAKE_COMMAND) -P CMakeFiles/puk.dir/cmake_clean.cmake
.PHONY : puk/CMakeFiles/puk.dir/clean

puk/CMakeFiles/puk.dir/depend:
	cd /Users/vioviooo/Desktop/numerical-methods/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/vioviooo/Desktop/numerical-methods /Users/vioviooo/Desktop/numerical-methods/puk /Users/vioviooo/Desktop/numerical-methods/build /Users/vioviooo/Desktop/numerical-methods/build/puk /Users/vioviooo/Desktop/numerical-methods/build/puk/CMakeFiles/puk.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : puk/CMakeFiles/puk.dir/depend

