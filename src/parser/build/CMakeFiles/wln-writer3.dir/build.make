# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mikey/Workspace/freycompchem/wisswesser/src/parser

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mikey/Workspace/freycompchem/wisswesser/src/parser/build

# Include any dependencies generated for this target.
include CMakeFiles/wln-writer3.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/wln-writer3.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/wln-writer3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/wln-writer3.dir/flags.make

CMakeFiles/wln-writer3.dir/writewln3.cpp.o: CMakeFiles/wln-writer3.dir/flags.make
CMakeFiles/wln-writer3.dir/writewln3.cpp.o: ../writewln3.cpp
CMakeFiles/wln-writer3.dir/writewln3.cpp.o: CMakeFiles/wln-writer3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikey/Workspace/freycompchem/wisswesser/src/parser/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/wln-writer3.dir/writewln3.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/wln-writer3.dir/writewln3.cpp.o -MF CMakeFiles/wln-writer3.dir/writewln3.cpp.o.d -o CMakeFiles/wln-writer3.dir/writewln3.cpp.o -c /home/mikey/Workspace/freycompchem/wisswesser/src/parser/writewln3.cpp

CMakeFiles/wln-writer3.dir/writewln3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wln-writer3.dir/writewln3.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikey/Workspace/freycompchem/wisswesser/src/parser/writewln3.cpp > CMakeFiles/wln-writer3.dir/writewln3.cpp.i

CMakeFiles/wln-writer3.dir/writewln3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wln-writer3.dir/writewln3.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikey/Workspace/freycompchem/wisswesser/src/parser/writewln3.cpp -o CMakeFiles/wln-writer3.dir/writewln3.cpp.s

# Object files for target wln-writer3
wln__writer3_OBJECTS = \
"CMakeFiles/wln-writer3.dir/writewln3.cpp.o"

# External object files for target wln-writer3
wln__writer3_EXTERNAL_OBJECTS =

wln-writer3: CMakeFiles/wln-writer3.dir/writewln3.cpp.o
wln-writer3: CMakeFiles/wln-writer3.dir/build.make
wln-writer3: ../../openbabel/build/lib/libopenbabel.so.7
wln-writer3: CMakeFiles/wln-writer3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mikey/Workspace/freycompchem/wisswesser/src/parser/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable wln-writer3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wln-writer3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/wln-writer3.dir/build: wln-writer3
.PHONY : CMakeFiles/wln-writer3.dir/build

CMakeFiles/wln-writer3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/wln-writer3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/wln-writer3.dir/clean

CMakeFiles/wln-writer3.dir/depend:
	cd /home/mikey/Workspace/freycompchem/wisswesser/src/parser/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mikey/Workspace/freycompchem/wisswesser/src/parser /home/mikey/Workspace/freycompchem/wisswesser/src/parser /home/mikey/Workspace/freycompchem/wisswesser/src/parser/build /home/mikey/Workspace/freycompchem/wisswesser/src/parser/build /home/mikey/Workspace/freycompchem/wisswesser/src/parser/build/CMakeFiles/wln-writer3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/wln-writer3.dir/depend

