# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/chronis/JoinsRevisited

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/chronis/JoinsRevisited

# Include any dependencies generated for this target.
include datagen/CMakeFiles/zipf.dir/depend.make

# Include the progress variables for this target.
include datagen/CMakeFiles/zipf.dir/progress.make

# Include the compile flags for this target's objects.
include datagen/CMakeFiles/zipf.dir/flags.make

datagen/CMakeFiles/zipf.dir/zipf.cpp.o: datagen/CMakeFiles/zipf.dir/flags.make
datagen/CMakeFiles/zipf.dir/zipf.cpp.o: datagen/zipf.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/chronis/JoinsRevisited/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object datagen/CMakeFiles/zipf.dir/zipf.cpp.o"
	cd /Users/chronis/JoinsRevisited/datagen && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/zipf.dir/zipf.cpp.o -c /Users/chronis/JoinsRevisited/datagen/zipf.cpp

datagen/CMakeFiles/zipf.dir/zipf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/zipf.dir/zipf.cpp.i"
	cd /Users/chronis/JoinsRevisited/datagen && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/chronis/JoinsRevisited/datagen/zipf.cpp > CMakeFiles/zipf.dir/zipf.cpp.i

datagen/CMakeFiles/zipf.dir/zipf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/zipf.dir/zipf.cpp.s"
	cd /Users/chronis/JoinsRevisited/datagen && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/chronis/JoinsRevisited/datagen/zipf.cpp -o CMakeFiles/zipf.dir/zipf.cpp.s

# Object files for target zipf
zipf_OBJECTS = \
"CMakeFiles/zipf.dir/zipf.cpp.o"

# External object files for target zipf
zipf_EXTERNAL_OBJECTS =

datagen/libzipf.a: datagen/CMakeFiles/zipf.dir/zipf.cpp.o
datagen/libzipf.a: datagen/CMakeFiles/zipf.dir/build.make
datagen/libzipf.a: datagen/CMakeFiles/zipf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/chronis/JoinsRevisited/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libzipf.a"
	cd /Users/chronis/JoinsRevisited/datagen && $(CMAKE_COMMAND) -P CMakeFiles/zipf.dir/cmake_clean_target.cmake
	cd /Users/chronis/JoinsRevisited/datagen && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zipf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
datagen/CMakeFiles/zipf.dir/build: datagen/libzipf.a

.PHONY : datagen/CMakeFiles/zipf.dir/build

datagen/CMakeFiles/zipf.dir/clean:
	cd /Users/chronis/JoinsRevisited/datagen && $(CMAKE_COMMAND) -P CMakeFiles/zipf.dir/cmake_clean.cmake
.PHONY : datagen/CMakeFiles/zipf.dir/clean

datagen/CMakeFiles/zipf.dir/depend:
	cd /Users/chronis/JoinsRevisited && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/chronis/JoinsRevisited /Users/chronis/JoinsRevisited/datagen /Users/chronis/JoinsRevisited /Users/chronis/JoinsRevisited/datagen /Users/chronis/JoinsRevisited/datagen/CMakeFiles/zipf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : datagen/CMakeFiles/zipf.dir/depend

