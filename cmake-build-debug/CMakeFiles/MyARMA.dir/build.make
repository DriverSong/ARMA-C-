# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /home/qiujiawei/Clion/clion-2017.3.4/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/qiujiawei/Clion/clion-2017.3.4/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/qiujiawei/CLionProjects/MyARMA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/MyARMA.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MyARMA.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MyARMA.dir/flags.make

CMakeFiles/MyARMA.dir/main.cpp.o: CMakeFiles/MyARMA.dir/flags.make
CMakeFiles/MyARMA.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MyARMA.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MyARMA.dir/main.cpp.o -c /home/qiujiawei/CLionProjects/MyARMA/main.cpp

CMakeFiles/MyARMA.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MyARMA.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qiujiawei/CLionProjects/MyARMA/main.cpp > CMakeFiles/MyARMA.dir/main.cpp.i

CMakeFiles/MyARMA.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MyARMA.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qiujiawei/CLionProjects/MyARMA/main.cpp -o CMakeFiles/MyARMA.dir/main.cpp.s

CMakeFiles/MyARMA.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/MyARMA.dir/main.cpp.o.requires

CMakeFiles/MyARMA.dir/main.cpp.o.provides: CMakeFiles/MyARMA.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/MyARMA.dir/build.make CMakeFiles/MyARMA.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/MyARMA.dir/main.cpp.o.provides

CMakeFiles/MyARMA.dir/main.cpp.o.provides.build: CMakeFiles/MyARMA.dir/main.cpp.o


CMakeFiles/MyARMA.dir/ARMACore.cpp.o: CMakeFiles/MyARMA.dir/flags.make
CMakeFiles/MyARMA.dir/ARMACore.cpp.o: ../ARMACore.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/MyARMA.dir/ARMACore.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MyARMA.dir/ARMACore.cpp.o -c /home/qiujiawei/CLionProjects/MyARMA/ARMACore.cpp

CMakeFiles/MyARMA.dir/ARMACore.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MyARMA.dir/ARMACore.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qiujiawei/CLionProjects/MyARMA/ARMACore.cpp > CMakeFiles/MyARMA.dir/ARMACore.cpp.i

CMakeFiles/MyARMA.dir/ARMACore.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MyARMA.dir/ARMACore.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qiujiawei/CLionProjects/MyARMA/ARMACore.cpp -o CMakeFiles/MyARMA.dir/ARMACore.cpp.s

CMakeFiles/MyARMA.dir/ARMACore.cpp.o.requires:

.PHONY : CMakeFiles/MyARMA.dir/ARMACore.cpp.o.requires

CMakeFiles/MyARMA.dir/ARMACore.cpp.o.provides: CMakeFiles/MyARMA.dir/ARMACore.cpp.o.requires
	$(MAKE) -f CMakeFiles/MyARMA.dir/build.make CMakeFiles/MyARMA.dir/ARMACore.cpp.o.provides.build
.PHONY : CMakeFiles/MyARMA.dir/ARMACore.cpp.o.provides

CMakeFiles/MyARMA.dir/ARMACore.cpp.o.provides.build: CMakeFiles/MyARMA.dir/ARMACore.cpp.o


CMakeFiles/MyARMA.dir/Predict.cpp.o: CMakeFiles/MyARMA.dir/flags.make
CMakeFiles/MyARMA.dir/Predict.cpp.o: ../Predict.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/MyARMA.dir/Predict.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MyARMA.dir/Predict.cpp.o -c /home/qiujiawei/CLionProjects/MyARMA/Predict.cpp

CMakeFiles/MyARMA.dir/Predict.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MyARMA.dir/Predict.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qiujiawei/CLionProjects/MyARMA/Predict.cpp > CMakeFiles/MyARMA.dir/Predict.cpp.i

CMakeFiles/MyARMA.dir/Predict.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MyARMA.dir/Predict.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qiujiawei/CLionProjects/MyARMA/Predict.cpp -o CMakeFiles/MyARMA.dir/Predict.cpp.s

CMakeFiles/MyARMA.dir/Predict.cpp.o.requires:

.PHONY : CMakeFiles/MyARMA.dir/Predict.cpp.o.requires

CMakeFiles/MyARMA.dir/Predict.cpp.o.provides: CMakeFiles/MyARMA.dir/Predict.cpp.o.requires
	$(MAKE) -f CMakeFiles/MyARMA.dir/build.make CMakeFiles/MyARMA.dir/Predict.cpp.o.provides.build
.PHONY : CMakeFiles/MyARMA.dir/Predict.cpp.o.provides

CMakeFiles/MyARMA.dir/Predict.cpp.o.provides.build: CMakeFiles/MyARMA.dir/Predict.cpp.o


CMakeFiles/MyARMA.dir/Control.cpp.o: CMakeFiles/MyARMA.dir/flags.make
CMakeFiles/MyARMA.dir/Control.cpp.o: ../Control.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/MyARMA.dir/Control.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MyARMA.dir/Control.cpp.o -c /home/qiujiawei/CLionProjects/MyARMA/Control.cpp

CMakeFiles/MyARMA.dir/Control.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MyARMA.dir/Control.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qiujiawei/CLionProjects/MyARMA/Control.cpp > CMakeFiles/MyARMA.dir/Control.cpp.i

CMakeFiles/MyARMA.dir/Control.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MyARMA.dir/Control.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qiujiawei/CLionProjects/MyARMA/Control.cpp -o CMakeFiles/MyARMA.dir/Control.cpp.s

CMakeFiles/MyARMA.dir/Control.cpp.o.requires:

.PHONY : CMakeFiles/MyARMA.dir/Control.cpp.o.requires

CMakeFiles/MyARMA.dir/Control.cpp.o.provides: CMakeFiles/MyARMA.dir/Control.cpp.o.requires
	$(MAKE) -f CMakeFiles/MyARMA.dir/build.make CMakeFiles/MyARMA.dir/Control.cpp.o.provides.build
.PHONY : CMakeFiles/MyARMA.dir/Control.cpp.o.provides

CMakeFiles/MyARMA.dir/Control.cpp.o.provides.build: CMakeFiles/MyARMA.dir/Control.cpp.o


# Object files for target MyARMA
MyARMA_OBJECTS = \
"CMakeFiles/MyARMA.dir/main.cpp.o" \
"CMakeFiles/MyARMA.dir/ARMACore.cpp.o" \
"CMakeFiles/MyARMA.dir/Predict.cpp.o" \
"CMakeFiles/MyARMA.dir/Control.cpp.o"

# External object files for target MyARMA
MyARMA_EXTERNAL_OBJECTS =

MyARMA: CMakeFiles/MyARMA.dir/main.cpp.o
MyARMA: CMakeFiles/MyARMA.dir/ARMACore.cpp.o
MyARMA: CMakeFiles/MyARMA.dir/Predict.cpp.o
MyARMA: CMakeFiles/MyARMA.dir/Control.cpp.o
MyARMA: CMakeFiles/MyARMA.dir/build.make
MyARMA: CMakeFiles/MyARMA.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable MyARMA"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MyARMA.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MyARMA.dir/build: MyARMA

.PHONY : CMakeFiles/MyARMA.dir/build

CMakeFiles/MyARMA.dir/requires: CMakeFiles/MyARMA.dir/main.cpp.o.requires
CMakeFiles/MyARMA.dir/requires: CMakeFiles/MyARMA.dir/ARMACore.cpp.o.requires
CMakeFiles/MyARMA.dir/requires: CMakeFiles/MyARMA.dir/Predict.cpp.o.requires
CMakeFiles/MyARMA.dir/requires: CMakeFiles/MyARMA.dir/Control.cpp.o.requires

.PHONY : CMakeFiles/MyARMA.dir/requires

CMakeFiles/MyARMA.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MyARMA.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MyARMA.dir/clean

CMakeFiles/MyARMA.dir/depend:
	cd /home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qiujiawei/CLionProjects/MyARMA /home/qiujiawei/CLionProjects/MyARMA /home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug /home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug /home/qiujiawei/CLionProjects/MyARMA/cmake-build-debug/CMakeFiles/MyARMA.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MyARMA.dir/depend

