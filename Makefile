# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/danielconde/Desktop/B3/B3a

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/danielconde/Desktop/B3/B3a

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components

.PHONY : list_install_components/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/danielconde/Desktop/B3/B3a/CMakeFiles /home/danielconde/Desktop/B3/B3a/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/danielconde/Desktop/B3/B3a/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named B3a

# Build rule for target.
B3a: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 B3a
.PHONY : B3a

# fast build rule for target.
B3a/fast:
	$(MAKE) -f CMakeFiles/B3a.dir/build.make CMakeFiles/B3a.dir/build
.PHONY : B3a/fast

#=============================================================================
# Target rules for targets named exampleB3a

# Build rule for target.
exampleB3a: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 exampleB3a
.PHONY : exampleB3a

# fast build rule for target.
exampleB3a/fast:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/build
.PHONY : exampleB3a/fast

exampleB3a.o: exampleB3a.cc.o

.PHONY : exampleB3a.o

# target to build an object file
exampleB3a.cc.o:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/exampleB3a.cc.o
.PHONY : exampleB3a.cc.o

exampleB3a.i: exampleB3a.cc.i

.PHONY : exampleB3a.i

# target to preprocess a source file
exampleB3a.cc.i:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/exampleB3a.cc.i
.PHONY : exampleB3a.cc.i

exampleB3a.s: exampleB3a.cc.s

.PHONY : exampleB3a.s

# target to generate assembly for a file
exampleB3a.cc.s:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/exampleB3a.cc.s
.PHONY : exampleB3a.cc.s

src/B3DetectorConstruction.o: src/B3DetectorConstruction.cc.o

.PHONY : src/B3DetectorConstruction.o

# target to build an object file
src/B3DetectorConstruction.cc.o:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o
.PHONY : src/B3DetectorConstruction.cc.o

src/B3DetectorConstruction.i: src/B3DetectorConstruction.cc.i

.PHONY : src/B3DetectorConstruction.i

# target to preprocess a source file
src/B3DetectorConstruction.cc.i:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.i
.PHONY : src/B3DetectorConstruction.cc.i

src/B3DetectorConstruction.s: src/B3DetectorConstruction.cc.s

.PHONY : src/B3DetectorConstruction.s

# target to generate assembly for a file
src/B3DetectorConstruction.cc.s:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.s
.PHONY : src/B3DetectorConstruction.cc.s

src/B3PhysicsList.o: src/B3PhysicsList.cc.o

.PHONY : src/B3PhysicsList.o

# target to build an object file
src/B3PhysicsList.cc.o:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o
.PHONY : src/B3PhysicsList.cc.o

src/B3PhysicsList.i: src/B3PhysicsList.cc.i

.PHONY : src/B3PhysicsList.i

# target to preprocess a source file
src/B3PhysicsList.cc.i:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.i
.PHONY : src/B3PhysicsList.cc.i

src/B3PhysicsList.s: src/B3PhysicsList.cc.s

.PHONY : src/B3PhysicsList.s

# target to generate assembly for a file
src/B3PhysicsList.cc.s:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.s
.PHONY : src/B3PhysicsList.cc.s

src/B3PrimaryGeneratorAction.o: src/B3PrimaryGeneratorAction.cc.o

.PHONY : src/B3PrimaryGeneratorAction.o

# target to build an object file
src/B3PrimaryGeneratorAction.cc.o:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o
.PHONY : src/B3PrimaryGeneratorAction.cc.o

src/B3PrimaryGeneratorAction.i: src/B3PrimaryGeneratorAction.cc.i

.PHONY : src/B3PrimaryGeneratorAction.i

# target to preprocess a source file
src/B3PrimaryGeneratorAction.cc.i:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.i
.PHONY : src/B3PrimaryGeneratorAction.cc.i

src/B3PrimaryGeneratorAction.s: src/B3PrimaryGeneratorAction.cc.s

.PHONY : src/B3PrimaryGeneratorAction.s

# target to generate assembly for a file
src/B3PrimaryGeneratorAction.cc.s:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.s
.PHONY : src/B3PrimaryGeneratorAction.cc.s

src/B3StackingAction.o: src/B3StackingAction.cc.o

.PHONY : src/B3StackingAction.o

# target to build an object file
src/B3StackingAction.cc.o:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o
.PHONY : src/B3StackingAction.cc.o

src/B3StackingAction.i: src/B3StackingAction.cc.i

.PHONY : src/B3StackingAction.i

# target to preprocess a source file
src/B3StackingAction.cc.i:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.i
.PHONY : src/B3StackingAction.cc.i

src/B3StackingAction.s: src/B3StackingAction.cc.s

.PHONY : src/B3StackingAction.s

# target to generate assembly for a file
src/B3StackingAction.cc.s:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.s
.PHONY : src/B3StackingAction.cc.s

src/B3aActionInitialization.o: src/B3aActionInitialization.cc.o

.PHONY : src/B3aActionInitialization.o

# target to build an object file
src/B3aActionInitialization.cc.o:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o
.PHONY : src/B3aActionInitialization.cc.o

src/B3aActionInitialization.i: src/B3aActionInitialization.cc.i

.PHONY : src/B3aActionInitialization.i

# target to preprocess a source file
src/B3aActionInitialization.cc.i:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.i
.PHONY : src/B3aActionInitialization.cc.i

src/B3aActionInitialization.s: src/B3aActionInitialization.cc.s

.PHONY : src/B3aActionInitialization.s

# target to generate assembly for a file
src/B3aActionInitialization.cc.s:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.s
.PHONY : src/B3aActionInitialization.cc.s

src/B3aEventAction.o: src/B3aEventAction.cc.o

.PHONY : src/B3aEventAction.o

# target to build an object file
src/B3aEventAction.cc.o:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o
.PHONY : src/B3aEventAction.cc.o

src/B3aEventAction.i: src/B3aEventAction.cc.i

.PHONY : src/B3aEventAction.i

# target to preprocess a source file
src/B3aEventAction.cc.i:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.i
.PHONY : src/B3aEventAction.cc.i

src/B3aEventAction.s: src/B3aEventAction.cc.s

.PHONY : src/B3aEventAction.s

# target to generate assembly for a file
src/B3aEventAction.cc.s:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.s
.PHONY : src/B3aEventAction.cc.s

src/B3aRunAction.o: src/B3aRunAction.cc.o

.PHONY : src/B3aRunAction.o

# target to build an object file
src/B3aRunAction.cc.o:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o
.PHONY : src/B3aRunAction.cc.o

src/B3aRunAction.i: src/B3aRunAction.cc.i

.PHONY : src/B3aRunAction.i

# target to preprocess a source file
src/B3aRunAction.cc.i:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.i
.PHONY : src/B3aRunAction.cc.i

src/B3aRunAction.s: src/B3aRunAction.cc.s

.PHONY : src/B3aRunAction.s

# target to generate assembly for a file
src/B3aRunAction.cc.s:
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.s
.PHONY : src/B3aRunAction.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... install/strip"
	@echo "... install/local"
	@echo "... B3a"
	@echo "... rebuild_cache"
	@echo "... list_install_components"
	@echo "... exampleB3a"
	@echo "... edit_cache"
	@echo "... install"
	@echo "... exampleB3a.o"
	@echo "... exampleB3a.i"
	@echo "... exampleB3a.s"
	@echo "... src/B3DetectorConstruction.o"
	@echo "... src/B3DetectorConstruction.i"
	@echo "... src/B3DetectorConstruction.s"
	@echo "... src/B3PhysicsList.o"
	@echo "... src/B3PhysicsList.i"
	@echo "... src/B3PhysicsList.s"
	@echo "... src/B3PrimaryGeneratorAction.o"
	@echo "... src/B3PrimaryGeneratorAction.i"
	@echo "... src/B3PrimaryGeneratorAction.s"
	@echo "... src/B3StackingAction.o"
	@echo "... src/B3StackingAction.i"
	@echo "... src/B3StackingAction.s"
	@echo "... src/B3aActionInitialization.o"
	@echo "... src/B3aActionInitialization.i"
	@echo "... src/B3aActionInitialization.s"
	@echo "... src/B3aEventAction.o"
	@echo "... src/B3aEventAction.i"
	@echo "... src/B3aEventAction.s"
	@echo "... src/B3aRunAction.o"
	@echo "... src/B3aRunAction.i"
	@echo "... src/B3aRunAction.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
