message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# set build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: Debug Test Release GProf GProfRel."
       FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
if(NOT ${CMAKE_BUILD_TYPE} AND NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Enable verbose make output.")
    message(STATUS "Enabling verbose make output.")
endif()

# general cmake settings
cmake_minimum_required(VERSION 2.8.12)
project(NANOCALL C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include(CheckLibraryExists)
include(CheckIncludeFileCXX)
include(ExternalProject)

# force out-of-source build
if(${PROJECT_BINARY_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
    message(FATAL_ERROR "In-source build not supported.")
endif()
# check that "nanocall.cpp" exists; if not, assume the source dir is wrong
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/nanocall/nanocall.cpp")
    message(FATAL_ERROR "${PROJECT_SOURCE_DIR}: source files not found")
endif()

# set project-related variables
set(PACKAGE_BUGREPORT "mdavid@oicr.on.ca")
set(PACKAGE_URL "https://github.com/jts/nanocall")
set(PACKAGE "${PROJECT_NAME}")
set(PACKAGE_NAME "${PROJECT_NAME}")
set(PACKAGE_TARNAME "${PROJECT_NAME}")

# directory where to build optional submodules
set(OPTIONAL_SUBMODULE_PREFIX
    ${PROJECT_SOURCE_DIR}/optional
    CACHE INTERNAL "Directory for installing optional submodules")

# header and source directories
set(SUBDIRS
    nanocall version
    CACHE INTERNAL "Subdirectories to descend into")
set(HEADER_SUBDIRS
    builtin_models fast5/src hpptools/include tclap/include
    nanocall version
    CACHE INTERNAL "Subdirectories containing header files")

### Resolve external dependencies
#
# prefer static libraries
#set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so)

# bake-in RPATH to prevent library search problems
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
set(CMAKE_BUILD_WITH_INSTALL_RPATH true)
set(CMAKE_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})

# find zlib
# => ZLIB_INCLUDE_DIRS, ZLIB_LIBRARIES
find_package(ZLIB REQUIRED)

# find HDF5
# => HDF5_INCLUDE_DIRS, HDF5_LIBRARIES
if(BUILD_HDF5)
    message(STATUS "Building HDF5 in: ${OPTIONAL_SUBMODULE_PREFIX}")
    if(UNIX)
        # use local copy if available
        if(EXISTS /tmp/hdf5-1.8.16.tar.bz2)
            set(HDF5_URL /tmp/hdf5-1.8.16.tar.bz2)
        else()
            set(HDF5_URL http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.16/src/hdf5-1.8.16.tar.bz2)
        endif()
        ExternalProject_Add(hdf5
            PREFIX ${OPTIONAL_SUBMODULE_PREFIX}
            URL ${HDF5_URL}
            URL_MD5 79c1593573ebddf734eee8d43ecfe483
            CONFIGURE_COMMAND ./configure --prefix=${OPTIONAL_SUBMODULE_PREFIX} --disable-hl --disable-shared --enable-threadsafe
            BUILD_IN_SOURCE 1
            BUILD_COMMAND make
            INSTALL_COMMAND make install
            )
        set(HDF5_INCLUDE_DIRS ${OPTIONAL_SUBMODULE_PREFIX}/include)
        set(HDF5_LIBRARIES ${OPTIONAL_SUBMODULE_PREFIX}/lib/libhdf5.a)
    elseif(WIN32)
        message(FATAL_ERROR "Building HDF5 not supported on Windows")
    endif()
else()
    find_package(HDF5 1.8.11 REQUIRED)
endif()

# Problem:
#
# By default, if using target_link_libraries(... ${HDF5_LIBRARIES}), cmake
# attempts to use the short name "-lhdf5" as an argument to ld, *without* setting
# "-L" appropriately. It's unclear why, or if this is related to LD_LIBRARY_PATH.
# The problem arises if an older libhdf5 is found instead. The solution below
# uses an IMPORTED library. In this case, cmake will use the full path during
# linking.
#
# Also, setting up RPATH apropriately using cmake is... elusive. If baking-in of
# linker paths is needed, be sure to use LD_RUN_PATH during make- this will be
# honoured by g++.
#
# Refs:
# https://cmake.org/pipermail/cmake/2013-December/056655.html
# https://cmake.org/Wiki/CMake/Tutorials/Exporting_and_Importing_Targets
#
add_library(libhdf5 UNKNOWN IMPORTED)
set_property(TARGET libhdf5 PROPERTY IMPORTED_LOCATION ${HDF5_LIBRARIES})
#get_filename_component(HDF5_LIBRARIES_DIR ${HDF5_LIBRARIES} DIRECTORY)
#set(CMAKE_INSTALL_RPATH ${HDF5_LIBRARIES_DIR})

message(STATUS "HDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5_LIBRARIES=${HDF5_LIBRARIES}")

# message(STATUS "CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH}")
# message(STATUS "CMAKE_SYSTEM_LIBRARY_PATH=${CMAKE_SYSTEM_LIBRARY_PATH}")
# message(STATUS "CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES=${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
# message(STATUS "CMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}")
# message(STATUS "CMAKE_INSTALL_RPATH_USE_LINK_PATH=${CMAKE_INSTALL_RPATH_USE_LINK_PATH}")
# message(STATUS "CMAKE_BUILD_WITH_INSTALL_RPATH=${CMAKE_BUILD_WITH_INSTALL_RPATH}")
# message(STATUS "CMAKE_SKIP_BUILD_RPATH=${CMAKE_SKIP_BUILD_RPATH}")

### Prepare default compilation flags
#
# headers
include_directories(SYSTEM
    ${ZLIB_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    )
include_directories(
    ${PROJECT_BINARY_DIR}
    ${HEADER_SUBDIRS}
    )
get_directory_property(include_directories INCLUDE_DIRECTORIES)
message(STATUS "INCLUDE_DIRECTORIES='${include_directories}'")

### general compile flags
set(EXTRA_FLAGS "-std=c++11 -pthread -Wall -Wextra -pedantic")

# compiler-specific flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -fmax-errors=1")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -ferror-limit=1")
endif()

# consolidate compile flags
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} ${EXTRA_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS}'")

### build-specific compile flags
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -ggdb -fno-inline -fno-eliminate-unused-debug-types")
set(CMAKE_CXX_FLAGS_TEST "-O3 -g3 -fno-eliminate-unused-debug-types")
set(CMAKE_CXX_FLAGS_TEST_O2 "-O2 -g3 -fno-eliminate-unused-debug-types")
set(CMAKE_CXX_FLAGS_TEST_O1 "-O1 -g3 -fno-eliminate-unused-debug-types")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DDISABLE_ASSERTS")
set(CMAKE_CXX_FLAGS_GPROF "-O3 -g3 -pg")
set(CMAKE_CXX_FLAGS_GPROFREL "-O3 -DNDEBUG -DDISABLE_ASSERTS -pg")

# link flags
#set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=$ENV{LD_RUN_PATH} -Wl,--as-needed")
if(APPLE)
	set(CMAKE_EXE_LINKER_FLAGS "-Wl")
else()
	set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed")
endif()
message(STATUS "CMAKE_EXE_LINKER_FLAGS='${CMAKE_EXE_LINKER_FLAGS}'")

### This target updates the package version
#
add_custom_target(package_version
    ${PROJECT_SOURCE_DIR}/get-dir-version -v -d ${PROJECT_SOURCE_DIR}/.. -r ${PROJECT_SOURCE_DIR}/.. -o ${PROJECT_BINARY_DIR}/package_version.h -t PACKAGE_VERSION
    )

### Enable testing targets
#enable_testing()

### Descend into subdirectories
#
foreach(dir ${SUBDIRS})
    add_subdirectory(${dir})
endforeach()
