# Find HDF5 (https://www.hdfgroup.org/)
# Uses hint:
#   HDF5_ROOT
# Sets:
#   HDF5_FOUND
#   HDF5_INCLUDE_DIRS
#   HDF5_LIBRARIES
# Saves:
#   HDF5_ROOT
#   HDF5_INCLUDE_DIRS_CACHED
#   HDF5_LIBRARIES_CACHED

if(NOT "${HDF5_ROOT}" STREQUAL "${OLD_HDF5_ROOT}")
    message(STATUS "Detecting HDF5: redetecing with new HDF5_ROOT=${HDF5_ROOT} (OLD_HDF5_ROOT=${OLD_HDF5_ROOT}).")
    unset(HDF5_INCLUDE_DIRS_CACHED CACHE)
    unset(HDF5_LIBRARIES_CACHED CACHE)
else()
    message(STATUS "Detecting HDF5: HDF5_ROOT=${HDF5_ROOT} is not new; using cached paths.")
    message(STATUS "HDF5_INCLUDE_DIRS_CACHED=${HDF5_INCLUDE_DIRS_CACHED}")
    message(STATUS "HDF5_LIBRARIES_CACHED=${HDF5_LIBRARIES_CACHED}")
endif()
set(OLD_HDF5_ROOT ${HDF5_ROOT} CACHE INTERNAL "Last used value of HDF5_ROOT")

# find headers
find_path(HDF5_INCLUDE_DIRS_CACHED H5pubconf.h PATHS ${HDF5_ROOT}/include NO_DEFAULT_PATH)
find_path(HDF5_INCLUDE_DIRS_CACHED H5pubconf.h)
if(HDF5_INCLUDE_DIRS_CACHED)
    execute_process(
        COMMAND grep H5_VERSION ${HDF5_INCLUDE_DIRS_CACHED}/H5pubconf.h
        COMMAND awk "{print \$3}"
        COMMAND tr -d "\"\n"
        OUTPUT_VARIABLE HDF5_INCLUDE_DIRS_VERSION
        )
    message(STATUS "Found HDF5 headers version ${HDF5_INCLUDE_DIRS_VERSION} in: ${HDF5_INCLUDE_DIRS_CACHED}")
endif()

# find library
find_library(HDF5_LIBRARIES_CACHED hdf5 PATHS ${HDF5_ROOT}/lib ${HDF5_ROOT}/lib64 NO_DEFAULT_PATH)
find_library(HDF5_LIBRARIES_CACHED hdf5)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5
    REQUIRED_VARS HDF5_INCLUDE_DIRS_CACHED HDF5_LIBRARIES_CACHED
    VERSION_VAR HDF5_INCLUDE_DIRS_VERSION
    #"HDF5 library (https://www.hdfgroup.org/) not found. Specify location with -DHDF5_ROOT=<path>"
    )
mark_as_advanced(HDF5_INCLUDE_DIRS_CACHED HDF5_LIBRARIES_CACHED)

if(HDF5_FOUND)
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS_CACHED})
    set(HDF5_LIBRARIES ${HDF5_LIBRARIES_CACHED})
endif()
