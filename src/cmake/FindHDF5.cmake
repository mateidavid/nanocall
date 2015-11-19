# Find HDF5 (https://www.hdfgroup.org/)
# Uses hint:
#   HDF5_ROOT
# Sets:
#   HDF5_FOUND
#   HDF5_INCLUDE_DIR
#   HDF5_LIBRARY
# Saves:
#   HDF5_ROOT
#   HDF5_INCLUDE_DIR_CACHED
#   HDF5_LIBRARY_CACHED

if(NOT HDF5_INCLUDE_DIR_CACHED OR NOT HDF5_LIBRARY_CACHED)
    # save HDF_ROOT
    set(HDF5_ROOT "$ENV{HDF5_ROOT}" CACHE PATH "Path to HDF5")

    # find headers
    find_path(HDF5_INCLUDE_DIR_CACHED H5pubconf.h PATHS ${HDF5_ROOT}/include)
    if(HDF5_INCLUDE_DIR_CACHED)
        execute_process(
            COMMAND grep H5_VERSION ${HDF5_INCLUDE_DIR_CACHED}/H5pubconf.h
            COMMAND awk "{print \$3}"
            COMMAND tr -d "\"\n"
            OUTPUT_VARIABLE HDF5_INCLUDE_DIR_VERSION
            )
        message(STATUS "Found HDF5 headers version ${HDF5_INCLUDE_DIR_VERSION} in: ${HDF5_INCLUDE_DIR_CACHED}")
    endif()

    # find library
    find_library(HDF5_LIBRARY_CACHED hdf5 PATHS ${HDF5_ROOT}/lib ${HDF5_ROOT}/lib64)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(HDF5
        REQUIRED_VARS HDF5_INCLUDE_DIR_CACHED HDF5_LIBRARY_CACHED
        VERSION_VAR HDF5_INCLUDE_DIR_VERSION
        #"HDF5 library (https://www.hdfgroup.org/) not found. Specify location with -DHDF5_ROOT=<path>"
        )
    mark_as_advanced(HDF5_INCLUDE_DIR_CACHED HDF5_LIBRARY_CACHED)
else()
    message(STATUS "Using HDF5_INCLUDE_DIR: ${HDF5_INCLUDE_DIR_CACHED}")
    message(STATUS "Using HDF5_LIBRARY: ${HDF5_LIBRARY_CACHED}")
    set(HDF5_FOUND TRUE)
endif()

if(HDF5_FOUND)
    set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR_CACHED})
    set(HDF5_LIBRARY ${HDF5_LIBRARY_CACHED})
endif()
