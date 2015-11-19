message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(version STATIC version.cpp)
add_dependencies(version package_version)
