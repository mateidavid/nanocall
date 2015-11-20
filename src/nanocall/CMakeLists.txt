message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(nanocall
    nanocall.cpp
    Builtin_Model.cpp
    )
target_link_libraries(nanocall
    version
    libhdf5
    ${CMAKE_DL_LIBS}
    ${ZLIB_LIBRARIES}
    )

if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_executable(compute-state-transitions compute-state-transitions.cpp)

    add_executable(compute-scaled-pore-model compute-scaled-pore-model.cpp)
    target_link_libraries(compute-scaled-pore-model libhdf5 ${CMAKE_DL_LIBS} ${ZLIB_LIBRARIES})

    add_executable(run-fwbw run-fwbw.cpp)
    target_link_libraries(run-fwbw ${ZLIB_LIBRARIES})

    add_executable(run-viterbi run-viterbi.cpp)
    target_link_libraries(run-viterbi ${ZLIB_LIBRARIES})

    add_executable(list-directory list-directory.cpp)
endif()
