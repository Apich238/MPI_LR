project(Boltzman)

set(SOURCE_FILES boltzmann.c)
find_package(MPI REQUIRED)

add_executable(Boltzman ${SOURCE_FILES})
target_link_libraries(Boltzman ${MPI_LIBRARIES})

if (MPI_COMPILE_FLAGS)
    set_target_properties(Boltzman PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(Boltzman PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()