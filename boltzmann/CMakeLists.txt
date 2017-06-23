project(Boltzman)

set(SOURCE_FILES boltzmann.c)
find_package(MPI REQUIRED)

add_executable(MPI_MainElementInColumn ${SOURCE_FILES})
target_link_libraries(MPI_MainElementInColumn ${MPI_LIBRARIES})

if (MPI_COMPILE_FLAGS)
    set_target_properties(MPI_MainElementInColumn PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(MPI_MainElementInColumn PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()