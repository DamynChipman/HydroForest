cmake_minimum_required(VERSION 3.21)
project("HydroForest")
set(CMAKE_CXX_STANDARD 20)

set(
    P4EST_BUILD
    /Users/damynchipman/packages/p4est/p4est_source_git/build/local
)
set(
    HOMEBREW_DIR
    /opt/homebrew
)

set(
    SOURCES
    main.cpp
    src/Options.cpp
)

add_executable(main ${SOURCES})
target_link_directories(main PUBLIC ${P4EST_BUILD}/lib ${HOMEBREW_DIR}/lib)
target_include_directories(main PUBLIC ${P4EST_BUILD}/include ${HOMEBREW_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(main petsc p4est sc z)