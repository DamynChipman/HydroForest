add_executable(p4main
    main.cpp
)

target_include_directories(p4main PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(p4main PUBLIC
    HydroForest
)