add_executable(p2main
    main.cpp
)

target_include_directories(p2main PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(p2main PUBLIC
    HydroForest
)