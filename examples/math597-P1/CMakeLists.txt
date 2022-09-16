add_executable(p1main
    main.cpp
)

target_include_directories(p1main PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(p1main PUBLIC
    HydroForest
)