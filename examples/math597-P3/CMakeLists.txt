add_executable(p3main
    main.cpp
)

target_include_directories(p3main PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(p3main PUBLIC
    HydroForest
)