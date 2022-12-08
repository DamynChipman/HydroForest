# --=== External Package ===--
FetchContent_GetProperties(matplotlibcpp)

# --=== Create HydroForest Library ===--
add_library(HydroForest)
set_target_properties(HydroForest PROPERTIES EXPORT_NAME HydroForest)

set(HYDRO_FOREST_SOURCE_FILES
    VTKBuilder.cpp
    P4est.cpp
)

# set(HYDRO_FOREST_HEADER_FILES
#     Element1DBase.hpp
#     GenericSingleton.hpp
#     Grid1DBase.hpp
#     HydroForestApp.hpp
#     Logger.hpp
#     Options.hpp
#     UniformGrid1D.hpp
#     VTKBuilder.hpp
#     XMLTree.hpp
# )

target_sources(HydroForest PRIVATE
    ${HYDRO_FOREST_SOURCE_FILES}
)

target_include_directories(HydroForest PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${EXTERNAL_PKGS_INCLUDE}
)

target_link_directories(HydroForest PUBLIC
    ${EXTERNAL_PKGS_LIB}
)

target_link_libraries(HydroForest PUBLIC
    petsc sc p4est python3.9 lapack blas
)

