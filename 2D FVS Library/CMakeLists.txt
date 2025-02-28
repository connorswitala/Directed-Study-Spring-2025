cmake_minimum_required(VERSION 3.10)
project(2DFVSLibrary)

add_library(FVSLibrary STATIC
    "2D FVS Library.cpp"
)

target_include_directories(FVSLibrary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(FVSLibrary PRIVATE LinearAlgebraLibrary GridGenerator)
