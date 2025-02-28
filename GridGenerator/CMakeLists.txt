cmake_minimum_required(VERSION 3.10)
project(GridGenerator)

add_library(GridGenerator STATIC
    "GridGenerator.cpp"
)

target_include_directories(GridGenerator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(GridGenerator PRIVATE LinearAlgebraLibrary)
