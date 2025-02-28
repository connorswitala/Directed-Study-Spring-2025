cmake_minimum_required(VERSION 3.10)
project(DPLRFiniteVolume)

add_executable(DPLRFiniteVolume "DPLR Finite Volume.cpp")

target_include_directories(DPLRFiniteVolume PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DPLRFiniteVolume PRIVATE LinearAlgebraLibrary GridGenerator FVSLibrary)
