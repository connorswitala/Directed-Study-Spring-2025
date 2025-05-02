cmake_minimum_required(VERSION 3.10)

# Define the project
project(DirectedStudy)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:500000000")  # 8 MB stack size

# Add subdirectories for each library and executable
add_subdirectory("Linear Algebra Library")
add_subdirectory("GridGenerator")
add_subdirectory("2D FVS Library")
add_subdirectory("DPLR Finite Volume")
