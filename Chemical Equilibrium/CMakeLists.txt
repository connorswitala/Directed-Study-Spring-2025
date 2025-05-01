cmake_minimum_required(VERSION 3.10)
project(ChemicalEquilibrium)

# Create executable
add_executable(ChemicalEquilibriumSolver
    "ChemicalEquilibrium.cpp"
    "C:/Users/Connor/source/repos/Directed Study Spring 2025/2D FVS Library/2D FVS Library.cpp"
    "C:/Users/Connor/source/repos/Directed Study Spring 2025/Linear Algebra Library/Linear Algebra Library.cpp"
    "C:/Users/Connor/source/repos/Directed Study Spring 2025/GridGenerator/GridGenerator.cpp"
)

# Add include directories
target_include_directories(ChemicalEquilibriumSolver PUBLIC
    "C:/Users/Connor/source/repos/Directed Study Spring 2025/2D FVS Library"
    "C:/Users/Connor/source/repos/Directed Study Spring 2025/Linear Algebra Library"
    "C:/Users/Connor/source/repos/Directed Study Spring 2025/GridGenerator"
)


