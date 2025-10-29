cmake_minimum_required(VERSION 3.10)
project(SystemEquation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

set(SOURCES
    src/SystemEquation.cpp
    src/Matrix.cpp
    src/main.cpp
    src/solver/SORSolver.cpp
)

set(CORE_SOURCES
    src/SystemEquation.cpp
    src/solver/SORSolver.cpp
    src/Matrix.cpp
)

# Header files
set(HEADERS
    include/SystemEquation.h
    include/solver/LinearSolverBase.h
    include/solver/SORSolver.h
    include/Matrix.h
)

# Create executable
add_executable(system_equation ${SOURCES} ${HEADERS})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(system_equation PRIVATE -Wall -Wextra -O3)
endif()


