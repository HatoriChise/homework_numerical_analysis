cmake_minimum_required(VERSION 3.10)
project(CubicSpline)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Source files
set(SOURCES
    src/cubic_spline.cpp
    src/main.cpp
)

# Library sources (without main)
set(CORE_SOURCES
    src/cubic_spline.cpp
)

# Header files
set(HEADERS
    include/cubic_spline.h
)

# Create executable
add_executable(cubic_spline ${SOURCES} ${HEADERS})

# Compiler options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cubic_spline PRIVATE -Wall -Wextra -O3)
endif()

# Create test executable (link only core sources, avoid main.cpp)
set(TEST_SOURCES tests/test_cubic_spline.cpp)
add_executable(test_cubic_spline ${TEST_SOURCES} ${CORE_SOURCES})
target_include_directories(test_cubic_spline PRIVATE include)