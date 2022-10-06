cmake_minimum_required(VERSION 3.0.0)   # Selects the minimum version of CMake required to run this file
project(bvh-trace VERSION 0.1.0)          # Here we select the project name and version

# Here we select C++17 with all the standards required and all compiler-specific extensions disabled
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# set(CMAKE_CXX_FLAGS_RELEASE "-flto -ffast-math -O3 -Ofast -ffloat-store -march=native -frename-registers -funroll-loops -fopenmp")

CMAKE_POLICY(SET CMP0074 NEW)

set(SOURCES 
  src/main.cpp
  src/bvh.cpp
  src/bvh_instance.cpp
  src/tlas.cpp

  src/vec_math.cpp
  src/util.cpp
)

# Define the directories in which to search for the included headers
include_directories(
  src/
)

set(CMAKE_DEBUG_POSTFIX _Debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY          ${PROJECT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG    ${PROJECT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

find_package(OpenMP REQUIRED)

# For each example, we add an executable target
# Each target compiles one example source file and the common & vendor source files
# Then we link GLFW with each target
add_executable(${PROJECT_NAME} src/main.cpp ${SOURCES})

target_link_libraries(
  ${PROJECT_NAME} 
  PRIVATE OpenMP::OpenMP_CXX
  -lraylib 
  -lpthread 
  -lGL 
  -lm 
  -lrt 
  -lX11 
  -ldl
  -msse2
)
