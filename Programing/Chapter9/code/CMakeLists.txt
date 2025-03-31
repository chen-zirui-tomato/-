cmake_minimum_required(VERSION 3.10)
project(NumericalAnalysisChapter9)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen (header-only library)
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(
    ${EIGEN3_INCLUDE_DIR}
    include
    external
)

# Add library
add_library(NumericalAnalysisLib src/Sparse.cpp)

# Add executable
add_executable(P test/P.cpp)
target_link_libraries(P PRIVATE NumericalAnalysisLib)

# Enhanced Python/matplotlib-cpp support
find_package(Python COMPONENTS Interpreter Development)
if(Python_FOUND)
    # Custom NumPy finder
    execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE NUMPY_FOUND
    )
    
    if(NUMPY_FOUND EQUAL 0)
        message(STATUS "Found NumPy at: ${NUMPY_INCLUDE_DIR}")
        include_directories(${NUMPY_INCLUDE_DIR})
        target_link_libraries(P PRIVATE Python::Python)
        add_definitions(-DMATPLOTLIB_CPP_ENABLED)
    else()
        message(WARNING "NumPy not found - matplotlib-cpp will be disabled")
    endif()
else()
    message(STATUS "Python not found - matplotlib-cpp will be disabled")
endif()

# Enable testing if needed
enable_testing()
