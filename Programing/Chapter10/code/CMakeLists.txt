cmake_minimum_required(VERSION 3.10)
project(NumericalAnalysisChapter10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# GDB debugging support
option(ENABLE_GDB "Enable GDB debugging support" OFF)
if(ENABLE_GDB)
    set(CMAKE_BUILD_TYPE Debug)
    add_compile_options(-ggdb3 -O0 -fno-omit-frame-pointer)
    message(STATUS "GDB debugging enabled with full debug symbols")
endif()
set(CMAKE_BUILD_TYPE Debug)

# Find Eigen (header-only library)
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(
    ${EIGEN3_INCLUDE_DIR}
    include
    external
)

# Add library (if needed)
# add_library(NumericalAnalysisLib src/...)

# Add main executable
add_executable(main test/main.cpp)
# target_link_libraries(Chapter10Test PRIVATE NumericalAnalysisLib)

# Python/matplotlib-cpp support for plotting
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
        target_link_libraries(main PRIVATE Python::Python)
        add_definitions(-DMATPLOTLIB_CPP_ENABLED)
    else()
        message(WARNING "NumPy not found - matplotlib-cpp will be disabled")
    endif()
else()
    message(STATUS "Python not found - matplotlib-cpp will be disabled")
endif()

# Enable testing if needed
enable_testing()
