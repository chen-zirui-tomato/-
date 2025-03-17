cmake_minimum_required(VERSION 3.10)
project(Chapter7)

# 设置默认构建类型为 Debug
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置工程目录变量
set(PROJECT_DIR ${CMAKE_SOURCE_DIR})
set(PROJECT_INCLUDE_DIR ${PROJECT_DIR}/include)
set(PROJECT_EXTERNAL_INCLUDE_DIR ${PROJECT_DIR}/external)
set(PROJECT_SOURCE_DIR ${PROJECT_DIR}/src)
set(PROJECT_TEST_DIR ${PROJECT_DIR}/test)

# 查找 OpenCV 库
find_package(OpenCV REQUIRED)

# 查找 eigen 库
find_package(Eigen3 REQUIRED)

# # 设置 lunasvg 的路径
# set(CMAKE_PREFIX_PATH ${PROJECT_EXTERNAL_INCLUDE_DIR}/lunasvg/cmake)

# # 添加 external 目录
# add_subdirectory(${PROJECT_EXTERNAL_INCLUDE_DIR}/lunasvg)

# # 查找 lunasvg 库
# find_package(lunasvg REQUIRED)

# 添加头文件搜索路
include_directories(
    ${PROJECT_EXTERNAL_INCLUDE_DIR} 
    ${PROJECT_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIRS}
)

# 得到所有的实现算法的 cpp 源文件
file(GLOB PROJECT_SOURCES 
    "${PROJECT_SOURCE_DIR}/*.cpp" 
    "${PROJECT_SOURCE_DIR}/*/*.cpp"
    )

# 将所有的源文件编译好，避免重复编译
add_library(CommonLib OBJECT ${PROJECT_SOURCES})

# 得到所有的用于测试的 cpp 文件
file(GLOB PROJECT_EXECUTABLES ${PROJECT_TEST_DIR}/*.cpp)

# 对每一个测试的 cpp 文件进行编译并链接生成可执行文件
foreach(PROJECT_EXECUTABLE ${PROJECT_EXECUTABLES})

    # 获取可执行文件的名称
    get_filename_component(PROJECT_EXECUTABLE_NAME ${PROJECT_EXECUTABLE} NAME_WE)

    # 添加可执行文件
    add_executable(${PROJECT_EXECUTABLE_NAME} ${PROJECT_EXECUTABLE})

    # 链接提前编译好的库文件
    target_link_libraries(${PROJECT_EXECUTABLE_NAME} CommonLib)

    # # 链接 lunasvg 库
    # target_link_libraries(${PROJECT_EXECUTABLE_NAME} PRIVATE ${OpenCV_LIBS} lunasvg)
    
endforeach(PROJECT_EXECUTABLE)
