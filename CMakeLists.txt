cmake_minimum_required(VERSION 3.10)
project(MyLTCApp)

set(CMAKE_CXX_STANDARD 17)


add_definitions(-Dcimg_display=0)  # 如果你用 X11 直接去掉=0
# add_subdirectory(
    # external/tinyexr
# )

# 设定包含路径
include_directories(
    ${CMAKE_SOURCE_DIR}/fit
    ${CMAKE_SOURCE_DIR}/external/glm
    ${CMAKE_SOURCE_DIR}/external/CImg
    # ${CMAKE_SOURCE_DIR}/external/tinyexr
)

# 定义宏，CImg 要求启用 X11 时需要定义这个（或你可以改用其他输出方式）


# 自动添加 fit 目录下所有 .cpp 和 .h 文件
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/fit/*.cpp  ${CMAKE_SOURCE_DIR}/fit/*.c)
file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/fit/*.h)

# 创建可执行程序
add_executable(MyLTCExecutable ${SOURCES} ${HEADERS})
