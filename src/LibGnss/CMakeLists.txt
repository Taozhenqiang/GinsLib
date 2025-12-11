#Minimum requirement of CMake version : 3.5.0
cmake_minimum_required(VERSION 3.10.0)	#项目执行所需最低版本cmake

#Project name and version number
project(${libGnss} VERSION 0.1.0)	#设置项目名称

# set include path
include_directories(${libInsSrc})
include_directories(${libFusionSrc})

# set include path
if (WIN32)
    link_libraries(-lWs2_32 -lwinmm)
ENDIF ()

# for all the header files and source_files
file(GLOB_RECURSE header_files ${libGnssSrc}/*.h)
file(GLOB_RECURSE source_files ${libGnssSrc}/*.c)

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    file(GLOB rcv_header_files rcv/*.h)
    file(GLOB rcv_source_files rcv/*.c)

    # add files to project
    source_group("CMake Files" FILES CMakeLists.txt)
    source_group("Header Files" FILES header_files)
    source_group("Soruce Files" FILES source_files)
    #source_group("Header Files\\rcv" FILES ${rcv_header_files})
    #source_group("Source Files\\rcv" FILES ${rcv_source_files})

else ()

endif ()

# 使用指定的源文件创建库
add_library(${PROJECT_NAME} ${header_files} ${source_files})

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    link_directories(${ROOT}/build/Lib/Debug)
    link_directories(${ROOT}/build/Lib/Release)
    link_directories(${ROOT}/build/Lib/RelWithDebInfo)
    link_directories(${ROOT}/build/Lib/MinSizeRel)
else ()
    link_directories(${ROOT}/build/Lib)
endif ()

# 将当前库/可执行文件与其他库进行链接，这意味着当前库/可执行文件可以调用其他库中的函数和变量
target_link_libraries(${PROJECT_NAME} ${libIns} ${libFusion})
