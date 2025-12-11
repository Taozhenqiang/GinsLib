#Minimum requirement of CMake version : 3.5.0
cmake_minimum_required(VERSION 3.10.0)	#项目执行所需最低版本cmake

#Project name and version number
project(${libFusion} VERSION 0.1.0)	#设置项目名称

# set include path
include_directories(${libGnssSrc})
include_directories(${libInsSrc})

# for all the header files and source_files
file(GLOB_RECURSE header_files ${libFusionSrc}/*.h)
file(GLOB_RECURSE source_files ${libFusionSrc}/*.c)

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    # add files to project
    source_group("CMake Files"  FILES CMakeLists.txt)
    source_group("Header Files" FILES header_files)
    source_group("Soruce Files" FILES source_files)
else ()

endif ()

add_library(${PROJECT_NAME} ${header_files} ${source_files})

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    link_directories(${ROOT}/build/Lib/Debug)
    link_directories(${ROOT}/build/Lib/Release)
    link_directories(${ROOT}/build/Lib/RelWithDebInfo)
    link_directories(${ROOT}/build/Lib/MinSizeRel)
else ()
    link_directories(${ROOT}/build/Lib)
endif ()

target_link_libraries(${PROJECT_NAME} ${libGnss} ${libIns})