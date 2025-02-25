#Minimum requirement of CMake version : 3.5.0
cmake_minimum_required(VERSION 3.5.0)
cmake_policy(SET CMP0074 NEW)

#Project name and version number
project(${App})
set(EXECUTABLE_OUTPUT_PATH ${ROOT}/build/Bin)

file(GLOB source_files *.c)
source_group("Cmake Files" FILES CMakeLists.txt)
source_group("Source Files" FILES source_files)

set(include_path ${PROJECT_SOURCE_DIR}/include ${libGnssSrc} ${libInsSrc} ${libFusionSrc})
include_directories(${include_path})
set(lib_list ${libGnss} ${libIns} ${libFusion})

add_executable(rnx2rtkp rnx2rtkp.c)
#add_executable(inspure inspure.c)

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    link_directories(${ROOT}/build/Lib/Debug)
    link_directories(${ROOT}/build/Lib/Release)
    link_directories(${ROOT}/build/Lib/RelWithDebInfo)
    link_directories(${ROOT}/build/Lib/MinSizeRel)
else ()
    link_directories(${ROOT}/build/Lib)
endif ()

target_link_libraries(rnx2rtkp ${lib_list})
#target_link_libraries(inspure  ${lib_list})
