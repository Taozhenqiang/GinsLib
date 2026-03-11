#Minimum requirement of CMake version : 3.5.0
cmake_minimum_required(VERSION 3.10.0)

#Project name and version number
project(${App})

file(GLOB source_files *.c)
source_group("Cmake Files" FILES CMakeLists.txt)
source_group("Source Files" FILES source_files)

# set(include_path ${PROJECT_SOURCE_DIR}/include ${libGnssSrc} ${libInsSrc} ${libFusionSrc})
# include_directories(${include_path})
set(lib_list ${libGnss} ${libIns} ${libFusion})

add_executable(GINSLIB rnx2rtkp.c)
add_executable(RTKRCV  rtkrcv.c)
add_executable(POS2KML pos2kml.c)
# add_executable(inspure inspure.c)

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    link_directories(${ROOT}/build/Lib/Debug)
    link_directories(${ROOT}/build/Lib/Release)
    link_directories(${ROOT}/build/Lib/RelWithDebInfo)
    link_directories(${ROOT}/build/Lib/MinSizeRel)
else ()
    link_directories(${ROOT}/build/Lib)
endif ()

# 链接库（添加数学库和系统库）
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    # Linux 需显式链接数学库(-lm)、zlib(-lz)和C++标准库(-lstdc++)
    target_link_libraries(GINSLIB ${lib_list} m z stdc++)
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    # macOS 可能需要链接系统框架
    target_link_libraries(GINSLIB ${lib_list} m)
elseif(WIN32)
    # Windows 数学库已包含在默认链接中，无需额外添加
    target_link_libraries(GINSLIB ${lib_list})
    # Windows 需要链接系统库，如 pthread（多线程）、ws2_32（网络通信）、winmm（定时器）、kernel32、advapi32（高级API服务）
    target_link_libraries(RTKRCV ${lib_list} pthread ws2_32 winmm kernel32 advapi32)
    target_link_libraries(POS2KML ${lib_list})
endif()
#target_link_libraries(inspure  ${lib_list})
