#Minimum requirement of CMake version : 3.5.0
cmake_minimum_required(VERSION 3.10.0) #项目执行所需最低版本cmake

#Project name and version number
project(RTKNav VERSION 0.1.0) #设置项目名称

# 根据构建类型设置优化选项
if(CMAKE_BUILD_TYPE MATCHES "Release")
    # Release模式：最高优化级别
    add_compile_options(-O3) 		#最高级别优化
    add_compile_options(-DNDEBUG) 	#禁用调试断言
    add_compile_options(-ffast-math) #快速数学运算
    if(CMAKE_SYSTEM_NAME MATCHES "Linux")
        add_compile_options(-march=native) #针对当前CPU架构优化
    endif()
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    # 带调试信息的发布版本
    add_compile_options(-O2) 		#平衡优化
    add_compile_options(-g) 		#保留调试信息
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    # Debug模式：不优化，便于调试
    add_compile_options(-O0) 		#不进行优化
    add_compile_options(-g) 		#生成调试信息
    add_compile_options(-Wno-unused-variable) #忽略未使用变量警告
elseif(CMAKE_BUILD_TYPE MATCHES "Asan")
    # ASAN模式：地址消毒器检测内存错误
    add_compile_options(-O0) 		#不进行优化
    # add_compile_options(-g) 		#生成调试信息
    add_compile_options(-fsanitize=address) 	#启用地址消毒器
    add_compile_options(-fno-omit-frame-pointer) #保留帧指针
    add_link_options(-fsanitize=address) 	#链接时也启用ASAN
else()
    # 默认模式：中等优化
    add_compile_options(-O2) 		#中等优化级别
endif()

# 平台特定的通用选项
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    add_compile_options(-pthread)   #启用多线程
    add_compile_options(-m64) 		#指定生成64位可执行文件
    add_compile_options(-Wno-unused-variable) #忽略未使用变量警告
    # add_compile_options(-lz) 		#链接zlib库
    # add_compile_options(-lstdc++) 	#链接标准C++库
    # add_compile_options(-lm)        # 添加数学库链接
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    add_compile_options(-pthread)
    add_compile_options(-Wall) 		#启用所有警告信息
elseif(CMAKE_BUILD_TYPE MATCHES "Asan") # 设置ASAN环境变量
    # 设置ASAN选项
    set(ENV{ASAN_OPTIONS} "detect_leaks=1:symbolize=1:abort_on_error=1:check_initialization_order=1:strict_init_order=1")
    # 在Windows上可能需要额外的设置
    if(WIN32)
        set(ENV{ASAN_OPTIONS} "$ENV{ASAN_OPTIONS}:windows_hook_rtl_allocators=1")
    endif()
endif()

if(WIN32)
    MATH(EXPR stack_size "600*1024*1024")   # 使用数学表达式计算变量stack_size的值，将其设置为600MB
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--stack,${stack_size}") # 将程序的栈大小设置为之前计算得到的stack_size值
endif()

include(CheckCXXCompilerFlag)#检查指定的C++编译器标志是否可用，并将结果存储在变量中供后续使用
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

# 检查编译器是否支持C++11和C++0x标准
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    # 链接时优化(LTO) - 如果编译器支持
    CHECK_CXX_COMPILER_FLAG("-flto" COMPILER_SUPPORTS_LTO)
    if(COMPILER_SUPPORTS_LTO AND CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -flto")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#Set the ROOT and subdirectory, you should put the CMakeList.txt in these file directories
set(ROOT ${PROJECT_SOURCE_DIR}) #将项目路径${PROJECT_SOURCE_DIR}存到变量ROOT中

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT}/build/Lib) #设置动态库输出路径
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT}/build/Lib) #设置静态库输出路径
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT}/build/Bin) #设置可执行文件输出路径

# 根据构建类型设置编译配置
if(CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_BUILD_POSTFIX "${CMAKE_RELEASE_POSTFIX}")
    # Release模式额外优化选项
    if(MSVC)
        add_compile_options(/Ox)  # MSVC最大优化
        add_compile_options(/Ob2) #内联函数扩展
        add_compile_options(/Ot)  #速度优先
    else()
        add_compile_options(-O3)  # GCC最大优化
        add_compile_options(-ffast-math) #快速数学运算
        if(NOT WIN32)
            add_compile_options(-march=native) #针对当前CPU架构优化
        endif()
    endif()
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_BUILD_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
    # Debug模式调试选项
    if(MSVC)
        add_compile_options(/Od)  #禁用优化
        add_compile_options(/Zi)  #生成调试信息
    else()
        add_compile_options(-O0)  # GCC禁用优化
        add_compile_options(-g)   #生成调试信息
    endif()
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    set(CMAKE_BUILD_POSTFIX "${CMAKE_RELWITHDEBINFO_POSTFIX}")
    # 带调试信息的发布版本
    if(MSVC)
        add_compile_options(/O2)  #最大速度优化
        add_compile_options(/Zi)  #生成调试信息
    else()
        add_compile_options(-O2)  # GCC平衡优化
        add_compile_options(-g)   #生成调试信息
    endif()
elseif(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
    set(CMAKE_BUILD_POSTFIX "${CMAKE_MINSIZEREL_POSTFIX}")
    # 最小尺寸优化
    if(MSVC)
        add_compile_options(/O1)  #最小尺寸优化
    else()
        add_compile_options(-Os)  # GCC最小尺寸优化
    endif()
else()
    set(CMAKE_BUILD_POSTFIX "")
    # 默认使用Release优化
    if(MSVC)
        add_compile_options(/O2)
    else()
        add_compile_options(-O2)
    endif()
endif()

set(CMAKE_DEBUG_POSTFIX "d") #在Debug配置下生成的目标文件名后添加"d"作为后缀
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_RELWITHDEBINFO_POSTFIX "rd")
set(CMAKE_MINSIZEREL_POSTFIX "s")

include_directories(include) #添加包含路径

# 仅在Windows平台添加WIN32宏定义，其他宏保持不变
# if(WIN32)
#     add_definitions(-DWIN32)
# endif()
# add_definitions(-DENAGLO -DENAGAL -DENACMP -DENAQZS -DENAIRN -DNEXOBS=3 -DNFREQ=7 -DTRACE)#添加预处理器定义，这些定义在编译时会被添加到源代码中的宏定义中
add_compile_definitions(WIN32 ENAGLO ENAGAL ENACMP ENAQZS ENAIRN NEXOBS=3 NFREQ=7 TRACE)#添加预处理器定义，这些定义在编译时会被添加到源代码中的宏定义中

# global value
# set lib name
set(libGnss libGnss)
set(libIns libIns)
set(libFusion libFusion)
# set lib path
set(libGnssSrc ${PROJECT_SOURCE_DIR}/src/LibGnss)
set(libInsSrc  ${PROJECT_SOURCE_DIR}/src/LibIns)
set(libFusionSrc  ${PROJECT_SOURCE_DIR}/src/LibFusion)

# For lib source files
#将 ${libGnssSrc}  目录下的源代码添加到当前构建过程中，并将构建结果输出到 ${ROOT}/build/${libGnss} 目录下
add_subdirectory(${libGnssSrc}   ${ROOT}/build/${libGnss}) 
add_subdirectory(${libInsSrc}    ${ROOT}/build/${libIns}) 
add_subdirectory(${libFusionSrc} ${ROOT}/build/${libFusion}) 

#CUI APP
set(App GINSLIB)
add_subdirectory(${ROOT}/src/app ${ROOT}/build/${App}) 

#Output Messages for debug the Cmake
message(STATUS "operation system is : ${CMAKE_SYSTEM}") #向控制台打印信息
message(STATUS "current platform is : ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMake version is    : ${CMAKE_SYSTEM_VERSION}")
message(STATUS "C compiler is : ${CMAKE_C_COMPILER}")
message(STATUS "The program main directory is : ${ROOT}")
message(STATUS "The program build type is : ${CMAKE_BUILD_TYPE}")