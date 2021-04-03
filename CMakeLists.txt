cmake_minimum_required(VERSION 3.10)

project(NanoporeCompression 
    VERSION 0.1
    LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Clear All Flags
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
set(CMAKE_CXX_FLAGS_RELEASE "")


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CHECKS=${CHECKS} (Use -DCHECKS=True to enable all checks")
message(STATUS "LOG=${LOG} (Use -LOG=True to output log files for analysis")
message(STATUS "Directory: ${CMAKE_CURRENT_SOURCE_DIR}")
find_package(OpenMP)

###############################
# boost
add_subdirectory(boost-cmake)

#######################################################
# libbsc
add_library(libbsc
    STATIC
        libbsc/libbsc/adler32/adler32.cpp
        libbsc/libbsc/bwt/divsufsort/divsufsort.c
        libbsc/libbsc/bwt/bwt.cpp
        libbsc/libbsc/coder/coder.cpp
        libbsc/libbsc/coder/qlfc/qlfc.cpp
        libbsc/libbsc/coder/qlfc/qlfc_model.cpp
        libbsc/libbsc/filters/detectors.cpp
        libbsc/libbsc/filters/preprocessing.cpp
        libbsc/libbsc/libbsc/libbsc.cpp
        libbsc/libbsc/lzp/lzp.cpp
        libbsc/libbsc/platform/platform.cpp
        libbsc/libbsc/st/st.cpp
)

target_include_directories(libbsc 
    PRIVATE libbsc)

target_compile_options(libbsc 
    PRIVATE 
    -g -Wall
    # Comment out CFLAGS line below for compatability mode for 32bit file sizes
    # (less than 2GB) and systems that have compilers that treat int as 64bit
    # natively (ie: modern AIX)
    -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
    # Comment out CFLAGS line below to disable optimizations
    -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math
    # Comment out CFLAGS line below to disable OpenMP optimizations
    -fopenmp -DLIBBSC_OPENMP_SUPPORT
    # Comment out CFLAGS line below to enable debug output
    -DNDEBUG
    )

if(OpenMP_CXX_FOUND)
    target_link_libraries(libbsc PUBLIC OpenMP::OpenMP_CXX -DLIBBSC_OPENMP_SUPPORT)
endif()

#######################################################
#minimap2

add_custom_target(
   MINIMAP2
   COMMAND make clean && make
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/minimap2/
)


#######################################################
# the NanoporeCompressionLib library
add_library(NanoporeCompressionLib 
    STATIC 
        src/AlignerTester.cpp
        src/bsc.cpp
        src/Compressor.cpp
        src/Consensus.cpp
        src/Decompressor.cpp
        src/Edits.cpp
        src/ReadAligner.cpp
        src/ReadData.cpp
        src/ReadFilter.cpp
        src/ConsensusGraph.cpp
        src/OmpMutex.cpp
        src/DirectoryUtils.cpp
        src/dnaToBits.cpp
        src/BBHashMap.cpp
        include/AlignerTester.h
        include/bsc_helper.h
        include/Compressor.h
        include/Consensus.h
        include/Decompressor.h
        include/Edits.h
        include/StringAligner.h
        include/ReadAligner.h
        include/ReadData.h
        include/ReadFilter.h
        include/testAligner.h
        include/ConsensusGraph.h
        include/OmpMutex.h
        include/DirectoryUtils.h
        include/Types.h
        include/dnaToBits.h
        include/BBHashMap.h
        include/BooPHF.h
        )

target_include_directories(NanoporeCompressionLib 
    PUBLIC include src
    PRIVATE libbsc
    PRIVATE minimap2)

target_compile_options(NanoporeCompressionLib 
    PRIVATE 
    -W -Wall -Wextra)

target_compile_options(NanoporeCompressionLib 
    PRIVATE 
    $<$<CONFIG:Release>:-O3 -DNDEBUG>)
target_compile_options(NanoporeCompressionLib 
    PRIVATE 
    $<$<CONFIG:RelWithDebInfo>:-g -O3>)
target_compile_options(NanoporeCompressionLib 
    PRIVATE 
    $<$<CONFIG:Debug>:-g -O3 -fsanitize=address -fno-omit-frame-pointer>)

if ( CHECKS )
    target_compile_options(NanoporeCompressionLib 
        PRIVATE 
        -DCHECKS)
endif()
if ( LOG )
    target_compile_options(NanoporeCompressionLib 
        PRIVATE 
        -DLOG)
endif()

target_compile_definitions(NanoporeCompressionLib PUBLIC "$<$<CONFIG:Debug>:DEBUG>")
target_compile_definitions(NanoporeCompressionLib PUBLIC "$<$<CONFIG:RelWithDebInfo>:DEBUG>")

add_dependencies(NanoporeCompressionLib MINIMAP2)

target_link_libraries(NanoporeCompressionLib 
    PRIVATE 
        libbsc
        ${CMAKE_CURRENT_SOURCE_DIR}/minimap2/libminimap2.a
        Boost::filesystem
        -lz -ldl -lm -lpthread
    )

target_link_libraries(NanoporeCompressionLib 
    PRIVATE 
    $<$<CONFIG:Debug>:-fsanitize=address -fno-omit-frame-pointer>)

# OpenMp support
if(OpenMP_CXX_FOUND)
    target_link_libraries(NanoporeCompressionLib PUBLIC OpenMP::OpenMP_CXX)
endif()

###############################################

add_executable(testCompressor src/testCompressor.cpp)
target_link_libraries(testCompressor PRIVATE NanoporeCompressionLib)

add_executable(testDecompressor src/testDecompressor.cpp)
target_link_libraries(testDecompressor PRIVATE NanoporeCompressionLib)

add_executable(testLoneReads src/testLoneReads.cpp)
target_include_directories(testLoneReads PRIVATE minimap2)
target_link_libraries(testLoneReads PRIVATE NanoporeCompressionLib)

