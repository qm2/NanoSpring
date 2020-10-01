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
# the NanoporeCompressionLib library
add_library(NanoporeCompressionLib 
    STATIC 
        src/AlignerTester.cpp
        src/bsc.cpp
        #src/Compressor.cpp
        #src/Consensus.cpp
        #src/Contig.cpp
        #src/Decompressor.cpp
        src/Edits.cpp
        src/LocalMyers_impl.h
        src/LocalMyersRollBack_impl.h
        # src/LocalMyersRollBackOld.cpp
        #src/ReadAligner.cpp
        src/ReadData.cpp
        src/ReadFilter.cpp
        #src/ConsensusGraph.cpp
        include/AlignerTester.h
        include/bsc_helper.h
        #include/Compressor.h
        #include/Consensus.h
        #include/Contig.h
        #include/Decompressor.h
        include/Edits.h
        include/StringAligner.h
        include/LocalMyers.h
        include/LocalMyersRollBack.h
        # include/LocalMyersRollBackOld.h
        include/ReadAligner.h
        include/ReadData.h
        include/ReadFilter.h
        #include/testAligner.h
        #include/ConsensusGraph.h
        )

target_include_directories(NanoporeCompressionLib 
    PUBLIC include src
    PRIVATE libbsc)

target_compile_options(NanoporeCompressionLib 
    PRIVATE 
    -W -Wall -Wextra)

target_compile_options(NanoporeCompressionLib 
    PRIVATE 
    $<$<CONFIG:Release>:-O3>)
target_compile_options(NanoporeCompressionLib 
    PRIVATE 
    $<$<CONFIG:RelWithDebInfo>:-g -O3>)
target_compile_options(NanoporeCompressionLib 
    PRIVATE 
    $<$<CONFIG:Debug>:-g -O3 -fsanitize=address -fno-omit-frame-pointer>)

target_compile_definitions(NanoporeCompressionLib PUBLIC "$<$<CONFIG:Debug>:DEBUG>")
target_compile_definitions(NanoporeCompressionLib PUBLIC "$<$<CONFIG:RelWithDebInfo>:DEBUG>")


target_link_libraries(NanoporeCompressionLib 
    PRIVATE 
        libbsc
        Boost::filesystem
    )

target_link_libraries(NanoporeCompressionLib 
    PRIVATE 
    $<$<CONFIG:Debug>:-fsanitize=address -fno-omit-frame-pointer>)

# OpenMp support
if(OpenMP_CXX_FOUND)
    target_link_libraries(NanoporeCompressionLib PUBLIC OpenMP::OpenMP_CXX)
endif()

###############################################

add_executable(testAligner src/testAligner.cpp)
target_link_libraries(testAligner PRIVATE NanoporeCompressionLib profiler)

add_executable(testAlignAlgs src/testAlignAlgs.cpp)
target_link_libraries(testAlignAlgs PRIVATE NanoporeCompressionLib profiler)

#add_executable(testCompressor src/testCompressor.cpp)
#target_link_libraries(testCompressor PRIVATE NanoporeCompressionLib profiler)

#add_executable(testContig src/testContigGenerator.cpp)
#target_link_libraries(testContig PRIVATE NanoporeCompressionLib profiler)

#add_executable(testDecompressor src/testDecompressor.cpp)
#target_link_libraries(testDecompressor PRIVATE NanoporeCompressionLib profiler)

add_executable(testMinHash src/testMinHash.cpp)
target_link_libraries(testMinHash PRIVATE NanoporeCompressionLib profiler)

add_executable(testMyers src/testMyers.cpp)
target_link_libraries(testMyers PRIVATE NanoporeCompressionLib profiler)

#add_executable(validateResult src/validateResult.cpp)
#target_link_libraries(validateResult PRIVATE NanoporeCompressionLib profiler)