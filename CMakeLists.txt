cmake_minimum_required(VERSION 3.10)

project(NanoSpring 
    VERSION 0.1
    LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
include(CheckIPOSupported)
include(CheckCXXSourceCompiles)

# Clear All Flags
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
set(CMAKE_CXX_FLAGS_RELEASE "")


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

check_ipo_supported(RESULT result)
if(result)
    message(STATUS "Support for Interprocedural Optimization detected.")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CHECKS=${CHECKS} (Use -DCHECKS=True to enable all checks")
message(STATUS "LOG=${LOG} (Use -LOG=True to output log files for analysis")
message(STATUS "Directory: ${CMAKE_CURRENT_SOURCE_DIR}")
find_package(OpenMP)

##############################
# malloc_trim
check_cxx_source_compiles(
    "#include <malloc.h>
     int main() {
         malloc_trim(0);
         return 0;
     }"
     USE_MALLOC_TRIM)
if(USE_MALLOC_TRIM)
    message(STATUS "Support for malloc_trim detected.")
else()
    message(STATUS "Support for malloc_trim not detected.")
endif()


###############################
# boost
set( Boost_NO_SYSTEM_PATHS ON ) # do not use system boost
add_subdirectory(boost-cmake)

#######################################################
# libbsc
add_library(libbsc
    STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/adler32/adler32.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/bwt/libsais/libsais.c
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/bwt/bwt.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/coder/coder.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/coder/qlfc/qlfc.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/coder/qlfc/qlfc_model.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/filters/detectors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/filters/preprocessing.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/libbsc/libbsc.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/lzp/lzp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/platform/platform.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/libbsc/st/st.cpp
)

target_include_directories(libbsc 
    PUBLIC libbsc)

target_compile_options(libbsc 
    PUBLIC
    -g -Wall

    # Comment out CFLAGS line below for compatability mode for 32bit file sizes
    # (less than 2GB) and systems that have compilers that treat int as 64bit
    # natively (ie: modern AIX)
    -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

    # Comment out CFLAGS line below to disable optimizations
    -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math

    # Comment out CFLAGS line below to disable AVX2 instruction set (performance will suffer)
    -mavx2

    # Comment out CFLAGS line below to disable OpenMP optimizations
    # We disable this here and parallelize at the top level
    # -fopenmp -DLIBBSC_OPENMP_SUPPORT

    # Comment out CFLAGS line below to enable debug output
    -DNDEBUG

    # Comment out CFLAGS line below to disable unaligned memory access
    -DLIBBSC_ALLOW_UNALIGNED_ACCESS
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
#fast-lzma2

add_custom_target(
   FASTLZMA2
   COMMAND make clean && make
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/fast-lzma2/
)

#######################################################
# the NanoSpringLib library
add_library(NanoSpringLib 
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
        src/lzma2.cpp
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
        include/lzma2_helper.h
        )

target_include_directories(NanoSpringLib 
    PUBLIC include src
    PUBLIC libbsc
    PUBLIC minimap2 fast-lzma2)

target_compile_options(NanoSpringLib 
    PUBLIC 
    -W -Wall -Wextra)

target_compile_options(NanoSpringLib 
    PUBLIC 
    $<$<CONFIG:Release>:-O3 -msse4.1>)
target_compile_options(NanoSpringLib 
    PUBLIC 
    $<$<CONFIG:RelWithDebInfo>:-g -O3 -msse4.1>)
target_compile_options(NanoSpringLib 
    PUBLIC 
    $<$<CONFIG:Debug>:-g -O3 -msse4.1 -fsanitize=address -fno-omit-frame-pointer>)

if ( CHECKS )
    target_compile_options(NanoSpringLib 
        PUBLIC 
        -DCHECKS)
endif()
if ( LOG )
    target_compile_options(NanoSpringLib 
        PUBLIC 
        -DLOG)
endif()
if ( USE_MALLOC_TRIM )
    target_compile_options(NanoSpringLib 
        PUBLIC 
        -DUSE_MALLOC_TRIM)
endif()

target_compile_definitions(NanoSpringLib PUBLIC "$<$<CONFIG:Debug>:DEBUG>")
target_compile_definitions(NanoSpringLib PUBLIC "$<$<CONFIG:RelWithDebInfo>:DEBUG>")

add_dependencies(NanoSpringLib MINIMAP2)
add_dependencies(NanoSpringLib FASTLZMA2)

target_link_libraries(NanoSpringLib 
    PUBLIC 
        libbsc
        ${CMAKE_CURRENT_SOURCE_DIR}/minimap2/libminimap2.a
        ${CMAKE_CURRENT_SOURCE_DIR}/fast-lzma2/libfast-lzma2.a
        Boost::filesystem
        Boost::iostreams
        Boost::program_options
        -lz -ldl -lm -lpthread
    )

target_link_libraries(NanoSpringLib 
    PUBLIC 
    $<$<CONFIG:Debug>:-fsanitize=address -fno-omit-frame-pointer>)

# OpenMp support
if(OpenMP_CXX_FOUND)
    target_link_libraries(NanoSpringLib PUBLIC OpenMP::OpenMP_CXX)
endif()

###############################################

add_executable(NanoSpring src/main.cpp)
target_link_libraries(NanoSpring PUBLIC NanoSpringLib)

add_executable(testLoneReads src/testLoneReads.cpp)
target_include_directories(testLoneReads PUBLIC minimap2)
target_link_libraries(testLoneReads PUBLIC NanoSpringLib)

