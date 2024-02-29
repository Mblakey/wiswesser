
project(wisswesser)
cmake_minimum_required(VERSION 3.21)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(
        ${PROJECT_SOURCE_DIR}/src/wlnparser
        ${PROJECT_SOURCE_DIR}/src/wlngrep
        ${PROJECT_SOURCE_DIR}/src/wlncompress
        ${PROJECT_SOURCE_DIR}/src/wlngen
)

add_executable(readwln 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/reader.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/readwln2.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/wln-nextmove.cpp
)

add_executable(writewln 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/writer.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/writewln2.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/modernwln.cpp
)

add_executable(obcomp ${PROJECT_SOURCE_DIR}/src/wlnparser/obcomp.cpp)
add_executable(wlngrep ${PROJECT_SOURCE_DIR}/src/wlngrep/wlngrep.cpp)

add_executable(wlnzip 
  ${PROJECT_SOURCE_DIR}/src/wlncompress/compress.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlncompress/wlnzip.cpp
  ${PROJECT_SOURCE_DIR}/src/wlncompress/lz.cpp
  ${PROJECT_SOURCE_DIR}/src/wlncompress/huffman.cpp
)

add_executable(wlnseed 
  ${PROJECT_SOURCE_DIR}/src/wlncluster/wlnseed.cpp 
)

add_executable(wlncluster 
  ${PROJECT_SOURCE_DIR}/src/wlncluster/wlncluster.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlncompress/wlnzip.cpp
  ${PROJECT_SOURCE_DIR}/src/wlncompress/lz.cpp
  ${PROJECT_SOURCE_DIR}/src/wlncompress/huffman.cpp
)

add_executable(wlnfp 
  ${PROJECT_SOURCE_DIR}/src/wlncluster/wlnFP.cpp
)

add_executable(wlngen ${PROJECT_SOURCE_DIR}/src/wlngen/wlnlearn.cpp ${PROJECT_SOURCE_DIR}/src/wlnparser/readwln2.cpp)

include_directories("${CMAKE_SOURCE_DIR}/external/openbabel/build/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/openbabel/include/")

if(WIN32)
  message(ERROR "Unsupported operating system - Windows")
elseif(APPLE)
  target_link_libraries(readwln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(writewln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(obcomp "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(wlngen "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
elseif(UNIX)
  target_link_libraries(readwln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(writewln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(obcomp "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(wlngen "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
endif()


target_compile_definitions(readwln PRIVATE ERRORS=1)
target_compile_definitions(wlngen PRIVATE ERRORS=0)

