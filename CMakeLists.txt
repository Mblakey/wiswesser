cmake_minimum_required(VERSION 3.21)
project(wisswesser)

message(STATUS "CMake System: ${CMAKE_SYSTEM_NAME}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(
        ${PROJECT_SOURCE_DIR}/src/wlnparser
        ${PROJECT_SOURCE_DIR}/src/wlngrep
        ${PROJECT_SOURCE_DIR}/src/wlncompress
)

add_executable(readwln 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/reader.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/readwln.cpp
  ${PROJECT_SOURCE_DIR}/src/wlnparser/readwln2.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/wln-nextmove.cpp
)

add_executable(writewln 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/writer.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlnparser/writewln2.cpp 
)


add_executable(wlngrep ${PROJECT_SOURCE_DIR}/src/wlngrep/wlngrep.cpp)

add_executable(wlnzip 
  ${PROJECT_SOURCE_DIR}/src/wlncompress/wlnzip.cpp
  ${PROJECT_SOURCE_DIR}/src/wlncompress/wlnppm.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlncompress/wlndeflate.cpp 
  ${PROJECT_SOURCE_DIR}/src/wlncompress/lempelz.cpp
  ${PROJECT_SOURCE_DIR}/src/wlncompress/huffman.cpp
  ${PROJECT_SOURCE_DIR}/src/wlncompress/context_trie.cpp
)


# this will not compile on Arm64 // fix at later date
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_executable(wlnpaq6
    ${PROJECT_SOURCE_DIR}/src/wlncompress/paq6.cpp
  )
  set_target_properties(wlnpaq6 PROPERTIES COMPILE_OPTIONS "-m32" LINK_FLAGS "-m32") 
endif() 

include_directories("${CMAKE_SOURCE_DIR}/external/openbabel/build/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/openbabel/include/")

if(WIN32)
  message(ERROR "Unsupported operating system - Windows")
elseif(MACOSX)
  target_link_libraries(readwln  "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(writewln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
elseif(UNIX)
  target_link_libraries(readwln  "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(writewln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
endif()

