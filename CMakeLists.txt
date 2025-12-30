cmake_minimum_required(VERSION 3.21)
project(wisswesser)
enable_testing()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
message(STATUS "CMake System: ${CMAKE_SYSTEM_NAME}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wpedantic -Wextra")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

find_package(OpenBabel3 REQUIRED)
if(OPENBABEL3_FOUND)
  set(OPENBABEL_INCLUDE_DIR ${OPENBABEL3_INCLUDE_DIR})
  set(OPENBABEL_LIBRARIES ${OPENBABEL3_LIBRARIES})
  include_directories(${OPENBABEL_INCLUDE_DIR})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_OPENBABEL")
else()
  message(FATAL_ERROR "OpenBabel library not found")
endif(OPENBABEL3_FOUND)

include_directories(${CMAKE_SOURCE_DIR}/src/)

add_executable(readwln 
  ${CMAKE_SOURCE_DIR}/src/readwln.cpp 
  ${CMAKE_SOURCE_DIR}/src/wlnreader.cpp
)

add_executable(writewln 
  ${CMAKE_SOURCE_DIR}/src/writewln.cpp 
  ${CMAKE_SOURCE_DIR}/src/wlnwriter.cpp
)

add_executable(benchwln 
  ${CMAKE_SOURCE_DIR}/src/benchwln.cpp 
  ${CMAKE_SOURCE_DIR}/src/wlnreader.cpp
)

add_executable(wlngrep ${CMAKE_SOURCE_DIR}/src/wlngrep.cpp)

target_link_libraries(readwln   ${OPENBABEL3_LIBRARIES})
target_link_libraries(benchwln  ${OPENBABEL3_LIBRARIES})
target_link_libraries(writewln  ${OPENBABEL3_LIBRARIES})

