cmake_minimum_required(VERSION 3.21)
project(wisswesser)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
message(STATUS "CMake System: ${CMAKE_SYSTEM_NAME}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fsanitize=address -Wall -Wpedantic -Wextra")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

find_package(OpenBabel3 REQUIRED)
if(OPENBABEL3_FOUND)
  set(USE_OPENBABEL true)
  set(OPENBABEL_INCLUDE_DIR ${OPENBABEL3_INCLUDE_DIR})
  set(OPENBABEL_LIBRARIES ${OPENBABEL3_LIBRARIES})
  include_directories(${OPENBABEL_INCLUDE_DIR})
else()
  message(FATAL_ERROR "OpenBabel library not found")
endif(OPENBABEL3_FOUND)

include_directories(
  ${CMAKE_SOURCE_DIR}/src/wlnparser
  ${CMAKE_SOURCE_DIR}/src/wlngrep
  ${CMAKE_SOURCE_DIR}/src/wlncompress
)

add_executable(readwln 
  ${CMAKE_SOURCE_DIR}/src/wlnparser/readwln.cpp 
  ${CMAKE_SOURCE_DIR}/src/wlnparser/wlnreader.cpp
)

add_executable(benchwln 
  ${CMAKE_SOURCE_DIR}/src/wlnparser/benchwln.cpp 
  ${CMAKE_SOURCE_DIR}/src/wlnparser/wlnreader.cpp
)

target_link_libraries(readwln   ${OPENBABEL3_LIBRARIES})
target_link_libraries(benchwln  ${OPENBABEL3_LIBRARIES})

add_executable(wlngrep ${CMAKE_SOURCE_DIR}/src/wlngrep/wlngrep.cpp)

# this will not compile on Arm64 // fix at later date
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_executable(wlnpaq6
    ${CMAKE_SOURCE_DIR}/src/wlncompress/paq6.cpp
  )
  set_target_properties(wlnpaq6 PROPERTIES COMPILE_OPTIONS "-m32" LINK_FLAGS "-m32") 
endif() 
