cmake_minimum_required(VERSION 3.21)
project(wisswesser)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")

if (DEBUGWLN)
  if(DEBUGWLN STREQUAL "on" OR DEBUGWLN STREQUAL "ON")
    message(STATUS "WLN Debugging mode - enabled\n")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUGWLN=1")
  elseif(DEBUGWLN STREQUAL "off" OR DEBUGWLN STREQUAL "OFF")
    message(STATUS "WLN Debugging mode - disabled\n")
  else() 
    message (WARNING "Invalid selection for DEBUGWLN, must be on|ON or off|OFF\n")
  endif()
else()
  message(STATUS "WLN Debugging mode - disabled\n")
endif()


include_directories(
        ${PROJECT_SOURCE_DIR}/src/parser
        ${PROJECT_SOURCE_DIR}/src/wlngrep
)

add_executable(readwln ${PROJECT_SOURCE_DIR}/src/parser/readwln2.cpp)
add_executable(OLD_readwln ${PROJECT_SOURCE_DIR}/src/parser/wln-nextmove.cpp)
add_executable(writewln ${PROJECT_SOURCE_DIR}/src/parser/writewln2.cpp)
add_executable(obcomp ${PROJECT_SOURCE_DIR}/src/parser/obcomp.cpp)

add_executable(wlngrep ${PROJECT_SOURCE_DIR}/src/wlngrep/wlngrep.cpp)

include_directories("${CMAKE_SOURCE_DIR}/external/openbabel/build/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/openbabel/include/")

if(WIN32)
  message(ERROR "Unsupported operating system - windows")
elseif(APPLE)
  target_link_libraries(readwln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(OLD_readwln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(writewln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(obcomp "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
elseif(UNIX)
  target_link_libraries(readwln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(OLD_readwln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(writewln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(obcomp "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
endif()
