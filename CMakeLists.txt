cmake_minimum_required(VERSION 3.21)
project(wisswesser)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
message(STATUS "CMake System: ${CMAKE_SYSTEM_NAME}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fsanitize=address -Wall -Wpedantic -Wextra")
endif()

option(UNIT_TESTS "Enable/disable unit testing" TRUE)

if(UNIT_TESTS)
  find_package(Catch2 3)
  if(Catch2_FOUND)
    include(CTest)
    include(Catch)
    message(STATUS "Catch2 found")
  else()
    message(STATUS "Catch2 not found! Testing is disabled")
  endif()
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
  ${CMAKE_SOURCE_DIR}/src/wlnparser/reader.cpp 
  ${CMAKE_SOURCE_DIR}/src/wlnparser/readwln.cpp
  ${CMAKE_SOURCE_DIR}/src/wlnparser/readwln2.cpp 
)

add_executable(writewln 
  ${CMAKE_SOURCE_DIR}/src/wlnparser/writer.cpp 
  ${CMAKE_SOURCE_DIR}/src/wlnparser/writewln2.cpp 
)

add_executable(wlngrep ${CMAKE_SOURCE_DIR}/src/wlngrep/wlngrep.cpp)

# add_executable(wlnzip 
#   ${CMAKE_SOURCE_DIR}/src/wlncompress/wlnzip.cpp
#   ${CMAKE_SOURCE_DIR}/src/wlncompress/wlnppm.cpp 
#   ${CMAKE_SOURCE_DIR}/src/wlncompress/wlndeflate.cpp 
#   ${CMAKE_SOURCE_DIR}/src/wlncompress/lempelz.cpp
#   ${CMAKE_SOURCE_DIR}/src/wlncompress/huffman.cpp
#   ${CMAKE_SOURCE_DIR}/src/wlncompress/context_trie.cpp
# )
#

# this will not compile on Arm64 // fix at later date
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_executable(wlnpaq6
    ${CMAKE_SOURCE_DIR}/src/wlncompress/paq6.cpp
  )
  set_target_properties(wlnpaq6 PROPERTIES COMPILE_OPTIONS "-m32" LINK_FLAGS "-m32") 
endif() 

target_link_libraries(readwln  ${OPENBABEL3_LIBRARIES})
target_link_libraries(writewln ${OPENBABEL3_LIBRARIES})

# if (Catch2_FOUND)
#   add_executable(test_read ${CMAKE_SOURCE_DIR}/test/wln_catch2.cpp
#                             ${CMAKE_SOURCE_DIR}/src/.cpp)
#   target_link_libraries(test_read PRIVATE Catch2::Catch2WithMain)
#   catch_discover_tests(test_nmidx)
#   target_compile_options(test_nmidx PUBLIC "-std=gnu++17") 
#   target_compile_definitions(test_nmidx PUBLIC TEST_DIR="${CMAKE_SOURCE_DIR}/tests" BIN_DIR="${CMAKE_CURRENT_BINARY_DIR}")
#   add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
#                     DEPENDS test_nmidx)
# endif()
#

