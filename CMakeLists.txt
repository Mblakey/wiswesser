
project(wisswesser)
cmake_minimum_required(VERSION 3.21)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(
        ${PROJECT_SOURCE_DIR}/src/parser
        ${PROJECT_SOURCE_DIR}/src/wlngrep
        ${PROJECT_SOURCE_DIR}/src/wlncompress
        ${PROJECT_SOURCE_DIR}/src/wlngen
)

add_executable(readwln ${PROJECT_SOURCE_DIR}/src/parser/reader.cpp ${PROJECT_SOURCE_DIR}/src/parser/readwln2.cpp ${PROJECT_SOURCE_DIR}/src/parser/wln-nextmove.cpp)
add_executable(writewln ${PROJECT_SOURCE_DIR}/src/parser/writer.cpp ${PROJECT_SOURCE_DIR}/src/parser/writewln2.cpp ${PROJECT_SOURCE_DIR}/src/parser/modernwln.cpp)
add_executable(obcomp ${PROJECT_SOURCE_DIR}/src/parser/obcomp.cpp)
add_executable(wlngrep ${PROJECT_SOURCE_DIR}/src/wlngrep/wlngrep.cpp)

# old execs
# add_executable(wlnarithmetic ${PROJECT_SOURCE_DIR}/src/wlncompress/wlnarithmetic.cpp)
# add_executable(wlnhuffman ${PROJECT_SOURCE_DIR}/src/wlncompress/wlnhuffman.cpp) # will remove, just for testing

add_executable(wlnzip ${PROJECT_SOURCE_DIR}/src/wlncompress/wlnzip.cpp)
add_executable(wlngen ${PROJECT_SOURCE_DIR}/src/wlngen/wlngen.cpp ${PROJECT_SOURCE_DIR}/src/parser/readwln2.cpp)
add_executable(wlnlearn ${PROJECT_SOURCE_DIR}/src/wlngen/wlnlearn.cpp ${PROJECT_SOURCE_DIR}/src/parser/readwln2.cpp)


include_directories("${CMAKE_SOURCE_DIR}/external/openbabel/build/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/openbabel/include/")

if(WIN32)
  message(ERROR "Unsupported operating system - Windows")
elseif(APPLE)
  target_link_libraries(readwln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(writewln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(obcomp "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(wlngen "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
  target_link_libraries(wlnlearn "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.7.dylib")
elseif(UNIX)
  target_link_libraries(readwln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(writewln "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(obcomp "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(wlngen "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
  target_link_libraries(wlnlearn "${CMAKE_SOURCE_DIR}/external/openbabel/build/lib/libopenbabel.so.7")
endif()


target_compile_definitions(readwln PRIVATE ERRORS=1)
target_compile_definitions(wlngen PRIVATE ERRORS=0)
target_compile_definitions(wlnlearn PRIVATE ERRORS=0)

