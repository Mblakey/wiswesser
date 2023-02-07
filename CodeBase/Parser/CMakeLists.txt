cmake_minimum_required(VERSION 3.21)
project(wisswesser)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")


include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
)

include_directories("${CMAKE_SOURCE_DIR}/../openbabel/build/include/")
include_directories("${CMAKE_SOURCE_DIR}/../openbabel/include/")


add_executable(wisswesser parser.cpp readwln.cpp writewln.cpp)

add_executable(new-writer writewln2.cpp)

target_link_libraries(wisswesser "${CMAKE_SOURCE_DIR}/../openbabel/build/lib/libopenbabel.7.0.0.dylib")
target_link_libraries(new-writer "${CMAKE_SOURCE_DIR}/../openbabel/build/lib/libopenbabel.7.0.0.dylib")