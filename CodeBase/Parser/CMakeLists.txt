cmake_minimum_required(VERSION 3.21)
project(wisswesser)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
include(GNUInstallDirs)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/api
)

file(GLOB PACKAGE_SRCS
        "${PROJECT_SOURCE_DIR}/include/*.h"
        "${PROJECT_SOURCE_DIR}/include/*.hpp"
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.c"
        "${PROJECT_SOURCE_DIR}/api/*.cpp"
        "${PROJECT_SOURCE_DIR}/api/*.c"
        "${PROJECT_SOURCE_DIR}/api/*.h"
        "${PROJECT_SOURCE_DIR}/api/*.hpp"
        )


file(GLOB EXEC_SRCS
        "${PROJECT_SOURCE_DIR}/include/*.h"
        "${PROJECT_SOURCE_DIR}/include/*.hpp"
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.c"
        )

set(PACKAGE_SOURCES ${PACKAGE_SRCS})
set(EXEC_SOURCES ${EXEC_SRCS})


find_package(OpenBabel3 REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
message("pybind11_VERSION:${PYBIND11_VERSION}")

add_executable(TEST ${EXEC_SOURCES})
pybind11_add_module(${PROJECT_NAME} ${PACKAGE_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC
        ${OpenBabel3_INCLUDE_DIRS}
        ${pybind11_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        PRIVATE src)

target_include_directories(TEST PUBLIC
        ${OpenBabel3_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        PRIVATE src)

target_link_libraries(${PROJECT_NAME} PUBLIC
        ${pybind11_INCLUDE_DIRS}
        ${OpenBabel3_LIBRARIES}
     )


target_link_libraries(TEST PUBLIC
        ${OpenBabel3_LIBRARIES}
        )

target_compile_definitions(${PROJECT_NAME}
        PUBLIC VERSION_INFO=${EXAMPLE_VERSION_INFO})

target_compile_options(TEST PRIVATE -Wall -Wextra -Wpedantic)

install(TARGETS TEST DESTINATION bin)