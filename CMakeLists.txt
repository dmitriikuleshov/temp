cmake_minimum_required(VERSION 3.2)
project(lw6)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES false ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_executable(main main.cpp)
add_library(${PROJECT_NAME} ${SOURCES})
add_executable(tests tests/tests.cpp)

target_link_libraries(main PRIVATE ${PROJECT_NAME})
target_include_directories(tests PRIVATE ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
target_link_libraries(tests PRIVATE gtest gtest_main ${PROJECT_NAME})