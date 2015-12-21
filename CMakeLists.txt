cmake_minimum_required(VERSION 3.3)

set (SRC ./src)

include_directories(${SRC})
add_executable(main ${SRC}/main.cpp ${SRC}/DBM.cpp)


set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_COMPILER_ID Clang)
add_definitions( -std=c++11 -Wall )
