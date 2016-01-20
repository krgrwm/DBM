cmake_minimum_required(VERSION 2.8)

set (SRC ./src)

include_directories(${SRC})
add_executable(main ${SRC}/main.cpp ${SRC}/DBM.cpp ${SRC}/plist.cpp)


set(CMAKE_BUILD_TYPE Release)
add_definitions( -std=c++11 -stdlib=libc++ -Wall )
