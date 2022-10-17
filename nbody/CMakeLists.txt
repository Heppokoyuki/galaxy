cmake_minimum_required(VERSION 3.5)

project(nbody LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(nbody main.cpp)

install(TARGETS nbody
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -save-temps")

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
