cmake_minimum_required(VERSION 3.5)

project(generate LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(generate main.cpp)

install(TARGETS generate
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
