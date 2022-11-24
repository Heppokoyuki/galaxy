cmake_minimum_required(VERSION 3.5)

project(cartwheel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(cartwheel main.cpp)

install(TARGETS cartwheel
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
