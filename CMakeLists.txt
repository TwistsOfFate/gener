cmake_minimum_required(VERSION 3.19)
project(gener)

set(CMAKE_CXX_STANDARD 14)

add_executable(gener
        main.cpp check_sv.cpp defines.h prepare.cpp declarations.h alignment.cpp commons.cpp concat.cpp)
