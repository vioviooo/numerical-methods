cmake_minimum_required(VERSION 3.10)
project(task_1_4)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

configure_file(matrix.txt matrix.txt COPYONLY)

add_executable(
        task_1_4
        main.cpp)