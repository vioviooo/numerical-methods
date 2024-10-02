cmake_minimum_required(VERSION 3.10)
project(task_1_2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

configure_file(test.txt test.txt COPYONLY)

add_executable(
        task_1_2
        main.cpp)