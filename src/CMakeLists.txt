cmake_minimum_required(VERSION 3.20)
project(MATRIX C)

add_library(matrix matrix.c)
add_executable(main main.c)

target_link_libraries(main PRIVATE matrix)