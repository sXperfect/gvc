cmake_minimum_required(VERSION 3.9)

project(gvc)

SET(CMAKE_C_FLAGS "-O3 -g -W -Wall -pedantic")

add_library(ds SHARED
    data_structure.c
)
