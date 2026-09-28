cmake_minimum_required(VERSION 3.15)

project(gvc_ds VERSION 1.0.0 LANGUAGES C)

add_library(ds SHARED data_structure.c)
target_include_directories(ds PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(ds PRIVATE c_std_99)

if(MSVC)
    target_compile_options(ds PRIVATE /W4)
else()
    target_compile_options(ds PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(UNIX)
    target_link_libraries(ds PRIVATE m)
endif()
