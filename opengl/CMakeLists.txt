cmake_minimum_required(VERSION 3.10)

project(glLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)


set(EX_INC "")
if(APPLE)
	set(EX_INC "/usr/local/include")
endif()

include_directories(
	../3rdparty/glad/include
	..
	${OPENGL_INCLUDE_DIRS}
	${glfw3_INCLUDE_DIRS}
	${EX_INC}
)

add_library(glLib STATIC
	line.cpp
	point.cpp
	vertarray.cpp
	visualize.cpp
	window.cpp
)
