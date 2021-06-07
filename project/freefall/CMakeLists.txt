cmake_minimum_required(VERSION 3.10)

project(freefall CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)


if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)


add_subdirectory(../../opengl ../../../opengl/build)
add_subdirectory(../../mathfunc ../../../mathfunc/build)

add_executable(Program main.cpp)

target_include_directories(Program PRIVATE
	../..
)



set(EX_GLFWLIB "")
if (APPLE)
	set(EX_GLFWLIB "-framework Cocoa -framework IOKit -framework CoreFoundation")
endif()

target_link_libraries(Program PRIVATE
	mathLib
	glLib
	${OPENGL_LIBRARY}
	glfw
	${EX_GLFWLIB}
	${GLEW_LIBRARIES}
)
