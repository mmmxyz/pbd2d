cmake_minimum_required(VERSION 3.10)

project(blfem CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)


if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)


add_subdirectory(../../opengl ../../../opengl/build)
add_subdirectory(../../mathfunc ../../../mathfunc/build)
add_subdirectory(../../3rdparty/glad ../../../3rdparty/glad/build)

add_executable(Program main.cpp)

target_include_directories(Program PRIVATE
	../..
)

if (APPLE)
	link_directories(Program
		"/usr/local/lib"
		)
endif()

set(EX_GLFWLIB "")
if (APPLE)
	set(EX_GLFWLIB "-framework Cocoa -framework OpenGL -framework IOKit")
endif()

target_link_libraries(Program PRIVATE
	mathLib
	glLib
	${EX_GLFWLIB}
	${OPENGL_LIBRARY}
	glfw
	gladLib
	${CMAKE_DL_LIBS}
)
