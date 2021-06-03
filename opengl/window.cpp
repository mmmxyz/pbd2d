#include <GL/glew.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <cstdint>
#include <iostream>

#include "window.hpp"

Window::Window(uint32_t width, uint32_t height) : window(glfwCreateWindow(width, height, "windoooooow", NULL, NULL))
{
		if (window == NULL)
		{
				std::cerr << "fail to create window" << std::endl;
				exit(1);
		}

		glfwMakeContextCurrent(window);
}

bool Window::is_shouldclose()
{
		return glfwWindowShouldClose(window);
}

void Window::swapbuf()
{
		glfwSwapBuffers(window);
}

void Window::pollevent()
{
		glfwPollEvents();
}

void Window::waitevent(double sec)
{
		glfwWaitEventsTimeout(sec);
}

void Window::clear()
{
		// ウィンドウを消去する
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//塗りつぶし
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}
