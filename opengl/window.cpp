#include <glad/glad.h>
#include <KHR/khrplatform.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <cstdint>
#include <iostream>

#include "window.hpp"

static void keyinput(GLFWwindow *window, int key, int scancode, int action, int mods);

Window::Window(uint32_t width, uint32_t height) : window(glfwCreateWindow(width, height, "windoooooow", NULL, NULL))
{
		if (window == NULL)
		{
				std::cerr << "fail to create window" << std::endl;
				exit(1);
		}

		glfwMakeContextCurrent(window);

		glfwSetWindowUserPointer(window, this);

		glfwSetKeyCallback(window, keyinput);
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

void Window::clearstatus()
{
		wait = false;
}

void Window::clear()
{
		// ウィンドウを消去する
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//塗りつぶし
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

bool Window::shouldwait()
{
		return wait;
}

static void keyinput(GLFWwindow *window, int key, int scancode, int action, int mods)
{
		Window *const instance(static_cast<Window *>(glfwGetWindowUserPointer(window)));

		if (GLFW_KEY_SPACE == 32 && (action == 1 || action == 2))
		{
				instance->wait = true;
		}
}
