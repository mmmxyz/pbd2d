#pragma once

//windowというか、ぶっちゃけ、OpenGLの状態そのものを管理する。

struct GLFWwindow;

class Window
{
		GLFWwindow *const window;

	  public:
		Window(uint32_t width = 640, uint32_t height = 480);

		bool is_shouldclose();

		void swapbuf();
		void pollevent();
		void waitevent(double sec);

		void clear();
};
