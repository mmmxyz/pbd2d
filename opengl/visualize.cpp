#include "3rdparty/glad/include/glad/glad.h"
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <iostream>
#include <fstream>
#include <vector>

#include "opengl/visualize.hpp"
#include "opengl/window.hpp"

Window visualizeinit()
{

		if (glfwInit() == GL_FALSE)
		{
				std::cout << "fail to initialize glfw" << std::endl;
				exit(1);
		}

		atexit(glfwTerminate);

		uint32_t width = 1024;
		uint32_t height = 1024;

		Window window = Window(width, height);

		if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
		{
				std::cout << "fail to create OpenGL context" << std::endl;
				exit(1);
		}

		//Version
		std::cout << "OpenGL Version : " << glGetString(GL_VERSION) << std::endl;
		std::cout << "OpenGL Version : " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
		std::cout << std::endl;

		//mics OpenGL setting
		glPointSize(10.0f);

		return window;
}

double getTime()
{
		return glfwGetTime();
}

//shader

void printcompilelog(GLuint shader)
{
		GLsizei bufSize;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &bufSize);

		if (bufSize > 1)
		{
				GLchar *logbuf = new GLchar[bufSize];
				GLsizei length;
				glGetShaderInfoLog(shader, bufSize, &length, logbuf);
				std::cerr << logbuf << std::endl;
				delete[] logbuf;
		}
}

void printlinklog(GLuint program)
{
		GLsizei bufSize;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &bufSize);

		if (bufSize > 1)
		{
				GLchar *logbuf = new GLchar[bufSize];
				GLsizei length;
				glGetProgramInfoLog(program, bufSize, &length, logbuf);
				std::cerr << logbuf << std::endl;
				delete[] logbuf;
		}
}

bool loadcode(const char *name, GLchar *&buffer)
{
		using namespace std;
		ifstream file(name, ios::binary);
		if (file.fail())
		{
				cout << "fail to open file: " << name << endl;
				return false;
		}

		file.seekg(0L, ios::end);
		GLsizei length = static_cast<GLsizei>(file.tellg());

		buffer = new GLchar[length + 1];
		file.seekg(0L, ios::beg);
		file.read(buffer, length);
		buffer[length] = '\0';

		file.close();

		////printcode
		//for (uint32_t i = 0; i <= length; i++)
		//{
		//		cout << buffer[i];
		//}
		//cout << endl;

		return true;
}

GLuint linkprogram(const GLchar *vsd, const GLchar *fsd)
{
		const GLuint program(glCreateProgram());

		//バーテクスシェーダ
		const GLuint vobj(glCreateShader(GL_VERTEX_SHADER));
		glShaderSource(vobj, 1, &vsd, nullptr);
		glCompileShader(vobj);

		//コンパイルの成否
		GLint status;
		glGetShaderiv(vobj, GL_COMPILE_STATUS, &status);
		if (status == GL_FALSE)
				std::cerr << "compile error in vertex shader" << std::endl;
		else
				std::cout << "vertex shader compilation is OK" << std::endl;

		printcompilelog(vobj);

		if (status)
				glAttachShader(program, vobj);
		glDeleteShader(vobj);

		//===================================

		//フラグメントシェーダ
		const GLuint fobj(glCreateShader(GL_FRAGMENT_SHADER));
		glShaderSource(fobj, 1, &fsd, nullptr);
		glCompileShader(fobj);

		glGetShaderiv(fobj, GL_COMPILE_STATUS, &status);
		if (status == GL_FALSE)
				std::cerr << "compile error in fragment shader" << std::endl;
		else
				std::cout << "fragment shader compilation is OK" << std::endl;

		printcompilelog(fobj);
		if (status)
				glAttachShader(program, fobj);
		glDeleteShader(fobj);

		//===================================

		//オブジェクトのリンク
		glBindAttribLocation(program, 0, "position");
		glBindAttribLocation(program, 1, "color");
		glBindFragDataLocation(program, 0, "fragment");
		glLinkProgram(program);

		glGetProgramiv(program, GL_LINK_STATUS, &status);
		if (status == GL_FALSE)
				std::cerr << "link error" << std::endl;
		else
				std::cout << "link is OK" << std::endl;

		printlinklog(program);

		if (status)
				return program;
		glDeleteProgram(program);

		return 0;
}

shader::shader(const char *vsname, const char *fsname)
{
		GLchar *vcode = nullptr;
		GLchar *fcode = nullptr;

		if (loadcode(vsname, vcode) && loadcode(fsname, fcode))
		{
				program = linkprogram(vcode, fcode);
		}
		else
		{
				program = 0;
		}

		delete[] vcode;
		delete[] fcode;
}

void shader::useprogram() const
{
		glUseProgram(program);
}
