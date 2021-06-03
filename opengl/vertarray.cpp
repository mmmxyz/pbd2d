#include <GL/glew.h>

#include "opengl/vertarray.hpp"

vertarray::vertarray(uint32_t size, vertex *data, uint32_t isize, uint32_t *ilist) : size(size), va(new vertex[size])
{

		if (data != nullptr)
		{
				for (uint32_t i = 0; i < size; i++)
				{
						va[i] = data[i];
				}
		}
		else
		{
				for (uint32_t i = 0; i < size; i++)
				{
						va[i].position[0] = 0.0f;
						va[i].position[1] = 0.0f;
						va[i].position[2] = 0.0f;
						va[i].color[0] = 0.0f;
						va[i].color[1] = 0.0f;
						va[i].color[2] = 0.0f;
				}
		}

		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_STATIC_DRAW);

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex *>(0)->position);
		glEnableVertexAttribArray(0);

		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex *>(0)->color);
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, 0);

		if (ilist != nullptr)
		{
				glGenBuffers(1, &ibo);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
				glBufferData(GL_ELEMENT_ARRAY_BUFFER, isize * sizeof(GLuint), ilist, GL_STATIC_DRAW);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		}

		glBindVertexArray(0);
}

void vertarray::setdata(vertex *data)
{
		for (uint32_t i = 0; i < size; i++)
		{
				va[i] = data[i];
		}

		glBindVertexArray(vao);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glBindVertexArray(0);
}

void vertarray::setposition(uint32_t index, float x, float y, float z)
{

		va[index].position[0] = x;
		va[index].position[1] = y;
		va[index].position[2] = z;

		glBindVertexArray(vao);

		//指定された位置だけの更新をするべきである。TODO
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glBindVertexArray(0);
}

void vertarray::setcolor(float r, float g, float b, float alpha)
{

		for (uint32_t i = 0; i < size; i++)
		{
				va[i].color[0] = r;
				va[i].color[1] = g;
				va[i].color[2] = b;
				va[i].color[3] = alpha;
		}

		glBindVertexArray(vao);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glBindVertexArray(0);
}

vertarray::~vertarray()
{
		delete[] va;
		glDeleteBuffers(1, &vao);
		glDeleteBuffers(1, &vbo);
		if (ibo != 0)
				glDeleteBuffers(1, &ibo);
}

void vertarray::bind() const
{
		glBindVertexArray(vao);
}
void vertarray::unbind() const
{
		glBindVertexArray(0);
}
