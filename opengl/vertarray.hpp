#pragma once

#include <cstdint>

struct vertex
{
		float position[3];
		float color[4];
};

class vertarray
{
		uint32_t size;
		vertex *va;

		uint32_t vao, vbo, ibo = 0;

	  public:
		vertarray(uint32_t size, vertex *data, uint32_t isize = 0, uint32_t *ilist = nullptr);
		~vertarray();

		void setdata(vertex *data);
		void setposition(uint32_t index, float x, float y, float z);
		void setcolor(float r, float g, float b, float alpha);
		void bind() const;
		void unbind() const;
};
