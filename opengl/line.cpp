#include <glad/glad.h>
#include <KHR/khrplatform.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <cstdint>

#include "opengl/line.hpp"
#include "opengl/vertarray.hpp"

#include "mathfunc/vec.hpp"

lineset::lineset(vertex *data, uint32_t vsize) : varray(vsize, data), vsize(vsize), isize(0)
{
		isindexed = false;
}

lineset::lineset(vertex *data, uint32_t vsize, uint32_t *idata, uint32_t isize)
	: varray(vsize, data, isize, idata), vsize(vsize), isize(isize)
{
		isindexed = true;
}

void lineset::setdata(vertex *data)
{
		varray.setdata(data);
}

void lineset::setposition(uint32_t index, float x, float y, float z)
{
		varray.setposition(index, x, y, z);
}

void lineset::setcolor(float r, float g, float b, float alpha)
{
		varray.setcolor(r, g, b, alpha);
}

void lineset::draw() const
{
		if (isindexed)
		{
				varray.bind();
				glDrawElements(GL_LINE_STRIP, isize, GL_UNSIGNED_INT, (void *)0);
				varray.unbind();
		}
		else
		{
				varray.bind();
				glDrawArrays(GL_LINE_STRIP, 0, vsize);
				varray.unbind();
		}
}

//=====

line2d::line2d(float p0[2], float p1[2], float color[4]) : lineset(nullptr, 2)
{
		varray.setcolor(color[0], color[1], color[2], color[3]);
		varray.setposition(0, p0[0], p0[1], 0.0);
		varray.setposition(1, p1[0], p1[1], 0.0);
}

line2d::line2d(float p0x, float p0y, float p1x, float p1y, float color[4]) : lineset(nullptr, 2)
{
		varray.setcolor(color[0], color[1], color[2], color[3]);
		varray.setposition(0, p0x, p0y, 0.0);
		varray.setposition(1, p1x, p1y, 0.0);
}

line2d::line2d(const fvec2 &v0, const fvec2 &v1, float color[4]) : lineset(nullptr, 2)
{

		varray.setcolor(color[0], color[1], color[2], color[3]);
		varray.setposition(0, v0.x, v0.y, 0.0);
		varray.setposition(1, v1.x, v1.y, 0.0);
}

void line2d::setposition(float p0x, float p0y, float p1x, float p1y)
{
		varray.setposition(0, p0x, p0y, 0.0);
		varray.setposition(1, p1x, p1y, 0.0);
}

void line2d::setposition(const fvec2 &v0, const fvec2 &v1)
{
		varray.setposition(0, v0.x, v0.y, 0.0);
		varray.setposition(1, v1.x, v1.y, 0.0);
}
