#include <glad/glad.h>
#include <KHR/khrplatform.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <cstdint>

#include "opengl/vertarray.hpp"
#include "opengl/point.hpp"

#include "mathfunc/vec.hpp"

pointset::pointset(vertex *data, uint32_t vsize) : varray(vsize, data), vsize(vsize)
{
}

pointset::pointset(vertex *data, uint32_t vsize, uint32_t *idata, uint32_t isize)
	: varray(vsize, data, isize, idata), vsize(vsize)
{
}

void pointset::setdata(vertex *data)
{
		varray.setdata(data);
}

void pointset::setposition(uint32_t index, float x, float y, float z)
{
		varray.setposition(index, x, y, z);
}

void pointset::setcolor(float r, float g, float b, float alpha)
{
		varray.setcolor(r, g, b, alpha);
}

void pointset::draw() const
{
		varray.bind();
		glDrawArrays(GL_POINTS, 0, vsize);
		varray.unbind();
}

//====

point2d::point2d(float p0[2], float p1[2], float color[4]) : pointset(nullptr, 2)
{
		varray.setcolor(color[0], color[1], color[2], color[3]);
		varray.setposition(0, p0[0], p0[1], 0.0);
		varray.setposition(1, p1[0], p1[1], 0.0);
}

point2d::point2d(float p0x, float p0y, float p1x, float p1y, float color[4]) : pointset(nullptr, 2)
{
		varray.setcolor(color[0], color[1], color[2], color[3]);
		varray.setposition(0, p0x, p0y, 0.0);
		varray.setposition(1, p1x, p1y, 0.0);
}

point2d::point2d(const fvec2 &a0, const fvec2 &a1, float color[4]) : pointset(nullptr, 2)
{
		varray.setcolor(color[0], color[1], color[2], color[3]);
		varray.setposition(0, a0.x, a0.y, 0.0);
		varray.setposition(1, a1.x, a1.y, 0.0);
}

void point2d::setposition(float p0x, float p0y, float p1x, float p1y)
{
		varray.setposition(0, p0x, p0y, 0.0);
		varray.setposition(1, p1x, p1y, 0.0);
}

void point2d::setposition(const fvec2 &p0, const fvec2 &p1)
{
		varray.setposition(0, p0.x, p0.y, 0.0);
		varray.setposition(1, p1.x, p1.y, 0.0);
}

//====

point1d::point1d(float p0[2], float color[4]) : pointset(nullptr, 1)
{
		varray.setcolor(color[0], color[1], color[2], color[3]);
		varray.setposition(0, p0[0], p0[1], 0.0);
}

point1d::point1d(const fvec2 &v, float color[4]) : pointset(nullptr, 1)
{

		varray.setcolor(color[0], color[1], color[2], color[3]);
		varray.setposition(0, v.x, v.y, 0.0);
}

void point1d::setposition(float p0x, float p0y)
{
		varray.setposition(0, p0x, p0y, 0.0);
}

void point1d::setposition(const fvec2 &v)
{
		varray.setposition(0, v.x, v.y, 0.0);
}
