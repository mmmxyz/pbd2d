#pragma once

#include <cstdint>

#include "opengl/vertarray.hpp"

class pointset
{
	  protected:
		vertarray varray;
		const uint32_t vsize;

	  public:
		pointset(vertex *data, uint32_t vsize);
		pointset(vertex *data, uint32_t vsize, uint32_t *idata, uint32_t isize);

		void setdata(vertex *data);
		void setposition(uint32_t index, float x, float y, float z);
		void setcolor(float r, float g, float b, float alpha);
		void draw() const;
};

//===

class point2d : public pointset
{
	  public:
		point2d(float p0[2], float p1[2], float color[4]);
		void setposition(float p0x, float p0y, float p1x, float p1y);
};

class point1d : public pointset
{

	  public:
		point1d(float p0[2], float color[4]);
		void setposition(float p0x, float p0y);
};
