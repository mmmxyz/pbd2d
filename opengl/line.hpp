#pragma once

#include <cstdint>

#include "opengl/vertarray.hpp"

#include "mathfunc/vec.hpp"

class lineset
{
	  protected:
		vertarray varray;
		const uint32_t vsize;

	  public:
		lineset(vertex *data, uint32_t vsize);
		lineset(vertex *data, uint32_t vsize, uint32_t *idata, uint32_t isize);

		void setdata(vertex *data);
		void setposition(uint32_t index, float x, float y, float z);
		void setcolor(float r, float g, float b, float alpha);
		void draw() const;
};

//===

class line2d : public lineset
{
	  public:
		line2d(float p0[2], float p1[2], float color[4]);
		line2d(float p0x, float p0y, float p1x, float p1y, float color[4]);
		line2d(const fvec2 &v0, const fvec2 &v1, float color[4]);
		void setposition(float p0x, float p0y, float p1x, float p1y);
		void setposition(const fvec2 &v0, const fvec2 &v1);
};
