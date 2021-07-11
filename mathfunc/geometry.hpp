#pragma once

#include <vector>

#include "mathfunc/vec.hpp"
#include "mathfunc/matrix.hpp"

struct closestdp
{
		float dist;
		fvec2 p;
};

struct closestdpp
{
		float dist;
		fvec2 p0;
		fvec2 p1;
};

struct closestdppp
{
		float dist;
		fvec2 p0;
		fvec2 p1;
		fvec2 p2;
};

closestdp distanceline2point(const fvec2 &a0, const fvec2 &a1, const fvec2 &b);

closestdp distancesegment2point(const fvec2 &a0, const fvec2 &a1, const fvec2 &b);

closestdpp distanceline2line(const fvec2 &a0, const fvec2 &a1, const fvec2 &b0, const fvec2 &b1);

closestdpp distancesegment2segment(const fvec2 &a0, const fvec2 &a1, const fvec2 &b0, const fvec2 &b1);

closestdppp distancetriangle2point(const fvec2 &a0, const fvec2 &a1, const fvec2 &a2, const fvec2 &b);

closestdppp distancetriangleedge2point(const fvec2 &a0, const fvec2 &a1, const fvec2 &a2, const fvec2 &b);
