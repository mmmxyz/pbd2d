#include <vector>

#include "mathfunc/vec.hpp"
#include "mathfunc/matrix.hpp"
#include "mathfunc/geometry.hpp"

static float epsilon = 0.00001;

closestdp distanceline2point(const fvec2 &a0, const fvec2 &a1, const fvec2 &b)
{
		fvec2 n = a1 - a0;
		float t = -(a0 - b).dot(n) / n.dot(n);

		fvec2 v = t * n + a0;
		float dist = (v - b).length();
		return {dist, v};
}

closestdp distancesegment2point(const fvec2 &a0, const fvec2 &a1, const fvec2 &b)
{
		fvec2 n = a1 - a0;
		float t = -(a0 - b).dot(n) / n.dot(n);

		if (0 <= t && t <= 1.0)
		{
				fvec2 v = t * n + a0;
				float dist = (v - b).length();

				return {dist, v};
		}
		else if (t < 0.0)
		{
				fvec2 v = a0;
				float dist = (v - b).length();

				return {dist, v};
		}
		else
		{
				fvec2 v = a1;
				float dist = (v - b).length();

				return {dist, v};
		}
}

closestdpp distanceline2line(const fvec2 &a0, const fvec2 &a1, const fvec2 &b0, const fvec2 &b1)
{

		fvec2 na = a1 - a0;
		fvec2 nb = b1 - b0;

		fvec2 a, b;
		float dist;

		if (na.sqlength() * nb.sqlength() - (na.dot(nb)) * (na.dot(nb)) < epsilon)
		{
				float s = (a0.dot(na) - b0.dot(na)) / na.dot(nb);
				a = a0;
				b = s * (a1 - a0) + a0;

				dist = (a - b).length();
				return {dist, a, b};
		}

		fmat2 hoge({na.sqlength(), -na.dot(nb), na.dot(nb), -nb.sqlength()});
		hoge = hoge.inverse();
		fvec2 ts = hoge * fvec2({-a0.dot(na) + b0.dot(na), -a0.dot(nb) + b0.dot(nb)});

		a = ts.x * (a1 - a0) + a0;
		b = ts.y * (b1 - b0) + b0;

		dist = (a - b).length();

		return {dist, a, b};
}

closestdpp distancesegment2segment(const fvec2 &a0, const fvec2 &a1, const fvec2 &b0, const fvec2 &b1)
{

		fvec2 na = a1 - a0;
		fvec2 nb = b1 - b0;

		fvec2 a, b;
		float dist;

		if (na.sqlength() * nb.sqlength() - (na.dot(nb)) * (na.dot(nb)) < epsilon)
		{

				auto [dist0, v0] = distancesegment2point(a0, a1, b0);
				auto [dist1, v1] = distancesegment2point(a0, a1, b1);
				auto [dist2, v2] = distancesegment2point(b0, b1, a0);
				auto [dist3, v3] = distancesegment2point(b0, b1, a1);

				if (dist0 <= dist1 + epsilon && dist0 <= dist2 + epsilon && dist0 <= dist3 + epsilon)
						return {dist0, v0, b0};
				else if (dist1 <= dist0 + epsilon && dist1 <= dist2 + epsilon && dist1 <= dist3 + epsilon)
						return {dist1, v1, b1};
				else if (dist2 <= dist0 + epsilon && dist2 <= dist1 + epsilon && dist2 <= dist3 + epsilon)
						return {dist2, a0, v2};
				else
						return {dist3, a1, v3};
		}

		fmat2 hoge({na.sqlength(), -na.dot(nb), na.dot(nb), -nb.sqlength()});
		hoge = hoge.inverse();
		fvec2 ts = hoge * fvec2({-a0.dot(na) + b0.dot(na), -a0.dot(nb) + b0.dot(nb)});
		//std::cout << ts << std::endl;

		if (0.0 <= ts.x && ts.x <= 1.0 && 0.0 <= ts.y && ts.y <= 1.0)
		{
				a = ts.x * (a1 - a0) + a0;
				b = ts.y * (b1 - b0) + b0;
		}
		else if (0.0 <= ts.x && ts.x <= 1.0)
		{
				if (ts.y < 0.0)
				{
						auto [dist, v] = distancesegment2point(a0, a1, b0);
						return {dist, v, b0};
				}
				else
				{
						auto [dist, v] = distancesegment2point(a0, a1, b1);
						return {dist, v, b0};
				}
		}
		else if (0.0 <= ts.y && ts.y <= 1.0)
		{
				if (ts.x < 0.0)
				{
						auto [dist, v] = distancesegment2point(b0, b1, a0);
						return {dist, a0, v};
				}
				else
				{
						auto [dist, v] = distancesegment2point(b0, b1, a1);
						return {dist, a1, v};
				}
		}
		else
		{
				if (ts.x < 0.0)
						a = a0;
				else
						a = a1;

				if (ts.y < 0.0)
						b = b0;
				else
						b = b1;
		}

		auto [dist1, v1] = distancesegment2point(a0, a1, b);
		auto [dist2, v2] = distancesegment2point(b0, b1, a);

		if (dist1 < dist2)
				return {dist1, v1, b};
		else
				return {dist2, a, v2};
}

closestdppp distancetriangle2point(const fvec2 &a0, const fvec2 &a1, const fvec2 &a2, const fvec2 &b)
{

		auto [dist0, v0] = distancesegment2point(a0, a1, b);
		auto [dist1, v1] = distancesegment2point(a1, a2, b);
		auto [dist2, v2] = distancesegment2point(a2, a0, b);

		if ((a0 - a1).rot().dot(b - a0) < -epsilon && (a1 - a2).rot().dot(b - a1) < -epsilon &&
			(a2 - a0).rot().dot(b - a2) < -epsilon)
		{
				if (dist0 <= dist1 + epsilon && dist0 <= dist2 + epsilon)
						return {-1.0, v0, a0, a1};
				else if (dist1 <= dist0 + epsilon && dist1 <= dist2 + epsilon)
						return {-1.0, v1, a1, a2};
				else
						return {-1.0, v2, a2, a0};
		}
		else
		{
				if (dist0 <= dist1 + epsilon && dist0 <= dist2 + epsilon)
						return {dist0, v0, a0, a1};
				else if (dist1 <= dist0 + epsilon && dist1 <= dist2 + epsilon)
						return {dist1, v1, a1, a2};
				else
						return {dist2, v2, a2, a0};
		}
}

closestdppp distancetriangleedge2point(const fvec2 &a0, const fvec2 &a1, const fvec2 &a2, const fvec2 &b)
{

		auto [dist0, v0] = distancesegment2point(a0, a1, b);
		auto [dist1, v1] = distancesegment2point(a1, a2, b);
		auto [dist2, v2] = distancesegment2point(a2, a0, b);

		if (dist0 <= dist1 + epsilon && dist0 <= dist2 + epsilon)
				return {dist0, v0, a0, a1};
		else if (dist1 <= dist0 + epsilon && dist1 <= dist2 + epsilon)
				return {dist1, v1, a1, a2};
		else
				return {dist2, v2, a2, a0};
}
