#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <memory>
#include <algorithm>

#include <unistd.h>

#include "mathfunc/vec.hpp"
#include "mathfunc/matrix.hpp"

#include "mathfunc/geometry.hpp"

#include "opengl/visualize.hpp"
#include "opengl/window.hpp"

#include "opengl/line.hpp"
#include "opengl/point.hpp"

constexpr uint32_t FPS = 120;
constexpr float dt = 1.0 / FPS;

constexpr float epsilon = 0.001;

class constraint
{
	  public:
		virtual void projection(float sn) = 0;
};

class collision : public constraint
{
		const fvec2 normal;
		const fvec2 q;
		const float ratio; //系に対する当該質点の質量比
		fvec2 &p;
		fvec2 &cn;

	  public:
		collision(fvec2 &x, fvec2 &cn, fvec2 q, float ratio, fvec2 n)
			: p(x), cn(cn), q(q), ratio(ratio), normal(n.normalize())
		{
		}

		void projection(float sn) override
		{
				float C = (p - q).dot(normal);
				if (C > epsilon) // not collide
						return;

				fvec2 dC = normal;

				fvec2 dposi = ratio * dC * (-(C / dC.sqlength()));

				p = p + sn * dposi;
				cn = cn + normal;
		}
};

class distance : public constraint
{
		const float L;
		const float m0, m1, w0, w1;
		fvec2 &p0;
		fvec2 &p1;

	  public:
		distance(float m0, float m1, fvec2 &p0, fvec2 &p1, float L)
			: m0(m0), m1(m1), w0(1.0 / m0), w1(1.0 / m1), p0(p0), p1(p1), L(L)
		{
		}

		void projection(float sn) override
		{
				float l = (p0 - p1).length();
				//if (std::abs(l - L) < epsilon)
				//		return;

				fvec2 dposi0 = (p1 - p0).normalize() * ((w0) / (w0 + w1)) * (l - L);
				fvec2 dposi1 = (p0 - p1).normalize() * ((w1) / (w0 + w1)) * (l - L);

				p0 = p0 + sn * dposi0;
				p1 = p1 + sn * dposi1;
		}
};

class bend : public constraint
{
		const float m0, m1, m2;
		const float w0, w1, w2;
		const float theta;
		fvec2 &p0;
		fvec2 &p1;
		fvec2 &p2;

	  public:
		bend(float m0, float m1, float m2, float theta, fvec2 &p0, fvec2 &p1, fvec2 &p2)
			: m0(m0), m1(m1), m2(m2), w0(1.0 / m0), w1(1.0 / m1), w2(1.0 / m2), theta(theta), p0(p0), p1(p1), p2(p2)
		{
		}

		void projection(float sn) override
		{
				fvec2 v0 = (p0 - p1).normalize();
				fvec2 v2 = (p2 - p1).normalize();

				float v0dv2 = v0.dot(v2);
				float C = std::acos(v0dv2) - theta;
				if (std::isnan(C))
						return;

				//std::cout << "hoge" << std::endl;

				//std::cout << std::acos(v0dv2) << std::endl;
				//std::cout << v0dv2 << std::endl;

				float darc = -1.0 / (std::sqrt(1 - v0dv2 * v0dv2));
				//std::cout << darc << std::endl;
				if (std::isnan(darc) || std::isinf(darc))
						return;

				fvec2 dC0 = (1.0 / (p0 - p1).length()) * (v2 - v0dv2 * v0);
				fvec2 dC2 = (1.0 / (p2 - p1).length()) * (v0 - v0dv2 * v2);
				fvec2 dC1 = -dC0 - dC2;

				//std::cout << dC0 << std::endl;
				//std::cout << dC2 << std::endl;
				//std::cout << dC1 << std::endl;

				dC0 = darc * dC0;
				dC2 = darc * dC2;
				dC1 = darc * dC1;

				//std::cout << dC0 << std::endl;
				//std::cout << dC2 << std::endl;
				//std::cout << dC1 << std::endl;

				fvec2 dp0 = -1.0 * C / (w0 * dC0.sqlength() + w2 * dC2.sqlength() + w1 * dC1.sqlength()) * w0 * dC0;
				fvec2 dp2 = -1.0 * C / (w0 * dC0.sqlength() + w2 * dC2.sqlength() + w1 * dC1.sqlength()) * w2 * dC2;
				fvec2 dp1 = -1.0 * C / (w0 * dC0.sqlength() + w2 * dC2.sqlength() + w1 * dC1.sqlength()) * w1 * dC1;

				//std::cout << dp0 << std::endl;

				p0 = p0 + sn * dp0;
				p2 = p2 + sn * dp2;
				p1 = p1 + sn * dp1;
		}
};

class bendsign : public constraint
{
		const float m0, m1, m2;
		const float w0, w1, w2;
		fvec2 &p0;
		fvec2 &p1;
		fvec2 &p2;

	  public:
		bendsign(float m0, float m1, float m2, fvec2 &p0, fvec2 &p1, fvec2 &p2)
			: m0(m0), m1(m1), m2(m2), w0(1.0 / m0), w1(1.0 / m1), w2(1.0 / m2), p0(p0), p1(p1), p2(p2)
		{
		}

		void projection(float sn) override
		{
				fvec2 v0 = (p0 - p1).normalize();
				fvec2 v2 = (p2 - p1).normalize();

				float C = v2.cross(v0);
				if (C > 0.0)
						return;

				fvec2 dC0(-v2.y, v2.x);
				fvec2 dC1(-v0.y + v2.y, v0.x - v2.x);
				fvec2 dC2(v0.y, -v0.x);

				fvec2 dp0 = -1.0 * C / (w0 * dC0.sqlength() + w2 * dC2.sqlength() + w1 * dC1.sqlength()) * w0 * dC0;
				fvec2 dp2 = -1.0 * C / (w0 * dC0.sqlength() + w2 * dC2.sqlength() + w1 * dC1.sqlength()) * w2 * dC2;
				fvec2 dp1 = -1.0 * C / (w0 * dC0.sqlength() + w2 * dC2.sqlength() + w1 * dC1.sqlength()) * w1 * dC1;

				p0 = p0 + sn * dp0;
				p2 = p2 + sn * dp2;
				p1 = p1 + sn * dp1;
		}
};

void wallcollisioncreate(fvec2 &tempposi, fvec2 &cn, std::vector<std::shared_ptr<constraint>> &v)
{

		if (tempposi.x < -0.8)
		{
				fvec2 q = tempposi;
				fvec2 normal(1.0, 0.0);
				q.x = -0.8;

				v.emplace_back(std::make_shared<collision>(tempposi, cn, q, 1.0, normal));
		}
		if (tempposi.x > 0.8)
		{
				fvec2 q = tempposi;
				fvec2 normal(-1.0, 0.0);
				q.x = 0.8;

				v.emplace_back(std::make_shared<collision>(tempposi, cn, q, 1.0, normal));
		}
		if (tempposi.y < -0.8)
		{
				fvec2 q = tempposi;
				fvec2 normal(0.0, 1.0);
				q.y = -0.8;

				v.emplace_back(std::make_shared<collision>(tempposi, cn, q, 1.0, normal));
		}
		if (tempposi.y > 0.8)
		{
				fvec2 q = tempposi;
				fvec2 normal(0.0, -1.0);
				q.y = 0.8;

				v.emplace_back(std::make_shared<collision>(tempposi, cn, q, 1.0, normal));
		}
}

class object
{
	  public:
		virtual void forwardEular() = 0;
		virtual void resetcn() = 0;
		virtual void refineposition() = 0;
		virtual void draw() = 0;
		virtual void createinternalconstraint(std::vector<std::shared_ptr<constraint>> &consset) = 0;
		virtual void createwallcollision(std::vector<std::shared_ptr<constraint>> &consset) = 0;
};

class rect : public object
{
	  public:
		const float m[4];
		const float width, height;

		fvec2 x[4];
		fvec2 p[4];
		fvec2 v[4];

		fvec2 cnormal[4];

		const float color[4];

		line2d lod[4]; //01,12,23,30
		point2d edge[4];

		rect(float m, fvec2 center, float width, float height, fvec2 gv, float omega, float *color)
			: m{m, m, m, m}, width(width), height(height), x{fvec2(center.x - width / 2.0, center.y - height / 2.0),
															 fvec2(center.x + width / 2.0, center.y - height / 2.0),
															 fvec2(center.x + width / 2.0, center.y + height / 2.0),
															 fvec2(center.x - width / 2.0, center.y + height / 2.0)},
			  color{color[0], color[1], color[2], color[3]}, lod{line2d(x[0].x, x[0].y, x[1].x, x[1].y, color),
																 line2d(x[1].x, x[1].y, x[2].x, x[2].y, color),
																 line2d(x[2].x, x[2].y, x[3].x, x[3].y, color),
																 line2d(x[3].x, x[3].y, x[0].x, x[0].y, color)},
			  edge{point2d(x[0].x, x[0].y, x[1].x, x[1].y, color), point2d(x[1].x, x[1].y, x[2].x, x[2].y, color),
				   point2d(x[2].x, x[2].y, x[3].x, x[3].y, color), point2d(x[3].x, x[3].y, x[0].x, x[0].y, color)}

		{
				float mat[4] = {0.0, -1.0, 1.0, 0.0};
				fmat2 rot90(mat);

				for (uint32_t i = 0; i < 4; i++)
				{
						fvec2 r = x[i] - center;
						r = rot90 * r;
						v[i] = gv + r * omega;
				}
		}

		void forwardEular() override
		{
				for (uint32_t i = 0; i < 4; i++)
				{
						v[i] = v[i] + fvec2(0.0, -9.8) * dt;

						p[i] = x[i] + v[i] * dt;
				}
		}

		void resetcn() override
		{
				for (uint32_t i = 0; i < 4; i++)
				{
						cnormal[i] = fvec2(0.0, 0.0);
				}
		}

		void refineposition() override
		{
				for (uint32_t i = 0; i < 4; i++)
				{

						v[i] = (p[i] - x[i]) / dt;
						x[i] = p[i];

						cnormal[i] = cnormal[i].normalize();

						if (cnormal[i].sqlength() > 0.0 && cnormal[i].dot(v[i]) < 0.0)
						{
								//restitution 0.8
								v[i] = v[i] - 0.9 * 2 * cnormal[i].dot(v[i]) * cnormal[i];
						}
				}
		}

		void draw() override
		{
				for (uint32_t i = 0; i < 4; i++)
				{
						uint32_t j = (i + 1) % 4;
						lod[i].setposition(x[i].x, x[i].y, x[j].x, x[j].y);
						edge[i].setposition(x[i].x, x[i].y, x[j].x, x[j].y);
						lod[i].draw();
						edge[i].draw();
				}
		}

		void createinternalconstraint(std::vector<std::shared_ptr<constraint>> &consset) override
		{
				for (uint32_t i = 0; i < 4; i++)
				{
						uint32_t j = (i + 1) % 4;
						float L = height;
						if (i % 2 == 0)
								L = width;

						consset.emplace_back(std::make_shared<distance>(m[i], m[j], p[i], p[j], L));

						uint32_t jj = (i + 2) % 4;
						consset.emplace_back(std::make_shared<bend>(m[i], m[j], m[jj], 3.14 * 0.5, p[i], p[j], p[jj]));
						consset.emplace_back(std::make_shared<bendsign>(m[i], m[j], m[jj], p[i], p[j], p[jj]));
				}
				//float L = std::sqrt(width * width + height * height);
				//consset.emplace_back(std::make_shared<distance>(m[0], m[2], p[0], p[2], L));
				//consset.emplace_back(std::make_shared<distance>(m[1], m[3], p[1], p[3], L));
		}

		void createwallcollision(std::vector<std::shared_ptr<constraint>> &consset) override
		{
				for (uint32_t i = 0; i < 4; i++)
				{
						wallcollisioncreate(p[i], cnormal[i], consset);
				}
		}
};

void drawminkowskidiff(rect &r0, rect &r1)
{
		float color[4] = {0.4, 0.4, 0.8, 1.0};

		for (uint32_t i = 0; i < 4; i++)
		{
				for (uint32_t j = 0; j < 4; j++)
				{
						point1d p(r0.p[i] - r1.p[j], color);
						p.draw();
				}
		}

		float zerocolor[4] = {0.0, 0.0, 0.0, 1.0};
		point1d zero(fvec2(), zerocolor);
		zero.draw();
}

struct support
{
		fvec2 vector;
		uint32_t index0, index1;
};

fvec2 supportfunc(const rect &r, const fvec2 d, uint32_t &index)
{
		fvec2 hoge = r.p[0];
		float max = hoge.dot(d);
		index = 0;
		for (uint32_t i = 1; i < 4; i++)
		{
				if (r.p[i].dot(d) >= max)
				{
						hoge = r.p[i];
						max = r.p[i].dot(d);
						index = i;
				}
		}

		return hoge;
}

support supportA_B(const rect &r0, const rect &r1, const fvec2 &d)
{
		fvec2 dn = d.normalize();
		uint32_t index0, index1;

		fvec2 vector = supportfunc(r0, dn, index0) - supportfunc(r1, -dn, index1);

		return {vector, index0, index1};
}

void collisionrect2rect(rect &r0, rect &r1, std::vector<std::shared_ptr<constraint>> &consset)
{

		fvec2 o = fvec2();

		fvec2 displace; //r0を動かすべき変位

		uint32_t l00, l01;
		uint32_t l10, l11;
		{
				//gjk
				support a0, a1, a2;

				{

						a0 = supportA_B(r0, r1, fvec2(1.0, 0.0));

						a1 = supportA_B(r0, r1, -a0.vector);

						if ((a1.vector - a0.vector).length() < epsilon)
								return;

						auto [dist, b] = distanceline2point(a0.vector, a1.vector, o);

						a2 = supportA_B(r0, r1, -b);

						if ((a2.vector - a0.vector).length() < epsilon || (a2.vector - a1.vector).length() < epsilon)
								return;

						if ((a1.vector - a0.vector).cross(a2.vector - a0.vector) < 0.0)
						{
								support temp = a0;
								a0 = a1;
								a1 = temp;
						}
				}

				struct supportset
				{
						float dist;
						fvec2 v;
						support p0;
						support p1;
				};

				std::vector<supportset> list;

				while (1)
				{

						auto [dist, c, l0, l1] = distancetriangle2point(a0.vector, a1.vector, a2.vector, o);

						if (dist <= 0.0)
						{

								auto [dist0, v0] = distancesegment2point(a0.vector, a1.vector, o);
								list.emplace_back(supportset{dist0, v0, a0, a1});
								auto [dist1, v1] = distancesegment2point(a1.vector, a2.vector, o);
								list.emplace_back(supportset{dist1, v1, a1, a2});
								auto [dist2, v2] = distancesegment2point(a2.vector, a0.vector, o);
								list.emplace_back(supportset{dist2, v2, a2, a0});

								break;
						}

						support b = supportA_B(r0, r1, -c);

						if ((a0.vector - b.vector).length() < epsilon || (a1.vector - b.vector).length() < epsilon ||
							(a2.vector - b.vector).length() < epsilon)
								return;

						if (epsilon < (a0.vector - l0).length() && epsilon < (a0.vector - l1).length())
						{
								a0 = b;
						}
						else if (epsilon < (a1.vector - l0).length() && epsilon < (a1.vector - l1).length())
						{
								a1 = b;
						}
						//else if (epsilon < (a2.vector - l0).length() < epsilon && epsilon < (a2.vector - l1).length())
						else
						{
								a2 = b;
						}

						//a0 = l0;
						//a1 = l1;
						//a2 = b;

						if ((a1.vector - a0.vector).cross(a2.vector - a0.vector) < 0.0)
						{
								support temp = a0;
								a0 = a1;
								a1 = temp;
						}
				}

				//std::sort(list.begin(), list.end(), [](supportset a, supportset b) -> bool { return a.dist < b.dist; });

				//for (auto x : list)
				//{
				//		std::cout << "!!!!!!" << std::endl;
				//		std::cout << x.dist << std::endl;
				//		std::cout << x.v << std::endl;
				//		std::cout << x.p0.vector << std::endl;
				//		std::cout << x.p1.vector << std::endl;
				//		std::cout << "!!!!!!" << std::endl;
				//}

				float tricolor[4] = {0.5, 0.0, 0.0, 1.0};
				float tricolorx[4] = {0.5, 0.4, 0.0, 1.0};
				line2d tri0(a0.vector, a1.vector, tricolorx);
				line2d tri1(a1.vector, a2.vector, tricolor);
				line2d tri2(a2.vector, a0.vector, tricolor);
				drawminkowskidiff(r0, r1);

				line2d cp(o, list[0].p0.vector, tricolorx);

				tri0.draw();
				tri1.draw();
				tri2.draw();
				cp.draw();

				//epa
				//uint32_t counter = 0;
				while (1)
				{
						//counter++;
						//if (counter > 100)
						//		return;

						std::sort(list.begin(), list.end(),
								  [](supportset a, supportset b) -> bool { return a.dist < b.dist; });

						fvec2 v = list[0].v;
						support l0 = list[0].p0;
						support l1 = list[0].p1;

						support b = supportA_B(r0, r1, v);

						//std::cout << "-----" << std::endl;
						//std::cout << v << std::endl;
						//std::cout << b.vector << std::endl;
						//std::cout << l0.vector << std::endl;
						//std::cout << l1.vector << std::endl;
						//std::cout << "-----" << std::endl;

						bool hoge = false;
						for (auto x : list)
						{

								if ((x.p0.index0 == b.index0 && x.p0.index1 == b.index1) ||
									(x.p1.index0 == b.index0 && x.p1.index1 == b.index1))
										hoge = true;
						}
						if (hoge)
						//if ((l0.vector - b.vector).length() < epsilon || (l1.vector - b.vector).length() < epsilon)
						{
								displace = -v;

								l00 = l0.index0;
								l01 = l0.index1;
								l10 = l1.index0;
								l11 = l1.index1;
								//fvec2 temp0 = findindex(r0, r1, l0, l00, l01);
								//if ((l0 - temp0).length() > epsilon)
								//		std::cout << "errRRRRRRRRRRRRRRR" << std::endl;

								//fvec2 temp1 = fetchsupportcomponent(r0, r1, l1, l10, l11);
								//if ((l1 - temp1).length() > epsilon)
								//		std::cout << "errRRRRRRRRRRRRRRR" << std::endl;

								//drawminkowskidiff(r0, r1);
								float color[4] = {0.0, 1.0, 0.0, 1.0};
								line2d disp(o, -displace, color);
								disp.draw();

								//std::cout << displace << std::endl;

								break;
						}

						list.erase(list.begin());

						auto [dist0, v0] = distancesegment2point(l0.vector, b.vector, o);
						auto [dist1, v1] = distancesegment2point(l1.vector, b.vector, o);

						list.emplace_back(supportset{dist0, v0, l0, b});
						list.emplace_back(supportset{dist1, v1, l1, b});

						line2d edge0(l0.vector, b.vector, tricolor);
						line2d edge1(l1.vector, b.vector, tricolor);
						edge0.draw();
						edge1.draw();
				}
		}

		float ratio = (r0.m[0] + r0.m[1] + r0.m[2] + r0.m[3]) /
					  (r0.m[0] + r0.m[1] + r0.m[2] + r0.m[3] + r1.m[0] + r1.m[1] + r1.m[2] + r1.m[3]);
		float lcolor[4] = {0.3, 0.3, 0.1, 1.0};

		{
				if (l00 == l10)
				{

						consset.emplace_back(std::make_shared<collision>(
							r0.p[l00], r0.cnormal[l00], r0.p[l00] + 1.0 * displace, 1.0f - ratio, displace));
				}
				else
				{
						consset.emplace_back(std::make_shared<collision>(
							r0.p[l00], r0.cnormal[l00], r0.p[l00] + 1.0 * displace, 0.5 * (1.0f - ratio), displace));
						consset.emplace_back(std::make_shared<collision>(
							r0.p[l10], r0.cnormal[l10], r0.p[l00] + 1.0 * displace, 0.5 * (1.0f - ratio), displace));
				}
				line2d hoge0(r0.p[l00], r0.p[l00] + 3.0 * displace, lcolor);
				point1d hogep0(r0.p[l00] + 3.0 * displace, lcolor);
				line2d hoge1(r0.p[l10], r0.p[l10] + 3.0 * displace, lcolor);
				point1d hogep1(r0.p[l10] + 3.0 * displace, lcolor);

				if (l01 == l11)
				{

						consset.emplace_back(std::make_shared<collision>(r1.p[l01], r1.cnormal[l01],
																		 r1.p[l01] - 1.0 * displace, ratio, -displace));
				}
				else
				{
						consset.emplace_back(std::make_shared<collision>(
							r1.p[l01], r1.cnormal[l01], r1.p[l01] - 1.0 * displace, 0.5 * (ratio), -displace));
						consset.emplace_back(std::make_shared<collision>(
							r1.p[l11], r1.cnormal[l11], r1.p[l01] - 1.0 * displace, 0.5 * (ratio), -displace));
				}
				line2d hoge2(r1.p[l01], r1.p[l01] - 3.0 * displace, lcolor);
				point1d hogep2(r1.p[l01] - 3.0 * displace, lcolor);
				line2d hoge3(r1.p[l11], r1.p[l11] - 3.0 * displace, lcolor);
				point1d hogep3(r1.p[l11] - 3.0 * displace, lcolor);

				hoge0.draw();
				hogep0.draw();
				hoge1.draw();
				hogep1.draw();
				hoge2.draw();
				hogep2.draw();
				hoge3.draw();
				hogep3.draw();
		}

		//for (uint32_t i = 0; i < 4; i++)
		//{
		//		consset.emplace_back(std::make_shared<collision>(r0.p[i], r0.cnormal[i], r0.p[i] + 1.0 * displace,
		//														 1.0f / ratio, displace));
		//		consset.emplace_back(
		//			std::make_shared<collision>(r1.p[i], r1.cnormal[i], r1.p[i] - 1.0 * displace, ratio, -displace));

		//		float lcolor[4] = {0.3, 0.3, 0.6, 0.0};

		//		line2d hoge0(r0.p[i], r0.p[i] + 1.0 * displace, lcolor);
		//		point1d hogep0(r0.p[i] + 1.0 * displace, lcolor);

		//		line2d hoge1(r1.p[i], r1.p[i] - 1.0 * displace, lcolor);
		//		point1d hogep1(r1.p[i] - 1.0 * displace, lcolor);

		//		hoge0.draw();
		//		hoge1.draw();
		//		hogep0.draw();
		//		hogep1.draw();
		//}
}

void timestep(rect &r0, rect &r1)
{

		//explicit eular
		//gravity acceleration = 9.8 m/ss
		r0.forwardEular();
		r1.forwardEular();

		//internal constraint
		std::vector<std::shared_ptr<constraint>> icv;
		r0.createinternalconstraint(icv);
		r1.createinternalconstraint(icv);

		//external constraint
		std::vector<std::shared_ptr<constraint>> ecv;

		//collision detection object and env
		r0.createwallcollision(ecv);
		r1.createwallcollision(ecv);

		//collision detection object and object
		//r0.createcollide(r1, ecv);
		//r1.createcollide(r0, ecv);
		collisionrect2rect(r0, r1, ecv);

		float sn = 1.0;

		//solver
		for (uint32_t loopcounter = 0; loopcounter < 10; loopcounter++)
		{
				for (auto &x : icv)
						x->projection(sn);

				for (auto &x : ecv)
						x->projection(sn);

				sn *= 0.8;
		}

		r0.refineposition();
		r0.resetcn();

		r1.refineposition();
		r1.resetcn();
}

int main(int argc, char const *argv[])
{
		using namespace std;

		Window mywindow = visualizeinit();

		shader mys = shader("../../../opengl/shadercode/vert.c", "../../../opengl/shadercode/frag.c");
		mys.useprogram();

		float wp0[2] = {-0.8, -0.8};
		float wp1[2] = {0.8, -0.8};
		float wp2[2] = {0.8, 0.8};
		float wp3[2] = {-0.8, 0.8};
		float wcolor[4] = {0.0, 0.0, 0.0, 1.0};
		line2d wl0(wp0, wp1, wcolor);
		line2d wl1(wp1, wp2, wcolor);
		line2d wl2(wp2, wp3, wcolor);
		line2d wl3(wp3, wp0, wcolor);

		float color0[4] = {1.0, 0.0, 0.0, 0.6};
		rect r0(1.0, fvec2(0.3, 0.3), 0.3, 0.4, fvec2(-1.0, 2.0), 10.0, color0);

		float color1[4] = {0.0, 1.0, 0.0, 0.6};
		rect r1(1.0, fvec2(-0.3, -0.20), 0.8, 0.1, fvec2(0.0, -3.0), 0.0, color1);

		double ctime = 0.0;
		double ltime = 0.0;

		uint32_t step = 0;
		while (!mywindow.is_shouldclose())
		{
				//clear buf
				mywindow.clear();

				ctime = getTime();
				if (ctime - ltime >= dt)
				{
						//physics
						timestep(r0, r1);
						ltime = getTime();
				}

				//rendering

				wl0.draw();
				wl1.draw();
				wl2.draw();
				wl3.draw();

				r0.draw();
				r1.draw();

				//swapbuff
				mywindow.swapbuf();
				//wait event
				mywindow.clearstatus();
				mywindow.pollevent();

				//mywindow.waitevent(100);

				step++;
				cout << "step: " << step;
				cout << "realtime : " << ctime << "  vtime : " << step * dt << endl;

				//while (mywindow.shouldwait())
				//{
				//		mywindow.clearstatus();
				//		mywindow.waitevent(100);
				//		std::cout << "waiting...." << std::endl;
				//}

				//for debug
				//mywindow.waitevent(100);
				//usleep(10000);
				//mywindow.waitevent(100);
		}

		return 0;
}
