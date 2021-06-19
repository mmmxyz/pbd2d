#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <memory>

#include <unistd.h>

#include "mathfunc/vec.hpp"
#include "mathfunc/matrix.hpp"

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

		bool is_collidebarline(const fvec2 a0, const fvec2 a1, const fvec2 b0, const fvec2 b1) const
		{
				//line a0,a1
				//bar b0,b1
				fvec2 v = b0 - a0;
				fvec2 w = b1 - a0;

				return v.cross(a1 - a0) * w.cross(a1 - a0) < 0.0;
		}

		bool is_collidebarbar(const fvec2 a0, const fvec2 a1, const fvec2 b0, const fvec2 b1) const
		{

				return is_collidebarline(a0, a1, b0, b1) && is_collidebarline(b0, b1, a0, a1);
		}

		bool is_collidebarray(const fvec2 a0, const fvec2 a1, const fvec2 b0, const fvec2 b1) const
		{
				//ray a0 to a1
				//bar b0,b1

				if (is_collidebarline(a0, a1, b0, b1))
				{
						fmat2 hoge(a1 - a0, b0 - b1);
						hoge = hoge.inverse();
						fvec2 ts = hoge * (b0 - a0);
						if (ts.x >= 0.0)
								return true;
				}
				return false;
		}

		bool is_inside(const fvec2 &v)
		{

				fvec2 center = (p[0] + p[1] + p[2] + p[3]) / 4.0;

				fvec2 v0 = v - p[0];
				fvec2 x = p[1] - p[0];
				fvec2 y = p[3] - p[0];
				float xc = x.dot(v0);
				float yc = y.dot(v0);

				if (epsilon < xc && xc < width * width - epsilon && epsilon < yc && yc < height * height - epsilon)
						return true;
				else
						return false;
		}

		void createcollide(rect &r, std::vector<std::shared_ptr<constraint>> &consset)
		{

				float ratio = (m[0] + m[1] + m[2] + m[3]) / (r.m[0] + r.m[1] + r.m[2] + r.m[3]);

				fvec2 rgv = (r.v[0] + r.v[1] + r.v[2] + r.v[3]) / 4;

				for (uint32_t i = 0; i < 4; i++)
				{
						if (r.is_inside(p[i]))
						{

								float ccolor[4] = {0.0, 1.0, 1.0, 0.0};
								float cposi[2] = {p[i].x - 0.001f, p[i].y + 0.001f};
								point1d contact(cposi, ccolor);
								contact.draw();

								float wcolor[4] = {1.0, 1.0, 0.0, 0.0};
								float hogep[2] = {p[i].x, p[i].y};
								float hogepv[2] = {p[i].x - (v[i].x - rgv.x), p[i].y - (v[i].y - rgv.y)};
								line2d contactlod(hogep, hogepv, wcolor);
								contactlod.draw();

								for (uint32_t j = 0; j < 4; j++)
								{
										uint32_t jj = (j + 1) % 4;

										if (is_collidebarray(p[i], p[i] - (v[i] - rgv), r.p[j], r.p[jj]))
										{

												fvec2 normal((r.p[jj] - r.p[j]).y, -(r.p[jj] - r.p[j]).x);
												normal = normal.normalize();

												fmat2 hoge(normal, r.p[j] - r.p[jj]);
												hoge = hoge.inverse();
												fvec2 ts = hoge * (r.p[j] - p[i]);
												////衝突ポイント collision point
												fvec2 cp = ts.x * normal + p[i];

												//test
												float ccolor[4] = {0.0, 0.0, 0.0, 0.0};
												float cposi[2] = {cp.x, cp.y};
												point1d contact(cposi, ccolor);
												contact.draw();

												consset.emplace_back(
													std::make_shared<collision>(p[i], cnormal[i], cp, ratio, normal));
												consset.emplace_back(std::make_shared<collision>(
													r.p[j], r.cnormal[j], cp, 0.5 / ratio, -normal));
												consset.emplace_back(std::make_shared<collision>(
													r.p[jj], r.cnormal[jj], cp, 0.5 / ratio, -normal));
										}
								}
						}
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
								v[i] = v[i] - 2 * cnormal[i].dot(v[i]) * cnormal[i];
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
		r0.createcollide(r1, ecv);
		r1.createcollide(r0, ecv);

		float sn = 1.0;

		//solver
		for (uint32_t loopcounter = 0; loopcounter < 10; loopcounter++)
		{
				for (auto &x : icv)
						x->projection(sn);

				for (auto &x : ecv)
						x->projection(sn);
				sn *= 0.9;
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
		float wcolor[4] = {0.0, 0.0, 0.0, 0.0};
		line2d wl0(wp0, wp1, wcolor);
		line2d wl1(wp1, wp2, wcolor);
		line2d wl2(wp2, wp3, wcolor);
		line2d wl3(wp3, wp0, wcolor);

		float color0[4] = {1.0, 0.0, 0.0, 0.0};
		rect r0(2.0, fvec2(-0.4, 0.5), 0.3, 0.4, fvec2(2.0, 5.0), -50.5, color0);

		float color1[4] = {0.0, 1.0, 0.0, 0.0};
		rect r1(1.0, fvec2(0.4, -0.6), 0.4, 0.1, fvec2(-2.0, 4.5), 50.8, color1);

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
