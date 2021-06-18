#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <memory>

#include "mathfunc/vec.hpp"
#include "mathfunc/matrix.hpp"

#include "opengl/visualize.hpp"
#include "opengl/window.hpp"

#include "opengl/line.hpp"
#include "opengl/point.hpp"

constexpr uint32_t FPS = 120;
constexpr float dt = 1.0 / FPS;

constexpr float epsilon = 0.0001;

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
				if (C > 0.0) // not collide
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
				fvec2 dposi0 = (p1 - p0).normalize() * ((w0) / (w0 + w1)) * (l - L);
				fvec2 dposi1 = (p0 - p1).normalize() * ((w1) / (w0 + w1)) * (l - L);

				p0 = p0 + sn * dposi0;
				p1 = p1 + sn * dposi1;
		}
};

class bend : public constraint
{
		const float m0, m1;
		const float w0, w1;
		const float theta;
		fvec2 &p0;
		fvec2 &p1;
		fvec2 &ref;

	  public:
		bend(float m0, float m1, float theta, fvec2 &p0, fvec2 &p1, fvec2 &ref)
			: m0(m0), m1(m1), w0(1.0 / m0), w1(1.0 / m1), theta(theta), p0(p0), p1(p1), ref(ref)
		{
		}

		void projection(float sn) override
		{
				fvec2 v0 = (p0 - ref).normalize();
				fvec2 v1 = (p1 - ref).normalize();

				float v0dv1 = v0.dot(v1);
				float C = std::acos(v0dv1) - theta;
				if (std::isnan(C))
						return;
				float darc = -1.0 / (std::sqrt(1 - v0dv1 * v0dv1));

				fvec2 dC0 = darc * (1.0 / (p0 - ref).length()) * (v1 - v0dv1 * v0);
				fvec2 dC1 = darc * (1.0 / (p1 - ref).length()) * (v0 - v0dv1 * v1);

				fvec2 dp0 = -1.0 * C / (w0 * dC0.sqlength() + w1 * dC1.sqlength()) * w0 * dC0;
				fvec2 dp1 = -1.0 * C / (w0 * dC0.sqlength() + w1 * dC1.sqlength()) * w1 * dC1;

				p0 = p0 + dp0;
				p1 = p1 + dp1;
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

		bool is_collidebar(const fvec2 a0, const fvec2 a1, const fvec2 b0, const fvec2 b1) const
		{
				fvec2 v = b0 - a0;
				fvec2 w = b1 - a0;

				return v.cross(a1 - a0) * w.cross(a1 - a0) < -0.0;
		}

		void createcollide(rect &r, std::vector<std::shared_ptr<constraint>> &consset)
		{

				for (uint32_t i = 0; i < 4; i++)
				{
						uint32_t ii = (i + 1) % 4;
						for (uint32_t j = 0; j < 4; j++)
						{
								uint32_t jj = (j + 1) % 4;

								//2つのrectについて、4*4=16通りの辺の組み合わせで衝突判定
								if (is_collidebar(p[i], p[ii], r.p[j], r.p[jj]) &&
									is_collidebar(r.p[jj], r.p[j], p[i], p[ii])

								)
								{
										fmat2 hoge(p[ii] - p[i], r.p[j] - r.p[jj]);
										hoge = hoge.inverse();
										fvec2 ts = hoge * (r.p[j] - p[i]);
										//衝突ポイント collision point
										fvec2 cp = ts.x * (p[ii] - p[i]) + p[i];

										//test cpを表示 debug用
										float ccolor[4] = {0.0, 0.0, 0.0, 0.0};
										float cposi[2] = {cp.x, cp.y};
										point1d contact(cposi, ccolor);
										contact.draw();
										float wcolor[4] = {0.0, 0.0, 1.0, 0.0};
										//

										//貫通している法の頂点を押し出す。
										if ((p[ii] - p[i]).cross(r.p[j] - p[i]) < 0.0)
										{
												consset.emplace_back(
													std::make_shared<collision>(p[i], cnormal[i], cp, 0.5, cp - p[i]));

												//debug用
												float hogep[2] = {p[i].x, p[i].y};
												line2d contactlod(hogep, cposi, wcolor);
												contactlod.draw();
										}
										else
										{
												consset.emplace_back(std::make_shared<collision>(p[ii], cnormal[ii], cp,
																								 0.5, cp - p[ii]));

												//debug用
												float hogep[2] = {p[ii].x, p[ii].y};
												line2d contactlod(hogep, cposi, wcolor);
												contactlod.draw();
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
						consset.emplace_back(std::make_shared<bend>(m[i], m[jj], 3.14 * 0.5, p[i], p[jj], p[j]));
				}
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
		for (uint32_t loopcounter = 0; loopcounter < 30; loopcounter++)
		{
				for (auto &x : icv)
						x->projection(sn);

				for (auto &x : ecv)
						x->projection(sn);
				sn *= 1.0;
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
		rect r0(1.0, fvec2(0.0, 0.5), 0.3, 0.2, fvec2(0.0, -3.0), -0.05, color0);

		float color1[4] = {0.0, 1.0, 0.0, 0.0};
		rect r1(1.0, fvec2(0.1, -0.7), 0.4, 0.1, fvec2(0.0, 0.5), 0.08, color1);

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
		}

		return 0;
}
