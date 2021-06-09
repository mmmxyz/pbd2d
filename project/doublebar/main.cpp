#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>

#include "mathfunc/vec.hpp"
#include "mathfunc/matrix.hpp"

#include "opengl/visualize.hpp"
#include "opengl/window.hpp"

#include "opengl/line.hpp"
#include "opengl/point.hpp"

constexpr uint32_t FPS = 120;
constexpr float dt = 1.0 / FPS;

constexpr float epsilon = 0.0001;

class collisionsolver
{
		const fvec2 normal;
		const fvec2 q;
		fvec2 &p;
		fvec2 &cn;

	  public:
		collisionsolver(fvec2 &x, fvec2 &cn, fvec2 q, fvec2 n) : p(x), cn(cn), q(q), normal(n.normalize())
		{
		}

		void projection(float sn)
		{
				float C = (p - q).dot(normal);
				if (C > 0.0) // not collide
						return;

				fvec2 dC = normal;

				fvec2 dposi = dC * (-(C / dC.sqlength()));

				p = p + sn * dposi;
				cn = cn + normal;
		}
};

class bar
{
	  public:
		const float m0, m1, L;
		const float w0, w1;
		fvec2 p0, p1, v0, v1;
		fvec2 temp0, temp1;
		fvec2 cnormal0, cnormal1;
		const float color[4];
		line2d lod;
		point2d edge;

		bar(float m0, float m1, float L, fvec2 p0, fvec2 p1, fvec2 v0, fvec2 v1, float *color)
			: m0(m0), m1(m1), w0(1 / m0), w1(1 / m1), L(L), p0(p0), p1(p1), temp0(p0), temp1(p1), v0(v0),
			  v1(v1), color{color[0], color[1], color[2], color[3]}, lod(p0.x, p0.y, p1.x, p1.y, color),
			  edge(p0.x, p0.y, p1.x, p1.y, color), cnormal0(0.0, 0.0), cnormal1(0.0, 0.0)
		{
		}

		void forwardEular()
		{
				v0 = v0 + fvec2(0.0, -9.8) * dt;
				v1 = v1 + fvec2(0.0, -9.8) * dt;

				temp0 = p0 + v0 * dt;
				temp1 = p1 + v1 * dt;
		}

		void resetcn()
		{
				cnormal0 = fvec2(0.0, 0.0);
				cnormal1 = fvec2(0.0, 0.0);
		}

		void refineposition()
		{
				v0 = (temp0 - p0) / dt;
				v1 = (temp1 - p1) / dt;

				p0 = temp0;
				p1 = temp1;

				cnormal0 = cnormal0.normalize();
				cnormal1 = cnormal1.normalize();

				if (cnormal0.sqlength() > 0.0 && cnormal0.dot(v0) < 0.0)
						v0 = v0 - 2 * cnormal0.dot(v0) * cnormal0;
				if (cnormal1.sqlength() > 0.0 && cnormal1.dot(v1) < 0.0)
						v1 = v1 - 2 * cnormal1.dot(v1) * cnormal1;
		}

		void draw()
		{
				lod.setposition(p0.x, p0.y, p1.x, p1.y);
				edge.setposition(p0.x, p0.y, p1.x, p1.y);
				lod.draw();
				edge.draw();
		}

		void projectdistance(float sn)
		{
				float l = (temp0 - temp1).length();
				fvec2 dposi0 = (temp1 - temp0).normalize() * ((w0) / (w0 + w1)) * (l - L);
				fvec2 dposi1 = (temp0 - temp1).normalize() * ((w1) / (w0 + w1)) * (l - L);

				temp0 = temp0 + sn * dposi0;
				temp1 = temp1 + sn * dposi1;
		}

		bool is_collide(const bar &a) const
		{
				fvec2 v = a.temp0 - temp0;
				fvec2 w = a.temp1 - temp0;

				return v.cross(temp1 - temp0) * w.cross(temp1 - temp0) < -epsilon;
		}

		void createcollision(bar &a, std::vector<collisionsolver> &v)
		{
				fmat2 hoge(a.temp1 - a.temp0, -temp1 + temp0);
				hoge = hoge.inverse();
				fvec2 ts = hoge * (temp0 - a.temp0);

				fvec2 cp = ts.y * (temp1 - temp0) + temp0;

				//if (std::abs(ts.y - 0.5) > 0.4)
				//		return;

				//test
				float ccolor[4] = {0.0, 0.0, 0.0, 0.0};
				float cposi[2] = {cp.x, cp.y};
				point1d contact(cposi, ccolor);
				contact.draw();
				////

				fvec2 gv = (m0 * v0 + m1 * v1) / (m0 + m1);
				fvec2 agv = (a.m0 * a.v0 + a.m1 * a.v1) / (a.m0 + a.m1);

				fvec2 gp = (m0 * temp0 + m1 * temp1) / (m0 + m1);

				v.emplace_back(collisionsolver(temp0, cnormal0, cp, (-gv)));
				v.emplace_back(collisionsolver(temp1, cnormal1, cp, (-gv)));
		}
};

void wallcollisionsolver(fvec2 &tempposi, fvec2 &cn, std::vector<collisionsolver> &v)
{

		if (tempposi.x < -0.8)
		{
				fvec2 q = tempposi;
				fvec2 normal(1.0, 0.0);
				q.x = -0.8;

				v.emplace_back(collisionsolver(tempposi, cn, q, normal));
		}
		if (tempposi.x > 0.8)
		{
				fvec2 q = tempposi;
				fvec2 normal(-1.0, 0.0);
				q.x = 0.8;

				v.emplace_back(collisionsolver(tempposi, cn, q, normal));
		}
		if (tempposi.y < -0.8)
		{
				fvec2 q = tempposi;
				fvec2 normal(0.0, 1.0);
				q.y = -0.8;

				v.emplace_back(collisionsolver(tempposi, cn, q, normal));
		}
		if (tempposi.y > 0.8)
		{
				fvec2 q = tempposi;
				fvec2 normal(0.0, -1.0);
				q.y = 0.8;

				v.emplace_back(collisionsolver(tempposi, cn, q, normal));
		}
}

void timestep(bar &b0, bar &b1)
{
		//explicit eular
		//gravity acceleration = 9.8 m/ss
		b0.forwardEular();
		b1.forwardEular();

		//create collisionsolver
		std::vector<collisionsolver> csolver;
		wallcollisionsolver(b0.temp0, b0.cnormal0, csolver);
		wallcollisionsolver(b0.temp1, b0.cnormal1, csolver);

		wallcollisionsolver(b1.temp0, b1.cnormal0, csolver);
		wallcollisionsolver(b1.temp1, b1.cnormal1, csolver);

		if (b0.is_collide(b1) && b1.is_collide(b0))
		{
				b0.createcollision(b1, csolver);
				b1.createcollision(b0, csolver);
		}

		float sn = 1.0;

		//solver
		for (uint32_t loopcounter = 0; loopcounter < 10; loopcounter++)
		{
				for (auto &x : csolver)
						x.projection(sn);

				b0.projectdistance(sn);
				b1.projectdistance(sn);
		}

		b0.refineposition();

		b1.refineposition();

		b0.resetcn();
		b1.resetcn();
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

		fvec2 p00(-0.2, 0.0);
		fvec2 p01(-0.3, 0.3);
		fvec2 v00(2.7, 1.4);
		fvec2 v01(2.4, 1.3);
		float barcolor0[4] = {1.0, 0.0, 0.0, 0.0};
		bar b0(1.0, 2.0, 0.3, p00, p01, v00, v01, barcolor0);

		fvec2 p10(0.4, 0.3);
		fvec2 p11(0.4, 0.0);
		fvec2 v10(-1.3, 1.4);
		fvec2 v11(-1.4, 2.3);
		float barcolor1[4] = {0.0, 1.0, 0.0, 0.0};
		bar b1(2.0, 3.0, 0.4, p10, p11, v10, v11, barcolor1);

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
						timestep(b0, b1);
						ltime = getTime();
				}

				//rendering

				wl0.draw();
				wl1.draw();
				wl2.draw();
				wl3.draw();

				b0.draw();
				b1.draw();

				//swapbuff
				mywindow.swapbuf();
				//wait event
				mywindow.waitevent(0.01);

				step++;
				cout << "step: " << step;
				cout << "realtime : " << ctime << "  vtime : " << step * dt << endl;
		}

		return 0;
}
