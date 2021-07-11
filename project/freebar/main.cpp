#include <iostream>
#include <cstdint>
#include <cmath>

#include "mathfunc/vec.hpp"
#include "mathfunc/matrix.hpp"

#include "opengl/visualize.hpp"
#include "opengl/window.hpp"

#include "opengl/line.hpp"
#include "opengl/point.hpp"

float dt = 0.005;

float m0 = 1.0;
float w0 = 1.0 / m0;
float m1 = 2.0;
float w1 = 1.0 / m1;

float L = 0.2;

float epsilon = 0.00001;

bool is_collide(const fvec2 &tempposi)
{
		return tempposi.x < -0.8 || tempposi.x > 0.8 || tempposi.y < -0.8 || tempposi.y > 0.8;
}

void projectcollision(fvec2 &tempposi, bool &isx_collide, bool &isy_collide)
{
		if (tempposi.x < -0.8)
		{
				fvec2 normal(1.0, 0.0);

				float C = tempposi.x - (-0.8);
				fvec2 dC = normal;

				fvec2 dposi = dC * (-1.0f * (C / dC.dot(dC)));

				tempposi = tempposi + dposi;

				isx_collide = true;
		}
		if (tempposi.x > 0.8)
		{
				fvec2 normal(-1.0, 0.0);

				float C = 0.8 - tempposi.x;
				fvec2 dC = normal;

				fvec2 dposi = dC * (-1.0 * (C / dC.dot(dC)));

				tempposi = tempposi + dposi;

				isx_collide = true;
		}
		if (tempposi.y < -0.8)
		{
				fvec2 normal(0.0, 1.0);

				float C = tempposi.y - (-0.8);
				fvec2 dC = normal;

				fvec2 dposi = dC * (-1.0 * (C / dC.dot(dC)));

				tempposi = tempposi + dposi;

				isy_collide = true;
		}
		if (tempposi.y > 0.8)
		{
				fvec2 normal(0.0, -1.0);

				float C = 0.8 - tempposi.y;
				fvec2 dC = normal;

				fvec2 dposi = dC * (-1.0 * (C / dC.dot(dC)));

				tempposi = tempposi + dposi;

				isy_collide = true;
		}
}

void projectdistance(fvec2 &tempposi0, fvec2 &tempposi1)
{

		float l = (tempposi0 - tempposi1).length();
		fvec2 dposi0 = (tempposi1 - tempposi0).normalize() * ((w0) / (w0 + w1)) * (l - L);
		fvec2 dposi1 = (tempposi0 - tempposi1).normalize() * ((w1) / (w0 + w1)) * (l - L);

		tempposi0 = tempposi0 + dposi0;
		tempposi1 = tempposi1 + dposi1;
}

void timestep(fvec2 &x0, fvec2 &v0, fvec2 &x1, fvec2 &v1)
{
		//explicit eular
		//gravity acceleration = 9.8 m/ss
		fvec2 velocity0 = v0 + fvec2(0.0, -9.8) * dt;
		fvec2 tempposi0 = x0 + velocity0 * dt;

		fvec2 velocity1 = v1 + fvec2(0.0, -9.8) * dt;
		fvec2 tempposi1 = x1 + velocity1 * dt;

		//solver
		bool isx0_collide = false;
		bool isy0_collide = false;

		bool isx1_collide = false;
		bool isy1_collide = false;

		uint32_t loopcounter = 0;
		while (is_collide(tempposi0) || is_collide(tempposi1) ||
			   std::abs((tempposi0 - tempposi1).length() - L) > epsilon)
		{
				projectcollision(tempposi0, isx0_collide, isy0_collide);
				projectcollision(tempposi1, isx1_collide, isy1_collide);
				projectdistance(tempposi0, tempposi1);

				loopcounter += 1;
				if (loopcounter > 10)
						break;
		}

		v0 = (tempposi0 - x0) / dt;
		v1 = (tempposi1 - x1) / dt;

		//flip velocity if collide happen
		if (isx0_collide)
				v0.x *= -1.0;
		if (isy0_collide)
				v0.y *= -1.0;

		if (isx1_collide)
				v1.x *= -1.0;
		if (isy1_collide)
				v1.y *= -1.0;

		x0 = tempposi0;
		x1 = tempposi1;
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

		float initp0[2] = {0.0, 0.0};
		float initp1[2] = {0.0, 0.2};
		float pcolor[4] = {1.0, 0.0, 0.0, 1.0};
		point1d po0(initp0, pcolor);
		point1d po1(initp1, pcolor);

		fvec2 position0(initp0[0], initp0[1]);
		fvec2 velocity0(0.5, 1.4);

		fvec2 position1(initp1[0], initp1[1]);
		fvec2 velocity1(-8.5, 8.4);

		float vcolor[4] = {0.0, 0.0, 1.0, 1.0};
		line2d visv0(initp0, initp0, vcolor);
		line2d visv1(initp1, initp1, vcolor);

		float lodcolor[4] = {0.0, 1.0, 0.0, 1.0};
		line2d lod(initp0, initp1, lodcolor);

		uint32_t step = 0;
		while (!mywindow.is_shouldclose())
		{
				//clear buf
				mywindow.clear();

				//physics

				timestep(position0, velocity0, position1, velocity1);

				//rendering

				wl0.draw();
				wl1.draw();
				wl2.draw();
				wl3.draw();

				po0.setposition(position0.x, position0.y);
				po0.draw();

				po1.setposition(position1.x, position1.y);
				po1.draw();

				visv0.setposition(position0.x, position0.y, position0.x + velocity0.x, position0.y + velocity0.y);
				visv0.draw();

				visv1.setposition(position1.x, position1.y, position1.x + velocity1.x, position1.y + velocity1.y);
				visv1.draw();

				lod.setposition(position0.x, position0.y, position1.x, position1.y);
				lod.draw();

				//swapbuff
				mywindow.swapbuf();
				//wait event
				mywindow.waitevent(0.01);

				step++;
				cout << "step: " << step << endl;
		}

		return 0;
}
