#include <iostream>
#include <cstdint>

#include "mathfunc/vec.hpp"
#include "mathfunc/matrix.hpp"

#include "opengl/visualize.hpp"
#include "opengl/window.hpp"

#include "opengl/line.hpp"
#include "opengl/point.hpp"

float dt = 0.005;

void timestep(fvec2 &x, fvec2 &v)
{
		//explicit eular
		//gravity acceleration = 9.8 m/ss
		fvec2 velocity = v + fvec2(0.0, -9.8) * dt;
		fvec2 tempposi = x + velocity * dt;

		//collision solver
		bool isx_collide = false;
		bool isy_collide = false;
		uint32_t loopcounter = 0;
		while (tempposi.x < -0.8 || tempposi.x > 0.8 || tempposi.y < -0.8 || tempposi.y > 0.8)
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

				loopcounter++;
				if (loopcounter > 10)
						break;
		}

		v = (tempposi - x) / dt;

		//flip velocity if collide happen
		if (isx_collide)
				v.x *= -1.0;
		if (isy_collide)
				v.y *= -1.0;

		x = tempposi;
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

		float initp[2] = {0.0, 0.0};
		float pcolor[4] = {1.0, 0.0, 0.0, 1.0};
		point1d po(initp, pcolor);

		fvec2 position(initp[0], initp[1]);
		fvec2 velocity(0.5, 1.4);

		float vcolor[4] = {0.0, 0.0, 1.0, 1.0};
		line2d visv(initp, initp, vcolor);

		uint32_t step = 0;
		while (!mywindow.is_shouldclose())
		{
				//clear buf
				mywindow.clear();

				//physics

				timestep(position, velocity);

				//rendering

				wl0.draw();
				wl1.draw();
				wl2.draw();
				wl3.draw();

				po.setposition(position.x, position.y);
				po.draw();

				visv.setposition(position.x, position.y, position.x + velocity.x, position.y + velocity.y);
				visv.draw();

				//swapbuff
				mywindow.swapbuf();
				//wait event
				mywindow.waitevent(0.01);

				step++;
				cout << "step: " << step << endl;
		}

		return 0;
}
