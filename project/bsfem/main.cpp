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

constexpr uint32_t N = 12;

constexpr float mu = 10.00;
constexpr float lambda = 0.500;
constexpr float m = 0.01;

fvec2 u[3 * (N + 1)];
fvec2 tempu[3 * (N + 1)];
fvec2 refposition[3 * (N + 1)];
fvec2 v[3 * (N + 1)];

uint32_t ilist[7 * N + 3];

void projection()
{

		float lambdalist[4 * N];
		for (uint32_t x = 0; x < 4 * N; x++)
				lambdalist[x] = 0.0;

		for (uint32_t x = 0; x < 3; x++)
		{
				//std::cout << "-------------------" << std::endl;
				for (uint32_t i = 0; i < 4 * N; i++)
				{
						//for each triangle
						uint32_t index0;
						uint32_t index1;
						uint32_t index2;
						if (i % 4 == 0)
						{
								index0 = (i / 4) * 3;
								index1 = (i / 4) * 3 + 1;
								index2 = (i / 4) * 3 + 4;
						}
						if (i % 4 == 1)
						{
								index0 = (i / 4) * 3;
								index1 = (i / 4) * 3 + 4;
								index2 = (i / 4) * 3 + 3;
						}
						if (i % 4 == 2)
						{
								index0 = (i / 4) * 3 + 1;
								index1 = (i / 4) * 3 + 2;
								index2 = (i / 4) * 3 + 5;
						}
						if (i % 4 == 3)
						{
								index0 = (i / 4) * 3 + 1;
								index1 = (i / 4) * 3 + 5;
								index2 = (i / 4) * 3 + 4;
						}

						fvec2 X0 = refposition[index0];
						fvec2 X1 = refposition[index1];
						fvec2 X2 = refposition[index2];

						fvec2 u0 = tempu[index0];
						fvec2 u1 = tempu[index1];
						fvec2 u2 = tempu[index2];

						fvec2 x0 = u0 + X0;
						fvec2 x1 = u1 + X1;
						fvec2 x2 = u2 + X2;

						fmat2 F = mat2(x1 - x0, x2 - x0) * (mat2(X1 - X0, X2 - X0).inverse());
						fmat2 E = F.transpose() * F - fmat2::identity();
						fmat2 T = 2.0 * mu * F * E * F.transpose() + lambda * E.trace() * F * F.transpose();

						float vol = mat2(x1 - x0, x2 - x0).det();
						float W = 2.0 * vol * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
						float C = std::sqrt(2.0 * vol * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace()));

						//std::cout << "F= " << F << std::endl;
						//std::cout << "E= " << E << std::endl;
						//std::cout << "T= " << T << std::endl;
						//std::cout << "W= " << W << std::endl;
						//std::cout << "C= " << C << std::endl;
						//std::cout << std::endl;

						if (C > epsilon)
						{
								fvec2 dC0 = -(2.0 / C) * T * fvec2(x2.y - x1.y, x1.x - x2.x);
								fvec2 dC1 = -(2.0 / C) * T * fvec2(x0.y - x2.y, x2.x - x0.x);
								fvec2 dC2 = -(2.0 / C) * T * fvec2(x1.y - x0.y, x0.x - x1.x);

								//std::cout << "dC0= " << dC0 << std::endl;
								//std::cout << "dC1= " << dC1 << std::endl;
								//std::cout << "dC2= " << dC2 << std::endl;

								float dtdlambda = (-C - lambdalist[i]) /
												  ((dC0.dot(dC0) + dC1.dot(dC1) + dC2.dot(dC2)) / m + 1.0 / (dt * dt));

								tempu[index0] = tempu[index0] + dtdlambda * (1.0 / m) * dC0;
								tempu[index1] = tempu[index1] + dtdlambda * (1.0 / m) * dC1;
								tempu[index2] = tempu[index2] + dtdlambda * (1.0 / m) * dC2;

								lambdalist[i] += dtdlambda / (dt * dt);
						}

						if (i < 4)
						{
								tempu[0] = fvec2(0.0);
								tempu[1] = fvec2(0.0);
								tempu[2] = fvec2(0.0);
						}
				}
		}
}

int main(int argc, char const *argv[])
{

		for (uint32_t i = 0; i <= N; i++)
		{
				double hx = 1.2 / N;
				u[3 * i + 0].x = 0.0;
				u[3 * i + 0].y = 0.0;
				refposition[3 * i + 0].x = hx * i - 0.8;
				refposition[3 * i + 0].y = 0.05;

				u[3 * i + 1].x = 0.0;
				u[3 * i + 1].y = 0.0;
				refposition[3 * i + 1].x = hx * i - 0.8;
				refposition[3 * i + 1].y = 0.0;

				u[3 * i + 2].x = 0.0;
				u[3 * i + 2].y = 0.0;
				refposition[3 * i + 2].x = hx * i - 0.8;
				refposition[3 * i + 2].y = -0.05;
		}

		for (uint32_t i = 0; i < N; i++)
		{
				ilist[7 * i + 0] = 3 * i + 0;
				ilist[7 * i + 1] = 3 * i + 1;
				ilist[7 * i + 2] = 3 * i + 2;
				ilist[7 * i + 3] = 3 * i + 5;
				ilist[7 * i + 4] = 3 * i + 1;
				ilist[7 * i + 5] = 3 * i + 4;
				ilist[7 * i + 6] = 3 * i + 0;
		}
		ilist[7 * N + 0] = 3 * N + 0;
		ilist[7 * N + 1] = 3 * N + 1;
		ilist[7 * N + 2] = 3 * N + 2;

		//
		using namespace std;

		Window mywindow = visualizeinit();

		shader mys = shader("../../../opengl/shadercode/vert.c", "../../../opengl/shadercode/frag.c");
		mys.useprogram();

		lineset rendererdline(nullptr, 3 * (N + 1), ilist, 7 * N + 3);
		rendererdline.setcolor(0.2, 0.2, 0.0, 1.0);

		uint32_t testilist[5];
		testilist[0] = 0;
		testilist[1] = 1;
		testilist[2] = 2;
		testilist[3] = 3;
		testilist[4] = 0;

		lineset test(nullptr, 4, testilist, 5);
		test.setposition(0, -0.8, 0.8, 0.0);
		test.setposition(1, -0.8, -0.8, 0.0);
		test.setposition(2, 0.8, -0.8, 0.0);
		test.setposition(3, 0.8, 0.8, 0.0);

		for (uint32_t i = 0; i < 3 * (N + 1); i++)
		{
				v[i] = fvec2(0.0, 0.0);
		}

		while (!mywindow.is_shouldclose())
		{
				//clear buf
				mywindow.clear();

				for (uint32_t i = 0; i < 3 * (N + 1); i++)
				{
						tempu[i] = u[i] + dt * v[i] + dt * dt * fvec2(0.0, -0.20);
				}

				//
				projection();

				for (uint32_t i = 0; i < 3 * (N + 1); i++)
				{
						v[i] = (tempu[i] - u[i]) / dt;
						u[i] = tempu[i];
				}

				for (uint32_t i = 0; i < 3 * (N + 1); i++)
				{
						rendererdline.setposition(i, u[i].x + refposition[i].x, u[i].y + refposition[i].y, 0.0);
				}

				rendererdline.draw();

				test.draw();

				//mywindow.waitevent(100);

				//swapbuff
				mywindow.swapbuf();
				//wait event
				mywindow.clearstatus();
				mywindow.pollevent();
		}

		return 0;
}
