#include "fastgl.hpp"
#include "Hermite.hpp"
#include <array>

#ifndef _Integrator_hpp
#define _Integrator_hpp

class Integrator
{
	public:
		Integrator(nlohmann::json params);
		~Integrator() = default;
		array<array<double, 4>, 4> b(Hermite u, Hermite v);
		array<double, 4> l(Hermite v);
	private:
		std::vector<fastgl::QuadPair> gl;
};

Integrator::Integrator(nlohmann::json params) : gl((int)params["N"])
{
	for (int i = 0; i < gl.size(); ++i)
	{
		gl[i] = fastgl::GLPair(gl.size(), i + 1);
	}
}

array<array<double, 4>, 4> Integrator::b()
{
	array<array<double, 4>, 4> matrix;
	Hermite u;
	for (int i = 0; i < gl.size(); ++i)
	{
		array<double, 4> ud = u.eval_dd(gl.x());
		matrix[0][0] += ud[0]*ud[0]*gl.weight;
		matrix[1][0] += ud[0]*ud[1]*gl.weight;
		matrix[2][0] += ud[0]*ud[2]*gl.weight;
		matrix[3][0] += ud[0]*ud[3]*gl.weight;
		matrix[0][1] += ud[1]*ud[0]*gl.weight;
		matrix[1][1] += ud[1]*ud[1]*gl.weight;
		matrix[2][1] += ud[1]*ud[2]*gl.weight;
		matrix[3][1] += ud[1]*ud[3]*gl.weight;
		matrix[0][2] += ud[2]*ud[0]*gl.weight;
		matrix[1][2] += ud[2]*ud[1]*gl.weight;
		matrix[2][2] += ud[2]*ud[2]*gl.weight;
		matrix[3][2] += ud[2]*ud[3]*gl.weight;
		matrix[0][3] += ud[3]*ud[0]*gl.weight;
		matrix[1][3] += ud[3]*ud[1]*gl.weight;
		matrix[2][3] += ud[3]*ud[2]*gl.weight;
		matrix[3][3] += ud[3]*ud[3]*gl.weight;
	}
	return matrix;
}

array<double, 4> Integrator::l()
{
	array<double, 4> rhs;
	Hermite v;
	for (int i = 0; i < gl.size(); ++i)
	{
		rhs += v.eval_dd(gl.x())*gl.weight;
	}
}

#endif