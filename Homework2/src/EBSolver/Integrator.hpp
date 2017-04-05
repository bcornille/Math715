#include "Basis1D.hpp"
#include "Hermite.hpp"
#include <array>
#include <cmath>
#include "Element1D.hpp"

#ifndef _Integrator_hpp
#define _Integrator_hpp

class Integrator
{
	public:
		Integrator(nlohmann::json params);
		~Integrator() = default;
		std::array<std::array<double, 4>, 4> b(Element1D e);
		std::array<double, 4> l(Element1D e);
	private:
		GaussLegendre gl;
};

Integrator::Integrator(nlohmann::json params) : gl((int)params["N"]) {}

std::array<std::array<double, 4>, 4> Integrator::b(Element1D e)
{
	std::array<std::array<double, 4>, 4> matrix = {0.0};
	Hermite u;
	for (int i = 0; i < gl.getN(); ++i)
	{
		std::array<double, 4> ud = u.eval_dd(gl.getNode(i)/2.0 + 0.5);
		double w = gl.getWeight(i)/2.0;
		double j = e.jacobian();
		matrix[0][0] += ud[0]*ud[0]*w/pow(j,3);
		matrix[1][0] += ud[0]*ud[1]*w/pow(j,2);
		matrix[2][0] += ud[0]*ud[2]*w/pow(j,3);
		matrix[3][0] += ud[0]*ud[3]*w/pow(j,2);
		matrix[0][1] += ud[1]*ud[0]*w/pow(j,2);
		matrix[1][1] += ud[1]*ud[1]*w/j;
		matrix[2][1] += ud[1]*ud[2]*w/pow(j,2);
		matrix[3][1] += ud[1]*ud[3]*w/j;
		matrix[0][2] += ud[2]*ud[0]*w/pow(j,3);
		matrix[1][2] += ud[2]*ud[1]*w/pow(j,2);
		matrix[2][2] += ud[2]*ud[2]*w/pow(j,3);
		matrix[3][2] += ud[2]*ud[3]*w/pow(j,2);
		matrix[0][3] += ud[3]*ud[0]*w/pow(j,2);
		matrix[1][3] += ud[3]*ud[1]*w/j;
		matrix[2][3] += ud[3]*ud[2]*w/pow(j,2);
		matrix[3][3] += ud[3]*ud[3]*w/j;
	}
	return matrix;
}

std::array<double, 4> Integrator::l(Element1D e)
{
	std::array<double, 4> rhs = {0.0};
	Hermite v;
	for (int i = 0; i < gl.getN(); ++i)
	{
		double w = gl.getWeight(i)/2.0;
		double x_hat = gl.getNode(i)/2.0 + 0.5;
		double j = e.jacobian();
		double x = e.f(x_hat);
		double c = exp(x)*(pow(x,4) + 14.0*pow(x,3) + 49.0*pow(x,2) + 32.0*x - 12.0);
		std::array<double, 4> vd = v.eval(x_hat);
		rhs[0] += c*vd[0]*w*j;
		rhs[1] += c*vd[1]*w*j;
		rhs[2] += c*vd[2]*w*j;
		rhs[3] += c*vd[3]*w*j;
	}
	return rhs;
}

#endif
