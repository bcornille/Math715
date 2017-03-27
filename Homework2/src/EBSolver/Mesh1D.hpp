#include "json.hpp"
#include <vector>
#include "Integrator.hpp"

#ifndef _Mesh1D_hpp
#define _Mesh1D_hpp

class Mesh1D
{
	public:
		Mesh1D(nlohmann::json params);
		~Mesh1D() = default;
		std::vector<double> nodes();
	private:
		std::vector<double> x;
		Integrator integral;
};

Mesh1D::Mesh1D(nlohmann::json params) : x((int)params["N"] + 2)
{
	double deltax = ((double)params["x_max"] - (double)params["x_min"])/(x.size() - 1);
	for (int i = 0; i < x.size(); ++i)
	{
		x[i] = (int)params["x_min"] + i*deltax;
	}
}

std::vector<double> Mesh1D::nodes()
{
	return x;
}

#endif