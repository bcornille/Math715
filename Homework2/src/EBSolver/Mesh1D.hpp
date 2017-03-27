#include "json.hpp"
#include <vector>
#include "Element1D.hpp"

#ifndef _Mesh1D_hpp
#define _Mesh1D_hpp

class Mesh1D
{
	public:
		Mesh1D(nlohmann::json params);
		~Mesh1D() = default;
		std::vector<double> nodes();
		Element1D getElement(int i);
	private:
		std::vector<double> x;
};

Mesh1D::Mesh1D(nlohmann::json params) : x((int)params["N"] + 2)
{
	double deltax = ((double)params["x_max"] - (double)params["x_min"])/(x.size() - 1);
	for (int i = 0; i < x.size(); ++i)
	{
		x[i] = (int)params["x_min"] + i*deltax;
	}
}

inline std::vector<double> Mesh1D::nodes()
{
	return x;
}

inline Element1D Mesh1D::getElement(int i)
{
	return Element1D(x[i], x[i+1]);
}

#endif