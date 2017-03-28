#include "json.hpp"
#include <iostream>
#include <fstream>
#include "Mesh1D.hpp"

int main(int argc, char const *argv[])
{
	// std::cout << argc << std::endl;
	if (argc != 2)
	{
		std::cerr << "Must provide one input file." << std::endl;
		return 1;
	}
	nlohmann::json input;
	std::ifstream infile(argv[1]);
	infile >> input;
	Mesh1D mesh(input["Mesh"]);
	Integrator integrate(input["Integrator"]);
	// for (auto const& value1 : integrate.b(mesh.getElement(1)))
	// {
	// 	for (auto const& value2 : value1)
	// 	{
	// 		std::cout << value2 << " ";
	// 	}
	// 	std::cout << std::endl;
	// }
	// for (auto const& value : integrate.l(mesh.getElement(1)))
	// {
	// 	std::cout << value << std::endl;
	// }
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
	Eigen::SparseMatrix<double> delta_sq = mesh.assembleMatrix(integrate);
	// std::cout << delta_sq << std::endl;
	solver.analyzePattern(delta_sq);
	solver.factorize(delta_sq);
	Eigen::VectorXd b = mesh.assembleRHS(integrate);
	// std::cout << b << std::endl;
	Eigen::VectorXd x =solver.solve(b);
	// std::cout << x << std::endl;
	std::vector<double> u((int)input["Mesh"]["N"]+2);
	std::vector<double> up((int)input["Mesh"]["N"]+2);
	u[0] = 0.0;
	up[0] = 0.0;
	for (int i = 1; i < u.size() - 1; ++i)
	{
		u[i] = x[2*i-2];
		up[i] = x[2*i-1];
	}
	u[u.size()-1] = 0.0;
	up[up.size()-1] = 0.0;
	std::ofstream outfile("test.json");
	nlohmann::json output;
	output["x"] = mesh.nodes();
	output["u"] = u;
	output["up"] = up;
	outfile << std::setw(4) << output << std::endl;
	return 0;
}