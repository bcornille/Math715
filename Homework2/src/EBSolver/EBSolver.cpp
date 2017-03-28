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
	// for (auto const& value1 : integrate.b(mesh.getElement(0)))
	// {
	// 	for (auto const& value2 : value1)
	// 	{
	// 		std::cout << value2 << " ";
	// 	}
	// 	std::cout << std::endl;
	// }
	// for (auto const& value : integrate.l(mesh.getElement(0)))
	// {
	// 	std::cout << value << std::endl;
	// }
	Eigen::SparseMatrix<double> delta_sq = mesh.assembleMatrix(integrate);
	std::ofstream outfile("test.json");
	nlohmann::json output;
	output["x"] = mesh.nodes();
	outfile << std::setw(4) << output << std::endl;
	return 0;
}