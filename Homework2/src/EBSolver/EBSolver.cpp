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
	std::ofstream outfile("test.json");
	nlohmann::json output;
	output["x"] = mesh.nodes();
	outfile << std::setw(4) << output << std::endl;
	return 0;
}