#include "json.hpp"
#include <vector>
#include "Element1D.hpp"
#include "Integrator.hpp"
#include "Eigen/Sparse"

#ifndef _Mesh1D_hpp
#define _Mesh1D_hpp

class Mesh1D
{
	public:
		Mesh1D(nlohmann::json params);
		~Mesh1D() = default;
		std::vector<double> nodes();
		Element1D getElement(int i);
		Eigen::SparseMatrix<double> assembleMatrix(Integrator integrate);
	private:
		std::vector<double> x;
		int N_el;
		int N_node;
};

Mesh1D::Mesh1D(nlohmann::json params) :
	N_node(params["N"]), N_el(N_node + 1), x(N_el + 1)
{
	double deltax = (((double)params["x_max"] - (double)params["x_min"])
		/(x.size() - 1));
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

Eigen::SparseMatrix<double> Mesh1D::assembleMatrix(Integrator integrate)
{
	Eigen::SparseMatrix<double> matrix(2*N_node,2*N_node);
	matrix.reserve(Eigen::VectorXi::Constant(2*N_node,6));
	matrix.setZero();
	std::array<std::array<double, 4>, 4> minimatrix = integrate.b(getElement(0));
	matrix.coeffRef(0,0) += minimatrix[2][2];
	matrix.coeffRef(1,0) += minimatrix[3][2];
	matrix.coeffRef(0,1) += minimatrix[2][3];
	matrix.coeffRef(1,1) += minimatrix[3][3];
	for (int i = 1; i < N_el - 1; ++i)
	{
		minimatrix = integrate.b(getElement(i));
		matrix.coeffRef(2*i-2,2*i-2) += minimatrix[0][0];
		matrix.coeffRef(2*i-1,2*i-2) += minimatrix[1][0];
		matrix.coeffRef(2*i,2*i-2) += minimatrix[2][0];
		matrix.coeffRef(2*i+1,2*i-2) += minimatrix[3][0];
		matrix.coeffRef(2*i-2,2*i-1) += minimatrix[0][1];
		matrix.coeffRef(2*i-1,2*i-1) += minimatrix[1][1];
		matrix.coeffRef(2*i,2*i-1) += minimatrix[2][1];
		matrix.coeffRef(2*i+1,2*i-1) += minimatrix[3][1];
		matrix.coeffRef(2*i-2,2*i) += minimatrix[0][2];
		matrix.coeffRef(2*i-1,2*i) += minimatrix[1][2];
		matrix.coeffRef(2*i,2*i) += minimatrix[2][2];
		matrix.coeffRef(2*i+1,2*i) += minimatrix[3][2];
		matrix.coeffRef(2*i-2,2*i+1) += minimatrix[0][3];
		matrix.coeffRef(2*i-1,2*i+1) += minimatrix[1][3];
		matrix.coeffRef(2*i,2*i+1) += minimatrix[2][3];
		matrix.coeffRef(2*i+1,2*i+1) += minimatrix[3][3];
	}
	minimatrix = integrate.b(getElement(N_el));
	matrix.coeffRef(2*N_node,2*N_node) += minimatrix[0][0];
	matrix.coeffRef(2*N_node,2*N_node+1) += minimatrix[0][1];
	matrix.coeffRef(2*N_node+1,2*N_node) += minimatrix[1][0];
	matrix.coeffRef(2*N_node+1,2*N_node+1) += minimatrix[1][1];
	matrix.makeCompressed();
	return matrix;
}

#endif