#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

struct InitialProblem
{


	int vn = 0;
	Eigen::SparseMatrix<double> A;
	Eigen::SparseMatrix<double> ATA;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	Eigen::MatrixXd S;
	Eigen::MatrixXd ATS;
	Eigen::MatrixXd X;
};

struct ReduceModel
{
	InitialProblem* ref = nullptr;
	Eigen::SparseMatrix<double> dATA;
	Eigen::MatrixXd S;
	Eigen::MatrixXd dATS;
	Eigen::MatrixXd r[3];
	Eigen::MatrixXd X;
	double error = 0.;
};
