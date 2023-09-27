#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

struct LinearSystem
{
	virtual void    test()  = 0;
	virtual void    Solve() = 0;
	bool            initialized = false;
	Eigen::MatrixXd X;
};

struct ReferenceSystem : public LinearSystem
{
	int                                                vn = 0;
	Eigen::SparseMatrix<double>                        AT;
	Eigen::SparseMatrix<double>                        ATA;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	Eigen::MatrixXd                                    S;
	Eigen::MatrixXd                                    ATS;

	void test() { std::cout << "ReferenceSystem" << std::endl; }
	void Solve();
};

struct ReducedSystem : public LinearSystem
{
	int                         rank = 0;
	ReferenceSystem*            ref = nullptr;
	Eigen::SparseMatrix<double> dATA;
	Eigen::MatrixXd             S;
	Eigen::MatrixXd             dATS;
	Eigen::MatrixXd             r[3];
	double                      error = 0.;

	void test() { std::cout << "ReducedSystem" << std::endl; }
	void                        Solve();
};

struct DirectSolver : LinearSystem
{
	int                                                vn = 0;
	ReferenceSystem*                                   ref = nullptr;
	Eigen::SparseMatrix<double>                        AT;
	Eigen::SparseMatrix<double>                        ATA;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	Eigen::MatrixXd                                    S;
	Eigen::MatrixXd                                    ATS;
	void                                               Solve();
};
