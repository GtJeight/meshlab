#include <FittingDocument.h>
#include <common/plugins/interfaces/meshlab_plugin_logger.h>


void ReferenceSystem::Solve()
{
	solver.compute(ATA);
	X = solver.solve(ATS);
	std::cout << "solving system" << std::endl;
}

void ReducedSystem::Solve()
{
	auto                        R     = ref->ATS + dATS;
	auto                        K     = ref->ATA + dATA;
	Eigen::SparseMatrix<double> B     = ref->solver.solve(dATA);
	Eigen::MatrixXd             ri    = ref->solver.solve(R);
	Eigen::MatrixXd             rB[3] = {
        Eigen::MatrixXd::Zero(dATS.rows(), rank),
        Eigen::MatrixXd::Zero(dATS.rows(), rank),
        Eigen::MatrixXd::Zero(dATS.rows(), rank)};

	Eigen::MatrixXd KB[3] = {
		Eigen::MatrixXd::Zero(rank, rank),
		Eigen::MatrixXd::Zero(rank, rank),
		Eigen::MatrixXd::Zero(rank, rank)};

	Eigen::MatrixXd RB[3] = {
		Eigen::VectorXd::Zero(rank), Eigen::VectorXd::Zero(rank), Eigen::VectorXd::Zero(rank)};

	for (int i = 0; i < rank; i++) {
		rB[0](Eigen::placeholders::all, i) = ri(Eigen::placeholders::all, 0);
		rB[1](Eigen::placeholders::all, i) = ri(Eigen::placeholders::all, 1);
		rB[2](Eigen::placeholders::all, i) = ri(Eigen::placeholders::all, 2);
		ri                                 = -B * ri;
	}

	KB[0] = rB[0].transpose() * K * rB[0];
	KB[1] = rB[1].transpose() * K * rB[1];
	KB[2] = rB[2].transpose() * K * rB[2];

	RB[0] = rB[0].transpose() * R(Eigen::placeholders::all, 0);
	RB[1] = rB[1].transpose() * R(Eigen::placeholders::all, 1);
	RB[2] = rB[2].transpose() * R(Eigen::placeholders::all, 2);

	X.resize(dATS.rows(), 3);
	X(Eigen::placeholders::all, 0) = rB[0] * (KB[0].colPivHouseholderQr().solve(RB[0]));
	X(Eigen::placeholders::all, 1) = rB[1] * (KB[1].colPivHouseholderQr().solve(RB[1]));
	X(Eigen::placeholders::all, 2) = rB[2] * (KB[2].colPivHouseholderQr().solve(RB[2]));

	std::cout << "solving system2" << std::endl;
}
