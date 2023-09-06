/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005-2021                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include "filter_subdivfitting.h"

using namespace vcg;

/**
 * @brief
 * Constructor usually performs only two simple tasks of filling the two lists
 *  - typeList: with all the possible id of the filtering actions
 *  - actionList with the corresponding actions.
 * If you want to add icons to your filtering actions you can do here by construction the QActions accordingly
 */
FilterSubdivFittingPlugin::FilterSubdivFittingPlugin()
{ 
	typeList = { FP_SUBDIV_FITTING};

	for(const ActionIDType& tt : typeList)
		actionList.push_back(new QAction(filterName(tt), this));
}

QString FilterSubdivFittingPlugin::pluginName() const
{
	return "FilterSubdivFitting";
}

QString FilterSubdivFittingPlugin::vendor() const
{
	return "CNR-ISTI-VCLab";
}

/**
 * @brief ST() must return the very short string describing each filtering action
 * (this string is used also to define the menu entry)
 * @param filterId: the id of the filter
 * @return the name of the filter
 */
QString FilterSubdivFittingPlugin::filterName(ActionIDType filterId) const
{
	switch(filterId) {
	case FP_SUBDIV_FITTING :
		return "Subdivision Surface Fitting";
	default :
		assert(0);
		return "";
	}
}


/**
 * @brief // Info() must return the longer string describing each filtering action
 * (this string is used in the About plugin dialog)
 * @param filterId: the id of the filter
 * @return an info string of the filter
 */
QString FilterSubdivFittingPlugin::filterInfo(ActionIDType filterId) const
{
	switch(filterId) {
	case FP_SUBDIV_FITTING :
		return "Fitting samples by a subdivision surface";
	default :
		assert(0);
		return "Unknown Filter";
	}
}

/**
 * @brief The FilterClass describes in which generic class of filters it fits.
 * This choice affect the submenu in which each filter will be placed
 * More than a single class can be chosen.
 * @param a: the action of the filter
 * @return the class od the filter
 */
FilterSubdivFittingPlugin::FilterClass FilterSubdivFittingPlugin::getClass(const QAction *a) const
{
	switch(ID(a)) {
	case FP_SUBDIV_FITTING :
		return FilterPlugin::Other;
	default :
		assert(0);
		return FilterPlugin::Generic;
	}
}

/**
 * @brief FilterSamplePlugin::filterArity
 * @return
 */
FilterPlugin::FilterArity FilterSubdivFittingPlugin::filterArity(const QAction*) const
{
	return SINGLE_MESH;
}

/**
 * @brief FilterSamplePlugin::getPreConditions
 * @return
 */
int FilterSubdivFittingPlugin::getPreConditions(const QAction*) const
{
	return MeshModel::MM_NONE;
}

/**
 * @brief FilterSamplePlugin::postCondition
 * @return
 */
int FilterSubdivFittingPlugin::postCondition(const QAction*) const
{
	return MeshModel::MM_VERTCOORD | MeshModel::MM_FACENORMAL | MeshModel::MM_VERTNORMAL;
}

/**
 * @brief This function returns a list of parameters needed by each filter.
 * For each parameter you need to define,
 * - the name of the parameter,
 * - the default value
 * - the string shown in the dialog
 * - a possibly long string describing the meaning of that parameter (shown as a popup help in the dialog)
 * @param action
 * @param m
 */
RichParameterList
FilterSubdivFittingPlugin::initParameterList(const QAction* action, const MeshDocument& md)
{
	RichParameterList parlst;
	switch(ID(action)) {
	case FP_SUBDIV_FITTING: {
		const MeshModel* target = md.mm();
		// looking for a second mesh different that the current one
		for (const MeshModel& t : md.meshIterator()) {
			if (&t != md.mm()) {
				target = &t;	// control mesh topo
				break;
			}
		}
		parlst.addParam(RichMesh(
			"samples",
			md.mm()->id(),
			&md,
			"Samples",
			"Samples to be fitted"));
		parlst.addParam(RichMesh(
			"control_mesh",
			target->id(),
			&md,
			"Control mesh",
			"Control mesh for subdivision"));
		}
		break;
	default :
		assert(0);
	}
	return parlst;
}

/**
 * @brief The Real Core Function doing the actual mesh processing.
 * @param action
 * @param md: an object containing all the meshes and rasters of MeshLab
 * @param par: the set of parameters of each filter, with the values set by the user
 * @param cb: callback object to tell MeshLab the percentage of execution of the filter
 * @return true if the filter has been applied correctly, false otherwise
 */
std::map<std::string, QVariant> FilterSubdivFittingPlugin::applyFilter(
		const QAction * action,
		const RichParameterList & par,
		MeshDocument &md,
		unsigned int& /*postConditionMask*/,
		vcg::CallBackPos *cb)
{
	switch(ID(action)) {
	case FP_SUBDIV_FITTING: {
		MeshModel* curMM = md.mm();
		CMeshO&    m     = curMM->cm;
		//for (size_t i = 0; i < m.fn; i++)
		//	log("%d: %f, %f, %f\n", i, m.face[i].N()[0], m.face[i].N()[1], m.face[i].N()[2]);
		solveFootPoints(
			md,
			*md.getMesh(par.getMeshId("samples")),
			*md.getMesh(par.getMeshId("control_mesh")),
			FootPointMode::MODE_SUBDIVISION);

		{
			//Eigen::MatrixXi P = matrixPickUP(7, 1);
			//std::cout << P;
			//int i = -1, N = 3;
			//std::cout << ((i % N)+N)%N;

			auto test = matrixPatchSubdiv(3, 2, true);
		}
	} break;
	default :
		wrongActionCalled(action);
	}
	return std::map<std::string, QVariant>();
}

void FilterSubdivFittingPlugin::solveFootPoints(
	MeshDocument&    md,
	MeshModel&       spl,
	const MeshModel& ctrlm,
	int              mode)
{
	CMeshO::PerVertexAttributeHandle<const CFaceO*> ftptrs;
	if (tri::HasPerVertexAttribute(spl.cm, "FootTriangle")) {
		ftptrs =
			tri::Allocator<CMeshO>::FindPerVertexAttribute<const CFaceO*>(spl.cm, "FootTriangle");
		if (!tri::Allocator<CMeshO>::IsValidHandle<const CFaceO*>(spl.cm, ftptrs)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		ftptrs =
			tri::Allocator<CMeshO>::AddPerVertexAttribute<const CFaceO*>(spl.cm, "FootTriangle");

	for (auto vi = spl.cm.vert.begin(); vi != spl.cm.vert.end(); vi++) {
		if (!(*vi).IsD()) {
			ftptrs[vi] = nullptr;
		}
	}


	CMeshO::PerVertexAttributeHandle<Point3f> ftbarycoords;
	if (tri::HasPerVertexAttribute(spl.cm, "BaryCoord")) {
		ftbarycoords = tri::Allocator<CMeshO>::FindPerVertexAttribute<Point3f>(spl.cm, "BaryCoord");
		if (!tri::Allocator<CMeshO>::IsValidHandle<Point3f>(spl.cm, ftbarycoords)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		ftbarycoords = tri::Allocator<CMeshO>::AddPerVertexAttribute<Point3f>(spl.cm, "BaryCoord");

	for (auto vi = spl.cm.vert.begin(); vi != spl.cm.vert.end(); vi++) {
		if (!(*vi).IsD()) {
			ftbarycoords[vi] = Point3f(0, 0, 0);
		}
	}


	switch (mode) {
	case FootPointMode::MODE_MESH: {
		// solve foot point by brute force
		for (auto si = spl.cm.vert.begin(); si != spl.cm.vert.end(); si++) {
			std::pair<float, vcg::Point3f> minpair(
				std::numeric_limits<float>::max(), Point3f(0.f, 0.f, 0.f));
			const CFaceO*                  footface = nullptr;
			for (auto fi = ctrlm.cm.face.begin(); fi != ctrlm.cm.face.end(); fi++) {
				auto d = distancePointTriangle(*si, *fi);
				if (d.first < minpair.first) {
					minpair = d;
					footface = &(*fi);
				}
			}
		}
	} break;
	case FootPointMode::MODE_SUBDIVISION: {
		// solve foot point by Newton method, need subdivision surface evaluation
		for (auto si = spl.cm.vert.begin(); si != spl.cm.vert.end(); si++) {
			std::pair<float, vcg::Point3f> minpair(
				std::numeric_limits<float>::max(), Point3f(0.f, 0.f, 0.f));
			const CFaceO* footface = nullptr;
			for (auto fi = ctrlm.cm.face.begin(); fi != ctrlm.cm.face.end(); fi++) {
				auto d = distancePointTriangle(*si, *fi);
				if (d.first < minpair.first) {
					minpair  = d;
					footface = &(*fi);
				}
			}

			ftptrs[si] = footface;
			ftbarycoords[si] = minpair.second;
		}

		for (auto si = spl.cm.vert.begin(); si != spl.cm.vert.end(); si++) {
			auto p = evaluateLimitPoint(ftptrs[si], ftbarycoords[si]);
		}
	} break;
	default: break;
	}

	// test code for vertex attribute
	//auto test1 = tri::Allocator<CMeshO>::FindPerVertexAttribute<const CFaceO*>(spl.cm, "FootTriangle");
	//if (!tri::Allocator<CMeshO>::IsValidHandle<const CFaceO*>(spl.cm, test1)) {
	//	throw MLException("attribute test failed");
	//}
	//for (auto vi = spl.cm.vert.begin(); vi != spl.cm.vert.end(); vi++) {
	//	if (!(*vi).IsD()) {
	//		log("%d\n", test1[vi] == nullptr);
	//	}
	//}

	//auto test2 = tri::Allocator<CMeshO>::FindPerVertexAttribute<Point3f>(spl.cm, "BaryCoord");
	//if (!tri::Allocator<CMeshO>::IsValidHandle<Point3f>(spl.cm, test2)) {
	//	throw MLException("attribute test failed");
	//}
	//for (auto vi = spl.cm.vert.begin(); vi != spl.cm.vert.end(); vi++) {
	//	if (!(*vi).IsD()) {
	//		log("%f, %f, %f\n", test2[vi][0], test2[vi][1], test2[vi][2]);
	//	}
	//}
}

std::pair<float,Point3f> FilterSubdivFittingPlugin::distancePointTriangle(const CVertexO& _p,const CFaceO& _f)
{
	auto n  = _f.N();
	auto p  = _p.P();
	auto v0 = _f.V(0)->P();
	auto v1 = _f.V(1)->P();
	auto v2 = _f.V(2)->P();
	float height                = (p - v0).dot(n);
	float squared_parallel_dist = 0.f, l0 = 0.f, l1 = 0.f, l2 = 0.f;
	auto  project = p - height * n;

	auto det = [](const Point3f& v0, const Point3f& v1, const Point3f& v2) -> float {
		return v0[0] * v1[1] * v2[2] + v0[1] * v1[2] * v2[0] + v0[2] * v1[0] * v2[1] -
			   v0[2] * v1[1] * v2[0] - v0[1] * v1[0] * v2[2] - v0[0] * v1[2] * v2[1];
	};

	// case corner
	if ((project - v0).dot(v0 - v1) > 0.f && (project - v0).dot(v0 - v2) > 0.f) {
		squared_parallel_dist = (project - v0).SquaredNorm();
		l0                    = 1.f;
		l1                    = 0.f;
		l2                    = 0.f;
	}
	else if ((project - v1).dot(v1 - v2) > 0.f && (project - v1).dot(v1 - v0) > 0.f) {
		squared_parallel_dist = (project - v1).SquaredNorm();
		l0                    = 0.f;
		l1                    = 1.f;
		l2                    = 0.f;
	}
	else if ((project - v2).dot(v2 - v0) > 0.f && (project - v2).dot(v2 - v1) > 0.f) {
		squared_parallel_dist = (project - v2).SquaredNorm();
		l0                    = 0.f;
		l1                    = 0.f;
		l2                    = 1.f;
	}
	else {
		float detf = det(n, v0-v1,v1-v2);
		l0         = det(n, project - v1, v1 - v2) / detf;
		l1         = det(n, v0 - project, project - v2) / detf;
		l2         = det(n, v0 - v1, v1 - project) / detf;

		if (l0 < 0.f) {
			l0 = 0.f;
			l1 = (project - v2).dot((v1 - v2).normalized()) / (v1 - v2).Norm();
			l2 = 1.f - l1;
			squared_parallel_dist = (project - v2).SquaredNorm() - l1 * l1;
		}
		else if (l1 < 0.f) {
			l1 = 0.f;
			l2 = (project - v0).dot((v2 - v0).normalized()) / (v2 - v0).Norm();
			l0 = 1.f - l2;
			squared_parallel_dist = (project - v0).SquaredNorm() - l2 * l2;
		}
		else if (l2 < 0.f) {
			l2 = 0.f;
			l0 = (project - v1).dot((v0 - v1).normalized()) / (v0 - v1).Norm();
			l1 = 1.f - l0;
			squared_parallel_dist = (project - v1).SquaredNorm() - l0 * l0;
		}
		else
			squared_parallel_dist = 0.f;
	}

	return std::pair<float, Point3f>(
		sqrtf(height * height + squared_parallel_dist), Point3f(l0, l1, l2));
}

Point3f
FilterSubdivFittingPlugin::evaluateLimitPoint(const CFaceO* ft, const vcg::Point3f& barycoord)
{
	auto b = weightsPatch(ft, barycoord[1], barycoord[2]);
	return Point3f(0.f, 0.f, 0.f);
}

Eigen::VectorXd FilterSubdivFittingPlugin::weightsPatch(const CFaceO* ft, float v, float w)
{
	//Eigen::VectorXf b(N + 6);
	//if (N == 6) {
	//	
	//}

	auto b = weightsIrregularPatch(7, v, w);


	return Eigen::VectorXd::Zero(1);
}

Eigen::VectorXd FilterSubdivFittingPlugin::weightsIrregularPatch(int N, float v, float w)
{
	if (N == 6)
		return weightsRegularPatch(1.f - v - w, v, w);

	if (v + w < eps) {
	}
	else {
		int k    = -1;
		int n    = 1 - log2f(v + w);
		int pow2 = powf(2.f, n - 1);
		v *= pow2;
		w *= pow2;

		// barycoord transform to regular patch Omega_n,k
		if (v > 0.5f) {
			k = 0;
			v = 2 * v - 1;
			w = 2 * w;
		}
		else if (w > 0.5f) {
			k = 2;
			v = 2 * v;
			w = 2 * w - 1;
		}
		else {
			k = 1;
			v = 1 - 2 * v;
			w = 1 - 2 * w;
		}
		float u = 1.f - v - w;

		auto B = weightsRegularPatch(u, v, w) * matrixPickUP(N, k) * matrixPatchSubdiv(N, n);

		//Eigen::RowVectorXd b  = weightsRegularPatch(u, v, w);
		//std::cout << "111111111111111" << std::endl;
		//Eigen::MatrixXi    P  = matrixPickUP(N, k);
		//std::cout << "222222222222222" << std::endl;
		//Eigen::MatrixXd    A_ = matrixPatchSubdiv(N);
		//std::cout << "333333333333333" << std::endl;

	}
	

	return Eigen::VectorXd::Zero(1);
}

Eigen::RowVectorXd FilterSubdivFittingPlugin::weightsRegularPatch(float u, float v, float w)
{
	Eigen::RowVectorXd b(12);
	b << pow(u, 4) + 2 * pow(u, 3) * v, pow(u, 4) + 2 * pow(u, 3) * w,
		pow(u, 4) + 2 * pow(u, 3) * w + 6 * pow(u, 3) * v + 6 * u * u * v * w + 12 * u * u * v * v +
			6 * u * v * v * w + 6 * u * pow(v, 3) + 2 * pow(v, 3) * w + pow(v, 4),
		6 * pow(u, 4) + 24 * pow(u, 3) * w + 24 * u * u * w * w + 8 * u * pow(w, 3) + pow(w, 4) +
			24 * pow(u, 3) * v + 60 * u * u * v * w + 36 * u * v * w * w + 6 * v * pow(w, 3) +
			24 * u * u * v * v + 36 * u * v * v * w + 12 * v * v * w * w + 8 * u * pow(v, 3) +
			6 * pow(v, 3) * w + pow(v, 4),
		pow(u, 4) + 6 * pow(u, 3) * w + 12 * u * u * w * w + 6 * u * pow(w, 3) + pow(w, 4) +
			2 * pow(u, 3) * v + 6 * u * u * v * w + 6 * u * v * w * w + 2 * v * pow(w, 3),
		2 * u * pow(v, 3) + pow(v, 4),
		pow(u, 4) + 6 * pow(u, 3) * w + 12 * u * u * w * w + 6 * u * pow(w, 3) + pow(w, 4) +
			8 * pow(u, 3) * v + 36 * u * u * v * w + 36 * u * v * w * w + 8 * v * pow(w, 3) +
			24 * u * u * v * v + 60 * u * v * v * w + 24 * v * v * w * w + 24 * u * pow(v, 3) +
			24 * pow(v, 3) * w + 6 * pow(v, 4),
		pow(u, 4) + 8 * pow(u, 3) * w + 24 * u * u * w * w + 24 * u * pow(w, 3) + 6 * pow(w, 4) +
			6 * pow(u, 3) * v + 36 * u * u * v * w + 60 * u * v * w * w + 24 * v * pow(w, 3) +
			12 * u * u * v * v + 36 * u * v * v * w + 24 * v * v * w * w + 6 * u * pow(v, 3) +
			8 * pow(v, 3) * w + pow(v, 4),
		2 * u * pow(w, 3) + pow(w, 4), 2 * pow(v, 3) * w + pow(v, 4),
		2 * u * pow(w, 3) + pow(w, 4) + 6 * u * v * w * w + 6 * v * pow(w, 3) + 6 * u * v * v * w +
			12 * v * v * w * w + 2 * u * pow(v, 3) + 6 * pow(v, 3) * w + pow(v, 4),
		pow(w, 4) + 2 * v * pow(w, 3);
	return b;
}

Eigen::MatrixXd FilterSubdivFittingPlugin::matrixPickUP(int N, int k)
{
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(12, N + 12);
	switch (k) {

	case 0: {
		std::vector<int> idx {2, 0, N + 3, 1, N, N + 8, N + 2, N + 1, N + 4, N + 7, N + 6, N + 9};
		P(Eigen::placeholders::all, idx) = Eigen::MatrixXd::Identity(12, 12);
	} break;

	case 1: {
		std::vector<int> idx {N + 9, N + 6, N + 4, N + 1, N + 2, N + 5, N, 1, N + 3, N - 1, 0, 2};
		P(Eigen::placeholders::all, idx) = Eigen::MatrixXd::Identity(12, 12);
	} break;

	case 2: {
		std::vector<int> idx {
			0, N - 1, 1, N, N + 5, N + 2, N + 1, N + 4, N + 11, N + 6, N + 9, N + 10};
		P(Eigen::placeholders::all, idx) = Eigen::MatrixXd::Identity(12, 12);
	} break;

	default: break;
	}
	return P;
}

static int id(int i, int N)
{
	return ((i % N) + N) % N;
}

Eigen::MatrixXd FilterSubdivFittingPlugin::matrixPatchSubdiv(int N, int n, bool test)
{
	assert(N >= 3 && n >= 1);
	double          alphaN = 5. / 8. - pow((3 + 2 * cos(2 * M_PI / N)), 2) / 64.;
	double          aN     = 1. - alphaN;
	double          bN     = alphaN / N;
	double          c      = 0.375;
	double          d      = 0.125;
	Eigen::MatrixXd A_ = Eigen::MatrixXd::Zero(N + 12, N + 6);

	// constructing S
	// 1-ring neighbor of EV
	A_(0, 0) = aN;
	A_(0,Eigen::seq(1,N)).array() = bN;
	Eigen::RowVector4d cd(c, c, d, d);
	for (int i = 1; i < N + 1; i++) {
		A_(i, {0, i, id(i - 1, N + 1), id(i + 1, N + 1)}) = cd;
	}

	// constructing S_11,S_12,S_21 and S_22
	// edge vertices
	A_(N + 1, {1, N, 0, N + 1}) = cd;
	A_(N + 3, {1, 2, 0, N + 3}) = cd;
	A_(N + 5, {N - 1, N, 0, N + 5}) = cd;
	A_(N + 6, {1, N + 1, N, N + 2}) = cd;
	A_(N + 7, {1, N + 2, N + 1, N + 3}) = cd;
	A_(N + 8, {1, N + 3, 2, N + 2})     = cd;
	A_(N + 9, {N, N + 1, 1, N + 4})     = cd;
	A_(N + 10, {N, N + 4, N + 1, N + 5}) = cd;
	A_(N + 11, {N, N + 5, N - 1, N + 4}) = cd;
	// vertex vertices
	A_(N + 2, {1, 0, 2, N + 3, N + 2, N + 1, N}) =
		(Eigen::ArrayXd(7) << 0.625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625).finished();
	A_(N + 4, {N, 0, 1, N + 1, N + 4, N + 5, N - 1}) =
		(Eigen::ArrayXd(7) << 0.625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625).finished();

	if (n == 1)
		return A_;

	// constructing eigen structure of A and submatrix S
	Eigen::MatrixXd LambdaPow = Eigen::MatrixXd::Zero(N + 6, N + 6);
	Eigen::MatrixXd V = Eigen::MatrixXd::Zero(N + 6, N + 6);

	if (N == 3) {
		LambdaPow.diagonal() << 1, pow(0.25, n - 1), pow(0.25, n - 1), pow(0.125, n - 1),
			pow(0.125, n - 1), pow(0.125, n - 1), pow(0.0625, n - 1), pow(0.0625, n - 1),
			pow(0.0625, n - 1);
		LambdaPow(N + 4, N + 5) = (n - 1) / pow(16., n - 2);
		V << 1, 0, 0, 0, 0, 0, 0, 0, 33, 1, 0, 1, 0, 0, 0, 0, 0, -22, 1, -1, -1, 0, 0, 0, 0, 0, -22,
			1, 1, 0, 0, 0, 0, 0, 0, -22, 1, 3, 3, 1, -1, 0, 0, 0, 198, 1, 0, 4, 1, 0, 0, 0,
			165. / 16., 473, 1, -3, 0, 0, 1, 0, 0, 0, 198, 1, 4, 0, 0, 0, 1, 1, 165. / 16., 438, 1, 0,
			-3, - 1, 1, 1, 0, 0, 198;
		Eigen::MatrixXd V_inv = Eigen::MatrixXd::Zero(N + 6, N + 6);
		V_inv << 0.4, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0, -1. / 3., -1. / 3., 2. / 3., 0, 0, 0, 0, 0, 0,
			2. / 3., -1. / 3., -1. / 3., 0, 0, 0, 0, 0, -8, 0, 3, 3, 1, 0, 1, 0, 0, -4, 0, 0, 3, 0, 0, 1,
			0, 0, -8, 3, 3, 0, 1, 0, 0, 0, 1, 7. / 11., 26. / 33., -7. / 33., -40. / 33., 0, -1, 1, 1, -1,
			-16. / 165., 0, 16. / 165., 16. / 165., -16. / 165., 16. / 165., -16. / 165., 0, 0, 1. / 55.,
			-1. / 165., -1. / 165., -1. / 165., 0, 0, 0, 0, 0;

		if (test) {
			std::cout << A_(Eigen::seq(0, N + 5), Eigen::placeholders::all) << std::endl;
			std::cout << V * LambdaPow * V_inv << std::endl;
			std::cout << V * V_inv << std::endl;
		}

		return A_ * V * LambdaPow * V_inv;
	}
	else {
		auto f = [&N](int k) -> double { return 0.375 + 0.25 * cos(2 * M_PI * k / N); };

		// constructing Sigma and U0
		LambdaPow(0, 0)  = 1;
		LambdaPow(1, 1)  = pow(0.625 - alphaN, n - 1);
		int  halfN       = (N % 2 == 1) ? ((N + 3) / 2) : (N / 2 + 1);
		int  eigid      = 2;
		auto idseq       = Eigen::ArrayXd(N);
		for (int i = 0; i < N; i++)
			idseq(i) = i;
		for (int i = 3; i <= halfN; i++) {
			LambdaPow(eigid, eigid)         = pow(f(i - 2), n - 1);
			LambdaPow(eigid + 1, eigid + 1) = pow(f(i - 2), n - 1);
			V(Eigen::seq(1, N), eigid).array()     = (2 * M_PI * (i - 2) / N * idseq).cos();
			V(Eigen::seq(1, N), eigid + 1).array() = (2 * M_PI * (i - 2) / N * idseq).sin();
			eigid += 2;
		}
		if (eigid == N) {
			LambdaPow(eigid + 1, eigid + 1)    = pow(0.125, n - 1);
			V(Eigen::seq(1, N), eigid + 1).array() = (2 * M_PI * (halfN - 1) / N * idseq).cos();
		}
		V(Eigen::placeholders::all, {0, 1}).array() = 1;
		V(0,1) = -8/3*alphaN;

		// constructing U1
		auto S_11U0 =
			A_(Eigen::seq(N + 1, N + 5), Eigen::seq(0, N)) * V.topLeftCorner(N + 1, N + 1);
		auto S_12 = A_(Eigen::seq(N + 1, N + 5), Eigen::seq(N + 1, N + 5));
		if (N % 2 == 1) {
			for (int col = 0; col < N + 1; col++) {
				V(Eigen::seq(N + 1, N + 5), col) =
					S_12.colPivHouseholderQr().solve(S_11U0(Eigen::placeholders::all, col));
			}
		}
		else {
			for (int col = 0; col < N; col++) {
				V(Eigen::seq(N + 1, N + 5), col) =
					S_12.colPivHouseholderQr().solve(S_11U0(Eigen::placeholders::all, col));
			}
			V(Eigen::seq(N + 1, N + 5), N) << 0, 8, 0, -8, 0;
		}


		// constructing W1
		V.bottomRightCorner<5, 5>() << 0, -1, 1, 0, 0, 1, -1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0,
			0, 1, 0, 0, 0;

		// constructing Delta
		LambdaPow.bottomRightCorner<5, 5>().diagonal() = (Eigen::ArrayXd(5) << pow(0.125, n - 1),
														  pow(0.125, n - 1),
														  pow(0.125, n - 1),
														  pow(0.0625, n - 1),
														  pow(0.0525, n - 1))
															 .finished();

		// constructing V_inv
	}


	if (test) {
		std::cout << "A:\n" << A_(Eigen::seq(0, N + 5), Eigen::placeholders::all) << std::endl;
		std::cout << "VLV:\n" << V * LambdaPow * V.inverse() << std::endl;
	}

	return A_*V*LambdaPow*V.inverse();
}
//
//Eigen::MatrixXd FilterSubdivFittingPlugin::matrixSubdivEigen(int N,int n)
//{
//	Eigen::MatrixXd LambdaPow = Eigen::MatrixXd::Zero(N + 6, N + 6);
//	if (N == 3) {
//
//	}
//	else {
//		double alphaN = 5 / 8 - pow((3 + 2 * cos(2 * M_PI / N)), 2) / 64;
//		double aN     = 1. - alphaN;
//		double bN     = alphaN / N;
//		double c      = 0.375;
//		double d      = 0.125;
//
//		auto f = [&N](int k) -> double { return 0.375 + 0.25 * cos(2 * M_PI * k / N); };
//
//		// constructing Sigma
//		LambdaPow(0, 0) = 1;
//		LambdaPow(1, 1) = pow(0.625 - alphaN, n);
//		int halfN       = (N % 2 == 1) ? ((N - 1) / 2) : (N / 2 - 1);
//		int eigid     = 2;
//		for (int i = 3; i <=halfN ; i++) {
//			LambdaPow(eigid, eigid)             = pow(f(i - 2), n);
//			LambdaPow(eigid + 1, eigid + 1) = pow(f(i - 2), n);
//			eigid += 2;
//		}
//		if (eigid == N)
//			LambdaPow(eigid, eigid) = pow(0.125, n);
//
//		// constructing Delta
//		LambdaPow.bottomRightCorner<5, 5>().diagonal() = (Eigen::ArrayXd(5) << pow(0.125, n),
//														  pow(0.125, n),
//														  pow(0.125, n),
//														  pow(0.0625, n),
//														  pow(0.0525, n))
//															 .finished();
//	}
//
//	return LambdaPow;
//}
//
//Eigen::MatrixXd FilterSubdivFittingPlugin::matrixSubdivEigenVector(int N)
//{
//	return Eigen::MatrixXd(0,0);
//}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterSubdivFittingPlugin)
