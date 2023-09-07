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
		if (!initflag) {
			solveFootPoints(
				*md.getMesh(par.getMeshId("samples")),
				*md.getMesh(par.getMeshId("control_mesh")),
				FootPointMode::MODE_SUBDIVISION);

			/*solvePickupVec(*md.getMesh(par.getMeshId("control_mesh")));*/
		}
		else
			log("already initialized!");

	} break;
	default :
		wrongActionCalled(action);
	}
	return std::map<std::string, QVariant>();
}

static int id(int i, int N)
{
	return ((i % N) + N) % N;
}

static int intPairHash(int m, int n)
{
	return (m + n) * (m + n + 1) / 2 + m + 1;
}

static std::pair<int, int> intPair(int h)
{
	int sum = -0.5 + 0.5 * sqrt(1. + 8. * h);
	int acc = sum * (sum + 1) / 2;
	int m   = h - acc - 1;
	int n   = sum - (h - acc) + 1;
	if (m < 0) {
		m = sum - 1;
		n = 0;
	}
	return std::make_pair(m, n);
}

void FilterSubdivFittingPlugin::solveFootPoints(
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
}

void FilterSubdivFittingPlugin::solvePickupVec(MeshModel& mm)
{
	if (!solveflag) {
		// store valence for each vertex
		CMeshO::PerVertexAttributeHandle<int> val;
		if (tri::HasPerVertexAttribute(mm.cm, "Valence")) {
			val = tri::Allocator<CMeshO>::FindPerVertexAttribute<int>(
				mm.cm, "Valence");
			if (!tri::Allocator<CMeshO>::IsValidHandle<int>(mm.cm, val)) {
				throw MLException("attribute already exists with a different type");
			}
		}
		else
			val = tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(
				mm.cm, "Valence");

		for (auto vi = mm.cm.vert.begin(); vi != mm.cm.vert.end(); vi++) {
			if (!vi->IsD()) {
				val[vi] = 0;
			}
		}

		for (auto fi = mm.cm.face.begin(); fi != mm.cm.face.end(); fi++) {
			if (!(*fi).IsD()) {
				std::cout << fi->Index() << " ";
				for (int vi = 0; vi < 3; vi++) {
					val[fi->V(vi)->Index()] += 1;
					std::cout << fi->V(vi)->Index() << " ";
				}
				std::cout << std::endl;
			}
		}

		// Construc patch pickup matrix
		CMeshO::PerFaceAttributeHandle<std::vector<Eigen::MatrixXd>> matPatchSubdiv;
		if (tri::HasPerFaceAttribute(mm.cm, "PatchSubdiv")) {
			matPatchSubdiv =
				tri::Allocator<CMeshO>::FindPerFaceAttribute<std::vector<Eigen::MatrixXd>>(
					mm.cm, "PatchSubdiv");
			if (!tri::Allocator<CMeshO>::IsValidHandle<std::vector<Eigen::MatrixXd>>(
					mm.cm, matPatchSubdiv)) {
				throw MLException("attribute already exists with a different type");
			}
		}
		else
			matPatchSubdiv =
				tri::Allocator<CMeshO>::AddPerFaceAttribute<std::vector<Eigen::MatrixXd>>(
					mm.cm, "PatchSubdiv");

		// 1.initialize to zero mat
		for (auto fi = mm.cm.face.begin(); fi != mm.cm.face.end(); fi++) {
			if (!(*fi).IsD()) {
				matPatchSubdiv[fi] = std::vector<Eigen::MatrixXd>(4);
				for (int vi = 0; vi < 3; vi++) {
					int vid           = fi->V(vi)->Index();
					matPatchSubdiv[fi][vi] = Eigen::MatrixXd::Zero(val[vid] + 6, mm.cm.vn);
				}
				matPatchSubdiv[fi][3] = Eigen::MatrixXd::Zero(12, mm.cm.vn);
			}
		}

		// 2.fill in patch subdiv mat
		for (auto fi = mm.cm.face.begin(); fi != mm.cm.face.end(); fi++) {
			if (!(*fi).IsD()) {
				matPatchSubdiv[fi] = std::vector<Eigen::MatrixXd>(4);
				for (int vi = 0; vi < 3; vi++) {
					int   vid  = fi->V(vi)->Index();
					int   pre  = fi->V(id(vi - 1, 3))->Index();
					int   next = fi->V(id(vi + 1, 3))->Index();
					auto& mats = matPatchSubdiv[fi];
					//mats[vi](0, {vid,pre,next}) = 
				}
			}
		}
	}
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
	// return dimension of 1*NV weight vector of all control points

	for (int i = 0; i < 3; i++) {
		const CVertexO* v = ft->V(i);
	}

	auto b = weightsIrregularPatch(9, v, w);
	weightsIrregularPatch(10, v, w);
	weightsIrregularPatch(11, v, w);


	return b;
}

Eigen::VectorXd FilterSubdivFittingPlugin::weightsIrregularPatch(int N, float v, float w)
{
	// return dimension of 1*(N+6) weight vector of control points in 1-ring of irregular patch
	if (N == 6)
		return weightsRegularPatch(1.f - v - w, v, w);

	if (v + w < eps) {
		// TODO: validate this formula
		double          alphaN = 5. / 8. - pow((3 + 2 * cos(2 * M_PI / N)), 2) / 64.;
		Eigen::VectorXd b = Eigen::VectorXd::Constant(N + 6, (8. / 3. * alphaN) / (1. + 8. / 3. * alphaN) / N);
		b(3) = 1. / (1. + 8. / 3. * alphaN);
		return b;
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

		return weightsRegularPatch(u, v, w) * matrixPickup(N, k) * matrixPatchSubdiv(N, n);
	}
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

Eigen::MatrixXd FilterSubdivFittingPlugin::matrixPickup(int N, int k)
{
	assert(N >= 3);
	auto cP = cacheP.find(intPairHash(N, k));
	if (cP != cacheP.end())
		return cP->second;

	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(12, N + 12);
	switch (k) {
	case 0: {
		std::vector<int> idx {
			2, 0, N + 3, 1, N, N + 8, N + 2, N + 1, N + 4, N + 7, N + 6, N + 9};
		P(Eigen::placeholders::all, idx) = Eigen::MatrixXd::Identity(12, 12);
	} break;

	case 1: {
		std::vector<int> idx {
			N + 9, N + 6, N + 4, N + 1, N + 2, N + 5, N, 1, N + 3, N - 1, 0, 2};
		P(Eigen::placeholders::all, idx) = Eigen::MatrixXd::Identity(12, 12);
	} break;

	case 2: {
		std::vector<int> idx {
			0, N - 1, 1, N, N + 5, N + 2, N + 1, N + 4, N + 11, N + 6, N + 9, N + 10};
		P(Eigen::placeholders::all, idx) = Eigen::MatrixXd::Identity(12, 12);
	} break;

	default: break;
	}
	cacheP[intPairHash(N, k)] = P;
	return P;

}

Eigen::MatrixXd FilterSubdivFittingPlugin::matrixPatchSubdiv(int N, int n)
{
	assert(N >= 3 && n >= 1);
	auto cabap = cacheAbarApow.find(intPairHash(N, n));
	if (cabap != cacheAbarApow.end())
		return cabap->second;

	double          alphaN = 5. / 8. - pow((3 + 2 * cos(2 * M_PI / N)), 2) / 64.;
	double          aN     = 1. - alphaN;
	double          bN     = alphaN / N;
	double          c      = 0.375;
	double          d      = 0.125;
	Eigen::MatrixXd A_;

	// constructing S
	auto cab = cacheAbar.find(N);
	if (cab != cacheAbar.end())
		A_ = cab->second;
	else {
		A_ = Eigen::MatrixXd::Zero(N + 12, N + 6);
		// 1-ring neighbor of EV
		A_(0, 0)                        = aN;
		A_(0, Eigen::seq(1, N)).array() = bN;
		Eigen::RowVector4d cd(c, c, d, d);
		for (int i = 1; i < N + 1; i++) {
			A_(i, {0, i, 1 + id(i - 2, N), 1 + id(i, N)}) = cd;
		}

		// constructing S_11,S_12,S_21 and S_22
		// edge vertices
		A_(N + 1, {1, N, 0, N + 1})          = cd;
		A_(N + 3, {1, 2, 0, N + 3})          = cd;
		A_(N + 5, {N - 1, N, 0, N + 5})      = cd;
		A_(N + 6, {1, N + 1, N, N + 2})      = cd;
		A_(N + 7, {1, N + 2, N + 1, N + 3})  = cd;
		A_(N + 8, {1, N + 3, 2, N + 2})      = cd;
		A_(N + 9, {N, N + 1, 1, N + 4})      = cd;
		A_(N + 10, {N, N + 4, N + 1, N + 5}) = cd;
		A_(N + 11, {N, N + 5, N - 1, N + 4}) = cd;
		// vertex vertices
		A_(N + 2, {1, 0, 2, N + 3, N + 2, N + 1, N}) =
			(Eigen::ArrayXd(7) << 0.625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625).finished();
		A_(N + 4, {N, 0, 1, N + 1, N + 4, N + 5, N - 1}) =
			(Eigen::ArrayXd(7) << 0.625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625).finished();

		cacheAbar[N] = A_;

		if (n == 1) {
			cacheAbarApow[intPairHash(N, n)] = A_;
			return A_;
		}
	}

	// constructing eigen structure of A and submatrix S
	Eigen::MatrixXd LambdaPow;
	Eigen::MatrixXd V;
	Eigen::MatrixXd V_inv;
	bool            computeV;
	auto            cV = cacheV.find(N);
	if (cV != cacheV.end()) {
		V        = cV->second;
		V_inv    = cacheVinv[N];
		computeV = false;
	}
	else {
		computeV = true;
		V        = Eigen::MatrixXd::Zero(N + 6, N + 6);
		V_inv    = Eigen::MatrixXd::Zero(N + 6, N + 6);
	}

	if (N == 3) {
		LambdaPow = Eigen::MatrixXd::Zero(9, 9);
		LambdaPow.diagonal() << 1, pow(0.25, n - 1), pow(0.25, n - 1), pow(0.125, n - 1),
			pow(0.125, n - 1), pow(0.125, n - 1), pow(0.0625, n - 1), pow(0.0625, n - 1),
			pow(0.0625, n - 1);
		LambdaPow(N + 4, N + 5) = (n - 1) / pow(16., n - 2);

		if (computeV) {
			V << 1, 0, 0, 0, 0, 0, 0, 0, 33, 1, 0, 1, 0, 0, 0, 0, 0, -22, 1, -1, -1, 0, 0, 0, 0, 0,
				-22, 1, 1, 0, 0, 0, 0, 0, 0, -22, 1, 3, 3, 1, -1, 0, 0, 0, 198, 1, 0, 4, 1, 0, 0, 0,
				165. / 16., 473, 1, -3, 0, 0, 1, 0, 0, 0, 198, 1, 4, 0, 0, 0, 1, 1, 165. / 16., 438,
				1, 0, -3, -1, 1, 1, 0, 0, 198;
			V_inv << 0.4, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0, -1. / 3., -1. / 3., 2. / 3., 0, 0, 0, 0,
				0, 0, 2. / 3., -1. / 3., -1. / 3., 0, 0, 0, 0, 0, -8, 0, 3, 3, 1, 0, 1, 0, 0, -4, 0,
				0, 3, 0, 0, 1, 0, 0, -8, 3, 3, 0, 1, 0, 0, 0, 1, 7. / 11., 26. / 33., -7. / 33.,
				-40. / 33., 0, -1, 1, 1, -1, -16. / 165., 0, 16. / 165., 16. / 165., -16. / 165.,
				16. / 165., -16. / 165., 0, 0, 1. / 55., -1. / 165., -1. / 165., -1. / 165., 0, 0,
				0, 0, 0;
			cacheV[N] = V;
			cacheVinv[N] = V_inv;
		}

		cacheAbarApow[intPairHash(N, n)] = A_ * V * LambdaPow * V_inv;

		return cacheAbarApow[intPairHash(N, n)];
	}
	else {
		LambdaPow      = Eigen::MatrixXd::Zero(N + 6, N + 6);
		auto f = [&](int k) -> double { return 0.375 + 0.25 * cos(2 * M_PI * k / N); };
		auto pairedEig = [&](int k) -> double {
			return (k == 0) ? 1 : (k == 1 ? (0.625 - alphaN) : f(k / 2));
		};

		// constructing Sigma and U0
		LambdaPow(0, 0)  = 1;
		LambdaPow(1, 1)  = pow(0.625 - alphaN, n - 1);
		int  halfN       = (N % 2 == 1) ? ((N + 3) / 2) : (N / 2 + 1);
		int  eigid      = 2;
		auto idseq       = Eigen::ArrayXd::LinSpaced(N, 0, N - 1);
		for (int i = 3; i <= halfN; i++) {
			LambdaPow(eigid, eigid)         = pow(f(i - 2), n - 1);
			LambdaPow(eigid + 1, eigid + 1) = pow(f(i - 2), n - 1);
			if (computeV) {
				V(Eigen::seq(1, N), eigid).array()     = (2 * M_PI * (i - 2.) / N * idseq).cos();
				V(Eigen::seq(1, N), eigid + 1).array() = (2 * M_PI * (i - 2.) / N * idseq).sin();
			}
			eigid += 2;
		}
		if (eigid == N) {
			LambdaPow(eigid, eigid)    = pow(0.125, n - 1);
			if (computeV) {
				V(Eigen::seq(1, N), eigid).array() = (2 * M_PI * (halfN - 1.) / N * idseq).cos();
			}
		}


		if (computeV) {
			V(Eigen::placeholders::all, {0, 1}).array() = 1;
			V(0, 1)                                     = -8. / 3. * alphaN;

					// constructing U1
			auto S_11U0 =
				A_(Eigen::seq(N + 1, N + 5), Eigen::seq(0, N)) * V.topLeftCorner(N + 1, N + 1);
			auto S_12 = A_(Eigen::seq(N + 1, N + 5), Eigen::seq(N + 1, N + 5));
			if (N % 2 == 1) {
				for (int col = 0; col < N + 1; col++) {
					V(Eigen::seq(N + 1, N + 5), col) =
						(pairedEig(col) * Eigen::MatrixXd::Identity(5, 5) - S_12)
							.colPivHouseholderQr()
							.solve(S_11U0(Eigen::placeholders::all, col));
					std::cout << pairedEig(col) << std::endl;
				}
			}
			else {
				for (int col = 0; col < N; col++) {
					V(Eigen::seq(N + 1, N + 5), col) =
						(pairedEig(col) * Eigen::MatrixXd::Identity(5, 5) - S_12)
							.colPivHouseholderQr()
							.solve(S_11U0(Eigen::placeholders::all, col));
				}
				V(Eigen::seq(N + 1, N + 5), N) << 0, 8, 0, -8, 0;
			}

			// constructing W1
			V.bottomRightCorner<5, 5>() << 0, -1, 1, 0, 0, 1, -1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1,
				1, 0, 0, 1, 0, 0, 0;

			cacheV[N]    = V;
			cacheVinv[N] = V.inverse();
		}

		// constructing Delta
		LambdaPow.bottomRightCorner<5, 5>().diagonal() = (Eigen::ArrayXd(5) << pow(0.125, n - 1),
														  pow(0.125, n - 1),
														  pow(0.125, n - 1),
														  pow(0.0625, n - 1),
														  pow(0.0625, n - 1))
															 .finished();
	}

	cacheAbarApow[intPairHash(N, n)] = A_ * V * LambdaPow * cacheVinv[N];
	return cacheAbarApow[intPairHash(N, n)];
}


MESHLAB_PLUGIN_NAME_EXPORTER(FilterSubdivFittingPlugin)
