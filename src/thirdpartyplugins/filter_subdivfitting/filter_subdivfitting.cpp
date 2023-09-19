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

#define MATEPS 1e-10


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

static double alpha(int N)
{
	return 0.625 - pow((3 + 2 * cos(2 * M_PI / N)), 2) / 64.;
}

/**
 * @brief
 * Constructor usually performs only two simple tasks of filling the two lists
 *  - typeList: with all the possible id of the filtering actions
 *  - actionList with the corresponding actions.
 * If you want to add icons to your filtering actions you can do here by construction the QActions accordingly
 */
FilterSubdivFittingPlugin::FilterSubdivFittingPlugin()
{ 
	typeList = {
		FP_INIT,
		FP_SUBDIV_FITTING,
		FP_REANALYSIS_CA,
		FP_FITTING_ERROR,
		FP_FITTING_CACHE_CLEAR,
		FP_SIMPLE_SAMPLE_DENSIFY,
		FP_QUALITY_TRANSFFER,
		FP_ADD_SAMPLES
	};

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
	case FP_INIT: return "Fitting: Initialize";
	case FP_SUBDIV_FITTING :
		return "Fitting: Subdivision Surface Fitting";
	case FP_FITTING_ERROR: return "Fitting: Render Distance Error";
	case FP_FITTING_CACHE_CLEAR: return "Fitting: Clear Cache";
	case FP_SIMPLE_SAMPLE_DENSIFY: return "Fitting: Densify Samples";
	case FP_QUALITY_TRANSFFER: return "Fitting: Transfer Samples Quality";
	case FP_ADD_SAMPLES: return "Fitting: Add Samples";
	case FP_REANALYSIS_CA: return "Fitting: Reanalysis";
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
	case FP_FITTING_ERROR: return "Coloring fitting error";
	case FP_INIT:
	case FP_FITTING_CACHE_CLEAR:
	case FP_SIMPLE_SAMPLE_DENSIFY:
	case FP_QUALITY_TRANSFFER:
	case FP_ADD_SAMPLES:
	case FP_REANALYSIS_CA:
		return "";
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
	case FP_INIT:
	case FP_SUBDIV_FITTING :
	case FP_REANALYSIS_CA:
	case FP_FITTING_ERROR:
	case FP_FITTING_CACHE_CLEAR:
	case FP_SIMPLE_SAMPLE_DENSIFY:
	case FP_QUALITY_TRANSFFER:
	case FP_ADD_SAMPLES:
		return FilterPlugin::SubdivFitting;
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
 * @brief FilterSamplePlugin::getRequirements
 * @return
 */
int FilterSubdivFittingPlugin::getRequirements(const QAction* act)
{
	switch (ID(act)) {
	case FP_SUBDIV_FITTING:
	case FP_REANALYSIS_CA:
		return MeshModel::MM_FACEFACETOPO | MeshModel::MM_VERTFACETOPO;
	case FP_INIT:
	case FP_FITTING_ERROR:
	case FP_FITTING_CACHE_CLEAR:
	case FP_SIMPLE_SAMPLE_DENSIFY:
	case FP_QUALITY_TRANSFFER:
	case FP_ADD_SAMPLES:
		return 0;
	default: assert(0); return 0;
	}
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
int FilterSubdivFittingPlugin::postCondition(const QAction* action) const
{
	switch (ID(action)) {
	case FP_FITTING_ERROR: return MeshModel::MM_VERTQUALITY + MeshModel::MM_VERTCOLOR;
	case FP_INIT:
	case FP_FITTING_CACHE_CLEAR:
	case FP_SIMPLE_SAMPLE_DENSIFY:
	case FP_QUALITY_TRANSFFER:
	case FP_ADD_SAMPLES:
	case FP_REANALYSIS_CA:
	default: break;
	}
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
	case FP_INIT: {
		const MeshModel* target = md.mm();
		// looking for a second mesh different that the current one
		for (const MeshModel& t : md.meshIterator()) {
			if (&t != md.mm()) {
				target = &t; // control mesh topo
				break;
			}
		}
		parlst.addParam(
			RichMesh("source_mesh", md.mm()->id(), &md, "Source Mesh", "Mesh to be sampled"));
		parlst.addParam(RichMesh("samples", md.mm()->id(), &md, "Samples", "Samples to be fitted"));
		parlst.addParam(RichMesh(
			"control_mesh", target->id(), &md, "Control mesh", "Control mesh for subdivision"));
	} break;

	case FP_SUBDIV_FITTING: {
		parlst.addParam(RichBool(
			"no_fitting",
			false,
			"Push to Limit Position only",
			" "
			" "));
		parlst.addParam(RichBool(
			"show_new_mesh",
			false,
			"Show New Control Mesh",
			" "
			" "));
		parlst.addParam(RichBool(
			"show_limit_samples",
			false,
			"Show Limit Samples",
			" "
			" "));
		parlst.addParam(RichBool(
			"check_bary_coord",
			false,
			"Check Barycenter",
			" "
			" "));
		}
		break;
	case FP_REANALYSIS_CA: {
			parlst.addParam(RichInt(
				"CA_rank",
				4,
				"Rank of CA",
				" "
				" "));
	} break;
	case FP_FITTING_ERROR: {
		parlst.addParam(RichMesh("fitting_samples", md.mm()->id(), &md, "Fitting Samples", ""));
		parlst.addParam(RichMesh("dest_mesh", md.mm()->id(), &md, "Destination Mesh", ""));
	} break;
	case FP_FITTING_CACHE_CLEAR: {
	} break;
	case FP_SIMPLE_SAMPLE_DENSIFY: {
		parlst.addParam(RichMesh("tobedensify", md.mm()->id(), &md, "Samples to be Densify", ""));
		parlst.addParam(RichMesh("ref_mesh", md.mm()->id(), &md, "Referrence Mesh", ""));
	} break;
	case FP_QUALITY_TRANSFFER: {
		parlst.addParam(RichMesh("from", md.mm()->id(), &md, "From which samples", ""));
		parlst.addParam(RichMesh("to", md.mm()->id(), &md, "To which samples", ""));
	} break;
	case FP_ADD_SAMPLES: {
		parlst.addParam(RichMesh("add_samples", md.mm()->id(), &md, "Samples to be added", ""));
	} break;
	default :
		assert(0);
	}
	return parlst;
}

void FilterSubdivFittingPlugin::clearFittingCache()
{
	if (!initflag)
		return;

	tri::Allocator<CMeshO>::DeletePerVertexAttribute(ptsample->cm, "SampleUpdate");
	tri::Allocator<CMeshO>::DeletePerVertexAttribute(ptsample->cm, "BaryCoord");
	tri::Allocator<CMeshO>::DeletePerVertexAttribute(ptsample->cm, "FootTriangle");
	tri::Allocator<CMeshO>::DeletePerVertexAttribute(ptsample->cm, "LimitStencil");
	tri::Allocator<CMeshO>::DeletePerVertexAttribute(ptctrlmesh->cm, "Valence");
	tri::Allocator<CMeshO>::DeletePerVertexAttribute(ptctrlmesh->cm, "ControlMeshUpdate");
	tri::Allocator<CMeshO>::DeletePerFaceAttribute(ptctrlmesh->cm, "PatchSubdiv");

	initflag     = false;
	topochange   = true;
	sampleupdate = true;
	solveflag    = false;

	cacheP.clear();
	cacheAbar.clear();
	cacheV.clear();
	cacheVinv.clear();
	cacheAbarApow.clear();

	mdptr           = nullptr;
	ptsource        = nullptr;
	ptsample        = nullptr;
	ptctrlmesh      = nullptr;
	fittingres      = nullptr;

	splpts          = Eigen::MatrixXd();
	projectedsplpts = Eigen::MatrixXd();
	controlmesh     = Eigen::MatrixXd();
	AT              = Eigen::SparseMatrix<double>();
	ATA             = Eigen::SparseMatrix<double>();
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
	MeshModel* curMM = md.mm();
	CMeshO&    m     = curMM->cm;

	mdptr      = &md;

	switch(ID(action)) {
	case FP_INIT: {
		if (initflag)
			clearFittingCache();

		ptsource   = md.getMesh(par.getMeshId("source_mesh"));
		ptsample   = md.getMesh(par.getMeshId("samples"));
		ptctrlmesh = md.getMesh(par.getMeshId("control_mesh"));
		oldvn      = ptsample->cm.VN();

		if (!initflag) {
			assignPerElementAtributes();
			initflag = true;
		}

	} break;
	case FP_SUBDIV_FITTING: {
		if (!initflag)
			throw MLException("Not initiated!");

		if (topochange){
			updateControlVertexAttribute();
			solvePickupVec();
			updateVertexComplete(ptctrlmesh, "ControlMeshUpdate");
			log("topochange");
			topochange = false;
		}

		if (sampleupdate) {
			parameterizeSamples(FootPointMode::MODE_MESH);
			updateLimitStencils(UpdateOptions::MODE_INIT);
			updateVertexComplete(ptsample, "SampleUpdate");
			log("sampleupdate");
			sampleupdate = false;
		}


		if (!solveflag) {
			assembleFittingQuery(par);
			solveflag = true;
		}
		else {
		}

		displayResults(par);
	} break;

	case FP_REANALYSIS_CA: {
		int rank = par.getInt("CA_rank");
		assembleIncrement(rank);
	} break;

	case FP_FITTING_ERROR: {
		auto fittingspls = md.getMesh(par.getMeshId("fitting_samples"));
		auto destmesh = md.getMesh(par.getMeshId("dest_mesh"));

		fittingspls->updateDataMask(MeshModel::MM_VERTQUALITY);

		for (auto& vi = fittingspls->cm.vert.begin(); vi != fittingspls->cm.vert.end(); vi++) {
			if (!vi->IsD()) {
				vi->Q() = std::numeric_limits<float>::max();
				for (auto& fi = destmesh->cm.face.begin(); fi != destmesh->cm.face.end(); fi++) {
					auto p = distancePointTriangle(*vi, *fi);
					if (p.first < vi->Q())
						vi->Q() = p.first;
				}
			}
		}

		//tri::UpdateQuality<CMeshO>::VertexNormalize(fittingspls->cm);
		tri::UpdateColor<CMeshO>::PerVertexQualityRamp(fittingspls->cm);
		fittingspls->updateDataMask(MeshModel::MM_VERTCOLOR);

	} break;

	case FP_FITTING_CACHE_CLEAR: {
		clearFittingCache();
	} break;

	case FP_SIMPLE_SAMPLE_DENSIFY: {
		auto tobedensify = md.getMesh(par.getMeshId("tobedensify"));
		auto refmesh     = md.getMesh(par.getMeshId("ref_mesh"));
	} break;

	case FP_QUALITY_TRANSFFER: {
		auto from_ = md.getMesh(par.getMeshId("from"));
		auto to_   = md.getMesh(par.getMeshId("to"));

		to_->updateDataMask(MeshModel::MM_VERTQUALITY);

		for (int vi = 0; vi < from_->cm.vert.size(); vi++) {
			to_->cm.vert[vi].Q() = from_->cm.vert[vi].Q();
		}

		tri::UpdateColor<CMeshO>::PerVertexQualityRamp(to_->cm);
		to_->updateDataMask(MeshModel::MM_VERTCOLOR);
	} break;
	case FP_ADD_SAMPLES: {
		if (!initflag)
			throw MLException("Not initiated!");

		sampleupdate = true;
		auto tobeadd = md.getMesh(par.getMeshId("add_samples"));
		auto oldVN   = ptsample->cm.VN();
		auto vi      = tri::Allocator<CMeshO>::AddVertices(ptsample->cm, tobeadd->cm.VN());

		for (int i = 0; i < tobeadd->cm.VN(); i++, vi++) {
			vi->P() = tobeadd->cm.vert[i].P();
		}

		auto tobeupdate =
			tri::Allocator<CMeshO>::FindPerVertexAttribute<bool>(ptsample->cm, "SampleUpdate");
		for (int i = oldVN; i < ptsample->cm.VN(); i++) {
			tobeupdate[i] = true;
		}

		auto ftptrs = tri::Allocator<CMeshO>::FindPerVertexAttribute<const CFaceO*>(
			ptsample->cm, "FootTriangle");
		for (int i = oldVN; i < ptsample->cm.VN(); i++) {
			ftptrs[i] = nullptr;
		}

		auto ftbarycoord =
			tri::Allocator<CMeshO>::FindPerVertexAttribute<Point3f>(ptsample->cm, "BaryCoord");
		for (int i = oldVN; i < ptsample->cm.VN(); i++) {
			ftbarycoord[i] = Point3f(0, 0, 0);
		}

		auto ls = tri::Allocator<CMeshO>::FindPerVertexAttribute<Eigen::SparseVector<double>>(
			ptsample->cm, "LimitStencil");

		parameterizeSamples(FootPointMode::MODE_MESH);
		updateLimitStencils(UpdateOptions::MODE_INIT);
		updateVertexComplete(ptsample, "SampleUpdate");
		sampleupdate = false;
	} break;
	default :
		wrongActionCalled(action);
	}
	return std::map<std::string, QVariant>();
}

void FilterSubdivFittingPlugin::assignPerElementAtributes()
{
	CMeshO::PerVertexAttributeHandle<const CFaceO*> ftptrs;
	if (tri::HasPerVertexAttribute(ptsample->cm, "FootTriangle")) {
		ftptrs = tri::Allocator<CMeshO>::FindPerVertexAttribute<const CFaceO*>(
			ptsample->cm, "FootTriangle");
		if (!tri::Allocator<CMeshO>::IsValidHandle<const CFaceO*>(ptsample->cm, ftptrs)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		ftptrs = tri::Allocator<CMeshO>::AddPerVertexAttribute<const CFaceO*>(
			ptsample->cm, "FootTriangle");

	for (auto vi = ptsample->cm.vert.begin(); vi != ptsample->cm.vert.end(); vi++) {
		if (!(*vi).IsD()) {
			ftptrs[vi] = nullptr;
		}
	}


	CMeshO::PerVertexAttributeHandle<Point3f> ftbarycoords;
	if (tri::HasPerVertexAttribute(ptsample->cm, "BaryCoord")) {
		ftbarycoords =
			tri::Allocator<CMeshO>::FindPerVertexAttribute<Point3f>(ptsample->cm, "BaryCoord");
		if (!tri::Allocator<CMeshO>::IsValidHandle<Point3f>(ptsample->cm, ftbarycoords)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		ftbarycoords =
			tri::Allocator<CMeshO>::AddPerVertexAttribute<Point3f>(ptsample->cm, "BaryCoord");

	for (auto vi = ptsample->cm.vert.begin(); vi != ptsample->cm.vert.end(); vi++) {
		if (!(*vi).IsD()) {
			ftbarycoords[vi] = Point3f(0, 0, 0);
		}
	}


	CMeshO::PerVertexAttributeHandle<int> val;
	if (tri::HasPerVertexAttribute(ptctrlmesh->cm, "Valence")) {
		val = tri::Allocator<CMeshO>::FindPerVertexAttribute<int>(ptctrlmesh->cm, "Valence");
		if (!tri::Allocator<CMeshO>::IsValidHandle<int>(ptctrlmesh->cm, val)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		val = tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(ptctrlmesh->cm, "Valence");

	for (auto vi = ptctrlmesh->cm.vert.begin(); vi != ptctrlmesh->cm.vert.end(); vi++) {
		if (!vi->IsD()) {
			val[vi] = 0;
		}
	}


	// Register patch pickup matrix attribute
	CMeshO::PerFaceAttributeHandle<std::vector<Eigen::MatrixXd>> matPatchSubdiv;
	if (tri::HasPerFaceAttribute(ptctrlmesh->cm, "PatchSubdiv")) {
		matPatchSubdiv = tri::Allocator<CMeshO>::FindPerFaceAttribute<std::vector<Eigen::MatrixXd>>(
			ptctrlmesh->cm, "PatchSubdiv");
		if (!tri::Allocator<CMeshO>::IsValidHandle<std::vector<Eigen::MatrixXd>>(
				ptctrlmesh->cm, matPatchSubdiv)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		matPatchSubdiv = tri::Allocator<CMeshO>::AddPerFaceAttribute<std::vector<Eigen::MatrixXd>>(
			ptctrlmesh->cm, "PatchSubdiv");


	CMeshO::PerVertexAttributeHandle<Eigen::SparseVector<double>> ls;
	if (tri::HasPerVertexAttribute(ptsample->cm, "LimitStencil")) {
		ls = tri::Allocator<CMeshO>::FindPerVertexAttribute<Eigen::SparseVector<double>>(
			ptsample->cm, "LimitStencil");
		if (!tri::Allocator<CMeshO>::IsValidHandle<Eigen::SparseVector<double>>(ptsample->cm, ls)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		ls = tri::Allocator<CMeshO>::AddPerVertexAttribute<Eigen::SparseVector<double>>(
			ptsample->cm, "LimitStencil");


	CMeshO::PerVertexAttributeHandle<bool> tobeupdate;
	if (tri::HasPerVertexAttribute(ptsample->cm, "SampleUpdate")) {
		tobeupdate =
			tri::Allocator<CMeshO>::FindPerVertexAttribute<bool>(ptsample->cm, "SampleUpdate");
		if (!tri::Allocator<CMeshO>::IsValidHandle<bool>(ptsample->cm, tobeupdate)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		tobeupdate = tri::Allocator<CMeshO>::AddPerVertexAttribute<bool>(ptsample->cm, "SampleUpdate");

	for (auto vi = ptsample->cm.vert.begin(); vi != ptsample->cm.vert.end(); vi++) {
		if (!vi->IsD()) {
			tobeupdate[vi] = true;
		}
	}

		CMeshO::PerVertexAttributeHandle<bool> tobeupdate2;
	if (tri::HasPerVertexAttribute(ptctrlmesh->cm, "ControlMeshUpdate")) {
		tobeupdate2 =
			tri::Allocator<CMeshO>::FindPerVertexAttribute<bool>(ptctrlmesh->cm, "ControlMeshUpdate");
		if (!tri::Allocator<CMeshO>::IsValidHandle<bool>(ptctrlmesh->cm, tobeupdate2)) {
			throw MLException("attribute already exists with a different type");
		}
	}
	else
		tobeupdate2 =
			tri::Allocator<CMeshO>::AddPerVertexAttribute<bool>(ptctrlmesh->cm, "ControlMeshUpdate");

	for (auto vi = ptctrlmesh->cm.vert.begin(); vi != ptctrlmesh->cm.vert.end(); vi++) {
		if (!vi->IsD()) {
			tobeupdate2[vi] = true;
		}
	}
}

void FilterSubdivFittingPlugin::parameterizeSamples(FootPointMode mode)
{
	auto tobeupdate =
		tri::Allocator<CMeshO>::FindPerVertexAttribute<bool>(ptsample->cm, "SampleUpdate");

	for (auto si = ptsample->cm.vert.begin(); si != ptsample->cm.vert.end(); si++) {
		if ((!si->IsD()) && tobeupdate[si])
			solveFootPoint(&(*si), mode);
	}
}

void FilterSubdivFittingPlugin::solveFootPoint(CVertexO* v, FootPointMode mode)
{
	auto ftptrs =
		tri::Allocator<CMeshO>::FindPerVertexAttribute<const CFaceO*>(ptsample->cm, "FootTriangle");

	auto ftbarycoords =
		tri::Allocator<CMeshO>::FindPerVertexAttribute<Point3f>(ptsample->cm, "BaryCoord");

	switch (mode) {
	case FilterSubdivFittingPlugin::MODE_MESH: {
		// solve foot point by brute force
		std::pair<float, vcg::Point3f> minpair(
			std::numeric_limits<float>::max(), Point3f(0.f, 0.f, 0.f));
		const CFaceO* footface = nullptr;
		for (auto fi = ptctrlmesh->cm.face.begin(); fi != ptctrlmesh->cm.face.end(); fi++) {
			auto d = distancePointTriangle(*v, *fi);
			if (d.first < minpair.first) {
				minpair  = d;
				footface = &(*fi);
			}
		}
		ftptrs[v->Index()]       = footface;
		ftbarycoords[v->Index()] = minpair.second;
	} break;
	case FilterSubdivFittingPlugin::MODE_SUBDIVISION: {
		std::pair<float, vcg::Point3f> minpair(
			std::numeric_limits<float>::max(), Point3f(0.f, 0.f, 0.f));
		const CFaceO* footface = nullptr;
		for (auto fi = ptctrlmesh->cm.face.begin(); fi != ptctrlmesh->cm.face.end(); fi++) {
			auto d = distancePointTriangle(*v, *fi);
			if (d.first < minpair.first) {
				minpair  = d;
				footface = &(*fi);
			}
		}
		ftptrs[v->Index()] = footface;
		ftbarycoords[v->Index()] = minpair.second;
	} break;
	default: break;
	}
}

void FilterSubdivFittingPlugin::updateControlVertexAttribute()
{
	ptctrlmesh->cm.face.EnableFFAdjacency();
	tri::UpdateTopology<CMeshO>::FaceFace(ptctrlmesh->cm);

	// store valence for each vertex
	auto val = tri::Allocator<CMeshO>::FindPerVertexAttribute<int>(ptctrlmesh->cm, "Valence");
	auto tobeupdate =
		tri::Allocator<CMeshO>::FindPerVertexAttribute<bool>(ptctrlmesh->cm, "ControlMeshUpdate");

	for (auto vi = ptctrlmesh->cm.vert.begin(); vi != ptctrlmesh->cm.vert.end(); vi++) {
		if ((!vi->IsD()) && tobeupdate[vi]) {
			val[vi] = 0;
		}
	}

	// TODO: not updating precisely
	for (auto fi = ptctrlmesh->cm.face.begin(); fi != ptctrlmesh->cm.face.end(); fi++) {
		if (!fi->IsD()) {
			val[fi->V(0)] += 1;
			val[fi->V(1)] += 1;
			val[fi->V(2)] += 1;
		}
	}

	splpts = Eigen::MatrixXd::Zero(ptsample->cm.VN(), 3);

	for (int si = 0; si != ptsample->cm.vert.size(); si++) {
		splpts(si, Eigen::placeholders::all) << ptsample->cm.vert[si].P()[0],
			ptsample->cm.vert[si].P()[1], ptsample->cm.vert[si].P()[2];
	}
}

void FilterSubdivFittingPlugin::solvePickupVec()
{
	auto  val = tri::Allocator<CMeshO>::FindPerVertexAttribute<int>(ptctrlmesh->cm, "Valence");

	// Register patch pickup matrix attribute
	auto matPatchSubdiv =
		tri::Allocator<CMeshO>::FindPerFaceAttribute<std::vector<Eigen::MatrixXd>>(
			ptctrlmesh->cm, "PatchSubdiv");

	// initialize subdiv mat
	for (auto fi = ptctrlmesh->cm.face.begin(); fi != ptctrlmesh->cm.face.end(); fi++) {
		if (!(*fi).IsD()) {
			matPatchSubdiv[fi].resize(4);
			for (int vi = 0; vi < 3; vi++) {
				// allocate zero mat
				int vid                = fi->V(vi)->Index();
				matPatchSubdiv[fi][vi] = Eigen::MatrixXd::Zero(val[vid] + 6, ptctrlmesh->cm.vn);
			}
			matPatchSubdiv[fi][3] = Eigen::MatrixXd::Zero(12, ptctrlmesh->cm.vn);
		}
	}

	Eigen::Array4d tempedge({0.375, 0.375, 0.125, 0.125});
	for (auto fi = ptctrlmesh->cm.face.begin(); fi != ptctrlmesh->cm.face.end(); fi++) {
		if (!(*fi).IsD()) {
			CFaceO*          startfp = &*fi;
			std::vector<int> ring1, ring2, ringNp;

			// construct irregular patch subdiv mat
			for (int vi = 0; vi < 3; vi++) {
				// pick vertice's id in rings
				int vid  = fi->V(vi)->Index();
				int prev = fi->V(fi->Prev(vi))->Index();
				int next = fi->V(fi->Next(vi))->Index();
				int id1  = vi;
				int id2  = fi->Next(id1);
				int idNp = fi->Prev(id1);
				ring1.clear();
				ring2.clear();
				ringNp.clear();
				ring1.reserve(val[fi->V0(id1)]);
				ring2.reserve(val[fi->V0(id2)]);
				ringNp.reserve(val[fi->V0(idNp)]);

				CFaceO*                fp = startfp;
				vcg::face::Pos<CFaceO> p(fp, id1, startfp->V0(id1));

				// ring1 circulate
				do {
					ring1.push_back(fp->V(fp->Next(p.z))->Index());
					p.FlipF();
					p.FlipE();
					fp = p.F();
				} while (fp != startfp);

				// ring2 criculate
				p.Set(fp, id1, startfp->V0(id2));
				do {
					ring2.push_back(fp->V(p.z)->Index());
					p.FlipE();
					p.FlipF();
					fp = p.F();
				} while (fp != startfp);

				// ringNp circulate
				p.Set(fp, idNp, startfp->V0(idNp));
				do {
					ringNp.push_back(fp->V(fp->Next(p.z))->Index());
					p.FlipF();
					p.FlipE();
					fp = p.F();
				} while (fp != startfp);

				// start filling mats
				int N                                    = val[vid];
				matPatchSubdiv[fi][vi](0, vid)           = 1 - alpha(N);
				matPatchSubdiv[fi][vi](0, ring1).array() = alpha(N) / N;
				for (int i = 0; i < N; i++) {
					matPatchSubdiv[fi][vi](
						i + 1, {vid, ring1[i], ring1[id(i - 1, N)], ring1[id(i + 1, N)]})
						.array() = tempedge;
				}
				matPatchSubdiv[fi][vi](N + 1, {ring1[0], ring1[id(-1, N)]}).array()    = 0.375;
				matPatchSubdiv[fi][vi](N + 1, {vid, ring2[2]}).array()                 = 0.125;
				matPatchSubdiv[fi][vi](N + 3, {ring1[0], ring1[1]}).array()            = 0.375;
				matPatchSubdiv[fi][vi](N + 3, {vid, ring2[id(-2, val[next])]}).array() = 0.125;
				matPatchSubdiv[fi][vi](N + 5, {ring1[id(-1, N)], ring1[id(-2, N)]}).array() =
					0.375;
				matPatchSubdiv[fi][vi](N + 5, {vid, ringNp[2]}).array() = 0.125;

				matPatchSubdiv[fi][vi](N + 2, next)           = 1 - alpha(val[next]);
				matPatchSubdiv[fi][vi](N + 2, ring2).array()  = alpha(val[next]) / val[next];
				matPatchSubdiv[fi][vi](N + 4, prev)           = 1 - alpha(val[prev]);
				matPatchSubdiv[fi][vi](N + 4, ringNp).array() = alpha(val[prev]) / val[prev];
			}
			// std::cout << "face: " << fi->Index() << " vert0: " << fi->V(0)->Index() <<
			// std::endl
			//		  << matPatchSubdiv[fi][0] << std::endl
			//		  << std::endl;

			////test
			// std::cout << "face: " << fi->Index()<<std::endl;
			// std::cout << "vert0: " << fi->V(0)->Index();
			// std::cout << std::endl << "ring1: ";
			// for (int i : ring1)
			//	std::cout << i << " ";
			// std::cout << std::endl << "ring2: ";
			// for (int i : ring2)
			//	std::cout << i << " ";
			// std::cout << std::endl << "ringNp: ";
			// for (int i : ringNp)
			//	std::cout << i << " ";
			// std::cout << std::endl << std::endl;

			// construct regular patch subdiv mat

			// vert 1
			matPatchSubdiv[fi][3](
				0,
				{fi->V2(0)->Index(),
					ring1[id(-2, ring1.size())],
					ring1[id(-1, ring1.size())],
					ring1[id(-3, ring1.size())]})
				.array() = tempedge;

			// vert 2
			matPatchSubdiv[fi][3](1, {fi->V2(2)->Index(), ringNp[1], ringNp[0], ringNp[2]})
				.array() = tempedge;

			// vert 3
			matPatchSubdiv[fi][3](2, fi->V2(0)->Index()) = 1 - alpha(ring1.size());
			matPatchSubdiv[fi][3](2, ring1).array()      = alpha(ring1.size()) / ring1.size();

			// vert 4
			matPatchSubdiv[fi][3](
				3, {fi->V2(0)->Index(), fi->V2(2)->Index(), ring1[0], ringNp[1]})
				.array() = tempedge;

			// vert 5
			matPatchSubdiv[fi][3](4, fi->V2(2)->Index()) = 1 - alpha(ringNp.size());
			matPatchSubdiv[fi][3](4, ringNp).array()     = alpha(ringNp.size()) / ringNp.size();

			// vert 6
			matPatchSubdiv[fi][3](5, {fi->V2(0)->Index(), ring1[1], ring1[0], ring1[2]})
				.array() = tempedge;

			// vert 7
			matPatchSubdiv[fi][3](
				6, {fi->V2(0)->Index(), ring1[0], ring1[1], ring1[id(-1, ring1.size())]})
				.array() = tempedge;

			// vert 8
			matPatchSubdiv[fi][3](
				7,
				{fi->V2(2)->Index(),
					ringNp[id(-1, ringNp.size())],
					ringNp[0],
					ringNp[id(-2, ringNp.size())]})
				.array() = tempedge;

			// vert 9
			matPatchSubdiv[fi][3](
				8,
				{fi->V2(2)->Index(),
					ringNp[id(-2, ringNp.size())],
					ringNp[id(-1, ringNp.size())],
					ringNp[id(-3, ringNp.size())]})
				.array() = tempedge;

			// vert 10dev
			matPatchSubdiv[fi][3](
				9,
				{fi->V2(1)->Index(),
					ring2[id(-1, ring2.size())],
					ring2[0],
					ring2[id(-2, ring2.size())]})
				.array() = tempedge;

			// vert 11
			matPatchSubdiv[fi][3](10, fi->V2(1)->Index()) = 1 - alpha(ring2.size());
			matPatchSubdiv[fi][3](10, ring2).array()      = alpha(ring2.size()) / ring2.size();

			// vert 12
			matPatchSubdiv[fi][3](
				11, {fi->V2(1)->Index(), ring2[2], ring2[1], ring2[id(3, ring2.size())]})
				.array() = tempedge;
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

	float p01 = (project - v0).dot((v1 - v0).normalized());
	float p10 = (project - v1).dot((v0 - v1).normalized());
	float p02 = (project - v0).dot((v2 - v0).normalized());
	float p20 = (project - v2).dot((v0 - v2).normalized());
	float p12 = (project - v1).dot((v2 - v1).normalized());
	float p21 = (project - v2).dot((v1 - v2).normalized());

	auto det = [](const Point3f& v0, const Point3f& v1, const Point3f& v2) -> float {
		return v0[0] * v1[1] * v2[2] + v0[1] * v1[2] * v2[0] + v0[2] * v1[0] * v2[1] -
			   v0[2] * v1[1] * v2[0] - v0[1] * v1[0] * v2[2] - v0[0] * v1[2] * v2[1];
	};

	// case corner
	if ((p01 < 0.f) && (p02 < 0.f)) {
		squared_parallel_dist = (project - v0).SquaredNorm();
		l0                    = 1.f;
		l1                    = 0.f;
		l2                    = 0.f;
	}
	else if ((p12 < 0.f) && (p10 < 0.f)) {
		squared_parallel_dist = (project - v1).SquaredNorm();
		l0                    = 0.f;
		l1                    = 1.f;
		l2                    = 0.f;
	}
	else if ((p20 < 0.f) && (p21 < 0.f)) {
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

		if ((l0 < 0.f) && (p12 * p21 >= 0.f)) {
			l0                    = 0.f;
			l1                    = p21 / (v1 - v2).Norm();
			l2                    = 1.f - l1;
			squared_parallel_dist = (project - v2).SquaredNorm() - p21 * p21;
		}
		else if ((l1 < 0.f) && (p02 * p20 >= 0.f)) {
			l1                    = 0.f;
			l2                    = p02 / (v2 - v0).Norm();
			l0                    = 1.f - l2;
			squared_parallel_dist = (project - v0).SquaredNorm() - p02 * p02;
		}
		else if ((l2 < 0.f) && (p01 * p10 >= 0.f)) {
			l2                    = 0.f;
			l0                    = p10 / (v0 - v1).Norm();
			l1                    = 1.f - l0;
			squared_parallel_dist = (project - v1).SquaredNorm() - p10 * p10;
		}
		else
			squared_parallel_dist = 0.f;
	}

	return std::pair<float, Point3f>(
		sqrtf(height * height + squared_parallel_dist), Point3f(l0, l1, l2));
}

void FilterSubdivFittingPlugin::updateVertexComplete(MeshModel* mm, std::string field)
{
	auto tobeupdate = tri::Allocator<CMeshO>::FindPerVertexAttribute<bool>(mm->cm, field);
	for (auto vi = mm->cm.vert.begin(); vi != mm->cm.vert.end(); vi++) {
		if (!vi->IsD()) {
			tobeupdate[vi] = false;
		}
	}
}

void FilterSubdivFittingPlugin::updateLimitStencils(UpdateOptions mode) {
	switch (mode) {
	case FilterSubdivFittingPlugin::MODE_INIT: {
		auto ls = tri::Allocator<CMeshO>::FindPerVertexAttribute<Eigen::SparseVector<double>>(
			ptsample->cm, "LimitStencil");

		auto ftptrs = tri::Allocator<CMeshO>::FindPerVertexAttribute<const CFaceO*>(
			ptsample->cm, "FootTriangle");

		auto ftbarycoords =
			tri::Allocator<CMeshO>::FindPerVertexAttribute<Point3f>(ptsample->cm, "BaryCoord");

		auto ControlMeshUpdate =
			tri::Allocator<CMeshO>::FindPerVertexAttribute<bool>(ptsample->cm, "SampleUpdate");

		// foot points should be solved now
		for (auto vi = ptsample->cm.vert.begin(); vi != ptsample->cm.vert.end(); vi++) {
			if ((!vi->IsD()) && ControlMeshUpdate[vi]) {
				Eigen::VectorXd densestencil =
					weightsPatch(ftptrs[vi], ftbarycoords[vi][1], ftbarycoords[vi][2]);
				ls[vi] = densestencil.sparseView(1., MATEPS);
			}
		}
	} break;
	case FilterSubdivFittingPlugin::MODE_UPDATE: break;
	default: break;
	}
}

void FilterSubdivFittingPlugin::assembleFittingQuery(const RichParameterList& par)
{
	projectedsplpts = Eigen::MatrixXd::Zero(ptctrlmesh->cm.VN(), 3);

	CMeshO::PerVertexAttributeHandle<Eigen::SparseVector<double>> ls;
	if (tri::HasPerVertexAttribute(ptsample->cm, "LimitStencil")) {
		ls = tri::Allocator<CMeshO>::FindPerVertexAttribute<Eigen::SparseVector<double>>(
			ptsample->cm, "LimitStencil");
		if (!tri::Allocator<CMeshO>::IsValidHandle<Eigen::SparseVector<double>>(ptsample->cm, ls)) {
			throw MLException("attribute already exists with a different type");
		}
	}

	AT.resize(ptctrlmesh->cm.vert.size(), ptsample->cm.vert.size());
	for (int si = 0; si != ptsample->cm.vert.size(); si++) {
		AT.col(si) = ls[si];
		projectedsplpts += ls[si] * splpts(si, Eigen::placeholders::all);
	}

	bool no_fitting = par.getBool("no_fitting");

	if (no_fitting) {
		controlmesh.resize(ptctrlmesh->cm.vert.size(), 3);
		for (int vi = 0; vi != ptctrlmesh->cm.vert.size(); vi++) {
			controlmesh(vi, Eigen::placeholders::all) << ptctrlmesh->cm.vert[vi].P()[0],
				ptctrlmesh->cm.vert[vi].P()[1], ptctrlmesh->cm.vert[vi].P()[2];
		}
	}
	else {
		ATA = AT * AT.transpose();

		solver.compute(ATA);
		controlmesh = solver.solve(projectedsplpts);
	}

}

void FilterSubdivFittingPlugin::assembleIncrement(int rank)
{
	auto ls = tri::Allocator<CMeshO>::FindPerVertexAttribute<Eigen::SparseVector<double>>(
		ptsample->cm, "LimitStencil");

	splpts.conservativeResize(ptsample->cm.VN(), 3);
	dATA.resize(ptctrlmesh->cm.VN(), ptctrlmesh->cm.VN());
	dATA.setZero();
	dps.resize(ptctrlmesh->cm.VN(),3);
	dps.setZero();
	for (int i = oldvn; i < ptsample->cm.VN(); i++) {
		splpts(i, Eigen::placeholders::all) << ptsample->cm.vert[i].P()[0],
			ptsample->cm.vert[i].P()[1], ptsample->cm.vert[i].P()[2];
		dATA += ls[i] * ls[i].transpose();
		dps += ls[i] * splpts(i, Eigen::placeholders::all);
	}

	auto R = splpts + dps;
	Eigen::SparseMatrix<double> B = solver.solve(dATA);
	Eigen::MatrixXd             rB[3];
	for (int dim = 0; dim < 3; dim++) {

	}
}

Point3f
FilterSubdivFittingPlugin::evaluateLimitPoint(int vi)
{
	auto ls = tri::Allocator<CMeshO>::FindPerVertexAttribute<Eigen::SparseVector<double>>(
		ptsample->cm, "LimitStencil");

	auto fitting_sample = ls[vi].transpose() * controlmesh;

	return Point3f(fitting_sample(0), fitting_sample(1), fitting_sample(2));
}

void FilterSubdivFittingPlugin::displayResults(const RichParameterList& par)
{
	bool show_ctrlmesh  = par.getBool("show_new_mesh");
	bool show_limitspls = par.getBool("show_limit_samples");
	bool check_bary_coord = par.getBool("check_bary_coord");

	if (show_ctrlmesh) {
		MeshModel* sourceCtrlMesh = ptctrlmesh; // source = current
		QString    newName        = sourceCtrlMesh->label() + "_fitting";
		MeshModel* destModel    = mdptr->addNewMesh(
            "",
            newName,
            true); // After Adding a mesh to a MeshDocument the new mesh is the current one
		destModel->updateDataMask(sourceCtrlMesh);
		tri::Append<CMeshO, CMeshO>::Mesh(destModel->cm, sourceCtrlMesh->cm);

		for (const std::string& tex : destModel->cm.textures) {
			destModel->addTexture(tex, sourceCtrlMesh->getTexture(tex));
		}

		log("Create control mesh solution %i", destModel->id());

		// init new layer
		destModel->updateBoxAndNormals();
		destModel->cm.Tr = sourceCtrlMesh->cm.Tr;

		for (int vi = 0; vi < destModel->cm.vert.size(); vi++) {
			destModel->cm.vert[vi].P()[0] = controlmesh(vi, 0);
			destModel->cm.vert[vi].P()[1] = controlmesh(vi, 1);
			destModel->cm.vert[vi].P()[2] = controlmesh(vi, 2);
		}
	}

	if (show_limitspls) {
		MeshModel* sourceCtrlMesh = ptsample; // source = current
		QString    newName        = sourceCtrlMesh->label() + "_fitting";
		MeshModel* destModel      = mdptr->addNewMesh(
            "",
            newName,
            true); // After Adding a mesh to a MeshDocument the new mesh is the current one
		destModel->updateDataMask(sourceCtrlMesh);
		tri::Append<CMeshO, CMeshO>::Mesh(destModel->cm, sourceCtrlMesh->cm);

		for (const std::string& tex : destModel->cm.textures) {
			destModel->addTexture(tex, sourceCtrlMesh->getTexture(tex));
		}

		log("Create fitting samples %i", destModel->id());

		// init new layer
		destModel->updateBoxAndNormals();
		destModel->cm.Tr = sourceCtrlMesh->cm.Tr;

		for (int vi = 0; vi < destModel->cm.vert.size(); vi++) {
			destModel->cm.vert[vi].P() = evaluateLimitPoint(vi);
		}

		fittingres = destModel;
	}

	if (check_bary_coord) {
		MeshModel* sourceCtrlMesh = ptsample; // source = current
		QString    newName        = sourceCtrlMesh->label() + "_barycoord";
		MeshModel* destModel      = mdptr->addNewMesh(
            "",
            newName,
            true); // After Adding a mesh to a MeshDocument the new mesh is the current one
		destModel->updateDataMask(sourceCtrlMesh);
		tri::Append<CMeshO, CMeshO>::Mesh(destModel->cm, sourceCtrlMesh->cm);

		for (const std::string& tex : destModel->cm.textures) {
			destModel->addTexture(tex, sourceCtrlMesh->getTexture(tex));
		}

		log("Create fitting samples %i", destModel->id());

		// init new layer
		destModel->updateBoxAndNormals();
		destModel->cm.Tr = sourceCtrlMesh->cm.Tr;

		CMeshO::PerVertexAttributeHandle<const CFaceO*> ftptrs;
		if (tri::HasPerVertexAttribute(ptsample->cm, "FootTriangle")) {
			ftptrs = tri::Allocator<CMeshO>::FindPerVertexAttribute<const CFaceO*>(
				ptsample->cm, "FootTriangle");
			if (!tri::Allocator<CMeshO>::IsValidHandle<const CFaceO*>(ptsample->cm, ftptrs)) {
				throw MLException("attribute already exists with a different type");
			}
		}

		CMeshO::PerVertexAttributeHandle<Point3f> ftbarycoords;
		if (tri::HasPerVertexAttribute(ptsample->cm, "BaryCoord")) {
			ftbarycoords =
				tri::Allocator<CMeshO>::FindPerVertexAttribute<Point3f>(ptsample->cm, "BaryCoord");
			if (!tri::Allocator<CMeshO>::IsValidHandle<Point3f>(ptsample->cm, ftbarycoords)) {
				throw MLException("attribute already exists with a different type");
			}
		}

		for (int vi = 0; vi < destModel->cm.vert.size(); vi++) {
			destModel->cm.vert[vi].P() = ftbarycoords[vi][0] * ftptrs[vi]->V(0)->P() +
										 ftbarycoords[vi][1] * ftptrs[vi]->V(1)->P() +
										 ftbarycoords[vi][2] * ftptrs[vi]->V(2)->P();
		}
	}
}

Eigen::VectorXd FilterSubdivFittingPlugin::weightsPatch(const CFaceO* ft, float v, float w)
{
	// return dimension of 1*NV weight vector of all control points
	// per vertex valence and per triangle subpatch pick up mat should be update

	// TODO: wrong check!
	CMeshO::PerVertexAttributeHandle<int> val;
	if (tri::HasPerVertexAttribute(ptctrlmesh->cm, "Valence")) {
		val = tri::Allocator<CMeshO>::FindPerVertexAttribute<int>(ptctrlmesh->cm, "Valence");
		if (!tri::Allocator<CMeshO>::IsValidHandle<int>(ptctrlmesh->cm, val)) {
			throw MLException("attribute already exists with a different type");
		}
	}

	CMeshO::PerFaceAttributeHandle<std::vector<Eigen::MatrixXd>> matPatchSubdiv;
	if (tri::HasPerFaceAttribute(ptctrlmesh->cm, "PatchSubdiv")) {
		matPatchSubdiv = tri::Allocator<CMeshO>::FindPerFaceAttribute<std::vector<Eigen::MatrixXd>>(
			ptctrlmesh->cm, "PatchSubdiv");
		if (!tri::Allocator<CMeshO>::IsValidHandle<std::vector<Eigen::MatrixXd>>(
				ptctrlmesh->cm, matPatchSubdiv)) {
			throw MLException("attribute already exists with a different type");
		}
	}

	if (v < 0.5 && w < 0.5 - v) {
		// patch correspond to vert 0
		return weightsIrregularPatch(val[ft->V(0)], 2 * v, 2 * w).transpose() *
			   matPatchSubdiv[ft][0];
	}
	else if (v > 0.5) {
		// patch correspond to vert 1
		return weightsIrregularPatch(val[ft->V(1)], 2 * w, 2 * (1 - v - w)).transpose() *
			   matPatchSubdiv[ft][1];
	}
	else if (w > 0.5) {
		// patch correspond to vert 2
		return weightsIrregularPatch(val[ft->V(2)], 2 * (1 - v - w), 2 * v).transpose() *
			   matPatchSubdiv[ft][2];
	}
	else {
		// patch correspond to vert 3
		return weightsIrregularPatch(6, 1 - 2 * v, 1 - 2 * w).transpose() * matPatchSubdiv[ft][3];
	}
}

Eigen::VectorXd FilterSubdivFittingPlugin::weightsIrregularPatch(int N, float v, float w)
{
	// return dimension of 1*(N+6) weight vector of control points in 1-ring of irregular patch
	if (N == 6) {
		return weightsRegularPatch(1.f - v - w, v, w);
	}

	if (v + w < eps) {
		// TODO: validate this formula
		double          alphaN_ = 3. / (11. - 8 * (0.375 + pow(0.375 + 0.25 * cos(2 * M_PI / N),
                                                      2))); // different from previous alphaN
		Eigen::VectorXd b       = Eigen::VectorXd::Zero(N + 6);
		b(Eigen::seq(0, N)).array() = (1 - alphaN_) / N;
		if (N == 6)
			b(3) = alphaN_;
		else
			b(0) = alphaN_;
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
	b /= 12.;
	return b;
}

Eigen::MatrixXd FilterSubdivFittingPlugin::matrixPickup(int N, int k)
{
	assert(N >= 3);
	auto cP = cacheP.find(4 * N + k);
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
	cacheP[4 * N + k] = P;
	return P;

}

Eigen::MatrixXd FilterSubdivFittingPlugin::matrixPatchSubdiv(int N, int n)
{
	// return dimension of (N+12)*(N+6) subdiv matrix irregular patch
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

void FilterSubdivFittingPlugin::testFVOutput(MeshModel& mm)
{
	for (auto fi = mm.cm.face.begin(); fi != mm.cm.face.end(); fi++) {
		std::cout << "face " << fi->Index() << std::endl << "vert ";
		for (int i = 0; i < 3; i++) {
			std::cout << fi->V(i)->Index() << " ";
		}
		std::cout << std::endl << std::endl;
	}
}


MESHLAB_PLUGIN_NAME_EXPORTER(FilterSubdivFittingPlugin)
