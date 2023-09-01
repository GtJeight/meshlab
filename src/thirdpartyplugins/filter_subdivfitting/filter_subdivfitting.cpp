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
		//for (size_t iv = 0; iv < m.vn; iv++)
		//	log("%d: %f, %f, %f\n", iv, m.vert[iv].P()[0], m.vert[iv].P()[1], m.vert[iv].P()[2]);
		for (size_t i = 0; i < m.fn; i++)
			log("%d: %f, %f, %f\n", i, m.face[i].N()[0], m.face[i].N()[1], m.face[i].N()[2]);
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
	switch (mode) {
	case FootPointMode::MODE_MESH: {

	} break;
	case FootPointMode::MODE_SUBDIVISION: {

	} break;
	default: break;
	}
}

std::pair<CFaceO*, vcg::Point3f>
FilterSubdivFittingPlugin::distancePointTriangle(const CVertexO& v, const CFaceO& f)
{
	auto n  = f.cN();
	auto p  = v.cP();

	return std::pair<CFaceO*, vcg::Point3f>(0, vcg::Point3f(0, 0, 0));
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterSamplePlugin)
