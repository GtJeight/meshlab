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

#include "filter_example.h"
#include <fstream>

/**
 * @brief
 * Constructor usually performs only two simple tasks of filling the two lists
 *  - typeList: with all the possible id of the filtering actions
 *  - actionList with the corresponding actions.
 * If you want to add icons to your filtering actions you can do here by construction the QActions accordingly
 */
FilterExamplePlugin::FilterExamplePlugin()
{ 
	typeList = { FP_MOVE_VERTEX};

	for(const ActionIDType& tt : typeList)
		actionList.push_back(new QAction(filterName(tt), this));
}

QString FilterExamplePlugin::pluginName() const
{
	return "FilterExample";
}

QString FilterExamplePlugin::vendor() const
{
	return "CNR-ISTI-VCLab";
}

/**
 * @brief ST() must return the very short string describing each filtering action
 * (this string is used also to define the menu entry)
 * @param filterId: the id of the filter
 * @return the name of the filter
 */
QString FilterExamplePlugin::filterName(ActionIDType filterId) const
{
	switch(filterId) {
	case FP_MOVE_VERTEX :
		return "Selected";
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
QString FilterExamplePlugin::filterInfo(ActionIDType filterId) const
{
	switch(filterId) {
	case FP_MOVE_VERTEX :
		return "Move the vertices of the mesh of a random quantity.";
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
FilterExamplePlugin::FilterClass FilterExamplePlugin::getClass(const QAction *a) const
{
	switch(ID(a)) {
	case FP_MOVE_VERTEX :
		return FilterPlugin::Smoothing;
	default :
		assert(0);
		return FilterPlugin::Generic;
	}
}

/**
 * @brief FilterSamplePlugin::filterArity
 * @return
 */
FilterPlugin::FilterArity FilterExamplePlugin::filterArity(const QAction*) const
{
	return SINGLE_MESH;
}

/**
 * @brief FilterSamplePlugin::getPreConditions
 * @return
 */
int FilterExamplePlugin::getPreConditions(const QAction*) const
{
	return MeshModel::MM_NONE;
}

int FilterExamplePlugin::getRequirements(const QAction* action)
{
	switch (ID(action)) {
	case FP_MOVE_VERTEX: return MeshModel::MM_FACECOLOR | MeshModel::MM_FACEQUALITY;
	default: assert(0);
	}
	return 0;
}

/**
 * @brief FilterSamplePlugin::postCondition
 * @return
 */
int FilterExamplePlugin::postCondition(const QAction*) const
{
	return MeshModel::MM_FACECOLOR + MeshModel::MM_FACEQUALITY;
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
RichParameterList FilterExamplePlugin::initParameterList(const QAction *action, const MeshModel &m)
{
	RichParameterList parlst;
	switch(ID(action)) {
	case FP_MOVE_VERTEX :
		parlst.addParam(RichInt(
			"aint",
			1,
			"Noise bits:",
			"Bits of noise added to each RGB channel. Example: 3 noise bits adds three random "
			"offsets in the [-4,+4] interval to each RGB channels."));

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
std::map<std::string, QVariant> FilterExamplePlugin::applyFilter(
		const QAction * action,
		const RichParameterList & par,
		MeshDocument &md,
		unsigned int& /*postConditionMask*/,
		vcg::CallBackPos *cb)
{
	switch(ID(action)) {
	case FP_MOVE_VERTEX : vertexDisplacement(md, cb, par.getInt("aint"));
		break;
	default :
		wrongActionCalled(action);
	}
	return std::map<std::string, QVariant>();
}

void FilterExamplePlugin::vertexDisplacement(
		MeshDocument &md,
		vcg::CallBackPos *cb, int i)
{
	CMeshO &m = md.mm()->cm;
	std::string filename;
	if (i == 1) {
		filename = "D:/source/udoboolean/zxdata/outputface1.txt";
	}
	else {
		filename = "D:/source/udoboolean/zxdata/outputface2.txt";
	}

	std::ifstream faces;
	faces.open(filename);
	size_t fi;
	while (!faces.eof()) {
		faces >> fi;
		m.face[fi].Q() = 10.;
	}
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterSamplePlugin)
