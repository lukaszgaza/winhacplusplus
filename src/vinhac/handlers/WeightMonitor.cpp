/*
 * \file WeightMonitor.cpp
 * \brief
 * \date Created on: 2009-07-08
 * \author     Author: siudmy
 */

#include "WeightMonitor.h"
#include <sstream>
#include <cstdlib>
#include <utility>
#include "../util/Commons.h"

VINHAC::WeightMonitor::WeightMonitor() :
	m_nbTot(0), m_wtSum(0), m_wt2Sum(0), m_WtMax(-1e10), m_WtMin(1e10) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: WeightMonitor::WeightMonitor() 1"<<std::endl;

#endif

}

VINHAC::WeightMonitor::WeightMonitor(const double& wtMax) :
	m_nbTot(0), m_wtSum(0), m_wt2Sum(0), m_WtMax(wtMax), m_WtMin(-wtMax) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: WeightMonitor::WeightMonitor() 2"<<std::endl;

#endif
}

bool VINHAC::WeightMonitor::addWeight(const double& MCwt) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: WeightMonitor::addWeight()"<<std::endl;

#endif
	// accumulation of statistics for the main MC weight
	m_wtSum += MCwt; // sum of Wt
	m_wt2Sum += MCwt * MCwt; // sum of Wt**2
	m_nbTot++; // sum of 1d0
	m_WtMax = std::max(m_WtMax, MCwt); // maximum wt
	m_WtMin = std::min(m_WtMin, MCwt); // minimum wt
	return true;
}

std::string VINHAC::WeightMonitor::makeSnapshot(
		std::string m_monitorName) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: WeightMonitor::makeSnapshot()"<<std::endl;

#endif
	std::stringstream ss;
	ss.precision(20);
	ss << m_monitorName << ":";
	ss.precision(20);
	ss << m_nbTot << ":";
	ss.precision(20);
	ss << m_wtSum << ":";
	ss.precision(20);
	ss << m_wt2Sum << ":";
	ss.precision(20);
	ss << m_WtMax << ":";
	ss.precision(20);
	ss << m_WtMin;
	return ss.str();
}

std::string VINHAC::WeightMonitor::restoreState(
		const std::vector<std::string>& state) {
	std::string m_monitorName;
	if (state.size() >= 1)
		m_monitorName = state[0];
	if (state.size() >= 2)
		m_nbTot = atof(state[1].c_str());
	if (state.size() >= 3)
		m_wtSum = atof(state[2].c_str());
	if (state.size() >= 4)
		m_wt2Sum = atof(state[3].c_str());
	if (state.size() >= 5)
		m_WtMax = atof(state[4].c_str());
	if (state.size() >= 6)
		m_WtMin = atof(state[5].c_str());

	return m_monitorName;
}
