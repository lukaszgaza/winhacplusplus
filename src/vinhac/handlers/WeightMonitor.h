/*
 * \file WeightMonitor.h
 * \brief
 * \date Created on: 2009-07-08
 * \author     Author: siudmy
 */
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include <utility>
#include "../core/enum/WeightsDefinition.h"
#include "../core/enum/PolarizationDefinition.h"

#ifndef WEIGHTMONITOR_H_
#define WEIGHTMONITOR_H_
namespace VINHAC {
class WeightMonitor;
}

namespace VINHAC {

/**
 * Class responsible for cummulation of statistic for some weight.
 * It is used by MonitorHandler
 */
class WeightMonitor {
private:
	double m_nbTot; //!<  Total number of weights
	double m_wtSum; //!<  Sum of weights
	double m_wt2Sum; //!<  Sum of weights^2

	double m_WtMax; //!<  Max weight
	double m_WtMin; //!<  Min weight

public:
	WeightMonitor();
	WeightMonitor(const double& wtMax); //!< set max weight in the constructor
	~WeightMonitor() {}

	//! add weight to the monitor and calculate all important quantities
	bool addWeight(const double& weight);

	//! get average weight
	inline double averageWeight() const{
		return averageWeight(m_nbTot);
	}
	//! get average weight which overrides stored number of entries
	inline double averageWeight(double m_nbTot) const{
		return m_wtSum / m_nbTot;
	}
	//! get error
	inline double error() const{
		return error(m_nbTot);
	}

	//! get error which overrides stored number of entries
	inline double error(double m_nbTot) const{
		if(m_wtSum!=0.0){
			return sqrt(fabs((m_wt2Sum / m_wtSum) / m_wtSum - 1.0 / m_nbTot));
		} else {
			return 0.0;
		}
	}

	inline const double& getSumWt2() const {
		return m_wt2Sum;
	}
	std::string makeSnapshot(std::string monitorName);
	std::string restoreState(const std::vector<std::string>& state);

	// getters and setters


	inline const double& getWtMax() const {
		return m_WtMax;
	}
	inline const double& getWtMin() const {
		return m_WtMin;
	}
};
}

#endif /* WEIGHTMONITOR_H_ */
