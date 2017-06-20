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
	std::string m_monitorName;//!<  Name of monitor handler = name of weight

	double m_nbTot; //!<  Total number of weights
	double m_wtSum; //!<  Sum of weights
	double m_wt2Sum; //!<  Sum of weights^2

	double m_WtMax; //!<  Max weight
	double m_WtMin; //!<  Min weight

	double dmax(double x, double y) {
		if (x > y)
			return x;
		else
			return y;
	}
	double dmin(double x, double y) {
		if (x > y)
			return y;
		else
			return x;
	}

public:
	WeightMonitor() {
	}
	;
	WeightMonitor(std::string name);
	WeightMonitor(std::string name, double wtMax); //!< set max weight in the constructor
	~WeightMonitor() {
	}

	//! add weight to the monitor and calculate all important quantities
	bool addWeight(double weight);

	//! get average weight
	inline double averageWeight() {
		return m_wtSum / m_nbTot;
	}
	//! get error
	inline double error() {
		return sqrt(fabs((m_wt2Sum / m_wtSum) / m_wtSum - 1.0 / m_nbTot));
	}

	inline double getSumWt2() {
		return m_wt2Sum;
	}
	std::string makeSnapshot();
	void restoreState(std::vector<std::string> state);

	// getters and setters


	inline double getWtMax() {
		return m_WtMax;
	}
	inline double getWtMin() {
		return m_WtMin;
	}
};
}

#endif /* WEIGHTMONITOR_H_ */
