/*
 * \file MonitorHandler.h
 * \brief
 * \date Created on: 2009-07-07
 * \author     Author: siudmy
 */

#ifndef MONITORHANDLER_H_
#define MONITORHANDLER_H_
#include "ProcessHandler.h"
#include "WeightMonitor.h"
namespace VINHAC {
class MonitorHandler;
class ProcessHandler;
class WeightMonitor;
}

namespace VINHAC {

//!  Class responsible for cumulation of weights
/*!
  This class cumulates weights statistics using WeightMonitor objects.
  It contains common (aka crude) weights monitors, model weights monitors, main weight monitor (calculated using ModelBest weight).
  It also contains cross sections (various combinantion of common and model weights) monitors.
*/
class MonitorHandler: public ProcessHandler {
private:
	int mode; //! mode of the Handler 0 - need to be initialise, 1 - collecting wt

	map<string, VINHAC::WeightMonitor> m_commonWtMonitors;
	map<string, VINHAC::WeightMonitor> m_modelWtMonitors;
	map<string, VINHAC::WeightMonitor> xsYFSMonitors;
	map<string, VINHAC::WeightMonitor> xsFixedOrderMonitors;
	map<string, VINHAC::WeightMonitor> xsAllMonitors;
	WeightMonitor mainWeightMonitor;
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
	MonitorHandler();

	/**
	* \brief evolves the event
	*
	* It does whole part of generation process which is to do in this handler
	* @see Event
	* @param event reference to current Event object
	* @return if weight was zero
	*/
	bool eventEvolve(VINHAC::Event& event);

	/**
	 * \brief saves monitors state when generator makes snapshot
	 *
	 * @param snapshotPath path to folder where snapshot is stored
	 */
	void makeSnapshot(string snapshotPath);

	/**
	 * \brief restores monitors state from snapshot
	 *
	 * @param snapshotPath path to folder where snapshot is stored
	 */
	void restoreState(string snapshotPath);

	~MonitorHandler();

	/**
	 * \brief prints average weigths
	 */
	void printAverageWeights();

	/**
	 * \brief prints cross sections
	 */
	void printCrossSections();
	// Getters and setters

	/**
	 * \brief returns map containing common weights monitors
	 *
	 * @return map of common weights monitors
	 */
	inline map<string, VINHAC::WeightMonitor> getCommonWtMonitors() {
		return m_commonWtMonitors;
	}

	/**
	 * \brief returns map containing model weights monitors
	 *
	 * @return map of model weights monitors
	 */
	inline map<string, VINHAC::WeightMonitor> getModelWtMonitors() {
		return m_modelWtMonitors;
	}

	/**
	 * \brief returns map containing cross sections YFS monitors
	 *
	 * @return map of cross sections monitors
	 */
	inline map<string, VINHAC::WeightMonitor> getCrossSectionsYFSMonitors() {
		return xsYFSMonitors;
	}

	/**
	 * \brief returns map containing cross sections fixed order monitors
	 *
	 * @return map of cross sections monitors
	 */
	inline map<string, VINHAC::WeightMonitor> getCrossSectionsFixedOrderMonitors() {
		return xsFixedOrderMonitors;
	}

	/**
	 * \brief returns map containing all cross sections fixed order monitors
	 *
	 * @return map of cross sections monitors
	 */
	inline map<string, VINHAC::WeightMonitor> getAllCrossSectionsMonitors() {
		return xsAllMonitors;
	}


	/**
	 * \brief returns main weight monitor
	 *
	 * @return main weight monitor
	 */
	inline WeightMonitor getMainWtMonitor() {
		return mainWeightMonitor;
	}
};

}
#endif /* MONITORHANDLER_H_ */
