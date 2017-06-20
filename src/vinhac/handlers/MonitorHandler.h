/*
 * \file MonitorHandler.h
 * \brief
 * \date Created on: 2009-07-07
 * \author     Author: siudmy
 */

#ifndef MONITORHANDLER_H_
#define MONITORHANDLER_H_

#include <map>

#include "WeightMonitor.h"
#include "ProcessHandler.h"
#include "../core/enum/WeightsDefinition.h"

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
	double numberOfEvents;

	map<WeightName::name, VINHAC::WeightMonitor> m_commonWtMonitors;
	map<pair<WeightName::name,PolarizationName::name>, VINHAC::WeightMonitor> m_modelWtMonitors;
	map<pair<WeightName::name,PolarizationName::name>, VINHAC::WeightMonitor> xsMonitors;
	WeightMonitor mainWeightMonitor;

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
	void printAverageWeights(ostream& out = cout);

	/**
	 * \brief prints cross sections
	 */
	void printCrossSections(ostream& out = cout);
	// Getters and setters

	/**
	 * \brief returns map containing common weights monitors
	 *
	 * @return map of common weights monitors
	 */
	inline const map<WeightName::name, VINHAC::WeightMonitor>& getCommonWtMonitors() const {
		return m_commonWtMonitors;
	}

	/**
	 * \brief returns map containing model weights monitors
	 *
	 * @return map of model weights monitors
	 */
	inline const map<pair<WeightName::name,PolarizationName::name>, VINHAC::WeightMonitor>& getModelWtMonitors() const {
		return m_modelWtMonitors;
	}

	/**
	 * \brief returns main weight monitor
	 *
	 * @return main weight monitor
	 */
	inline const WeightMonitor& getMainWtMonitor()const {
		return mainWeightMonitor;
	}

	/**
	 * \brief returns number of registered events
	 *
	 * @return number of registered events
	 */
	inline double getNumberOfEvents(){
		return numberOfEvents;
	}
};

}
#endif /* MONITORHANDLER_H_ */
