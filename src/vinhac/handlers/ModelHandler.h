/*
 * ModelHandler.h
 *
 *  Created on: 26 Aug 2009
 *      Author: siudmy
 */

#ifndef MODELHANDLER_H_
#define MODELHANDLER_H_
#include "../../utils/SpecialMatrix.h"
#include "../../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include <cstdlib>
#include <map>
#include "ProcessHandler.h"
using namespace std;

namespace VINHAC {
class Manager;
class Event;
class EventStep;
class Particle;
class ModelHandler;
class HardProcessHandler;
class ProcessHandler;
class DataSource;
}

namespace VINHAC {

//!  Class responsible for model weights calculation
/*!
 This class calculates various model weights using proper matrix elements
 */
class ModelHandler: public ProcessHandler {
private:

	VINHAC::Manager* p_Manager;
	VINHAC::HardProcessHandler* p_HardMatrix;

	double ckmElement;
	const double normalizationFactor0;
	const double normalizationFactor1;
	const double alfpi;

	map<int, double> deltaWeak;

	CLHEP::HepLorentzVector quark;
	CLHEP::HepLorentzVector antiquark;
	CLHEP::HepLorentzVector boson;
	CLHEP::HepLorentzVector lepton;
	CLHEP::HepLorentzVector antilepton;
	int chargeOfBoson;
	int overWeightedEvents;

	void prepareParticles(VINHAC::Event& event);

	double differentialXsectionCrude(const double& s, Event&);
	void printHighWeightWaring(Event& event, string name,
			double beta00, double DsiCru, double beta01_EW, double beta01_QED,
			double beta11, double beta01_LL, double WtCrud, double FacWcr);

public:


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
	 * \brief constructor
	 *
	 * @param ds pointer to structure with input data
	 * @param p_Hard pointer to HadrProcessHandler object
	 */
	ModelHandler( VINHAC::HardProcessHandler* p_Hard);
	~ModelHandler() {
	}

};
}

#endif /* MODELHANDLER_H_ */
