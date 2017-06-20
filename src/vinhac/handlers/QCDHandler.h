/*
 * QCDHandler.h
 *
 *  Created on: 03-08-2011
 *      Author: sobol
 */

#ifndef QCDHANDLER_H_
#define QCDHANDLER_H_

#include "ProcessHandler.h"
#include "Pythia.h"
#include "../input/DataSource.h"
#include "../../utils/FOAM/TRND.h"
#include "../../utils/pythia/RandomEngineAdapter.h"
#include "../../utils/pythia/PythiaThread.h"
#include "../interfaces/LHEWriter.h"
#include "../../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include <boost/thread/thread.hpp>

namespace VINHAC {

class Event;

struct PhotonComarator {
	bool operator()(const Particle& photon1, const Particle& photon2) {
		return photon1.compareTo(photon2) == 1;
	}
};

class QCDHandler: public VINHAC::ProcessHandler {
private:
	Pythia8::Pythia pythia;
	RandomEngineAdapter *randomEngine;
	string fifoname;
	Manager& manager;
	LHEWriter *lheWriter;
	fstream lhefile;
	boost::mutex *winhacLock;
	boost::mutex *pythiaLock;
	bool stopFlag;
	PythiaThread *pythiaThread;
	boost::thread *thread;
	PhotonComarator photonComparator;

	void getAllSisters(set<int>& target, Pythia8::Event &pythiaEvent, int i);
	void updatePythiaSettings(Pythia8::Pythia& pythia, DataSource& ds);
	void updateQuarksMomentum(
			CLHEP::HepLorentzVector &targetQuark,
			CLHEP::HepLorentzVector &targetAntiQuark, Pythia8::Vec4 sourceQuark,
			Pythia8::Vec4 sourceAntiQuark, Pythia8::Vec4 sourceBoson);
	void updateEventMomenta(Event& event, Pythia8::Event &pythiaEvent);
	void updateMomentum(CLHEP::HepLorentzVector &target, Pythia8::Vec4 source);
	inline bool comparePhotons(Particle& photon1, Particle& photon2) {
		return (abs(photon1.getFourMomentum().e()) < abs(
				photon2.getFourMomentum().e()));
	}
public:
	QCDHandler(TRND *p_randomGenerator, Manager& manager,
			std::string inputPath, std::string pythiaXMLdocdir);
	virtual ~QCDHandler();

	/**
	 * \brief evolves the event
	 *
	 * It does whole part of generation process which is to do in this handler
	 * @see Event
	 * @param event reference to current Event object
	 * @return if weight was zero
	 */
	bool eventEvolve(VINHAC::Event& event);

};

}

#endif /* QCDHANDLER_H_ */
