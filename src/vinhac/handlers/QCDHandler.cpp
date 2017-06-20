/*
 * QCDHandler.cpp
 *
 *  Created on: 03-08-2011
 *      Author: sobol
 */

#include "QCDHandler.h"
#include "../core/VinhacException.h"
#include "../util/Transformation.h"
#include <iostream>
#include <cstdlib>
using namespace std;

VINHAC::QCDHandler::QCDHandler(TRND *p_randomGenerator, Manager& manager,
		std::string inputPath, std::string pythiaXMLdocdir) :
	pythia(pythiaXMLdocdir), randomEngine(0), manager(manager),
			winhacLock(new boost::mutex()), pythiaLock(new boost::mutex()),
			stopFlag(false), pythiaThread(0), thread(0) {
	try {
		randomEngine = new RandomEngineAdapter(*p_randomGenerator);
		pythia.setRndmEnginePtr(randomEngine);
		updatePythiaSettings(pythia, DataSource::get());

		fifoname = tmpnam(0);
		if (system(("mkfifo " + fifoname).c_str())) {
			throw VinhacException("Unable to create FIFO PIPE: " + fifoname);
		}

		pythiaThread = new PythiaThread(pythia, *pythiaLock, *winhacLock,
				stopFlag, fifoname);

		thread = new boost::thread(*pythiaThread);

		lhefile.open(fifoname.c_str(), fstream::out);
		lheWriter = new LHEWriter(manager, lhefile);
		lheWriter->writeOpening();
		lheWriter->writeInit();
		winhacLock->lock();
	} catch (...) {
		throw VinhacException("Initialization of Pythia failed");
	}
}

VINHAC::QCDHandler::~QCDHandler() {

	lheWriter->writeClosing();
	lhefile.close();

	if (lheWriter) {
		delete lheWriter;
		lheWriter = 0;
	}

	if (thread) {
		stopFlag = true;
		pythiaLock->unlock();
		thread->join();
		thread = 0;
	}
	if (pythiaThread) {
		pythiaThread = 0;
	}
	if (randomEngine) {
		delete randomEngine;
		randomEngine = 0;
	}
	system(("rm -f " + fifoname).c_str());

	if (winhacLock) {
		winhacLock = 0;
	}
	if (pythiaLock) {
		pythiaLock = 0;
	}
}

bool VINHAC::QCDHandler::eventEvolve(VINHAC::Event& event) {

	lheWriter->writeEvent(event);
	pythiaLock->unlock();
	winhacLock->lock();
	//cout << endl;
	//cout << event.printInLabFrame() << endl;
	updateEventMomenta(event, pythia.event);
	//cout << event.printInLabFrame() << endl;

	return true;
}

void VINHAC::QCDHandler::updateEventMomenta(VINHAC::Event& event,
		Pythia8::Event &pythiaEvent) {
	// Extracting data from Event object
	Pythia8::Particle quark;
	Pythia8::Particle antiQuark;
	Pythia8::Particle boson;
	Pythia8::Particle lepton;
	Pythia8::Particle neutrino;
	std::vector < Pythia8::Particle > photons;

	Pythia8::Vec4 pq;
	Pythia8::Vec4 paq;

	for (int i = 0; i < pythiaEvent.size(); ++i) {
		int id = pythiaEvent[i].id();
		int absid = abs(pythiaEvent[i].id());
		if (absid > 0 && absid <= 6) {
			if (id > 0) {
				quark = pythiaEvent[pythiaEvent.iBotCopyId(i)];
				set<int> sisters;
				getAllSisters(sisters, pythiaEvent, i);
				for (set<int>::iterator it = sisters.begin(); it
						!= sisters.end(); ++it) {
					pq += pythiaEvent[*it].p();
				}
				break;
			}
		}
	}

	for (int i = 0; i < pythiaEvent.size(); ++i) {
		int id = pythiaEvent[i].id();
		int absid = abs(pythiaEvent[i].id());
		if (absid > 0 && absid <= 6) {
			if (id < 0) {
				antiQuark = pythiaEvent[pythiaEvent.iBotCopyId(i)];
				set<int> sisters;
				getAllSisters(sisters, pythiaEvent, i);
				for (set<int>::iterator it = sisters.begin(); it
						!= sisters.end(); ++it) {
					paq += pythiaEvent[*it].p();
				}
				break;
			}
		}
	}

	for (int i = 0; i < pythiaEvent.size(); ++i) {
		int absid = abs(pythiaEvent[i].id());

		if (absid == 24) {
			boson = pythiaEvent[pythiaEvent.iBotCopyId(i)];
			break;
		}
	}

	for (int i = 0; i < pythiaEvent.size(); ++i) {
		int absid = abs(pythiaEvent[i].id());

		if (absid == 11 || absid == 13 || absid == 15) {
			lepton = pythiaEvent[pythiaEvent.iBotCopyId(i)];
			break;
		}
	}

	for (int i = 0; i < pythiaEvent.size(); ++i) {
		int absid = abs(pythiaEvent[i].id());

		if (absid == 12 || absid == 14 || absid == 16) {
			neutrino = pythiaEvent[pythiaEvent.iBotCopyId(i)];
			break;
		}
	}

	for (int i = 0; i < pythiaEvent.size(); ++i) {
		int absid = abs(pythiaEvent[i].id());
		if (absid == 22) {
			photons.push_back(pythiaEvent[i]);
		}
	}

	updateQuarksMomentum(event.getQuark().getFourMomentum(),
			event.getAntiQuark().getFourMomentum(), quark.p() - pq,
			antiQuark.p() - paq, boson.p());
	updateMomentum(event.getBoson().getFourMomentum(), boson.p());
	updateMomentum(event.getLepton().getFourMomentum(), lepton.p());
	updateMomentum(event.getNeutrino().getFourMomentum(), neutrino.p());

	event.getPhotons().clear();
	event.getPhotons().resize(photons.size());
	for (unsigned i = 0; i < photons.size(); ++i) {
		updateMomentum(event.getPhotons()[i].getFourMomentum(), photons[i].p());
		int motherIdx = 0;
		if (photons[i].mother1() != 0) {
			motherIdx = photons[i].mother1();
		} else if (photons[i].mother2() != 0) {
			motherIdx = photons[i].mother2();
		}
		if (motherIdx > 0) {
			Particle* parent = 0;
			int motherId = pythiaEvent[motherIdx].id();
			int motherAbsid = pythiaEvent[motherIdx].idAbs();
			if (motherAbsid > 0 && motherAbsid <= 6) {
				if (motherId > 0) {
					parent = &event.getQuark();
				} else {
					parent = &event.getAntiQuark();
				}
			}
			if (motherAbsid == 24)
				parent = &event.getBoson();
			if (motherAbsid == 11 || motherAbsid == 13 || motherAbsid == 15)
				parent = &event.getLepton();
			if (motherAbsid == 12 || motherAbsid == 14 || motherAbsid == 16)
				parent = &event.getNeutrino();
			if (parent != 0) {
				event.getPhotons()[i].addParent(parent);
			}
		}
	}

	//sort photons in the descending order in energy
	sort(event.getPhotons().begin(), event.getPhotons().end(), photonComparator);

}

void VINHAC::QCDHandler::getAllSisters(set<int>& target,
		Pythia8::Event &pythiaEvent, int i) {
	vector<int> sisters = pythiaEvent.sisterListTopBot(i);
	for (unsigned k = 0; k < sisters.size(); ++k) {
		target.insert(sisters[k]);
	}

	vector<int> mothers = pythiaEvent.motherList(i);
	for (unsigned k = 0; k < mothers.size(); ++k) {
		if (mothers[k] > 2) {
			getAllSisters(target, pythiaEvent, mothers[k]);
		}
	}
}

void VINHAC::QCDHandler::updateMomentum(CLHEP::HepLorentzVector &target,
		Pythia8::Vec4 source) {
	target.set(source.px(), source.py(), source.pz(), source.e());
}

void VINHAC::QCDHandler::updateQuarksMomentum(
		CLHEP::HepLorentzVector &targetQuark,
		CLHEP::HepLorentzVector &targetAntiQuark, Pythia8::Vec4 sourceQuark,
		Pythia8::Vec4 sourceAntiQuark, Pythia8::Vec4 sourceBoson) {

	double qT2 = sourceQuark.px()*sourceQuark.px() + sourceQuark.py()*sourceQuark.py();
	double aqT2 = sourceAntiQuark.px()*sourceAntiQuark.px() + sourceAntiQuark.py()*sourceAntiQuark.py();
	double mq2 = targetQuark.m2();
	double maq2 = targetAntiQuark.m2();
	double Q4 = sourceBoson.e();
	double Q3 = sourceBoson.pz();

	double C = (aqT2-qT2+maq2-mq2+Q4*Q4-Q3*Q3)/2.0/Q4;
	double D = Q3/Q4;

	double delta = 4*D*D*C*C-4*(D*D-1.0)*(C*C-aqT2-maq2);

	double aq3 = 0.0;
	double aq4 = 0.0;
	double q3 = 0.0;
	double q4 = 0.0;

	double aq3_2 = 0.0;
	double aq4_2 = 0.0;
	double q3_2 = 0.0;
	double q4_2 = 0.0;

	aq3 = (-2.0*D*C+sqrt(delta))/2.0/(D*D-1);
	aq4 = aq3*D+C;
	q3 = Q3 - aq3;
	q4 = Q4 - aq4;
	aq3_2 = (-2.0*D*C-sqrt(delta))/2.0/(D*D-1);
	aq4_2 = aq3_2*D+C;
	q3_2 = Q3 - aq3_2;
	q4_2 = Q4 - aq4_2;

	if(fabs(targetQuark.pz()) > fabs(targetAntiQuark.pz())){
		if(fabs(q3) > fabs(aq3)) {
			targetQuark.set(sourceQuark.px(), sourceQuark.py(), q3, q4);
			targetAntiQuark.set(sourceAntiQuark.px(), sourceAntiQuark.py(), aq3, aq4);
		} else {
			targetQuark.set(sourceQuark.px(), sourceQuark.py(), q3_2, q4_2);
			targetAntiQuark.set(sourceAntiQuark.px(), sourceAntiQuark.py(), aq3_2, aq4_2);
		}
	} else {
		if(fabs(q3) < fabs(aq3)) {
			targetQuark.set(sourceQuark.px(), sourceQuark.py(), q3, q4);
			targetAntiQuark.set(sourceAntiQuark.px(), sourceAntiQuark.py(), aq3, aq4);
		} else {
			targetQuark.set(sourceQuark.px(), sourceQuark.py(), q3_2, q4_2);
			targetAntiQuark.set(sourceAntiQuark.px(), sourceAntiQuark.py(), aq3_2, aq4_2);
		}
	}
}

void VINHAC::QCDHandler::updatePythiaSettings(Pythia8::Pythia& pythia,
		DataSource& ds) {
	pythia.readString("PhaseSpace:mHatMin = 0.0");
	pythia.readString("ProcessLevel:resonanceDecays = off");
	pythia.readString("PartonLevel:MI = off");
	pythia.readString("PartonLevel:FSR = off");
	if (ds.isPartonShower) {
		pythia.readString("PartonLevel:ISR = on");
	} else {
		pythia.readString("PartonLevel:ISR = off");
	}
	pythia.readString("HadronLevel:Decay = off");
	pythia.readString("HadronLevel:BoseEinstein = off");
	if (ds.isHadronization) {
		pythia.readString("HadronLevel:Hadronize = on");
		pythia.readString("PartonLevel:Remnants = on");
	} else {
		pythia.readString("HadronLevel:Hadronize = off");
		pythia.readString("PartonLevel:Remnants = off");
	}
	{
		stringstream ss;
		ss << "StandardModel:alphaEM0 = " << ds.alphaQED;
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "SigmaProcess:alphaSvalue = " << ds.alphaQCD;
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:sin2thetaW = " << ds.sinThetaW2;
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vcb = " << ds.ckm[1][2];
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vcd = " << ds.ckm[1][0];
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vcs = " << ds.ckm[1][1];
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vtb = " << ds.ckm[2][2];
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vtd = " << ds.ckm[2][0];
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vts = " << ds.ckm[2][1];
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vub " << ds.ckm[0][2];
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vud = " << ds.ckm[0][0];
		pythia.readString(ss.str());
	}
	{
		stringstream ss;
		ss << "StandardModel:Vus = " << ds.ckm[0][1];
		pythia.readString(ss.str());
	}

}
