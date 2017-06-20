#include <string>
#include <vector>
#include <exception>
#include <ostream>
#include <fstream>
using namespace std;

#include "Manager.h"
#include "VinhacException.h"
#include "../handlers/InitializationHandler.h"
#include "../input/DataSource.h"
#include "Event.h"
#include "../handlers/BeamHandler.h"
#include "../handlers/MonitorHandler.h"
#include "../handlers/HardProcessHandler.h"
#include "../handlers/ModelHandler.h"
#include "../handlers/QCDHandler.h"
#include "enum/PolarizationDefinition.h"
#include "../util/Commons.h"
//// Constructors
VINHAC::Manager::Manager(const string& aInputStream) :
	m_userFilePath(aInputStream), p_InitializationHandler(0),
			p_randomGenerator(0) {
#ifdef DEBUG_FLAGw

	std::cout<<"DEBUG: Manager::Manag`er() 1"<<std::endl;

#endif
	std::cout << ">>>>>>>> Reading user file: " << m_userFilePath << endl;
	m_snapshotNb = 50;
}

VINHAC::Manager::Manager() :
	m_userFilePath(), p_InitializationHandler(0), p_randomGenerator(0) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::Manager() 2"<<std::endl;

#endif
	m_userFilePath.clear();
	m_snapshotNb = 1;
}

VINHAC::Manager::~Manager() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::~Manager()"<<std::endl;

#endif

	if (p_randomGenerator) {
		delete p_randomGenerator;
		p_randomGenerator = 0;
	}

	if (p_InitializationHandler) {
		delete p_InitializationHandler;
		p_InitializationHandler = 0;
	}
}
void VINHAC::Manager::initializeRandomNumberGenerator() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::InitializeRandomNumberGenerator()"<<std::endl;

#endif

	//// read seed from the seed file
	string datfile_comp_path = input_Path + "/seed.dat";
	ifstream seed_file(datfile_comp_path.c_str());

	if (seed_file.is_open()) {
		int N(0);
		long seed;
		while (seed_file >> seed) {
			N++;
		}

		seed_file.close();
		//// end of reading a seed file

		if (N != 1) { // should be one and only one seed in the file
			cout << "There is somthing wrong with the seed from file: "
					<< datfile_comp_path << endl;
			throw VinhacException(
					"There is somthing wrong with the seed from file: "
							+ datfile_comp_path);
		} else {
			if (seed == 0) { // if seed = 0 then use time for seed
				seed = time(0);
			}
			cout << "$$$$$$$$$$ SEED = " << seed << endl;
			p_randomGenerator = new CLHEP::RanluxEngine(seed, 3);
			return;
		}
	} else {
		cout << "Unable to open seed file!!" << endl;
		throw VinhacException("Unable to open seed file!!");
	}

}

///// Z VINHAC'a
void VINHAC::Manager::initializeGenerator(int argc, char* argv[]) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::initializeGenerator()"<<std::endl;

#endif

	p_InitializationHandler = new InitializationHandler();

	p_InitializationHandler->setUserFilePath(m_userFilePath);
	p_InitializationHandler->setInputPath(input_Path);
	p_InitializationHandler->setOutputPath(output_Path);

	p_InitializationHandler->setManager(this);

	initializeRandomNumberGenerator();

	//// Initialise handlers
	p_InitializationHandler->initializeAllHandlers();

	//// add them to the handler vector. Becareful have to be added in the proper order
	_handlers.push_back(p_InitializationHandler->getBeamHandler());
	_handlers.push_back(p_InitializationHandler->getHardProcessHandler());
	if (DataSource::get().keyPol > 0 && (DataSource::get().isPartonShower
			|| DataSource::get().isHadronization)) {
		if (p_InitializationHandler->getQCDHandler() == 0) {
			throw VinhacException(
					"The polarization and QCD effects are set but QCD Handler is not present ...");
		} else {
			_handlers.push_back(p_InitializationHandler->getQCDHandler());
		}
	}
	_handlers.push_back(
			p_InitializationHandler->getFrameTransformationHandler());
	_handlers.push_back(p_InitializationHandler->getModelHandler());
	//_handlers.push_back(p_InitializationHandler->getAnalisysHandler());
	_handlers.push_back(p_InitializationHandler->getMonitorHandler());
	if (DataSource::get().keyPol == 0 && (DataSource::get().isPartonShower
			|| DataSource::get().isHadronization)) {
		if (p_InitializationHandler->getQCDHandler() == 0) {
			throw VinhacException(
					"The polarization is off and QCD effects are set but QCD Handler is not present ...");
		} else {
			_handlers.push_back(p_InitializationHandler->getQCDHandler());
		}
	}

	p_MonitorHandler = p_InitializationHandler->getMonitorHandler();
	p_BeamHandler = p_InitializationHandler->getBeamHandler();

	fstream crudefile((output_Path + "/XSectionCrude.dat").c_str(),
			fstream::out);
	crudefile << p_BeamHandler->getXsCrud() << ":"
			<< p_BeamHandler->getXsCrudError();
	crudefile.close();

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::initializeGenerator() # handlers.size(): "<<_handlers.size()<<std::endl;

#endif

}

void VINHAC::Manager::printSummary(ostream& out) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::printSummary()"<<std::endl;

#endif
	this->p_MonitorHandler->printAverageWeights(out);
	this->p_MonitorHandler->printCrossSections(out);
}

//END

bool VINHAC::Manager::checkSemafor(string semaforFile, string snapshotPath,
		double& i, double& nbEvents) {
	string semafor;
	fstream fs(semaforFile.c_str(), ifstream::in);
	fs >> semafor;
	fs.close();

	if (semafor == "stop") {
		return false;
	}

	if (semafor == "resume") {
		p_MonitorHandler->restoreState(snapshotPath);
		fstream fs2(semaforFile.c_str(), ifstream::out);
		fs2.clear();
		fs2 << "continue";
		fs2.close();

		fstream nEv((output_Path + "/snapshots/numberOfEvents.dat").c_str(),
				fstream::in);
		double tmp;
		nEv >> tmp;
		i = tmp + 1;
		nbEvents += tmp;
		nEv.close();

		p_randomGenerator->restoreStatus(
				(output_Path + "/snapshots/generator.dat").c_str());
	}

	return true;
}

bool VINHAC::Manager::checkSemafor(string semaforFile, double evntNb) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::checkSemafor()"<<std::endl;

#endif
	string semafor;
	ifstream ifs(semaforFile.c_str(), ifstream::in);
	ifs >> semafor;
	ifs.close();

	if (semafor == "continue") {
		return true;
	} else if (semafor == "stop")
		return false;

	else
		return true;

}

bool VINHAC::Manager::generateEvents(double nbEvents) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::GenerateEvents()"<<std::endl;

#endif
	double ratio = nbEvents / 100;
	if (ratio == 0)
		ratio = 100;
	double i = 1;

	for (i = 1; i < nbEvents + 1; i++) {
		if (i == 1) {
			string snapshotPath = output_Path + "/snapshots/";
			string semafor = input_Path + "/semaphore/semaphore";
			if (checkSemafor(semafor, snapshotPath, i, nbEvents) == false)
				return false;

		}

		Event evt = generateEvent();

		if (i == 1) {
			cout << evt.printInLabFrame() << endl;
			cout << evt.printInQuarksCMSFrame() << endl;
		}
		if (fmod(i, ratio) == 0) {
			cout << "####################" << i << "####################"
					<< endl;
		}
		if (fmod(i, 100000) == 0) {
			string snapshotPath = output_Path + "/snapshots/";
			makeSnapshot(snapshotPath, i);
			string semafor = input_Path + "/semaphore/semaphore";
			if (checkSemafor(semafor, i) == false)
				break;
		}

	}
	string output = output_Path;
	makeSnapshot(output, i - 1);
	return true;
}

void VINHAC::Manager::makeSnapshot(string snapshotPath, double eventNumber) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::makeSnapshot()"<<std::endl;

#endif
	//// Convert event number to string

	//// Prepare paths

	string generatorFile = snapshotPath + "/generator" + ".dat";
	string generatorFileBak = snapshotPath + "/generator" + ".bak";

	string mv_generator_comand = "mv " + generatorFile + " " + generatorFileBak;
	// Snapshot of the random number generator
	system(mv_generator_comand.c_str());

	string number_of_eventsFile = snapshotPath + "/numberOfEvents" + ".dat";
	string number_of_eventsFileBak = snapshotPath + "/numberOfEvents" + ".bak";

	string mv_number_of_events_comand = "mv " + number_of_eventsFile + " "
			+ number_of_eventsFileBak;
	// Snapshot of the random number generator
	system(mv_number_of_events_comand.c_str());

	fstream nEv(number_of_eventsFile.c_str(), fstream::out);
	nEv << eventNumber;
	nEv.close();

	p_randomGenerator->saveStatus(generatorFile.c_str());
	// Snapshot of the analysis handler

	p_MonitorHandler->makeSnapshot(snapshotPath);

	string summaryFile = snapshotPath + "/summary" + ".txt";
	string summaryFileBak = snapshotPath + "/sumary" + ".bak";

	string summarycomand = "mv " + summaryFile + " "
			+ summaryFileBak;
	// Snapshot of the summary
	system(summarycomand.c_str());
	fstream summary(summaryFile.c_str(), fstream::out);
	printSummary(summary);

}
void VINHAC::Manager::removeProcessHandler(const string& aName) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::removeProcessHandler()"<<std::endl;

#endif
	vector<ProcessHandler*>::iterator it;
	for (it = _handlers.begin(); it != _handlers.end(); it++) {
		if ((*it)->getHandlerName() == aName) {
			_handlers.erase(it);
		}
	}
}

void VINHAC::Manager::addProcessHandler(VINHAC::ProcessHandler* pointerProcess) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::addProcessHandler()"<<std::endl;

#endif

	_handlers.push_back(pointerProcess);
}

VINHAC::Event VINHAC::Manager::generateEvent() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::GenerateEvent() 2"<<std::endl;

#endif
	Event event;

	while (true) {
		event.clean();

		vector<ProcessHandler*>::iterator it;

		// show content:
		for (it = _handlers.begin(); it != _handlers.end(); ++it) {
			(*it)->eventEvolve(event);
		}

		if (DataSource::get().keyWgt == 0) {
			double randomNumber = p_randomGenerator->Flat();
			double weight = event.mainWt();
			if (randomNumber * DataSource::get().maxWeightRej > abs(weight)) {
				continue;
			} else {
				event.setMainWt(copysign(1.0, weight));
				break;
			}
		} else {
			break;
		}
	}
#ifdef DEBUG_FLAG

	std::cout<<event.printInLabFrame()<<std::endl;
	std::cout<<event.printInQuarksCMSFrame()<<std::endl;

#endif
	return event;

}

double VINHAC::Manager::getFinalXsection() {
	double result = p_MonitorHandler->getMainWtMonitor().averageWeight(
			p_MonitorHandler->getNumberOfEvents());
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::getFinalXsection() # "<<result<<std::endl;

#endif
	result *= p_BeamHandler->getXsCrud();

	return result;
}
double VINHAC::Manager::getFinalXsectionError() {
	double result = p_MonitorHandler->getMainWtMonitor().averageWeight(
			p_MonitorHandler->getNumberOfEvents());
	result *= p_MonitorHandler->getMainWtMonitor().error(
			p_MonitorHandler->getNumberOfEvents());
	result *= p_BeamHandler->getXsCrud();

	return result;
}

void VINHAC::Manager::printCrudeXsectionTimesWt() const {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::printCrudeXsectionTimesWt()"<<std::endl;

#endif

	const map<WeightName::name, VINHAC::WeightMonitor>& weightMonitors =
			getMonitorHandler()->getCommonWtMonitors();
	map<WeightName::name, VINHAC::WeightMonitor>::const_iterator it;

	for (it = weightMonitors.begin(); it != weightMonitors.end(); it++)
		cout << (*it).first << " = " << (*it).second.averageWeight() << " +/- "
				<< (*it).second.averageWeight() * (*it).second.error()
				<< ", Sum of wt2 = " << (*it).second.getSumWt2() << endl;
}

double VINHAC::Manager::getCrudeXsection() const {
	return p_BeamHandler->getXsCrud();
}

#ifdef HEPMC_ENABLED
HepMC::GenEvent VINHAC::Manager::toHepMC(Event& event) const {
	HepMC::GenEvent result;
	result.use_units(HepMC::Units::GEV, HepMC::Units::MM);

	/************** SETUP PARTICLES AND VERTICES ****************/
	HepMC::GenParticle* beamParticle1;
	HepMC::GenParticle* beamParticle2;

	HepMC::GenParticle* quark = new HepMC::GenParticle(
			event.getQuark().getFourMomentum(), event.getQuark().getPDGid(), 4);
	beamParticle1 = quark;
	if (event.getQuark().getParents().size() > 0) {
		HepMC::GenVertex* proton_quark = new HepMC::GenVertex();
		result.add_vertex(proton_quark);
		Particle* proton = event.getQuark().getParents()[0];
		HepMC::GenParticle* proton_in = new HepMC::GenParticle(
				proton->getFourMomentum(), proton->getPDGid(), 4);
		beamParticle1 = proton_in;
		proton_quark->add_particle_in(proton_in);
		proton_quark->add_particle_out(quark);
	}

	HepMC::GenParticle* antiquark = new HepMC::GenParticle(
			event.getAntiQuark().getFourMomentum(),
			event.getAntiQuark().getPDGid(), 4);
	beamParticle2 = antiquark;
	if (event.getAntiQuark().getParents().size() > 0) {
		HepMC::GenVertex* proton_antiquark = new HepMC::GenVertex();
		result.add_vertex(proton_antiquark);
		Particle* proton = event.getAntiQuark().getParents()[0];
		HepMC::GenParticle* proton_in = new HepMC::GenParticle(
				proton->getFourMomentum(), proton->getPDGid(), 4);
		beamParticle2 = proton_in;
		proton_antiquark->add_particle_in(proton_in);
		proton_antiquark->add_particle_out(antiquark);
	}
	result.set_beam_particles(beamParticle1, beamParticle2);

	HepMC::GenVertex* quarks_boson = new HepMC::GenVertex();
	result.add_vertex(quarks_boson);
	quarks_boson->add_particle_in(quark);
	quarks_boson->add_particle_in(antiquark);
	HepMC::GenParticle* boson =
	new HepMC::GenParticle(event.getBoson().getFourMomentum(),
			event.getBoson().getPDGid(), 11);
	quarks_boson->add_particle_out(boson);

	HepMC::GenVertex* boson_output = new HepMC::GenVertex();
	result.add_vertex(boson_output);
	boson_output->add_particle_in(boson);
	HepMC::GenParticle* lepton = new HepMC::GenParticle(
			event.getLepton().getFourMomentum(), event.getLepton().getPDGid(),
			1);
	boson_output->add_particle_out(lepton);
	HepMC::GenParticle* neutrino = new HepMC::GenParticle(
			event.getNeutrino().getFourMomentum(),
			event.getNeutrino().getPDGid(), 1);
	boson_output->add_particle_out(neutrino);
	for (vector<Particle>::iterator it = event.getPhotons().begin(); it
			!= event.getPhotons().end(); ++it) {
		HepMC::GenParticle* photon = new HepMC::GenParticle(
				(*it).getFourMomentum(), (*it).getPDGid(), 1);
		boson_output->add_particle_out(photon);
	}
	/************** SETUP PARTICLES AND VERTICES - END **********/

	/************** SETUP PDF INFO ******************************/
	if (event.isPDFUsed()) {
		HepMC::PdfInfo pdfInfo;
		pdfInfo.set_id1(event.getQuark().getPDGid());
		pdfInfo.set_id2(event.getAntiQuark().getPDGid());
		pdfInfo.set_pdf_id1(event.getQuarkPDFid());
		pdfInfo.set_pdf_id2(event.getAntiquarkPDFid());
		pdfInfo.set_x1(event.getQuarkX());
		pdfInfo.set_x2(event.getAntiquarkX());
		pdfInfo.set_scalePDF(event.getPDFScale());
		pdfInfo.set_pdf1(event.getQuarkXpdf());
		pdfInfo.set_pdf2(event.getAntiquarkXpdf());
		result.set_pdf_info(pdfInfo);
	}
	/************** SETUP PDF INFO - END ************************/

	/************** SETUP EVENT WEIGHT **************************/
	HepMC::WeightContainer& weights = result.weights();
	double xscrude = getCrudeXsection();
	xscrude *= event.currentCommonWtMultiplication();
	weights[Commons::mapWeightPolarizationToString(DataSource::get().mainWeight)] = event.getModelWeightByName(DataSource::get().mainWeight);
	for (map<pair<WeightName::name,PolarizationName::name>, double>::const_iterator it =
			event.getModelWeightsMap().begin(); it
			!= event.getModelWeightsMap().end(); ++it) {
		if(it->first != DataSource::get().mainWeight) {
			weights[Commons::mapWeightPolarizationToString(it->first)] = it->second;
		}
	}
	/************** SETUP EVENT WEIGHT - END*********************/
	return result;
}
#endif
