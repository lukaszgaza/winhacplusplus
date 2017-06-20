#include <string>
#include <vector>
#include <exception>
#include <ostream>
using namespace std;

#include "Manager.h"
#include "VinhacException.h"
//// Constructors
VINHAC::Manager::Manager(string aInputStream) :
	m_userFilePath(aInputStream), p_InitializationHandler(0), ds() {
#ifdef DEBUG_FLAGw

	std::cout<<"DEBUG: Manager::Manager() 1"<<std::endl;

#endif
	std::cout << ">>>>>>>> Reading user file: " << m_userFilePath << endl;
	m_snapshotNb = 50;
}

VINHAC::Manager::Manager() :
	m_userFilePath(), p_InitializationHandler(0), ds() {
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

	p_InitializationHandler = new InitializationHandler(&ds);

	p_InitializationHandler->setUserFilePath(m_userFilePath);
	p_InitializationHandler->setInputPath(input_Path);
	p_InitializationHandler->setOutputPath(output_Path);

	p_InitializationHandler->setManager(this);

	initializeRandomNumberGenerator();

	//// Initialise handlers
	p_InitializationHandler->initializeAllHandlers();

	//// add them to the handler vector. Becareful have to be added in the proper order

	_handlers.push_back(p_InitializationHandler->getBeamHandler());
	//_handlers.push_back(p_InitializationHandler->getBornHandler());
	_handlers.push_back(p_InitializationHandler->getHardProcessHandler());
	_handlers.push_back(p_InitializationHandler->getModelHandler());
	//_handlers.push_back(p_InitializationHandler->getAnalisysHandler());
	_handlers.push_back(p_InitializationHandler->getMonitorHandler());
	p_MonitorHandler = p_InitializationHandler->getMonitorHandler();
	p_BeamHandler = p_InitializationHandler->getBeamHandler();

	//// Print out used input date
	p_InitializationHandler->inputEcho();

	fstream crudefile((output_Path + "/XSectionCrude.dat").c_str(),
			fstream::out);
	crudefile << p_BeamHandler->getXsCrud() << ":"
			<< p_BeamHandler->getXsCrudError();
	crudefile.close();

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::initializeGenerator() # handlers.size(): "<<_handlers.size()<<std::endl;

#endif

}

void VINHAC::Manager::printSummary() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::printSummary()"<<std::endl;

#endif
	this->p_MonitorHandler->printAverageWeights();
	this->p_MonitorHandler->printCrossSections();
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

		p_randomGenerator->restoreStatus((output_Path
				+ "/snapshots/generator.dat").c_str());
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

void VINHAC::Manager::makeSnapshot(string snapshotPath, int eventNumber) {
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

}
void VINHAC::Manager::removeProcessHandler(string aName) {
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
	event.clean();

	vector<ProcessHandler*>::iterator it;

	// show content:
	for (it = _handlers.begin(); it != _handlers.end(); ++it) {
		(*it)->eventEvolve(event);
	}
#ifdef DEBUG_FLAG

	std::cout<<event.printInLabFrame()<<std::endl;
	std::cout<<event.printInQuarksCMSFrame()<<std::endl;

#endif
	return event;

}

double VINHAC::Manager::getFinalXsection() {
	double result = p_MonitorHandler->getMainWtMonitor().averageWeight();
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::getFinalXsection() # "<<result<<std::endl;

#endif
	result *= p_BeamHandler->getXsCrud();

	return result;
}
double VINHAC::Manager::getFinalXsectionError() {
	double result = p_MonitorHandler->getMainWtMonitor().averageWeight();
	result *= p_MonitorHandler->getMainWtMonitor().error();
	result *= p_BeamHandler->getXsCrud();

	return result;
}

void VINHAC::Manager::printCrudeXsectionTimesWt() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: Manager::printCrudeXsectionTimesWt()"<<std::endl;

#endif

	map<string, VINHAC::WeightMonitor> weightMonitors =
			getMonitorHandler()->getCommonWtMonitors();
	map<string, VINHAC::WeightMonitor>::iterator it;

	for (it = weightMonitors.begin(); it != weightMonitors.end(); it++)
		cout << (*it).first << " = " << (*it).second.averageWeight() << " +/- "
				<< (*it).second.averageWeight() * (*it).second.error()
				<< ", Sum of wt2 = " << (*it).second.getSumWt2() << endl;
}

double VINHAC::Manager::getCrudeXsection() {
	return p_BeamHandler->getXsCrud();
}
