#include <string>
#include <vector>
#include <exception>
#include <map>
using namespace std;

#include "InitializationHandler.h"

VINHAC::InitializationHandler::InitializationHandler(DataSource *ds) :
		p_Manager(0),
		p_HardProcessHandler(0),
		p_BeamHandler(0),
		p_InputFileReader(0),
		ds(ds) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::InitializationHandler()"<<std::endl;

#endif

}

void VINHAC::InitializationHandler::initializeDS(){
	ds->MW = p_InputFileReader->getParticleMass(24);
	ds->MZ = p_InputFileReader->getParticleMass(23);
	string Zwidth = p_InputFileReader->getParticleAttribute(23, "mWidth");
	ds->GZ = strtod(Zwidth.c_str(), NULL);
	string Wwidth = p_InputFileReader->getParticleAttribute(24, "mWidth");
	ds->GW = strtod(Wwidth.c_str(), NULL);

	//TODO Higgs mass from input
	ds->MH = 150.0;
	ds->GH = 1.0;

	double sinW2 = p_InputFileReader->getConstant("sinThetaW2");
	if( sinW2 < 0 ){
			sinW2 = 1.0 - (ds->MW/ds->MZ)*(ds->MW/ds->MZ);
	}
	ds->sinThetaW2 = sinW2;
	// read fermion masses

	vector<double> fermionMasses = readFermionMasses();
	for(unsigned int i = 0 ;i<fermionMasses.size(); ++i){
		ds->masses[static_cast<int>(i)] = fermionMasses[i];
	}
	ds->masses[2212] = p_InputFileReader->getParticleMass(2212);

	ds->quarksCharges = p_InputFileReader->readQuarksCharges();
	ds->pi = p_InputFileReader->getConstant("pi");
	ds->invAlphaQED = p_InputFileReader->getConstant("invAlphaQED");
	ds->alphaQED = 1.0/ds->invAlphaQED;

	ds->initialParticles = p_InputFileReader->readInitialParticles();
	// read bosons list
	ds->intermediateParticles =
			p_InputFileReader->readIntermediateParticles();
	// read final state particles
	vector<int> finalParticlesTmp = p_InputFileReader->readFinalParticles();
	vector<int> finalParticles;

	for(unsigned i = 0 ; i< finalParticlesTmp.size(); ++i){
		bool wasE = false; bool wasMu = false; bool wasTau = false;
		if((finalParticlesTmp[i] == e || finalParticlesTmp[i] == -e)&&!wasE){
			finalParticles.push_back(e);
			wasE=true;
		}
		if((finalParticlesTmp[i] == mu || finalParticlesTmp[i] == -mu)&&!wasMu){
			finalParticles.push_back(mu);
			wasMu=true;
		}
		if((finalParticlesTmp[i] == tau || finalParticlesTmp[i] == -tau)&&!wasTau){
			finalParticles.push_back(tau);
			wasTau=true;
		}
	}

	ds->finalParticles = finalParticles;

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::initializeDS() initial particles : ";
	for(vector<int>::iterator it = ds->initialParticles.begin(); it!=ds->initialParticles.end(); ++it){
		std::cout<<(*it)<<" ";
	}
	std::cout<<std::endl;

#endif

	string polarization = p_InputFileReader->getPolarization();
	if(polarization=="Unpolarized"){
		ds->keyPol = 0;
	} else if(polarization=="Transversely"){
		ds->keyPol = 1;
	} else if(polarization=="Longitudinally"){
		ds->keyPol = 2;
	} else {
		std::cout<<"Unparsable choise of polarization !!"<<std::endl;
		throw VinhacException("Unparsable choise of polarization !!");
	}


	ds->fermiConst = p_InputFileReader->getConstant("FermiConst");
	ds->invGeV2toNb = p_InputFileReader->getConstant("invGeV2toNb");
	ds->geV2toNb = 1.0/ds->invGeV2toNb;

	ds->q2MinCut = p_InputFileReader->getQ2minCut();

	ds->beamA = initializeBeam("beamA");
	ds->beamB = initializeBeam("beamB");

	ds->EgMin = 1.0E-3; //TODO

	ds->keyGmu = p_InputFileReader->getIPSscheme();

	ds->keyRad = p_InputFileReader->getRadiativeCorrection();
	if(ds->keyRad<0) throw VinhacException("Bad settings of radiative corrections");

	ds->keyFac = 0; //todo
	ds->keyWid = 0; //todo
	ds->keyCol = 0; //todo

	ds->ckm = p_InputFileReader->readCKM();
	ds->ckm2 = p_InputFileReader->readCKM();

	for(unsigned i = 0; i<ds->ckm2.size(); ++i){
		for(unsigned j =0 ; j<ds->ckm2[i].size(); ++j){
			ds->ckm2[i][j] = ds->ckm2[i][j]*ds->ckm2[i][j];
		}
	}
#ifdef DEBUG_FLAG
	for(unsigned i = 0; i<ds->ckm2.size(); ++i){
		std::cout<<"DEBUG: InitializationHandler:: ckm^2 ";
		for(unsigned j =0 ; j<ds->ckm2[i].size(); ++j){
			cout<<ds->ckm2[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
#endif
}

void VINHAC::InitializationHandler::initializeAllHandlers() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::InitializeAllHandlers()"<<std::endl;

#endif

	if(ds==0){
		throw VinhacException("DataSource not properly initialized");
	}

	//initialize Reader
	initializeReader();
	//reads settings to DataSource structure
	initializeDS();



	//initialize BeamHandler
	//Create Beam Handler
	p_BeamHandler = new BeamHandler(ds,p_Manager->getRandomNumberGenerator());

	p_HardProcessHandler = new HardProcessHandler(ds);
	p_HardProcessHandler->setRandomNumberGenerator(p_Manager->getRandomNumberGenerator());


	//Create Beam Handler
	p_ModelHandler = new ModelHandler(ds, p_HardProcessHandler);




	//initialize MonitorHandler
	p_MonitorHandler = new MonitorHandler();
}

std::vector<double> VINHAC::InitializationHandler::readFermionMasses(){
		std::vector<double> fermionMasses;
		for (int i = 0; i<16; ++i) fermionMasses.push_back(0.0);
		for (int i = 1; i < 7; i++){
			fermionMasses[i] = p_InputFileReader->getParticleMass(i);
		}
			fermionMasses[11] = p_InputFileReader->getParticleMass(11);
			fermionMasses[13] = p_InputFileReader->getParticleMass(13);
			fermionMasses[15] = p_InputFileReader->getParticleMass(15);

		return fermionMasses;
}

bool VINHAC::InitializationHandler::initializeReader() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::InitializeReader()"<<std::endl;

#endif

	// Create and initialize XercescReader
	p_InputFileReader = new XercescReader(inputPath+"/InputData.xml",
			inputPath+"/ParticleDataBase.xml");
	// Load default file
	p_InputFileReader->loadDefaultFile();

	// Load user's file
	if (m_userFilePath.size() != 0) {
		cout << "Reading user file: " << m_userFilePath << endl;
		p_InputFileReader->loadUserFile(m_userFilePath);
	} else {

	}

	return true;
}

bool VINHAC::InitializationHandler::inputEcho() {
	p_InputFileReader->inputEcho(outputPath);

	return true;
}

std::map<string, string> VINHAC::InitializationHandler::initializeBeam(
		string beam_name) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::initializeBeam()"<<std::endl;

#endif
	std::cout << "====================Initialize " << beam_name
			<< "===================" << std::endl;
	try {

		std::map<string, string> beamSettings;
		beamSettings = p_InputFileReader->readBeamSettings(beam_name);

		map<string, string>::iterator it;

#ifndef NOTEST // For tests
		for (it = beamSettings.begin(); it != beamSettings.end(); it++)
		cout << (*it).first << " => " << (*it).second << endl;
#endif

		return beamSettings;

	} catch (const std::runtime_error e) {
		std::cerr << "Error loading config: " << e.what() << std::endl;
		throw VinhacException(string("Error loading config: ")+string(e.what()));
	}

}


VINHAC::InitializationHandler::~InitializationHandler() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::~InitializationHandler()"<<std::endl;

#endif
	if (p_InputFileReader) {
		delete p_InputFileReader;
		p_InputFileReader = 0;
	}

	if (p_HardProcessHandler) {
		delete p_HardProcessHandler;
		p_HardProcessHandler = 0;
	}

	if (p_BeamHandler) {
		delete p_BeamHandler;
		p_BeamHandler = 0;
	}


	if (p_MonitorHandler){
		delete p_MonitorHandler;
		p_MonitorHandler = 0;
	}

	if (p_ModelHandler) {
			delete p_ModelHandler;
			p_ModelHandler = 0;
		}

}
