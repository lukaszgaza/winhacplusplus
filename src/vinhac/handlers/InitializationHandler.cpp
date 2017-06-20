#include <string>
#include <vector>
#include <exception>
#include <map>
#include <cstdlib>
using namespace std;

#include "InitializationHandler.h"
#include "../input/DataSource.h"
#include "../input/Reader.h"
#include "HardProcessHandler.h"
#include "../core/VinhacException.h"
#include "BeamHandler.h"
#include "MonitorHandler.h"
#include "../core/Manager.h"
#include "ModelHandler.h"
#include "../input/XercescReader.h"
#include "QCDHandler.h"

VINHAC::InitializationHandler::InitializationHandler() :
		p_Manager(0),
		p_HardProcessHandler(0),
		p_BeamHandler(0),
		p_ModelHandler(0),
		p_QCDHandler(0),
		p_InputFileReader(0) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::InitializationHandler()"<<std::endl;

#endif

}

void VINHAC::InitializationHandler::initializeDS(){
	DataSource::get().printOutLevel = p_InputFileReader->getPrintOutLevel();

	DataSource::get().MW = p_InputFileReader->getParticleMass(24);
	DataSource::get().MZ = p_InputFileReader->getParticleMass(23);
	string Zwidth = p_InputFileReader->getParticleAttribute(23, "mWidth");
	DataSource::get().GZ = strtod(Zwidth.c_str(), NULL);
	string Wwidth = p_InputFileReader->getParticleAttribute(24, "mWidth");
	DataSource::get().GW = strtod(Wwidth.c_str(), NULL);

	DataSource::get().MH = p_InputFileReader->getConstant("HiggsMass");
	DataSource::get().GH = p_InputFileReader->getConstant("HiggsWidth");

	double sinW2 = p_InputFileReader->getConstant("sinThetaW2");
	if( sinW2 < 0 ){
			sinW2 = 1.0 - (DataSource::get().MW/DataSource::get().MZ)*(DataSource::get().MW/DataSource::get().MZ);
	}
	DataSource::get().sinThetaW2 = sinW2;
	// read fermion masses

	vector<double> fermionMasses = readFermionMasses();
	for(unsigned int i = 0 ;i<fermionMasses.size(); ++i){
		DataSource::get().masses[static_cast<int>(i)] = fermionMasses[i];
	}
	DataSource::get().masses[2212] = p_InputFileReader->getParticleMass(2212);

	DataSource::get().quarksCharges = p_InputFileReader->readQuarksCharges();
	DataSource::get().pi = p_InputFileReader->getConstant("pi");
	DataSource::get().invAlphaQED = p_InputFileReader->getConstant("invAlphaQED");
	DataSource::get().alphaQED = 1.0/DataSource::get().invAlphaQED;
	DataSource::get().alphaQCD = p_InputFileReader->getConstant("alphaQCD");

	DataSource::get().initialParticles = p_InputFileReader->readInitialParticles();
	// read bosons list
	DataSource::get().intermediateParticles =
			p_InputFileReader->readIntermediateParticles();
	// read final state particles
	vector<int> finalParticlesTmp = p_InputFileReader->readFinalParticles();
	vector<int> finalParticles;

	for(unsigned i = 0 ; i< finalParticlesTmp.size(); ++i){
		bool wasE = false; bool wasMu = false; bool wasTau = false;
		if((finalParticlesTmp[i] == electron || finalParticlesTmp[i] == -electron)&&!wasE){
			finalParticles.push_back(electron);
			wasE=true;
		}
		if((finalParticlesTmp[i] == muon || finalParticlesTmp[i] == -muon)&&!wasMu){
			finalParticles.push_back(muon);
			wasMu=true;
		}
		if((finalParticlesTmp[i] == tau || finalParticlesTmp[i] == -tau)&&!wasTau){
			finalParticles.push_back(tau);
			wasTau=true;
		}
	}

	DataSource::get().finalParticles = finalParticles;

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::initializeDS() initial particles : ";
	for(vector<int>::iterator it = DataSource::get().initialParticles.begin(); it!=DataSource::get().initialParticles.end(); ++it){
		std::cout<<(*it)<<" ";
	}
	std::cout<<std::endl;

#endif

	DataSource::get().polarizations = p_InputFileReader->getPolarizations();
	if(DataSource::get().polarizations.size()>1 ||
			DataSource::get().polarizations.count(PolarizationName::unpolarized) == 0){
		DataSource::get().keyPol = 1;
	} else {
		DataSource::get().keyPol = 0;
	}


	DataSource::get().pythiaXMLdocDir = p_InputFileReader->getPythiaXMLdocDir();
	DataSource::get().isPartonShower = p_InputFileReader->isPartonShower();
	DataSource::get().isHadronization = p_InputFileReader->isHadronization();

	DataSource::get().fermiConst = p_InputFileReader->getConstant("FermiConst");
	DataSource::get().invGeV2toNb = p_InputFileReader->getConstant("invGeV2toNb");
	DataSource::get().geV2toNb = 1.0/DataSource::get().invGeV2toNb;

	DataSource::get().isPDFOn = p_InputFileReader->isPdfOn();
	DataSource::get().pdfInterface = p_InputFileReader->getPdfInterface();
	DataSource::get().pdfName = p_InputFileReader->getPdfName();
	DataSource::get().pdfSubset = p_InputFileReader->getPdfSubset();
	DataSource::get().pdfXMin = p_InputFileReader->getPdfXMin();
	DataSource::get().pdfXMax = p_InputFileReader->getPdfXMax();
	DataSource::get().pdfQ2min = p_InputFileReader->getPdfQ2min();

	DataSource::get().beamA = initializeBeam("A");
	DataSource::get().beamB = initializeBeam("B");

	DataSource::get().EgMin = p_InputFileReader->getSoftPhotonCutOff();

	DataSource::get().keyGmu = p_InputFileReader->getAlphaScheme();

	DataSource::get().keyRad = p_InputFileReader->getRadiativeCorrection();
	DataSource::get().mainWeight = p_InputFileReader->getMainWeight();

	DataSource::get().modelReferenceFrame = p_InputFileReader->getModelReferenceFrame();

	DataSource::get().keyFac = p_InputFileReader->getFactorizationScale();
	DataSource::get().keyWid = p_InputFileReader->getWidthScheme();
	DataSource::get().keyCol = 0; //todo


	DataSource::get().ckm = p_InputFileReader->readCKM();
	DataSource::get().ckm2 = p_InputFileReader->readCKM();

	for(unsigned i = 0; i<DataSource::get().ckm2.size(); ++i){
		for(unsigned j =0 ; j<DataSource::get().ckm2[i].size(); ++j){
			DataSource::get().ckm2[i][j] = DataSource::get().ckm2[i][j]*DataSource::get().ckm2[i][j];
		}
	}

	DataSource::get().keyWgt = p_InputFileReader->isWeightedEvent();
	DataSource::get().maxWeightRej = p_InputFileReader->getMaxWeightRejection();

#ifdef DEBUG_FLAG
	for(unsigned i = 0; i<DataSource::get().ckm2.size(); ++i){
		std::cout<<"DEBUG: InitializationHandler:: ckm^2 ";
		for(unsigned j =0 ; j<DataSource::get().ckm2[i].size(); ++j){
			cout<<DataSource::get().ckm2[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
#endif
}

void VINHAC::InitializationHandler::initializeAllHandlers() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: InitializationHandler::InitializeAllHandlers()"<<std::endl;

#endif

	//initialize Reader
	initializeReader();

	inputEcho();
	//reads settings to DataSource structure
	initializeDS();



	//initialize BeamHandler
	//Create Beam Handler
	p_BeamHandler = new BeamHandler(p_Manager->getRandomNumberGenerator());
	p_Manager->p_BeamHandler = p_BeamHandler;

	p_HardProcessHandler = new HardProcessHandler();
	p_HardProcessHandler->setRandomNumberGenerator(p_Manager->getRandomNumberGenerator());


	//Create QCDHandler
	if((DataSource::get().isPartonShower
			|| DataSource::get().isHadronization)) {
		p_QCDHandler = new QCDHandler(p_Manager->getRandomNumberGenerator(), *p_Manager,
				inputPath, DataSource::get().pythiaXMLdocDir);
	}

	//Create Frame Transformation Handler
	p_FrameTransformationHandler = new FrameTransformationHandler(p_Manager->getRandomNumberGenerator());


	//Create Model Handler
	p_ModelHandler = new ModelHandler( p_HardProcessHandler);




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
	if(DataSource::get().printOutLevel > 0){
		std::cout << "====================Initialize Beam " << beam_name
			<< "===================" << std::endl;
	}
	try {

		std::map<string, string> beamSettings;
		beamSettings = p_InputFileReader->readBeamSettings(beam_name);

		map<string, string>::iterator it;

		if(DataSource::get().printOutLevel > 0){
			for (it = beamSettings.begin(); it != beamSettings.end(); it++)
			cout << (*it).first << " => " << (*it).second << endl;
		}

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

	if(p_FrameTransformationHandler){
		delete p_FrameTransformationHandler;
		p_FrameTransformationHandler = 0;
	}

	if (p_QCDHandler) {
		delete p_QCDHandler;
		p_QCDHandler = 0;
	}
}
