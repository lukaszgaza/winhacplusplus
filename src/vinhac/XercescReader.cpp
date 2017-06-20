/*
 * XercescReader.cpp
 *
 *  Created on: 2009-05-13
 *      Author: siudmy
 */
#include "XercescReader.h"

VINHAC::XercescReader::XercescReader(const string& fromFile,
		const string& particle) :
	xmlFile_(fromFile), xmlFileParticle_(particle) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::XercescReader()"<<std::endl;

#endif

	try {

		xercesc::XMLPlatformUtils::Initialize();

	} catch (xercesc::XMLException& e) {

		std::cerr << "XML toolkit initialization error: " << DualString(
				e.getMessage()) << std::endl;

		throw(ERROR_XERCES_INIT);

	}

	parser_ = new XercesDOMParser();
	finder_ = new ElementFinder();
	finderParticle_ = new ElementFinder();
	finderUserFile_ = new ElementFinder();

}

void VINHAC::XercescReader::loadDefaultFile() throw (runtime_error) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::loadDefaultFile()"<<std::endl;

#endif

	parser_->setValidationScheme(XercesDOMParser::Val_Never);
	parser_->setDoNamespaces(false);
	parser_->setDoSchema(false);
	parser_->setLoadExternalDTD(false);

	try {

		parser_->parse(xmlFile_.c_str());

		// there's no need to free this pointer -- it's
		// owned by the parent parser object
		DOMDocument* xmlDoc = parser_->getDocument();

		if (NULL == xmlDoc->getDocumentElement()) {
			throw(runtime_error("empty XML document"));
		}

		finder_->setDocument(xmlDoc);

		parser_->parse(xmlFileParticle_.c_str());
		DOMDocument* xmlParticle = parser_->getDocument();

		if (NULL == xmlParticle->getDocumentElement()) {
			throw(runtime_error("empty XML document"));
		}

		finderParticle_->setDocument(xmlParticle);

	} catch (XMLException& e) {

		// believe it or not, XMLException is not
		// a parent class of DOMException

		ostringstream buf;
		buf << "Error parsing file: " << DualString(e.getMessage()) << flush;

		throw(runtime_error(buf.str()));

	} catch (const DOMException& e) {

		ostringstream buf;
		buf << "Encountered DOM Exception: " << DualString(e.getMessage())
				<< flush;
		throw(runtime_error(buf.str()));
	}
	return;
} // load()

void VINHAC::XercescReader::loadUserFile(string userFile) throw (runtime_error) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::loadUserFile()"<<std::endl;

#endif

	parser_->setValidationScheme(XercesDOMParser::Val_Never);
	parser_->setDoNamespaces(false);
	parser_->setDoSchema(false);
	parser_->setLoadExternalDTD(false);

	try {

		parser_->parse(userFile.c_str());

		// there's no need to free this pointer -- it's
		// owned by the parent parser object
		DOMDocument* xmlUserDoc = parser_->getDocument();

		if (NULL == xmlUserDoc->getDocumentElement()) {
			throw(runtime_error("empty XML document"));
		}

		finder_->setUserDocument(xmlUserDoc);

	} catch (XMLException& e) {

		// believe it or not, XMLException is not
		// a parent class of DOMException

		ostringstream buf;
		buf << "Error parsing file: " << DualString(e.getMessage()) << flush;

		throw(runtime_error(buf.str()));

	} catch (const DOMException& e) {

		ostringstream buf;
		buf << "Encountered DOM Exception: " << DualString(e.getMessage())
				<< flush;
		throw(runtime_error(buf.str()));
	}
	return;
} // load()


bool VINHAC::XercescReader::initialize() {

	return true;
}

void VINHAC::XercescReader::inputEcho(string path) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::inputEcho()"<<std::endl;

#endif

	finder_->inputEcho(path);
}
/*! \brief Method which read from input a beam settings and returns it as
 * 	a map.
 *  \param beamName - name of the beam
 *
 *
 */
map<string, string> VINHAC::XercescReader::readBeamSettings(
		const string& beamName) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::readBeamSettings()"<<std::endl;

#endif

	// Map with beam settings map <keyName, content>
	map<string, string> beamSettings;
	StringManager smBeam;

	const DualString TAG_BEAM(beamName);

	//Beam element
	DOMElement* element = finder_->getElement(TAG_BEAM.asXMLString());

	//Energy of the beam
	const char* beamEnergy = smBeam.convert(element->getAttribute(
			finder_->ATTR_ENERGY_.asXMLString()));
	beamSettings[finder_->ATTR_ENERGY_.asCString()] = beamEnergy;

	//idPDG of the beam
	DOMElement* idPDGelement = finder_->getSubElement(
			finder_->TAG_idPDG_.asXMLString(), element); // Find idPDGelement
	beamSettings[finder_->TAG_idPDG_.asCString()] = smBeam.convert(
			idPDGelement->getTextContent());
	//idPDGcontent.asCString();

	////PDF content
	//keyPDFelement
	DOMElement* keyPDFelement = finder_->getSubElement(
			finder_->TAG_keyPDF_.asXMLString(), element);
	//PDF name
	const char* PDFname = smBeam.convert(keyPDFelement->getAttribute(
			finder_->ATTR_name_.asXMLString()));
	beamSettings[finder_->ATTR_name_.asCString()] = PDFname;

	const char* PDFinterface = smBeam.convert(keyPDFelement->getAttribute(
			finder_->ATTR_PDFinterface_.asXMLString()));
	beamSettings[finder_->ATTR_PDFinterface_.asCString()] = PDFinterface;

	const DualString keyPDFcontent(keyPDFelement->getTextContent());

	const DualString xMinPDFcontent(keyPDFelement->getTextContent());
	//xMinPDF
	beamSettings[finder_->TAG_xMinPDF_.asCString()]
			= finder_->getSubElementContent(
					finder_->TAG_xMinPDF_.asXMLString(), element);
	//xMaxPDF
	beamSettings[finder_->TAG_xMaxPDF_.asCString()]
			= finder_->getSubElementContent(
					finder_->TAG_xMaxPDF_.asXMLString(), element);;

	beamSettings[finder_->TAG_PDFsubset_.asCString()]
			= finder_->getSubElementContent(
					finder_->TAG_PDFsubset_.asXMLString(), element);
	////HeavyIon settings

	DOMElement* heavyIonElement = finder_->getSubElement(
			finder_->TAG_heavyIon_.asXMLString(), element);

	const char* heavyIonSwitch = smBeam.convert(heavyIonElement->getAttribute(
			finder_->ATTR_switchOn_.asXMLString()));

	beamSettings[finder_->ATTR_switchOn_.asCString()] = heavyIonSwitch;

	//atomicNb
	beamSettings[finder_->TAG_atomicNb_.asCString()]
			= finder_->getSubElementContent(
					finder_->TAG_atomicNb_.asXMLString(), heavyIonElement);

	beamSettings[finder_->TAG_chargeNb_.asCString()]
			= finder_->getSubElementContent(
					finder_->TAG_chargeNb_.asXMLString(), heavyIonElement);

	beamSettings[finder_->TAG_nuclearCorre_.asCString()]
			= finder_->getSubElementContent(
					finder_->TAG_nuclearCorre_.asXMLString(), heavyIonElement);

	return beamSettings;
}

int VINHAC::XercescReader::PDGid(string particleName) {
	if (particleName == "d")
		return 1;
	else if (particleName == "dbar")
		return -1;
	else if (particleName == "u")
		return 2;
	else if (particleName == "ubar")
		return -2;
	else if (particleName == "s")
		return 3;
	else if (particleName == "sbar")
		return -3;
	else if (particleName == "c")
		return 4;
	else if (particleName == "cbar")
		return -4;
	else if (particleName == "b")
		return 5;
	else if (particleName == "bbar")
		return -5;
	else if (particleName == "t")
		return 6;
	else if (particleName == "tbar")
		return -6;
	else if (particleName == "e-" || particleName == "e")
		return 11;
	else if (particleName == "e+")
		return -11;
	else if (particleName == "mu-" || particleName == "mu")
		return 13;
	else if (particleName == "mu+")
		return -13;
	else if (particleName == "tau-" || particleName == "tau")
		return 15;
	else if (particleName == "tau+")
		return -15;
	else if (particleName == "gamma")
		return 22;
	else if (particleName == "Z0")
		return 23;
	else if (particleName == "W+")
		return 24;
	else if (particleName == "W-")
		return -24;
	else if (particleName == "nu_e")
		return 12;
	else if (particleName == "nu_ebar")
		return -12;
	else if (particleName == "nu_mu")
		return 14;
	else if (particleName == "nu_mubar")
		return -14;
	else if (particleName == "nu_tau")
		return 16;
	else if (particleName == "nu_taubar")
		return -16;
	else
		throw VinhacException("PDGid: wrong particle name");

}

vector<int> VINHAC::XercescReader::readInitialParticles() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::readInitialParticles()"<<std::endl;

#endif

	//get initial particle element
	DOMElement* initialParticlesElement = finder_->getElement(
			finder_->TAG_initialParticles_.asXMLString());
	//get particle list

	vector<int> particleVector;

	DOMNodeList* particleList = initialParticlesElement->getElementsByTagName(
			finder_->TAG_particleName_.asXMLString());

	for (unsigned i = 0; i < particleList->getLength(); i++) {
		DualString particleName = particleList->item(i)->getTextContent();
		int PDGid = this->PDGid(particleName.asCString());
		particleVector.push_back(PDGid);
	}

	return particleVector;
}

vector<int> VINHAC::XercescReader::readFinalParticles() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::readFinalParticles()"<<std::endl;

#endif

	//get initial particle element
	DOMElement* finalParticlesElement = finder_->getElement(
			finder_->TAG_finialParticles_.asXMLString());
	//get particle list

	vector<int> particleVector;

	DOMNodeList* particleList = finalParticlesElement->getElementsByTagName(
			finder_->TAG_particleName_.asXMLString());

	for (unsigned i = 0; i < particleList->getLength(); i++) {
		DualString particleName = particleList->item(i)->getTextContent();
		int PDGid = this->PDGid(particleName.asCString());
		particleVector.push_back(PDGid);
	}

	return particleVector;
}

vector<int> VINHAC::XercescReader::readIntermediateParticles() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::readIntermediateParticles()"<<std::endl;

#endif

	//get initial particle element
	DOMElement* intermediateParticlesElement = finder_->getElement(
			finder_->TAG_intermidiateBosons_.asXMLString());
	//get particle list

	vector<int> particleVector;

	DOMNodeList* particleList =
			intermediateParticlesElement->getElementsByTagName(
					finder_->TAG_particleName_.asXMLString());

	for (unsigned i = 0; i < particleList->getLength(); i++) {
		DualString particleName = particleList->item(i)->getTextContent();
		int PDGid = this->PDGid(particleName.asCString());
		particleVector.push_back(PDGid);
	}

	return particleVector;
}

double VINHAC::XercescReader::getParticleMass(int PDGid) {

	string mass = this->getParticleAttribute(PDGid, "m0");

	return strtod(mass.c_str(), NULL);

}

vector<int> VINHAC::XercescReader::readQuarksCharges() {
	vector<int> charges;
	charges.push_back(0);

	for (int i = 1; i <= 6; ++i) {
		charges.push_back(atoi(getParticleAttribute(i, "chargeType").c_str()));
	}

	return charges;
}

string VINHAC::XercescReader::getParticleAttribute(int PDGid,
		const char* propertyName) {

	// get Element particleDataBase
	DOMElement* particleDataBase = finderParticle_->getConfigElement();
	// get all particles as a node list
	DOMNodeList* nodelist = particleDataBase->getElementsByTagName(
			finderParticle_->TAG_particle_.asXMLString());
	////find particle with PDGid
	// loop over particles
	for (unsigned i = 0; i < nodelist->getLength(); i++) {

		// get particle id
		DualString id = nodelist->item(i)->getAttributes()->getNamedItem(
				DualString("id").asXMLString())->getTextContent();
		//cout << id.asCString() << endl;

		// return its property if PDGid
		if (atoi(id.asCString()) == PDGid) {
			DualString
					property =
							nodelist->item(i)->getAttributes()->getNamedItem(
									DualString(propertyName).asXMLString())->getTextContent();
			//	cout << property.asCString() << endl;

			string result(property.asCString());
			return result;
		}

	}
	std::cout << "There is no particle with id = " << PDGid
			<< "or required property in particleDataBase.xml file" << std::endl;
	return 0;

}

double VINHAC::XercescReader::getQ2minCut() {
	DOMElement* Q2maxCutelement = finder_->getElement(
			finder_->TAG_Q2minCut_.asXMLString());
	DualString Q2maxCut(Q2maxCutelement->getTextContent());
	char const * strc = Q2maxCut.asCString();//
	return strtol(strc, NULL, 10);
}

double VINHAC::XercescReader::getConstant(const char* constantName) {
	DOMElement* constMathElement = finder_->getElement(
			DualString(constantName).asXMLString());
	DualString constantContent(constMathElement->getTextContent());
	char const * strc = constantContent.asCString();//
	return strtod(strc, NULL);
}

/// Read IPS scheme
int VINHAC::XercescReader::getIPSscheme() {
	DOMElement* IPSelement = finder_->getElement(
			finder_->TAG_IPS_.asXMLString());
	DualString IPSscheme(IPSelement->getTextContent());
	char const * strc = IPSscheme.asCString();//
	return strtol(strc, NULL, 10);
}

int VINHAC::XercescReader::getWidthScheme() {
	DOMElement* IPSelement = finder_->getElement(
			finder_->TAG_widthSchem_.asXMLString());
	DualString IPSscheme(IPSelement->getTextContent());
	char const * strc = IPSscheme.asCString();//
	return strtol(strc, NULL, 10);
}

string VINHAC::XercescReader::getPolarization() {
	DOMElement* IPSelement = finder_->getElement(
			finder_->TAG_POLARIZATION.asXMLString());
	DualString IPSscheme(IPSelement->getTextContent());
	char const * strc = IPSscheme.asCString();//
	return string(strc);
}

int VINHAC::XercescReader::getFactorizationScale() {
	cout << "FAC2" << endl;
	DOMElement* Facelement = finder_->getElement(
			finder_->TAG_factorizationScale_.asXMLString());
	DualString Facscheme(Facelement->getTextContent());
	char const * strc = Facscheme.asCString();//
	cout << Facscheme << " FAC2 " << strc << endl;
	return strtol(strc, NULL, 10);
}

vector<vector<double> > VINHAC::XercescReader::readCKM() {
	vector<vector<double> > ckm;

	for (int i = 1; i < 4; ++i) {
		vector<double> tmp;
		for (int j = 1; j < 4; ++j) {
			stringstream ss;
			ss << i << j;
			DOMElement* ckmelement = finder_->getElement(DualString(string(
					finder_->TAG_CKM.asCString()) + ss.str()).asXMLString());
			DualString ckmel(ckmelement->getTextContent());
			tmp.push_back(atof(ckmel.asCString()));

		}
		ckm.push_back(tmp);
	}

#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::readCKM()"<<std::endl;
	for(unsigned i = 0;i<ckm.size(); ++i) {
		std::cout<<"DEBUG: XercescReader::readCKM() # ";
		for(unsigned j= 0; j<ckm[i].size(); ++j)
		std::cout<<ckm[i][j]<<" ";
		std::cout<<std::endl;

	}

#endif

	return ckm;
}

int VINHAC::XercescReader::getRadiativeCorrection() {
	DOMElement* ewCorr = finder_->getElement(
			finder_->TAG_ewCorrections_.asXMLString());
	string on(DualString(ewCorr->getAttribute(
			finder_->ATTR_switchOn_.asXMLString())).asCString());
	if (on != "on")
		return 0;

	DOMElement* correctionType = finder_->getElement(
			finder_->TAG_correctionType_.asXMLString());
	DualString corrType(correctionType->getTextContent());

	if (string(corrType.asCString()) == "1") {
		return 1;
	}

	if (string(corrType.asCString()) == "2") {
		return 2;
	}

	throw VinhacException("Not implemented radiative correction type");
}

VINHAC::XercescReader::~XercescReader() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::~XercescReader()"<<std::endl;

#endif
	delete parser_;
	delete finder_;
	delete finderParticle_;
	delete finderUserFile_;
	try {
		xercesc::XMLPlatformUtils::Terminate();
	} catch (xercesc::XMLException& e) {

		std::cerr << "XML toolkit teardown error: " << DualString(
				e.getMessage()) << std::endl;

	}
}
