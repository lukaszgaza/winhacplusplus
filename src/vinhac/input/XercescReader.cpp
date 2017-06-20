/*
 * XercescReader.cpp
 *
 *  Created on: 2009-05-13
 *      Author: siudmy
 */
#include <iostream>
#include <sstream>
#include "XercescReader.h"
#include "XercescErrorHandler.h"
#include "../core/VinhacException.h"
#include <xercesc/framework/psvi/XSValue.hpp>

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
	parser_->setErrorHandler(new XercescErrorHandler());
	finder_ = new ElementFinder();
	finderParticle_ = new ElementFinder();
	finderUserFile_ = new ElementFinder();

}

void VINHAC::XercescReader::loadDefaultFile()  {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::loadDefaultFile()"<<std::endl;

#endif

	parser_->setValidationScheme(XercesDOMParser::Val_Always);
	parser_->setDoNamespaces(true);
	parser_->setDoSchema(true);
	parser_->setValidationSchemaFullChecking(true);

	try {
		parser_->parse(xmlFile_.c_str());
		if(parser_->getErrorCount()>0){
			throw VinhacException(xmlFile_+" does not match InputData.xsd");
		}

		// there's no need to free this pointer -- it's
		// owned by the parent parser object
		DOMDocument* xmlDoc = parser_->getDocument();

		if (NULL == xmlDoc->getDocumentElement()) {
			throw VinhacException("empty XML document");
		}

		finder_->setDocument(xmlDoc);

		parser_->setValidationScheme(XercesDOMParser::Val_Never);
		parser_->parse(xmlFileParticle_.c_str());

		DOMDocument* xmlParticle = parser_->getDocument();

		if (NULL == xmlParticle->getDocumentElement()) {
			throw(runtime_error("empty XML document"));
		}

		finderParticle_->setDocument(xmlParticle);


	} catch (const VinhacException& e){
		throw e;
	} catch (const SAXException& e) {

			// believe it or not, XMLException is not
			// a parent class of DOMException

			cerr << "Error parsing file: " << DualString(e.getMessage()) << endl;

			throw VinhacException(DualString(e.getMessage()).asCString());

	} catch (const XMLException& e) {

		// believe it or not, XMLException is not
		// a parent class of DOMException

		cerr << "Error parsing file: " << DualString(e.getMessage()) << endl;

		throw VinhacException(DualString(e.getMessage()).asCString());

	} catch (const DOMException& e) {

		cerr << "Encountered DOM Exception: " << DualString(e.getMessage())
				<< endl;
		throw VinhacException(DualString(e.getMessage()).asCString());
	} catch (...) {
		throw VinhacException("Unexpected error during parsing input data, probably xml does not match InputData.xsd");
	}
	return;
} // load()

void VINHAC::XercescReader::loadUserFile(string userFile)  {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: XercescReader::loadUserFile()"<<std::endl;

#endif

	parser_->setValidationScheme(XercesDOMParser::Val_Always);
	parser_->setDoNamespaces(true);
	parser_->setDoSchema(true);
	parser_->setValidationSchemaFullChecking(true);

	try {

		parser_->parse(userFile.c_str());
		if(parser_->getErrorCount()>0){
			throw VinhacException(userFile+" does not match InputData.xsd");
		}

		// there's no need to free this pointer -- it's
		// owned by the parent parser object
		DOMDocument* xmlUserDoc = parser_->getDocument();

		if (NULL == xmlUserDoc->getDocumentElement()) {
			throw(runtime_error("empty XML document"));
		}

		finder_->setUserDocument(xmlUserDoc);
	} catch (const VinhacException& e){
		throw e;
	}  catch (XMLException& e) {

		// believe it or not, XMLException is not
		// a parent class of DOMException

		cerr << "Error parsing file: " << DualString(e.getMessage()) << endl;

		throw VinhacException(DualString(e.getMessage()).asCString());

	} catch (const DOMException& e) {

		cerr << "Encountered DOM Exception: " << DualString(e.getMessage())
				<< endl;
		throw VinhacException(DualString(e.getMessage()).asCString());
	} catch (...){
		throw VinhacException("Unexpected error during parsing input data, probably xml does not match InputData.xsd");
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


	//Beam element
	DOMNodeList* elements = finder_->getElements(finder_->TAG_beam.asXMLString());

	DOMElement* element = 0;
	for(unsigned i = 0; i< elements->getLength(); ++i){
		if (DOMNode::ELEMENT_NODE == elements->item(i)->getNodeType()) {
			element = dynamic_cast<DOMElement*> (elements->item(i));
		}
		if(element!=0){
			if(smBeam.convert(element->getAttribute(finder_->ATTR_beam_id.asXMLString())) == beamName){
				break;
			} else {
				element = 0;
			}
		}
	}
	if(element == 0){
		throw VinhacException("Beam settings for beam id="+beamName+" not found please check the input");
	}

	//Energy of the beam
	beamSettings[finder_->TAG_energy.asCString()] = finder_->getSubElementContent(finder_->TAG_energy.asXMLString(),element);


	beamSettings[finder_->TAG_PDGid.asCString()] = finder_->getSubElementContent(finder_->TAG_PDGid.asXMLString(),element);

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
			finder_->TAG_quarks.asXMLString());
	//get particle list

	vector<int> particleVector;

	DOMNodeList* particleList = initialParticlesElement->getElementsByTagName(
			finder_->TAG_name.asXMLString());

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
			finder_->TAG_leptons.asXMLString());
	//get particle list

	vector<int> particleVector;

	DOMNodeList* particleList = finalParticlesElement->getElementsByTagName(
			finder_->TAG_name.asXMLString());

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
			finder_->TAG_bosons.asXMLString());
	//get particle list

	vector<int> particleVector;

	DOMNodeList* particleList =
			intermediateParticlesElement->getElementsByTagName(
					finder_->TAG_name.asXMLString());

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
			finderParticle_->TAG_particle.asXMLString());
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


double VINHAC::XercescReader::getConstant(const char* constantName) {
	StringManager sm;
	DOMElement* constElement = 0;
	DOMNodeList* constList = finder_->getElements(
			finder_->TAG_constant.asXMLString());
	for(unsigned i = 0 ; i < constList->getLength(); ++i){
		if (DOMNode::ELEMENT_NODE == constList->item(i)->getNodeType()) {
			constElement = dynamic_cast<DOMElement*> (constList->item(i));
		}
		if(constElement!=0){
			if(string(sm.convert(constElement->getAttribute(finder_->ATTR_constant_name.asXMLString()))) == string(constantName)){
				break;
			} else {
				constElement = 0;
			}
		}
	}
	if(constElement!=0){
		DualString constantContent(constElement->getTextContent());
		char const * strc = constantContent.asCString();//
		return strtod(strc, NULL);
	} else {
		throw VinhacException("constant with name="+string(constantName)+" not found");
	}
}

/// Read Alpha Scheme
VINHAC::AlphaScheme::name VINHAC::XercescReader::getAlphaScheme() {
	DOMElement* IPSelement = finder_->getElement(
			finder_->TAG_alphaScheme.asXMLString());
	DualString IPSscheme(IPSelement->getTextContent());
	string strc(IPSscheme.asCString());//
	return AlphaScheme::mapStringToEnum(strc);
}

VINHAC::WidthScheme::name VINHAC::XercescReader::getWidthScheme() {
	DOMElement* element = finder_->getElement(
			finder_->TAG_widthScheme.asXMLString());
	DualString scheme(element->getTextContent());
	string strc(scheme.asCString());//
	return WidthScheme::mapStringToEnum(strc);
}

set<VINHAC::PolarizationName::name> VINHAC::XercescReader::getPolarizations() {
	set<VINHAC::PolarizationName::name> result;
	StringManager sm;
	DOMElement* polarizationsElement = 0;
	DOMElement* polarizationElement = 0;

	polarizationsElement = finder_->getElement(finder_->TAG_polarizations.asXMLString());
	DOMNodeList* polarizationList = finder_->getSubElements(
			finder_->TAG_polarization.asXMLString(),polarizationsElement);
	for(unsigned i = 0 ; i < polarizationList->getLength(); ++i){
		polarizationElement = 0;
		if (DOMNode::ELEMENT_NODE == polarizationList->item(i)->getNodeType()) {
			polarizationElement = dynamic_cast<DOMElement*> (polarizationList->item(i));
		}
		if(polarizationElement!=0){
			DualString polarization(polarizationElement->getTextContent());
			string pol(polarization.asCString());
			result.insert(PolarizationName::mapStringToEnum(pol));
		}
	}

	if(result.size()==0){
		throw VINHAC::VinhacException("No polarization has been set ! check the input files");
	}
	return result;
}

VINHAC::FactorizationScale::name VINHAC::XercescReader::getFactorizationScale() {
	DOMElement* element = finder_->getElement(
			finder_->TAG_factorizationScale.asXMLString());
	DualString scheme(element->getTextContent());
	string strc(scheme.asCString());//
	return FactorizationScale::mapStringToEnum(strc);
}

vector<vector<double> > VINHAC::XercescReader::readCKM() {
	StringManager sm;
	DOMElement* ckmElement = 0;
	DOMNodeList* ckmList = finder_->getElements(
			finder_->TAG_CKM.asXMLString());


	vector<vector<double> > ckm;

	for (int i = 1; i < 4; ++i) {
		vector<double> tmp;
		for (int j = 1; j < 4; ++j) {
			ostringstream ss_i;
			ss_i << i;
			string str_i = ss_i.str();
			ostringstream ss_j;
			ss_j << j;
			string str_j = ss_j.str();
			for(unsigned k = 0 ; k < ckmList->getLength(); ++k){
					ckmElement = 0;
					if (DOMNode::ELEMENT_NODE == ckmList->item(k)->getNodeType()) {
						ckmElement = dynamic_cast<DOMElement*> (ckmList->item(k));
					}
					if(ckmElement!=0){
						if(string(sm.convert(ckmElement->getAttribute(finder_->ATTR_CKM_i.asXMLString()))) == str_i
								&& string(sm.convert(ckmElement->getAttribute(finder_->ATTR_CKM_j.asXMLString()))) == str_j){
							break;
						} else {
							ckmElement = 0;
						}
					}
			}
			if(ckmElement!=0){
				DualString ckmel(ckmElement->getTextContent());
				tmp.push_back(atof(ckmel.asCString()));
			} else {
				throw VINHAC::VinhacException("Missing CKM element i= "+str_i+" j="+str_j);
			}



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

VINHAC::RadiativeCorrection::name VINHAC::XercescReader::getRadiativeCorrection() {
	DOMElement* radiation = finder_->getElement(
			finder_->TAG_radiation.asXMLString());
	string level = finder_->getSubElementContent(finder_->TAG_level.asXMLString(),radiation);

	return RadiativeCorrection::mapStringToEnum(level);
}

VINHAC::FrameName::name VINHAC::XercescReader::getModelReferenceFrame(){
	DOMElement* element = finder_->getElement(
			finder_->TAG_modelReferenceFrame.asXMLString());
	DualString scheme(element->getTextContent());
	string strc(scheme.asCString());//
	return FrameName::mapStringToEnum(strc);
}

bool VINHAC::XercescReader::isPdfOn(){
	const XMLCh *text = finder_->getElement(finder_->TAG_PDFs.asXMLString())->getAttribute(finder_->ATTR_PDFs_on.asXMLString());
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_boolean,
	status);
	bool result = val -> fData.fValue.f_bool;

	delete val;
	return result;
}
VINHAC::PDFInterface::name VINHAC::XercescReader::getPdfInterface(){
	DOMElement* element = finder_->getElement(
			finder_->TAG_PDFs.asXMLString());
	string strc = finder_->getSubElementContent(finder_->TAG_interface.asXMLString(),element);
	return PDFInterface::mapStringToEnum(strc);
}
string VINHAC::XercescReader::getPdfName(){
	DOMElement* element = finder_->getElement(
			finder_->TAG_PDFs.asXMLString());
	string strc = finder_->getSubElementContent(finder_->TAG_name.asXMLString(),element);
	return strc;
}
int VINHAC::XercescReader::getPdfSubset(){
	DOMElement* element = finder_->getElement(finder_->TAG_PDFs.asXMLString());
	element = finder_->getSubElement(finder_->TAG_subset.asXMLString(),element);
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_integer,
	status);
	int result = val -> fData.fValue.f_int;

	delete val;
	return result;
}
double VINHAC::XercescReader::getPdfXMin(){
	DOMElement* element = finder_->getElement(finder_->TAG_PDFs.asXMLString());
	element = finder_->getSubElement(finder_->TAG_xMin.asXMLString(),element);
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_double,
	status);
	double result = val -> fData.fValue.f_doubleType.f_double;

	delete val;
	return result;
}
double VINHAC::XercescReader::getPdfXMax(){
	DOMElement* element = finder_->getElement(finder_->TAG_PDFs.asXMLString());
	element = finder_->getSubElement(finder_->TAG_xMax.asXMLString(),element);
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_double,
	status);
	double result = val -> fData.fValue.f_doubleType.f_double;

	delete val;
	return result;
}
double VINHAC::XercescReader::getPdfQ2min(){
	DOMElement* element = finder_->getElement(finder_->TAG_PDFs.asXMLString());
	element = finder_->getSubElement(finder_->TAG_Q2min.asXMLString(),element);
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_double,
	status);
	double result = val -> fData.fValue.f_doubleType.f_double;

	delete val;
	return result;
}

double VINHAC::XercescReader::getSoftPhotonCutOff(){
	DOMElement* element = finder_->getElement(finder_->TAG_softPhotonCutOff.asXMLString());
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_double,
	status);
	double result = val -> fData.fValue.f_doubleType.f_double;

	delete val;
	return result;
}

bool VINHAC::XercescReader::isPartonShower(){
	DOMElement* element = finder_->getElement(finder_->TAG_partonShower.asXMLString());
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_boolean,
	status);
	bool result = val -> fData.fValue.f_bool;

	delete val;
	return result;
}

bool VINHAC::XercescReader::isHadronization(){
	DOMElement* element = finder_->getElement(finder_->TAG_hadronization.asXMLString());
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_boolean,
	status);
	bool result = val -> fData.fValue.f_bool;

	delete val;
	return result;
}

pair<VINHAC::WeightName::name,VINHAC::PolarizationName::name> VINHAC::XercescReader::getMainWeight(){
	DOMElement* element = finder_->getElement(
			finder_->TAG_main.asXMLString());
	string weight = finder_->getSubElementContent(finder_->TAG_modelWeight.asXMLString(),element);
	string polarization = finder_->getSubElementContent(finder_->TAG_polarization.asXMLString(),element);
	return make_pair(WeightName::mapStringToEnum(weight),PolarizationName::mapStringToEnum(polarization));
}

bool VINHAC::XercescReader::isWeightedEvent(){
	DOMElement* element = finder_->getElement(finder_->TAG_weightedEvents.asXMLString());
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_boolean,
	status);
	bool result = val -> fData.fValue.f_bool;

	delete val;
	return result;
}

double VINHAC::XercescReader::getMaxWeightRejection(){
	DOMElement* element = finder_->getElement(finder_->TAG_maxWeightRejection.asXMLString());
	const XMLCh *text = element->getTextContent();
	XSValue::Status status;
	  // Extract the XMLCh from the node
	XSValue *val = XSValue::getActualValue(text, XSValue::dt_double,
	status);
	double result = val -> fData.fValue.f_doubleType.f_double;

	delete val;
	return result;
}

VINHAC::PrintOutLevel::name VINHAC::XercescReader::getPrintOutLevel(){
	DOMElement* element = finder_->getElement(
			finder_->TAG_technical.asXMLString());
	string strc = finder_->getSubElementContent(finder_->TAG_printOutLevel.asXMLString(),element);
	return PrintOutLevel::mapStringToEnum(strc);
}

string VINHAC::XercescReader::getPythiaXMLdocDir(){
	DOMElement* element = finder_->getElement(
			finder_->TAG_pythiaXMLdocDir.asXMLString());
	DualString scheme(element->getTextContent());
	string strc(scheme.asCString());//
	return strc;
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
