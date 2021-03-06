/*
 * XercescReader.cpp
 *
 *  Created on: 2009-05-13
 *      Author: siudmy
 */
#include <vector>
#include <exception>
#include <map>
#include <vector>
using namespace std;

#ifndef __XercescReader_h__
#define __XercescReader_h__

#include "Reader.h"
#include<xercesc/util/PlatformUtils.hpp>

#include<xercesc/dom/DOM.hpp>
#include<xercesc/dom/DOMDocument.hpp>
#include<xercesc/dom/DOMDocumentType.hpp>
#include<xercesc/dom/DOMElement.hpp>
#include<xercesc/dom/DOMImplementation.hpp>
#include<xercesc/dom/DOMImplementationLS.hpp>
#include<xercesc/dom/DOMNodeIterator.hpp>
#include<xercesc/dom/DOMNodeList.hpp>
#include<xercesc/dom/DOMText.hpp>
#include<xercesc/util/XMLString.hpp>

#include<xercesc/parsers/XercesDOMParser.hpp>
#include<xercesc/util/XMLUni.hpp>

#include "helper-classes.h"
#include "ElementFinder.h"

#include "HardProcessHandler.h"
using namespace xercesc;

namespace VINHAC {
class Reader;
class XercescReader;
class ElementFinder;
}

namespace VINHAC {
/**
 * \brief class handling input provided by XML files.
 */
class XercescReader: public VINHAC::Reader {
private:
	string xmlFile_;
	string xmlFileParticle_;
	XercesDOMParser *parser_;
	ElementFinder *finder_;
	ElementFinder *finderParticle_;
	ElementFinder *finderUserFile_;

	void handleElement(const DOMElement* element) {
	}

public:
	XercescReader(const string& fromFile, const string& particle);
	void loadDefaultFile() throw (runtime_error);
	void loadUserFile(string path) throw (runtime_error);
	bool initialize();
	void inputEcho(string path);

	/*! Read beam setings
	 * \param beamName - beamA or beamB
	 * \return map with parameters (Energy, PDF settings ect) needed for
	 *  BeamParticle class setting
	 *********** \author Andrzej Siodmok \date Created on: 2009-06-10 ******/
	map<string, string> readBeamSettings(const string& beamName);

	/*! Read initial particles
	 * \return vector of PDG codes of initial particles
	 *********** \author Andrzej Siodmok \date Created on: 2009-06-10 ******/
	vector<int> readInitialParticles();

	/*! Read intermediate particles
	 * \return vector of PDG codes of intermediate particles
	 *********** \author Andrzej Siodmok \date Created on: 2009-06-10 ******/
	vector<int> readIntermediateParticles();

	/*! Read final particles
	 * \return vector of PDG codes of final particles
	 *********** \author Andrzej Siodmok \date Created on: 2009-06-10 ******/
	vector<int> readFinalParticles();

	vector<int> readQuarksCharges();

	/*! \return a PDGid of particle with name \param particleName
	 *********** \author Andrzej Siodmok \date Created on: 2009-06-10 ******/
	int PDGid(string particleName);

	double getParticleMass(int PDGid);
	string getParticleAttribute(int PDGid, const char* propertyName);

	double getConstant(const char* constantName);

	vector<vector<double> > readCKM();

	int getIPSscheme();
	int getWidthScheme();
	int getFactorizationScale();
	int getRadiativeCorrection();
	string getPolarization();
	double getQ2minCut();
	virtual ~XercescReader();
};

}

#endif
