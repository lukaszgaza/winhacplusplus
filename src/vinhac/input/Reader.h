/*
 * FileReader.h
 *
 *  Created on: 2009-05-13
 *      Author: siudmy
 */

#ifndef FILEREADER_H_
#define FILEREADER_H_

#include<string>
#include<iostream>
#include<sstream>
#include<stdexcept>
#include<list>
#include<set>
#include<utility>

#include <vector>
#include <map>
#include "../core/enum/AlphaSchemeDefinition.h"
#include "../core/enum/WidthSchemeDefinition.h"
#include "../core/enum/PolarizationDefinition.h"
#include "../core/enum/FactorizationScaleDefinition.h"
#include "../core/enum/RadiativeCorrectionDefinition.h"
#include "../core/enum/FrameDefinition.h"
#include "../core/enum/PDFInterfaceDefinition.h"
#include "../core/enum/WeightsDefinition.h"
#include "../core/enum/PrintOutLevelDefinition.h"
using namespace std;

namespace VINHAC {
class Reader;
}

namespace VINHAC {
/**
 * \brief base class for classes handling input
 */
class Reader {

public:
	virtual void loadDefaultFile() = 0;
	virtual void loadUserFile(string path) = 0;
	virtual bool initialize() = 0;
	virtual map<string, string> readBeamSettings(const string& in) = 0;
	virtual bool isPdfOn() = 0;
	virtual PDFInterface::name getPdfInterface() = 0;
	virtual string getPdfName() = 0;
	virtual int getPdfSubset() = 0;
	virtual double getPdfXMin() = 0;
	virtual double getPdfXMax() = 0;
	virtual double getPdfQ2min() = 0;

	virtual vector<int> readInitialParticles() = 0;
	virtual void inputEcho(string path) = 0;
	virtual vector<int> readIntermediateParticles() = 0;
	virtual vector<int> readFinalParticles() = 0;
	virtual vector<int> readQuarksCharges() = 0;
	virtual double getParticleMass(int PDGid) = 0;
	virtual string getParticleAttribute(int PDGid, const char* propertyName) = 0;
	virtual double getConstant(const char* constantName) = 0;
	virtual vector<vector<double> > readCKM()=0;

	virtual AlphaScheme::name getAlphaScheme() = 0;
	virtual WidthScheme::name getWidthScheme() = 0;
	virtual FactorizationScale::name getFactorizationScale() = 0;
	virtual RadiativeCorrection::name getRadiativeCorrection() = 0;
	virtual double getSoftPhotonCutOff() = 0;
	virtual set<PolarizationName::name> getPolarizations() = 0;
	virtual FrameName::name getModelReferenceFrame() = 0;

	virtual pair<WeightName::name,PolarizationName::name> getMainWeight() = 0;
	virtual bool isWeightedEvent() = 0;
	virtual double getMaxWeightRejection() = 0;

	virtual string getPythiaXMLdocDir() = 0;
	virtual bool isPartonShower() = 0;
	virtual bool isHadronization() = 0;

	virtual PrintOutLevel::name getPrintOutLevel() = 0;

	virtual ~Reader() {
	}
	;
};
}
#endif /* FILEREADER_H_ */
