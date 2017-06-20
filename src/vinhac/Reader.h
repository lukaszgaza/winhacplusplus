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

#include <vector>
#include <map>
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
	virtual void loadDefaultFile() throw (runtime_error) = 0;
	virtual void loadUserFile(string path) throw (runtime_error) = 0;
	virtual bool initialize() = 0;
	virtual map<string, string> readBeamSettings(const string& in) = 0;
	virtual vector<int> readInitialParticles() = 0;
	virtual void inputEcho(string path) = 0;
	virtual vector<int> readIntermediateParticles() = 0;
	virtual vector<int> readFinalParticles() = 0;
	virtual vector<int> readQuarksCharges() = 0;
	virtual double getParticleMass(int PDGid) = 0;
	virtual string getParticleAttribute(int PDGid, const char* propertyName) = 0;
	virtual double getConstant(const char* constantName) = 0;
	virtual vector<vector<double> > readCKM()=0;

	virtual int getIPSscheme() = 0;
	virtual int getWidthScheme() = 0;
	virtual int getFactorizationScale() = 0;
	virtual int getRadiativeCorrection() = 0;
	virtual string getPolarization() = 0;
	virtual double getQ2minCut() = 0;
	virtual ~Reader() {
	}
	;
};
}
#endif /* FILEREADER_H_ */
