/*
 * XercescErrorHandler.h
 *
 *  Created on: 01-09-2011
 *      Author: sobol
 */

#ifndef XERCESCERRORHANDLER_H_
#define XERCESCERRORHANDLER_H_
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <iostream>
#include "helper-classes.h"
using namespace xercesc;
using namespace std;

namespace VINHAC {

class XercescErrorHandler: public ErrorHandler {
public:
	XercescErrorHandler();

	void warning(const SAXParseException& exc){
		cerr << "Error parsing file: " << DualString(exc.getMessage()) << endl;
	}

	void error(const SAXParseException& exc){
		cerr << "Error parsing file: " << DualString(exc.getMessage()) << endl;
	}

	void fatalError(const SAXParseException& exc){
		cerr << "Error parsing file: " << DualString(exc.getMessage()) << endl;
	}

	void resetErrors(){}

	virtual ~XercescErrorHandler();
};

}

#endif /* XERCESCERRORHANDLER_H_ */
