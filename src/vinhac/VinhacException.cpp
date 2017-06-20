/*
 * VinhacException.cpp
 *
 *  Created on: 2009-11-11
 *      Author: kamil
 */

#include "VinhacException.h"

namespace VINHAC {

VinhacException::VinhacException() {
	// TODO Auto-generated constructor stub

}



VinhacException::~VinhacException() {
	// TODO Auto-generated destructor stub
}

std::string VinhacException::getMessage(){
	return std::string(message);
}

}
