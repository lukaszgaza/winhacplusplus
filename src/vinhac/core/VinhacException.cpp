/*
 * VinhacException.cpp
 *
 *  Created on: 2009-11-11
 *      Author: kamil
 */

#include "VinhacException.h"

namespace VINHAC {

VinhacException::VinhacException() {

}



VinhacException::~VinhacException() {
}

std::string VinhacException::getMessage(){
	return std::string(message);
}

}
