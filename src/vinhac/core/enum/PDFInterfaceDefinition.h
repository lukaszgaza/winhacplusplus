/*
 * FrameDefinition.h
 *
 *  Created on: 2010-12-09
 *      Author: kamil
 */

#ifndef PDFINTERFACEDEFINITION_H_
#define PDFINTERFACEDEFINITION_H_

#include "../VinhacException.h"

namespace VINHAC {


struct PDFInterface {
	enum name {
		LHAPDF = 0,
		Simple = 1
	};

	static name mapStringToEnum(const std::string& target){
		if(target=="LHAPDF"){
			return LHAPDF;
		} else if(target=="Simple"){
			return Simple;
		}
		throw VinhacException("Unknown PDFInterface name");
	}

	static std::string mapEnumToString(const name& target){
		switch(target){
			case LHAPDF:						return "LHAPDF";
			case Simple:						return "Simple";
		}
		throw VinhacException("Unknown PDFInterface name");
	}
};

}

#endif /* PDFINTERFACEDEFINITION_H_ */
