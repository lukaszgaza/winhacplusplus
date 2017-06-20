/*
 * Commons.cpp
 *
 *  Created on: 25-08-2011
 *      Author: sobol
 */

#include "Commons.h"

namespace VINHAC {

Commons::Commons() {
}

Commons::~Commons() {
}


	std::string Commons::mapWeightPolarizationToString(const std::pair<WeightName::name,PolarizationName::name>& target){
		return WeightName::mapEnumToString(target.first)+ "%" + PolarizationName::mapEnumToString(target.second);
	}

	std::pair<WeightName::name,PolarizationName::name> Commons::mapStringToWeightPolarization(const std::string& target){
		size_t found;

		found = target.find("%");

		if (found==std::string::npos) throw VinhacException("Tried to parse Weight-Polarization string, but not found delimiter %");

		std::string weight = target.substr(0,found);
		std::string polarization = target.substr(found+1);

		return std::make_pair(WeightName::mapStringToEnum(weight),PolarizationName::mapStringToEnum(polarization));
	}

}
