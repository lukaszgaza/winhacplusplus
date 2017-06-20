/*
 * Commons.h
 *
 *  Created on: 25-08-2011
 *      Author: sobol
 */

#ifndef COMMONS_H_
#define COMMONS_H_

#include <string>
#include <utility>
#include "../core/enum/PolarizationDefinition.h"
#include "../core/enum/WeightsDefinition.h"

namespace VINHAC {

class Commons {
public:
	Commons();
	virtual ~Commons();
	static std::string mapWeightPolarizationToString(
			const std::pair<WeightName::name, PolarizationName::name>& target);
	static std::pair<WeightName::name,PolarizationName::name> mapStringToWeightPolarization(const std::string& target);
};

}

#endif /* COMMONS_H_ */
