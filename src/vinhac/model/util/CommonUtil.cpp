/*
 * CommonUtil.cpp
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#include "CommonUtil.h"

namespace VINHAC {

CommonUtil::CommonUtil() {
}

CommonUtil::~CommonUtil() {
}


double VINHAC::CommonUtil::inverseWPropagator(const double& s) {
	double result;
	if (DataSource::get().keyWid == 0) {
		result = (s - (DataSource::get().MW) * (DataSource::get().MW)) * (s - (DataSource::get().MW) * (DataSource::get().MW))
				+ ((DataSource::get().MW) * (DataSource::get().GW)) * ((DataSource::get().MW) * (DataSource::get().GW));
	} else {
		result = (s - (DataSource::get().MW) * (DataSource::get().MW)) * (s - (DataSource::get().MW) * (DataSource::get().MW)) + (s * s
				* (DataSource::get().GW) * (DataSource::get().GW)) / ((DataSource::get().MW) * (DataSource::get().MW));
	}

	return result;
}

}
