/*
 * CommonUtil.h
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#ifndef COMMONUTIL_H_
#define COMMONUTIL_H_
#include "../../input/DataSource.h"

namespace VINHAC {

class CommonUtil {
public:
	CommonUtil();
	virtual ~CommonUtil();

	/**
	 * \brief provides value of inverse propagator of W boson
	 *
	 * Inverse of propagator is given by \f$ (s-m_W^2 )^2+(m_W \Gamma_W)^2 \f$
	 * @param s Mandelstam variable
	 * @param ds pointer to structure with input data
	 * @return value of inverse propagator
	 */
	static double inverseWPropagator(const double& s);
};

}

#endif /* COMMONUTIL_H_ */
