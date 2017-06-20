/*
 * VectorUtil.h
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#ifndef VECTORUTIL_H_
#define VECTORUTIL_H_

#include <complex>
#include "../../../utils/CLHEP/Vector/Vector/LorentzVector.h"

using namespace std;
namespace VINHAC {

class VectorUtil {
public:
	VectorUtil();
	virtual ~VectorUtil();

	static void VecPol(CLHEP::HepLorentzVector& epsilonW,
			const CLHEP::HepLorentzVector& p, const int& lambda);
	static void HelEig(const CLHEP::HepLorentzVector& p, const int& hel, complex<
			double> chip[2]);
	static void MxV2dC(complex<double> matrix[2][2], complex<double> vec[2], complex<
			double> res[2]);
	static void
	MvtoWm(const CLHEP::HepLorentzVector& p, complex<double> as[2][2],
			const int& alpha);

	//*********
	/*! \brief A method for calculating \f$\omega_{\pm}(p)=(E\pm |\vec{p}|)^{1/2}\f$, where
	 * \f$p^{\mu}=(E, \vec{p}) = (E, p_x, p_y, p_z)\f$.
	 ***
	 * \param E - energy
	 * \param pp - \f$|\vec{p}|^2 \f$
	 * \param imp = \f$\pm\f$
	 * \return \f$\omega\f$
	 *********** \author Andrzej Siodmok \date Created on: 2009-05-23 ******/
	static inline double omega(const int& imp, const double& En, const double& pp) {
		return (sqrt(fabs(En + imp * pp)));
	}
};

}

#endif /* VECTORUTIL_H_ */
