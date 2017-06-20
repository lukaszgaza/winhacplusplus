/*
 * YFSUtil.h
 *
 *  Created on: 22-08-2011
 *      Author: sobol
 */

#ifndef YFSUTIL_H_
#define YFSUTIL_H_

#include "../../../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include <cstdlib>

namespace VINHAC {

class YFSUtil {
public:
	YFSUtil();
	virtual ~YFSUtil();

	/**
	 * \brief This function provides a value of YFS infrared function Btilde+ReB
	 * for leptonic W-decays.
	 *
	 * \f$ W(Q) \rightarrow l(q_l) + \nu(q_{\nu}) \f$
	 * @param Q 4-momentum of W
	 * @param ql 4-momenta of charged lepton
	 * @param aKmax soft photon cut-off in GeV
	 * @return YFS infrared function
	 */
	static double YWldec(const CLHEP::HepLorentzVector& Q,
			const CLHEP::HepLorentzVector& ql, const double& aKmax);

	/**
	 * \brief The t or u-channel virtual infrared YFS function ReB corresponding
	 * to any pair of massive charged particles participating in scattering,
	 * production or neutral-particle decay processes.
	 *
	 * @param t any \f$ (p_1 - p_2)^2 \f$, where \f$p_1\f$, \f$p_2\f$ are four-momenta of the initial and final state particles, respectively
	 * @param xm1 mass of particle
	 * @param xm2 mass of particle
	 * @param amg dummy "photon mass" (IR regulator)
	 * @return The t or u-channel virtual infrared YFS function ReB
	 */
	static double ReBtug(const double& t, const double& xm1, const double& xm2,
			const double& amg);

	/**
	 * \brief Dilogarithm function.
	 *
	 * This is the CERNLIB version.
	 *
	 * @param x integration limit
	 * @return \f$ -\int_0^x \! \frac{\ln(1-z)}{z} \, dz \f$
	 */
	static double DILOGY(const double& x);

	/**
	 * \brief This function provides a value of YFS real photon IR function
	 * B-tilde for any pair of charged particles.
	 *
	 * @param p1 particle 4-momenta
	 * @param p2 particle 4-momenta
	 * @param am1 particles masses
	 * @param am2 particles masses
	 * @param amgam "photon mass" (IR regulator)
	 * @param aKmax maximum soft photon energy [GeV]
	 * @return function value
	 */
	static double Btilde(const CLHEP::HepLorentzVector& p1,
			const CLHEP::HepLorentzVector& p2, const double& am1,
			const double& am2, const double& aKmax, const double& amgam);

	/**
	 * \brief This function provides an analytical result for the integral
	 * \f$A_4(p_1,p_2)\f$ being a part of the YFS IR function B-tilde.
	 *
	 * This is a general case without any approximation!\n
	 * Added option for \f$Q^2=0\f$ (e.g. leptonic W-decay)\n
	 * Added option for calculation in particle rest frame\n
	 * Added version for \f$Q^2<0\f$ (numerically safe)\n
	 * @param p1p2 scalar product of the 4-momenta \f$p_1\f$ and \f$p_2\f$
	 * @param E1 particle energy
	 * @param E2 particle energy
	 * @param am1 particle mass
	 * @param am2 particle mass
	 * @return function value
	 */
	static double A4anal(const double& p1p2, const double& E1,
			const double& E2, const double& am1, const double& am2);

	/**
	 * \brief This function provides an analytical result for the integral
	 * \f$A_4(p_1,p_2)\f$ being a part of the YFS IR function B-tilde for the usual
	 * case, i.e. \f$Q^2 > 0\f$ and NOT in \f$p_1\f$/\f$p_2\f$ Rest Frame.
	 *
	 * @param p1p2 scalar product of the 4-momenta \f$p_1\f$ and \f$p_2\f$
	 * @param En1 particle energy
	 * @param En2 particle energy
	 * @param xm1 particle mass
	 * @param xm2 particle mass
	 * @return function value
	 */
	static double A4Q2pos(const double& p1p2, double En1, double En2,
			double xm1, double xm2);

	/**
	 * \brief Some auxiliary function (combination of Logs and Dilogs) used in
	 * the function A4anal for \f$A_4(p_1,p_2)\f$.
	 *
	 */
	static double Yijeta(const double& yi, const double& yj, const double& eta);

	/**
	 * Function
	 * @return \f$ \ln| \frac{(z-x_1)(z-x_2)}{(z-x_3)(z-x_4)} |\f$
	 */
	inline static double etaln(const double& x1, const double& x2,
			const double& x3, const double& x4, const double& z) {
		return log(fabs((z - x1) * (z - x2) / (z - x3) / (z - x4)));
	}

	/**
	 * \brief This function provides an analytical result for the integral
	 * \f$A_4(p_1,p_2)\f$ being a part of the YFS IR function B-tilde for the usual
	 * case, i.e. \f$Q^2 < 0\f$ and NOT in \f$p_1\f$/\f$p_2\f$ Rest Frame.
	 * @param p1p2 scalar product of the 4-momenta \f$p_1\f$ and \f$p_2\f$
	 * @param En1 particle energy
	 * @param En2 particle energy
	 * @param xm1 particle mass
	 * @param xm2 particle mass
	 * @return function value
	 *
	 */
	static double A4Q2neg(const double& p1p2, double En1, double En2,
			double xm1, double xm2);

	/**
	 * \brief This function provides an analytical result for the integral
	 * \f$A_4(p_1,p_2)\f$ being a part of the YFS IR function B-tilde for the case:
	 * \f$Q^2 = 0\f$ (e.g. leptonic W deacys) but NOT in \f$p_1\f$/\f$p_2\f$ Rest Frame.
	 * @param E1 particle energy
	 * @param E2 particle energy
	 * @param am1 particle mass
	 * @param am2 particle mass
	 * @return function value
	 */
	static double A4Q2zer(double E1, double E2, double am1, double am2);

	/**
	 * \brief This function provides an analytical result for the integral
	 * \f$A_4(p_1,p_2)\f$ being a part of the YFS IR function B-tilde for \f$Q^2 \neq 0\f$.
	 *
	 * NOTE: Calculation is performed in the \f$p_2\f$ Rest Frame!
	 * @param E1 energy of \f$p_1\f$;
	 * @param am1 particle mass
	 * @param am2 particle mass
	 * @return function value
	 */
	static double A4Q2fRF(const double& E1, const double& am1,
			const double& am2);

	/**
	 * \brief This function provides an analytical result for the integral
	 * \f$A_4(p_1,p_2)\f$ being a part of the YFS IR function B-tilde
	 * for the case: \f$Q^2 = 0\f$ (e.g. leptonic W deacys).
	 *
	 * NOTE: Calculation is performed in the \f$p_2\f$ Rest Frame!
	 *
	 * @param am1 particle mass
	 * @param am2 particle mass
	 * @return function variable
	 */
	static double A4Q2zRF(const double& am1, const double& am2);

	/**
	 * \brief Virtual QED radiative corrections to the YFS beta0 function
	 * for leptonic W-decays.
	 *
	 * @param amW W-boson mass
	 * @param aml lepton mass
	 * @return function value
	 */
	static double deltaQED(const double& amW, const double& aml);

	/**
	 * \brief S-tilde factor for any pair of charged particles.
	 *
	 * @param am1s squared mass of particle
	 * @param am2s squared mass of particle
	 * @param p1 particle's 4-momenta
	 * @param p2 particle's 4-momenta
	 * @param pk photon 4-momentum
	 * @return function value
	 */
	static double StPair(const double& am1s, const double& am2s,
			const CLHEP::HepLorentzVector& p1,
			const CLHEP::HepLorentzVector& p2,
			const CLHEP::HepLorentzVector& pk);

	/**
	 * \brief Virtual + soft-real QED radiative corrections for leptonic W-decays
	 * in W rest frame
	 *
	 * (Fleischer & Jegerlehner, Z. Phys. C26 (1985) 629).
	 *
	 * @param amW W-boson mass
	 * @param aml lepton mass
	 * @param aKmax maximum soft real-photon energy in GeV
	 * @return function value
	 */
	static double VirSofQED(const double& amW, const double& aml,
			const double& aKmax);
};

}

#endif /* YFSUTIL_H_ */
