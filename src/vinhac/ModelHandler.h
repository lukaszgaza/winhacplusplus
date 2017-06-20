/*
 * ModelHandler.h
 *
 *  Created on: 26 Aug 2009
 *      Author: siudmy
 */

#ifndef MODELHANDLER_H_
#define MODELHANDLER_H_
#include "ProcessHandler.h"
#include "HardProcessHandler.h"
#include "../utils/SpecialMatrix.h"
#include "../utils/CLHEP/Vector/Vector/LorentzVector.h"

namespace VINHAC {
class Manager;
class Event;
class EventStep;
class Particle;
class ModelHandler;
class HardProcessHandler;
}

namespace VINHAC {

//!  Class responsible for model weights calculation
/*!
  This class calculates various model weights using proper matrix elements
*/
class ModelHandler: public VINHAC::ProcessHandler {
private:

	DataSource *ds; //!< pointer to object with input parameters

	VINHAC::Manager* p_Manager;
	VINHAC::HardProcessHandler* p_HardMatrix;

	double ckmElement;
	double e;
	double colourFactor;
	double normalizationFactor0;
	double normalizationFactor1;

	map<int, double> deltaWeak;

	void prepareParticles(VINHAC::Event& event);
	double matrixElement0(Event& event);
	double matrixElementFSR1(Event& event, Particle& photon);
	SpecialMatrix bornWproduction(Event& event);
	SpecialMatrix bornWdecay(Event& event);
	SpecialMatrix radiativeWdecay(Event& event, Particle& photon);

	double differentialXsectionCrude(double s, Event&);
	CLHEP::HepLorentzVector VecPol(CLHEP::HepLorentzVector p, int lambda);
	void HelEig(CLHEP::HepLorentzVector p, int hel, complex<double> chip[2]);
	void MxV2dC(complex<double> matrix[2][2], complex<double> vec[2], complex<
			double> res[2]);
	void
	MvtoWm(CLHEP::HepLorentzVector& p, complex<double> as[2][2], int alpha);
	complex<double> spinorialFunction(complex<double> hip1[2],
			CLHEP::HepLorentzVector a, complex<double> hip2[2], int ialfa);
	complex<double> spinorialFunction(complex<double> hip1[2],
			CLHEP::HepLorentzVector a, CLHEP::HepLorentzVector b,
			CLHEP::HepLorentzVector c, complex<double> hip2[2], int ialfa);

	//*********
	/*! \brief A method for calculating \f$\omega_{\pm}(p)=(E\pm |\vec{p}|)^{1/2}\f$, where
	 * \f$p^{\mu}=(E, \vec{p}) = (E, p_x, p_y, p_z)\f$.
	 ***
	 * \param E - energy
	 * \param pp - \f$|\vec{p}|^2 \f$
	 * \param imp = \f$\pm\f$
	 * \return \f$\omega\f$
	 *********** \author Andrzej Siodmok \date Created on: 2009-05-23 ******/
	inline double omega(int imp, double En, double pp) {
		return (sqrt(fabs(En + imp * pp)));
	}

public:
	/**
	 * \brief provides value of inverse propagator of W boson
	 *
	 * Inverse of propagator is given by \f$ (s-m_W^2 )^2+(m_W \Gamma_W)^2 \f$
	 * @param s Mandelstam variable
	 * @param ds pointer to structure with input data
	 * @return value of inverse propagator
	 */
	static double inverseWPropagator(double s, DataSource *ds);

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
	static double YWldec(CLHEP::HepLorentzVector Q, CLHEP::HepLorentzVector ql,
			double aKmax);


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
	static double ReBtug(double t, double xm1, double xm2, double amg);

	/**
	 * \brief Dilogarithm function.
	 *
	 * This is the CERNLIB version.
	 *
	 * @param x integration limit
	 * @return \f$ -\int_0^x \! \frac{\ln(1-z)}{z} \, dz \f$
	 */
	static double DILOGY(double x);


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
	static double Btilde(CLHEP::HepLorentzVector p1,
			CLHEP::HepLorentzVector p2, double am1, double am2, double aKmax,
			double amgam);

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
	static double A4anal(double p1p2, double E1, double E2, double am1,
			double am2);

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
	static double A4Q2pos(double p1p2, double En1, double En2, double xm1,
			double xm2);

	/**
	 * \brief Some auxiliary function (combination of Logs and Dilogs) used in
	 * the function A4anal for \f$A_4(p_1,p_2)\f$.
	 *
	 */
	static double Yijeta(double yi, double yj, double eta);

	/**
	 * Function
     * @return \f$ \ln| \frac{(z-x_1)(z-x_2)}{(z-x_3)(z-x_4)} |\f$
	 */
	inline static double etaln(double x1, double x2, double x3, double x4,
			double z) {
		return log(abs((z - x1) * (z - x2) / (z - x3) / (z - x4)));
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
	static double A4Q2neg(double p1p2, double En1, double En2, double xm1,
			double xm2);

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
	static double A4Q2fRF(double E1, double am1, double am2);

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
	static double A4Q2zRF(double am1, double am2);

	/**
	 * \brief Virtual QED radiative corrections to the YFS beta0 function
	 * for leptonic W-decays.
	 *
	 * @param amW W-boson mass
	 * @param aml lepton mass
	 * @return function value
	 */
	static double deltaQED(double amW, double aml);

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
	static double StPair(double am1s, double am2s, CLHEP::HepLorentzVector p1,
			CLHEP::HepLorentzVector p2, CLHEP::HepLorentzVector pk);

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
	static double VirSofQED(double amW, double aml, double aKmax);

	/**
	* \brief evolves the event
	*
	* It does whole part of generation process which is to do in this handler
	* @see Event
	* @param event reference to current Event object
	* @return if weight was zero
	*/
	bool eventEvolve(VINHAC::Event& event);

	/**
	 * \brief constructor
	 *
	 * @param ds pointer to structure with input data
	 * @param p_Hard pointer to HadrProcessHandler object
	 */
	ModelHandler(DataSource *ds, VINHAC::HardProcessHandler* p_Hard);
	~ModelHandler() {
	}


};
}

#endif /* MODELHANDLER_H_ */
