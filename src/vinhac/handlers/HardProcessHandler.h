#include <string>
using namespace std;

#ifndef __HardProcessHandler_h__
#define __HardProcessHandler_h__

#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include "../../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include "../../utils/CLHEP/Random/Random/RanluxEngine.h"
#include <complex>
#include "ProcessHandler.h"

enum Fermions {
	quark_d = 1,
	quark_u = 2,
	quark_s = 3,
	quark_c = 4,
	quark_b = 5,
	quark_t = 6,
	electron = 11,
	nu_e = 12,
	muon = 13,
	nu_mu = 14,
	tau = 15,
	nu_tau = 16,
};

enum Bosons {
	photon = 22, Z = 23, W = 24,
};

namespace VINHAC {
class DataSource;
class Manager;
class Event;
class HardProcessHandler;
class ModelHandler;
class ProcessHandler;
}

namespace VINHAC {

//!  Class resposible for parton level process handling.
/*!
 This class generates leptons and photons. It calculates also kinematics of the process.
 */
class HardProcessHandler: public ProcessHandler {
private:
	VINHAC::Manager* _unnamed_Manager;

	TRND * p_randomGenerator; //! < pointer to the random number generator


	CLHEP::HepLorentzVector sumOfPhotonsMomentum;

	double phi;
	double theta;

	double Qe;

	double gW; //!< Coupling constant \f$g_W = \frac{1}{\sin{\Theta_W}\cos{\Theta_W}}\f$  where \f$\Theta_W\f$ is the Weinberg angle or weak mixing angle
	double Ecm; //!< CMS energy of incoming hadron beams 14e3 Ecm for LHC sum of two beams

	int KeyWmp; //!< Intermediate bosons choice (Z, gamma, Z+gamma, other)
	int KeyQua; //!< Quark doublet choice
	int KeyWid; //!< Typ of Z width: fixed or running


	double x_min, x_max; //!< x range
	double Q2_min; //!< minimum of Q^2 for PDF

	double epsk;
	double delta;

	//sudakov variables
	vector<double> ygr;
	vector<double> zet;

	double prepareParticles(VINHAC::Event&);
	double generateLepton(VINHAC::Event&);
	void generateLeptonAngles(Event&);
	double generateLeptonMomenta(VINHAC::Event&);
	void transformationToLAB(VINHAC::Event&);

	double generateMultiphotonRadiation(Event& event);
	double YWlRem(Event& event);
	void generatePhotonAngularDistribution(const double& aM, const double& aml,
			const double& amx2, const double& beta, double& costhg,
			double& sinthg, double& dist);
	double setLeptonMomentaInRadiativeDecay(Event& event, const double& sprim,
			const double& aM, const CLHEP::HepLorentzVector& phSum);
	double yfsCut(const double& am1, const double& am2,
			const CLHEP::HepLorentzVector& p1,
			const CLHEP::HepLorentzVector& p2, const double& epsg);
	/*!
	 // \brief  This generates photon multipl. nphot according to poisson distr.
	 // \param  average = average multiplicity
	 // \param        nmax  = maximum multiplicity
	 // \return  mult = generated multiplicity
	 // \return  vector rr list of ordered uniform random numbers,
	 //           a byproduct result, to be eventually used for some further
	 //           purpose (i.e.  generation of photon energies).				*/

	int PoissGenerator(double average, int nmax, vector<double> &rr);

public:

	/**
	 * \brief constructor
	 *
	 * It does the initialization of HardProcessHandler
	 * @param ds pointer to structure with input data
	 */
	HardProcessHandler();

	~HardProcessHandler();

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
	 * \brief returns gW
	 *
	 * @return gW
	 */
	inline double getGw() {
		return this->gW;
	}

	/**
	 * \brief sets poniter to random number generator
	 *
	 * @see TRND
	 * @param p_Generator pointer to random number generator
	 */
	inline void setRandomNumberGenerator(TRND *p_Generator) {
		p_randomGenerator = p_Generator;
	}

};
}

#endif
