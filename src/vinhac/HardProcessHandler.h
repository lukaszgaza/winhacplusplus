#include <string>
using namespace std;

#ifndef __HardProcessHandler_h__
#define __HardProcessHandler_h__

#include "VinhacException.h"
#include "Event.h"
#include "ProcessHandler.h"
#include "DataSource.h"
#include "ModelHandler.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include "../utils/CLHEP/Vector/Vector/LorentzVector.h"
#include "../utils/CLHEP/Random/Random/RanluxEngine.h"
#include <complex>

enum Fermions {
	d = 1,
	u = 2,
	s = 3,
	c = 4,
	b = 5,
	t = 6,
	e = 11,
	nu_e = 12,
	mu = 13,
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
}

namespace VINHAC {

//!  Class resposible for parton level process handling.
/*!
  This class generates leptons and photons. It calculates also kinematics of the process.
*/
class HardProcessHandler: public ProcessHandler {
private:
	VINHAC::Manager* _unnamed_Manager;

	DataSource *ds; //!< pointer to object with input parameters

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
	CLHEP::HepLorentzVector BOSTDQ(int MODE, CLHEP::HepLorentzVector QQ,
			CLHEP::HepLorentzVector PP);

	double generateMultiphotonRadiation(Event& event);
	double YWlRem(Event& event);
	void generatePhotonAngularDistribution(double aM, double aml, double amx2,
			double beta, double& costhg, double& sinthg, double& dist);
	double setLeptonMomentaInRadiativeDecay(Event& event, double sprim,
			double aM, CLHEP::HepLorentzVector phSum);
	double yfsCut(double am1, double am2, CLHEP::HepLorentzVector p1,
			CLHEP::HepLorentzVector p2, double epsg);
	CLHEP::HepLorentzVector Rot2le(double the, double phi,
			CLHEP::HepLorentzVector pvec);
	CLHEP::HepLorentzVector RXTOD2(double PHI, CLHEP::HepLorentzVector PVEC);
	CLHEP::HepLorentzVector RXTOD3(double PHI, CLHEP::HepLorentzVector PVEC);
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
	HardProcessHandler(DataSource *ds);

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
