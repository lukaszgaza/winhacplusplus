#include <string>
#include <vector>
#include <map>
#include <exception>
#include <sstream>

#ifndef _EVENT_H_
#define _EVENT_H_

#include "Particle.h"

using namespace std;

namespace VINHAC {
class Particle;
class Event;
class EventStep;
}

namespace VINHAC {
/**
 * \brief class representing one generated event.
 *
 * This class stores all particles involved in one event.
 * There is also stored differential crude cross section, lepton angles
 * and maps of common and model weights
 */
class Event {
private:
	Particle beamParticleA;
	Particle beamParticleB;
	Particle beamRestA;
	Particle beamRestB;
	Particle quark;
	Particle antiQuark;
	Particle boson;
	Particle neutrino;
	Particle lepton;
	vector<Particle> photons;

	friend ostream& operator<<(ostream& out, const Event& event);
	map<string, double> m_commonWeightsMap; //!< map of weights, key name of weight
	map<string, double> m_modelWeightsMap; //!< map of weights, key name of weight
	double m_crudeXsection;
	double differentialXsectionCrude;
	double m_cosTheta;
	double m_phi;

	void printParticle(stringstream& ss, Particle& p);
	void printParticleWrf(stringstream& ss, Particle& p);
public:
	Event();

	/**
	 * \brief returns value of common weight specified by name
	 *
	 * @param wtName name of weight
	 * @return value of weight
	 */
	inline double getCommonWeightByName(string wtName) {
		return m_commonWeightsMap.find(wtName)->second;
	}

	/**
	 * \brief adds value of common weight specified by name
	 *
	 * @param wtName name of weight
	 * @param wt value of weight
	 */
	inline void addCommonWeight(string wtName, double wt) {
		m_commonWeightsMap.insert(pair<string, double> (wtName, wt));
	}

	/**
	 * \brief returns common weights map
	 *
	 * @return common weights map
	 */
	inline map<string, double> getCommonWeightsMap() {
		return m_commonWeightsMap;
	}

	/**
	 * \brief returns value of model weight specified by name
	 *
	 * @param wtName name of weight
	 * @return value of weight
	 */
	inline double getModelWeightByName(string wtName) {
		return m_modelWeightsMap.find(wtName)->second;
	}

	/**
	 * \brief adds value of common weight specified by name
	 *
	 * @param wtName name of weight
	 * @param wt value of weight
	 */
	inline void addModelWeight(string wtName, double wt) {
		m_modelWeightsMap.insert(pair<string, double> (wtName, wt));
	}

	/**
	 * \brief returns model weights map
	 *
	 * @return model weights map
	 */
	inline map<string, double> getModelWeightsMap() {
		return m_modelWeightsMap;
	}

	/**
	 * \brief return beam particle A
	 *
	 * @return beam particle A
	 */
	inline Particle& getBeamParticleA() {
		return beamParticleA;
	}

	/**
	 * \brief return beam particle B
	 *
	 * @return beam particle B
	 */
	inline Particle& getBeamParticleB() {
		return beamParticleB;
	}

	/**
	 * \brief return beam rest A
	 *
	 * @return beam rest A
	 */
	inline Particle& getBeamRestA() {
		return beamRestA;
	}

	/**
	 * \brief return beam rest B
	 *
	 * @return beam rest B
	 */
	inline Particle& getBeamRestB() {
		return beamRestB;
	}

	/**
	 * \brief return quark
	 *
	 * @return quark
	 */
	inline Particle& getQuark() {
		return quark;
	}

	/**
	 * \brief return antiquark
	 *
	 * @return antiquark
	 */
	inline Particle& getAntiQuark() {
		return antiQuark;
	}

	/**
	 * \brief return boson
	 *
	 * @return boson
	 */
	inline Particle& getBoson() {
		return boson;
	}

	/**
	 * \brief return lepton
	 *
	 * @return lepton
	 */
	inline Particle& getLepton() {
		return lepton;
	}

	/**
	 * \brief return neutrino
	 *
	 * @return neutrino
	 */
	inline Particle& getNeutrino() {
		return neutrino;
	}

	/**
	 * \brief return vector of photons
	 *
	 * @return vector of photons
	 */
	inline vector<Particle>& getPhotons() {
		return photons;
	}


	/**
	 * \brief computes main weight
	 *
	 * 	Main weight is multiplication of common weights and best model weight
	 *
	 * @return main weight
	 */
	inline double mainWt() {
		double result = 1.0;
		for (map<string, double>::iterator it = m_commonWeightsMap.begin(); it
				!= m_commonWeightsMap.end(); ++it) {
			result *= (*it).second;
		}
		result *= m_modelWeightsMap["ModelBest"];
		return result;
	}

	/**
	 * \brief computes whole common
	 *
	 * @return multiplication of all common weights
	 */
	inline double currentCommonWtMultiplication() {
		map<string, double>::iterator it;
		double wtMultiplication = 1.0;
		for (it = m_commonWeightsMap.begin(); it != m_commonWeightsMap.end(); it++) {
			string name = (*it).first;
			double wt = (*it).second;
			wtMultiplication = wtMultiplication * wt;
		}
		return wtMultiplication;
	}

	/**
	 * \brief cleans event data
	 */
	bool clean();

	/**
	 * \brief sets crude cross section
	 *
	 * @param crude crude cross section to set
	 */
	inline void setXsection(double crude) {
		m_crudeXsection = crude;
	}

	/**
	 * \brief sets differential crude cross section
	 *
	 * @param dsigcru differential crude cross section to set
	 */
	inline void setDifferentialXsectionCrude(double dsigcru) {
		differentialXsectionCrude = dsigcru;
	}

	/**
	 * \brief returns differential crude cross section
	 *
	 * @return differential crude cross section
	 */
	inline double getDifferentialXsectionCrude() {
		return differentialXsectionCrude;
	}

	/**
	 * \brief returns \f$cos(\theta)\f$, \f$ \theta\f$ - polar angle of lepton
	 *
	 * @return \f$cos(\theta)\f$
	 */
	inline double getCosTheta() {
		return m_cosTheta;
	}

	/**
	 * \brief sets \f$cos(\theta)\f$, \f$ \theta\f$ - polar angle of lepton
	 *
	 * @param cos \f$cos(\theta)\f$ to set
	 */
	inline void setCosTheta(double cos) {
		m_cosTheta = cos;
	}

	/**
	 * \brief sets \f$\phi\f$ - azimuthal angle of lepton
	 *
	 * @param phi \f$\phi\f$ to set
	 */
	inline void setPhi(double phi) {
		m_phi = phi;
	}

	/**
	 * \brief returns \f$\phi\f$ - azimuthal angle of lepton
	 *
	 * @return phi \f$\phi\f$
	 */
	inline double getPhi() {
		return m_phi;
	}

	/**
	 * \brief generates string containing event parameters in laboratory frame
	 *
	 * @return string which can be printed
	 */
	string printInLabFrame();

	/**
	 * \brief generates string containing event parameters in quarks cms frame
	 *
	 * @return string which can be printed
	 */
	string printInQuarksCMSFrame();

};

}

#endif
