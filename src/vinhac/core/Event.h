#include <string>
#include <vector>
#include <map>
#include <exception>
#include <sstream>
#include "Particle.h"
#include "enum/WeightsDefinition.h"
#include "enum/PolarizationDefinition.h"

#ifndef _EVENT_H_
#define _EVENT_H_

using namespace std;

namespace VINHAC {
class Particle;
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

	bool pdfUsed;
	double quarkX;
	double antiquarkX;
	double quarkXpdf;
	double antiquarkXpdf;
	double pdfScale;
	int quarkPDFid;
	int antiquarkPDFid;

	map<WeightName::name, double> m_commonWeightsMap; //!< map of weights, key name of weight
	map<pair<WeightName::name,PolarizationName::name>, double> m_modelWeightsMap; //!< map of weights, key name of weight
	double m_crudeXsection;
	double differentialXsectionCrude;
	double m_cosTheta;
	double m_phi;

	void printParticle(stringstream& ss, Particle& p);
	void printParticle(stringstream& ss, Particle& p, FrameName::name name);
	void printParticleWrf(stringstream& ss, Particle& p);

	friend ostream& operator<<(ostream& out, const Event& event);
	friend class BeamHandler;
public:
	Event();

	/**
	 * \brief returns value of common weight specified by name
	 *
	 * @param wtName name of weight
	 * @return value of weight
	 */
	inline double getCommonWeightByName(const WeightName::name& wtName) const {
		return m_commonWeightsMap.find(wtName)->second;
	}

	/**
	 * \brief adds value of common weight specified by name
	 *
	 * @param wtName name of weight
	 * @param wt value of weight
	 */
	inline void addCommonWeight(const WeightName::name& wtName,const double& wt) {
		m_commonWeightsMap.insert(pair<WeightName::name, double> (wtName, isnan(wt)?0.0:wt));
	}

	/**
	 * \brief returns common weights map
	 *
	 * @return common weights map
	 */
	inline const map<WeightName::name, double>& getCommonWeightsMap() const{
		return m_commonWeightsMap;
	}

	/**
	 * \brief returns value of model weight specified by name
	 *
	 * @param wtName name of weight
	 * @return value of weight
	 */
	inline double getModelWeightByName(const pair<WeightName::name,PolarizationName::name>& wtName) const{
		return m_modelWeightsMap.find(wtName)->second;
	}

	/**
	 * \brief returns value of model weight specified by name
	 *
	 * @param wtName name of weight
	 * @return value of weight
	 */
	inline double getModelWeightByName(const WeightName::name& wtName) const{
		return m_modelWeightsMap.find(make_pair(wtName,PolarizationName::unpolarized))->second;
	}

	/**
	 * \brief returns value of model weight specified by name and polarization
	 *
	 * @param wtName name of weight
	 * @return value of weight
	 */
	inline double getModelWeightByName(const WeightName::name& wtName, const PolarizationName::name& polName) const{
		return m_modelWeightsMap.find(make_pair(wtName,polName))->second;
	}

	/**
	 * \brief adds value of common weight specified by name
	 *
	 * @param wtName name of weight
	 * @param wt value of weight
	 */
	inline void addModelWeight(const pair<WeightName::name,PolarizationName::name>& wtName,const double& wt) {
		m_modelWeightsMap.insert(make_pair(wtName, isnan(wt)?0.0:wt));
	}

	/**
	 * \brief adds value of common weight specified by name
	 *
	 * @param wtName name of weight
	 * @param wt value of weight
	 */
	inline void addModelWeight(const WeightName::name& wtName,const double& wt) {
		m_modelWeightsMap.insert(make_pair(make_pair(wtName,PolarizationName::unpolarized), isnan(wt)?0.0:wt));
	}

	/**
	 * \brief returns model weights map
	 *
	 * @return model weights map
	 */
	inline const map<pair<WeightName::name,PolarizationName::name>, double>& getModelWeightsMap() const {
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
	 * \brief return quark's Bjorken's x
	 *
	 * @return quark's Bjorken's x
	 */
	inline double getQuarkX() const {
		return quarkX;
	}

	/**
	 * \brief return antiquark's Bjorken's x
	 *
	 * @return antiquark's Bjorken's x
	 */
	inline double getAntiquarkX() const {
		return antiquarkX;
	}

	/**
	 * \brief return quark's  x times pdf
	 *
	 * @return quark's x times pdf
	 */
	inline double getQuarkXpdf() const {
		return quarkXpdf;
	}

	/**
	 * \brief return antiquark's  x times pdf
	 *
	 * @return antiquark's x times pdf
	 */
	inline double getAntiquarkXpdf() const {
		return antiquarkXpdf;
	}

	/**
	 * \brief return pdf scale
	 *
	 * @return  pdf scale
	 */
	inline double getPDFScale() const {
		return pdfScale;
	}

	/**
	 * \brief return if pdf was used
	 *
	 * @return  if pdf was used
	 */
	inline bool isPDFUsed() const {
		return pdfUsed;
	}

	/**
	 * \brief return id of pdf used for quark
	 *
	 * @return  id of pdf used for quark
	 */
	inline int getQuarkPDFid() const{
		return quarkPDFid;
	}

	/**
	 * \brief return id of pdf used for antiquark
	 *
	 * @return  id of pdf used for antiquark
	 */
	inline int getAntiquarkPDFid() const{
		return antiquarkPDFid;
	}

	/**
	 * \brief computes main weight
	 *
	 * 	Main weight is multiplication of common weights and best model weight
	 *
	 * @return main weight
	 */
	double mainWt() const;

	/**
	 *	\brief overrides event's weights and sets custom weight
	 */
	void setMainWt(double mainWt);

	/**
	 * \brief computes whole common
	 *
	 * @return multiplication of all common weights
	 */
	inline double currentCommonWtMultiplication() const {
		map<WeightName::name, double>::const_iterator it;
		double wtMultiplication = 1.0;
		for (it = m_commonWeightsMap.begin(); it != m_commonWeightsMap.end(); it++) {
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
	inline void setXsection(const double& crude) {
		m_crudeXsection = crude;
	}

	/**
	 * \brief sets differential crude cross section
	 *
	 * @param dsigcru differential crude cross section to set
	 */
	inline void setDifferentialXsectionCrude(const double& dsigcru) {
		differentialXsectionCrude = dsigcru;
	}

	/**
	 * \brief returns differential crude cross section
	 *
	 * @return differential crude cross section
	 */
	inline double getDifferentialXsectionCrude() const {
		return differentialXsectionCrude;
	}

	/**
	 * \brief returns \f$cos(\theta)\f$, \f$ \theta\f$ - polar angle of lepton
	 *
	 * @return \f$cos(\theta)\f$
	 */
	inline double getCosTheta() const {
		return m_cosTheta;
	}

	/**
	 * \brief sets \f$cos(\theta)\f$, \f$ \theta\f$ - polar angle of lepton
	 *
	 * @param cos \f$cos(\theta)\f$ to set
	 */
	inline void setCosTheta(const double& cos) {
		m_cosTheta = cos;
	}

	/**
	 * \brief sets \f$\phi\f$ - azimuthal angle of lepton
	 *
	 * @param phi \f$\phi\f$ to set
	 */
	inline void setPhi(const double& phi) {
		m_phi = phi;
	}

	/**
	 * \brief returns \f$\phi\f$ - azimuthal angle of lepton
	 *
	 * @return phi \f$\phi\f$
	 */
	inline double getPhi() const {
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

	/**
	 * \brief generates string containing event parameters in some frame
	 *
	 * @return string which can be printed
	 */
	string printInFrame(FrameName::name name);

};

}

#endif
