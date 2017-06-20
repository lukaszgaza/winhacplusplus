/*
 * DataSource.h
 *
 *  Created on: 2009-11-27
 *      Author: kamil
 */

#ifndef DATASOURCE_H_
#define DATASOURCE_H_

#include <map>
#include <vector>
#include <set>
#include <utility>
#include <string>

#include "../core/enum/PolarizationDefinition.h"
#include "../core/enum/WeightsDefinition.h"
#include "../core/enum/FrameDefinition.h"
#include "../core/enum/PDFInterfaceDefinition.h"
#include "../core/enum/PrintOutLevelDefinition.h"

namespace VINHAC {
//!  Class used for storage of input parameters. This is a singleton !.
/*!
 * It is filled by InitializationHandler and used by other components
 * @see InitializationHandler
*/
class DataSource {
private:

	DataSource();
	DataSource(const DataSource &);
	DataSource& operator=(const DataSource&);
public:
	static DataSource& get(){
		static DataSource instance;
		return instance;
	}

	double pi; //!< \f$\pi\f$ constant
	double alphaQED; //!<   QED coupling constant
	double alphaQCD; //!<   QCD coupling constant
	double invAlphaQED; //!<  inversion of QED coupling constant
	double invGeV2toNb; //!< inversion of \f$GeV^2\f$ to \f$nb\f$ conversion constant
	double geV2toNb; //!< \f$GeV^2\f$ to \f$nb\f$ conversion constant
	double fermiConst; //!< Fermi's constant

	double MW; //!< W boson mass
	double MZ; //!< Z boson mass
	double GW; //!< W boson width
	double GZ; //!< Z boson width
	double MH; //!< Higg's boson mass
	double GH; //!< Higg's boson width

	double sinThetaW2; //!< \f$ sin^2(\theta_W) \f$ (\f$ \theta_W \f$ - Weinberg's angle)

	std::map<int, double> masses; //!< map of particle masses, key - absolute value of PDGid
	std::vector<int> quarksCharges; //!< vector of quark's charges (indexed by PDGid)

	std::vector<int> initialParticles; //!< vector of PDGids of initial particles
	std::vector<int> intermediateParticles; //!< vector of PDGids of intermediate particles
	std::vector<int> finalParticles; //!< vector of PDGids of final particles


	bool isPDFOn; //!< if PDF is active
	PDFInterface::name pdfInterface; //!< PDF interface
	std::string pdfName;
	int pdfSubset;
	double pdfXMin;
	double pdfXMax;
	double pdfQ2min;
	std::map<std::string, std::string> beamA; //!< map of beam A settings
	std::map<std::string, std::string> beamB; //!< map of beam B settings

	std::vector<std::vector<double> > ckm; //!< CKM matrix
	std::vector<std::vector<double> > ckm2; //!< CKM matrix squared
	double q2MinCut; //!< minimum of quarks \f$Q^2 [GeV^2]\f$
	double EgMin; //!< Soft-photon cut-off in GeV
	double Ecm; //!< energy in center of mass of beams

	std::string pythiaXMLdocDir;//!< pythiaXMLdocDir
	bool isPartonShower; //!< QCD partonshower
	bool isHadronization; //!< Hadronization

	int keyPol; //!< W boson polarization switch
	int keyRad; //!< radiation type switch
	int keyWid; //!< width type switch
	int keyGmu; //!< input-parameter scheme (IPS) switch
	int keyFac; //!< factorization scale switch
	int keyCol; //!< collider type switch

	int keyWgt; //!< 0 - unweighted events, 1 - weighted events
	double maxWeightRej; //!<

	FrameName::name modelReferenceFrame;
	std::set<WeightName::name> weights;
	std::set<PolarizationName::name> polarizations;
	std::pair<WeightName::name,PolarizationName::name> mainWeight;

	PrintOutLevel::name printOutLevel;

};

}

#endif /* DATASOURCE_H_ */
