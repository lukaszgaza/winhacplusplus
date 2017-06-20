/*
 * ElementFinder.h
 *
 *  Created on: 2009-05-15
 *      Author: siudmy
 */

#ifndef ELEMENTFINDER_H_
#define ELEMENTFINDER_H_

#include<xercesc/util/PlatformUtils.hpp>

#include<xercesc/dom/DOM.hpp>
#include<xercesc/dom/DOMDocument.hpp>
#include<xercesc/dom/DOMDocumentType.hpp>
#include<xercesc/dom/DOMElement.hpp>
#include<xercesc/dom/DOMImplementation.hpp>
#include<xercesc/dom/DOMImplementationLS.hpp>
#include<xercesc/dom/DOMNodeIterator.hpp>
#include<xercesc/dom/DOMNodeList.hpp>
#include<xercesc/dom/DOMText.hpp>
#include<xercesc/util/XMLString.hpp>

#include<xercesc/parsers/XercesDOMParser.hpp>
#include<xercesc/util/XMLUni.hpp>
#include<string>
#include "helper-classes.h"
using namespace xercesc;
using namespace std;

namespace VINHAC {
//!  helper class for Xerces XML parser
/*!
 * It helps to find tags, attributes and values in input file.
*/
class ElementFinder {

public:
	//XML input elements tags
	const DualString TAG_idPDG_; //!< name for tag idPDG
	//PDF
	const DualString TAG_Q2minCut_;//!< name for tag Q2minCut
	const DualString TAG_keyPDF_;//!< name for tag keyPDF
	const DualString ATTR_name_;//!< name for attribute name
	const DualString TAG_xMinPDF_;//!< name for tag xMinPDF
	const DualString TAG_xMaxPDF_;//!< name for tag xMaxPDF
	const DualString TAG_Q2minPDF_;//!< name for tag Q2minPDF
	const DualString ATTR_PDFinterface_;//!< name for attribute PDFinterface
	const DualString TAG_PDFsubset_;//!< name for tag PDFsubset
	//HeavyIon
	const DualString TAG_heavyIon_;//!< name for tag heavyIon
	const DualString ATTR_switchOn_;//!< name for attribute switchOn
	const DualString TAG_atomicNb_;//!< name for tag atomicNb
	const DualString TAG_chargeNb_;//!< name for tag chargeNb
	const DualString TAG_nuclearCorre_;//!< name for tag nuclearCorre
	const DualString ATTR_ENERGY_;//!< name for attribute ENERGY

	//Initial particle
	const DualString TAG_initialParticles_;//!< name for tag initialParticles
	const DualString TAG_finialParticles_;//!< name for tag finialParticles
	const DualString TAG_intermidiateBosons_;//!< name for tag intermidiateBosons
	const DualString TAG_particleName_;//!< name for tag particleName
	const DualString TAG_particle_;//!< name for tag particle

	//const
	const DualString TAG_constMath_;//!< name for tag constMath
	const DualString TAG_constPhys_;//!< name for tag constPhys

	//Schemes
	const DualString TAG_IPS_;//!< name for tag IPS
	const DualString TAG_widthSchem_;//!< name for tag widthSchem
	const DualString TAG_factorizationScale_;//!< name for tag factorizationScale

	const DualString TAG_POLARIZATION;//!< name for tag POLARIZATION
	const DualString TAG_CKM;//!< name for tag CKM

	const DualString TAG_ewCorrections_;//!< name for tag ewCorrections
	const DualString TAG_correctionType_;//!< name for tag correctionType

	ElementFinder() :
				//XML input elements tags
				TAG_idPDG_("idPDG"),
				TAG_Q2minCut_("Q2minCut"),
				TAG_keyPDF_("PDF"),
				ATTR_name_("name"),
				TAG_xMinPDF_("xMin"),
				TAG_xMaxPDF_("xMax"),
				TAG_Q2minPDF_("Q2min"),
				ATTR_PDFinterface_("PDFinterface"),
				TAG_PDFsubset_("PDFsubset"),
				//HeavyIon
				TAG_heavyIon_("heavyIon"),
				ATTR_switchOn_("switchOn"),
				TAG_atomicNb_("atomicNb"),
				TAG_chargeNb_("chargeNb"),
				TAG_nuclearCorre_("nuclearCorre"),

				//XML input Attributes
				ATTR_ENERGY_("energy"),

				//Initial particle
				TAG_initialParticles_("initialParticles"),
				TAG_finialParticles_("finialParticles"),
				TAG_intermidiateBosons_("intermidiateBosons"),

				TAG_particleName_("particleName"),
				TAG_particle_("particle"),

				//const
				TAG_constMath_("constMath"),
				TAG_constPhys_("constPhys"),

				//Schemes
				TAG_IPS_("ips"), TAG_widthSchem_("widthScheme"),
				TAG_factorizationScale_("factorizationScale"),

				TAG_POLARIZATION("polarization"), TAG_CKM("VCKM_"),
				TAG_ewCorrections_("ewCorrections"), TAG_correctionType_(
						"correctionType"),

				xmlDoc_(0), xmlUserDoc_(0)

	{
		return;
	}

	/**
	 * \brief returns first document element
	 *
	 * @return first document element
	 */
	DOMElement* getConfigElement() {
		DOMElement* result = xmlDoc_->getDocumentElement();
		return (result);
	}

	/**
	 * \brief returns element specified by name
	 *
	 * @param name name of wanted tag
	 * @return wanted tag
	 */
	DOMElement* getElement(const XMLCh* name) {
		DOMElement* result = NULL;
		DOMNodeList* list = xmlDoc_->getElementsByTagName(name);
		DOMNode* node = list->item(0);
		if (DOMNode::ELEMENT_NODE == node->getNodeType()) {
			result = dynamic_cast<DOMElement*> (node);
		}
		return (result);
	}

	/**
	 * \brief returns subelement of some element specified by name
	 *
	 * @param name name of wanted tag
	 * @param element element to search in
	 * @return wanted tag
	 */
	DOMElement* getSubElement(const XMLCh* name, DOMElement* element) {
		DOMElement* result = NULL;
		DOMNodeList* list = element ->getElementsByTagName(name);
		DOMNode* node = list->item(0);
		if (DOMNode::ELEMENT_NODE == node->getNodeType()) {
			result = dynamic_cast<DOMElement*> (node);
		}
		return (result);
	}

	/**
	 * \brief returns text content of subelement specified by name
	 *
	 * @param subElementName name of wanted subelement
	 * @param element element to search in
	 * @return text content of found subelement
	 */
	string getSubElementContent(const XMLCh* subElementName,
			DOMElement* element) {
		DOMElement* subElement = getSubElement(subElementName, element);
		const DualString subElementContent = subElement->getTextContent();
		return subElementContent.asCString();
	}

	/**
	 * \brief sets xml document
	 *
	 * @param doc xml document to set
	 */
	void setDocument(DOMDocument* const doc) {
		xmlDoc_ = doc;
	}

	/**
	 * \brief sets user's xml document
	 *
	 * @param doc user's xml document to set
	 */
	void setUserDocument(DOMDocument* const doc);

	/**
	 * \brief produces echo of input
	 *
	 * @param path path where input echo should be placed
	 */
	void inputEcho(string path);

private:
	void updateBeamSettings(string beamTagName);
	DOMDocument* xmlDoc_;
	DOMDocument* xmlUserDoc_;

}; // class ElementFinder
}
#endif /* ELEMENTFINDER_H_ */
