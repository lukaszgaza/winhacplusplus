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
	//Constants
	const DualString TAG_constant; //!< name for tag constant
	const DualString ATTR_constant_name; //!< name for attribute name
	const DualString TAG_CKM; //!< name for tag CKM
	const DualString ATTR_CKM_i; //!< name for attribute i
	const DualString ATTR_CKM_j; //!< name for attribute j
	//Beams
	const DualString TAG_beam;//!< name for tag beam
	const DualString ATTR_beam_id;//!< name for attribute id
	const DualString TAG_PDGid;//!< name for tag PDGid
	const DualString TAG_energy;//!< name for tag energy
	const DualString TAG_PDFs;//!< name for tag PDFs
	const DualString ATTR_PDFs_on;//!< name for attribute on
	const DualString TAG_interface;//!< name for tag interface
	const DualString TAG_name;//!< name for tag name
	const DualString TAG_subset;//!< name for tag subset
	const DualString TAG_xMin;//!< name for tag xMin
	const DualString TAG_xMax;//!< name for tag xMax
	const DualString TAG_Q2min;//!< name for tag Q2min


	//HardProcess
	const DualString TAG_quarks;//!< name for tag quarks
	const DualString TAG_bosons;//!< name for tag quarks
	const DualString TAG_leptons;//!< name for tag quarks
	const DualString TAG_polarizations;//!< name for tag polarizations
	const DualString TAG_polarization;//!< name for tag polarization
	const DualString TAG_radiation;//!< name for tag radiation
	const DualString TAG_level;//!< name for tag level
	const DualString TAG_softPhotonCutOff;//!< name for tag softPhotonCutOff
	const DualString TAG_modelReferenceFrame;//!< name for tag modelReferenceFrame
	const DualString TAG_factorizationScale;//!< name for tag factorizationScale
	const DualString TAG_alphaScheme;//!< name for tag alphaScheme
	const DualString TAG_widthScheme;//!< name for tag widthScheme

	//External
	const DualString TAG_pythiaXMLdocDir;//!< name for tag pythiaXMLdocDir
	const DualString TAG_partonShower;//!< name for tag partonShower
	const DualString TAG_hadronization;//!< name for tag hadronization

	//weights
	const DualString TAG_weightedEvents;//!< name for tag weightedEvents
	const DualString TAG_maxWeightRejection;//!< name for tag maxWeightRejection
	const DualString TAG_main;//!< name for tag main
	const DualString TAG_crudeWeight;//!< name for tag crudeWeight
	const DualString TAG_modelWeight;//!< name for tag modelWeight
	const DualString TAG_crude;//!< name for tag crude
	const DualString TAG_model;//!< name for tag model

	//others
	const DualString TAG_technical;//!< name for tag technical
	const DualString TAG_printOutLevel;//!< name for tag printOutLevel

	//particle data base
	const DualString TAG_particle;//!< name for tag particle

	ElementFinder() :
		TAG_constant("constant"),
		ATTR_constant_name("name"),
		TAG_CKM("CKM"),
		ATTR_CKM_i("i"),
		ATTR_CKM_j("j"),

		TAG_beam("beam"),
		ATTR_beam_id("id"),
		TAG_PDGid("PDGid"),
		TAG_energy("energy"),
		TAG_PDFs("PDFs"),
		ATTR_PDFs_on("on"),
		TAG_interface("interface"),
		TAG_name("name"),
		TAG_subset("subset"),
		TAG_xMin("xMin"),
		TAG_xMax("xMax"),
		TAG_Q2min("Q2min"),

		TAG_quarks("quarks"),
		TAG_bosons("bosons"),
		TAG_leptons("leptons"),
		TAG_polarizations("polarizations"),
		TAG_polarization("polarization"),
		TAG_radiation("radiation"),
		TAG_level("level"),
		TAG_softPhotonCutOff("softPhotonCutOff"),
		TAG_modelReferenceFrame("modelReferenceFrame"),
		TAG_factorizationScale("factorizationScale"),
		TAG_alphaScheme("alphaScheme"),
		TAG_widthScheme("widthScheme"),

		TAG_pythiaXMLdocDir("pythiaXMLdocDir"),
		TAG_partonShower("partonShower"),
		TAG_hadronization("hadronization"),

		TAG_weightedEvents("weightedEvents"),
		TAG_maxWeightRejection("maxWeightRejection"),
		TAG_main("main"),
		TAG_crudeWeight("crudeWeight"),
		TAG_modelWeight("modelWeight"),
		TAG_crude("crude"),
		TAG_model("model"),

		TAG_technical("technical"),
		TAG_printOutLevel("printOutLevel"),

		TAG_particle("particle")

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
	 * \brief returns elements specified by name
	 *
	 * @param name name of wanted tag
	 * @return wanted tag
	 */
	DOMNodeList* getElements(const XMLCh* name) {
		DOMNodeList* list = xmlDoc_->getElementsByTagName(name);
		return list;
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
	 * \brief returns subelements of some element specified by name
	 *
	 * @param name name of wanted tag
	 * @param elements element to search in
	 * @return wanted tag
	 */
	DOMNodeList* getSubElements(const XMLCh* name, DOMElement* element) {
		DOMNodeList* list = element ->getElementsByTagName(name);
		return list;
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
	DOMDocument* xmlDoc_;
	DOMDocument* xmlUserDoc_;

}; // class ElementFinder
}
#endif /* ELEMENTFINDER_H_ */
