/*
 * ElementFinder.cpp
 *
 *  Created on: 2009-05-31
 *      Author: kamil
 */

#include "ElementFinder.h"
#include <iostream>
#include <xercesc/util/XMLString.hpp>
#include <vector>

namespace VINHAC {

void ElementFinder::setUserDocument(DOMDocument* const doc) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: ElementFinder::setUserDocument()"<<std::endl;

#endif
	try {
		xmlUserDoc_ = doc;

		vector<string> nodesToSimplyReplace;
		nodesToSimplyReplace.push_back("pi");
		nodesToSimplyReplace.push_back("invGeV2toNb");
		nodesToSimplyReplace.push_back("sinThetaW2");
		nodesToSimplyReplace.push_back("invAlphaQED");
		nodesToSimplyReplace.push_back("FermiConst");
		nodesToSimplyReplace.push_back("alphaQCD");
		nodesToSimplyReplace.push_back("VCKM_11");
		nodesToSimplyReplace.push_back("VCKM_12");
		nodesToSimplyReplace.push_back("VCKM_13");
		nodesToSimplyReplace.push_back("VCKM_21");
		nodesToSimplyReplace.push_back("VCKM_22");
		nodesToSimplyReplace.push_back("VCKM_23");
		nodesToSimplyReplace.push_back("VCKM_31");
		nodesToSimplyReplace.push_back("VCKM_32");
		nodesToSimplyReplace.push_back("VCKM_33");
		nodesToSimplyReplace.push_back("maxWeightRej");
		nodesToSimplyReplace.push_back("printOut");
		nodesToSimplyReplace.push_back("intermidiateBosons");
		nodesToSimplyReplace.push_back("polarization");
		nodesToSimplyReplace.push_back("initialParticles");
		nodesToSimplyReplace.push_back("finialParticles");
		nodesToSimplyReplace.push_back("widthScheme");
		nodesToSimplyReplace.push_back("ips");
		nodesToSimplyReplace.push_back("factorizationScale");
		nodesToSimplyReplace.push_back("radiativeCorrections");
		nodesToSimplyReplace.push_back("Q2minCut");

		for (unsigned i = 0; i < nodesToSimplyReplace.size(); i++) {
			if (xmlUserDoc_->getElementsByTagName(DualString(
					nodesToSimplyReplace[i]).asXMLString())->getLength() == 1
					&& xmlDoc_->getElementsByTagName(DualString(
							nodesToSimplyReplace[i]).asXMLString())->getLength()
							== 1) {
				DOMNode
						*fromBaseFile =
								xmlDoc_->getElementsByTagName(DualString(
										nodesToSimplyReplace[i]).asXMLString())->item(
										0);
				DOMNode
						*fromUserFile =
								xmlUserDoc_->getElementsByTagName(DualString(
										nodesToSimplyReplace[i]).asXMLString())->item(
										0);

				fromUserFile = xmlDoc_->importNode(fromUserFile, true);
				fromBaseFile->getParentNode()->replaceChild(fromUserFile,
						fromBaseFile);

			}
		}

		vector<string> nodesToAttributesReplace;
		nodesToAttributesReplace.push_back("randomGenrator");
		nodesToAttributesReplace.push_back("beamA");
		nodesToAttributesReplace.push_back("beamB");

		for (unsigned i = 0; i < nodesToAttributesReplace.size(); i++) {
			if (xmlUserDoc_->getElementsByTagName(DualString(
					nodesToAttributesReplace[i]).asXMLString())->getLength()
					== 1 && xmlDoc_->getElementsByTagName(DualString(
					nodesToAttributesReplace[i]).asXMLString())->getLength()
					== 1) {
				DOMNode
						*fromBaseFile =
								xmlDoc_->getElementsByTagName(
										DualString(nodesToAttributesReplace[i]).asXMLString())->item(
										0);
				DOMNode
						*fromUserFile =
								xmlUserDoc_->getElementsByTagName(
										DualString(nodesToAttributesReplace[i]).asXMLString())->item(
										0);

				DOMNamedNodeMap *userAttr = fromUserFile->getAttributes();
				if (!userAttr)
					continue;
				DOMNamedNodeMap *baseAttr = fromBaseFile->getAttributes();

				for (unsigned i = 0; i < userAttr->getLength(); i++) {
					DOMNode *userAttrItem = userAttr->item(i);
					DOMNode *baseAttrItem = baseAttr->getNamedItem(
							userAttrItem->getNodeName());
					userAttrItem = xmlDoc_->importNode(userAttrItem, true);
					if (baseAttrItem) {
						baseAttr->removeNamedItem(baseAttrItem->getNodeName());
						baseAttr->setNamedItem(userAttrItem);
					} else {
						baseAttr->setNamedItem(userAttrItem);
					}

				}

			}
		}

		updateBeamSettings("beamA");
		updateBeamSettings("beamB");

	} catch (const XMLException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		std::cout << "Exception message is: \n" << message << "\n";
		XMLString::release(&message);

	} catch (const DOMException& toCatch) {
		char* message = XMLString::transcode(toCatch.msg);
		cout << "Exception message is: \n" << message << "\n";
		XMLString::release(&message);

	} catch (...) {
		cout << "Unexpected Exception \n";

	}
}

void ElementFinder::updateBeamSettings(string beamTagName) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: ElementFinder::updateBeamSettings()"<<std::endl;

#endif
	try {

		vector<string> nodesToSimplyReplace;
		nodesToSimplyReplace.push_back("idPDG");
		nodesToSimplyReplace.push_back("PDFsubset");
		nodesToSimplyReplace.push_back("xMin");
		nodesToSimplyReplace.push_back("xMax");
		nodesToSimplyReplace.push_back("atomicNb");
		nodesToSimplyReplace.push_back("chargeNb");
		nodesToSimplyReplace.push_back("energy");
		nodesToSimplyReplace.push_back("nuclearCorre");

		for (unsigned i = 0; i < nodesToSimplyReplace.size(); i++) {
			if (xmlUserDoc_->getElementsByTagName(DualString(
					nodesToSimplyReplace[i]).asXMLString())->getLength() >= 1
					&& xmlDoc_->getElementsByTagName(DualString(
							nodesToSimplyReplace[i]).asXMLString())->getLength()
							>= 1) {
				DOMNodeList *fromBaseFileList = xmlDoc_->getElementsByTagName(
						DualString(nodesToSimplyReplace[i]).asXMLString());
				DOMNodeList *fromUserFileList =
						xmlUserDoc_->getElementsByTagName(DualString(
								nodesToSimplyReplace[i]).asXMLString());

				DOMNode *fromUserFile = 0;
				DOMNode *fromBaseFile = 0;

				DOMNode *tmpNode = 0;

				for (unsigned i = 0; i < fromBaseFileList->getLength(); i++) {
					tmpNode = fromBaseFileList->item(i);

					DOMNode *tmp = tmpNode;

					while (tmp) {
						if (string(DualString(tmp->getNodeName()).asCString())
								== beamTagName)
							break;
						tmp = tmp->getParentNode();
					}

					if (tmp != 0) {
						fromBaseFile = tmpNode;
						break;
					}

				}

				for (unsigned i = 0; i < fromUserFileList->getLength(); i++) {
					tmpNode = fromUserFileList->item(i);

					DOMNode *tmp = tmpNode;

					while (tmp) {
						if (string(DualString(tmp->getNodeName()).asCString())
								== beamTagName)
							break;
						tmp = tmp->getParentNode();
					}

					if (tmp != 0) {
						fromUserFile = tmpNode;
						break;
					}

				}

				if (fromUserFile == 0 || fromBaseFile == 0)
					continue;

				fromUserFile = xmlDoc_->importNode(fromUserFile, true);
				fromBaseFile->getParentNode()->replaceChild(fromUserFile,
						fromBaseFile);

			}
		}

		vector<string> nodesToAttributesReplace;
		nodesToAttributesReplace.push_back("heavyIon");
		nodesToAttributesReplace.push_back("PDF");

		for (unsigned i = 0; i < nodesToAttributesReplace.size(); i++) {
			if (xmlUserDoc_->getElementsByTagName(DualString(
					nodesToAttributesReplace[i]).asXMLString())->getLength()
					>= 1 && xmlDoc_->getElementsByTagName(DualString(
					nodesToAttributesReplace[i]).asXMLString())->getLength()
					>= 1) {
				DOMNodeList *fromBaseFileList = xmlDoc_->getElementsByTagName(
						DualString(nodesToAttributesReplace[i]).asXMLString());
				DOMNodeList *fromUserFileList =
						xmlUserDoc_->getElementsByTagName(DualString(
								nodesToAttributesReplace[i]).asXMLString());

				DOMNode *fromUserFile = 0;
				DOMNode *fromBaseFile = 0;

				DOMNode *tmpNode = 0;

				for (unsigned i = 0; i < fromBaseFileList->getLength(); i++) {
					tmpNode = fromBaseFileList->item(i);

					DOMNode *tmp = tmpNode;

					while (tmp) {
						if (string(DualString(tmp->getNodeName()).asCString())
								== beamTagName)
							break;
						tmp = tmp->getParentNode();
					}

					if (tmp != 0) {
						fromBaseFile = tmpNode;
						break;
					}
				}

				for (unsigned i = 0; i < fromUserFileList->getLength(); i++) {
					tmpNode = fromUserFileList->item(i);

					DOMNode *tmp = tmpNode;

					while (tmp) {
						if (string(DualString(tmp->getNodeName()).asCString())
								== beamTagName)
							break;
						tmp = tmp->getParentNode();
					}

					if (tmp != 0) {
						fromUserFile = tmpNode;
						break;
					}

				}

				if (fromUserFile == 0 || fromBaseFile == 0)
					continue;

				DOMNamedNodeMap *userAttr = fromUserFile->getAttributes();
				if (!userAttr)
					continue;
				DOMNamedNodeMap *baseAttr = fromBaseFile->getAttributes();

				for (unsigned i = 0; i < userAttr->getLength(); i++) {
					DOMNode *userAttrItem = userAttr->item(i);
					DOMNode *baseAttrItem = baseAttr->getNamedItem(
							userAttrItem->getNodeName());
					userAttrItem = xmlDoc_->importNode(userAttrItem, true);
					if (baseAttrItem) {
						baseAttr->removeNamedItem(baseAttrItem->getNodeName());
						baseAttr->setNamedItem(userAttrItem);
					} else {
						baseAttr->setNamedItem(userAttrItem);
					}

				}

			}
		}

	} catch (const XMLException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		std::cout << "Exception message is: \n" << message << "\n";
		XMLString::release(&message);

	} catch (const DOMException& toCatch) {
		char* message = XMLString::transcode(toCatch.msg);
		cout << "Exception message is: \n" << message << "\n";
		XMLString::release(&message);

	} catch (...) {
		cout << "Unexpected Exception \n";

	}
}

void ElementFinder::inputEcho(string path) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: ElementFinder::inputEcho()"<<std::endl;

#endif

	XMLCh tempStr[100];
	XMLString::transcode("LS", tempStr, 99);
	DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(
			tempStr);
	DOMLSSerializer* theSerializer =
			((DOMImplementationLS*) impl)->createLSSerializer();

	try {
		// do the serialization through DOMLSSerializer::write();
		theSerializer->writeToURI(xmlDoc_,
				DualString(path + "/InputEcho.xml").asXMLString());
	} catch (const XMLException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		std::cout << "Exception message is: \n" << message << "\n";
		XMLString::release(&message);

	} catch (const DOMException& toCatch) {
		char* message = XMLString::transcode(toCatch.msg);
		cout << "Exception message is: \n" << message << "\n";
		XMLString::release(&message);

	} catch (...) {
		cout << "Unexpected Exception \n";

	}

	theSerializer->release();

}

}

