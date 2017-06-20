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
#include "../core/VinhacException.h"

namespace VINHAC {

void ElementFinder::setUserDocument(DOMDocument* const doc) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: ElementFinder::setUserDocument()"<<std::endl;

#endif
	try {
		xmlUserDoc_ = doc;

		vector < string > nodesToSimplyReplace;
		nodesToSimplyReplace.push_back("beams");
		nodesToSimplyReplace.push_back("quarks");
		nodesToSimplyReplace.push_back("bosons");
		nodesToSimplyReplace.push_back("leptons");
		nodesToSimplyReplace.push_back("polarizations");
		nodesToSimplyReplace.push_back("radiation");
		nodesToSimplyReplace.push_back("modelReferenceFrame");
		nodesToSimplyReplace.push_back("factorizationScale");
		nodesToSimplyReplace.push_back("alphaScheme");
		nodesToSimplyReplace.push_back("widthScheme");
		nodesToSimplyReplace.push_back("pythiaXMLdocDir");
		nodesToSimplyReplace.push_back("partonShower");
		nodesToSimplyReplace.push_back("hadronization");
		nodesToSimplyReplace.push_back("weightedEvents");
		nodesToSimplyReplace.push_back("main");
		nodesToSimplyReplace.push_back("crude");
		nodesToSimplyReplace.push_back("model");
		nodesToSimplyReplace.push_back("printOutLevel");

		for (unsigned i = 0; i < nodesToSimplyReplace.size(); i++) {
			if (xmlUserDoc_->getElementsByTagName(
					DualString(nodesToSimplyReplace[i]).asXMLString())->getLength()
					== 1
					&& xmlDoc_->getElementsByTagName(
							DualString(nodesToSimplyReplace[i]).asXMLString())->getLength()
							== 1) {
				DOMNode
						*fromBaseFile =
								xmlDoc_->getElementsByTagName(
										DualString(nodesToSimplyReplace[i]).asXMLString())->item(
										0);
				DOMNode
						*fromUserFile =
								xmlUserDoc_->getElementsByTagName(
										DualString(nodesToSimplyReplace[i]).asXMLString())->item(
										0);

				fromUserFile = xmlDoc_->importNode(fromUserFile, true);
				fromBaseFile->getParentNode()->replaceChild(fromUserFile,
						fromBaseFile);
			}
		}

		StringManager sm;
		//replace constants
		if (xmlUserDoc_->getElementsByTagName(TAG_constant.asXMLString())->getLength()
				> 0) {

			DOMNodeList* constantsFromBaseFile = xmlDoc_->getElementsByTagName(
					TAG_constant.asXMLString());

			DOMNodeList* constantsFromUserFile =
					xmlUserDoc_->getElementsByTagName(
							TAG_constant.asXMLString());

			DOMElement *fromBaseFile = 0;
			DOMElement *fromUserFile = 0;
			for (unsigned i = 0; i < constantsFromUserFile->getLength(); ++i) {
				fromUserFile = 0;
				if (DOMNode::ELEMENT_NODE
						== constantsFromUserFile->item(i)->getNodeType()) {
					fromUserFile
							= dynamic_cast<DOMElement*> (constantsFromUserFile->item(
									i));
				}

				if (fromUserFile != 0) {
					for (unsigned j = 0; j < constantsFromBaseFile->getLength(); ++j) {
						fromBaseFile = 0;
						if (DOMNode::ELEMENT_NODE
								== constantsFromBaseFile->item(j)->getNodeType()) {
							fromBaseFile
									= dynamic_cast<DOMElement*> (constantsFromBaseFile->item(
											j));
						}
						if (fromBaseFile != 0) {
							if (string(sm.convert(
									fromUserFile->getAttribute(
											ATTR_constant_name.asXMLString())))
									== string(sm.convert(
											fromBaseFile->getAttribute(
													ATTR_constant_name.asXMLString())))) {
								break;
							} else {
								fromBaseFile = 0;
							}
						}
					}
				}

				if (fromUserFile != 0 && fromBaseFile != 0) {
					DOMNode* fromUserFileImp = xmlDoc_->importNode(
							fromUserFile, true);
					fromBaseFile->getParentNode()->replaceChild(
							fromUserFileImp, fromBaseFile);
				}
			}
		}

		//replace CKM
		if (xmlUserDoc_->getElementsByTagName(TAG_CKM.asXMLString())->getLength()
				> 0) {

			DOMNodeList* constantsFromBaseFile = xmlDoc_->getElementsByTagName(
					TAG_CKM.asXMLString());

			DOMNodeList* constantsFromUserFile =
					xmlUserDoc_->getElementsByTagName(TAG_CKM.asXMLString());

			DOMElement *fromBaseFile = 0;
			DOMElement *fromUserFile = 0;
			for (unsigned i = 0; i < constantsFromUserFile->getLength(); ++i) {
				fromUserFile = 0;
				if (DOMNode::ELEMENT_NODE
						== constantsFromUserFile->item(i)->getNodeType()) {
					fromUserFile
							= dynamic_cast<DOMElement*> (constantsFromUserFile->item(
									i));
				}

				if (fromUserFile != 0) {
					for (unsigned j = 0; j < constantsFromBaseFile->getLength(); ++j) {
						fromBaseFile = 0;
						if (DOMNode::ELEMENT_NODE
								== constantsFromBaseFile->item(j)->getNodeType()) {
							fromBaseFile
									= dynamic_cast<DOMElement*> (constantsFromBaseFile->item(
											j));
						}
						if (fromBaseFile != 0) {
							if (string(sm.convert(
									fromUserFile->getAttribute(
											ATTR_CKM_i.asXMLString())))
									== string(sm.convert(
											fromBaseFile->getAttribute(
													ATTR_CKM_i.asXMLString())))
									&& string(sm.convert(
											fromUserFile->getAttribute(
													ATTR_CKM_j.asXMLString())))
											== string(sm.convert(
													fromBaseFile->getAttribute(
															ATTR_CKM_j.asXMLString())))) {
								break;
							} else {
								fromBaseFile = 0;
							}
						}
					}
				}

				if (fromUserFile != 0 && fromBaseFile != 0) {
					DOMNode* fromUserFileImp = xmlDoc_->importNode(
							fromUserFile, true);
					fromBaseFile->getParentNode()->replaceChild(
							fromUserFileImp, fromBaseFile);
				}
			}
		}

	} catch (const XMLException& toCatch) {
		char* message = XMLString::transcode(toCatch.getMessage());
		std::cout << "Exception message is: \n" << message << "\n";
		string msg(message);
		XMLString::release(&message);
		throw VinhacException(msg);

	} catch (const DOMException& toCatch) {
		char* message = XMLString::transcode(toCatch.msg);
		cout << "Exception message is: \n" << message << "\n";
		string msg(message);
		XMLString::release(&message);
		throw VinhacException(msg);

	} catch (...) {
		throw VinhacException("Unexpected exception during UserFile reading");

	}
}

/*void ElementFinder::updateBeamSettings(string beamTagName) {
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
 }*/

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
		string msg(message);
		XMLString::release(&message);
		throw VinhacException(msg);

	} catch (const DOMException& toCatch) {
		char* message = XMLString::transcode(toCatch.msg);
		cout << "Exception message is: \n" << message << "\n";
		string msg(message);
		XMLString::release(&message);
		throw VinhacException(msg);

	} catch (...) {
		throw VinhacException("Unexpected error during echoing input");
	}

	theSerializer->release();

}

}

