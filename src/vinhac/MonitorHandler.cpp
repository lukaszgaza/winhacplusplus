/*
 * \file MonitorHandler.cpp
 * \brief
 * \date Created on: 2009-07-07
 * \author     Author: siudmy
 */
#include "MonitorHandler.h"
#include <fstream>
#include "../utils/Functions.h"
#include "VinhacException.h"

VINHAC::MonitorHandler::MonitorHandler() :
	mainWeightMonitor("MainWeight") {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::MonitorHandler()"<<std::endl;

#endif

	setHandlerName("MonitorHandler");
	m_commonWtMonitors.clear();
	m_modelWtMonitors.clear();
	xsYFSMonitors.clear();
	xsFixedOrderMonitors.clear();
	xsAllMonitors.clear();
	mode = 0;
}

bool VINHAC::MonitorHandler::eventEvolve(VINHAC::Event& newEvent) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::eventEvolve()"<<std::endl;

#endif

	//	vector<double> weights = newEvent.getWeights();
	map<string, double> commonWeights = newEvent.getCommonWeightsMap();
	map<string, double> modelWeights = newEvent.getModelWeightsMap();

	if (mode == 0) {
		for (map<string, double>::iterator it = commonWeights.begin(); it
				!= commonWeights.end(); ++it) {
			VINHAC::WeightMonitor monitor((*it).first);
			m_commonWtMonitors[(*it).first] = monitor;
		}
		for (map<string, double>::iterator it = modelWeights.begin(); it
				!= modelWeights.end(); ++it) {
			VINHAC::WeightMonitor monitor((*it).first);
			m_modelWtMonitors[(*it).first] = monitor;

			VINHAC::WeightMonitor monitorXS((*it).first);
			xsAllMonitors[(*it).first] = monitorXS;
		}

		if (m_modelWtMonitors.count("beta00") > 0) {
			VINHAC::WeightMonitor monitor("XsTot     O(alf0)");
			xsYFSMonitors["XsTot     O(alf0)"] = monitor;
		}

		if (m_modelWtMonitors.count("beta01_LL + beta11_LL") > 0) {
			VINHAC::WeightMonitor monitor("XsTot_LL  O(alf1)_LL");
			xsYFSMonitors["XsTot_LL  O(alf1)_LL"] = monitor;
		}

		if (m_modelWtMonitors.count("beta01_QED + beta11") > 0) {
			VINHAC::WeightMonitor monitor("XsTot_QED O(alf1)_QED");
			xsYFSMonitors["XsTot_QED O(alf1)_QED"] = monitor;
		}

		if (m_modelWtMonitors.count("beta01_EW + beta11") > 0) {
			VINHAC::WeightMonitor monitor("XsTot_EW  O(alf1)_EW");
			xsYFSMonitors["XsTot_EW  O(alf1)_EW"] = monitor;
		}

		if (m_modelWtMonitors.count("xs00") > 0) {
			VINHAC::WeightMonitor monitor("Born");
			xsFixedOrderMonitors["Born"] = monitor;
		}

		if (m_modelWtMonitors.count("xs0gQED + xs1g") > 0) {
			VINHAC::WeightMonitor monitor("O(alpha)_QED");
			xsFixedOrderMonitors["O(alpha)_QED"] = monitor;
		}

		if (m_modelWtMonitors.count("xs0gEW + xs1g") > 0) {
			VINHAC::WeightMonitor monitor("O(alpha)_EW");
			xsFixedOrderMonitors["O(alpha)_EW"] = monitor;
		}

		if (m_modelWtMonitors.count("xs0gQED") > 0) {
			VINHAC::WeightMonitor monitor("Born + virtual_QED + soft");
			xsFixedOrderMonitors["Born + virtual_QED + soft"] = monitor;
		}

		if (m_modelWtMonitors.count("xs0gEW") > 0) {
			VINHAC::WeightMonitor monitor("Born + virtual_EW + soft");
			xsFixedOrderMonitors["Born + virtual_EW + soft"] = monitor;
		}

		if (m_modelWtMonitors.count("xs1g") > 0) {
			VINHAC::WeightMonitor monitor("One real hard photon");
			xsFixedOrderMonitors["One real hard photon"] = monitor;
		}

		mode = 1;
	}

	double mainWt = 1.0;
	for (map<string, double>::iterator it = commonWeights.begin(); it
			!= commonWeights.end(); it++) {
		string name = (*it).first;
		double wt = (*it).second;
		mainWt *= wt;
		m_commonWtMonitors[name].addWeight(wt);
	}
	for (map<string, double>::iterator it = modelWeights.begin(); it
			!= modelWeights.end(); it++) {
		string name = (*it).first;
		double wt = (*it).second;
		m_modelWtMonitors[name].addWeight(wt);

		xsAllMonitors[name].addWeight(wt*mainWt);
	}

	if (m_modelWtMonitors.count("beta00") > 0) {
		xsYFSMonitors["XsTot     O(alf0)"].addWeight(
				newEvent.getModelWeightByName("beta00") * mainWt);
	}

	if (m_modelWtMonitors.count("beta01_LL + beta11_LL") > 0) {
		xsYFSMonitors["XsTot_LL  O(alf1)_LL"].addWeight(
				newEvent.getModelWeightByName("beta01_LL + beta11_LL") * mainWt);
	}

	if (m_modelWtMonitors.count("beta01_QED + beta11") > 0) {
		xsYFSMonitors["XsTot_QED O(alf1)_QED"].addWeight(
				newEvent.getModelWeightByName("beta01_QED + beta11") * mainWt);
	}

	if (m_modelWtMonitors.count("beta01_EW + beta11") > 0) {
		xsYFSMonitors["XsTot_EW  O(alf1)_EW"].addWeight(
				newEvent.getModelWeightByName("beta01_EW + beta11") * mainWt);
	}


	if (m_modelWtMonitors.count("xs00") > 0) {
		xsFixedOrderMonitors["Born"].addWeight(
				newEvent.getModelWeightByName("xs00") * mainWt);
	}

	if (m_modelWtMonitors.count("xs0gQED + xs1g") > 0) {
		xsFixedOrderMonitors["O(alpha)_QED"].addWeight(
				newEvent.getModelWeightByName("xs0gQED + xs1g") * mainWt);
	}

	if (m_modelWtMonitors.count("xs0gEW + xs1g") > 0) {
		xsFixedOrderMonitors["O(alpha)_EW"].addWeight(
				newEvent.getModelWeightByName("xs0gEW + xs1g") * mainWt);
	}

	if (m_modelWtMonitors.count("xs0gQED") > 0) {
		xsFixedOrderMonitors["Born + virtual_QED + soft"].addWeight(
				newEvent.getModelWeightByName("xs0gQED") * mainWt);
	}

	if (m_modelWtMonitors.count("xs0gEW") > 0) {
		xsFixedOrderMonitors["Born + virtual_EW + soft"].addWeight(
				newEvent.getModelWeightByName("xs0gEW") * mainWt);
	}

	if (m_modelWtMonitors.count("xs1g") > 0) {
		xsFixedOrderMonitors["One real hard photon"].addWeight(
				newEvent.getModelWeightByName("xs1g") * mainWt);
	}



	mainWt *= modelWeights["ModelBest"];
	mainWeightMonitor.addWeight(mainWt);
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::eventEvolve() # mainWt: "<<mainWt<<std::endl;

#endif

	return true;
}

void VINHAC::MonitorHandler::makeSnapshot(string snapshotPath) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::makeSnapshot()"<<std::endl;

#endif

	string filePath = snapshotPath + "/XSections.dat";
	string filePathBak = snapshotPath + "/XSections.bak";
	string mv_file_comand = "mv " + filePath + " " + filePathBak;
	system(mv_file_comand.c_str());

	fstream output(filePath.c_str(), fstream::out);

	output << "Best:" << mainWeightMonitor.averageWeight() << ":"
			<< mainWeightMonitor.error() * mainWeightMonitor.averageWeight()
			<< endl;

	output << "# YFS Exponentiation" << endl;
	for (map<string, VINHAC::WeightMonitor>::iterator it = xsYFSMonitors.begin(); it
			!= xsYFSMonitors.end(); ++it) {
		output << (*it).first << ":" << (*it).second.averageWeight() << ":"
				<< (*it).second.error() * (*it).second.averageWeight() << endl;

	}

	output << "# Fixed order" << endl;
	for (map<string, VINHAC::WeightMonitor>::iterator it = xsFixedOrderMonitors.begin(); it
			!= xsFixedOrderMonitors.end(); ++it) {
		output << (*it).first << ":" << (*it).second.averageWeight() << ":"
				<< (*it).second.error() * (*it).second.averageWeight() << endl;

	}

	output << "# All" << endl;

	for (map<string, VINHAC::WeightMonitor>::iterator it = xsAllMonitors.begin(); it
			!= xsAllMonitors.end(); ++it) {
		output << (*it).first << ":" << (*it).second.averageWeight() << ":"
				<< (*it).second.error() * (*it).second.averageWeight() << endl;

	}
	output.close();

	string monfilePath = snapshotPath + "/MonitorsBackup.dat";
	string monfilePathBak = snapshotPath + "/MonitorsBackup.bak";
	string mv_monfile_comand = "mv " + monfilePath + " " + monfilePathBak;
	system(mv_monfile_comand.c_str());

	fstream monoutput(monfilePath.c_str(), fstream::out);

	for (map<string, VINHAC::WeightMonitor>::iterator it =
			m_commonWtMonitors.begin(); it != m_commonWtMonitors.end(); ++it) {
		monoutput << (*it).second.makeSnapshot() << endl;
	}
	monoutput << "#" << endl;
	for (map<string, VINHAC::WeightMonitor>::iterator it =
			m_modelWtMonitors.begin(); it != m_modelWtMonitors.end(); ++it) {
		monoutput << (*it).second.makeSnapshot() << endl;
	}
	monoutput << "#" << endl;
	for (map<string, VINHAC::WeightMonitor>::iterator it = xsYFSMonitors.begin(); it
			!= xsYFSMonitors.end(); ++it) {
		monoutput << (*it).second.makeSnapshot() << endl;
	}

	monoutput << "#" << endl;
	for (map<string, VINHAC::WeightMonitor>::iterator it = xsFixedOrderMonitors.begin(); it
			!= xsFixedOrderMonitors.end(); ++it) {
		monoutput << (*it).second.makeSnapshot() << endl;
	}

	monoutput << "#" << endl;
	for (map<string, VINHAC::WeightMonitor>::iterator it = xsAllMonitors.begin(); it
			!= xsAllMonitors.end(); ++it) {
		monoutput << (*it).second.makeSnapshot() << endl;
	}
	monoutput << "#" << endl;
	monoutput << mainWeightMonitor.makeSnapshot() << endl;

	monoutput.close();
}

void VINHAC::MonitorHandler::restoreState(string snapshotPath) {
	string monfilePath = snapshotPath + "/MonitorsBackup.dat";
	fstream monoutput(monfilePath.c_str(), fstream::in);
	if (monoutput.fail())
		throw VinhacException(
				"Failed to restore monitors state, check the [outputPath]/snapshots/MonitorsBackup.dat file!");

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector<string> line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(line[0]);
		monitor.restoreState(line);
		m_commonWtMonitors[line[0]] = monitor;
	}

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector<string> line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(line[0]);
		monitor.restoreState(line);
		m_modelWtMonitors[line[0]] = monitor;
	}

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector<string> line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(line[0]);
		monitor.restoreState(line);
		xsYFSMonitors[line[0]] = monitor;
	}

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector<string> line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(line[0]);
		monitor.restoreState(line);
		xsFixedOrderMonitors[line[0]] = monitor;
	}

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector<string> line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(line[0]);
		monitor.restoreState(line);
		xsAllMonitors[line[0]] = monitor;
	}

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector<string> line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(line[0]);
		monitor.restoreState(line);
		mainWeightMonitor = monitor;
		break;
	}

	monoutput.close();
}

void VINHAC::MonitorHandler::printCrossSections() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::printCrossSections()"<<std::endl;

#endif

	int width1 = 27;
	int width2 = 20;
	int width3 = 20;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "#";
	}
	cout << endl;
	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		cout << " ";
	}
	cout << "Cross Sections [nb]";
	cout << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "=";
	}
	cout << endl;

	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		cout << " ";
	}
	cout << "Best";
	cout << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "-";
	}
	cout << endl;

	int precision = 13;

	cout.width(width1);
	cout.flags(ios_base::left);
	cout << "Best";

	cout << " ";
	cout.width(width2);
	cout.precision(precision);
	cout << mainWeightMonitor.averageWeight();
	cout << " +/- ";
	cout.width(width3);
	cout.precision(precision);
	cout << mainWeightMonitor.error() * mainWeightMonitor.averageWeight();
	cout << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "=";
	}
	cout << endl;

	if (xsYFSMonitors.size() > 0) {
		for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
			cout << " ";
		}
		cout << "YFS exponentiation ";
		cout << endl;

		for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
			cout << "-";
		}
		cout << endl;

		int precision = 13;
		for (map<string, VINHAC::WeightMonitor>::iterator it =
				xsYFSMonitors.begin(); it != xsYFSMonitors.end(); it++) {

			cout.width(width1);
			cout.flags(ios_base::left);
			cout << (*it).first;

			cout << " ";
			cout.width(width2);
			cout.precision(precision);
			cout << (*it).second.averageWeight();
			cout << " +/- ";
			cout.width(width3);
			cout.precision(precision);
			cout << (*it).second.error() * (*it).second.averageWeight();
			cout << endl;

		}

		for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
			cout << "=";
		}
		cout << endl;
	}

	if (xsFixedOrderMonitors.size() > 0) {
		for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
			cout << " ";
		}
		cout << "Fixed Order ";
		cout << endl;

		for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
			cout << "-";
		}
		cout << endl;

		int precision = 13;
		for (map<string, VINHAC::WeightMonitor>::iterator it =
				xsFixedOrderMonitors.begin(); it != xsFixedOrderMonitors.end(); it++) {

			cout.width(width1);
			cout.flags(ios_base::left);
			cout << (*it).first;

			cout << " ";
			cout.width(width2);
			cout.precision(precision);
			cout << (*it).second.averageWeight();
			cout << " +/- ";
			cout.width(width3);
			cout.precision(precision);
			cout << (*it).second.error() * (*it).second.averageWeight();
			cout << endl;

		}

		for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
			cout << "=";
		}
		cout << endl;
	}

	if (xsAllMonitors.size() > 0) {
		for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
			cout << " ";
		}
		cout << "All ";
		cout << endl;

		for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
			cout << "-";
		}
		cout << endl;

		int precision = 13;
		for (map<string, VINHAC::WeightMonitor>::iterator it =
				xsAllMonitors.begin(); it != xsAllMonitors.end(); it++) {

			cout.width(width1);
			cout.flags(ios_base::left);
			cout << (*it).first;

			cout << " ";
			cout.width(width2);
			cout.precision(precision);
			cout << (*it).second.averageWeight();
			cout << " +/- ";
			cout.width(width3);
			cout.precision(precision);
			cout << (*it).second.error() * (*it).second.averageWeight();
			cout << endl;

		}

		for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
			cout << "=";
		}
		cout << endl;
	}

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "#";
	}
	cout << endl;

}

void VINHAC::MonitorHandler::printAverageWeights() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::printAverageWeights()"<<std::endl;

#endif

	int width1 = 27;
	int width2 = 20;
	int width3 = 20;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "#";
	}
	cout << endl;
	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		cout << " ";
	}
	cout << "Weights";
	cout << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "=";
	}
	cout << endl;

	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		cout << " ";
	}
	cout << "Crude";
	cout << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "-";
	}
	cout << endl;

	int precision = 13;
	for (map<string, VINHAC::WeightMonitor>::iterator it =
			m_commonWtMonitors.begin(); it != m_commonWtMonitors.end(); it++) {

		cout.width(width1);
		cout.flags(ios_base::left);
		cout << (*it).first;

		cout << " ";
		cout.width(width2);
		cout.precision(precision);
		cout << (*it).second.averageWeight();
		cout << " +/- ";
		cout.width(width3);
		cout.precision(precision);
		cout << (*it).second.error() * (*it).second.averageWeight();
		cout << endl;

	}

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "=";
	}
	cout << endl;

	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		cout << " ";
	}
	cout << "Model";
	cout << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "-";
	}
	cout << endl;

	for (map<string, VINHAC::WeightMonitor>::iterator it =
			m_modelWtMonitors.begin(); it != m_modelWtMonitors.end(); it++) {

		cout.width(width1);
		cout << (*it).first;

		cout << " ";
		cout.width(width2);
		cout.precision(precision);
		cout << (*it).second.averageWeight();
		cout << " +/- ";
		cout.width(width3);
		cout.precision(precision);
		cout << (*it).second.error() * (*it).second.averageWeight();
		cout << endl;

	}

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "=";
	}
	cout << endl;

	cout.width(width1);
	cout << "Main Weight";

	cout << " ";
	cout.width(width2);
	cout.precision(precision);
	cout << mainWeightMonitor.averageWeight();
	cout << " +/- ";
	cout.width(width3);
	cout.precision(precision);
	cout << mainWeightMonitor.error() * mainWeightMonitor.averageWeight();
	cout << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		cout << "#";
	}
	cout << endl;

}
VINHAC::MonitorHandler::~MonitorHandler() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::~MonitorHandler()"<<std::endl;

#endif
}

