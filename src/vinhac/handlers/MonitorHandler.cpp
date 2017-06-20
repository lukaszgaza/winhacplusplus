/*
 * \file MonitorHandler.cpp
 * \brief
 * \date Created on: 2009-07-07
 * \author     Author: siudmy
 */
#include "MonitorHandler.h"
#include <fstream>
#include "../../utils/Functions.h"
#include "../core/VinhacException.h"
#include "../core/Event.h"
#include "../core/enum/PolarizationDefinition.h"
#include "../util/Commons.h"
#include "../input/DataSource.h"

VINHAC::MonitorHandler::MonitorHandler() :
	numberOfEvents(0),
	mainWeightMonitor() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::MonitorHandler()"<<std::endl;

#endif

	setHandlerName("MonitorHandler");
	m_commonWtMonitors.clear();
	m_modelWtMonitors.clear();
	xsMonitors.clear();
	mode = 0;
}

bool VINHAC::MonitorHandler::eventEvolve(VINHAC::Event& newEvent) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::eventEvolve()"<<std::endl;

#endif

	const map<WeightName::name, double>& commonWeights =
			newEvent.getCommonWeightsMap();
	const map<pair<WeightName::name, PolarizationName::name> , double>
			& modelWeights = newEvent.getModelWeightsMap();

	double mainWt = 1.0;
	for (map<WeightName::name, double>::const_iterator it =
			commonWeights.begin(); it != commonWeights.end(); it++) {
		WeightName::name name = it->first;
		double wt = it->second;
		mainWt *= wt;
		m_commonWtMonitors[name].addWeight(wt);
	}
	for (map<pair<WeightName::name, PolarizationName::name> , double>::const_iterator
			it = modelWeights.begin(); it != modelWeights.end(); it++) {
		pair<WeightName::name, PolarizationName::name> name = it->first;
		double wt = it->second;
		m_modelWtMonitors[name].addWeight(wt);

		if (name.first > 100) {
			xsMonitors[name].addWeight(wt * mainWt);
		}
	}

	mainWt *= modelWeights.find(DataSource::get().mainWeight)->second;
	mainWeightMonitor.addWeight(mainWt);

	numberOfEvents++;
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

	output << "Best:" << mainWeightMonitor.averageWeight(numberOfEvents) << ":"
			<< mainWeightMonitor.error(numberOfEvents) * mainWeightMonitor.averageWeight(numberOfEvents)
			<< endl;

	output << "# All" << endl;

	for (map<pair<WeightName::name, PolarizationName::name> ,
			VINHAC::WeightMonitor>::iterator it = xsMonitors.begin(); it
			!= xsMonitors.end(); ++it) {
		output << Commons::mapWeightPolarizationToString((*it).first) << ":"
				<< (*it).second.averageWeight() << ":" << (*it).second.error()
				* (*it).second.averageWeight() << endl;

	}
	output.close();

	string monfilePath = snapshotPath + "/MonitorsBackup.dat";
	string monfilePathBak = snapshotPath + "/MonitorsBackup.bak";
	string mv_monfile_comand = "mv " + monfilePath + " " + monfilePathBak;
	system(mv_monfile_comand.c_str());

	fstream monoutput(monfilePath.c_str(), fstream::out);

	for (map<WeightName::name, VINHAC::WeightMonitor>::iterator it =
			m_commonWtMonitors.begin(); it != m_commonWtMonitors.end(); ++it) {
		monoutput << it->second.makeSnapshot(
				WeightName::mapEnumToString(it->first)) << endl;
	}
	monoutput << "#" << endl;
	for (map<pair<WeightName::name, PolarizationName::name> ,
			VINHAC::WeightMonitor>::iterator it = m_modelWtMonitors.begin(); it
			!= m_modelWtMonitors.end(); ++it) {
		monoutput << it->second.makeSnapshot(
				Commons::mapWeightPolarizationToString(it->first)) << endl;
	}
	monoutput << "#" << endl;
	for (map<pair<WeightName::name, PolarizationName::name> ,
			VINHAC::WeightMonitor>::iterator it = xsMonitors.begin(); it
			!= xsMonitors.end(); ++it) {
		monoutput << it->second.makeSnapshot(
				Commons::mapWeightPolarizationToString(it->first)) << endl;
	}
	monoutput << "#" << endl;
	monoutput
			<< mainWeightMonitor.makeSnapshot(
					Commons::mapWeightPolarizationToString(
							DataSource::get().mainWeight)) << endl;

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
		vector < string > line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(WeightName::mapStringToEnum( line[0]));
		monitor.restoreState(line);
		m_commonWtMonitors[WeightName::mapStringToEnum(line[0])] = monitor;
	}

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector < string > line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(WeightName::mapStringToEnum( line[0]));
		monitor.restoreState(line);
		m_modelWtMonitors[Commons::mapStringToWeightPolarization(line[0])]
				= monitor;
	}

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector < string > line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(WeightName::mapStringToEnum( line[0]));
		monitor.restoreState(line);
		xsMonitors[Commons::mapStringToWeightPolarization(line[0])] = monitor;
	}

	while (!monoutput.eof()) {
		string tmp;
		monoutput >> tmp;
		if (tmp.at(0) == '#')
			break;
		vector < string > line;
		Tokenize(tmp, line, ":");
		if (line.size() != 6)
			continue;
		VINHAC::WeightMonitor monitor(WeightName::mapStringToEnum( line[0]));
		monitor.restoreState(line);
		mainWeightMonitor = monitor;
		break;
	}

	monoutput.close();
}

void VINHAC::MonitorHandler::printCrossSections(ostream& out) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::printCrossSections()"<<std::endl;

#endif

	int width1 = 50;
	int width2 = 20;
	int width3 = 20;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "#";
	}
	out << endl;
	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		out << " ";
	}
	out << "Cross Sections [nb]";
	out << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "=";
	}
	out << endl;

	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		out << " ";
	}
	out << "Best";
	out << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "-";
	}
	out << endl;

	int precision = 13;

	out.width(width1);
	out.flags(ios_base::left);
	out
			<< Commons::mapWeightPolarizationToString(
					DataSource::get().mainWeight);

	out << " ";
	out.width(width2);
	out.precision(precision);
	out << mainWeightMonitor.averageWeight(numberOfEvents);
	out << " +/- ";
	out.width(width3);
	out.precision(precision);
	out << mainWeightMonitor.error(numberOfEvents) * mainWeightMonitor.averageWeight(numberOfEvents);
	out << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "=";
	}
	out << endl;

	if (xsMonitors.size() > 0) {
		for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
			out << " ";
		}
		out << "All ";
		out << endl;

		for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
			out << "-";
		}
		out << endl;

		int precision = 13;
		for (map<pair<WeightName::name, PolarizationName::name> ,
				VINHAC::WeightMonitor>::iterator it = xsMonitors.begin(); it
				!= xsMonitors.end(); it++) {

			out.width(width1);
			out.flags(ios_base::left);
			out << Commons::mapWeightPolarizationToString((*it).first);

			out << " ";
			out.width(width2);
			out.precision(precision);
			out << (*it).second.averageWeight();
			out << " +/- ";
			out.width(width3);
			out.precision(precision);
			out << (*it).second.error() * (*it).second.averageWeight();
			out << endl;

		}

		for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
			out << "=";
		}
		out << endl;
	}

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "#";
	}
	out << endl;

}

void VINHAC::MonitorHandler::printAverageWeights(ostream& out) {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::printAverageWeights()"<<std::endl;

#endif

	int width1 = 50;
	int width2 = 20;
	int width3 = 20;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "#";
	}
	out << endl;
	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		out << " ";
	}
	out << "Weights";
	out << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "=";
	}
	out << endl;

	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		out << " ";
	}
	out << "Crude";
	out << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "-";
	}
	out << endl;

	int precision = 13;
	for (map<WeightName::name, VINHAC::WeightMonitor>::iterator it =
			m_commonWtMonitors.begin(); it != m_commonWtMonitors.end(); it++) {

		out.width(width1);
		out.flags(ios_base::left);
		out << WeightName::mapEnumToString((*it).first);

		out << " ";
		out.width(width2);
		out.precision(precision);
		out << (*it).second.averageWeight();
		out << " +/- ";
		out.width(width3);
		out.precision(precision);
		out << (*it).second.error() * (*it).second.averageWeight();
		out << endl;

	}

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "=";
	}
	out << endl;

	for (int i = 0; i < (width1 + width2 + width3) / 2; ++i) {
		out << " ";
	}
	out << "Model";
	out << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "-";
	}
	out << endl;

	for (map<pair<WeightName::name, PolarizationName::name> ,
			VINHAC::WeightMonitor>::iterator it = m_modelWtMonitors.begin(); it
			!= m_modelWtMonitors.end(); it++) {

		out.width(width1);
		out << Commons::mapWeightPolarizationToString((*it).first);

		out << " ";
		out.width(width2);
		out.precision(precision);
		out << (*it).second.averageWeight();
		out << " +/- ";
		out.width(width3);
		out.precision(precision);
		out << (*it).second.error() * (*it).second.averageWeight();
		out << endl;

	}

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "=";
	}
	out << endl;

	out.width(width1);
	out << "Main Weight";

	out << " ";
	out.width(width2);
	out.precision(precision);
	out << mainWeightMonitor.averageWeight(numberOfEvents);
	out << " +/- ";
	out.width(width3);
	out.precision(precision);
	out << mainWeightMonitor.error(numberOfEvents) * mainWeightMonitor.averageWeight(numberOfEvents);
	out << endl;

	for (int i = 0; i < width1 + width2 + width3 + 5; ++i) {
		out << "#";
	}
	out << endl;

}
VINHAC::MonitorHandler::~MonitorHandler() {
#ifdef DEBUG_FLAG

	std::cout<<"DEBUG: MonitorHandler::~MonitorHandler()"<<std::endl;

#endif
}

