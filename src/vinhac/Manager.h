#include <string>
#include <vector>
#include <exception>
#include <iostream>
using namespace std;

#ifndef __Manager_h__
#define __Manager_h__

#include "DataSource.h"
#include "HardProcessHandler.h"
#include "BeamHandler.h"
#include "MonitorHandler.h"
#include "InitializationHandler.h"
#include "Event.h"
#include "ProcessHandler.h"
#include <map>
#include "../utils/CLHEP/Random/Random/RanluxEngine.h"
#include "../utils/CLHEP/Random/Random/RandomEngine.h"

#include <fstream>
#include <ctime>

namespace VINHAC {
class DataSource;
class HardProcessHandler;
class BeamHandler;
class InitializationHandler;
class MonitorHandler;
class Event;
class Manager;
class BornHandler;
}

namespace VINHAC {

/*! \brief Manager is a base class for the construction of Drell-Yan Monte Carlo
 *  generators such as ZINHAC++ and WINHAC++. The aim of this class is to manage
 *  the whole generation of Drell-Yan events.
 *  \li read data from default and user files
 *  \li initialise all handlers
 */

class Manager {
private:
	string m_userFilePath; //!<  user file path needed for the initialisation phase
	string input_Path;
	string output_Path;

	virtual void PrintLogo() = 0; //!<  method which prints the logo of generator

	int m_snapshotNb; //!< number of snapshot in the run

	std::vector<VINHAC::ProcessHandler*> _handlers;
	VINHAC::HardProcessHandler* p_HardProcessHandler;
	VINHAC::BeamHandler* p_BeamHandler;
	VINHAC::InitializationHandler* p_InitializationHandler;
	VINHAC::MonitorHandler* p_MonitorHandler;
	VINHAC::BornHandler* p_BornHandler;
	CLHEP::RanluxEngine * p_randomGenerator;

	DataSource ds; //!< object with input parameters

public:
	/**
	 * \brief constructor
	 *
	 * @param aInputStream path to folder with input files
	 */
	Manager(string aInputStream);

	Manager();

	~Manager();

	//! Initialisation of the random number generator used by all handlers
	void initializeRandomNumberGenerator();

	/*! Initialisation of the generator, creating handlers and putting
	 *  them into handlers vector */
	void initializeGenerator(int argc, char* argv[]);

	//! Adds \param handler to handlers vector
	void
	addProcessHandler(VINHAC::ProcessHandler* handler);

	//! Removes a ProcessHandler from handlers vector
	void removeProcessHandler(string Name);

	//! Generates one event
	/**
	 * @return generated event object
	 */
	Event generateEvent();

	//! Generates \param  NbEvents events
	bool generateEvents(double NbEvents);



	//! Save the most important information from the run
	/**
	 * @param snapshotPath path where snapshot should be saved
	 * @param eventNumber current event number
	 */
	void makeSnapshot(string snapshotPath, int eventNumber);

	/**
	 * \brief checks semaphore flag in "semaphore" file
	 *
	 * Semaphore file should contain: continue - when generator should just generate,
	 * resume - when generator should be resumed from some point of previous generation, to do this
	 * the snapshot shoud be placed in proper folder, stop - stops generation
	 * @param semaforFile path to file
	 * @param eventNb events generated
	 */
	bool checkSemafor(string semaforFile, double eventNb);

	/**
	 * \brief checks semaphore flag in "semaphore" file
	 *
	 * This method is used when generator is resumed
	 * @param semaforFile path to file
	 * @param snapshotPath path containing snapshot
	 * @param i reference to variable which stores events number
	 * @param nbEvents reference to variable which stores events to be generated
	 */
	bool checkSemafor(string semaforFile, string snapshotPath, double& i,
			double& nbEvents);

	//! Prints output
	void printSummary();

	/**
	 * \brief prints crude cross section multipied by weights
	 */
	void printCrudeXsectionTimesWt();

	/**
	 * \brief returns final cross section
	 *
	 * @return final cross section
	 */
	double getFinalXsection();

	/**
	 * \brief returns final cross section error
	 *
	 * @return final cross section error
	 */
	double getFinalXsectionError();

	/**
	 * \brief returns crude cross section
	 *
	 * @return crude cross section
	 */
	double getCrudeXsection();

	//Getters and setters

	/**
	 * \brief sets factor which decides how many times is snapshot done
	 *
	 *	It is currently unused
	 * @param nb factor
	 */
	inline void setSnapshotNb(int nb) {
		m_snapshotNb = nb;
	}

	/**
	 * \brief sets path to user file
	 *
	 * @param path path to user file
	 */
	inline void setUserFilePath(string path) {
		m_userFilePath = path;
	}

	/**
	 * \brief returns user file path
	 *
	 * @return user file path
	 */

	inline string getUserFilePath() {
		return m_userFilePath;
	}

	/**
	 * \brief sets input folder path
	 *
	 * @param name input folder path
	 */
	inline void setInputPath(string name) {
		input_Path = name;
	}

	/**
	 * \brief returns input folder path
	 *
	 * @return input folder path
	 */
	inline string getInputPath() {
		return input_Path;
	}

	/**
	 * \brief sets path to output folder
	 *
	 * @param name path to output folder
	 */
	inline void setOutputPath(string name) {
		output_Path = name;
	}

	/**
	 * \brief returns path to output folder
	 *
	 * @return path to output folder
	 */
	inline string getOutputPath() {
		return output_Path;
	}

	/**
	 * \brief returns pointer to InitializationHandler
	 *
	 * @return pointer to InitializationHandler
	 */
	inline VINHAC::InitializationHandler * const GetInitHandler() const {
		return p_InitializationHandler;
	}

	/**
	 * \brief returns pointer to HardProcessHandler
	 *
	 * @return pointer to HardProcessHandler
	 */
	inline VINHAC::HardProcessHandler* getHardProcessHandler() {
		return this->p_HardProcessHandler;
	}

	/**
	 * \brief returns pointer to MonitorHandler
	 *
	 * @return pointer to MonitorHandler
	 */
	inline VINHAC::MonitorHandler* getMonitorHandler() {
		return this->p_MonitorHandler;
	}

	/**
	 * \brief returns pointer to random number generator
	 *
	 * @return pointer to random number generator
	 */
	inline TRND * getRandomNumberGenerator() {
		return this->p_randomGenerator;
	}

};
}

#endif
