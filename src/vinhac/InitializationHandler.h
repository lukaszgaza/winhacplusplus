#include <string>
#include <vector>
#include <exception>
using namespace std;

#ifndef __InitializationHandler_h__
#define __InitializationHandler_h__

#include "BeamHandler.h"
#include "HardProcessHandler.h"
#include "MonitorHandler.h"
#include "ModelHandler.h"
#include "Reader.h"
#include "XercescReader.h"
#include "Manager.h"

namespace VINHAC {
class Manager;
class InitializationHandler;
class HardProcessHandler;
class BeamHandler;
class XercescReader;
class ModelHandler;

}
namespace VINHAC {

//!  Class resposible for initialization process
/*!
	This class do the whole initialization and creates all internall process handlers.
*/
class InitializationHandler {
private:

	Manager* p_Manager;
	HardProcessHandler* p_HardProcessHandler;
	BeamHandler* p_BeamHandler;
	MonitorHandler* p_MonitorHandler;
	ModelHandler* p_ModelHandler;

	Reader* p_InputFileReader;
	DataSource *ds; //!< pointer to object with input parameters

	string m_userFilePath;
	string inputPath;
	string outputPath;

	std::vector<double> readFermionMasses();
	void initializeDS();
	std::map<string, string> initializeBeam(string beam_name);

	bool initializeReader();

public:
	/**
	 * \brief constructor
	 *
	 *	@param ds pointer to data structure which will be filled
	 */
	InitializationHandler(DataSource *ds);

	/**
	 * \brief initializes all internal handlers
	 */
	void initializeAllHandlers();

	/**
	 * \brief prints echo if readed input
	 */
	bool inputEcho();

	/**
	 * \brief sets path to user file
	 *
	 * @param path path to set
	 */
	inline void setUserFilePath(string path) {
		m_userFilePath = path;
	}

	/**
	 * \brief sets path where it should look for input files
	 *
	 * @param path path to set
	 */
	inline void setInputPath(string path) {
		inputPath = path;
	}

	/**
	 *	\brief sets path where should be placed output from generator
	 *
	 * @param path path to set
	 */
	inline void setOutputPath(string path) {
		outputPath = path;
	}
	~InitializationHandler();

	/**
	 * \brief returns pointer to HardProcessHandler object
	 *
	 * @return pointer to HardProcessHandler
	 */
	inline HardProcessHandler* getHardProcessHandler() {
		return p_HardProcessHandler;
	}

	/**
	 * \brief returns pointer to BeamHandler object
	 *
	 * @return pointer to BeamHandler
	 */
	inline BeamHandler* getBeamHandler() {
		return p_BeamHandler;
	}

	/**
	 * \brief returns pointer to MonitorHandler object
	 *
	 * @return pointer to MonitorHandler
	 */
	inline MonitorHandler* getMonitorHandler() {
		return p_MonitorHandler;
	}

	/**
	 * \brief returns pointer to ModelHandler object
	 *
	 * @return pointer to ModelHandler
	 */
	inline ModelHandler* getModelHandler() {
		return p_ModelHandler;
	}

	/**
	 * \brief sets pointer to Manager
	 *
	 * @param manager pointer to set
	 */
	inline void setManager(Manager* manager) {
		p_Manager = manager;
	}
};
}

#endif
