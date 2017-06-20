#include <string>
#include <vector>
#include <exception>
#include "../handlers/FrameTransformationHandler.h"
using namespace std;

#ifndef __InitializationHandler_h__
#define __InitializationHandler_h__


namespace VINHAC {
class Manager;
class InitializationHandler;
class HardProcessHandler;
class BeamHandler;
class XercescReader;
class ModelHandler;
class MonitorHandler;
class Reader;
class DataSource;
class QCDHandler;
class FrameTransformationHandler;
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
	QCDHandler* p_QCDHandler;
	FrameTransformationHandler* p_FrameTransformationHandler;

	Reader* p_InputFileReader;

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
	InitializationHandler();

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
	 * \brief returns pointer to FrameTransformationHandler object
	 *
	 * @return pointer to FrameTransformationHandler
	 */
	inline FrameTransformationHandler* getFrameTransformationHandler() {
		return p_FrameTransformationHandler;
	}

	/**
	 * \brief returns pointer to QCDHandler object
	 *
	 * @return pointer to QCDHandler
	 */
	inline QCDHandler* getQCDHandler() {
		return p_QCDHandler;
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
