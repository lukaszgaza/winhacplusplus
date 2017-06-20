#include <iostream>
#include <fstream>
#include <cstdlib>
#include <boost/thread/thread.hpp>

using namespace std;

#include "../vinhac/core/WINHAC.h"
#include "../vinhac/core/Event.h"
#include "../vinhac/core/VinhacException.h"
#include "Pythia.h"
#include "../utils/pythia/EventAnalizer.h"
#include "../utils/pythia/RandomEngineAdapter.h"
#include "../utils/pythia/PythiaThread.h"

//CMake automatically checks if there is ROOT on your system and sets proper flag
#ifdef ROOT_HIST
#include "TCanvas.h"
#include "TApplication.h"
#endif

//#define PYTHIA_ONLY

class HistogramThread {
private:
#ifndef PYTHIA_ONLY
	VINHAC::WINHAC& winhac;
#else
	Pythia8::Pythia& pythiaHard;
#endif
	Pythia8::Pythia& pythia;
	EventAnalizer& eventAnalizer;
	boost::mutex& histogramLock;
	boost::mutex& winhacLock;
	bool& stopFlag;
public:
	HistogramThread(
#ifndef PYTHIA_ONLY
			VINHAC::WINHAC& winhac
#else
			Pythia8::Pythia& pythiaHard
#endif
			, Pythia8::Pythia& pythia, EventAnalizer& eventAnalizer,
			boost::mutex& histogramLock, boost::mutex& winhacLock,
			bool& stopFlag) :
#ifndef PYTHIA_ONLY
		winhac(winhac)
#else
		pythiaHard(pythiaHard)
#endif
	, pythia(pythia), eventAnalizer(eventAnalizer), histogramLock(histogramLock),
				winhacLock(winhacLock), stopFlag(stopFlag) {
	}

	void operator()() {
		while (true) {
			histogramLock.lock();
			if (stopFlag)
							break;
			double weight = pythia.info.weight();

			eventAnalizer.analize(pythia.event,weight);


			winhacLock.unlock();
		}
		winhacLock.unlock();
		histogramLock.unlock();
	}
};

class WinhacThread {
private:
#ifndef PYTHIA_ONLY
	VINHAC::WINHAC& winhac;
#else
	Pythia8::Pythia& pythiaHard;
#endif
	HistogramThread& histogramThread;
	boost::mutex& winhacLock;
	boost::mutex& pythiaLock;
	string outputPath;
	double eventsNb;
	bool& stopFlag;
public:
	WinhacThread(
#ifndef PYTHIA_ONLY
			VINHAC::WINHAC& winhac
#else
		Pythia8::Pythia& pythiaHard
#endif
			, HistogramThread& histogramThread,
			boost::mutex& winhacLock, boost::mutex& pythiaLock, bool& stopFlag) :
#ifndef PYTHIA_ONLY
		winhac(winhac)
#else
		pythiaHard(pythiaHard)
#endif
		, histogramThread(histogramThread),
				winhacLock(winhacLock), pythiaLock(pythiaLock), stopFlag(
						stopFlag) {
	}

	void operator()() {
#ifndef PYTHIA_ONLY
		fstream lhefile;
		lhefile.open("SamplePipe.lhe", fstream::out);
		VINHAC::LHEWriter lheWriter = winhac.getLHEWriter(lhefile);
					lheWriter.writeOpening();
					lheWriter.writeInit();
#else
		// Create an LHAup object that can access relevant information in pythia.
		Pythia8::LHAupFromPYTHIA8 myLHA(&pythiaHard.process, &pythiaHard.info);

		// Open a file on which LHEF events should be stored, and write header.
		myLHA.openLHEF("SamplePipe.lhe");
		pythiaHard.init();
		// Store initialization info in the LHAup object.
		myLHA.setInit();
		// Write out this initialization info on the file.
		myLHA.initLHEF();
#endif

		double i = 0;
		for (i = 0; i <= eventsNb; i++) {

			winhacLock.lock();
			if (i == eventsNb) {
				stopFlag = true;
#ifndef PYTHIA_ONLY
				lheWriter.writeClosing();
				lhefile.close();
#else
				// Write endtag. Overwrite initialization info with new cross sections.
				myLHA.closeLHEF(false);
#endif
				pythiaLock.unlock();
				break;
			}
#ifndef PYTHIA_ONLY
			VINHAC::Event event = winhac.generateEvent();
			//std::cout<<event.printInLabFrame()<<endl;
			lheWriter.writeEvent(event);
#else
			pythiaHard.next();
			// Store event info in the LHAup object.
			myLHA.setEvent();

			    // Write out this event info on the file.
			myLHA.eventLHEF();
#endif
			// Generator makes snapshot
			if (fmod(i + 1, 1000) == 0) {
#ifndef PYTHIA_ONLY
				winhac.makeSnapshot(outputPath + "/snapshots", i + 1);
#endif
				cout << i + 1 << endl;
			}


			pythiaLock.unlock();
		}
		pythiaLock.unlock();
		winhacLock.unlock();
	}

	inline void setEventsNumber(double eventsNb) {
		this->eventsNb = eventsNb;
	}
	inline void setOutputPath(string outputPath) {
		this->outputPath = outputPath;
	}
};



int main(int argc, char* argv[]) {

	string pythiaXMLdocdir;
	string userPath;
	string inputPath;
	string outputPath;
	double eventsNb;

	/* Reading command line parameters */
	if (argc < 4) {
		cout
				<< "Usage : DemoPythia [Pythia xmldoc dir] [input path] [output path] <number of events> <user file path>"
				<< endl;
		return 1;
	}
	if (argc >= 4) {
		/* reading path to input and output directory */
		pythiaXMLdocdir = argv[1];
		inputPath = argv[2];
		outputPath = argv[3];
	}
	if (argc > 4) {
		/* reading number of events to generate */
		eventsNb = atof(argv[4]);
	} else {
		/* default number of events */
		eventsNb = 10000;
	}
	if (argc > 5) {
		/* path to UserFile.xml-like file */
		userPath = argv[5];
	}

	try {
#ifndef PYTHIA_ONLY
		// Creation of WINHAC instance
		VINHAC::WINHAC generator;
		// Setting path to input, output and UserFile.xml-like file
		generator.setInputPath(inputPath);
		generator.setOutputPath(outputPath);
		generator.setUserFilePath(userPath);
		// Initialization of WINHAC
		generator.initializeGenerator(argc, argv);
		RandomEngineAdapter randomEngine(*generator.getRandomNumberGenerator());
#else
		Pythia8::Pythia pythiaHard(pythiaXMLdocdir);
#endif
		// Generator
		Pythia8::Pythia pythia(pythiaXMLdocdir);
#ifndef PYTHIA_ONLY
		pythia.setRndmEnginePtr(&randomEngine);
#endif
		pythia.readFile(inputPath+"/pythia/pythiaCommon.cmnd");
		pythia.readFile(inputPath+"/pythia/pythiaSetup.cmnd");
#ifdef PYTHIA_ONLY
		pythiaHard.readFile(inputPath+"/pythia/pythiaCommon.cmnd");
		pythiaHard.readFile(inputPath+"/pythia/pythiaFullCommon.cmnd");
		pythiaHard.readFile(inputPath+"/pythia/pythiaFullSetup.cmnd");
#endif

		system("mkfifo SamplePipe.lhe");

		bool stopFlag = false;
		boost::mutex winhacLock;
		winhacLock.lock();
		boost::mutex pythiaLock;
		pythiaLock.lock();
		boost::mutex histogramLock;
		histogramLock.lock();

		EventAnalizer eventAnalizer;

		HistogramThread histogramThread(
#ifndef PYTHIA_ONLY
				generator
#else
				pythiaHard
#endif
				, pythia, eventAnalizer, histogramLock,
				winhacLock, stopFlag);
		WinhacThread winhacThread(
#ifndef PYTHIA_ONLY
				generator
#else
				pythiaHard
#endif
				, histogramThread, winhacLock,
				pythiaLock, stopFlag);
		winhacThread.setEventsNumber(eventsNb);
		winhacThread.setOutputPath(outputPath);
		PythiaThread pythiaThread(pythia, pythiaLock,
				histogramLock, stopFlag, "SamplePipe.lhe");

		// The main events generation loop
		cout << "######## GENERATION #########" << endl;
		boost::thread thread1(winhacThread);
		boost::thread thread2(pythiaThread);
		boost::thread thread3(histogramThread);

		winhacLock.unlock();

		thread1.join();
		thread2.join();
		thread3.join();

		cout << "############################" << endl;
		cout << "Generated: " << eventsNb << " events" << endl;
		system("rm SamplePipe.lhe");
#ifdef PYTHIA_ONLY
		eventAnalizer.scaleHistograms(pythiaHard.info.sigmaGen());
		cout.precision(15);
		cout<<"## Pythia cross section: "<< pythiaHard.info.sigmaGen()
				<<" +/- "<<pythiaHard.info.sigmaErr()<< " ##"<<endl;
#else
		// Generator makes snapshot
		generator.makeSnapshot(outputPath, eventsNb);
		//eventAnalizer.scaleHistograms(generator.getFinalXsection());
		cout.precision(15);
		// The main cross section with its error
		cout << "Xsection = " << generator.getFinalXsection() << " +/- "
				<< generator.getFinalXsectionError() << endl;
		// This method prints weights and cross sections
		generator.printSummary();
#endif
		eventAnalizer.makeSnapshot(outputPath);


		// Below you can find ROOT application drawing collected histograms
	#ifdef ROOT_HIST

	/*	TApplication theApp("theApp",&argc,argv);
		TCanvas kanwa("kanwa","Histograms",5,5,1024,768);
		kanwa.Divide(2,3);
		kanwa.cd(1);
		eventAnalizer.getWrapidityHist().Draw();
		kanwa.cd(2);
		eventAnalizer.getCosThetaHist().Draw();
		kanwa.cd(3);
		eventAnalizer.getLeptonPTHist().Draw();
		kanwa.cd(4);
		eventAnalizer.getLeptonPseudorapidityHist().Draw();
		kanwa.cd(5);
		eventAnalizer.getWTransverseMass().Draw();
		kanwa.cd(6);
		eventAnalizer.getPhotonMultiplicity().Draw();

		TCanvas kanwa2("kanwa2","Histograms2",5,5,1024,768);
		kanwa2.cd(1);
		eventAnalizer.getPhotonsEnergy().Draw();

		TCanvas kanwa4("kanwa4","Histograms4",5,5,1024,768);
		kanwa4.cd(1);
		eventAnalizer.getHardPhotonPt().Draw();

		TCanvas kanwa3("kanwa3","Histograms3",5,5,1024,768);
		kanwa3.cd(1);
		eventAnalizer.getHardPhotonEnergy().Draw();

		TCanvas kanwa5("kanwa5","Histograms5",5,5,1024,768);
		kanwa5.cd(1);
		eventAnalizer.getHardPhotonEta().Draw();

		TCanvas kanwa6("kanwa6","Histograms6",5,5,1024,768);
		kanwa6.cd(1);
		eventAnalizer.getWPT().Draw();

		TCanvas kanwa7("kanwa7","Histograms7",5,5,1024,768);
		kanwa7.cd(1);
		eventAnalizer.getWMass().Draw();


		theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
		theApp.Run();*/

	#endif
		return 0;

	} catch (VINHAC::VinhacException e) {
		//In case of problems WINHAC throws VinhacException with message what's wrong.
		cout << "EXCEPTION : " << e.getMessage() << endl;
	}

#ifndef PYTHIA_ONLY
	system("rm SamplePipe.lhe");
#endif

	return 0;
}

