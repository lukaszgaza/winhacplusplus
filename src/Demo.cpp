
#include <iostream>
#include "LHAPDF/LHAPDF.h"

using namespace std;

//CMake automatically checks if there is ROOT on your system and sets proper flag
#ifdef ROOT_HIST
#include "TCanvas.h"
#include "TApplication.h"
#endif

#include "vinhac/WINHAC.h"
#include "vinhac/VinhacException.h"
#include "vinhac/DemoAnalysisHandler.h"
#include "vinhac/Event.h"

int main(int argc, char* argv[]) {

	string userPath;
	string inputPath;
	string outputPath;
	double eventsNb;

	/* Reading command line parameters */
	if(argc <3 ){
		cout<<"Usage : Demo [input path] [output path] <number of events> <user file path>"<<endl;
		return 1;
	}
	if (argc >= 3) {
		/* reading path to input and output directory */
		inputPath = argv[1];
		outputPath = argv[2];
	}
	if(argc >3 ){
		/* reading number of events to generate */
		eventsNb = atof(argv[3]);
	} else {
		/* default number of events */
		eventsNb = 10000;
	}
	if(argc >4){
		/* path to UserFile.xml-like file */
		userPath = argv[4];
	}


	try{
	// Creation of WINHAC instance
	VINHAC::WINHAC generator;

	// Creation of DemoAnalysisHandler instance
	VINHAC::DemoAnalysisHandler demoAnalysisHandler;

	// Setting path to input, output and UserFile.xml-like file
	generator.setInputPath(inputPath);
	generator.setOutputPath(outputPath);
	generator.setUserFilePath(userPath);


	// Initialization of WINHAC
	generator.initializeGenerator(argc, argv);

	// Adding own DemoAnalysisHandler to generator's ProcessHandlers vector
	// This step should be done after initialization
	generator.addProcessHandler(&demoAnalysisHandler);

	// The main events generation loop
	cout<<"######## GENERATION #########"<<endl;
	double i=0;
	for(i=0; i<eventsNb; i++){
		VINHAC::Event event = generator.generateEvent();


		// Generator makes snapshot
		if(fmod(i+1,100000)==0){
			generator.makeSnapshot(outputPath+"/snapshots",i+1);
			cout<<i+1<<endl;
		}
		// Here could be some user's code handling event

	}
	cout<<"############################"<<endl;
	cout<<"Generated: "<<i<<" events"<<endl;

	// Generator makes snapshot
	generator.makeSnapshot(outputPath,i);
	// Here are some outputs after generation
	cout.precision(15);
	// The main cross section with its error
	cout<<"Xsection = "<<generator.getFinalXsection()<<" +/- "<<generator.getFinalXsectionError()<<endl;

	// This method prints weights and cross sections
	generator.printSummary();

	// This method causes that demoAnalysisHandler writes ROOT's histograms to output directory
	demoAnalysisHandler.makeSnapshot(outputPath);

	// Below you can find ROOT application drawing collected histograms
#ifdef ROOT_HIST
	TApplication theApp("theApp",&argc,argv);
	TCanvas kanwa("kanwa","Histograms",5,5,1024,768);
	kanwa.Divide(2,3);
	kanwa.cd(1);
	demoAnalysisHandler.getWrapidityHist().Draw();
	kanwa.cd(2);
	demoAnalysisHandler.getCosThetaHist().Draw();
	kanwa.cd(3);
	demoAnalysisHandler.getLeptonPTHist().Draw();
	kanwa.cd(4);
	demoAnalysisHandler.getLeptonPseudorapidityHist().Draw();
	kanwa.cd(5);
	demoAnalysisHandler.getWTransverseMass().Draw();
	kanwa.cd(6);
	demoAnalysisHandler.getPhotonMultiplicity().Draw();

	TCanvas kanwa2("kanwa2","Histograms2",5,5,1024,768);
	kanwa2.cd(1);
	demoAnalysisHandler.getPhotonsEnergy().Draw();

	TCanvas kanwa4("kanwa4","Histograms4",5,5,1024,768);
	kanwa4.cd(1);
	demoAnalysisHandler.getHardPhotonPt().Draw();

	TCanvas kanwa3("kanwa3","Histograms3",5,5,1024,768);
	kanwa3.cd(1);
	demoAnalysisHandler.getHardPhotonEnergy().Draw();

	TCanvas kanwa5("kanwa5","Histograms5",5,5,1024,768);
	kanwa5.cd(1);
	demoAnalysisHandler.getHardPhotonEta().Draw();


	theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
	theApp.Run();

#endif
	} catch(VINHAC::VinhacException e){
		//In case of problems WINHAC throws VinhacException with message what's wrong.
		cout<<"EXCEPTION : "<<e.getMessage()<<endl;
	}
}

