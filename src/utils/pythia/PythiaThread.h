/*
 * PythiaThread.h
 *
 *  Created on: 04-08-2011
 *      Author: sobol
 */

#ifndef PYTHIATHREAD_H_
#define PYTHIATHREAD_H_

#include <boost/thread/thread.hpp>

class PythiaThread {
private:
	Pythia8::Pythia& pythia;
	boost::mutex& myLock;
	boost::mutex& concurentLock;
	bool& stopFlag;
	string fifoPath;
public:
	PythiaThread(Pythia8::Pythia& pythia,
			boost::mutex& myLock, boost::mutex& concurentLock,
			bool& stopFlag,
			string fifoPath) :
		pythia(pythia),
				myLock(myLock), concurentLock(concurentLock), stopFlag(
						stopFlag), fifoPath(fifoPath) {
	}

	void operator()() {

		int iEvent = 0;

		Pythia8::LHAupLHEF lhaUp(fifoPath.c_str());
		// Initialize Les Houches Event File run.
		pythia.init(&lhaUp);
		// List initialization information.
		//pythia.settings.listChanged();
		//pythia.settings.listChanged();
		//cout<<"#################"<<endl;
		//pythia.settings.listAll();

		while (true) {
			myLock.lock();
			if (stopFlag)
				break;
			//cout << "Begin event " << iEvent << endl;

			if(!pythia.next()) {
				// If failure because reached end of file then exit event loop
				if (pythia.info.atEndOfFile()) {
					//cout << "Info: end of Les Houches file reached" << endl;
					break;
				}

				//cout << "Warning: event " << iEvent << "failed" << endl;
				//continue;
			}

			// List first event: Les Houches, hard process and complete
			//pythia.LHAeventList();
			//pythia.info.list();
			//pythia.process.list();
			//pythia.event.list();
			//cout<<"##"<<pythia.info.weight()<<"##"<<endl;
			//cout<<pythia.event.size()<<endl;

			++iEvent;

			concurentLock.unlock();
		}
		concurentLock.unlock();
	}

	~PythiaThread(){
	}
};

#endif /* PYTHIATHREAD_H_ */
