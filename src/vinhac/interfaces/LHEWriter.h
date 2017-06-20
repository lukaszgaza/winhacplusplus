/*
 * LHEWriter.h
 *
 *  Created on: 2010-11-04
 *      Author: kamil
 */

#ifndef LHEWRITER_H_
#define LHEWRITER_H_

#include <iostream>

namespace VINHAC {

class Manager;
class Event;
class DataSource;
class Particle;

class LHEWriter {
private:
	const Manager& manager;
	std::ostream& os;
	LHEWriter(const Manager& manager,std::ostream& os);

	int getNumberOfParticles(Event& event);
	double getEventWeight(Event& event);
	void printParticle(Particle& particle, int statusPart,
			int mother1Part, int mother2Part,int colorTag=0,int antiColorTag=0);
	friend class Manager;
	friend class QCDHandler;
public:

	virtual ~LHEWriter(){};

	void writeOpening();
	void writeInit();
	void writeEvent(Event& event);
	void writeClosing();
};

}

#endif /* LHEWRITER_H_ */
