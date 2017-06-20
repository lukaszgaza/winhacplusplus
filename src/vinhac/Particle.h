#ifndef __Particle_h__
#define __Particle_h__
#include <string>
#include <vector>
#include <exception>
#include "../utils/CLHEP/Vector/Vector/LorentzVector.h"

using namespace std;

namespace VINHAC {
class Particle;
}

namespace VINHAC {

/**
 * \brief class representing particle
 *
 * This class represents particle, it contains needed properties of particle used in generation.
 */
class Particle {

private:
	vector<VINHAC::Particle*> _parents;
	vector<VINHAC::Particle*> _children;

	int _PDGid; //!< PDG of the particle
	CLHEP::HepLorentzVector p; //!< Four momentum of a particle
	CLHEP::HepLorentzVector p_brf2;
	string m_tag; //!< Tags initial/final particles

public:

	/**
	 * \brief default contructor
	 */
	Particle();

	/**
	 * \brief copy constructor
	 *
	 * @param target Particle to copy
	 */
	Particle(const Particle& target);

	/**
	 * \brief contructor
	 *
	 * @param PDGid PDGid to set
	 */
	Particle(int PDGid);

	/**
	 * \brief contructor
	 *
	 * @param PDGid PDGid to set
	 * @param momentum fourmomentum to set
	 */
	Particle(int PDGid, CLHEP::HepLorentzVector momentum);

	/**
	 * \brief contructor
	 *
	 * @param PDGid PDGid to set
	 * @param parents parents of the particle
	 */
	Particle(int PDGid, vector<VINHAC::Particle*> parents);

	/**
		 * \brief contructor
		 *
		 * @param PDGid PDGid to set
		 * @param parents parents of the particle
		 * @param momentum fourmomentum to set
		 */
	Particle(int PDGid, vector<VINHAC::Particle*> parents,
			CLHEP::HepLorentzVector momentum);

	//setters

	/**
	 * \brief sets PDGid
	 *
	 * @param PDGcode PDGid to set
	 */
	inline void setPDGid(int PDGcode) {
		this->_PDGid = PDGcode;
	}

	/**
	 * \brief adds child to the particle
	 *
	 * It also adds current particle to child's parents
	 * @param newChild child to add
	 */
	inline void addChild(Particle* newChild) {
		_children.push_back(newChild);
		newChild->_parents.push_back(this);
	}

	/**
	 * \brief adds parent for this particle
	 *
	 * It also adds current particle to parent's childrens
	 * @param newParent parent to add
	 */
	inline void addParent(Particle* newParent) {
		_parents.push_back(newParent);
		newParent->_children.push_back(this);
	}


	/**
	 * \brief sets four-momentum for this particle
	 *
	 * @param momentum four-momentum to set
	 */
	inline void setFourMomentum(CLHEP::HepLorentzVector momentum) {
		p = momentum;
	}

	/**
	 * \brief sets four-momentum for this particle in boson rest frame
	 *
	 * @param momentum four-momentum to set
	 */
	inline void setFourMomentumBosonRestFrame(CLHEP::HepLorentzVector momentum) {

		this->p_brf2 = momentum;

	}
	// getters

	/**
	 * \brief returns PDGid of particle
	 *
	 * @return PDGid
	 */
	inline int getPDGid() {
		return this->_PDGid;
	}

	/**
	 * \brief returns four-momentum of particle
	 *
	 * @return four-momentum
	 */
	inline CLHEP::HepLorentzVector getFourMomentum() {
		return this->p;
	}

	/**
	 * \brief returns four-momentum of particle in boson rest frame
	 *
	 * @return four-momentum
	 */
	inline CLHEP::HepLorentzVector getFourMomentumBosonRestFrame() {
		return this->p_brf2;
	}

	/**
	 * \brief returns vector of parents
	 *
	 * @return vector of Particle objects
	 */
	inline vector<VINHAC::Particle*> getParents() {
		return _parents;
	}

	/**
	 * \brief returns vector of childrens
	 *
	 * @return vector of Particle objects
	 */
	inline vector<VINHAC::Particle*> getChildren() {
		return _children;
	}

	/**
	 * \brief removes all children
	 */
	inline void clearChildren() {
		_children.clear();
	}

	/**
	 * \brief prints particle
	 *
	 * In prints PDGid, four-momentum and four-momentum in boson rest frame
	 */
	inline void print() {
		cout << "PDGid: " << _PDGid << " " << p << " " << p_brf2 << endl;
	}

	/**
	 * \brief sets tag
	 *
	 * Tag is a string which is used to recognize particle by handlers
	 * @param newTag tag to set
	 */
	inline void setTag(string newTag) {
		m_tag = newTag;
	}

	/**
	 * \brief returns tag
	 *
	 * Tag is a string which is used to recognize particle by handlers
	 * @return tag
	 */
	inline string getTag() {
		return m_tag;
	}

	/**
	 * \brief clears particle data
	 */
	void clear();
};

}

#endif
