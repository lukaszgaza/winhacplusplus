#ifndef __PDF_h__
#define __PDF_h__

#include <string>
#include <vector>
#include <exception>
using namespace std;

namespace VINHAC {
class PDF;
}

namespace VINHAC {

/**
 * \brief base class for classes handling Parton Distribution Function
 */
class PDF {

private:
	string m_name; // PDF name
public:
	//*********
	/*! \brief Nucleon PDF: returns a vector \f$ x f_i(x, Q) \f$ with index
	 * \f$ 0 < i < 12 \f$.
	 * \f$ 0..5 = \bar{t} , ..., \bar{u}, \bar{d};
	 * 6 = g;
	 * 7..12 = d, u, ..., t \f$.
	 *********** \author Andrzej Siodmok \date Created on: 2009-05-23 ******/

	virtual vector<double> xfx(double x, double q) = 0;
	//*********
	/*! \brief Nucleon PDF: returns \f$ x f(x, Q) \f$ for flavour fl -
	 *   this time the flavour encoding is as in the LHAPDF manual.
	 *  \f$-6..-1 = \bar{t}, ..., \bar{u}, \bar{d};
	 *   0 = g
	 *   1..6 = d, u, ..., t \f$.
	 *********** \author Andrzej Siodmok \date Created on: 2009-05-23 ******/
	virtual double xfx(double x, double q, int fl) = 0;
	/*
	 public: virtual void initPDF() = 0;
	 */

	/**
	 * @param m set number
	 * @return minimal x
	 */
	virtual double getXmin(int m) = 0;

	/**
	 * @param m set number
	 * @return maximal x
	 */
	virtual double getXmax(int m) = 0;

	/**
	 * @param m set number
	 * @return minimal \f$Q^2\f$
	 */
	virtual double getQ2min(int m) = 0;

	/**
	 * @param m set number
	 * @return maximal \f$Q^2\f$
	 */
	virtual double getQ2max(int m) = 0;

	/**
	 * sets name of the pdf
	 *
	 * @param name to set
	 */
	inline void setName(const string& name) {
		m_name = name;
	}

	/**
	 * returns name of pdf
	 *
	 * @return name
	 */
	inline const string& getName() const {
		return m_name;
	}
	virtual ~PDF() {
	}
	;
};
}

#endif
