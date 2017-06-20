/*
 * Functions.h
 *
 *  Created on: 2010-04-14
 *      Author: kamil
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <string>
#include <vector>

void Tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters);

#endif /* FUNCTIONS_H_ */
