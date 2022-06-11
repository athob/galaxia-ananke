/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Author: Sanjib Sharma                                                     *
 * Copyright (c) 2012 Sanjib Sharma                                          *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of Galaxia. The full Galaxia copyright notice, including*
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and COPYRIGHT, which can be found at the root           *
 * of the source code distribution tree.                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SURVEYERROR_H_
#define SURVEYERROR_H_
#include "StarParticle.h"
#include "Matrix.h"

class SurveyError
{
public:
	SurveyError(int errorOption1):errorOption(errorOption1),gauss(0.0,1.0,33){initialize();}
	vector<double> appMag;
	vector<double> vr;
	double sigma_w,sigma_mu,sigma_r,sigma_vr,sigma_vlb,sigma_fe,sigma_al;
	int errorOption;
	Normaldev gauss;
	void initialize();
	void add(StarParticle &Star);
	~SurveyError();
};




#endif /* SURVEYERROR_H_ */
