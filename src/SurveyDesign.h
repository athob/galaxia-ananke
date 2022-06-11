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

#ifndef SURVEYDESIGN_H_
#define SURVEYDESIGN_H_
#include "Geometry.h"
#include "Parameters.h"
#include "SurveyError.h"
#include "ebfvector.hpp"


class SurveyDesign
{
public:

	SurveyDesign(const string fname, double appMagLimits[2],double absMagLimits[2],double colorLimits[2],Parameters* All1)
	{
		All=All1;
//		sunCentricOn=sunCentricOn1;
//		colorMagCutOption=colorMagCutOption1;
		appMag[0]=appMagLimits[0];
		appMag[1]=appMagLimits[1];
		absMag[0]=absMagLimits[0];
		absMag[1]=absMagLimits[1];
		color[0]=colorLimits[0];
		color[1]=colorLimits[1];
		r_max=All->r_max;
		r_min=All->r_min;
		nstart=All->nstart;
//		age_min=All->age_min;
//		age_max=All->age_max;
		// if(sunCentricOn==1)
		// {
		// 	posC[0]=-8.5; posC[1]=0.0; posC[2]=0.0;
		// 	posC[3]=0.0; posC[4]=0.0; posC[5]=0.0;
		// 	cout<<"Sun centering set"<<endl;
		// }
		// else
		// {
		// 	posC[0]=0.0; posC[1]=0.0; posC[2]=0.0;
		// 	posC[3]=0.0; posC[4]=0.0; posC[5]=0.0;
		// }
		outputFileR=fname;
		outputFile=outputFileR;		outputFile+=".tmp";

		printFlag.clear();
		printFlag.resize(18,1);

//changed to output int64 and doubles rather than int32 and floats

		ebfA.clear();
		ebfA.resize(30);
		ebfA[0].Open(outputFileR+".tmp0","/px","w",5,"1000 parsec");
		ebfA[1].Open(outputFileR+".tmp1","/py","w",5,"1000 parsec");
		ebfA[2].Open(outputFileR+".tmp2","/pz","w",5,"1000 parsec");
		ebfA[3].Open(outputFileR+".tmp3","/vx","w",5,"1000 meter/second");
		ebfA[4].Open(outputFileR+".tmp4","/vy","w",5,"1000 meter/second");
		ebfA[5].Open(outputFileR+".tmp5","/vz","w",5,"1000 meter/second");
		ebfA[6].Open(outputFileR+".tmp6","/feh","w",5);
		ebfA[7].Open(outputFileR+".tmp7","/alpha","w",5);
		ebfA[8].Open(outputFileR+".tmp8","/smass","w",5,"solar_mass");
		ebfA[9].Open(outputFileR+".tmp9","/age","w",5,"log year");
		ebfA[10].Open(outputFileR+".tmp10","/rad","w",5,"1000 parsec");
		ebfA[11].Open(outputFileR+".tmp11","/mag0","w",5);
		ebfA[12].Open(outputFileR+".tmp12","/mag1","w",5);
		ebfA[13].Open(outputFileR+".tmp13","/mag2","w",5);
		ebfA[14].Open(outputFileR+".tmp14","/popid","w",13);
		ebfA[15].Open(outputFileR+".tmp15","/satid","w",3);
		ebfA[16].Open(outputFileR+".tmp16","/fieldid","w",3);
		ebfA[17].Open(outputFileR+".tmp17","/partid","w",3);
		ebfA[18].Open(outputFileR+".tmp18","/parentid","w",3);

//added formation distance
		ebfA[19].Open(outputFileR+".tmp19","/dform","w",5,"1000 parsec");

//added extra abundances		
		ebfA[20].Open(outputFileR+".tmp20","/helium","w",5);
		ebfA[21].Open(outputFileR+".tmp21","/carbon","w",5);
		ebfA[22].Open(outputFileR+".tmp22","/nitrogen","w",5);
		ebfA[23].Open(outputFileR+".tmp23","/oxygen","w",5);
		ebfA[24].Open(outputFileR+".tmp24","/neon","w",5);
		ebfA[25].Open(outputFileR+".tmp25","/magnesium","w",5);
		ebfA[26].Open(outputFileR+".tmp26","/silicon","w",5);
		ebfA[27].Open(outputFileR+".tmp27","/sulphur","w",5);
		ebfA[28].Open(outputFileR+".tmp28","/calcium","w",5);
//added distance modulus
		ebfA[29].Open(outputFileR+".tmp29","/dmod","w",5);


		if(All->fieldTable.col.size() == 0)
			printFlag[12]=0;

		outputFile=outputFileR;		outputFile+=".tmp";

		Stars.reserve(100000);
		nstars=0;
		fieldID=0;
		geo=NULL;
	}

	~SurveyDesign();
	
	Geometry* geo;
	int force_push_back(StarParticle &Star);
	int push_back(StarParticle &Star);
	int push_check(StarParticle &Star);
	void flush();
	void close();
	bool checkColMag(StarParticle &Star);

	void writeStars(vector<StarParticle> &Stars1);
	string outputFileR;
	string outputFile;
	vector<ebf::EbfFile> ebfA;

	vector<int> printFlag;
	double appMag[2];
	double absMag[2];
	double color[2];
	double r_min,r_max;//,age_min,age_max;
	Parameters* All;
//	int sunCentricOn;//colorMagCutOption;
	double posC[6];
//	double velC[3];
	vector<StarParticle> Stars;
	uint64_t nstars, nstart;
	vector<uint64_t> npart,npartc;
	int fieldID;
	void setGeometry(int option,double l,double b,double area,int fieldNo=0);
	void setError(int option,double sigma_r,double sigma_vr,double sigma_mu,double sigma_fe,double sigma_al);
	void setCenter(double *pos);
};



#endif /* SURVEYDESIGN_H_ */
