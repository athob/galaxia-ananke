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



#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "TableInterpolator.h"

class ParMem
{
public:
	ParMem(const char *s,void *addr1,int id1,const char *buf,int readOn1):id(id1),status(0),addr(addr1),readOn(readOn1)
	{
		tag=s;
		read(buf);
	}
	string tag;
	int id;
	int status;
	void *addr;
	int readOn;
	void read(const char *buf2)
	{
		switch (id)
		{
		case 1:
			*((double*) addr) = atof(buf2);
			//			fprintf(fdout, "%-35s%g\n", buf1, *((double*) addr[j]));
			break;
		case 2:
			*((string *) addr) = buf2;
			//				strcpy((char *) addr, buf2);
			//			fprintf(fdout, "%-35s%s\n", buf1, buf2);
			break;
		case 3:
			*((int*) addr) = atoi(buf2);
			//			fprintf(fdout, "%-35s%d\n", buf1, *((int*) addr[j]));
			break;
		case 4:
			*((uint64_t*) addr) = strtoull(buf2, NULL, 10);
			//			fprintf(fdout, "%-35s%d\n", buf1, *((int*) addr[j]));
			break;
		}
	}

	void write(stringstream &out)
	{
		switch (id)
		{
		case 1:
			if(readOn==1)
				out<<setw(35)<<left<<tag.c_str()<<" "<<(*((double*) addr))<<endl;
			else
				out<<"# "<<setw(35)<<left<<tag.c_str()<<" "<<(*((double*) addr))<<endl;
			//				fprintf(fdout, "%-35s%g\n", tag.c_str(), *((double*) addr));
			break;
		case 2:
			//				fprintf(fdout, "%-35s%s\n", tag, (char *) addr);
			if(readOn==1)
				out<<setw(35)<<left<<tag.c_str()<<" "<<(*(string *) addr)<<endl;
			else
				out<<"# "<<setw(35)<<left<<tag.c_str()<<" "<<(*(string *) addr)<<endl;
			break;
		case 3:
			if(readOn==1)
				out<<setw(35)<<left<<tag.c_str()<<" "<<(*((int*) addr))<<endl;
			else
				out<<"# "<<setw(35)<<left<<tag.c_str()<<" "<<(*((int*) addr))<<endl;
			break;
		case 4:
			if(readOn==1)
				out<<setw(35)<<left<<tag.c_str()<<" "<<(*((uint64_t*) addr))<<endl;
			else
				out<<"# "<<setw(35)<<left<<tag.c_str()<<" "<<(*((uint64_t*) addr))<<endl;
			break;

		}
	}

};




class Parameters
{
public:
	Parameters();
	~Parameters();
	void setFromParameterFile(const string fname);
	void saveParameterFile(const string fnameToSave);
	string outParameterFile();
	void setFromArguments(int argc,char **argv);
	void checkCompilation( );
	void print( );
	void usage( );
	void load_sat_list();
	void copyright();
	//---General  options--------------------------------------------
	int Debug;
	string outputFile,inputDir,SuSuffix,parameterFile,outputDir,fieldTableFile,nbodyFile,addstring;
	//---Survey  options--------------------------------------------
	int geometryOption;
	double surveyArea,fSample;
	double absMagLimits[2];
	double appMagLimits[2];
	double colorLimits[2];
	double posC[6];
	double r_min,r_max;//,age_min,age_max;
	string photoSys,magcolorNames;

	int starType,photoError,hdim,nres;
	uint64_t popID, nstart;
	int warpFlareOn,option,seed;
	double sigma_r,sigma_vr,sigma_mu,sigma_fe,sigma_al;
	double latitude,longitude;
//	int sunCentricOn;
//	int colorMagCutOption;
	Table fieldTable;
	std::vector<ParMem>  par_list;
	std::vector<std::pair<string,int> > sat_list;
};

string version();

#endif /*PARAMETERS_H_*/
