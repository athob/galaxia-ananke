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

#include "SurveyDesign.h"
#include "Functions.h"

void stringstreamTovchar(stringstream &sout,vector<char> &cv)
{
	int64_t pos = sout.tellg();
	sout.seekg (0, ios::end);
	int length = sout.tellg();
	cv.resize(length+1);
	sout.seekg (0, ios::beg);
	sout.read (&cv[0],length);
	cv[length]=0;
	sout.seekg (pos, ios::beg);
}



SurveyDesign::~SurveyDesign()
{
	if(geo!=NULL)
		delete geo;
// TODO Auto-generated destructor stub
}

void SurveyDesign::setCenter(double *pos)
{
	for (int i=0;i<6;++i)
		posC[i]=pos[i];
}

void SurveyDesign::writeStars(vector<StarParticle> &Stars1)
///changed to export doubles rather than floats
{
	for (size_t i = 0; i < Stars1.size(); ++i)
	{


		ebfA[0].Write(&(Stars1[i].pos(0)),1);
		ebfA[1].Write(&(Stars1[i].pos(1)),1);
		ebfA[2].Write(&(Stars1[i].pos(2)),1);
		ebfA[3].Write(&(Stars1[i].pos(3)),1);
		ebfA[4].Write(&(Stars1[i].pos(4)),1);
		ebfA[5].Write(&(Stars1[i].pos(5)),1);
		ebfA[6].Write(&(Stars1[i].feh()),1);
		ebfA[7].Write(&(Stars1[i].alpha()),1);
		ebfA[8].Write(&(Stars1[i].smass()),1);
		ebfA[9].Write(&(Stars1[i].age()),1);
		ebfA[10].Write(&(Stars1[i].rad()),1);
		ebfA[11].Write(&(Stars1[i].mag(0)),1);
		ebfA[12].Write(&(Stars1[i].mag(1)),1);
		ebfA[13].Write(&(Stars1[i].mag(2)),1);
		ebfA[14].Write(&(Stars1[i].popID()),1);
		ebfA[15].Write(&(Stars1[i].satID()),1);
		ebfA[16].Write(&fieldID,1);
		ebfA[17].Write(&(Stars1[i].partID()),1);
		//ID of parent star particle
		ebfA[18].Write(&(Stars1[i].parentID()),1);
		//ID of parent particle partition
		ebfA[19].Write(&(Stars1[i].partitionID()),1);
		//formation distance
		ebfA[20].Write(&(Stars1[i].dform()),1);
		//extra abundances
/*		if(i==0){
			cout<<"extra abundances: "<<Stars1[i].helium()<<" "<<Stars1[i].calcium()<<endl;
		}*/
		ebfA[21].Write(&(Stars1[i].helium()),1);
		ebfA[22].Write(&(Stars1[i].carbon()),1);
		ebfA[23].Write(&(Stars1[i].nitrogen()),1);
		ebfA[24].Write(&(Stars1[i].oxygen()),1);
		ebfA[25].Write(&(Stars1[i].neon()),1);
		ebfA[26].Write(&(Stars1[i].magnesium()),1);
		ebfA[27].Write(&(Stars1[i].silicon()),1);
		ebfA[28].Write(&(Stars1[i].sulphur()),1);
		ebfA[29].Write(&(Stars1[i].calcium()),1);
		//export distance modulus
		ebfA[30].Write(&(Stars1[i].dmod()),1);
	}
}


int SurveyDesign::push_back(StarParticle &Star)
{

	if (checkColMag(Star))
	{
		if (geo->checkP(&Star.pos(0)))
		{
			//Star.convert(Starf);
			Stars.push_back(Star);
			nstars++;
			if(nstars==numeric_limits<uint64_t>::max())
			{
				cout<<"ERROR: number of stars larger than maximum int limit- options are"<<endl;
				cout<<"1) Reduce magnitude limits"<<endl;
				cout<<"2) Reduce fSample"<<endl;
				cout<<"3) Make changes in code in file 'SurveyDesign.cpp' to handle larger integers"<<endl;
				exit(1);
			}
			if (Stars.size() == 100000)
			{
				writeStars(Stars);
				Stars.clear();
			}
			return 1;
		}
		else
			return 0;
	}
	return 0;
}

int SurveyDesign::force_push_back(StarParticle &Star)
{
	//Star.convert(Starf);
	Stars.push_back(Star);
	nstars++;
	if(nstars==numeric_limits<uint64_t>::max())
	{
		cout<<"ERROR: number of stars larger than maximum int limit- options are"<<endl;
		cout<<"1) Reduce magnitude limits"<<endl;
		cout<<"2) Reduce fSample"<<endl;
		cout<<"3) Make changes in code in file 'SurveyDesign.cpp' to handle larger integers"<<endl;
		exit(1);
	}
	if (Stars.size() == 100000)
	{
		writeStars(Stars);
		Stars.clear();
	}
	return 1;
}

int SurveyDesign::push_check(StarParticle &Star)
{
	if (checkColMag(Star))
	{
		if (geo->checkP(&Star.pos(0)))
			return 1;
		else
			return 0;
	}
	return 0;
}


void SurveyDesign::flush( )
{
	if(Stars.size()>0)
	{
		writeStars(Stars);
		Stars.clear();
	}
	npart.push_back(nstars);
	nstars=0;

}

void SurveyDesign::close()
{
	certify(Stars.size()==0);
	nstars=total(npart);
	if(nstars==numeric_limits<uint64_t>::max())
	{
		cout<<"ERROR: number of stars larger than maximum int limit- options are"<<endl;
		cout<<"1) Reduce magnitude limits"<<endl;
		cout<<"2) Reduce fSample"<<endl;
		cout<<"3) Make changes in code in file 'SurveyDesign.cpp' to handle larger integers"<<endl;
		exit(1);
	}

	cout<<left<<setw(36)<<"Total stars written "<<setw(24)<<nstars<<endl;
	{
		ebf::EbfFile ebf3;


		stringstream sout;
		sout<<"# File generated by "<<version()<<endl;
		sout<<"# <parameterfile>"<<endl;
		vector<char> buf;
		sout<<All->outParameterFile();
		sout<<"# </parameterfile>"<<endl;
		stringstreamTovchar(sout,buf);
		ebf::Write(outputFileR,"/Log",&buf[0],"w","",buf.size());

		cout << ebfA.size() << endl;

		for(size_t i=0;i<ebfA.size();++i)
		{
			ebfA[i].SaveTo(outputFileR,"a");
		}

		ebf::Write(outputFileR,"/Center",&posC[0],"a","1000 parsec",6);

		if(All->fieldTable.col.size() != 0)
		{
			ebf::Write(outputFileR,"/Field/Longitude",&(All->fieldTable.col[0][0]),"a","degree",All->fieldTable.col[0].size());
			ebf::Write(outputFileR,"/Field/Latitude",&(All->fieldTable.col[1][0]),"a","degree",All->fieldTable.col[1].size());
		}


	}

	cout<<left<<setw(36)<<"File written- "<<outputFileR<<endl;
}

void SurveyDesign::setError(int option,double sigma_r,double sigma_vr,double sigma_mu,double sigma_fe,double sigma_al)
{
//	stub
}

void SurveyDesign::setGeometry(int option,double l,double b,double area, int fieldNo)
{
	if(geo!=NULL)
		delete geo;
	fieldID=fieldNo;

	double temp,th,dth;
	double rad2deg=180.0/PI;
	double deg2rad=PI/180.0;
	switch (option)
	{
	case 0:
		geo=new AllSky();
		cout<<left<<setw(36)<<"Using geometry:"<<"All Sky"<<endl;
		break;
	case 1:
		temp=area*(deg2rad*deg2rad)/(2*PI);
		assert((temp>=0.0)&&(temp<=2));
		dth=rad2deg*min(acos(1-temp),PI);
		cout<<left<<setw(36)<<"Using geometry:"<<"Patch at l ,b : ("<<l<<" "<<b<<") d_theta="<<dth<<endl;
		geo=new GPatch(l,b,dth);
		break;
	case 2:
		temp=60.0;
		th=rad2deg*max(acos(cos(temp*deg2rad)+area*(deg2rad*deg2rad)/(2*PI)),0.0);
		dth=temp-th;
		cout<<left<<setw(36)<<"Using geometry:"<<"Angle Cone Strip th dth: "<<th<<" "<<th+dth<<endl;
		geo=new Cone(th,dth);
		break;
	case 3:
		temp=rad2deg*min(acos(area*(3.0/2.0)*(deg2rad*deg2rad)/(2*PI)),PI/2);
		th=temp;
		dth=2*(90.0-th);
		cout<<left<<setw(36)<<"Using geometry:"<<"Angle Wedge X th dth: "<<th<<" "<<th+dth<<endl;
		geo=new Wedge(th,dth,0);
		break;
	case 4:
		temp=rad2deg*min(acos(area*(3.0/2.0)*(deg2rad*deg2rad)/(2*PI)),PI/2);
		th=temp;
		dth=2*(90.0-th);
		cout<<left<<setw(36)<<"Using geometry:"<<"Angle Wedge Y th dth: "<<th<<" "<<th+dth<<endl;
		geo=new Wedge(th,dth,1);
		break;
	case 5:
		th=0.0;
		dth=rad2deg*min(acos(1-area*(deg2rad*deg2rad)/(2*PI)),PI);
		cout<<left<<setw(36)<<"Using geometry:"<<"Angle Cone Patch th dth:"<<th<<" "<<th+dth<<endl;
		geo=new Cone(th,dth);
		break;
	default:
		break;
	}
}

bool SurveyDesign::checkColMag(StarParticle & Star)
{
	if(All->starType==1)
		if((Star.teff()>(4.00279 -0.079*Star.lum()))||(Star.teff()<(4.00279 -0.079*Star.lum()-0.09))||(Star.lum()<1.4)||(Star.lum()>2.0))
			return 0;

	if(All->starType==2)
		if((Star.lum()<1.4)||(Star.lum()>2.0)||(Star.teff()<3.8)||(Star.smass()<0.5)||(Star.smass()>1.0))
			return 0;


	if((Star.mag(0) >= appMag[0])&&(Star.mag(0) <= appMag[1])&&((Star.mag(1)-Star.mag(2)) >= color[0])&&((Star.mag(1)-Star.mag(2)) <= color[1])&&((Star.mag(0)-Star.dmod()) >= absMag[0])&&((Star.mag(0)-Star.dmod()) <= absMag[1])&&(Star.rad() <=r_max)&&(Star.rad()>r_min))
		return 1;
	else
		return 0;


}









