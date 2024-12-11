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

#include"Parameters.h"
#include "Timer.h"
#include "ebfvector.hpp"
#include "Satellite.h"
#include"IsochroneDB.h"
#include"SurveyDesign.h"
#include "Population.h"
#include "ebf.hpp"

double Chabrier_exp_imf1(double x)
{
	return 22.8978*(exp(-pow((716.4/x),0.25)))*pow(x,-3.3);
}

double Kroupa_imf(double m)
{
	if(m >0.5)
		return 0.04*pow(m,-2.3);
	else if (m <0.08)
		return pow(m,-0.3);
	else
		return 0.08*pow(m,-1.3);
}


void append_ext(const string& filename,const string dir)
{
//	double PI=3.14159;
//	double DTOR=0.0174533;
	ebf::EbfDataInfo dinfo;
	double RADEG=57.2958;
	GInterpolator gi;
	gi.readFromFile(dir+"Extinction/ExMap3d_1024.ebf","/ExMap3d");
	GInterpolator ge_solar;
	ge_solar.readFromFile(dir+"Extinction/ExMap3d_1024.ebf","/ExMap2d");
	GInterpolator ge_schlegel;
	ge_schlegel.readFromFile(dir+"Extinction/Schlegel_4096.ebf","/ExMap2d");
	double PosC[3],PosS[3];
//	fbvector<float> x(filename+".ebf","r","/Pos3");
	ebf::EbfVector<float> x(filename+".ebf","/px");
	ebf::EbfVector<float> y(filename+".ebf","/py");
	ebf::EbfVector<float> z(filename+".ebf","/pz");
	vector<float> exbv1(x.dim(0),0.0);
	vector<float> exbv2(x.dim(0),0.0);
	vector<float> exbv3(x.dim(0),0.0);
	vector<float> exbv4(x.dim(0),0.0);
	vector<float> exbv5(x.dim(0),0.0);
//	certify(x.dim(1)==3,"Check map dim");
	for(int64_t i=0;i<x.dim(0);++i)
	{
//		cout<<i<<" "<<x.dim(0)<<" "<<x.size()<<" "<<i*x.dim(0)+2<<endl;
//		PosC[0]=x[i*x.dim(1)+0];	PosC[1]=x[i*x.dim(1)+1];	PosC[2]=x[i*x.dim(1)+2];
		PosC[0]=x[i];	PosC[1]=y[i];	PosC[2]=z[i];
		Cot::xyz_to_lbr(PosC,PosS,3);
		PosS[0]*=RADEG;		PosS[1]*=RADEG;
		if(PosS[0]<0)
			PosS[0]=360.0+PosS[0];
		PosS[2]=log10(PosS[2]);
		exbv1[i]=gi.interpol(PosS);
		exbv2[i]=exbv1[i];
		exbv1[i]*=ge_schlegel.interpol(PosS);
		exbv2[i]*=ge_solar.interpol(PosS);
		exbv3[i]=ge_schlegel.interpol(PosS);
		exbv4[i]=PosS[0];
		exbv5[i]=PosS[1];
	}

	if(x.dim(0) >0)
	{
		std::string filename1=filename+".ebf";
		if(ebf::ContainsKey(filename1,"/ExBV_Schlegel",dinfo)==0)
		{
			ebf::Write(filename1,"/ExBV_Schlegel",&exbv1[0],"a","",exbv1.size());
			ebf::Write(filename1,"/ExBV_Solar",&exbv2[0],"a","",exbv2.size());
			ebf::Write(filename1,"/ExBV_Schlegel_Inf",&exbv3[0],"a","",exbv3.size());
		}
		else
		{
			cout<<"Skipping Extinction as data already exists"<<endl;
		}

		if(ebf::ContainsKey(filename1,"/glon",dinfo)==0) {
			ebf::Write(filename1,"/glon",&exbv4[0],"a","degree",exbv4.size());
			ebf::Write(filename1,"/glat",&exbv5[0],"a","degree",exbv5.size());
		}
		else
		{
			cout<<"Skipping glat/glon as data already exists"<<endl;
		}
	}

}

void append_lb(const string& filename,const string dir)
{
//calculates and appends galactic lat and longitude to data file
	ebf::EbfDataInfo dinfo;

	std::string filename1=filename+".ebf";

	//check if file already contains this info
	if(ebf::ContainsKey(filename1,"/glat",dinfo)) {
		cout<<"File already contains galactic coordinates"<<endl;
		return;
	}

	//if not, load one coordinate to check size
	ebf::EbfVector<float> x(filename+".ebf","/px");

	if(x.dim(0) >0) {  //if file is not empty

		double RADEG=57.2958;
		double PosC[3],PosS[3];
		ebf::EbfVector<float> y(filename+".ebf","/py"); //read other coordinates
		ebf::EbfVector<float> z(filename+".ebf","/pz");
		vector<float> glon(x.dim(0),0.0);
		vector<float> glat(x.dim(0),0.0);

		for(int64_t i=0;i<x.dim(0);++i)
		{
			PosC[0]=x[i];	PosC[1]=y[i];	PosC[2]=z[i];
			Cot::xyz_to_lbr(PosC,PosS,3); //convert to l,b in radians

			PosS[0]*=RADEG;		PosS[1]*=RADEG;  //convert rads to degs

			if(PosS[0]<0) //re-wrap longitude from 0 to 360 degrees
				PosS[0]=360.0+PosS[0];

			glon[i]=PosS[0];
			glat[i]=PosS[1];
		}
		ebf::Write(filename1,"/glon",&glon[0],"a","degree",glon.size());
		ebf::Write(filename1,"/glat",&glat[0],"a","degree",glat.size());

	} else {
		cout<<"File contains no data"<<endl;
	}
}

void append_radec(const string& filename,const string dir)
{
//	double PI=3.14159;
//	double DTOR=0.0174533;
	ebf::EbfDataInfo dinfo;
	double RADEG=57.2958;
	double PosC[3],PosS[3],alpha,delta;
//	fbvector<float> x(filename+".ebf","r","/Pos3");
	ebf::EbfVector<float> x(filename+".ebf","/px");
	ebf::EbfVector<float> y(filename+".ebf","/py");
	ebf::EbfVector<float> z(filename+".ebf","/pz");
	vector<float> exbv1(x.dim(0),0.0);
	vector<float> exbv2(x.dim(0),0.0);
//	certify(x.dim(1)==3,"Check map dim");
	for(int64_t i=0;i<x.dim(0);++i)
	{
//		cout<<i<<" "<<x.dim(0)<<" "<<x.size()<<" "<<i*x.dim(0)+2<<endl;
//		PosC[0]=x[i*x.dim(1)+0];	PosC[1]=x[i*x.dim(1)+1];	PosC[2]=x[i*x.dim(1)+2];
		PosC[0]=x[i];	PosC[1]=y[i];	PosC[2]=z[i];
		Cot::xyz_to_lbr(PosC,PosS,3);
		PosS[0]*=RADEG;		PosS[1]*=RADEG;
		if(PosS[0]<0)
			PosS[0]=360.0+PosS[0];
		Cot::lb_to_radec(PosS[0],PosS[1],alpha,delta);
		exbv1[i]=alpha;
		exbv2[i]=delta;
	}

	if(x.dim(0) >0)
	{
		std::string filename1=filename+".ebf";
		if(ebf::ContainsKey(filename1,"/ra",dinfo)==0)
		{
			ebf::Write(filename1,"/ra",&exbv1[0],"a","degree",exbv1.size());
			ebf::Write(filename1,"/dec",&exbv2[0],"a","degree",exbv1.size());
		}
		else
		{
			cout<<"Skipping ra,dec addition as data already exists"<<endl;
		}
	}
}


void append1(const string& filename,IsochroneDB& ic,const string& photodir,const string& photosys,const string&magcolorNames)
{
		ebf::EbfDataInfo dinfo;	
		IsoFileDescriptor isofile_info(photodir+"/"+photosys+"/IsoFileDescriptor.txt",photosys,magcolorNames);
		vector<double> x(isofile_info.magnames.size()+5,0.0);
		ebf::EbfVector<float>   age(filename,"/Age");
		ebf::EbfVector<float>   feh(filename,"/FeH");
		ebf::EbfVector<float> smass(filename,"/Minit");  // used to be stored under "/smass"

	//	for(size_t i=0;i<isofile_info.magnames.size();++i)
	//		cout<<i<<" "<<isofile_info.magnames[i]<<endl;

		//check that there are stars in the file
		if(age.size()>0) {


		// initialize fbvectors
		cout<<"initializing isochrone data"<<endl;

		ebf::EbfFile lum,teff,grav,mtip,mact;
		lum.Open(filename+".Atmp0","/lum","w",4,"log solar_luminosity");
		teff.Open(filename+".Atmp1","/teff","w",4,"log kelvin");
		grav.Open(filename+".Atmp2","/grav","w",4,"log 0.01 meter/second/second");
		mact.Open(filename+".Atmp3","/mact","w",4,"solar_mass");
		mtip.Open(filename+".Atmp4","/mtip","w",4,"solar_mass");

		vector<ebf::EbfFile>  mags;
		mags.resize(isofile_info.magnames.size());
		for(size_t k=0;k<isofile_info.magnames.size();++k)
		{
			stringstream sout;
			sout<<filename<<".Atmp"<<k+5;
			mags[k].Open(sout.str(),"/"+photosys+"_"+isofile_info.magnames[k],"w",4,"magnitude");
		}
		// interpolate

		cout<<"interpolating on isochrone tables"<<endl;
		for(size_t i=0;i<age.size();++i)
		{
			ic.interpolateTGM(double(age[i]),double(feh[i]),double(smass[i]),x);
			lum.Write(&x[0],1);
			teff.Write(&x[1],1);
			grav.Write(&x[2],1);
			mact.Write(&x[3],1);
			mtip.Write(&x[4],1);
			for(size_t k=0;k<isofile_info.magnames.size();++k)
			{
				mags[k].Write(&x[k+5],1);
			}
		}
		if(ebf::ContainsKey(filename,"/teff",dinfo)==0)
		{
			lum.SaveTo(filename,"a");
			teff.SaveTo(filename,"a");
			grav.SaveTo(filename,"a");
		}
		else
		{
			cout<<"Skipping"<<" "<<"/teff"<<" "<<"/lum"<<" "<<"/grav"<<endl;
			lum.Remove();
			teff.Remove();
			grav.Remove();
		}
		if(ebf::ContainsKey(filename,"/mact",dinfo)==0)
		{
			mact.SaveTo(filename,"a");
			mtip.SaveTo(filename,"a");
		}
		else
		{
			cout<<"Skipping"<<" "<<"/mact"<<" "<<"/mtip"<<endl;
			mact.Remove();
			mtip.Remove();
		}

		if(ebf::ContainsKey(filename,mags[0].getDataName(),dinfo)==0)
		{
			for(size_t k=0;k<isofile_info.magnames.size();++k)
				mags[k].SaveTo(filename,"a");
		}
		else
		{
			cout<<"Skipping"<<" "<<mags[0].getDataName()<<" ---- "<<mags[isofile_info.magnames.size()-1].getDataName()<<endl;
			for(size_t k=0;k<isofile_info.magnames.size();++k)
				mags[k].Remove();
		}
	} else {
		cout<<"Error: survey contains no stars"<<endl;
	}


}




void runModel(SurveyDesign &sur, IsochroneDB& ic, Parameters &All)
{		
	Sampler imf(10000, 0.01, 120.0, Kroupa_imf, 1);

	for (unsigned int i = 0; i < All.sat_list.size(); ++i)
	{
		cout << "------------------------------" << endl;
		cout << All.sat_list[i].first << "  Sat No="
				<< All.sat_list[i].second << endl;

		Satellite Sat(All.inputDir+All.sat_list[i].first, 0, All.sat_list[i].second, All.hdim, All.nres);
		Sat.spawn1(sur, imf, ic, All.fSample, All.seed + i + 100);
		cout << "-----------Done---------------" << endl;
	}
	sur.flush();


}

int main(int argc, char **argv)
{
	Timer timer1, timer2;
	srand48(13);

	Parameters All;
	All.setFromArguments(argc, argv);
	nrRan = Ran(All.seed+4);
	nrGauss = Normaldev(0.0, 1.0, All.seed);
	

	if (All.option == 0) //generate density files
	{
		cout << "Sorry, density estimation not available yet" << endl;
	}

	if ((All.fSample > 0.0) && (All.option == 1)) //main catalog generation
	{
		cout << "Generating catalog. Max allowed stars is " << numeric_limits<uint64_t>::max() << endl;
		cout << "Source numbering will start at " << All.nstart <<endl;
		string fname = All.outputDir + All.outputFile + ".ebf";
		cout << "Writing to " << fname << endl;
		SurveyDesign sur(fname, All.appMagLimits, All.absMagLimits,
				All.colorLimits, &All);
		sur.setGeometry(All.geometryOption, All.longitude, All.latitude, All.surveyArea);
		cout<<"setting center to";
		for(size_t i=0;i<6;i++){
			cout<<" "<<All.posC[i];
		}
		cout<<endl;
		sur.setCenter(All.posC);
		srand48(7);
		//----------------------------------------------------
		timer2.start();
		IsochroneDB ic(All.inputDir + "Isochrones/", All.photoCateg + "/", All.photoSys, All.magcolorNames, 1);
		ic.print();
		timer2.print("Time Isochrone Reading");
		//-----------------------------------------------------
		timer2.start();


		if(All.fieldTableFile.size()==0)
		{
			runModel(sur,ic,All);
		}
		else
		{
			for(size_t i=0;i<All.fieldTable.col[0].size();++i)
			{
				cout<<i<<" long="<< All.fieldTable.col[0][i]<<" lat="<<All.fieldTable.col[1][i]<<endl;
				sur.setGeometry(All.geometryOption, All.fieldTable.col[0][i], All.fieldTable.col[1][i], All.surveyArea,i);
				runModel(sur,ic,All);
			}

		}

		sur.close();

		cout << "Calculating magnitudes in "<<All.photoSys<<" system................" << endl;
		append1(All.outputDir + All.outputFile + ".ebf", ic,
						All.inputDir + "Isochrones/" + All.photoCateg + "/", All.photoSys,
						All.magcolorNames);

		cout << "Appending spherical coordinates................" << endl;
		append_radec(All.outputDir + All.outputFile, All.inputDir);
		append_lb(All.outputDir + All.outputFile, All.inputDir);


	}

	if (All.option == 2) //command line argument -a --psys="name_of_photo_system"
	{
		if (All.addstring == "")
		{
			cout << "Appending magnitudes in "<<All.photoSys<<" system................" << endl;
			//----------------------------------------------------
			timer2.start();
			IsochroneDB ic(All.inputDir + "Isochrones/", All.photoCateg + "/",
					All.photoSys, All.magcolorNames, 1);
			ic.print();
			timer2.print("Time Isochrone Reading");
			//----------------------------------------------------
			append1(All.outputFile + ".ebf", ic,
						All.inputDir + "Isochrones/" + All.photoCateg + "/", All.photoSys,
						All.magcolorNames);

		}
		else if (All.addstring == "radec")
		{
			append_radec(All.outputFile, All.inputDir);
		}
		else if (All.addstring == "lb")
		{
			append_lb(All.outputFile, All.inputDir);
		}
	}

	if (All.option == 3) //command line argument -e
	{
		cout << "Calulating Extinction................" << endl;
		timer2.start();
		append_ext(All.outputFile, All.inputDir);
		timer2.print("Time for extinction calculation");
	}

	timer1.print("Total Time=");

	return 0;
}






