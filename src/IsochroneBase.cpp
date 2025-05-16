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

#include<algorithm>
#include "IsochroneBase.h"
#include "TableInterpolator.h"
#include "ebf.hpp"



IsoFileDescriptor::IsoFileDescriptor(const string &fname,const string& photoSys,const string& magcolorNames)
{

	ifstream fd;
//	cout<<fname<<endl;
	fd.open(fname.c_str());
	if (fd.is_open())
	{
		string s;
		int nmags1;

		while(fd.good()==1)
		{
			getline(fd,s);
			stringstream ss(s);
			ss>>photoSysName>>fields>>startid>>nmags1;
			{
				vector<string> sv1;
				stringSplit(s," ",sv1);
				if(int(sv1.size()) != (nmags1+4))
				{
					cout << "scan problem in IsoFileDescriptor no of params, expected " << nmags1+4 << ", received " << int(sv1.size()) << endl;
					exit(1);
				}
			}
			magnames.clear();
			for(int i=0;i<nmags1;++i)
			{
				ss>>s;
				magnames.push_back(s);
//							cout<<i<<" "<<s<<" "<<ss.good()<<" "<<ss.eof()<<endl;
			}
//			cout<<photoSysName<<" "<<fields<<" "<<startid<<" "<<nmags1<<" "<<magnames.size()<<" "<<ss.good()<<endl;
			if((photoSysName==photoSys)||(photoSysName.find(photoSys)))
				break;
		}

	    if(fd.good()!=1)
		{
			cout << "scan problem in IsoFileDescriptor: Incorrect no of params" << endl;
			exit(1);
		}
	    fd.close();

	    vector<string> sv;
	    stringSplit(magcolorNames," ,-",sv);
	    certify(sv.size()==3,"mag color string not correct");
//	    cout<<sv[0]<<" "<<sv[1]<<" "<<sv[2]<<endl;

	    if(sv[0]=="?")
	    	sv[0]=magnames[0];
	    if(sv[1]=="?")
	    	sv[1]=magnames[1];
	    if(sv[2]=="?")
	    	sv[2]=magnames[0];


	    for(size_t i=0;i<sv.size();++i)
	    {
	    	vector<string>::iterator it=find(magnames.begin(),magnames.end(),sv[i]);
	    	certify(it!=magnames.end(),"color or magnitude name not in list");
	    	magid[i]=startid+it-magnames.begin();
//	    	cout<<i<<" "<<magid[i]<<endl;
	    }
//    	cout<<sv[0]<<"("<<magid[0]<<"), "<<sv[1]<<"("<<magid[1]<<")-"<<sv[2]<<"("<<magid[2]<<")"<<endl;

		int fields_r = *max_element(magid, magid + 3);
		if ((fields_r < 0) || (fields_r >= fields) || (fields_r >= FIELDS_MAX))
		{
			cout << "Something wrong with field descriptor " << fields << " "
					<< fields_r << endl;
			exit(1);
		}
	}
	else
	{
		cout << "Isochrone File Descriptor \"" << fname << "\" not found" << endl;
		exit(1);
	}


}



void IsochroneBase::print()
{
	cout<<left<<setw(36)<<"Isochrone Grid Size:"<<"(Age bins="<<Age.size()<<",Feh bins="<<FeH.size()<<",Alpha bins="<<Alpha.size()<<")"<<endl;
//	printv(Age,"Age");
//	printv(FeH,"FeH");
//	printv(Alpha,"Alpha");
}




IsochroneBase::IsochroneBase(const string& inputDir,const string& dirname,const string& photoSys,int extraFieldsOn,const string &magcolorNames)
{
	magswap=0;
	im.resize(300);
	for(int i=0;i<300;++i)
		im[i]=i;
	readIsochrones(inputDir,dirname,photoSys,extraFieldsOn,magcolorNames);
}


IsochroneBase::~IsochroneBase()
{
	// TODO Auto-generated destructor stub
}

void IsochroneBase:: readIsochrones(const string& inputDir,const string& dirname,const string& photoSys,int extraFieldsOn,const string &magcolorNames)
{
//	int dwarfOn=1;
	int dwarfOn=0;
	cout<<left<<setw(36)<<"Reading Isochrones from dir- "<<inputDir+dirname<<photoSys<<endl;
	if (ebf::ebfutils::FileExists(inputDir+"BolCorr/"+photoSys+"/interp_keys.txt") == 0)
	{
		dwarfOn=0;
	}

	double  al=0.0;

	double fe1[]={0.0001,0.0002,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.0012,0.0014,0.0016,0.0018,0.002,0.0022,0.0024,0.0026,0.003,0.0034,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.012,0.014,0.016,0.018,0.02,0.024,0.028,0.03};
	vector<double> fe(fe1,fe1+sizeof(fe1)/sizeof(fe1[0]));
	for(size_t i=0;i<fe.size();++i)
	{
		if((dirname.compare("py_custom/")==0)||(photoSys.compare("WFIRST+LSST")==0)||(photoSys.compare("WFIRST")==0)||(photoSys.compare("WFIRST+HST__WFC3")==0)||(photoSys.compare("HST__WFC3")==0)||(photoSys.compare("GAIA__0")==0)||(photoSys.compare("GAIA__DR2")==0)||(photoSys.compare("CTIO__DECam")==0)||(photoSys.compare("LSST_DP0")==0)) {
			//cout<<"Using new Zsun";
			FeH.push_back(log10(fe[i]/0.0152));
		}
		else
			FeH.push_back(log10(fe[i]/0.019)); //changes for PARSEC isochrones - then Zsun = 0.0152
	}

	Alpha.push_back(0.0);
	string buf=inputDir+dirname+"/"+photoSys+"/IsoFileDescriptor.txt";



	IsoFileDescriptor iso_fileinfo(buf,photoSys,magcolorNames);
	iso_fileinfo.setExtraFields(extraFieldsOn);

	if(iso_fileinfo.magid[0]==iso_fileinfo.magid[1])
		magswap=1;
	if(iso_fileinfo.magid[0]==iso_fileinfo.magid[2])
		magswap=2;

	if(extraFieldsOn>0)
	{
		certify(iso_fileinfo.extraid.size()==(iso_fileinfo.magnames.size()),"extraid vs nmags match");
	}

	int n_ages=0;
	for(unsigned int i=0;i<FeH.size();i++)
	{
		stringstream sout;
		sout<<inputDir<<dirname<<photoSys<<"/output_"<<fixed<<setprecision(6)<<fe[i]<<".dat";
		//cout<<"Reading file "<<sout.str()<<endl;
		if(i==0)
			n_ages=readfile(sout.str(),al,fe[i],iso_fileinfo,dwarfOn);
		else
		{
			assert(n_ages==readfile(sout.str(),al,fe[i],iso_fileinfo,dwarfOn));
		}
	}


	//---------------------------------------------------------
	// changed from with readfile (last one might be missing)
	for(size_t i=0;i<icv.size();++i)
	{
		icv[i].setMon();
		icv[i].setTip();
//		cout<<"l302 "<<icv[i].m[0]<<" "<<icv[i].Vmon[0]<<endl;
	}

// few checks
	//cout<<icv.size()<<" "<<n_ages<<" "<<FeH.size()<<endl;
	assert(icv.size()==n_ages*FeH.size());
	for(size_t i=1;i<FeH.size();i++)
	{
		//cout<<"l241 "<<icv[0].age<<" "<<icv[i*n_ages].age<<endl;
		assert(icv[0].age==icv[i*n_ages].age);
		assert(icv[n_ages-1].age==icv[(i+1)*n_ages-1].age);
	}
// set Age vector
	assert(icv[0].FeH==icv[1].FeH);
	//double temp=icv[1].age-icv[0].age;
	for(int i=0;i<n_ages;i++)
		Age.push_back(icv[i].age);
		//Age.push_back(icv[0].age+i*temp);


}

int IsochroneBase:: readfile(const string& fname,double alpha1,double feH1,const IsoFileDescriptor &iso_fileinfo,int dwarfOn)
{
	char buf[512];
	char fmt[512];
	char* cptr;
    char* token;
    const char* delimiters = " \t";
	FILE* fd=NULL;
	char buf_check[512];
	vector<float> x(iso_fileinfo.fields,0.0);
	certify(int(x.size())>=16,"number of fields greater than equal to 16");
	certify(int(iso_fileinfo.extraid.size())<=(iso_fileinfo.fields-iso_fileinfo.startid),"number of magnames greater than number of available fields");
	Isochrone icData;
	icData.age=-1.0;

	if((iso_fileinfo.photoSysName.find("Python")!=string::npos)||(iso_fileinfo.photoSysName.compare("WFIRST+LSST")==0)||(iso_fileinfo.photoSysName.compare("WFIRST")==0)||(iso_fileinfo.photoSysName.compare("WFIRST+HST__WFC3")==0)||(iso_fileinfo.photoSysName.compare("HST__WFC3")==0)||(iso_fileinfo.photoSysName.compare("GAIA__0")==0)||(iso_fileinfo.photoSysName.compare("GAIA__DR2")==0)||(iso_fileinfo.photoSysName.compare("CTIO__DECam")==0)||(iso_fileinfo.photoSysName.compare("LSST_DP0")==0)) { 
		//cout<<"Using new Zsun";
		icData.FeH=log10(feH1/0.0152);
	}
	else
		icData.FeH=log10(feH1/0.019);  //changes for PARSEC isochrones - then Zsun = 0.0152

	icData.alpha=alpha1;
	icData.addDwarfs(int(iso_fileinfo.extraid.size()),dwarfOn);


	//cout<<"magids are "<<iso_fileinfo.magid[0]<<" "<<iso_fileinfo.magid[1]<<" "<<iso_fileinfo.magid[2]<<"\n";

	int j,i=0;
	int k=0;
	int iage;
	bool linage;
	if((fd=fopen(fname.c_str(),"r")))
	{
		while(1)
		{
			cptr=fgets(buf,512,fd);
			if(feof(fd))
				break;
			if(buf[0]=='#'){
				//cout<<"Skipping comments"<<endl;
			    continue;
			}
			else if(iso_fileinfo.photoSysName.find("Python")!=string::npos){
				// cout<<"using py-custom Isochrones"<<endl;
				iage=0; // index of age column
				linage=0; // flag if ages are linear instead of log
                k = 0;
                token = strtok(buf,delimiters);
                while(k<iso_fileinfo.fields)
                {
                    sscanf(token,"%f", &x[k]);
                    token = strtok(NULL,delimiters);
                    k++;
                }
			}
			else if((iso_fileinfo.photoSysName.compare("Roman")==0)||(iso_fileinfo.photoSysName.compare("Euclid")==0)||(iso_fileinfo.photoSysName.find("JWST")!=string::npos)){
				// cout<<"using py-custom Isochrones"<<endl;
				iage=2; // index of age column
				linage=0; // flag if ages are linear instead of log
                k = 0;
                token = strtok(buf,delimiters);
                while(k<iso_fileinfo.fields)
                {
                    sscanf(token,"%f", &x[k]);
                    token = strtok(NULL,delimiters);
                    k++;
                }
			}
			else if(iso_fileinfo.photoSysName.compare("WFIRST")==0){
				//cout<<"using new WFIRST Isochrones"<<endl;
				iage=1;
				linage=1;
				sscanf(buf,"%f%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20],&x[21],&x[22],&x[23],&x[24],&x[25],&x[26],&x[27],&x[28],&x[29]);
			}
			else if(iso_fileinfo.photoSysName.compare("WFIRST+HST__WFC3")==0){
				//cout<<"using WFIRST-HST Isochrones"<<endl;
				iage=1;
				linage=0;
				sscanf(buf,"%f%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20],&x[21],&x[22],&x[23],&x[24],&x[25],&x[26]);
			}
			else if(iso_fileinfo.photoSysName.compare("LSST")==0){
				//cout<<"using LSST Isochrones"<<endl;
				iage=1;
				linage=1;
				sscanf(buf,"%f%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20],&x[21],&x[22],&x[23],&x[24],&x[25],&x[26],&x[27],&x[28]);
			}
			else if(iso_fileinfo.photoSysName.compare("CTIO__DECam")==0){
				//cout<<"using DECAM Isochrones"<<endl;
				iage=2;
				linage=0;
				sscanf(buf,"%f%f%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%e%f%e%e%e%e%e%f%e%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20],&x[21],&x[22],&x[23],&x[24],&x[25],&x[26],&x[27],&x[28],&x[29],&x[30],&x[31],&x[32],&x[33]);
			}
			else if(iso_fileinfo.photoSysName.compare("LSST_DP0")==0){
				//cout<<"using LSST_DP0 Isochrones"<<endl;
				iage=2;
				linage=0;
				sscanf(buf,"%f%f%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%e%f%e%e%e%e%e%f%e%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20],&x[21],&x[22],&x[23],&x[24],&x[25],&x[26],&x[27],&x[28],&x[29],&x[30],&x[31],&x[32],&x[33]);
			}
			else if(iso_fileinfo.photoSysName.compare("HST__WFC3")==0){
				//cout<<"using HST__WFC3 Isochrones"<<endl;
				iage=1;
				linage=1;
                k = 0;
                token = strtok(buf,delimiters);
                while(k<iso_fileinfo.fields)
                {
                    sscanf(token,"%f", &x[k]);
                    token = strtok(NULL,delimiters);
                    k++;
                } // TODO testing here if new implementation works well and can be generalized
				// sscanf(buf,"%f%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20],&x[21],&x[22],&x[23],&x[24],&x[25]);
			}
			else if(iso_fileinfo.photoSysName.compare("GAIA__0")==0){
				//cout<<"using Gaia Isochrones"<<endl;
				iage=1;
				linage=1;
				sscanf(buf,"%f%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20],&x[21],&x[22],&x[23],&x[24],&x[25]);
			}
			else if(iso_fileinfo.photoSysName.compare("GAIA__DR2")==0){
				//cout<<"using Gaia Isochrones"<<endl;
				iage=1;
				linage=1;
				sscanf(buf,"%f%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20],&x[21],&x[22],&x[23],&x[24],&x[25]);
			} else{
				iage=0;
				linage=0;
				sscanf(buf,"%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15]);
			}
			// if((i==0)){
			// 	cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<endl;
			// }

			//cout<<"iage "<<iage<<" linage "<<linage<<endl;
			sprintf(buf_check,"%f ",x[2]);
			if(strncmp(buf_check,"nan",3)==0)
				continue;

			//condition added because some isochrones use linear rather than log age in years
			if(linage){
				x[iage]=log10(x[iage]);
			}

			//check if we are on a new age in the table
			if(double(x[iage])!=icData.age)
			{
				icData.age=x[iage];
				//cout<<"Finished with age "<<x[iage]<<" index "<<iage<<endl;
				icv.push_back(icData);
				i++;
			}

			if(iso_fileinfo.photoSysName.find("Python")!=string::npos) {
				icv.back().m.push_back(x[1]);
				icv.back().Mact.push_back(x[2]);
				icv.back().Lum.push_back(x[3]);
				icv.back().Teff.push_back(x[4]);
				icv.back().Grav.push_back(x[5]);
			}
			else if((iso_fileinfo.photoSysName.compare("Roman")==0)||(iso_fileinfo.photoSysName.compare("Euclid")==0)||(iso_fileinfo.photoSysName.find("JWST")!=string::npos)) {
				icv.back().m.push_back(x[3]);
				icv.back().Mact.push_back(x[5]);
				icv.back().Lum.push_back(x[6]);
				icv.back().Teff.push_back(x[7]);
				icv.back().Grav.push_back(x[8]);
			}
			else if((iso_fileinfo.photoSysName.compare("WFIRST+HST__WFC3")==0)||(iso_fileinfo.photoSysName.compare("WFIRST")==0)||(iso_fileinfo.photoSysName.compare("HST__WFC3")==0)||(iso_fileinfo.photoSysName.compare("LSST")==0)||(iso_fileinfo.photoSysName.compare("GAIA__0")==0)||(iso_fileinfo.photoSysName.compare("GAIA__DR2")==0)) {
				icv.back().m.push_back(x[2]);
				icv.back().Mact.push_back(x[3]);
				icv.back().Lum.push_back(x[4]);
				icv.back().Teff.push_back(x[5]);
				icv.back().Grav.push_back(x[6]);
			}
			else if ((iso_fileinfo.photoSysName.compare("CTIO__DECam")==0)||(iso_fileinfo.photoSysName.compare("LSST_DP0")==0)) {
				icv.back().m.push_back(x[3]);
				icv.back().Mact.push_back(x[5]); // M_act does not exist in the new isochrone
				icv.back().Lum.push_back(x[6]);
				icv.back().Teff.push_back(x[7]);
				icv.back().Grav.push_back(x[8]);
			}
			else {
				icv.back().m.push_back(x[1]);
				icv.back().Mact.push_back(x[2]);
				icv.back().Lum.push_back(x[3]);
				icv.back().Teff.push_back(x[4]);
				icv.back().Grav.push_back(x[5]);
			}

			icv.back().Mag0.push_back(x[iso_fileinfo.magid[0]]);
			icv.back().Mag1.push_back(x[iso_fileinfo.magid[1]]);
			icv.back().Mag2.push_back(x[iso_fileinfo.magid[2]]);

			for(size_t k=0;k<iso_fileinfo.extraid.size();++k)
			{
				icv.back().Mags[k].push_back(x[iso_fileinfo.extraid[k]]);
			}

		}
		fclose(fd);
	}
	else
	{
		cout<<"file not found "<<fname<<endl;
		exit(1);
	}
	return i;
}


