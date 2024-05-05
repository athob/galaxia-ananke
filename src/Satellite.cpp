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

#include <utility>
#include <iostream>
#include "Satellite.h"
#include "Functions.h"
#include "ebfvector.hpp"

#define ALL(A) A.begin,A.end()
using namespace std;

vector<int> lindgen(int N)
{
	vector<int> x(N);
	for(int i=0;i<N;++i)
		x[i]=i;
	return x;
}



void Satellite::print()
{
	cout << "Satellite Info" << endl;
//	cout << "M_tot " << m_tot << endl;
//	cout << "Feh min max array size" << (*min_element(feh.begin(), feh.end()))
//			<< " " << (*max_element(feh.begin(), feh.end())) << " "
//			<< feh.size() << endl;
//	cout << "Feh min max array size" << (*min_element(feh_a.begin(),
//			feh_a.end())) << " " << (*max_element(feh_a.begin(), feh_a.end()))
//			<< " " << feh_a.size() << endl;
//	cout << "Age min max array size" << (*min_element(age.begin(), age.end()))
//			<< " " << (*max_element(age.begin(), age.end())) << " "
//			<< age.size() << endl;
//	cout << "Age min max array size" << (*min_element(age_a.begin(),
//			age_a.end())) << " " << (*max_element(age_a.begin(), age_a.end()))
//			<< " " << age_a.size() << endl;
//	cout<<"hdim"<<hdim<<endl;
//	printv(age_a.begin(), age_a.end(), "age_a");
//	printv(faca_a.begin(), faca_a.end(), "faca_a");
//	printv(feh_a.begin(), feh_a.end(), "feh_a");
	//	printv(m1_a.begin(), m1_a.end(), "m1_a");
	//	printv(m2_a.begin(), m2_a.end(), "m2_a");
	//	printv(fac_a.begin(), fac_a.end(), "fac_a");
}

void Satellite::initialize_stellar_data(Sampler& imf,IsochroneDB &ic)
{
	if(debug) cout <<"Satellite Initializing ....... ";
	age.resize(nsize);
	feh.resize(nsize);
	alpha.resize(nsize);
	int k = 0;
	for (ParticleStar* Part = pBegin; Part < pEnd; Part++)
	{
		age[k] = Part->Age();
//		if(k==0) cout<<"Age="<<age[k];
//		age[k] = log10((Part->Age()) * 1.0e9);
//		feh[k] = 0.019 * pow(10.0, Part->FeH());
		feh[k] = Part->FeH();
		alpha[k] = Part->Alpha();
		m_tot += Part->Mass();
		k++;
	}
//	m_tot *= 2.0;

	alpha2 = alpha;
	feh2 = feh;
	//sort feh and alpha by feh and make an iterator that returns feh and alpha in order. leftover, what is it used for now?
	sort2(feh2.begin(), feh2.end(), alpha2.begin());
	vector<double>::iterator it1 = unique2(feh2.begin(), feh2.end(), alpha2.begin());
	alpha2.resize(it1 - feh2.begin());
	feh2.resize(it1 - feh2.begin());

	//this part sorts feh and alpha by age and creates a sampler that I don't think is used anywhere now...but might be
	vector<int> ind=lindgen(age.size());
	sorti2(age.begin(), age.end(), ind.begin());
	orderv(feh,ind);
	orderv(alpha,ind);
	vector<double> age1 = age;
	vector<double> feh1 = feh;
	ageSampler = new Sampler(age);

 
	//another couple of sorting and iteration on unique values
	ind=lindgen(age.size());
	vector<double>::iterator it = unique2(age.begin(), age.end(), ind.begin());
	orderv(feh,ind);
	orderv(alpha,ind);
	feh.resize(it - age.begin());
	alpha.resize(it - age.begin());
	age.resize(it - age.begin());

	ind=lindgen(age.size());
	it = unique2(feh.begin(), feh.end(), ind.begin());
	orderv(age,ind);
	orderv(alpha,ind);
	age.resize(it - feh.begin());
	alpha.resize(it - feh.begin());
	feh.resize(it - feh.begin());


	//linear (in log space) age range spaced by 0.02
	cout <<"linear generator from "<<(*min_element(age.begin(), age.end()))<<" to "<<(*max_element(age.begin(), age.end()))<<" (log yr)"<<endl;
	linspace(age_a, (*min_element(age.begin(), age.end())), (*max_element(age.begin(), age.end())), 0.02);

	//interpolate feh on age range in linear generator, using the age and feh arrays built above
	interpolate(age, feh, age_a, feh_a);

	//construct interpolation of the cumulative distribution of ages -- aagh this is the problem
	interpolate(ageSampler->x, ageSampler->cpd, age_a, faca_a);

	if(debug) cout << "Done" << endl;
}

Satellite::~Satellite()
{
	delete ageSampler;
	ageSampler = NULL;
	delete[] data1;
	delete[] part;
}

double kernel_func(double x)
{
	double x2=x*x;
	return (1-x2)*x2*x2*x;
}

void Satellite::spawn1(SurveyDesign &sur, Sampler& imf,IsochroneDB &ic,double fSample,int seed1)
{
	print();
	double distMod,distModMax;
	double mmin, mmax, frac;
	//initialize_stellar_data(imf, ic);  //not needed for simulated star particles
	Sampler kerSampler(1000,0.0,1.0,kernel_func,0);
	double pr[6];
	Normaldev randomn(0.0,1.0,seed1);
	Ran randomu(seed1+12);

	StarParticle Star,Star1;
	vector<double> m3_a;
	m3_a=m2_a;
	nstars = 0;
	nstars_rejected = 0;
	int nparts_processed = 0;
	int nstars_generated = 0;
	l_tot_v = 0.0;
	cout<<"Particles="<< pEnd - pBegin<<" Mass="<<m_tot <<" "<<imf.meanx<< endl;

//	cout<<"#index mass mmin mmax frac nstars"<<endl;
	for (ParticleStar* Part = pBegin; Part < pEnd; Part++)
	{
		for(int i=0;i<6;++i)
			Star.pos(i)=Part->Pos[i]-sur.posC[i];
		distMod = distance(Part->Pos, sur.posC, 3) - Part->h(0); //closest possible distance within kernel
		distModMax = distance(Part->Pos, sur.posC, 3) + Part->h(0); //furthest possible distance within kernel
		//checks if star *particle* is within 1 kernel size of distance limits
		if((sur.geo->checkSphere(&(Star.pos(0)),Part->h(0)) == 1)&&(distMod<=sur.r_max)&&(distModMax>sur.r_min))
		{
			nparts_processed++;
			distMod = 5.0 * log10(max(distMod,0.01) * 100.0); //largest possible distance modulus given kernel size

//			cout<<Part->parentID()<<" "<<Part->Mass()<<" ";

			//calculates the min and max stellar masses for the magnitude range allowed, as a function of age and metallicity
			//this defines the fraction of the IMF to be sampled ...but I don't know that we need this any more. We should just be sampling the *whole* IMF
			//ic.min_max_m(age_a, feh_a, feh_a, m1_a, m2_a, min(sur.absMag[1],sur.appMag[1]-distMod)+0.6,sur.All->starType); 

			//calculate the min and max stellar masses for the magnitude range allowed, given age and metallicity of star *particle*
			//to avoid wasting time sampling stars with masses that don't exist anymore (b/c stellar remnant) or will be too faint to fall in the survey
			ic.min_max_m_new(Part->Age(), Part->FeH(), Part->Alpha(), &mmin, &mmax, min(sur.absMag[1],sur.appMag[1]-distMod)+0.6,sur.All->starType); 

			//calculate the fraction of the total IMF that we are sampling, given faint mag limit
			frac = imf.getFac(mmin, mmax);

			//use this to estimate how many stars we should sample given the mass of the parent particle
			//assuming they are all of the median stellar mass for the IMF
			double nstars_float = (Part->Mass() / imf.meanx) * frac * sur.All->fSample;

			//use stochastic rounding to determine the integer number of stars to include
			int stars = int(nstars_float);
			if (randomu.doub() <= (nstars_float - stars))
				stars++;

//			cout<<mmin<<" "<<mmax<<" "<<frac<<" "<<stars<<endl;

			//if some stars will be included, set up the sampler to sample only this range of the IMF
			if (stars > 0) {
				imf.setRange(mmin,mmax);
			}

			nstars_generated+=stars;

			//sample star particles from parent particle
			int status=0;
			for (int l = 0; l < stars; ++l)
			{
				for(int i=0;i<6;++i)
					Star.pos(i)=Part->Pos[i]-sur.posC[i];  //start by placing each particle at the location of the parent


				//compute a random unit vector in 6D by sampling the normal distribution six times
				double temp=0.0,temp1=0.0;
				for(int i=0;i<hdim;++i)
				{
					pr[i]=randomn.dev();
					temp1+=pr[i]*pr[i];
				}
				//vector norm such that ehat = pr/temp1
				temp1=sqrt(temp1);
				
				//sample radius from Epanechnikov kernel with radius 1
				temp=kerSampler.rand();


				//sample x and v using 2x3D kernels h(0) and h(1)
				for(int i=0;i<3;++i)
					Star.pos(i)+=Part->h(0)*pr[i]*temp/temp1;
				if(hdim==6)
				{
					for(int i=0;i<3;++i)
						Star.pos(i+3)+=Part->h(1)*pr[i+3]*temp/temp1;
				}


				//in original code the age, fe/h and alpha were being sampled/interpolated. 
				//changed to make them exactly equal to the parent particle
				//Star.age() = ageSampler->rand();  //age from random sample in 1D
				Star.age() = Part->Age();
				//interpolate(age, feh, Star.age(), Star.feh());  //feh from 2D interp with given age
				Star.feh() = Part->FeH();
				Star.smass() = imf.rand(); //initial mass of star, sampled from truncated IMF
				//Star.alpha() = 0.0;
				Star.alpha() = Part->Alpha();
				Star.mag(0) = 0.0;
				Star.mag(1) = 0.0;
				Star.mag(2) = 0.0;
				Star.satID() = satID;
				Star.popID() = nstars + total(sur.npart) + sur.nstart;
				Star.partID() = 1;
				Star.parentID() = Part->parentID(); //keep track of which simulation particle spawned this star
				Star.partitionID() = Part->partitionID(); //keep track of which partition the particle is part of
				Star.dform() = Part->dform(); //formation distance of parent particle
				ic.interpolateStar(Star);

//	carry forward abundances from parent particle
				Star.helium() = Part->helium();
				Star.carbon() = Part->carbon();
				Star.nitrogen() = Part->nitrogen();
				Star.oxygen() = Part->oxygen();
				Star.neon() = Part->neon();
				Star.silicon() = Part->silicon();
				Star.sulphur() = Part->sulphur();
				Star.calcium() = Part->calcium();
				Star.magnesium() = Part->magnesium();

//					cout << Star.parentID() << " ";

				//a status of zero means this is the first star sampled from this parent particle, 
				//in which case we place it exactly at the position of the parent
				if(status==0)
				{
					Star1=Star;
					for(int i=0;i<6;++i)
						Star1.pos(i)=Part->Pos[i]-sur.posC[i];
					Star1.AbsToAppMag();
					Star1.AddPhotoError(randomn,sur.All->photoError,ic.magswap);
					Star1.partID() = 0;
				}

				Star.AbsToAppMag();
				Star.AddPhotoError(randomn,sur.All->photoError,ic.magswap);


				if(status==0)
				{
					// star placed at position of parent is only included if *both* the parent position 
					// *and* what would have been the perturbed position of star
					// are in the survey volume ... if not, try again to place a star at the parent particle position
					if (sur.push_check(Star))
					{
						if(sur.push_back(Star1))
						{
							nstars++;
							status=1;
						}
						else if(sur.push_back(Star))
						{
							nstars++;
						}
						else
							nstars_rejected++;
					}
					else
						nstars_rejected++;

				}
				else
				{
					if (sur.push_back(Star))
					{
						nstars++;
					}
					else
						nstars_rejected++;
				}

				if (((nstars + nstars_rejected) % 10000000) == 0)
				{
					cout<< nstars<< " accepted  " << nstars_rejected<<" rejected Parts="<<nparts_processed<<" outof "<<Part-pBegin<<endl;
				}
			}
		}
	}

	cout << "Total Stars="<<nstars_generated<<" accepted="<< nstars << " rejected=" << nstars_rejected << endl;
	sur.flush();
}



void Satellite::initialize(int size1)
{
	nsize = size1;
	items = 0;
	for (size_t i = 0; i < tags.size(); ++i)
	{
		items += tags[i].dim;
	}
	data1 = new double[nsize * items];
	data = data1;
	part = new ParticleStar[nsize];
	for (int i = 0; i < nsize; i++)
		part[i].Pos = &data[i * items];

	vector<ParticleTag>::iterator it;
	it = find_if(tags.begin(), tags.end(), ParticleTag_eq("/H_smooth", 0));
	if (it != tags.end())
		for (long i = 0; i < long(nsize); i++)
		{
			part[i].Pos[it->id] = 1.0;
			part[i].Pos[it->id + 1] = 1.0;
		}
	pBegin = part;
	pEnd = part + nsize;
	if(debug) cout<<"Particles="<<nsize<<endl;
}

Satellite::Satellite(const string& fname, int sat_no, int satID1,int hdim1,int nres1)
{
	hdim=hdim1;
	nres=nres1;
	satID = satID1;
	pBegin = NULL;
	pEnd = NULL;
	items = 0;
	nsize = 0;
	m_tot = 0.0;
	nstars = 0;
	debug=1;
	ageSampler = NULL;
	ParticleStar s;
	s.getTags(tags, hdim);
	readEbfFile(fname, sat_no);
}




void Satellite::readEbfFile(const string &fname, int sat_no)
{
	ebf::EbfVector<float> fb_pos3(fname,"/Pos3");
	ebf::EbfVector<float> fb_vel3(fname,"/Vel3");
	ebf::EbfVector<float> fb_feh(fname,"/FeH");
	ebf::EbfVector<float> fb_alpha(fname,"/Alpha");
	ebf::EbfVector<float> fb_age(fname,"/Age");
	ebf::EbfVector<float> fb_mass(fname,"/Mass");
	ebf::EbfVector<int> fb_parentID(fname,"/ParentID");
	ebf::EbfVector<int> fb_partitionID(fname,"/PartitionID");
	ebf::EbfVector<float> fb_dform(fname,"/Dform");

	//eventually replace this method with e.g. readGizmoFile(fname)

	ebf::EbfVector<float> fb_helium(fname,"/Helium");
	ebf::EbfVector<float> fb_carbon(fname,"/Carbon");
	ebf::EbfVector<float> fb_nitrogen(fname,"/Nitrogen");
	ebf::EbfVector<float> fb_oxygen(fname,"/Oxygen");
	ebf::EbfVector<float> fb_neon(fname,"/Neon");
	ebf::EbfVector<float> fb_magnesium(fname,"/Magnesium");
	ebf::EbfVector<float> fb_silicon(fname,"/Silicon");
	ebf::EbfVector<float> fb_sulphur(fname,"/Sulphur");
	ebf::EbfVector<float> fb_calcium(fname,"/Calcium");

	initialize(fb_mass.size());


	{
		ebf::EbfVector<int> fb_satid(fname,"/id");
		satID=fb_satid[0];
	}


	m_tot =0.0;
	for(int i=0;i<nsize;++i)
	{
		part[i].Pos[0]=fb_pos3[i*3+0];
		part[i].Pos[1]=fb_pos3[i*3+1];
		part[i].Pos[2]=fb_pos3[i*3+2];
		part[i].Pos[3]=fb_vel3[i*3+0];
		part[i].Pos[4]=fb_vel3[i*3+1];
		part[i].Pos[5]=fb_vel3[i*3+2];
		part[i].FeH()=fb_feh[i];
		part[i].Alpha()=fb_alpha[i];
		part[i].Mass()=fb_mass[i];
		part[i].Age()=fb_age[i];
		part[i].parentID()=fb_parentID[i];
		part[i].partitionID()=fb_partitionID[i];
		part[i].dform() = fb_dform[i];
		part[i].helium() = fb_helium[i];
		part[i].carbon() = fb_carbon[i];
		part[i].nitrogen() = fb_nitrogen[i];
		part[i].oxygen() = fb_oxygen[i];
		part[i].neon() = fb_neon[i];
		part[i].magnesium() = fb_magnesium[i];
		part[i].silicon() = fb_silicon[i];
		part[i].sulphur() = fb_sulphur[i];
		part[i].calcium() = fb_calcium[i];

		m_tot+=fb_mass[i];
	}

/*	for(int i=0;i<10;++i) {
		cout << part[i].parentID() << " ";
	}
	cout << endl;
*/
	string s(fname);
	if(hdim==3)
		if(nres==32)
			s.insert(s.find(".ebf"),"_d3n32_den");
		else if(nres==64)
			s.insert(s.find(".ebf"),"_d3n64_den");
		else if(nres==128)
			s.insert(s.find(".ebf"),"_d3n128_den");
		else
			s.insert(s.find(".ebf"),"_d3n8_den");
	else
		if(nres==32)
			s.insert(s.find(".ebf"),"_d6n32_den");
		else if(nres==64)
			s.insert(s.find(".ebf"),"_d6n64_den");
		else if(nres==128)
			s.insert(s.find(".ebf"),"_d6n128_den");
		else
			s.insert(s.find(".ebf"),"_d6n8_den");

	ebf::EbfVector<float> fb_hcubic(s,"/H_cubic");
	if(fb_hcubic.rank()==1)
	{
		for(int i=0;i<nsize;++i)
		{
			part[i].h(0)=fb_hcubic[i];
			part[i].h(1)=0.0;
		}
	}
	if(fb_hcubic.rank()==2)
	{
		for(int i=0;i<nsize;++i)
		{
			part[i].h(0)=fb_hcubic[i*2+0];
			part[i].h(1)=fb_hcubic[i*2+1];
		}
	}

}

