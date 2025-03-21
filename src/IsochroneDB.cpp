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

#include "IsochroneDB.h"
#include "Functions.h"




IsochroneDB::~IsochroneDB()
{
	// TODO Auto-generated destructor stub
}

void IsochroneDB::calculateAgeProb(double Age2,double dAge2)  // seems unused
{
	double Age1[2];
	Age1[0]=log10(Age2-dAge2);
	Age1[1]=log10(Age2+dAge2);
	double temp_up,temp_lo,temp;
	temp=0.0;

	if(dAge2>0)
	{
	for(size_t i=0;i<Age.size();++i)
	{
		W_Age[i]=0;
		if(i==0)
			temp_lo=min(Age1[0],Age[i]);
		else
			temp_lo=(Age[i]+Age[i-1])*0.5;
		if(temp_lo > Age1[1])
			continue;

		if(i==(Age.size()-1))
			temp_up=max(Age1[1],Age[i]);
		else
			temp_up=(Age[i]+Age[i+1])*0.5;
		if(temp_up < Age1[0])
			continue;

		temp_lo=max(temp_lo,Age1[0]);
		temp_up=min(temp_up,Age1[1]);
		W_Age[i]=(pow(10.0,temp_up)-pow(10.0,temp_lo))/(2*dAge2);
//		cout<<Age1[0]<<" "<<Age1[1]<<" "<<Age[i]<<" "<<W_Age[i]<<" "<<temp_up<<" "<<temp_lo<<endl;
		temp+=W_Age[i];
	}
	// corrected from Feh.size()
	for(size_t i=0;i<Age.size();++i)
		W_Age[i]/=temp;
	}
	else
	{
		for(size_t i=0;i<Age.size();++i)
			W_Age[i]=0;
		size_t i=int((Age1[0]-Age[0])/dAge);
		if(i>=(Age.size()-1))
		{
			W_Age[i]=1.0;
		}
		else
		{
			assert((i>0)&&(i<(Age.size()-1)));
			W_Age[i]=1-(Age1[0]-Age[i])/(Age[i+1]-Age[i]);
			W_Age[i+1]=1-W_Age[i];
		}

	}

}

void IsochroneDB::calculateFeHProb(double FeH2,double dFeH2)  // seems unused
{
	double FeH1[2];
	FeH1[0]=FeH2-3*dFeH2;
	FeH1[1]=FeH2+3*dFeH2;
	double temp_up,temp_lo,temp;
	temp=0;
	if(dFeH2>0)
	{
	for(size_t i=0;i<FeH.size();++i)
	{
		W_FeH[i]=0;
		if(i==0)
			temp_lo=min(FeH1[0],FeH[i]);
		else
			temp_lo=(FeH[i]+FeH[i-1])*0.5;
		if(temp_lo > FeH1[1])
			continue;
		if(i==(FeH.size()-1))
			temp_up=max(FeH1[1],FeH[i]);
		else
			temp_up=(FeH[i]+FeH[i+1])*0.5;
		if(temp_up < FeH1[0])
			continue;
		temp_lo=max(temp_lo,FeH1[0]);
		temp_up=min(temp_up,FeH1[1]);
//		W_FeH[i]=(temp_up-temp_lo)/(FeH1[1]-FeH1[0]);
		W_FeH[i]=(gaussian.cdf((temp_up-FeH2)/dFeH2)-gaussian.cdf((temp_lo-FeH2)/dFeH2));
		temp+=W_FeH[i];
//		cout<<FeH1[0]<<" a "<<FeH1[1]<<" "<<FeH[i]<<" "<<W_FeH[i]<<" "<<FeH2<<endl;
	}
	for(size_t i=0;i<FeH.size();++i)
		W_FeH[i]/=temp;
	}
	else
	{
		for(size_t i=0;i<FeH.size();++i)
			W_FeH[i]=0;
		size_t i=locate(FeH,FeH2);
		W_FeH[i]=1-(FeH1[0]-FeH[i])/(FeH[i+1]-FeH[i]);
		W_FeH[i+1]=1-W_FeH[i];

	}
}




void IsochroneDB::setIsochrones(double Age2,double dAge2,double FeH2,double dFeH2,double maxV1,double mass)  // this is called nowhere
{
	box_mmin=1000.0;
	box_mmax=0.0;

	calculateAgeProb(Age2,dAge2);
	calculateFeHProb(FeH2,dFeH2);
	icNeighbor icN;
	icQ.clear();
	ic_p.clear();
	double temp;
	n_tot=0;
	if(fabs(maxV-maxV1)>delta_maxV)
		maxV=maxV1;
	temp=0.0;
	for(size_t p=0;p<Age.size();++p)
	{
//		cout<<"here1 "<<p<<" "<<W_Age[p]<<endl;
		if(W_Age[p]>0)
		for(size_t q=0;q<FeH.size();++q)
		{
//			cout<<"here2 "<<q<<" "<<W_FeH[q]<<endl;
			if(W_FeH[q]>0)
			{
				ic_p.push_back(temp);
				icN.ic=&icv[index(p,q,0)];
				if(icN.ic->maxV!=maxV)
					icN.ic->setmaxV(maxV);
				icN.prob=W_Age[p]*W_FeH[q]*imfP->getFac(icN.ic->m_min,icN.ic->m_max);

				if(icN.ic->m_min < box_mmin)
					box_mmin=icN.ic->m_min;
				if(icN.ic->m_max > box_mmax)
					box_mmax=icN.ic->m_max;

				icQ.push_back(icN);
				temp+=icN.prob;
			}
		}
	}
	ic_p.push_back(temp*1.00001);
	p_tot=temp;
	temp=temp*mass/imfP->meanx;
	n_tot=int(temp);
	if(grandomu()<(temp-n_tot))
		n_tot++;
}

void IsochroneDB::generateStar(StarParticle &Star)  // this is called nowhere
{

	double temp;
	temp=p_tot*grandomu();
	int i=upper_bound(ic_p.begin(),ic_p.end(),temp)-ic_p.begin()-1;

	assert((i>=0)&&(i<int(ic_p.size()-1)));
	imfP->setRange(icQ[i].ic->m_min,icQ[i].ic->m_max);
	Star.smass() = imfP->rand();
	int i_age,i_feh,i_al;
	reverse_index(i_age,i_feh,i_al,(icQ[i].ic-&icv[0]));

	if(i_age<=0)
		Star.age() = Age[0]+(Age[1]-Age[0])*grandomu()*0.5;
	else if(i_age>=int(Age.size()-1))
		Star.age() = Age[Age.size()-2]+(Age[Age.size()-1]-Age[Age.size()-2])*(0.5+grandomu()*0.5);
	else
		Star.age() = (Age[i_age-1]+Age[i_age])*0.5+(Age[i_age+1]-Age[i_age-1])*grandomu()*0.5;

	if(i_feh<=0)
		Star.feh() = FeH[0]+(FeH[1]-FeH[0])*grandomu()*0.5;
	else if(i_feh>=int(FeH.size()-1))
		Star.feh() = FeH[FeH.size()-2]+(FeH[FeH.size()-1]-FeH[FeH.size()-2])*(0.5+grandomu()*0.5);
	else
		Star.feh() = (FeH[i_feh-1]+FeH[i_feh])*0.5+(FeH[i_feh+1]-FeH[i_feh-1])*grandomu()*0.5;

	Star.alpha() = Star.feh();
	icQ[i].ic->interpolate(Star);

}


void IsochroneDB::setIsochrones1(double Age2,double dAge2,double FeH2,double dFeH2,double maxV1,double mass,int starType)  // this is called nowhere
{
	box_mmin=1000.0;
	box_mmax=0.0;
	if((last_age!=Age2)||(last_dage!=dAge2))
	{
		calculateAgeProb(Age2,dAge2);
		last_age=Age2;
		last_dage=dAge2;
	}
	calculateFeHProb(FeH2,dFeH2);



	icNeighbor icN;
	double temp;
	n_tot=0;
	if(fabs(maxV-maxV1)>delta_maxV)
		maxV=maxV1;
	temp=0.0;
	for(size_t p=0;p<Age.size();++p)
	{
//		cout<<"here1 "<<p<<" "<<W_Age[p]<<endl;
		if(W_Age[p]>0)
		for(size_t q=0;q<FeH.size();++q)
		{
//			cout<<"here2 "<<q<<" "<<W_FeH[q]<<endl;
			if(W_FeH[q]>0)
			{

				icN.ic=&icv[index(p,q,0)];
				if(icN.ic->maxV!=maxV)
					icN.ic->setmaxV(maxV);

				if(icN.ic->m_min < box_mmin)
					box_mmin=icN.ic->m_min;
				if(icN.ic->m_max > box_mmax)
					box_mmax=icN.ic->m_max;


			}
		}
	}

	box_mmin=max(imfP->x_min,box_mmin);

	if((starType==1)||(starType==2))
	{
		box_mmin=0.75;
		box_mmax=1.0;
	}


	temp=mass*imfP->getFac(box_mmin,box_mmax)/imfP->meanx;
	n_tot=int(temp);
	if(grandomu()<(temp-n_tot))
		n_tot++;
	imfP->setRange(box_mmin,box_mmax);

//	cout<<box_mmin<<" "<<box_mmax<<endl;

}


void IsochroneDB::interpolateStar(StarParticle &Star)  // used when initializing
{

	int i,j;
	float tempf=Star.age();
	// if(tempf<=Age.front())
	// 	i=0;
	// else if(tempf>=Age.back())
	// 	i=int(Age.size()-1);
	// else
	// 	i=int((tempf-Age[0])/dAge+0.5);
	i=locate_nearest(Age,double(tempf));

	tempf=Star.feh();
	j=locate_nearest(FeH,double(tempf));


	size_t k=index(i,j,0);
	if((k<0)||(k>=icv.size()))
	{
		cout<<"A grid index="<<k<<" is out of range max_size="<<icv.size()<<" i="<<i<<" j="<<j<<endl;
		cout<<"Age-->"<<Star.age()<<" "<<Age.front()<<" "<<Age.back()<<endl;
		cout<<"FeH-->"<<Star.feh()<<" "<<FeH.front()<<" "<<FeH.back()<<endl;
		exit(1);
	}
	icv[k].interpolate(Star);

}



void IsochroneDB::interpolateTGM(float age1,float feh1,float smass1,vector<double> &x)  // used when appending
{
	// int i;
	// if(age1<=Age.front())
	// 	i=0;
	// else if(age1>=Age.back())
	// 	i=int(Age.size()-1);
	// else
	// 	i=int((age1-Age[0])/dAge+0.5);
	int i=locate_nearest(Age,double(age1));

	int j=locate_nearest(FeH,double(feh1));
	size_t k=index(i,j,0);
	if((k<0)||(k>=icv.size()))
	{
		cout<<"A grid index="<<k<<" is out of range max_size="<<icv.size()<<" i="<<i<<" j="<<j<<endl;
		cout<<"Age-->"<<age1<<" "<<Age.front()<<" "<<Age.back()<<endl;
		cout<<"FeH-->"<<feh1<<" "<<FeH.front()<<" "<<FeH.back()<<endl;
		exit(1);
	}
	icv[k].interpolateTGM(smass1,x);
}

void IsochroneDB::min_max_m_new(double age, double feh, double alpha, double* m_min, double* m_max,
		double maxV, int OptRR)
{	
	int iage_up, iage_lo, iage;
	int ife_up, ife_lo, ifeh;
	int ii, imag;

	double dage1, dage2, dfe1, dfe2;

	if ((OptRR == 1)||(OptRR == 2))
		{
			*m_min = 0.75;
			*m_max = 1.0;

	} else {

		//find index of isochrones closest to age and feh in table

	// 	//this finds the index of the first element larger than the age and feh values
	// 	iage_up = int(upper_bound(Age.begin(), Age.end(), age) - Age.begin() - 1);
	// 	ife_up = int(upper_bound(FeH.begin(), FeH.end(), feh) - FeH.begin() - 1);
	// 	iage_lo = (iage_up>0 ? iage_up-1 : 0);
	// 	ife_lo = (ife_up>0 ? ife_up-1 : 0);
	// //	iage_up = max(iage_up, 0);
	// //	ife_up = max(iage_up, 0);
	// //	iage_lo = min(iage_lo,int(Age.end()-Age.begin()-1));
	// //	ife_lo = min(ife_lo, int(FeH.end()-FeH.begin()-1));
	// 	//determine which index points to the closer one
	// 	dage1 = fabs(Age[iage_up] - age);
	// 	dage2 = fabs(age - Age[iage_up]);
	// 	dfe1 = fabs(FeH[ife_up] - feh);
	// 	dfe2 = fabs(feh - FeH[ife_lo]);
	// 	iage = (dage1<dage2 ? iage_up : iage_lo);
	// 	ifeh = (dfe1<dfe2 ? ife_up : ife_lo); 

		iage = locate_nearest(Age, age);
		ifeh = locate_nearest(FeH, feh);

		//cout<<age<<" "<<feh<<" "<<iage<<" "<<ifeh<<" ";

		//find faint limit in magnitude sequence for this isochrone
		ii = index(iage,ifeh,0);
		//cout<<ii<<" "<<icv[ii].age<<" "<<icv[ii].FeH<<endl;
		imag = int(upper_bound(icv[ii].Vmon.begin(), icv[ii].Vmon.end(),maxV, greater<double> ()) - icv[ii].Vmon.begin() - 1);
		imag = max(min(int(icv[ii].Vmon.size()),imag),0);

		//set minimum stellar mass to the corresponding value
		*m_min = icv[ii].m[imag];

		//set maximum mass to the max for this isochrone
		*m_max = icv[ii].m.back();

		//cout<<*m_min<<" "<<*m_max<<" ";
		
	}
} 



void IsochroneDB::min_max_m(vector<double> &age_a, vector<double> &FeH_a,
		vector<double> &Alpha_a, vector<double> &m1_a, vector<double> &m2_a,
		double maxV, int OptRR)  // this is called nowhere
{
	int p1, q1, p2, q2, ii, i1;
	m1_a.resize(age_a.size(), 0.0);
	m2_a.resize(age_a.size(), 0.0);
	//	m3_a.resize(age_a.size(), 0.0);
	for (unsigned int j = 1; j < age_a.size(); ++j)
	{
		p1 = int(upper_bound(Age.begin() + 1, Age.end() - 1, age_a[j])
				- Age.begin() - 1);
		q1 = int(upper_bound(FeH.begin() + 1, FeH.end() - 1, FeH_a[j])
				- FeH.begin() - 1);
		p2 = int(upper_bound(Age.begin() + 1, Age.end() - 1, age_a[j - 1])
				- Age.begin() - 1);
		q2 = int(upper_bound(FeH.begin() + 1, FeH.end() - 1, FeH_a[j - 1])
				- FeH.begin() - 1);
		//		cout<<p1<<" "<<q1<<" "<<p2<<" "<<q2<<endl;
		if ((OptRR == 1)||(OptRR == 2))
		{
			m1_a[j] = 0.75;
			m2_a[j] = 1.0;

		}
		else
		{

			i1 = 0;
			for (int j1 = p1; j1 < p1 + 2; ++j1)
				for (int j2 = q1; j2 < q1 + 2; ++j2)
				{
					ii = index(j1, j2, 0);
					if (maxV < icv[ii].Vmon.back())
						i1++;
				}
			for (int j1 = p2; j1 < p2 + 2; ++j1)
				for (int j2 = q2; j2 < q2 + 2; ++j2)
				{
					ii = index(j1, j2, 0);
					if (maxV < icv[ii].Vmon.back())
						i1++;
				}

			if (i1 == 8)
			{
				ii = index(p1 + 1, q1, 0);
				m1_a[j] = icv[ii].m.back();
				m2_a[j] = m1_a[j];
			}
			else
			{
				// max age min Fe
				ii = index(p1 + 1, q1, 0);
				//safe method
				//			i1 = int(upper_bound(icv[ii].V.begin() + 1, (icv[ii].V.begin()
				//					+ icv[ii].tip1), maxV, greater<double> ())
				//					- icv[ii].V.begin() - 1);
				// risky
				i1 = int(upper_bound(icv[ii].Vmon.begin(), icv[ii].Vmon.end(),
						maxV, greater<double> ()) - icv[ii].Vmon.begin() - 1);
				if(i1<0)
				{
					m1_a[j] = icv[ii].m[0];
				}
				else
				{
//					assert((i1>=0)&&(i1<icv[ii].m.size()));
					m1_a[j] = icv[ii].m[i1];
				}

				// semi safe
				//						m1_a[j]=1000.0;
				//			for (unsigned int j1 = p1; j1 < p1 + 2; ++j1)
				//				for (unsigned int j2 = q1; j2 < q1 + 2; ++j2)
				//				{
				//					ii = index(j1, j2, 0);
				//					i1 = int(upper_bound(icv[ii].Vmon.begin(), icv[ii].Vmon.end(), maxV, greater<double> ())
				//							- icv[ii].Vmon.begin() - 1);
				//					if(m1_a[j] >icv[ii].m[i1]) 			m1_a[j] = icv[ii].m[i1];
				//				}
				//			m3_a[j] = icv[ii].m[i1];

				// min age max Fe
				ii = index(p2, q2 + 1, 0);
				m2_a[j] = icv[ii].m.back();
				//			cout<<m1_a[j]<<" "<<m2_a[j]<<endl;
			}
		}



	}
	m1_a[0] = m1_a[1];
	m2_a[0] = m2_a[1];

}


