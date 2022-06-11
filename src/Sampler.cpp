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



#include "Sampler.h"
#include "Functions.h"


Sampler::Sampler(int nsize,double xmin1,double xmax1,double (* func) (double x),int optionlinlog)
{
	 function=func;
	 if(optionlinlog==0)
		 linspace(x,xmin1,xmax1,nsize);
	 else
	 {
		 certify(xmax1>xmin1,"Smapler test failed xmin1>xmax1");
		 if(xmin1<=0)
		 {
			 cout<<"Sampler:log spacing not allowed for xmin<=0"<<endl;
			 exit(1);
		 }
		 logspace(x,xmin1,xmax1,nsize);
	 }
	px.resize(x.size());
	for(int i=0;i<int(x.size());++i)
		px[i]=function(x[i]);
	normalize();
}

void Sampler::print( )
{
	cout<<"Printing Sampler Data "<<endl;
	for(unsigned int i=0;i<x.size();++i)
		cout<<i<<" x="<<x[i]<<" px="<<px[i]<<" cpd="<<cpd[i]<<endl;
}

void Sampler::plot( )
{
//stub
}

Sampler::Sampler(const string fname)
{
	   float temp;
	    ifstream fd;
	    int i,no;

	    fd.open(fname.c_str());
	    if (fd.is_open())
	    {
	    	cout<<"Reading luminosity function file: "<<fname<<endl;
		fd>>no;
		for(i=0;i<no;++i)
		{
		    fd>>temp;
		    x.push_back(temp);
		}
		for(i=0;i<no;++i)
		{
		    fd>>temp;
		    px.push_back(temp);
		}
		fd>>temp;

		fd.close();
	    }
	    else
		cout<<"Error opening file"<<endl;
		   cout<<"Int P(l)dl "<<int_tabulated(x,px)<<endl;
		normalize();
}

Sampler::Sampler(vector<double> &x1,vector<double> &px1)
{
	x=x1;
	px=px1;
	normalize();
	calculateCpd();
}

Sampler::Sampler(vector<double> &x1)
{
	x=x1;
	sort(x.begin(),x.end());
	double temp;
	cpd.resize(x.size());
	for(unsigned int i=0;i<x.size();++i)
    {
	    cpd[i]=i*1.0/(x.size()-1);
    }
	vector<double>::iterator it=unique2(x.begin(),x.end(),cpd.begin());
	cpd.resize(it-x.begin());
	x.resize(it-x.begin());
	cpd.back()=1.0;
	px.resize(x.size());
	for(unsigned int i=1;i<(x.size()-1);++i)
    {
	    temp=(cpd[i]-cpd[i-1])/(x[i]-x[i-1])+(cpd[i+1]-cpd[i])/(x[i+1]-x[i]);
	    temp*=0.5;
	    px[i]=temp;
    }
    px[0]=0.0;
    px[x.size()-1]=0.0;
	x_min=*min_element(x.begin(),x.end());
	x_max=*max_element(x.begin(),x.end());
}


Sampler::~Sampler()
{
}


void Sampler::normalize( )
{
	double norm=int_tabulated(x,px)	;
//	cout<<"Normalization factor "<<norm<<" "<<x.size()<<endl;
//	norm=1.0;
	for(unsigned int i=0;i<x.size();++i)
	px[i]/=norm;
	meanx=0.0;
	vector<double> temp(x.size(),0.0);
	for(unsigned int i=0;i<x.size();++i)
		temp[i]=px[i]*x[i];
	meanx=int_tabulated(x,temp);
	x_min=*min_element(x.begin(),x.end());
	x_max=*max_element(x.begin(),x.end());
	calculateCpd();
	//	cout<<"Mean x "<<meanx<<endl;
}




void Sampler::setSeed(int64_t seed)
{
	srand48(seed)	;
}

void Sampler::setRange(double xmin1,double xmax1)
{
	unsigned int	up=0;

	//swap arguments if they were input backward
	if(xmax1 < xmin1)
	{
		double temp=xmin1;
		xmin1=xmax1;
		xmax1=temp;
	}

	//parse the upper limit of the requested range
	if(xmax1>=x_max) //above upper end of full range
	{
		cpd_max=cpd.back();
	}
	else if(xmax1<=x_min) //below lower end of full range
	{
		cpd_max=cpd.front();
	}
	else
	{
		up=(upper_bound(x.begin(),x.end(),xmax1)-x.begin())-1; //find closest tabulated point to upper end


		cpd_max=cpd[up]+(cpd[up+1]-cpd[up])*(xmax1-x[up])/(x[up+1]-x[up]); //interpolate upper end 
	}


	// parse lower limit of requested range
	if(xmin1>=x_max)
	{
		cpd_min=cpd.back();
	}
	else	if(xmin1<=x_min)
	{
		cpd_min=cpd.front();
	}
	else
	{
		up=(upper_bound(x.begin(),x.end(),xmin1)-x.begin())-1;
		cpd_min=cpd[up]+(cpd[up+1]-cpd[up])*(xmin1-x[up])/(x[up+1]-x[up]);
	}

	// changed Jul 28 ,2009
	// resets the range?
	if(xmax1<=x_min)
		cpd_min=cpd[1];
	if(xmin1>=x_max)
		cpd_max=cpd[cpd.size()-2];

	//check that the requested range isn't zero
	if(cpd_min==cpd_max)
	{
		cout<<"Warning:Range may be too small "<<xmin1<<" "<<xmax1<<" "<<x.front()<<" "<<x.back()<<endl;
		exit(1);
	}

}

void Sampler::getFacv(vector<double> &x1_a,vector<double> &x2_a,vector<double> &fac_a)
{
	fac_a.resize(x1_a.size());
	for(unsigned int j=0;j<(x1_a.size());++j)
		fac_a[j]=fabs(interpolate(x,cpd,x2_a[j])-interpolate(x,cpd,x1_a[j]));
}




double Sampler::getFac(double xmin1,double xmax1)
{
	return fabs(interpolate(x,cpd,xmin1)-interpolate(x,cpd,xmax1));
}


double Sampler::rand( )
{
	double u=cpd_min+grandomu()*(cpd_max-cpd_min);
	int up;
	if(cpd.front()==0.0)
		up=int(upper_bound(cpd.begin(),cpd.end(),u)-cpd.begin())-1;
	else
		up=int(upper_bound(cpd.begin(),cpd.end(),u,greater<double>())-cpd.begin())-1;
	//	cout<<cpd.front()<<" a "<<cpd.back()<<" "<<up<<" "<<cpd[up+1]<<" "<<cpd[up]<<endl;
	return x[up]+(x[up+1]-x[up])*(u-cpd[up])/(cpd[up+1]-cpd[up]);
}

vector<double> Sampler::randv(int nsize1)
{
	vector<double> y(nsize1);
	for(int i=0;i<nsize1;++i)
		y[i]=rand();
	return y;
}

void Sampler::calculateCpd( )
{
	cpd.resize(px.size());
	if(px.front()<px.back())
	{
		cpd.front()=0.0;
		double temp=0.0;
		for(unsigned int i=1;i<x.size();++i)
		{
			temp+=(x[i]-x[i-1])*(px[i]+px[i-1])*0.5;
			cpd[i]=temp;
//			cout<<i<<" "<<cpd[i]<<" "<<(cpd[i]-cpd[i-1])<<endl;
		}
	}
	else
	{
		cpd.back()=0.0;
		double temp=0.0;
		for(int i=x.size()-2;i>=0;--i)
		{
			temp+=(x[i+1]-x[i])*(px[i+1]+px[i])*0.5;
			cpd[i]=temp;
//			cout<<i<<" "<<cpd[i]<<" "<<(cpd[i]-cpd[i+1])<<endl;
		}
	}

	cpd_min=cpd.front();
	cpd_max=cpd.back();
}











