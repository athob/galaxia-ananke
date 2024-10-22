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


#ifndef STARPARTICLE_H_
#define STARPARTICLE_H_
#include "Functions.h"
#include<string>
#include<vector>


/*
class StarParticlef
{
public:
	float Pos[29];
    inline float& feh() {return Pos[6];}
    inline float& alpha() {return Pos[7];}
	inline float& age() {return Pos[8];}
    inline float& smass() {return Pos[9];}
    inline float& rad() {return Pos[10];}
    inline float& mag(int k) {return Pos[11+k];}
   inline int& satID() {return *(int *)(Pos+14);}
   inline int& popID() {return *(int *)(Pos+15);}
   inline int& partID() {return *(int *)(Pos+16);}
   inline int& parentID() {return *(int *)(Pos+17);}
//added formation distance
   inline float& dform() {return Pos[18];}
//added extra abundances
    inline float& helium() {return Pos[19];}
    inline float& carbon() {return Pos[20];}
    inline float& nitrogen() {return Pos[21];}
    inline float& oxygen() {return Pos[22];}
    inline float& neon() {return Pos[23];}
    inline float& magnesium() {return Pos[24];}
    inline float& silicon() {return Pos[25];}
    inline float& sulphur()  {return Pos[26];}
    inline float& calcium() {return Pos[27];}
    inline float& dmod() {return Pos[28];}
};
*/


class StarParticle
{
private:
	double Pos[33];
public:
    inline double& pos(int i) {return Pos[i];}
    inline double& feh() {return Pos[6];}
    inline double& alpha() {return Pos[7];}
	inline double& age() {return Pos[8];}
    inline double& smass() {return Pos[9];}
    inline double& rad() {return Pos[10];}
    inline double& mag(int k) {return Pos[11+k];}
    inline double& lum() {return Pos[14];}
    inline double& teff() {return Pos[15];}
    inline double& grav() {return Pos[16];}
    inline double& dmod() {return Pos[17];}
    inline int64_t& satID() {return *(int64_t *)(Pos+18);}
    inline int64_t& popID() {return *(int64_t *)(Pos+19);}
    inline int64_t& partID() {return *(int64_t *)(Pos+20);}
    inline int64_t& parentID() {return *(int64_t *)(Pos+21);}
    inline int64_t& partitionID() {return *(int64_t *)(Pos+22);}
//added formation distance
    inline double& dform() {return Pos[23];}
//added extra abundances
    inline double& helium() {return Pos[24];}
    inline double& carbon() {return Pos[25];}
    inline double& nitrogen() {return Pos[26];}
    inline double& oxygen() {return Pos[27];}
    inline double& neon() {return Pos[28];}
    inline double& magnesium() {return Pos[29];}
    inline double& silicon() {return Pos[30];}
    inline double& sulphur()  {return Pos[31];}
    inline double& calcium() {return Pos[32];}
    
 /*   void convert(StarParticlef &Starf)
    {
    	for(int i=0;i<14;++i)
    		Starf.Pos[i]=float(Pos[i]);
        for(int i=18;i<32;++i)
            Starf.Pos[i]=float(Pos[i+4]);
        Starf.Pos[28]=float(Pos[17]);

    	Starf.satID()=int(satID());
    	Starf.popID()=int(popID());
    	Starf.partID()=int(partID());
        Starf.parentID()=int(parentID());

    }
*/
    void print()
    {
    	cout<<setw(15)<<Pos[0]<<setw(15)<<Pos[3]<<setw(15)<<feh()<<setw(15)<<age()<<setw(15)<<smass()<<" "<<mag(0)<<" "<<mag(1)<<" "<<mag(2)<<" "<<dmod()<<endl;
    }

    inline void AbsToAppMag()
    {
    	rad()=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]+Pos[2]*Pos[2]);
    	dmod()=5*log10(rad()*100.0);
    	mag(0)+=dmod();
    	mag(1)+=dmod();
    	mag(2)+=dmod();
    }

    inline void AddPhotoError(Normaldev &gauss,int photoError,int magswap)
    {
    	double sigma_mag=0.0;
    	if(photoError==2)
    	{
    		double temp=pow(10.0,0.4*(min(25.0,mag(0))-24.5));
    		sigma_mag=sqrt((0.04-0.039)*temp+0.039*temp*temp);
    	}
    	else if(photoError==3)
    	{
    		double temp=pow(10.0,0.4*(min(28.0,mag(0))-27.5));
    		sigma_mag=sqrt((0.04-0.039)*temp+0.039*temp*temp);
    	}
    	else if(photoError==4)
    	{
    		double temp=pow(10.0,0.4*(min(23.0,mag(0))-22.6));
    		sigma_mag=sqrt((0.04-0.039)*temp+0.039*temp*temp);
    	}

		if(sigma_mag>0.0)
		{
			mag(0)+=gauss.dev()*sigma_mag;
			mag(1)+=gauss.dev()*sigma_mag;
			mag(2)+=gauss.dev()*sigma_mag;
			if(magswap==2)
				mag(2)=mag(0);
			else if (magswap==1)
				mag(1)=mag(0);
		}

    }

};


class ParticleTag
{
public:
	ParticleTag(string name1,int dim1,int id1,int datatype1,int status1):name(name1),dim(dim1),id(id1),status(status1), datatype(datatype1){}
	string name;
	int dim;
	int id;
	int status;
	int datatype;
	bool operator== (const ParticleTag& a)
	{
		return (a.name==name);
	}
};

class ParticleTag_eq: public unary_function<ParticleTag,bool>
{
	string s;
	int status;
public:
	explicit ParticleTag_eq(const string& ss,const int &st): s(ss), status(st){}
	bool operator() (const ParticleTag &t) const {return ((t.name==s)&&(t.status==status));}
};



class ParticleStar
{
public:
//	ParticleStar();
//	~ParticleStar();
    double* Pos;       // particle position
    inline double& FeH() {return Pos[6];}
    inline double& Alpha() {return Pos[7];}
    inline double& Age() {return Pos[8];}
    inline double& Mass() {return Pos[9];}
    inline double& parentID() {return Pos[10];}
    inline double& partitionID() {return Pos[11];}
//added formation distance
    inline double& dform() {return Pos[12];}
//added extra abundances
    inline double& helium() {return Pos[13];}
    inline double& carbon() {return Pos[14];}
    inline double& nitrogen() {return Pos[15];}
    inline double& oxygen() {return Pos[16];}
    inline double& neon() {return Pos[17];}
    inline double& magnesium() {return Pos[18];}
    inline double& silicon() {return Pos[19];}
    inline double& sulphur()  {return Pos[20];}
    inline double& calcium() {return Pos[21];}
    inline double& h(int k) {return Pos[22+k];};    

//    inline double& Density() {return Pos[11];} //replaced in new version

     void getTags(vector<ParticleTag> &tags,int hdim1)
    {
    	int c=0;
    	int dtype=5;
    	tags.clear();
//    	tags.push_back(ParticleTag("/Pos",pos_size,c,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Pos3",3,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Vel3",3,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/FeH",1,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Alpha",1,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Age",1,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Mass",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/ParentID",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/PartitionID",1,c,dtype,1));c+=tags.back().dim;

//added formation distance
        tags.push_back(ParticleTag("/Dform",1,c,dtype,1));c+=tags.back().dim;

// added extra abundances
        tags.push_back(ParticleTag("/Helium",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/Carbon",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/Nitrogen",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/Oxygen",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/Neon",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/Magnesum",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/Silicon",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/Sulphur",1,c,dtype,1));c+=tags.back().dim;
        tags.push_back(ParticleTag("/Calcium",1,c,dtype,1));c+=tags.back().dim;


    	if(hdim1==3)
			{tags.push_back(ParticleTag("/H_smooth",2,c,dtype,1));c+=tags.back().dim;}
    	else if(hdim1==6)
    		{tags.push_back(ParticleTag("/H_cubic",2,c,dtype,1));c+=tags.back().dim;}
    	else
    	{cout<<"hdim must be 3 or 6"<<endl; exit(1);}




 //   	tags.push_back(ParticleTag("/Density",1,c,0));c+=tags.back().dim; //replaced in new version
    }

};


#endif /* PARTICLE2_H_ */
