
#include "IsochroneDB.h"
#include "Functions.h"
#include <utility>
#include <iostream>
#include "Satellite.h"
#include "ebfvector.hpp"
#include"Parameters.h"

int main(int argc, char **argv) {

	double mmin, mmax, age, feh;

	Parameters All;
	All.setFromArguments(argc, argv);

	IsochroneDB ic(All.inputDir + "Isochrones/", "padova/", All.photoSys, All.magcolorNames, 1);

	cout<<"iage "<<"ifeh "<<"Mmin "<<"Mmax "<<"age "<<"feh"<<endl;
	for(age=6.5;age<10.2;age+=0.5) {
		for(feh=-2.25;feh<0.2;feh+=0.5) {
			ic.min_max_m_new(age,feh,0.0,&mmin,&mmax,11.0,0);
			cout<<age<<" "<<feh<<endl;
		}
	}
	return 0;
}