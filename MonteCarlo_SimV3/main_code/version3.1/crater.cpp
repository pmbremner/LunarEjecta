#include "crater.h"

using namespace std;


crater::crater(impactor* imp, target* tar){
	currentImpactor = imp;
	currentTarget   = tar;

	computeRadius();
	computeSaveVmin();
	computeVmax();
	computeMtot();
}




double crater::getRadius(){
	return radius;
}

double crater::getVmin(){
	return vmin;
}

double crater::getVmax(){
	return vmax;
}

double crater::getMtot(){
	return Mtot;
}

double crater::getRegime(){
	return regime;
}
