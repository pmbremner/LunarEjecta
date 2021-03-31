#ifndef LUNAREJECTA_SECONDARYEJECTA_H
#define LUNAREJECTA_SECONDARYEJECTA_H

#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"
#include "LunarEjecta_scalinglaw.h"
#include "LunarEjecta_bearingdistmap.h"

using namespace std;



iglooSet* compute_ejecta(input *p, vector<iglooSet*> &primaryFluxes, scalingLaw *ejectaFactors, hist3DSet* bearingDistMap, int ejectaDistType);


#endif 