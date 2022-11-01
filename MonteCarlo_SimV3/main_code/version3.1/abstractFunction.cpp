#include "crater.h"
#include "abstractFunction.h"

using namespace std;


double abstractFunction::getXfromY(double y){
	return rootSolver(y);
}


double abstractFunction::rootSolver(double x){
	// use getYfromX as the function to solve
}



VvsX::abstractFunction(crater* ctr){
	initializeAbstractFunction(ctr);
}

void VvsX::initializeAbstractFunction(crater* ctr){
	n2     = (ctr->get_regime() == ctr->strengthRegime ? currentTarget->get_n2s() : currentTarget->get_n2g());
	n1     = currentTarget->get_n1();
	k      = currentTarget->get_k();
	C1     = currentTarget->get_C1();
	mu     = currentTarget->get_mu();
	nu     = currentTarget->get_nu();
	p      = currentTarget->get_p();
	rho    = currentTarget->get_density_rho();
	a      = currentImpactor->get_radius_a();
	m      = currentImpactor->get_mass_m();
	delta  = currentImpactor->get_density_delta();
	Unorm = computeUnorm(currentImpactor->get_speed_U(), currentImpactor->get_elevationAngle_alpha());
}