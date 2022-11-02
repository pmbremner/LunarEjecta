#include "crater.h"

using namespace std;


crater::crater(impactor* imp, target* tar){
	reinitializeCrater(impactor* imp, target* tar);
}

void crater::reinitializeCrater(impactor* imp, target* tar) {
	currentImpactor = imp;
	currentTarget   = tar;

	computeRadius();
	computeDiameterToDepthRatio();
	computeVmin();
	computeVmax();
	computeMtot();
}



double crater::get_radius_R(){
	return radius_R;
}

double crater::get_diameterToDepthRatio_K(){
	return diameterToDepthRatio_K;
}

double crater::get_vmin(){
	return vmin;
}

double crater::get_vmax(){
	return vmax;
}

double crater::get_Mtot(){
	return Mtot;
}

double crater::get_regime(){
	return regime;
}

double crater::computeMfromV(double v){

}

double crater::computedMdVfromV(double v){

}


void crater::h_computeRegime(double radiusS, double radiusG){
	if (radiusS < radiusG)
	{
		regime = strengthRegime;
	}
	else
	{
		regime = gravityRegime;
	}
}


void crater::computeRadiusAndRegime(){
	double radiusStrength = h_computeRadiusStrength();
	double radiusGravity  = h_computeRadiusGravity();

	h_computeRegime(radiusStrength, radiusGravity);

	if (regime == strengthRegime){
		radius_R = radiusStrength;
	}
	else if (regime == gravityRegime){
		radius_R = radiusGravity;
	}
}

// Note, in general, can depend on impactor properties
void crater::computeDiameterToDepthRatio()
{
	return 5.;
}

// This is the minimum speed at which a particle can leave the bottom
// of the crater and reach the crater rim
void crater::computeVmin(){
	vmin = 2. * sqrt(currentTarget->get_lunarGravity_g() * radius_R / diameterToDepthRatio_K);
}

void crater::computeVmax(){
	vmin = h_speedVsPositionInCrater(currentTarget->get_n1() * currentImpactor->getRadius());
}

void crater::computeMtot(){
	double n2 = (regime == strengthRegime ? currentTarget->get_n2s() : currentTarget->get_n2g());
	Mtot = h_ejectedMassVsPositionInCrater(n2 * radius_R);
}

// There are two ways to do this, O'Keefe & Ahrens 1985, or Koschny & Grun 2001
// The first gives larger fragments up to a total ejected mass of 0.1 kg
void crater::computeMmaxFragment(){
	MmaxFragment_mb = h_MmaxFragmentOkeefeAhrens1985();
}

double crater::h_rootSolver();

double crater::h_computeRadiusStrength(){
	double H2     = currentTarget->get_H2();
	double mu     = currentTarget->get_mu();
	double nu     = currentTarget->get_nu();
	double rho    = currentTarget->get_density_rho();
	double Y      = currentTarget->get_strength_Y();
	double m      = currentImpactor->get_mass_m();
	double delta  = currentImpactor->get_density_delta();
	double U      = currentImpactor->get_speed_U();
	double alpha  = currentImpactor->get_elevationAngle_alpha();

	double Unorm = computeUnorm(U, alpha);

	return H2 * pow(rho / m, -1./3.)
		* pow(rho / delta, (1. - 3.*nu) / 3.)
		* pow(Y / (rho * sqr(Unorm)), -mu/2.);
}


double crater::h_computeRadiusGravity(){
	double H1     = currentTarget->get_H1();
	double mu     = currentTarget->get_mu();
	double nu     = currentTarget->get_nu();
	double rho    = currentTarget->get_density_rho();
	double g      = currentTarget->get_lunarGravity_g();
	double a      = currentImpactor->get_radius_a();
	double m      = currentImpactor->get_mass_m();
	double delta  = currentImpactor->get_density_delta();
	double U      = currentImpactor->get_speed_U();
	double alpha  = currentImpactor->get_elevationAngle_alpha();

	double Unorm = computeUnorm(U, alpha);

	return H1 * pow(rho / m, -1./3.)
		* pow(rho / delta, (2. + mu - 6.*nu) / (3.*(2. + mu)))
		* pow(g * a / sqr(Unorm), -mu / (2. + mu));
}

double crater::computeUnorm(double U, double alpha){
	return U * cos(alpha);
}

double crater::h_speedVsPositionInCrater(double x){
	double n2     = (regime == strengthRegime ? currentTarget->get_n2s() : currentTarget->get_n2g());
	double n1     = currentTarget->get_n1();
	double C1     = currentTarget->get_C1();
	double mu     = currentTarget->get_mu();
	double nu     = currentTarget->get_nu();
	double p      = currentTarget->get_p();
	double rho    = currentTarget->get_density_rho();
	double a      = currentImpactor->get_radius_a();
	double delta  = currentImpactor->get_density_delta();
	double U      = currentImpactor->get_speed_U();
	double alpha  = currentImpactor->get_elevationAngle_alpha();

	double Unorm = computeUnorm(U, alpha);

	x = h_forceXBounds(n1*a, n2*radius_R, x);

	return Unorm * C1
		* pow((x / a) * pow(rho / delta, nu), -1./mu)
		* pow(1. - x / (n2*radius_R), p);
}

double crater::h_ejectedMassVsPositionInCrater(double x){
	double n2     = (regime == strengthRegime ? currentTarget->get_n2s() : currentTarget->get_n2g());
	double n1     = currentTarget->get_n1();
	double k      = currentTarget->get_k();
	double rho    = currentTarget->get_density_rho();
	double a      = currentImpactor->get_radius_a();
	double m      = currentImpactor->get_mass_m();
	double delta  = currentImpactor->get_density_delta();

	x = h_forceXBounds(n1*a, n2*radius_R, x);

	return m * (3.*k / (4.*PI)) * (rho / delta)
		* (pow(x / a, 3) - pow(n1, 3));
}

double crater::h_MmaxFragmentOkeefeAhrens1985(){
	return 0.2 * pow(Mtot * 1000./*kg to g*/, 0.8) / 1000. /*g to kg*/;
}
double crater::h_MmaxFragmentKoschnyGrun2001(){
	return 0.01 * Mtot;
}



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////




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
	Unorm  = computeUnorm(currentImpactor->get_speed_U(), currentImpactor->get_elevationAngle_alpha());
	R      = ctr->radius_R;

	xmin = n1*a;
	xmax = n2*R;

}