#ifndef LUNAREJECTA_ASSEMBLY_H
#define LUNAREJECTA_ASSEMBLY_H

#include "lunarEjecta_FractalIntegration.h"
#include "lunarEjecta_AdaptiveMesh.h"
#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_Regolith.h"
#include "lunarEjecta_SecondaryFluxData.h"
#include "lunarEjecta_MeteoroidFlux.h"
#include "lunarEjecta_NearEarthObjectFlux.h"

//#include "lunarEjecta_SpeedZenithIntegration.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>

using namespace std;

// Note: a template class must be definied in the header file
//  There is another solution...
//  https://stackoverflow.com/questions/495021/why-can-templates-only-be-implemented-in-the-header-file?rq=1
template <class genMEMdataHi, class genMEMdataLo, class genOutput>
class lunarEjecta_Assembly
{
public:
	lunarEjecta_Assembly(
		/*  For lunarEjecta_Regolith */
		int HH11_targetMaterial,
		int regolithDensType,
		double new_lowDensity,
		double new_avgDensity,
		double new_highDensity,
	    double new_escapeSpeed,

		/*  For ImpactSites_and_ROI */
		double new_ND,     // total number of distance increments
        double new_Nazm,   // total number of azimuth increments
        double new_radius, // radius of Moon
        latLon& new_ROI,   // lat-lon location of Region-Of-Interest

		/* For MEM_LatData */
		string dn,   // directory name of lat data
		double lMin, // Minimum lat (most likely -90)
		double lMax, // Maximum lat (most likely +90)
		int NL,      // Assumes lat data equally spaced with NL # of lat files

		/* For lunarEjecta_NearEarthObjectFlux */
		string NEOfn,
		double new_m_min,
		double new_m_max,
		int densType,
		double userDefDens,
		string dn_NEO,
		double lMin_NEO,
		double lMax_NEO,
		int NL_NEO,

		/* For SecondaryFluxData */
		string fn, // file name
		double new_xMin, // min of x-axis of integral flux
		double new_xMax, // max of x-axis of integral flux
		int new_Nx,          // number of spacings on x-axis
		int new_xScale,      // xScaleType = linear or log10
		int new_NSetsXY,           // number of sets of x-y data, if 0 will ignore setMin and setMax
		vector<double> new_setMin, // minimum of set range i
		vector<double> new_setMax, // maximum of set range i

		/* For lunarEjecta_AdaptiveMesh */
		//int new_Nx_mesh,  // number of x cell edges
		//int new_Ny_mesh,  // number of y cell edges
		int new_maxLevelMesh,    // the division level of the integration mesh
		int new_maxLevelFractal) // the division level of the integrand-domain probing
	{
		cout << "----------------------------------------\n";
		cout << "lunarEjecta_Assembly template class init \n";
		cout << "----------------------------------------\n";

		// init regolith
		RegolithProperties = new lunarEjecta_Regolith(HH11_targetMaterial, regolithDensType, new_lowDensity, new_avgDensity, new_highDensity, new_radius, new_escapeSpeed);

		// init site and ROI locations
		ImpactSitesROILoc = new ImpactSites_and_ROI(new_ND, new_Nazm, new_radius, new_ROI);

		// establish high density and low density MEM input data
		MEMLatDataHi = new MEM_LatData<genMEMdataHi>(dn, lMin, lMax, NL);
		MEMLatDataLo = new MEM_LatData<genMEMdataLo>(dn, lMin, lMax, NL);

		// establish secondary flux data
		SecFluxOutputData = new SecondaryFluxData<genOutput>(fn, new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax);
	
		NEOLatData = new lunarEjecta_NearEarthObjectFlux(NEOfn, new_m_min, new_m_max, densType, userDefDens, dn_NEO, lMin_NEO, lMax_NEO, NL_NEO);

		cout << "H_compH11GrunMassFactor = \n" << this->H_compH11GrunMassFactor(100)*1E6 << " E-6" << endl;

		cout << "H_compH11HiDensFactor = \n" << this->H_compH11HiDensFactor() << endl;
		cout << "H_compH11LoDensFactor = \n" << this->H_compH11LoDensFactor() << endl << endl;

		cout << "H_C03RegolithParticleDiameterCDF = (> 1 cm) \n" << this->H_C03RegolithParticleDiameterCDF(1.E-2) << endl;
		cout << "H_C03RegolithParticleDiameterCDF = (> 0.1 mm) \n" << this->H_C03RegolithParticleDiameterCDF(1.E-4) << endl;
		cout << "H_C03RegolithParticleDiameterCDF = (> 0.1 um) \n" << this->H_C03RegolithParticleDiameterCDF(1.E-7) << endl << endl;


		cout << "H_C03RegolithParticleMassCDF = (> 10 g, low dens) \n" << this->H_C03RegolithParticleMassCDF(1.E-2, 0) << endl;
		cout << "H_C03RegolithParticleMassCDF = (> 10 g, high dens) \n" << this->H_C03RegolithParticleMassCDF(1.E-2, 1) << endl << endl;

		cout << "H_C03RegolithParticleMassCDF = (> 1 ug, low dens) \n" << this->H_C03RegolithParticleMassCDF(1.E-9, 0) << endl;
		cout << "H_C03RegolithParticleMassCDF = (> 1 ug, high dens) \n" << this->H_C03RegolithParticleMassCDF(1.E-9, 1) << endl << endl;

		cout << "H_C03RegolithParticleMassCDF = (> 1 ng, low dens) \n" << this->H_C03RegolithParticleMassCDF(1.E-12, 0) << endl;
		cout << "H_C03RegolithParticleMassCDF = (> 1 ng, high dens) \n" << this->H_C03RegolithParticleMassCDF(1.E-12, 1) << endl << endl;
	
		cout << "H_C03RegolithParticleMassCDF = (> 1 pg, low dens) \n" << this->H_C03RegolithParticleMassCDF(1.E-15, 0) << endl;
		cout << "H_C03RegolithParticleMassCDF = (> 1 pg, high dens) \n" << this->H_C03RegolithParticleMassCDF(1.E-15, 1) << endl << endl;
	
		this->H_init_normalization();

		AdaptiveMesh = new lunarEjecta_AdaptiveMesh(SecFluxOutputData->getNalt(), SecFluxOutputData->getNvel(), new_maxLevelMesh, new_maxLevelFractal);
		
	}

	~lunarEjecta_Assembly() {
		delete RegolithProperties;
		delete ImpactSitesROILoc;
		delete MEMLatDataHi;
		delete MEMLatDataLo;
		delete SecFluxOutputData;
		delete AdaptiveMesh;
		delete NEOLatData;
	}

	void computeSecondaryFlux() { // All the magic happens here!
		double runtime, percentFinished;

		int i_siteDist, j_siteAzm, k_impactHorzAngle, l_impactAzm, m_impactSpeed, count = 1, count2 = 1;
		int ii_ejectaAzm, jj_ejectaAlt, kk_ejectaSpeed;
		int Ni, Nj, Nk, Nl, Nm;
		int Nii, Njj, Nkk;
		double D0, D1; // units of circumference
		double siteLat; // in degrees, used for getting the right impact latitude MEM file
		double siteSA; // the site surface area, in m^2
		double incomingAzm_at_ROI;  // units of rads
		double outgoingAzm_at_site; // units of rads
		double impactHorzAngle; // units of rads
		double impactAzm; // units of rads
		double impactSpeed; // units of km/s
		double Dbeta; // swath of azm the secondaries reach the ROI, units of rads

		double MEM_fluxLo;   // #/m^2/yr
		double MEM_fluxHi;   // #/m^2/yr
		double NEO_massflux; // kg/m^2/yr (already takes into acount integrating over the HH11 mass term wrt the NEO mass spectrum)

		cout << "------------------------------------\n";
		cout << "---- Computing Secondary Fluxes ----\n";
		cout << "------------------------------------\n";

		Ni = ImpactSitesROILoc->getND();
		Nj = ImpactSitesROILoc->getNazm();
		Nk = MEMLatDataHi->getNphi()/2; // ignoring below the horizon
		Nl = MEMLatDataHi->getNtheta();
		Nm = MEMLatDataHi->getNvel();

		Nii = SecFluxOutputData->getNazm();
		Njj = SecFluxOutputData->getNalt(); // the same as our integration mesh
		Nkk = SecFluxOutputData->getNvel(); // ''

		cout << " Ni (i_siteDist) = " << Ni << endl;
		cout << " Nj (j_siteAzm) = " << Nj << endl;
		cout << " Nk (k_impactHorzAngle) = " << Nk << endl;
		cout << " Nl (l_impactAzm) = " << Nl << endl;
		cout << " Nm (m_impactSpeed) = " << Nm << endl;
		cout << " Nii (ii_ejectaAzm) = " << Nii << endl;
		cout << " Njj (jj_ejectaAlt) = " << Njj << endl;
		cout << " Nkk (kk_ejectaSpeed) = " << Nkk << endl;

		auto t1 = std::chrono::high_resolution_clock::now();

		// For each location, outer loop distance, inner loop azm
		// Defines D0 and D1, given a ROI
		for (i_siteDist = 0; i_siteDist < Ni; ++i_siteDist)
		{
			// all distance dependent only terms should be computed here
			// units of circumference (so all distances range from 0 to 1)
			D0 = (ImpactSitesROILoc->getD(i_siteDist) - ImpactSitesROILoc-> getROI_radius()/ ImpactSitesROILoc->getradius()) / (2.*PI); 
			D1 = (ImpactSitesROILoc->getD(i_siteDist) + ImpactSitesROILoc-> getROI_radius()/ ImpactSitesROILoc->getradius()) / (2.*PI); 

			Dbeta = ImpactSitesROILoc->getDbeta(D0, D1);

			siteSA = ImpactSitesROILoc->getsite_SA(i_siteDist); // units m^2

			for (j_siteAzm = 0; j_siteAzm < Nj; ++j_siteAzm)
			{
				// units of rads
				/// used for output binning
				incomingAzm_at_ROI  = ImpactSitesROILoc->getsiteAzm(j_siteAzm, i_siteDist);
				// used as beta (in x = beta - beta_i, where beta_i is the primary flux azimuth)
				outgoingAzm_at_site = ImpactSitesROILoc->getROIAzm(j_siteAzm);

				siteLat = ImpactSitesROILoc->getsiteLatDeg(j_siteAzm, i_siteDist);

				cout << " Site Latitude = " << siteLat << endl;
				cout << "  Lat idx = " << MEMLatDataHi->getLatIdx(siteLat) << endl;

				// timing stuff
				// see: https://stackoverflow.com/questions/12231166/timing-algorithm-clock-vs-time-in-c
				auto t2 = std::chrono::high_resolution_clock::now();
				runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000. / 60.;
				percentFinished = count / double(Ni*Nj*Nk*Nl);

				//cout << "*****************************************************\n";
				cout << " Percent finished = " << 100.*percentFinished << endl;
				//cout << "*****************************************************\n";
				cout << "  run time = " << runtime << " minutes |  time remaining = " << runtime * (1./percentFinished - 1.) << endl;

				// Next, we will loop over the impact angle, impact azm, and impact speeds
				// for the MEM low and high density populations
				// We will be ignoring negative horizon angles that MEM has
				//
				// The x-v integration depends on the impact angle and impact azm

				
				for (k_impactHorzAngle = 0; k_impactHorzAngle < Nk; ++k_impactHorzAngle)
				{
					impactHorzAngle = PI/2. * k_impactHorzAngle / double(Nk);

					for (l_impactAzm = 0; l_impactAzm < Nl; ++l_impactAzm)
					{
						impactAzm = 2.*PI * l_impactAzm / double(Nl);

						// First, compute the integration grid (which gives the secondary
						// flux speed and horizon bins, the azm is given by the incomingAzm_at_ROI)
						AdaptiveMesh->restartBins();
						AdaptiveMesh->evalBins(D0, D1, outgoingAzm_at_site - impactAzm, Dbeta, RegolithProperties->getHH11_mu(), PI/2. - impactHorzAngle);

						count++;
						
						// The normalization terms precompute the impact speed depence,
						// so all we need to do here is loop over the normalization terms
						// to finish the computation
						// The norm terms also integrate out the impactor density distribution
						// as well as the impactor mass distribution, assuming a single density
						// for the regolith
						for (m_impactSpeed = 0; m_impactSpeed < Nm; ++m_impactSpeed)
						{ // speed in km/s
							impactSpeed = MEMLatDataHi->getvMin() + (MEMLatDataHi->getvMax() - MEMLatDataHi->getvMin()) * (m_impactSpeed+0.5) / double(Nm);
							//cout << impactspeed << endl;

							// at this point, we can get specific MEM fluxes for high and low densities, and NEO's
							// units of (#/yr)
							MEM_fluxLo = siteSA * MEMLatDataLo->getFlux_atAngleVelLat(impactHorzAngle/DtoR, impactAzm/DtoR, impactSpeed, siteLat);
							MEM_fluxHi = siteSA * MEMLatDataHi->getFlux_atAngleVelLat(impactHorzAngle/DtoR, impactAzm/DtoR, impactSpeed, siteLat);
							NEO_massflux   = siteSA * NEOLatData->getMassFluxNEO_atAngleVelLat(impactHorzAngle/DtoR, impactAzm/DtoR, impactSpeed, siteLat);

							cout << " Alt, Azm, Vel = " << impactHorzAngle/DtoR << ' ' << impactAzm/DtoR << ' ';
							cout << impactSpeed << " | Lo, Hi, and NEO fluxes = " << scientific << MEM_fluxLo << ", " << MEM_fluxHi << ", " << NEO_massflux << endl;
							
							// Next, we loop over the secondary ejecta bins for 
							// ejecta azimuth, altitude, and speed
							for (ii_ejectaAzm = 0; ii_ejectaAzm < Nii; ++ii_ejectaAzm)
							{
								for (jj_ejectaAlt = 0; jj_ejectaAlt < Njj; ++jj_ejectaAlt)
								{
									for (kk_ejectaSpeed = 0; kk_ejectaSpeed < Nkk; ++kk_ejectaSpeed)
									{
										count2++;
									}
								}
							}


						} // END FOR, impact speed
					} // END FOR, impact azimuth angle
				} // END FOR, impact horizon angle
				
			}// END FOR, impact outgoing secondary azimuth
		} // END FOR, impact distance
		cout << " Total count = " << count << endl;
		cout << " Total count2 = " << count2 << endl;
	}



// once these are initialized, they shouldn't need to change since we will always have
	//  the same resolution throughout the simulation
	vector<double> x; // x = 1 - cos(zenith), for integration
	vector<double> y; // y = v/v_esc, for integration
	vector<vector<double>> z; // z[Nx][Ny]


private:

	// Modeled fit to mean of Fig. 1 of Carrier 2003 using bi-exponential function, useful range 0.001 mm to 10 mm
	// For size-limiting of secondary ejecta
	double H_C03RegolithParticleDiameterCDF(double d_min) { // units m, outputs % greater than d_min
		const double a =  0.05480447;
		const double b = -1.01472478;
		const double c =  0.33749929;
		const double d = -0.25180816;
		d_min *= 1000.; // convert from m to mm

		return exp(-1. / (a*pow(d_min, b) + c*pow(d_min, d)));
	}

	// For mass-limiting of secondary ejecta, based on Carrier 2003 regolith size distribution
	double H_C03RegolithParticleMassCDF(double m_min, int lowHighDens) { // units of kg
		double dens = 0.0, d_min;
		switch(lowHighDens){
			case 0: // low density regolith
				dens = RegolithProperties->getlowDensity();
				break;
			case 1: // high density regolith
				dens = RegolithProperties->gethighDensity();
				break;
			case 2: // average density regolith
				dens = RegolithProperties->getavgDensity();
				break;
			default:
				cerr << "ERROR: H_C03RegolithParticleMassCDF invalid density type selection\n";
		}

		d_min = 2. * pow(3.*m_min / (4.*PI*dens), 1./3.);

		return H_C03RegolithParticleDiameterCDF(d_min);
	}

	// eta defaults to 0 (to match HH11), but really should be 1
	// speed in units of km/s, alt in units rad
	double H_compH11ProjSpeedFactor(double speed, double alt, double eta = 0) {
		double sin_term;
		if(eta < 0.01) {
			sin_term = 1.;
		} else {
			sin_term = pow(sin(alt), eta); 
		}

		return pow(speed * sin_term, 3.*RegolithProperties->getHH11_mu());
	}

	inline double HH_Grun(double m) { // m = units of g, output in units of 1/(m^2*yr)
		return pow(2.2E3 * pow(m, 0.306) + 15., -4.38)
		     + 1.3E-9 * pow(m + 1.E11 * pow(m, 2.) + 1.E27 * pow(m, 4.), -0.36)
		     + 1.3E-16 * pow(m + 1.E6 * pow(m, 2.), -0.85); 
	}

	// fit made in SciDAVis, already normalized to 10^-6 g
	inline double HH_MDGrun(double m) {
		return 1. / (1.587E7 * pow(m, 2.3235-1.) + 7.056E4 * pow(m, 1.8189-1.));
	}


	double H_compH11GrunMassFactor(int N = 14) { // N ~100 is usually sufficient
		int i;
		double massMin, massMax = 10.; // mass limit in MEM is 10 g
		double Ai, Bi, sum = 0.;
		vector<double> m;
		vector<double> f;
		f.resize(N);

		massMin = MEMLatDataLo->getGrunMinMass(); // should be the same as Hi

		// prepare x and y vectors for power law fits
		logspace(m, massMin, massMax, N, 0, N);
		for (i = 0; i < N; ++i)
			f[i] = HH_MDGrun(m[i]) * (HH_Grun(massMin) / HH_Grun(1.E-6));

		// compute integral by breaking up in log spacings and approx by power law
		for (i = N-2; i >= 0; i--) // start at low end of mass to sum smaller #'s first (less floating-point error)
		{
			Bi = log10(f[i+1] / f[i]) / log10(m[i+1] / m[i]); 
			Ai = (f[i] + f[i+1]) / (pow(m[i], Bi) + pow(m[i+1], Bi));

			sum += Ai/(Bi+1.) * (pow(m[i+1], Bi+1) - pow(m[i], Bi+1));
			//cout << i << ' ' << Bi << endl;
		}
		return sum;	
	}

	// NEED TO CHANGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double H_compH11NEAMassFactor() {
		return 1.;
	}

	double H_compH11RegDensFactor(int lowHighDens){
		double dens = 0.0;
		switch(lowHighDens){
			case 0: // low density regolith
				dens = RegolithProperties->getlowDensity();
				break;
			case 1: // high density regolith
				dens = RegolithProperties->gethighDensity();
				break;
			case 2: // average density regolith
				dens = RegolithProperties->getavgDensity();
				break;
			default:
				cerr << "ERROR: H_compH11RegDensFactor invalid density type selection\n";
		}
		return pow(dens, -(3.*RegolithProperties->getHH11_nu()-1.));
	}

	// Note: checked against integrating in Mathematica, found out don't need the factor of 50 to make work.
	double H_compH11HiDensFactor() {
		double mass_sum = 0.0, cur_dens;

		for (int i = 0; i < MEMLatDataHi->getNdens(); ++i)
		{
			cur_dens = (MEMLatDataHi->getdensLEdge(i) + MEMLatDataHi->getdensREdge(i)) / 2.0; // kg/m^3
			mass_sum += MEMLatDataHi->getdensFraction(i) * pow(cur_dens, 3.*RegolithProperties->getHH11_nu()-1.) ; // bins are always 50 per kg/m^3 large
		}
		return mass_sum;
	}

	double H_compH11LoDensFactor() {
		double mass_sum = 0.0, cur_dens;

		for (int i = 0; i < MEMLatDataLo->getNdens(); ++i)
		{
			cur_dens = (MEMLatDataLo->getdensLEdge(i) + MEMLatDataLo->getdensREdge(i)) / 2.0; // kg/m^3
			mass_sum += MEMLatDataLo->getdensFraction(i) * pow(cur_dens, 3.*RegolithProperties->getHH11_nu()-1.) ; // bins are always 50 per kg/m^3 large
		}
		return mass_sum;
	}

	// NEED TO CHANGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double H_compH11NEADensFactor() {
		return 1.;
	}

	// beta - beta_i = pi
	inline double HH_zenithDownstream(double impactZenith) { // input and output in degrees
		return 20. + impactZenith * (1.5206 + impactZenith * (-0.036 + impactZenith * 0.0003));
	}

	// beta - beta_i = 0
	inline double HH_zenithUpstream(double impactZenith) { // input and output in degrees
		return 20. + impactZenith * (0.129 + impactZenith * (0.0236 + impactZenith * -0.00042));
	}

	// x = beta - beta_i
	// Note: beta_i is the impact azimuth (direction seen from impact point), which is defined to be due East (direction or Moon's rotation)
	//       
	//   beta is the secondary ejecta direction leaving the impact point (not direction seen from secondary impact point, ie ROI)
	//
	// impactZenith in degrees, x in radians, output in radians
	inline double HH_zenithGeneral(double impactZenith, double x) {
		return (HH_zenithDownstream(impactZenith)
			+ (HH_zenithUpstream(impactZenith) - HH_zenithDownstream(impactZenith)) / PI * fabs(x)) * DtoR;
	}

	// a = cos(a_max) / (1 - cos(amax)), converted to half sin to avoid subtraction of possible close #'s
	double a_power(double impactZenith, double x) {
		double alpha_max = HH_zenithGeneral(impactZenith, x); // in units of radians
		return cos(alpha_max) / (2. * sqr(sin(alpha_max / 2.)));
	}

	// exlusion zone around upstream direction wrt to beta_i (impact azm, ie East), in radians
	//  used in integration bounds of d(beta-beta_i) integral, impactZenith from rad to Degrees
	double exclusion_zone(double impactZenith) {
		double A = HH_zenithDownstream(impactZenith); // A and B are in degrees, but they cancel out their units
		double B = HH_zenithUpstream(impactZenith);

		//cout << A << ' ' << B << endl;

		if (B >= 0.)
		{
			return 0.;
		} else {
			return -PI * B / (A - B);
		} 
	}

	void H_init_normalization() {
		// assuming all MEM data has the same vel and angle resolution, separately
		/// First, we need to compute the integral term for each impact angle
		///  The impact angle implicitly controls the "exclusion zone", \Delta\beta
		cout << " Initializing normalization terms...\n";

		int i, j, idx;
		int Nalt = MEMLatDataHi->getNphi()/2; // ignoring below the horizon angles
		int Nvel = MEMLatDataHi->getNvel();

		double vMin = MEMLatDataHi->getvMin();
		double vMax = MEMLatDataHi->getvMax();

		cout << "Nalt = " << Nalt << " | Nvel = " << Nvel << endl;
		cout << "vMin = " << vMin << " | vMax = " << vMax << endl;

		double Gint_alt_i, vel, alt;
		double denseMEMLo, densMEMHi, densNEA;
		double massMEM, massNEA;
		double densRegolith, speedTerm, C4;

		normalizationHiDens.resize(Nalt * Nvel);
		normalizationLoDens.resize(Nalt * Nvel);
		normalizationNEA.resize(Nalt * Nvel);

		// compute density terms
		//// compute regolith dens term (currently assuming the same over whole Moon)
		densRegolith = H_compH11RegDensFactor(2/* Average Density */);
		denseMEMLo   = H_compH11LoDensFactor();
		densMEMHi    = H_compH11HiDensFactor();
		densNEA      = H_compH11NEADensFactor(); // NEED TO FINISH FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		cout << " densRegolith = " << densRegolith << endl;
		cout << " denseMEMLo   = " << denseMEMLo << endl;
		cout << " densMEMHi    = " << densMEMHi << endl;
		cout << " densNEA      = " << densNEA << endl;

		// compute mass term, not dependent on speed or impact angle
		massMEM = H_compH11GrunMassFactor(200); // units of kg
		massNEA = H_compH11NEAMassFactor(); // NEED TO FINISH FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

		cout << " massMEM = " << massMEM << endl;
		cout << " massNEA = " << massNEA << endl;

		C4 = RegolithProperties->getHH11_C4();

		cout << " C4 = " << C4 << endl;

		for (i = 0; i < Nalt; ++i)
		{
			alt = PI/2. * double(i) / double(Nalt-1.); // rad
			// compute integral G term
			Gint_alt_i = HH_compGint(alt);
			//cout << " i = " << i << " | Gint_alt_i = " << Gint_alt_i << endl;
			cout << alt << ' ' << Gint_alt_i << endl;

			for (j = 0; j < Nvel; ++j) 
			{
				idx = j + i * Nvel;
				vel = vMin + (vMax - vMin) * double(j + 0.5) / double(Nvel);

				speedTerm = H_compH11ProjSpeedFactor(vel, alt, 1);

				// build norm terms (the velocity term cancles out from the HH11 scaling laws, by definition)
				normalizationHiDens[idx] = C4 * massMEM * speedTerm * densRegolith * densMEMHi  / Gint_alt_i;
				normalizationLoDens[idx] = C4 * massMEM * speedTerm * densRegolith * denseMEMLo / Gint_alt_i;
				normalizationNEA[idx]    = C4 * massNEA * speedTerm * densRegolith * densNEA    / Gint_alt_i; // place holder #'s for now, need to finish the functions~~!'

				// cout << "alt = " << alt << " | vel = " << vel << " | Hi, Lo, NEA: ";
				// cout << scientific << normalizationHiDens[idx] << ' ' << normalizationLoDens[idx] << ' ';
				// cout << scientific << normalizationNEA[idx] << endl;
			}
		}
	}

	// derived in Mathematica, checked with plotting
	double HHH_BetaApproxLarge(double a) {
	/* order 1/a^1 */ return (1.
	/* order 1/a^2 */ - (1. + EulerGamma - log(a) 
	/* order 1/a^3 */ + (-6. + 6.*EulerGamma*(2. + EulerGamma) + sqr(PI) + 6.*log(a)*(2. + 2.*EulerGamma + log(a)) ) / (12.*a)) / a) / a;
	}

	double HHH_beta(double a) {
		// helps to avoid overflow errors doing it this way
		return exp(lgamma(1./a + 1.) + lgamma(a + 1.) - lgamma(a + 1./a + 2.));
	}

	// altitude distribution integral
	// x = \beta - \beta_i
	double HH_AltInt(double zenith, double x) {

		double a = a_power(zenith, x);

		if(a > 15) // large 'a', worst 0.46% error
		{
			return HHH_BetaApproxLarge(a);

		} else if(a < 1./15.) { // small 'a' worst 0.46% error

			return HHH_BetaApproxLarge(1./a);

		} else { // intermediate 'a'

			return HHH_beta(a);
		}
	}

	// zenith in units of rad, x = \beta - \beta_i
	double HH_AzmDist(double zenith, double x) {


		if(zenith < PI/3.)
		{
			return (1. + cos(x) * 3.*zenith / (2.*PI - 3.*zenith)) / (2.*PI);

		} else { // zenith > PI/3

			double b = (0.05 - 1.) /(PI/2. - PI/3.) * (zenith - PI/3.) + 1.;

			return exp(-(x - 2.*(zenith - PI/3.) ) / (PI*b)) + exp(-(x + 2.*(zenith - PI/3.) ) / (PI*b)); // not normalized to itself, it will get normalized with everything later
		}
	}


	// Integrand of cal(G), here, x = \beta - \beta_i
	double HH_GInt(double alt, double x) {
		// both functions need zenith angle and not altitude
		return HH_AzmDist(PI/2. - alt, x) * HH_AltInt(PI/2. - alt, x);
	}

	// using Romberg integration, see Appendix C of my dissertation
	double HH_compGint(double alt) {

		double delta_beta = exclusion_zone((PI/2. - alt)/DtoR /* zenith angle in deg */); // units of rad
		//cout << (PI/2. - alt)/DtoR << ' ' << delta_beta/DtoR << endl;
		double a = 0.; // left integration bounds
		double b = PI - delta_beta; // right integration bounds
		double eps = 1.E-3; // epsilon error in numerical integration

		double pow4m, hn, err = 10.*eps, fsum;
		int n = 1, m, k;
		vector<double> R;

		// base case, R(0,0), trapezoid rule
		hn = (b-a) / 2.;
		R.push_back(hn * (HH_GInt(alt, a) + HH_GInt(alt, b)));

		while(n < 5 || err > eps) { // do at least 2^4=16 function evals

			R.push_back(0.); // extend R array
			fsum = 0.;

			// fill in func eval gaps for next level
			for (k = 1; k <= (1 << (n-1)); k++)
				fsum += HH_GInt(alt, a + (2.*k - 1.)*hn);

			R[tri_idx(n,0)] = R[tri_idx(n-1,0)] / 2. + hn * fsum;

			pow4m = 1.;

			// compute n-th row of m's
			for (m = 1; m <= n; m++)
			{
				R.push_back(0.); // extend R array
				pow4m *= 4.;
				R[tri_idx(n,m)] = R[tri_idx(n,m-1)]
				                + (R[tri_idx(n,m-1)] - R[tri_idx(n-1,m-1)])
				                / (pow4m - 1.);

			}

			// compute relative error
			err = fabs((R[tri_idx(n,n)] - R[tri_idx(n-1,n-1)]) / R[tri_idx(n,n)]);
			n++;
			hn /= 2.;
		}
		return 2. * R[tri_idx(n-1, n-1)]; // multiply by 2 since we cut the domain in half at zero (even function)
	}


	// void H_init_integrationGrid(int Nx, int Ny)
	// {
	// 	int i;

	// 	// init the x bin edges

	// }




	lunarEjecta_Regolith*               RegolithProperties;
	ImpactSites_and_ROI*                ImpactSitesROILoc;
	lunarEjecta_NearEarthObjectFlux*    NEOLatData;
	MEM_LatData<genMEMdataHi>*          MEMLatDataHi;
	MEM_LatData<genMEMdataLo>*          MEMLatDataLo;
	SecondaryFluxData<genOutput>*       SecFluxOutputData;
	lunarEjecta_AdaptiveMesh*           AdaptiveMesh;
	//lunarEjecta_SpeedZenithIntegration* SpeedZenithIntegration;


	// A function of impact angle and impact speed, after integrating out impactor density and mass
	//  for each low and high density populations in MEM
	vector<double> normalizationHiDens; // will be size Nphi[alt]/2 * (Nvel), from meteoriodFlux
	vector<double> normalizationLoDens; // will be size Nphi[alt]/2 * (Nvel), from meteoriodFlux
	vector<double> normalizationNEA;    // will be size Nphi[alt]/2 * (Nvel), from meteoriodFlux
};


#endif 