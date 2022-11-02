#ifndef CRATER_H
#define CRATER_H


using namespace std;


class crater
{
public:
	crater(impactor* imp, target* tar);
	~crater();

	void reinitializeCrater(impactor* imp, target* tar);

	double get_radius_R();
	double get_diameterToDepthRatio_K();
	double get_vmin();
	double get_vmax();
	double get_Mtot();
	double get_regime();

	double computeMfromV(double v);
	double computedMdVfromV(double v);
	double computeUnorm(double U, double alpha);


	enum regimeType {strengthRegime, gravityRegime};

private:
	void h_computeRegime(double radiusS, double radiusG);
	void computeRadiusAndRegime();
	void computeDiameterToDepthRatio();
	void computeVmin();
	void computeVmax();
	void computeMtot();
	void computeMmaxFragment();

	double h_rootSolver();
	double h_computeRadiusStrength();
	double h_computeRadiusGravity();
	
	double h_speedVsPositionInCrater(double x);
	double h_ejectedMassVsPositionInCrater(double x);

	double h_MmaxFragmentOkeefeAhrens1985();
	double h_MmaxFragmentKoschnyGrun2001();



	impactor* currentImpactor;
	target*   currentTarget;

	abstractFunction* ejectedSpeed;
	abstractFunction* ejectedMass;

	regimeType regime; // strengthRegime or gravityRegime
	double radius_R; // [m]
	double diameterToDepthRatio_K;
	double vmin;   // [m/s]
	double vmax;   // [m/s]
	double Mtot;   // [kg]
	double MmaxFragment_mb; // [kg]
	

	//const double lunarGravity = 1.625; // [m/s^2]
};


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


// pure virtual, need to use derived class
class abstractFunction
{
public:
	abstractFunction();
	~abstractFunction();

	virtual double getYfromX(double x) = 0;
	double getXfromY(double y);

protected:
	double rootSolver(double x);

	double xmin;
	double xmax;
	double ymin;
	double ymax;
};



class VvsX: public abstractFunction
{
public:
	abstractFunction(crater* ctr);
	~abstractFunction();

	void initializeAbstractFunction(crater* ctr);
	
	double getYfromX(double x);
	//double getXfromY(double y);

private:
	double n2;
	double n1;
	double k;
	double C1;
	double mu;
	double nu;
	double p;
	double rho;
	double a;
	double m;
	double delta;
	double Unorm;
	double R;

};


class MvsX: public abstractFunction
{
public:
	abstractFunction(crater* ctr);
	~abstractFunction();

	void initializeAbstractFunction(crater* ctr);
	
	double getYfromX(double x);
	//double getXfromY(double y);


private:
	double n2;
	double n1;
	double k;
	double C1;
	double mu;
	double nu;
	double p;
	double rho;
	double a;
	double m;
	double delta;
	double Unorm;
	double R;

};


#endif 