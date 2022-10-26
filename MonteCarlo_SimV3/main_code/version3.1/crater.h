#ifndef CRATER_H
#define CRATER_H


using namespace std;


class crater
{
public:
	crater(impactor* imp, target* tar);
	~crater();

	double getRadius();
	double getVmin();
	double getVmax();
	double getMtot();
	double getRegime();

	double computeMfromV(double v);
	double computedMdVfromV(double v);

private:
	void computeRadius();
	void computeSaveVmin();
	void computeVmax();
	void computeMtot();

	double rootSolver();

	impactor* currentImpactor;
	target*   currentTarget;

	double radius;
	double vmin;
	double vmax;
	double Mtot;
	int    regime; // strength or gravity
};



#endif 