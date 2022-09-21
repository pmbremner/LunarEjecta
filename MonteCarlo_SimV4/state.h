#ifndef STATE_H
#define STATE_H

#include "dimension.h"

using namespace std;

// virtual class, must create derived classes
class state
{
public:
	state();
	~state();

	// pure virtual function, derived class must implement this function
	virtual void nextSampleN(vector<state*> state_in, int Nsample) = 0;

	double getSample_byname(string& dim_name);
	double getSample_byidx(int dim_idx);

	//void get_dim_sample_all();


protected:
	int NDimTot;
	vector<dimension*> dim; // each dimension is responsible for knowing how to choose its sample(s)

	vector<string> sampleName; // size NDimTot
	vector<double> sample;  // Nsample x NDimTot, flattened array
};

#endif 