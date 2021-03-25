#include "LunarEjecta_params.h"
#include "LunarEjecta_speedzenith.h"

using namespace std;

// g++ -O2 test_speedzenith.cpp LunarEjecta_speedzenith.cpp LunarEjecta_params.cpp -o test_speedzenith.exe

int main(int argc, char const *argv[])
{

	double D0 = atof(argv[1]);
	double D1 = atof(argv[2]);
	int max_level = atoi(argv[3]);
	
	double qsize = 0., ratio = 10., qsum_romb;

	vector<quad> q;
	vector<double> qsum(max_level, 0.);

	for (int i = 0; i < max_level; ++i)
	{
		getEvalQuads(q, D0, D1, 0.0, PI/4., 0.1/2.4, 1., i+1, 0);

		ratio = 200. * qsize / q.size();
		qsize = q.size();

		qsum[i] = getSumQuad(q);

		if (i > 0)
			qsum_romb = 2.*qsum[i] - qsum[i-1];


		cout << i+1 << ' ' << ratio << ' ' << qsize << ' ' << qsum[i] << ' ' << qsum_romb << endl;
		q.resize(0);
	}

	//getEvalQuads(q, D0, D1, 0.0, PI/4., 0.1/2.4, 1., max_level, 0);
	//getEvalQuads(q, D0, D1, 0.0, PI/2., 0.1/2.4, 0.4, max_level, 0);

	//cout << q.size() << endl;

	//printQuad(q);

	return 0;
}