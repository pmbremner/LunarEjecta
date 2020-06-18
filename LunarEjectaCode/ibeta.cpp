using namespace std;

#include <iostream>
#include <cmath>

double d_r(int n, double a, double b) {
	return (n % 2 == 0 ? n/2.*(b-n/2.)/((a+n-1.)*(a+n)) : -(a+(n-1.)/2.)*(a+b+(n-1.)/2.)/((a+n-1.)*(a+n)) );
}

double beta(double a, double b) {
	// helps to avoid overflow errors doing it this way
	return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}
// 
//		cout << "n = " << n << " | d = " << d_r(n-1, a, b) <<  " | a0 = " << a0 << "| hn, kn = " << hn << ", " << kn << " | f = " << fn << endl; 
double ibeta(double x, double a, double b)
{
	if (x > (a+1.) / (a+b+2.))
		return beta(a,b) - ibeta(1.-x, b, a);

	if(x == 0.)
		return 0.;

	double f0 = 1., fn = 1., rel_err;
	double a0 = 1.;

	double h0 = 0.;
	double h1 = 1., hn;

	double k0 = 1.;
	double k1 = 1., kn;

	int n = 2;
	int nMax = 100;

	// avoids having to check every loop
	nMax = (fmod(b, 1.) < 0.0001 ? int(2*b + 1) : nMax);
	nMax = (fmod(a, 1.) < 0.0001 ? (nMax < a + 2 ? nMax : int(a + 2) ) : nMax);

	do
	{
		f0 = fn;

		a0 = 1. / (a0 * d_r(n-1, a, b) * x);
		
		hn = a0 * h1 + h0;
		kn = a0 * k1 + k0;

		fn = hn / kn;

		rel_err = fabs((f0-fn)/fn);
		n++;

		// swap
		h0 = h1;
		h1 = hn;
		k0 = k1;
		k1 = kn;

		// // check terminating series
		// if (fabs((n-1.)/2. - b) < 0.0001 || fabs(n-2. -a) < 0.0001)
		// 	n = 1000;

	} while (n < nMax && rel_err > 0.00001);

	return fn * pow(x, a) * pow(1.-x, b) / a;
}


int main(int argc, char const *argv[])
{
	double x = atof(argv[1]);
	double a = atof(argv[2]);
	double b = atof(argv[3]);

	cout << "x = " << x << endl;
	cout << "a = " << a << endl;
	cout << "b = " << b << endl;
	cout << "ibeta(" << x << "; " << a << ", " << b << ") = " << ibeta(x, a, b) << endl;


	return 0;
}