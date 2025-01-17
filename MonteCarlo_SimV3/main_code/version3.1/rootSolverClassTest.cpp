#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

template<class T>
class inverseFunction
{
public:
	inverseFunction(T);
	~inverseFunction();

	virtual double fx(double x) = 0;
	double fx_inv(double y);

protected:
	
	double zbrent(double y);

	T f_params;
	double xmin;
	double xmax;
	int i_iter_max;
	int eps;
	int tol;

};


template<class T>
inverseFunction<T>::inverseFunction(T Tparams){
	f_params = Tparams;
	i_iter_max = 100;
	eps = 3.E-8;
	tol = 1.E-3;
}

template<class T>
double inverseFunction<T>::zbrent(double y){
	int i_iter;
	double a = xmin, b = xmax, c = b, d, e, min1, min2;
	double fa = fx(a)-y, fb = fx(b)-y, fc = fb, p, q, r, s, tol1, xm;

	for (i_iter = 0; i_iter < i_iter_max; ++i_iter){
		// rename a, b, c and adjust bounding interval
		if (fb*fc > 0.){
			c = a;
			fc = fa;
			e = d = b-a;
		}
		if (fabs(fc) < fabs(fb)){
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol1 = 2.*eps*fabs(b) + 0.5*tol; // check convergence
		xm = 0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0){
			return b;
		}
		// attempt inverse quadratic interplolation
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)){
			s = fb/fa;
			if (a == c){
				p = 2.*xm*s;
				q = 1. - s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*(2.*xm*q*(q - r) - (b - a)*(r - 1.));
				q = (q - 1.)*(r - 1.)*(s - 1.);
			}
			// check whether in bounds
			if (p > 0.){
				q *= -1.;
			}
			p = fabs(p);

			min1 = 3.*xm*q - fabs(tol1*q);
			min2 = fabs(e*q);

			if (2.*p < (min1 < min2 ? min1 : min2)){
				// Accept interpolation
				e = d;
				d = p/q;
			} else {
				// interpolation failed, use bisection
				d = xm;
				e = d;
			}
		} else {
			// Bounds decreasing too slowly, use bisection
			d = xm;
			e = d;
		}
		// move last best guess to a
		a = b;
		fa = fb;

		// evaluate new trial root
		if (fabs(d) > tol1){
			b += d;
		} else {
			b += (xm > 0. ? fabs(tol1) : -fabs(tol1));
		}
		fb = fx(b)-y;
		cout << setprecision(20) << i_iter << ' ' << b << endl;
	}
	cout << "Max iterations reached!\n";
	return b;
}


template<class T>
double inverseFunction<T>::fx_inv(double y) {
	return zbrent(y);
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


template<class T>
class myFunc: public inverseFunction<T>
{
public:
	myFunc(T);
	~myFunc();

	double fx(double x);
	
};

template<class T>
myFunc<T>::myFunc(T Tparams){
	//inverseFunction<T>(Tparams);

	xmin = 0.;
	xmax = 6.28;
}

template<class T>
double myFunc<T>::fx(double x) {
	return sin(x);
}



int main(int argc, char const *argv[])
{
	int a = 8;
	myFunc<int*> myActualFunc = new myFunc<int*>(&a);

	cout << myActualFunc.fx_inv(0.5) << endl;

	return 0;
}