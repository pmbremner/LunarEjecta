#include "functionStrategy.h"
#include "rootSolverStrategy.h"

#include <cmath>
#include <iostream>
//#include <iomanip>
//#include <limits>

using namespace std;

rootSolverStrategy::rootSolverStrategy(functionStrategy* specificFunction){
	initFunction(specificFunction);
}

void rootSolverStrategy::initFunction(functionStrategy* specificFunction){
	cout << "rootSolverStrategy::initFunction\n";
	function = specificFunction;
}

void rootSolverStrategy::updateMaximumIterations(int n) {
	i_iter_max = n;
}


double zbrent::solve(double x_lhs, double xmin, double xmax){
	// from Section 9.3 of Numerical Recipes in C 2nd ed by Press et al 1992
	//  see also Chapter 4 of Algorithms for Minimization without Derivatives by Brent 1973
	// typedef std::numeric_limits< double > dbl;
	// std::cout.precision(dbl::max_digits10);

	int i_iter;
	double a = xmin, b = xmax, c = b, d, e, min1, min2;
	double fa = function->execute(a) - x_lhs, fb = function->execute(b) - x_lhs, fc = fb, p, q, r, s, tol1, xm;

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
		fb = function->execute(b) - x_lhs;
		
		//cout << i_iter << ' ' << b << endl;
	}
	cout << "Max iterations reached!\n";
	return b;
}


