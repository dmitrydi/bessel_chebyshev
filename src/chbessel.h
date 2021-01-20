/*
 * chbessel.h
 *
 *  Created on: 15 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#pragma once

#include <vector>
#include <quadmath.h>
#include <cmath>

static const double PI = 3.141592653589793;
static const double PI2 = 1.570796326794897;
static const __float128 PI2Q = 1.5707963267948966192313216916397514q;
static const __float128 SQRT_PI2Q = 1.2533141373155002512078826424055226q;

template <typename T>
T ChebyshevT_(const int n, const T x) {
	if (n==0) return 1.;
	if (n==1) return 2.*x-1.;
	return 2.*(2.*x-1.)*ChebyshevT_(n-1, x) - ChebyshevT_(n-2, x);
}

struct Bess {
	// calculates modified bessel function of second kind and its integral using Chebyshev polynumial expansion
	Bess(); // default constructor using 21 (up to x^20) members in  series expansion and d=2 as a parameter in Chebyshev approximation
	Bess(const int n, const double dd); // constructor using custom number n of members and parameter d in Chebyshev series expansion
	const double d;
	double k0f(const double x);
	double ik0f(const double x);
	double k0(const double x); // returns BesselK0(x)
	double ik0(const double x); // returns Integral[BesselK0(t), {t, x, inf}]
	double _k0(const double x); // returns scaled function: x^1/2 * Exp[x] * BesselK0(x)
	double _ik0(const double x); // returns scaled function: x^1/2 * Exp[x] * Integral[BesselK0(t), {t, x, inf}]
	std::vector<double> fak(const int m, const int n, double d); // calculates coefficients of Chebyshev series for _k0(x)
	std::vector<double> fck(const int m, const int n, double d, double mu); // calculates coefficients of Chebyshev series for _ik0(x)
	std::vector<__float128> fakq(const int m, const int n, const __float128 d, const __float128 multt = 1.q); // same as fak but with quadruple presicion
	std::vector<__float128> fckq(const int m, const int n, const __float128 d, const __float128 mu, const __float128 multt = 1.q); // same as fck but with quadruple precision
	const std::vector<double> ak, ck; // coefficients for constructor Bess()
	const int sNa, sNc, Na, Nc;
	std::vector<double> yya, yyc;
	std::vector<double> ch;
};


