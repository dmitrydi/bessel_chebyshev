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
static const double PI2 = 1.5707'9632'6794'8966'2;
static const double EUL_GAMMA_D = 0.5772'1566'4901'5329;
static const long double EUL_GAMMA_LD = 0.5772'1566'4901'5328'6;
static const __float128 PI2Q = 1.5707963267948966192313216916397514q;
static const __float128 SQRT_PI2Q = 1.2533141373155002512078826424055226q;


double ikbess0_d(double x);

struct Bess {
	// calculates modified bessel function of second kind and its integral using Chebyshev polynomial expansion
	Bess(const int n, const int m, const double dd);
	const int MAXIT_IKBESS=20;
	const double d;
	const int N, M;
	double _k0(const double x);
	double _ik0(const double x);
	double k0(const double x);
	double ik0(const double x);
	double ik02(const double x);

	long double _ik02(const double x);
	double ikab(const double x1, const double x2);
	std::vector<long double>  _coef(), _ns();
	std::vector<double> fak(const int m, const int n, double d); // calculates coefficients of Chebyshev series for _k0(x)
	std::vector<double> fck(const int m, const int n, double d, double mu); // calculates coefficients of Chebyshev series for _ik0(x)
	std::vector<__float128> fakq(const int m, const int n, const __float128 d, const __float128 multt = 1.q); // same as fak but with quadruple presicion
	std::vector<__float128> fckq(const int m, const int n, const __float128 d, const __float128 mu, const __float128 multt = 1.q); // same as fck but with quadruple precision
	const std::vector<double> ak, ck; // coefficients for constructor Bess()
	std::vector<double> cha, chc;
	const std::vector<long double> coef, ns;
};


