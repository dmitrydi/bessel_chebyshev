/*
 * chbessel.h
 *
 *  Created on: 15 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#pragma once

#include <vector>
#include <quadmath.h>
#include <limits>
#include "bessel.h"
#include "qgaus.h"
#include <cmath>
#include <cstdint>
#include <iostream>

#define SQR(x) \
	x*x;

#define CYCLE(x) \
	for (int i = 0; i < N; ++i) { \
		sum += x; \
	};


#define CACL_CHEB_SUM(z, sum) \
	{						  \
		double z2 = 2.*(2*z-1.); \
		double ch_0 = 1.; \
		double ch_1 = (2.*z-1.); \
		double ch_2; \
		sum = ch_0*ck[0] + ch_1*ck[1]; \
		for (int i = 2; i <= M; ++i) { \
			ch_2 = z2*ch_1 - ch_0; \
			sum += ch_2*ck[i]; \
			ch_0 = ch_1; \
			ch_1 = ch_2; \
		} \
	};

#define STUPID_BESS(z, sum) \
		sum = 1.25331413731550e+00+ \
			z*(-3.91660667890412e-01+ \
			z*(3.15776408289850e-01+ \
			z*(-4.06194424195295e-01+ \
			z*(7.19597732003498e-01+ \
			z*(-1.62708622337386e+00+ \
			z*(4.46581915059184e+00+ \
			z*(-1.42188765206382e+01+ \
			z*(4.97110608920295e+01+ \
			z*(-1.78790459190384e+02+ \
			z*(6.21256867237806e+02+ \
			z*(-1.98596426994130e+03+ \
			z*(5.65054653532254e+03+ \
			z*(-1.40148864771337e+04+ \
			z*(2.99068912288121e+04+ \
			z*(-5.44194703279499e+04+ \
			z*(8.38397009322609e+04+ \
			z*(-1.08611516820013e+05+ \
			z*(1.17395578367076e+05+ \
			z*(-1.04843244388150e+05+ \
			z*(7.63708826316075e+04+ \
			z*(-4.45748586206830e+04+ \
			z*(2.03233348896703e+04+ \
			z*(-6.96738114384029e+03+ \
			z*(1.68804707483831e+03+ \
			z*(-2.57529748166367e+02+ \
			z*1.85968591391546e+01)))))))))))))))))))))))));

/*
 *double Bess::_ik0(const double x) const {
	if (x < d) throw;
	double z = d/x;
	double z2 = 2.*(2*z-1.);
	double ch_0 = 1.;
	double ch_1 = (2.*z-1.);
	double ch_2;
	double sum = ch_0*ck[0] + ch_1*ck[1];
	for (int i = 2; i <=M; ++i) {
		ch_2 = z2*ch_1 - ch_0;
		sum += ch_2*ck[i];
		ch_0 = ch_1;
		ch_1 = ch_2;
	}
	return sum;
}
 */

namespace Chebyshev {

template <typename T>
struct ChPoly {
	ChPoly(): _coefs({}), _pwr(0) {};
	ChPoly(const std::vector<T> coefs): _coefs(coefs), _pwr(_coefs.size()) {};
	ChPoly(const ChPoly& other): _coefs(other._coefs), _pwr(other._pwr) {};
	ChPoly multx(const int n) const {
		std::vector<T> ans(_pwr+n, 0);
		for (int i = 0; i < _pwr; ++i) {
			ans[i+n] = this->_coefs[i];
		}
		return ChPoly(ans);
	};
	template <typename U>
	ChPoly operator*(const U val) const {
		ChPoly ans(*this);
		for (auto& cf: ans._coefs) {
			cf *= val;
		}
		return ans;
	};
	ChPoly operator-(const T val) const;
	ChPoly operator+(const T val) const;
	ChPoly operator-(const ChPoly& other) const {
		if (_pwr > other._pwr) {
			ChPoly ans(*this);
			for (int i = 0; i < other._pwr; ++i) {
				ans._coefs[i] -= other._coefs[i];
			}
			return ans;
		} else {
			ChPoly ans(other);
			for (int i = 0; i < _pwr; ++i) {
				ans._coefs[i] -= _coefs[i];
			}
			return ans;
		}
	};
	ChPoly operator+(const ChPoly& other) const {
		if (_pwr > other._pwr) {
			ChPoly ans(*this);
			for (int i = 0; i < other._pwr; ++i) {
				ans._coefs[i] += other._coefs[i];
			}
			return ans;
		} else {
			ChPoly ans(other);
			for (int i = 0; i < _pwr; ++i) {
				ans._coefs[i] += _coefs[i];
			}
			return ans;
		}
	};
	ChPoly operator*(const ChPoly& other) const {
		ChPoly ans;
		for (int i = 0; i < other._pwr; ++i) {
			ChPoly dum = (this->multx(i))*other._coefs[i];
			ans = ans + dum;
		}
		return ans;
	}
	std::vector<T> _coefs;
	int _pwr;
};

struct Poly {
	Poly();
	Poly(const std::vector<int64_t> coefs);
	Poly(const Poly& other);
	std::vector<int64_t> GetCoefs() const;
	Poly multx(const int n) const;
	Poly operator*(const int64_t val) const;
	Poly operator-(const int64_t val) const;
	Poly operator+(const int64_t val) const;
	Poly operator*(const Poly& other) const;
	Poly operator-(const Poly& other) const;
	Poly operator+(const Poly& other) const;
	int power() const;
	std::vector<int64_t> _coefs;
	int _pwr;
};


}

std::ostream& operator<<(std::ostream& os, __float128 x);

template <typename T>
std::ostream& operator << (std::ostream& os, const Chebyshev::ChPoly<T>& p) {
	int n = p._pwr - 1;
	for (auto it = p._coefs.rbegin(); it != p._coefs.rend(); ++it, n--) {
		if (*it != 0) {
			if (it == p._coefs.rbegin()) {
				os << *it << "x^" << n;
			} else {
				os << ' ';
				if (*it > 0) {
					os << '+';
				} else {
					os << '-';
				}
				os << ' ';
				os << fabsq(*it);
				if (n>0) os << "x^" << n;
			}
		}
	}
	return os;
}

std::ostream& operator << (std::ostream& os, const Chebyshev::Poly& p);

namespace FastBessel {

static const double PI = 3.141592653589793;
static const long double PI2 = 1.5707'9632'6794'8966'2;
static const double EUL_GAMMA_D = 0.5772'1566'4901'5329;
static const long double EUL_GAMMA_LD = 0.5772'1566'4901'5328'6;
static const __float128 PI2Q = 1.5707963267948966192313216916397514q;
static const __float128 SQRT_PI2Q = 1.2533141373155002512078826424055226q;
static const double TINY = std::numeric_limits<double>::min();
static const double s_cf[27] = {
		1.25331413731550e+00,
		-3.91660667890412e-01,
		3.15776408289850e-01,
		-4.06194424195295e-01,
		7.19597732003498e-01  ,
		-1.62708622337386e+00 ,
		4.46581915059184e+00  ,
		-1.42188765206382e+01 ,
		4.97110608920295e+01  ,
		-1.78790459190384e+02 ,
		6.21256867237806e+02  ,
		-1.98596426994130e+03 ,
		5.65054653532254e+03  ,
		-1.40148864771337e+04 ,
		2.99068912288121e+04  ,
		-5.44194703279499e+04 ,
		8.38397009322609e+04  ,
		-1.08611516820013e+05 ,
		1.17395578367076e+05  ,
		-1.04843244388150e+05 ,
		7.63708826316075e+04  ,
		-4.45748586206830e+04 ,
		2.03233348896703e+04  ,
		-6.96738114384029e+03 ,
		1.68804707483831e+03  ,
		-2.57529748166367e+02 ,
		1.85968591391546e+01
};

double Stupid_ik0(const double x);
double Stupid_ik02(const double x);
double ik0ab_stupid(const double x1, const double x2);


struct Bess {
	// calculates modified bessel function of second kind and its integral using Chebyshev polynomial expansion
	Bess(const int n = 34, const int m = 34, const double dd = 2.);
	const int MAXIT_IKBESS=20;
	const int GAUSS_POINTS = 20;
	const double d;
	const int N, M;
	/*
	 * double _ik0_inf(const double x) // x>d
	 * double ik0_inf(const double x)  // x>d
	 * double ik0_x2(const double x)   // x<2
	 *
	 */
	double _k0(const double x) const; // sqrt(x)*exp(x)*BesselK[0,x], x >= d;
	double _ik0(const double x) const; // sqrt(x)*exp(x)*Integral[BesselK[x, inf]], x >= d;
	double k0(const double x) const; // BesselK[0,x], x >= d
	double ik0(const double x) const; // Integral BessK0[x, inf], x >= d
	double ik02(const double x) const; // Integral BessK0[x, inf], x < 2
	long double _ik02(const double x) const; // Integral BessK0[0, x], x < 2
	double analytic_ik0ab(const double x1, const double x2) const; // integral BessK0[x1, x2]
	double ik0ab(const double x1, const double x2) const;
	std::vector<long double>  _coef(), _ns();
	std::vector<double> fak(const int m, const int n, double d); // calculates coefficients of Chebyshev series for _k0(x)
	std::vector<double> fck(const int m, const int n, double d, double mu); // calculates coefficients of Chebyshev series for _ik0(x)
	std::vector<__float128> fakq(const int m, const int n, const __float128 d, const __float128 multt = 1.q); // same as fak but with quadruple presicion
	std::vector<__float128> fckq(const int m, const int n, const __float128 d, const __float128 mu, const __float128 multt = 1.q); // same as fck but with quadruple precision
	const std::vector<double> ak, ck; // coefficients for constructor Bess()
	std::vector<double> cha, chc;
	const std::vector<long double> coef, ns;
	BessK0 bess;
	GaussIntegrator gs;
	double num_ik0ab(const double x1, const double x2) const;
	double dum(const double x);
	double dum2(const double x);
	double ik0ab_macro(const double x1, const double x2) const;

};

}
