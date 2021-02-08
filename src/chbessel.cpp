/*
 * chbessel.cpp
 *
 *  Created on: 15 џэт. 2021 у.
 *      Author: Dmitry_Di
 */
#include "chbessel.h"

using namespace std;

namespace FastBessel {

//---PUBLIC---------

Bess::Bess(const bool fast, const int n, const int m, const double dd): is_fast(fast), d(dd), N(n), M(m),
		ak(fak(N, N, d)), ck(fck(M, M, d, 0.)),
		coef(_coef()), ns(_ns()), bess(), gs(GAUSS_POINTS) {};

double Bess::ik00x_ch(const double x) const {
	if (x < d) throw invalid_argument("in Bess::ik00x_ch(x): x < d\n");
	double sum = 0.;
	CALC_CHEB_SUM(x, sum);
	return PI2D - exp(-x)/sqrt(x)*sum;
};

long double Bess::ik00x_pwr(const double x) const {
	if (x < 0.) throw invalid_argument("in Bess::ik00x_pwr(x): x < 0");
	long double sum = 0.;
	POWER_SUM(x, sum);
	return sum;
};

double Bess::ik00x(const double x) const {
	if (x < d) {
		long double sum = 0.;
		POWER_SUM(x, sum);
		return static_cast<double>(sum);
	} else {
		double sum = 0.;
		CALC_CHEB_SUM(x, sum);
		return PI2D - exp(-x)/sqrt(x)*sum;
	}
};

double Bess::ik0ab_ch(const double x1, const double x2) const {
	if (x1 < d || x2 < d) throw invalid_argument("in Bess::ik0ab_ch(x): x < d\n");
	double sum1 = 0., sum2 = 0.;
	CALC_CHEB_SUM(x1, sum1);
	CALC_CHEB_SUM(x2, sum2);
	return exp(-x1)/sqrt(x1)*sum1 - exp(-x2)/sqrt(x2)*sum2;
}

long double Bess::ik0ab_pwr(const double x1, const double x2) const {
	if (x1 < 0. || x2 < 0.) throw invalid_argument("in Bess::ik0ab_pwr(x): x < 0");
	long double sum1 = 0., sum2 = 0.;
	POWER_SUM(x1, sum1);
	POWER_SUM(x2, sum2);
	return sum2 - sum1;
};

double Bess::ik0ab_num(const double x1, const double x2) const {
	return gs.Integrate(bess, x1, x2);
}

double Bess::ik0ab(const double x1, const double x2) const {
	if (x1 <= 0.1 || (x2-x1) > 0.5) {
		if (x2 < d) {
			if (x1 < 0. || x2 < 0.) throw invalid_argument("in Bess::ik0ab(x): x < 0");
			long double sum1 = 0., sum2 = 0.;
			POWER_SUM(x1, sum1);
			POWER_SUM(x2, sum2);
			return sum2 - sum1;
		} else if (x1 < d && x2 >= d) {
			long double sum1 = 0.;
			POWER_SUM(x1, sum1);
			double sum2 = 0.;
			if (is_fast) {
				FAST_CHEB_SUM(x2, sum2);
			} else {
				CALC_CHEB_SUM(x2, sum2);
			}
			return (PI2 - sum1) - exp(-x2)/sqrt(x2)*sum2;
		} else {
			double sum1 = 0., sum2 = 0.;
			if (is_fast) {
				FAST_CHEB_SUM(x1, sum1);
				FAST_CHEB_SUM(x2, sum2);
			} else {
				CALC_CHEB_SUM(x1, sum1);
				CALC_CHEB_SUM(x2, sum2);
			}
			return exp(-x1)/sqrt(x1)*sum1 - exp(-x2)/sqrt(x2)*sum2;
		}
	} else {
		return gs.Integrate(bess, x1, x2);
	}
}

double Bess::_k0(const double x) const {
	if (x < d) throw;
	double z = d/x;
	double z2 = 2.*(2*z-1.);
	double ch_0 = 1.;
	double ch_1 = (2.*z-1.);
	double ch_2;
	double sum = ch_0*ak[0] + ch_1*ak[1];
	for (int i = 2; i <=N; ++i) {
		ch_2 = z2*ch_1 - ch_0;
		sum += ch_2*ak[i];
		ch_0 = ch_1;
		ch_1 = ch_2;
	}
	return sum;
}

double Bess::k0(const double x) const {
	return exp(-x)/sqrt(x)*_k0(x);
}

//---PRIVATE---------

vector<long double> Bess::_coef() {
	vector<long double> ans(MAXIT_IKBESS);
	long double dum;
	long double inv_k_fact = 1.;
	for (int i = 1; i < MAXIT_IKBESS; ++i) {
		dum = 1./(long double)i;
		inv_k_fact *= dum;
		ans[i] = 1./(2.*i + 1.)*inv_k_fact*inv_k_fact;
	}
	return ans;
}

vector<long double> Bess::_ns() {
	vector<long double> ans(MAXIT_IKBESS);
	long double inv_n_sum = 0.;
	long double dum;
	for (int i = 1; i < MAXIT_IKBESS; ++i) {
		dum = 1./(long double)i;
		inv_n_sum += dum;
		ans[i] = inv_n_sum;
	}
	return ans;
}



vector<double> Bess::fak(const int m, const int n, const double d) {
	vector<__float128> fq = fakq(m,n,(__float128)d, SQRT_PI2Q); // here scaling is used SQRT_PI2Q
	vector<double> ans(fq.size());
	for (size_t i = 0; i< fq.size(); ++i) {
		ans[i] = (double)fq[i];
	}
	return ans;
}

vector<double> Bess::fck(const int m, const int n, const double d, const double mu) {
	vector<__float128> fq = fckq(m,n,(__float128)d, (__float128)mu, SQRT_PI2Q); // here scaling is used
	vector<double> ans(fq.size());
		for (size_t i = 0; i< fq.size(); ++i) {
			ans[i] = (double)fq[i];
		}
	return ans;
}

vector<__float128> Bess::fakq(const int m, const int n, const __float128 d, const __float128 multt) {
	vector<__float128> cfa(n+4, 0.q);
	vector<__float128> ans(m+1);
	cfa[n] = 1.q;
	__float128 cf1, cf2, cf3;
	__float128 ek;
	for (int k = n-1; k>=0; k--) {
		cf1 = 2.q*(k+1.q)*(1.q-(2.q*k+3.q)*(k+1.5q)*(k+1.5q)/(2.q*(k+2.q)*(k+0.5q)*(k+0.5q))-4.q*d/((k+0.5q)*(k+0.5q)));
		cf2 = (1.q-2.q*(k+1.q)*(2.q*k+3.q-4.q*d)/((k+0.5q)*(k+0.5q)));
		cf3 = -(k+1.q)*(k+2.5q)*(k+2.5q)/((k+2.q)*(k+0.5q)*(k+0.5q));
		ek = k==0 ? 1.q: 2.q;
		cfa[k] = ek/2.q*(cfa[k+1]*cf1+cfa[k+2]*cf2+cfa[k+3]*cf3);
	}
	__float128 un=0.q;
	__float128 mult = n%2==0? 1.0q : -1.0q;
	for (int k=n; k>=0; k--) {
		un += mult*cfa[k];
		mult *= -1.q;
	}
	for (int k=0; k<=m; k++) {
		ans[k] = multt*cfa[k]/un;
	}
	return ans;
}

vector<__float128> Bess::fckq(const int m, const int n, const __float128 d, const __float128 mu, const __float128 multt) {
	vector<__float128> ak = fakq(n+3, n+3, d);
	vector<__float128> dk(n+4, 0.q);
	vector<__float128> eek(n+4, 0.q);
	__float128 ek = 2.q;
	dk[n] = -ek/2.q/(n+0.5q-mu)*4.q*d*(ak[n+2]-ak[n+1]); // ?
	for (int k=n-1; k>=0; k--) {
		if (k==0) ek=1.q;
		dk[k] = -ek/2.q/(k+0.5q-mu)*(4.q*d*(ak[k+2]-ak[k+1]) +
				(3.q*k+3.5q-mu+4.q*d)*dk[k+1] +
				(3.q*k+5.5q+mu-4.q*d)*dk[k+2] +
				(k+2.5q+mu)*dk[k+3]);
		if (k==n-1) dk[n] = 0.; // ?
	}
	ek = 2.0q;
	eek[n] = 1.q;
	for (int k=n-1; k>=0; k--) {
		if (k==0) ek=1.q;
		eek[k] = -ek/2.q/(k+0.5q-mu)*( (3.q*k+3.5q-mu+4.q*d)*eek[k+1] + (3.q*k+5.5q+mu-4.q*d)*eek[k+2] +
				(k+2.5q+mu)*eek[k+3]);
	}
	__float128 sumd=0.q, sume=0.q;
	__float128 mult = n%2==0? 1.0q : -1.0q;
	for (int k=n; k>=0; k--) {
		sume += mult*eek[k];
		sumd += mult*dk[k];
		mult *= -1.q;
	}
	__float128 An = (1.q-sumd)/sume;
	vector<__float128> ans(m+1);
	for (int k=0; k<=m; ++k) {
		ans[k] = multt*(dk[k]+An*eek[k]);
	}
	return ans;
}

}


