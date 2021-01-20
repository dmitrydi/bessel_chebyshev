/*
 * chbessel.cpp
 *
 *  Created on: 15 џэт. 2021 у.
 *      Author: Dmitry_Di
 */
#include "chbessel.h"

using namespace std;

Bess::Bess(): d(2.),  ak({1.2201515410329777e+00,
	   -3.1448101311964500e-02,
	    1.5698838857300533e-03,
	   -1.2849549581627802e-04,
	    1.3949813718876499e-05,
	   -1.8317555227191195e-06,
	    2.7668136394450150e-07,
	   -4.6604898976879476e-08,
	    8.5740340174142259e-09,
	   -1.6975345093890615e-09,
	    3.5773972814003284e-10,
	   -7.9574892444773969e-11,
	    1.8559491149549265e-11,
	   -4.5145978833745191e-12,
	    1.1403405882073442e-12,
	   -2.9800969231481783e-13,
	    8.0328907750683742e-14,
	   -2.2275133267462963e-14,
	    6.3400764762766458e-15,
	   -1.8485933779209071e-15,
	    5.5120559994043333e-16,
	   -1.6782311257549006e-16
}), ck({1.1206780274986525e+00,
	   -1.1715459069850548e-01,
	    1.3036970236823027e-02,
	   -1.9806761821353283e-03,
	    3.6376878778932641e-04,
	   -7.6299030718489181e-05,
	    1.7691407526436762e-05,
	   -4.4412599633740433e-06,
	    1.1899212195550536e-06,
	   -3.3672827874851224e-07,
	    9.9858504718735408e-08,
	   -3.0845382476077599e-08,
	    9.8762264407071931e-09,
	   -3.2649783182163052e-09,
	    1.1108335763913329e-09,
	   -3.8789825589294646e-10,
	    1.3870353745126009e-10,
	   -5.0687673062241771e-11,
	    1.8898287875308313e-11,
	   -7.1779913917234740e-12,
	    2.7738142687119739e-12,
	   -1.0892937934539412e-12,
	    4.3427063743150941e-13,
	   -1.7560032645303451e-13,
	    7.1958402502578853e-14,
	   -2.9861094323682945e-14,
	    1.2540208296211895e-14,
	   -5.3261294841264028e-15,
	    2.2865611125307969e-15,
	   -9.9173250041067941e-16,
	    4.3435161400121426e-16,
	   -1.9201449922763557e-16

}), sNa(ak.size()), sNc(ck.size()), Na(sNa-1), Nc(sNc-1), yya(sNa+2, 0.), yyc(sNc+2, 0.) {};

Bess::Bess(const int n, const double dd):  d(dd), ak(fak(n, n, d)),
		ck(fck(n, n, d, 0.)),  sNa(ak.size()), sNc(ck.size()), Na(sNa-1),
		Nc(sNc-1), yya(sNa+2, 0.), yyc(sNc+2, 0.), ch(sNc+1, 1.){};

double Bess::k0f(const double x) {
	double z = d/x;
	ch[0] = 1.;
	ch[1] = 2.*z-1.;
	for (int i=2; i<=sNa; ++i) {
		ch[i] = 2.*(2.*z-1)*ch[i-1] - ch[i-2];
	}
	double sum=0.;
	for (int i=sNa; i>=0; --i) {
		sum += ch[i]*ak[i];
	}
	return exp(-x)/sqrt(x)*sum;
}

double Bess::ik0f(const double x) {
	double z = d/x;
	ch[0] = 1.;
	ch[1] = 2.*z-1.;
		for (int i=2; i<=sNc; ++i) {
			ch[i] = 2.*(2.*z-1)*ch[i-1] - ch[i-2];
	}
	double sum = 0.;
	for (int i=sNc; i>=0; --i) {
			sum += ch[i]*ck[i];
	}
	return exp(-x)/sqrt(x)*sum;
}

double Bess::_k0(const double x) {
	if (x < d) throw;
	double z = x/d;
	double z2 = 2.*(2.*z-1.);
	for(int k=Na; k<=1; --k) {
		yya[k] = z2*yya[k+1] - yya[k+2] + ak[k];
	}
	return z2*yya[1] - yya[2] + ak[0];
}

double Bess::_ik0(const double x) {
	if (x<d) throw;
	double z = x/d;
	double z2 = 2.*(2.*z-1.);
	for(int k=Nc; k<=1; --k) {
		yyc[k] = z2*yyc[k+1] - yyc[k+2] + ck[k];
	}
	return z2*yyc[1] - yyc[2] + ck[0];
}

double Bess::k0(const double x) {
	return exp(-x)/sqrt(x)*_k0(x);
}

double Bess::ik0(const double x) {
	return exp(-x)/sqrt(x)*_ik0(x);
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
	vector<__float128> ak = fakq(n+4, 5*n, d);
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




