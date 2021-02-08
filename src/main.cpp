//============================================================================
// Name        : bessel_cheb.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include "chbessel.h"
#include "auxillary.h"
#include "profile.h"
#include "tests.h"
#include "bessel.h"
#include "qgaus.h"

using namespace std;

const int MAXIT_IKBESS = 20;

using namespace Chebyshev;
using namespace FastBessel;

int main() {
	//Test_accuracy();
	//Test_speed();
//	const int m = 26;
//	ChPoly<__float128> mult({-2.q, 4.q});
//	vector<ChPoly<__float128>> polys = {ChPoly<__float128>({1.q}), ChPoly<__float128>({-1.q, 2.q})};
//	for (int i = 2; i <= m; ++i) {
//		polys.push_back(polys[i-1]*mult - polys[i-2]);
//	}
//	Bess bess(m, m);
//	auto cks = bess.fckq(m, m, 2.q, 0.q, SQRT_PI2Q);
//	//cout << polys.size() << ' ' << cks.size();
//	ChPoly<__float128> res;
//	for (int i = 0; i < cks.size(); ++i) {
//		res = res + polys[i]*cks[i];
//	}
//	cout << res._coefs[0] << "L+" << endl;
//	for (int i = 1; i < m; ++i) {
//		cout << "z*(" << res._coefs[i] << "L+" << endl;
//	}
//	cout << "z*" << res._coefs[m];
//	for (int i = 0; i < m-1; ++i) {
//		cout <<')';
//	}
//	for (auto x: res._coefs) {
//		cout << x << endl;
//	}
//	Bess bess;
//	double x = 5.;
//	double ans1, ans2, ans3;
//	{
//		LOG_DURATION("precise");
//		for (int i = 0; i < 10'000'000; ++i) {
//			ans1 = bess._ik0(x);
//		}
//	}
//	{
//		LOG_DURATION("fast 1");
//		for (int i = 0; i < 10'000'000; ++i) {
//			ans2 = Stupid_ik0(x);
//		}
//	}
//	{
//		LOG_DURATION("fast 2");
//		for (int i = 0; i < 10'000'000; ++i) {
//			ans3 = Stupid_ik02(x);
//		}
//	}
//	cout << ans1 - ans2 << endl;
//	cout << ans1 - ans3;
	Bess bess;
	double x1 = 3.;
	double ans1, ans2, ans3;
	int n = 20'000'000;
	{ LOG_DURATION("function");
		for (int i = 0; i < n; ++i)
			ans1 = bess.analytic_ik0ab(x1, x1+0.5);
	}
	{ LOG_DURATION("macro");
		for (int i = 0; i < n; ++i)
			ans2 = bess.ik0ab_macro(x1, x1+0.5);
	}
	{ LOG_DURATION("stupid");
		for (int i = 0; i < n; ++i)
			ans3 = ik0ab_stupid(x1, x1+0.5);
	}
	cout << ans1 << endl;
	cout << ans2 << endl;
	cout << ans3 << endl;
	cout << ans1 - ans2 << endl;
	cout << ans1 - ans3 << endl;
}
