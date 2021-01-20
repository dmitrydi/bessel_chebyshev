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
using namespace std;

double ChebyshevT_(const int n, const double x) {
	if (n==0) return 1.;
	if (n==1) return 2.*x-1.;
	return 2.*(2.*x-1.)*ChebyshevT_(n-1, x) - ChebyshevT_(n-2, x);
}

int main() {
//	const int N = 32;
//	const double d = 2.;
//	Bess bess(N, d);
//	vector<double> xs(19);
//	for (int i =0; i<19; ++i) {
//		xs[i] = i+2.;
//	}
//	vector<double> expected_ik0 = {0.09712059247806794,0.03084823672701536,0.01014675618972248,
//			0.003408936066530570,0.001161773525522248,0.0004000522249843623,
//			0.0001388415054529972,0.00004848397564095034,0.00001701517891775940,
//			5.995906150138932e-6,2.120148319899660e-6,7.518798030529142e-7,
//			2.673163078650657e-7,9.524797739457035e-8,3.400350281083561e-8,
//			1.215999458227736e-8,4.355166486884810e-9,1.561960309952770e-9,5.608816316551497e-10};
//
//	vector<double> expected_k0 = {0.1138938727495334,0.03473950438627925,0.01115967608585302,
//			0.003691098334042594,0.001243994328013123,0.0004247957418692318,
//			0.0001464707052228154,0.00005088131295645925,0.00001778006231616765,
//			6.243020547653677e-6,2.200825397311491e-6,7.784543861420496e-7,
//			2.761370823981620e-7,9.819536482396435e-8,3.499411663936499e-8,
//			1.249466402631773e-8,4.468753337309383e-9,1.600671286929361e-9,5.741237815336524e-10};
//	double max_eps_ik0 = 0.;
//	for (int i=0; i< 19; ++i) {
//		double eps = abs(bess.ik0f(xs[i])-expected_ik0[i])/expected_ik0[i];
//		if (eps > max_eps_ik0) max_eps_ik0 = eps;
//	}
//	cout << "max eps ik0 " << max_eps_ik0 << endl;
//	cout << "------\n";
//	double max_eps_k0 = 0.;
//	for (int i=0; i< 19; ++i) {
//		double eps = abs(bess.k0f(xs[i])-expected_k0[i])/expected_k0[i];
//		if (eps > max_eps_k0) max_eps_k0 = eps;
//	}
//	cout << "max eps k0 " << max_eps_k0 << endl;
//	double ans;
//	int NITER = 5'000'000;
//	vector<double> xxs = LinSpaced(2.,40.,NITER);
//	{
//		LOG_DURATION("bess");
//		for (auto x: xxs) {
//			ans = bess.ik0f(x);
//		}
//	}
//	cout << setprecision(16) << ans;
	Test_ik0_2_40();
}
