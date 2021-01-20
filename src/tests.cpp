/*
 * tests.cpp
 *
 *  Created on: 20 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#include "tests.h"

using namespace std;

void Test_ik0_2_40() {
	const int N = 32;
	const double d = 2.;
	const double ROUND = 100.;

	Bess bess(N, d);

	ifstream dfile("./data/data_ik0_2_40.txt");

	pair<double, double> max_eps_pair = {0., -1};

	if (dfile) {
	string dline;
		while (getline(dfile, dline)) {
			istringstream is(dline);
			string x_, val_;
			getline(is, x_, ',');
			getline(is, val_);
			double x = round(stod(x_)*ROUND)/ROUND;
			double expected = stod(val_);
			double y = bess.ik0f(x);
			double eps = abs(y-expected)/expected;
			if (eps > max_eps_pair.first) {
				max_eps_pair = {eps, x};
			}
		}
	}
	cout << "max eps = " << max_eps_pair.first << ", x = " << max_eps_pair.second;
}


