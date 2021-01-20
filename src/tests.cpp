/*
 * tests.cpp
 *
 *  Created on: 20 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#include "tests.h"

using namespace std;

void Test_ik0_2_40() {
	const double d = 2.;
	for (int k = 20; k<= 40; k+=2){
		Bess bess(k, d);
		ifstream dfile("./data/data_ik0_2_40.txt");
		pair<double, double> max_eps_pair = {0., -1};
		if (dfile) {
		string dline;
			while (getline(dfile, dline)) {
				istringstream is(dline);
				string x_, val_;
				getline(is, x_, ',');
				getline(is, val_);
				double x = stod(x_);
				double expected = stod(val_);
				double y = bess.ik0f(x);
				double eps = abs(y-expected)/expected;
				if (eps > max_eps_pair.first) {
					max_eps_pair = {eps, x};
				}
			}
		}
		cout <<"k = " << k << ", x = " << max_eps_pair.second  << ", max eps = " << max_eps_pair.first  << "\n";
	}
	cout << numeric_limits<double>::epsilon();
}


