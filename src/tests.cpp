/*
 * tests.cpp
 *
 *  Created on: 20 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#include "tests.h"

using namespace std;

void Tests_ik0_0_2() {
	ifstream dfile("./data/data_ik0_0_2_0.001.txt");
	pair<double, double> max_eps_pair = {0., -1};
	if (dfile) {
	Bess bess(34, 34, 2.);
	string dline;
		while (getline(dfile, dline)) {
			istringstream is(dline);
			string x_, val_;
			getline(is, x_, ',');
			getline(is, val_);
			double x = stod(x_);
			double expected = stod(val_);
			double y = bess.ik0(x);
			double eps = abs(y-expected)/expected;
			if (eps > max_eps_pair.first) {
				max_eps_pair = {eps, x};
			}
		}
	}
	cout << "x = " << max_eps_pair.second  << ", max eps = " << max_eps_pair.first  << "\n";
}

void Test_ik0_2_40(const double d, const double xmax, const int kmin) {
	for (int k = kmin; k<= 40; k+=2){
		Bess bess(k, k, d);
		ifstream dfile("./data/data_ik0_2_40_0.001.txt");
		pair<double, double> max_eps_pair = {0., -1};
		if (dfile) {
		string dline;
			while (getline(dfile, dline)) {
				istringstream is(dline);
				string x_, val_;
				getline(is, x_, ',');
				getline(is, val_);
				double x = stod(x_);
				if (x >= d && x <= xmax) {
					double expected = stod(val_);
					double y = bess.ik0(x);
					double eps = abs(y-expected)/expected;
					if (eps > max_eps_pair.first) {
						max_eps_pair = {eps, x};
					}
				}
			}
		}
		cout <<"k = " << k << ", x = " << max_eps_pair.second  << ", max eps = " << max_eps_pair.first  << "\n";
	}
}


