/*
 * tests.cpp
 *
 *  Created on: 20 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#include "tests.h"
#include "auxillary.h"
#include "profile.h"

using namespace std;

using namespace FastBessel;

struct Ans {
	double x0;
	double dx;
	double eps = 0.;
};

ostream& operator<<(ostream& os, const DataPoint p) {
	os << "x = " << p.x << " dx = " << p.dx << " y = " << p.val;
	return os;
}

ostream& operator<<(ostream& os, const ResultPoint& point) {
	os << "{";
	os << point.method_name +" "+ to_string(point.napprox);
	os << " niter: " << point.niter;
	os << " x: " << point.x;
	os << " dx: " << point.dx;
	os << " maxeps: " << point.eps;
	os << " time: " << point.time << " microseconds";
	os << "}";
	return os;
}

vector<DataPoint> ParseDataFile(const std::string& fname,
		const double xmin,
		const double xmax,
		const double dxmin,
		const double dxmax) {
	vector<DataPoint> ans;
	ifstream dfile(fname);
	if (dfile) {
		string dline;
		while(getline(dfile, dline)) {
			istringstream is(dline);
			string x_, dx_, val_;
			getline(is, x_, ',');
			getline(is, dx_, ',');
			getline(is, val_);
			double x = stod(x_);
			double dx = stod(dx_);
			double val = stod(val_);
			if ((xmin<0.||x>=xmin) && (xmax<0.||x<=xmax) && \
					(dxmin<0.||dx>=dxmin) && (dxmax<0.||dx<=dxmax)) {
				ans.push_back({x, dx, val});
			}
		}
	}
	return ans;
}

void Test_accuracy() {
	string fname = "./data/data_ikab_0_40_dx_0.01_0.0001_100.txt";
	auto data = ParseDataFile(fname, -1, -1, -1, -1);
	{
		double maxeps = 0.;
		Bess bess;
		{ LOG_DURATION("fast = false");
			for (auto p: data) {
				double y = bess.ik0ab(p.x, p.x+p.dx);
				double eps = abs(p.val-y)/p.val;
				if (eps > maxeps) maxeps = eps;
			}
		}
		cout << "max eps: " << maxeps << endl;
	}
	{
		double maxeps = 0.;
		Bess bess(true);
		{ LOG_DURATION("fast = true");
			for (auto p: data) {
				double y = bess.ik0ab(p.x, p.x+p.dx);
				double eps = abs(p.val-y)/p.val;
				if (eps > maxeps) maxeps = eps;
			}
		}
		cout << "max eps: " << maxeps << endl;
	}
}


void Test_speed() {
	{
		Bess bess;
		vector<double> xs = LinSpaced(0., 40., 1000);
		vector<double> dxs = LogSpaced(0.0001, 10., 1000);
		double ans;
		std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
		for (auto x: xs) {
			for (auto dx: dxs) {
				ans = bess.ik0ab(x, x+dx);
			}
		}
		std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
		auto dur = finish - start;
		double t = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
		size_t niter = xs.size()*dxs.size();
		cout << "Precise solution:: points: " << niter << " Total time: " << t << " us " << " Average per point: " << t/niter << " us\n";
		cout << ans << endl;
	}
	{
		Bess bess(true);
		vector<double> xs = LinSpaced(0., 40., 1000);
		vector<double> dxs = LogSpaced(0.0001, 10., 1000);
		double ans;
		std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
		for (auto x: xs) {
			for (auto dx: dxs) {
				ans = bess.ik0ab(x, x+dx);
			}
		}
		std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
		auto dur = finish - start;
		double t = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
		size_t niter = xs.size()*dxs.size();
		cout << "Fast solution:: points: " << niter << " Total time: " << t << " us " << " Average per point: " << t/niter << " us\n";
		cout << ans << endl;
	}
}
