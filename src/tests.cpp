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
	double maxeps = 0.;
	Bess bess;
	{ LOG_DURATION(to_string(data.size()));
	for (auto p: data) {
		double y = bess.ik0ab(p.x, p.x+p.dx);
		double eps = abs(p.val-y)/p.val;
		if (eps > maxeps) maxeps = eps;
	}
	}
	cout << "max eps: " << maxeps << endl;
}

std::vector<ResultPoint> Test_cheb_gauss(const std::string& fname,
		const int ncheb_min,
		const int ncheb_max,
		const int ngauss_min,
		const int ngauss_max,
		const double xmin,
		double xmax,
		double dxmin,
		double dxmax) {
	auto data = ParseDataFile(fname, xmin, xmax, dxmin, dxmax);
	const int ndata = data.size();
	vector<ResultPoint> res;
	if (ndata>0) {
		// chebyshev
		for (int k=ncheb_min; k <= ncheb_max; k+=2) {
			Bess bess(k,k,2.);
			vector<double> vexpected(ndata);
			// compute time
			std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
			for (int i = 0; i < ndata; ++i) {
				vexpected[i] = bess.analytic_ik0ab(data[i].x, data[i].x+data[i].dx);
			}
			std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
			auto dur = finish - start;
			double t = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
			//
			double maxeps = 0., xeps = 0., dxeps = 0.;
			for (int i = 0; i < ndata; ++i) {
				double eps = abs(vexpected[i] - data[i].val)/data[i].val;
				if (eps > maxeps) {
					maxeps = eps;
					xeps = data[i].x;
					dxeps = data[i].dx;
				}
			}
			res.push_back({"cheb", ndata, k, xeps, dxeps, maxeps, t});
		}
		// gauss
		for (int k=ngauss_min; k <= ngauss_max; k+=10) {
			Bessik bess;
			auto func = [&bess](double x) {return bess.k0(x);};
			GaussIntegrator gs(k);
			vector<double> vexpected(ndata);
			// compute time
			std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
			for (int i = 0; i < ndata; ++i) {
				vexpected[i] = gs.Integrate(func, data[i].x, data[i].x+data[i].dx);
			}
			std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
			auto dur = finish - start;
			double t = std::chrono::duration_cast<std::chrono::microseconds>(dur).count();
			//
			double maxeps = 0., xeps = 0., dxeps = 0.;
			for (int i = 0; i < ndata; ++i) {
				double eps = abs(vexpected[i] - data[i].val)/data[i].val;
				if (eps > maxeps) {
					maxeps = eps;
					xeps = data[i].x;
					dxeps = data[i].dx;
				}
			}
			res.push_back({"gauss", ndata, k, xeps, dxeps, maxeps, t});
		}
	}
	return res;
}


void Test_speed() {
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
	cout << "Number of points: " << niter << " Total time: " << t << " us " << " Average per point: " << t/niter << " us\n";
	cout << ans << endl;
}

//void Test_gauss (const double xmin, const double xmax, const int ngaus) {
//	//ifstream dfile("./data/data_ikab_0.1_20_0.0001_1.txt");
//	ifstream dfile("./data/data_ikab_0.001_0.1_0.0001_1.txt");
//	if (dfile) {
//		Ans ans;
//		Bessik bess;
//		auto func = [&bess](double x) {return bess.k0(x);};
//		GaussIntegrator gs(ngaus);
//		string dline;
//		while(getline(dfile, dline)) {
//			istringstream is(dline);
//			string x_, dx_, val_;
//			getline(is, x_, ',');
//			getline(is, dx_, ',');
//			getline(is, val_);
//			double x = stod(x_);
//			double dx = stod(dx_);
//			double expected = stod(val_);
//			if (dx < 1.) {
//				double y = gs.Integrate(func, x, x+dx);
//				double eps = abs(y-expected)/expected;
//				if (eps > ans.eps && x >= xmin && x <= xmax) {
//					ans.eps = eps;
//					ans.x0 = x;
//					ans.dx = dx;
//				}
//			}
//		}
//		cout << "x = " << ans.x0 << ", dx = " << ans.dx << ", eps = " << ans.eps << endl;
//	}
//}
//
//void Test_ikab() {
//	//ifstream dfile("./data/data_ikab_0_20_0.001_100.txt");
//	ifstream dfile("./data/data_ikab_0.001_0.1_0.0001_1.txt");
//	if (dfile) {
//		Ans ans;
//		Bess bess(36,36,2.);
//		string dline;
//		while(getline(dfile, dline)) {
//			istringstream is(dline);
//			string x_, dx_, val_;
//			getline(is, x_, ',');
//			getline(is, dx_, ',');
//			getline(is, val_);
//			double x = stod(x_);
//			double dx = stod(dx_);
//			double expected = stod(val_);
//			double y = bess.ik0ab(x, x+dx);
//			double eps = abs(y-expected)/expected;
//			if (eps > ans.eps) {
//				ans.eps = eps;
//				ans.x0 = x;
//				ans.dx = dx;
//			}
//		}
//		cout << "x = " << ans.x0 << ", dx = " << ans.dx << ", eps = " << ans.eps << endl;
//	}
//
//}
//
//void Tests_ik0_0_2() {
//	ifstream dfile("./data/data_ik0_0_2_0.001.txt");
//	pair<double, double> max_eps_pair = {0., -1};
//	if (dfile) {
//	Bess bess(34, 34, 2.);
//	string dline;
//		while (getline(dfile, dline)) {
//			istringstream is(dline);
//			string x_, val_;
//			getline(is, x_, ',');
//			getline(is, val_);
//			double x = stod(x_);
//			double expected = stod(val_);
//			double y = bess.ik0(x);
//			double eps = abs(y-expected)/expected;
//			if (eps > max_eps_pair.first) {
//				max_eps_pair = {eps, x};
//			}
//		}
//	}
//	cout << "x = " << max_eps_pair.second  << ", max eps = " << max_eps_pair.first  << "\n";
//}
//
//void Test_ik0_2_40(const double d, const double xmax, const int kmin) {
//	for (int k = kmin; k<= 40; k+=2){
//		Bess bess(k, k, d);
//		ifstream dfile("./data/data_ik0_2_40_0.001.txt");
//		pair<double, double> max_eps_pair = {0., -1};
//		if (dfile) {
//		string dline;
//			while (getline(dfile, dline)) {
//				istringstream is(dline);
//				string x_, val_;
//				getline(is, x_, ',');
//				getline(is, val_);
//				double x = stod(x_);
//				if (x >= d && x <= xmax) {
//					double expected = stod(val_);
//					double y = bess.ik0(x);
//					double eps = abs(y-expected)/expected;
//					if (eps > max_eps_pair.first) {
//						max_eps_pair = {eps, x};
//					}
//				}
//			}
//		}
//		cout <<"k = " << k << ", x = " << max_eps_pair.second  << ", max eps = " << max_eps_pair.first  << "\n";
//	}
//}
//
//
