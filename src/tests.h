/*
 * tests.h
 *
 *  Created on: 20 џэт. 2021 у.
 *      Author: Dmitry_Di
 */

#pragma once
#include "chbessel.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <utility>
#include <limits>
#include "qgaus.h"
#include <chrono>

struct DataPoint {
	double x;
	double dx;
	double val;
};

struct ResultPoint {
	std::string method_name;
	int niter;
	int napprox;
	double x;
	double dx;
	double eps = 0.;
	double time;
};

std::ostream& operator<<(std::ostream& os, const DataPoint p);

std::ostream& operator<<(std::ostream& os, const ResultPoint& point);

std::vector<DataPoint> ParseDataFile(const std::string& fname,
		const double xmin=-1.,
		const double xmax=-1.,
		const double dxmin=-1.,
		const double dxmax=-1.);

void Test_accuracy();
void Test_speed();
void Test_speed_bessel();

