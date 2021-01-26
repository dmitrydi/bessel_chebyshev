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

struct Dum {
	vector<double> coefs;
	vector<double> ns;
	double ans;
};


int main() {
	Test_accuracy();
	Test_speed();
}
