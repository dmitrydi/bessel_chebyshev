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

int main() {
	Tests_ik0_0_2();
	Test_ik0_2_40(2., 40., 30);
}
