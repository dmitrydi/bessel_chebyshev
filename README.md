# bessel_cheb

Provides FastBessel::Bess class for fast and accurate calculation of Modified Bessel function of second kind (K0) integral using Chebyshev polynomial expansion. 
Uses quadruple-precision calculations internally to provide at least 1e-14 relative precision of  K0(t)dt \[0,x\] integral calculation for 0 <= x < Inf.   

## Installation

Put __chbelssel__ .h and .cpp files along with __qgaus__ .h and .cpp in your working directory.

### Requirements:

C++11 or later version compiler

### Required compile options:
```
-fext-numeric-literals
```
### Required linker options:
```
-lquadmath
```

## Usage

Create a `FastBessel::Bess` default-constructed instance for calculations. Although constructor takes some parameters, don't change them unless you know how it affects precision.

```C++
#include "chbessel.h"

FastBessel::Bess bess;

// make necessary calculations

```

## Member functions
function | description
-------|------------
[__(constructor)__](#constructor) | creates an instance
[__ik00x__](#ik00x) | calculates K0 integral in \[0, x\] for any x
[__ik00x\_ch__](#ik00x_ch) | calculates K0 integral in \[0, x\] , using Chebyshev expansion only
[__ik00x\_pwr__](#ik00x_pwr) | calculates K0 integral in \[0, x\] using power series expansion only
[__ik0ab__](#ik0ab)| calculates K0 integral in \[a, b\], 0 <= a < b < Inf
[__ik0ab\_ch__](#ik0ab_ch) | calculates K0 integral in \[a, b\], 2. <= a < b < Inf using Chebyshev expansion only (see notes)
[__ik0ab\_pwr__](#ik0ab_pwr) | calculates K0 integral in \[a, b\], 0 <= a < b <= 20. using power series expansion only (see notes)
[__ik0ab\_num__](#ik0ab_num) | calculates K0 integral in \[a, b\] using Gauss quadrature formula (see notes)

## Bess::Bess<a name=constructor></a>
Constructs an instance for calculations.

```C++
Bess(const int n = 34, const int m = 34, const double d = 2.);
```
Paramterers `n` and `m` define the number of members in Chebyshev expansion for K0, parameter `d` defines lower bound of the expansion, as Chebyshev expansion of K0 (and its integral) is defined only for x >= d, d > 0. Default values are chosen as a tradeoff between accuracy and speed of calculations. Is is checked that with default parameters for 2 <= x < Inf relative error for integral calculation does not exceed 1e-14. Decreasing `n` and `m` speeds up calculation, but affects accuracy (note that `n >= m`). Decreasing `d` expands Chebyshev approximation domain but leads to a lower accuracy with `n` and `m` fixed. For x < d approximations other than Chebyshev are used (polynomial and Gauss quadrature integration).  
It is recommended to use only one instance of `Bess` for calculations as all member fucntions have `const` qualifiers and are thread-safe.  

## Bess::ik00x<a name=ik00x></a>

Calculates an integral of K0 in \[0,x\] with relative error < 1e-14.

```C++
double Bess::ik00x(const double x) const;
```

## Bess::ik00x_ch<a name=ik00x_ch></a>

Calculates an integral of K0 in \[0,x\], x >= d (see constructor) with relative error < 1e-14 using Chebyshev expansion.

```C++
double Bess::ik00x_ch(const double x) const;
```

## Bess::ik00x_pwr<a name=ik00x_pwr></a>

Calculates an integral of K0 in \[0,x\], x < 20 with relative error < 1e-14 using power series expansion.

```C++
double Bess::ik00x_pwr(const double x) const;
```

## Bess::ik0ab<a name=ik0ab></a>

Calculates an integral of K0 in \[a,b\], 0 <= a < b < Inf with relative error < 1e-14 using either Chebyshev or power appriximations or Gauss quardature formula.

```C++
double Bess::ik0ab(const double a, const double b) const;
```

## Bess::ik0ab_ch<a name=ik0ab_ch></a>

Calculates an integral of K0 in \[a,b\], [d](#constructor) <= a < b < Inf with relative error < 1e-14 using either Chebyshev approximation.

```C++
double Bess::ik0ab_ch(const double a, const double b) const;
```


## Bess::ik0ab_pwr<a name=ik0ab_pwr></a>

Calculates an integral of K0 in \[a,b\], 0 <= a < b < 20 with relative error < 1e-14 using power appriximations.

```C++
long double Bess::ik0ab_pwr(const double a, const double b) const;
```

## Bess::ik0ab_num<a name=ik0ab_num></a>

Calculates an integral of K0 in \[a,b\], 0 <= a < b < Inf with Gauss 10-point quardature formula. Note that high precision is only acquired whith a > 0.1, (b-a) <= 0.5 

```C++
double Bess::ik0ab_num(const double a, const double b) const;
```


## License
[MIT](https://choosealicense.com/licenses/mit/)