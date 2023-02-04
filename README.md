# als-math
A mathematical library for C++. The file "Math.ipynb" contains good examples of how to use the library. In order to run it, you need to install [als-xeus-cling](https://github.com/a-lain/als-xeus-cling). For the moment, the library supports
- vector arithmetic (Vector.hpp and Vector.cpp),
- automatic differentiation (DualNumbers.hpp and DualNumbers.hpp),
- some matrix operations (Matrix.hpp and Matrix.cpp),
- bisection method, secant method, Newton-Raphson method, LU decomposition to solve linear systems (AlgebraicSolvers.hpp and AlgebraicSolvers.cpp),
- rectangle rule, trapezoide rule, Gauss-Konrad G7-K15 integration methods (Integration.hpp and Integration.cpp),
- average linear and quadratic interpolations, linear interpolation (Interpolation.hpp and Interpolation.cpp),
- euler explicit mehtod, Runge-Kutta method of order 4, Runge-Kutta-Fehlberg method to solve ODEs (ODESolvers.hpp and ODESolvers.cpp).

## Dependencies:
- [als-basic-utilities](https://github.com/a-lain/als-basic-utilities)
- (Optional, but recommended) [als-xeus-cling](https://github.com/a-lain/als-xeus-cling)

## Installation:
- Linux: adapt the contents of the Makefile to match the configuration of your system.
- Windows: you are on your own.
