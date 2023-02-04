#ifndef ALS_MATH_INTEGRATION_HPP
#define ALS_MATH_INTEGRATION_HPP

#include <functional>

namespace als::math::integrators
{
    double rectangle_rule(const std::function<double(double x)>& f, const double a, const double b,
        const unsigned int N);

    double trapezoidal_rule(const std::function<double(double x)>& f, const double a, const double b,
        const unsigned int N);

    double Gauss_Konrad_G7_K15(const std::function<double(double x)>& f, const double a, const double b,
        const double tol = 1e-6);
}

#endif // ALS_MATH_INTEGRATION_HPP