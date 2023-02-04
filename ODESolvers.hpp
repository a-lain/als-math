#ifndef ALS_MATH_ODE_SOLVERS_HPP
#define ALS_MATH_ODE_SOLVERS_HPP

#include "Vector.hpp"
#include <functional>

namespace als::math::ode_solvers
{
    Vector<double> euler_explicit(const Vector<double>& x, const double t,
        std::function<Vector<double>(const Vector<double>&, const double)> f, const double dt);
    Vector<double> euler_explicit(const Vector<double>& x,
        std::function<Vector<double>(const Vector<double>&)> f, const double dt);

    Vector<double> runge_kutta_order_4(const Vector<double>& x, const double t,
        std::function<Vector<double>(const Vector<double>&, const double)> f, const double dt);
    Vector<double> runge_kutta_order_4(const Vector<double>& x,
        std::function<Vector<double>(const Vector<double>&)> f, const double dt);

    void runge_kutta_fehlberg(Vector<double>& x, double& t, double& dt, const std::function<Vector<double>(const Vector<double>&, const double)>& f,
        const Vector<double>& abs_tol);
    void runge_kutta_fehlberg(Vector<double>& x, double& t, double& dt, const std::function<Vector<double>(const Vector<double>&)>& f,
        const Vector<double>& abs_tol);
}

#endif // ALS_MATH_ODE_SOLVERS_HPP