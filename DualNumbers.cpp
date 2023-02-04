#ifndef ALS_MATH_DUALNUMBERS_CPP
#define ALS_MATH_DUALNUMBERS_CPP

#include "DualNumbers.hpp"
#include <cmath>

template <typename K>
const als::math::DualNumber<K> als::math::DualNumber<K>::epsilon = DualNumber(0, 1);

template <typename K>
als::math::DualNumber<K>::DualNumber() : a(0), b(0)
{

}

template <typename K>
als::math::DualNumber<K>::DualNumber(const K a) : a(a), b(0)
{
    
}

template <typename K>
als::math::DualNumber<K>::DualNumber(const K a, const K b) : a(a), b(b)
{
    
}

template <typename K>
als::math::DualNumber<K> als::math::operator+(const als::math::DualNumber<K>& x, const als::math::DualNumber<K>& y)
{
    return als::math::DualNumber<K>(x.a + y.a, x.b + y.b);
}

template <typename K>
als::math::DualNumber<K> als::math::operator+(const als::math::DualNumber<K>& x, const K y)
{
    return x + als::math::DualNumber<K>(y);
}

template <typename K>
als::math::DualNumber<K> als::math::operator+(const K x, const als::math::DualNumber<K>& y)
{
    return als::math::DualNumber<K>(x) + y;
}

template <typename K>
als::math::DualNumber<K> als::math::operator-(const als::math::DualNumber<K>& x, const als::math::DualNumber<K>& y)
{
    return als::math::DualNumber<K>(x.a - y.a, x.b - y.b);
}

template <typename K>
als::math::DualNumber<K> als::math::operator-(const als::math::DualNumber<K>& x, const K y)
{
    return x - als::math::DualNumber<K>(y);
}

template <typename K>
als::math::DualNumber<K> als::math::operator-(const K x, const als::math::DualNumber<K>& y)
{
    return als::math::DualNumber<K>(x) - y;
}

template <typename K>
als::math::DualNumber<K> als::math::operator*(const als::math::DualNumber<K>& x, const als::math::DualNumber<K>& y)
{
    return als::math::DualNumber<K>(x.a * y.a, x.a * y.b + x.b * y.a);
}

template <typename K>
als::math::DualNumber<K> als::math::operator*(const als::math::DualNumber<K>& x, const K y)
{
    return x * als::math::DualNumber<K>(y);
}

template <typename K>
als::math::DualNumber<K> als::math::operator*(const K x, const als::math::DualNumber<K>& y)
{
    return als::math::DualNumber<K>(x) * y;
}

template <typename K>
als::math::DualNumber<K> als::math::operator/(const als::math::DualNumber<K>& x, const als::math::DualNumber<K>& y)
{
    return als::math::DualNumber<K>(x.a / y.a, (x.b * y.a - x.a * y.b) / y.a / y.a);
}

template <typename K>
als::math::DualNumber<K> als::math::operator/(const als::math::DualNumber<K>& x, const K y)
{
    return x / als::math::DualNumber<K>(y);
}

template <typename K>
als::math::DualNumber<K> als::math::operator/(const K x, const als::math::DualNumber<K>& y)
{
    return als::math::DualNumber<K>(x) / y;
}

template <typename K>
template <class... Args>
std::string als::math::DualNumber<K>::to_string(
    const als::utilities::RepresentationType rt, Args... args) const
{
    std::string epsilon = 
        rt == als::utilities::RepresentationType::PLAIN ? "ε" : "\\varepsilon";    
    if (b < 0)
    {
        return als::utilities::to_string(a, rt, args...) + " - " +
            als::utilities::to_string(-b, rt, args...) + epsilon;
    }
    else
    {
        return als::utilities::to_string(a, rt, args...) + " + " +
            als::utilities::to_string(b, rt, args...) + epsilon;
    }
}

template <typename K>
als::math::DualNumber<K> als::math::cos(const als::math::DualNumber<K>& x)
{
    using ::cos;
    using ::sin;
    return als::math::DualNumber<K>(cos(x.a), - sin(x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::sin(const als::math::DualNumber<K>& x)
{
    using ::sin;
    using ::cos;
    return als::math::DualNumber<K>(sin(x.a), cos(x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::tan(const als::math::DualNumber<K>& x)
{
    using ::cos;
    using ::tan;
    // sec(a)
    K seca = 1./cos(x.a);
    return als::math::DualNumber<K>(tan(x.a), seca*seca * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::acos(const als::math::DualNumber<K>& x)
{
    using ::acos;
    using ::sqrt;
    return als::math::DualNumber<K>(acos(x.a), - 1./sqrt(1 - x.a*x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::asin(const als::math::DualNumber<K>& x)
{
    using ::asin;
    using ::sqrt;
    return als::math::DualNumber<K>(asin(x.a), 1./sqrt(1 - x.a*x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::atan(const als::math::DualNumber<K>& x)
{
    using ::atan;
    return als::math::DualNumber<K>(atan(x.a), 1./(1 + x.a*x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::cosh(const als::math::DualNumber<K>& x)
{
    using ::cosh;
    using ::sinh;
    return als::math::DualNumber<K>(cosh(x.a), sinh(x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::sinh(const als::math::DualNumber<K>& x)
{
    using ::sinh;
    using ::cosh;
    return als::math::DualNumber<K>(sinh(x.a), cosh(x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::tanh(const als::math::DualNumber<K>& x)
{
    // sech(a)
    using ::cosh;
    using ::tanh;
    K secha = 1./cosh(x.a);
    return als::math::DualNumber<K>(tanh(x.a), secha*secha * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::acosh(const als::math::DualNumber<K>& x)
{
    using ::acosh;
    using ::sqrt;
    return als::math::DualNumber<K>(acosh(x.a), 1. / sqrt(x.a*x.a - 1) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::asinh(const als::math::DualNumber<K>& x)
{
    using ::asinh;
    using ::sqrt;
    return als::math::DualNumber<K>(asinh(x.a), 1. / sqrt(x.a*x.a + 1) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::atanh(const als::math::DualNumber<K>& x)
{
    using ::atanh;
    return als::math::DualNumber<K>(atanh(x.a), 1. / (1 - x.a*x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::exp(const als::math::DualNumber<K>& x)
{
    using ::exp;
    // exp(a)
    K expa = ::exp(x.a);
    return als::math::DualNumber<K>(expa, expa * x.b);
}
template <typename K>
als::math::DualNumber<K> als::math::log(const als::math::DualNumber<K>& x)
{
    using ::log;
    return als::math::DualNumber<K>(log(x.a), 1. / x.a * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::log10(const als::math::DualNumber<K>& x)
{
    using ::log10;
    using ::log;
    return als::math::DualNumber<K>(log10(x.a), 1. / (x.a * log(10)) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::exp2(const als::math::DualNumber<K>& x)
{
    using ::exp2;
    using ::log;
    // exp2(a)
    K exp2a = exp2(x.a);
    return als::math::DualNumber<K>(exp2a, log(2) * exp2a * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::expm1(const als::math::DualNumber<K>& x)
{
    using ::expm1;
    // expm1(a)
    K expm1a = expm1(x.a);
    return als::math::DualNumber<K>(expm1a, (1 + expm1a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::log1p(const als::math::DualNumber<K>& x)
{
    using ::log1p;
    return als::math::DualNumber<K>(log1p(x.a), 1. / (1 + x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::log2(const als::math::DualNumber<K>& x)
{
    using ::log2;
    using ::log;
    return als::math::DualNumber<K>(log2(x.a), 1. / (x.a * log(2)) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::logb(const als::math::DualNumber<K>& x)
{
    using ::logb;
    using ::log;
    return als::math::DualNumber<K>(logb(x.a), 1. / (x.a * log(2)) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::pow(const als::math::DualNumber<K>& x, const als::math::DualNumber<K>& y)
{
    using ::pow;
    using ::log;
    K powxy = pow(x.a, y.a);
    return als::math::DualNumber<K>(powxy, powxy * (y.b*log(x.a) + y.a*x.b/x.a));
}

template <typename K>
als::math::DualNumber<K> als::math::pow(const als::math::DualNumber<K>& x, const K& y)
{
    using ::pow;
    K powxy = pow(x.a, y);
    return als::math::DualNumber<K>(powxy, powxy * y*x.b/x.a);
}

template <typename K>
als::math::DualNumber<K> als::math::pow(const K& x, const als::math::DualNumber<K>& y)
{
    using ::pow;
    using ::log;
    K powxy = pow(x, y.a);
    return als::math::DualNumber<K>(powxy, powxy * y.b*log(x));
}

template <typename K>
als::math::DualNumber<K> als::math::sqrt(const als::math::DualNumber<K>& x)
{
    using ::sqrt;
    // sqrt(a)
    K sqrta = sqrt(x.a);
    return als::math::DualNumber<K>(sqrta, 1. / (2 * sqrta) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::cbrt(const als::math::DualNumber<K>& x)
{
    using ::cbrt;
    return als::math::DualNumber<K>(cbrt(x.a), 1. / (3 * cbrt(x.a*x.a)) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::erf(const als::math::DualNumber<K>& x)
{
    using ::erf;
    using ::sqrt;
    using ::exp;
    return als::math::DualNumber<K>(erf(x.a), 2./sqrt(M_PI)*exp(-x.a*x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::erfc(const als::math::DualNumber<K>& x)
{
    using ::erfc;
    using ::sqrt;
    using ::exp;
    return als::math::DualNumber<K>(erfc(x.a), -2./sqrt(M_PI)*exp(-x.a*x.a) * x.b);
}

template <typename K>
als::math::DualNumber<K> als::math::fabs(const als::math::DualNumber<K>& x)
{
    using ::fabs;
    return als::math::DualNumber<K>(fabs(x.a), (x.a > 0) ? x.b : -x.b);
}


#endif // ALS_MATH_DUALNUMBERS_CPP