#include "../Vector.cpp"
#include "../Matrix.cpp"
#include "../DualNumbers.cpp"

namespace als::math
{
    // Vector.
    template class Vector<double>;

    template Vector<double> operator+(const Vector<double>& v, const Vector<double>& w);
    template Vector<double> operator+(const Vector<double>& v, const double& alpha);
    template Vector<double> operator+(const double& alpha, const Vector<double>& v);
    template Vector<double> operator-(const Vector<double>& v, const Vector<double>& w);
    template Vector<double> operator-(const Vector<double>& v, const double& alpha);
    template Vector<double> operator-(const double& alpha, const Vector<double>& v);
    template Vector<double> operator*(const Vector<double>& v, const Vector<double>& w);
    template Vector<double> operator*(const Vector<double>& v, const double& alpha);
    template Vector<double> operator*(const double& alpha, const Vector<double>& v);
    template Vector<double> operator/(const Vector<double>& v, const Vector<double>& w);
    template Vector<double> operator/(const Vector<double>& v, const double& alpha);
    template Vector<double> operator/(const double& alpha, const Vector<double>& v);
    template bool operator<(const Vector<double>& v, const double& alpha);
    template bool operator<(const double& alpha, const Vector<double>& v);
    template bool operator<(const Vector<double>& v, const Vector<double>& w); 
    template bool operator<=(const Vector<double>& v, const double& alpha);
    template bool operator<=(const double& alpha, const Vector<double>& v);
    template bool operator<=(const Vector<double>& v, const Vector<double>& w);
    template bool operator>(const Vector<double>& v, const double& alpha);
    template bool operator>(const double& alpha, const Vector<double>& v);
    template bool operator>(const Vector<double>& v, const Vector<double>& w);
    template bool operator>=(const Vector<double>& v, const double& alpha);
    template bool operator>=(const double& alpha, const Vector<double>& v);
    template bool operator>=(const Vector<double>& v, const Vector<double>& w);
    template bool operator==(const Vector<double>& v, const Vector<double>& w);
    template bool operator==(const Vector<double>& v, const double& alpha);
    template bool operator==(const double& alpha, const Vector<double>& v);
    template bool operator!=(const Vector<double>& v, const Vector<double>& w);
    template bool operator!=(const Vector<double>& v, const double& alpha);
    template bool operator!=(const double& alpha, const Vector<double>& v);
    template Vector<double> operator&(const Vector<double>& v, const Vector<double>& w);
    template Vector<double> operator&(const Vector<double>& v, const double& alpha);
    template Vector<double> operator&(const double& alpha, const Vector<double>& v);
    template double operator|(const Vector<double>& v, const Vector<double>& w);    
    template double vector_product_2d(const Vector<double>& v, const Vector<double>& w);
    template Vector<double> vector_product_3d(const Vector<double>& v, const Vector<double>& w);
    template std::string Vector<double>::to_string(const als::utilities::RepresentationType rt) const;
    template std::string Vector<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision) const;
    template std::string Vector<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign) const;
    template std::string Vector<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign, const int lim_inf) const;
    template std::string Vector<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign, const int lim_inf, const int lim_sup) const;
    template std::ostream& als::math::operator<< <double>(std::ostream&os, const Vector<double>& v);    

    template double min(const Vector<double>& v);
    template double max(const Vector<double>& v);
    template double sum(const Vector<double>& v);
    template double multiply(const Vector<double>& v);

    template Vector<double> cos(const Vector<double>& v);    
    template Vector<double> sin(const Vector<double>& v);    
    template Vector<double> tan(const Vector<double>& v);
    template Vector<double> acos(const Vector<double>& v);    
    template Vector<double> asin(const Vector<double>& v);    
    template Vector<double> atan(const Vector<double>& v);    
    template Vector<double> atan2(const double v, const Vector<double>& w);    
    template Vector<double> atan2(const Vector<double>& v, const double w);    
    template Vector<double> atan2(const Vector<double>& v, const Vector<double>& w);    
    template Vector<double> cosh(const Vector<double>& v);    
    template Vector<double> sinh(const Vector<double>& v);    
    template Vector<double> tanh(const Vector<double>& v);    
    template Vector<double> acosh(const Vector<double>& v);    
    template Vector<double> asinh(const Vector<double>& v);    
    template Vector<double> atanh(const Vector<double>& v);    
    template Vector<double> exp(const Vector<double>& v);    
    template Vector<double> frexp(const Vector<double>& v, Vector<int>* exp);    
    template Vector<double> ldexp(const Vector<double>& v, const int exp);    
    template Vector<double> ldexp(const Vector<double>& v, const Vector<int>& exp);    
    template Vector<double> log(const Vector<double>& v);    
    template Vector<double> log10(const Vector<double>& v);    
    template Vector<double> modf(const Vector<double>& v, Vector<double>* intpart);    
    template Vector<double> exp2(const Vector<double>& v);    
    template Vector<double> expm1(const Vector<double>& v);    
    template Vector<double> ilogb(const Vector<double>& v);    
    template Vector<double> log1p(const Vector<double>& v);    
    template Vector<double> log2(const Vector<double>& v);    
    template Vector<double> logb(const Vector<double>& v);    
    template Vector<double> scalbn(const Vector<double>& v, const int n);    
    template Vector<double> scalbn(const Vector<double>& v, const Vector<int>& n);    
    template Vector<double> scalbln(const Vector<double>& v, const long int n);    
    template Vector<double> pow(const double v, const Vector<double>& exponent);    
    template Vector<double> pow(const Vector<double>& v, const double exponent);    
    template Vector<double> pow(const Vector<double>& v, const Vector<double>& exponent);    
    template Vector<double> sqrt(const Vector<double>& v);    
    template Vector<double> cbrt(const Vector<double>& v);    
    template Vector<double> hypot(const double x, const Vector<double>& y);    
    template Vector<double> hypot(const Vector<double>& x, const double y);    
    template Vector<double> hypot(const Vector<double>& x, const Vector<double>& y);    
    template Vector<double> erf(const Vector<double>& v);    
    template Vector<double> erfc(const Vector<double>& v);    
    template Vector<double> tgamma(const Vector<double>& v);    
    template Vector<double> lgamma(const Vector<double>& v);    
    template Vector<double> ceil(const Vector<double>& v);    
    template Vector<double> floor(const Vector<double>& v);    
    template Vector<double> fmod(const double numer, const Vector<double>& denom);    
    template Vector<double> fmod(const Vector<double>& numer, const double denom);    
    template Vector<double> fmod(const Vector<double>& numer, const Vector<double>& denom);    
    template Vector<double> trunc(const Vector<double>& v);    
    template Vector<double> round(const Vector<double>& v);    
    template Vector<long int> lround(const Vector<double>& v);    
    template Vector<long long int> llround(const Vector<double>& v);    
    template Vector<double> rint(const Vector<double>& v);    
    template Vector<long int> lrint(const Vector<double>& v);
    template Vector<long long int> llrint(const Vector<double>& v);    
    template Vector<double> nearbyint(const Vector<double>& v);    
    template Vector<double> remainder(const double numer, const Vector<double>& denom);    
    template Vector<double> remainder(const Vector<double>& numer, const double denom);    
    template Vector<double> remainder(const Vector<double>& numer, const Vector<double>& denom);    
    template Vector<double> remquo(const double numer, const Vector<double>& denom, Vector<int>* quot);    
    template Vector<double> remquo(const Vector<double>& numer, const double denom, Vector<int>* quot);    
    template Vector<double> remquo(const Vector<double>& numer, const Vector<double>& denom, Vector<int>* quot);    
    template Vector<double> copysign(const Vector<double>& x, const double y);    
    template Vector<double> copysign(const Vector<double>& x, const Vector<double>& y);    
    template Vector<double> nan(const unsigned int N, const char* tagp);
    template Vector<double> nextafter(const Vector<double>& x, const Vector<double>& y);
    template Vector<double> fdim(const double x, const Vector<double>& y);    
    template Vector<double> fdim(const Vector<double>& x, double y);    
    template Vector<double> fdim(const Vector<double>& x, const Vector<double>& y);    
    template Vector<double> fabs(const Vector<double>& v);    
    template Vector<double> abs(const Vector<double>& v);    
    template Vector<double> fma(const double x, const Vector<double>& y, const Vector<double>& z);    
    template Vector<double> fma(const Vector<double>& x, const double y, const Vector<double>& z);    
    template Vector<double> fma(const Vector<double>& x, const Vector<double>& y, const double z);    
    template Vector<double> fma(const double x, const double y, const Vector<double>& z);    
    template Vector<double> fma(const double x, const Vector<double>& y, const double z);    
    template Vector<double> fma(const Vector<double>& x, const double y, const double z);    
    template Vector<double> fma(const Vector<double>& x, const Vector<double>& y, const Vector<double>& z);    
    template Vector<int> fpclassify(const Vector<double>& v);    
    template bool isfinite(const Vector<double>& v);    
    template bool isinf(const Vector<double>& v);    
    template bool isnan(const Vector<double>& v);    
    template bool isnormal(const Vector<double>& v);    

    // Matrix.
    template class Matrix<double>;
    template Matrix<double> operator+(const Matrix<double>& A, const Matrix<double>& B);
    template Matrix<double> operator+(const Matrix<double>& A, const double& alpha);
    template Matrix<double> operator+(const double& alpha, const Matrix<double>& A);
    template Matrix<double> operator-(const Matrix<double>& A, const Matrix<double>& B);
    template Matrix<double> operator-(const Matrix<double>& A, const double& alpha);
    template Matrix<double> operator-(const double& alpha, const Matrix<double>& A);
    template Matrix<double> operator*(const Matrix<double>& A, const double& alpha);
    template Matrix<double> operator*(const double& alpha, const Matrix<double>& A);
    template Matrix<double> operator/(const Matrix<double>& A, const double& alpha);
    template Matrix<double> operator|(const Matrix<double>& A, const Matrix<double>& B);
    template Vector<double> operator|(const Matrix<double>& A, const Vector<double>& v);
    template Vector<double> operator|(const Vector<double>& v, const Matrix<double>& A);

    template std::string Matrix<double>::to_string(const als::utilities::RepresentationType rt) const;
    template std::string Matrix<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision) const;
    template std::string Matrix<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign) const;
    template std::string Matrix<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign, const int lim_inf) const;
    template std::string Matrix<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign, const int lim_inf, const int lim_sup) const;

    // Dual numbers.
    template class DualNumber<double>;

    template std::string DualNumber<double>::to_string(const als::utilities::RepresentationType rt) const;
    template std::string DualNumber<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision) const;
    template std::string DualNumber<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign) const;
    template std::string DualNumber<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign, const int lim_inf) const;
    template std::string DualNumber<double>::to_string(const als::utilities::RepresentationType rt,
        const unsigned int precision, const bool show_sign, const int lim_inf, const int lim_sup) const;

    template DualNumber<double> operator+(const DualNumber<double>& x, const DualNumber<double>& y);
    template DualNumber<double> operator+(const DualNumber<double>& x, const double y);
    template DualNumber<double> operator+(const double x, const DualNumber<double>& y);
    template DualNumber<double> operator-(const DualNumber<double>& x, const DualNumber<double>& y);
    template DualNumber<double> operator-(const DualNumber<double>& x, const double y);
    template DualNumber<double> operator-(const double x, const DualNumber<double>& y);
    template DualNumber<double> operator*(const DualNumber<double>& x, const DualNumber<double>& y);
    template DualNumber<double> operator*(const DualNumber<double>& x, const double y);
    template DualNumber<double> operator*(const double x, const DualNumber<double>& y);
    template DualNumber<double> operator/(const DualNumber<double>& x, const DualNumber<double>& y);
    template DualNumber<double> operator/(const DualNumber<double>& x, const double y);
    template DualNumber<double> operator/(const double x, const DualNumber<double>& y);

    template DualNumber<double> cos(const DualNumber<double>& x);
    template DualNumber<double> sin(const DualNumber<double>& x);
    template DualNumber<double> tan(const DualNumber<double>& x);
    template DualNumber<double> acos(const DualNumber<double>& x);
    template DualNumber<double> asin(const DualNumber<double>& x);
    template DualNumber<double> atan(const DualNumber<double>& x);
    template DualNumber<double> cosh(const DualNumber<double>& x);
    template DualNumber<double> sinh(const DualNumber<double>& x);
    template DualNumber<double> tanh(const DualNumber<double>& x);
    template DualNumber<double> acosh(const DualNumber<double>& x);
    template DualNumber<double> asinh(const DualNumber<double>& x);
    template DualNumber<double> atanh(const DualNumber<double>& x);
    template DualNumber<double> exp(const DualNumber<double>& x);
    template DualNumber<double> log(const DualNumber<double>& x);
    template DualNumber<double> log10(const DualNumber<double>& x);
    template DualNumber<double> exp2(const DualNumber<double>& x);
    template DualNumber<double> expm1(const DualNumber<double>& x);
    template DualNumber<double> log1p(const DualNumber<double>& x);
    template DualNumber<double> log2(const DualNumber<double>& x);
    template DualNumber<double> logb(const DualNumber<double>& x);
    template DualNumber<double> pow(const DualNumber<double>& x, const DualNumber<double>& y);
    template DualNumber<double> pow(const DualNumber<double>& x, const double& y);
    template DualNumber<double> pow(const double& x, const DualNumber<double>& y);
    template DualNumber<double> sqrt(const DualNumber<double>& x);
    template DualNumber<double> cbrt(const DualNumber<double>& x);
    template DualNumber<double> erf(const DualNumber<double>& x);
    template DualNumber<double> erfc(const DualNumber<double>& x);
    // tgamma and lgamma missing because they are not implemented.
    template DualNumber<double> fabs(const DualNumber<double>& x);
}
