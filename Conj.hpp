#ifndef ALS_MATH_CONJ_HPP
#define ALS_MATH_CONJ_HPP

#include <complex>

namespace als::math
{
    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return char 
     */
    char inline conj(const char& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned char 
     */
    unsigned char inline conj(const unsigned char& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return short int 
     */
    short int inline conj(const short int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned short int 
     */
    unsigned short int inline conj(const unsigned short int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return int 
     */
    int inline conj(const int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned int 
     */
    unsigned int inline conj(const unsigned int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return long int 
     */
    long int inline conj(const long int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned long int 
     */
    unsigned long int inline conj(const unsigned long int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return long long int 
     */
    long long int inline conj(const long long int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned long long int 
     */
    unsigned long long int inline conj(const unsigned long long int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return float 
     */
    float inline conj(const float& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return double 
     */
    double inline conj(const double& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return long double 
     */
    long double inline conj(const long double& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return wchar_t 
     */
    wchar_t inline conj(const wchar_t& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns the complex conjugate of @a val .
     * 
     * @tparam K 
     * @param val 
     * @return std::complex<K> 
     */
    template<typename K>
    std::complex<K> inline conj(const std::complex<K>& val)
    {
        return std::conj(val);
    }
    
    template<typename K>
    std::string inline to_string(const std::complex<K>& z)
    {
        using std::to_string;
        return to_string(z.real()) + " + " + to_string(z.imag()) + "i"; 
    }

}

#endif // ALS_MATH_CONJ_HPP