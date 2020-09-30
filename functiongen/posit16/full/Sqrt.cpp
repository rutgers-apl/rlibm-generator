#include "PolySynth.h"
#include "Elementary.h"

template <>
posit16 Elementary<posit16>::MpfrCalculateFunction(posit16 x) {
    mpfr_set_d(mval, x.toDouble(), MPFR_RNDN);
    mpfr_sqrt(mval, mval, MPFR_RNDN);
    return Elementary<posit16>::FromMPFR(mval);
}

template <>
bool Elementary<posit16>::ComputeSpecialCase(posit16 x, posit16& res) {
    if (x.value == 0) {
        res.value = 0x0;
        return true;
    }
    if (x.value >= 0x8000) {
        // If x <= 0 or x is NaR, it should be NaR
        res.value = 0x8000;
        return true;
    }
    
    return false;
}

template <>
double Elementary<posit16>::RangeReduction(posit16 x, double& modifier) {
    // Extract exponent and mantissa
    int m;
    float fx = frexpf(x.toDouble(), &m);
    
    if ((m & 0x1) == 0) {
        fx *= 4.0;
        m -= 2;
    } else {
        fx *= 2.0;
        m--;
    }
    
    modifier = (m >> 1);
    
    return fx;
}

template <class T>
double Elementary<T>::RangePropagation(double y, double modifier) {
    return ldexp(y, modifier);
}

template <class T>
double Elementary<T>::ReverseRangePropagation(double y, double modifier) {
    return ldexp(y, -modifier);
}

template <class T>
bool Elementary<T>::FlipLbAndUb(double modifier) {
    return false;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 2000);
    
    PolySynth<posit16, Elementary<posit16>> p16sqrt;
    p16sqrt.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    p16sqrt.FindPolynomials({0, 1, 2, 3, 4, 5, 6});
    printf("\tCOMPLETED\n\n");
    
    p16sqrt.poly->PrintPiecewiseInfo();
    p16sqrt.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}
