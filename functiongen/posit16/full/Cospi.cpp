#include "PolySynth.h"
#include "Elementary.h"

template <>
posit16 Elementary<posit16>::MpfrCalculateFunction(posit16 x) {
    // Special case that MPFR can't handle:
    double dx = x.toDouble();
    double integral, frac;
    frac = modf(dx, &integral);
    if (frac == 0.5 || frac == -0.5) {
        return 0.0;
    }
    mpfr_const_pi(mval, MPFR_RNDN);
    mpfr_mul_d(mval, mval, dx, MPFR_RNDN);
    mpfr_cos(mval, mval, MPFR_RNDN);
    return Elementary<posit16>::FromMPFR(mval);
}

template <>
bool Elementary<posit16>::ComputeSpecialCase(posit16 x, posit16& res) {
    if (x.value == 0x8000) {
        // Take care of NaN
        res.value = 0x8000;
        return true;
    }
    
    return false;
}

template <>
double Elementary<posit16>::RangeReduction(posit16 x, double& modifier) {
    modifier = 1;

    // cos(-x) = cos(x)
    if (x < 0.0f) x = -1 * x.toDouble();
    
    // How do we reduce range of x?
    // Reduce x to [0, 1)
    double intPart;
    double frac = modf(x.toDouble(), &intPart);
    int iIntPart = intPart;
    
    // if iIntPart is odd, then flip modifier
    if (iIntPart % 2 == 1) modifier *= -1;
    
    // cos(pi - x) = -cos(x)
    if (frac >= 0.5) {
        frac = 1.0 - frac;
        modifier *= -1;
    }
    
    return frac;
}

template <class T>
double Elementary<T>::RangePropagation(double y, double modifier) {
    return y * modifier;
}

template <class T>
double Elementary<T>::ReverseRangePropagation(double y, double modifier) {
    return y * modifier;
}

template <class T>
bool Elementary<T>::FlipLbAndUb(double modifier) {
    return modifier != 1.0;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 2000);
    
    PolySynth<posit16, Elementary<posit16>> p16cospi;
    p16cospi.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    p16cospi.FindPolynomialOnce({0});
    p16cospi.FindPolynomialsUntil({0, 2, 4, 6, 8}, 12287);
    p16cospi.FindPolynomials({0});
    printf("\tCOMPLETED\n\n");
    
    p16cospi.poly->PrintPiecewiseInfo();
    p16cospi.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}
