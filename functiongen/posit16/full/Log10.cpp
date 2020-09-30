#include "PolySynth.h"
#include "Elementary.h"

template <>
posit16 Elementary<posit16>::MpfrCalculateFunction(posit16 x) {
    mpfr_set_d(mval, x.toDouble(), MPFR_RNDN);
    mpfr_log10(mval, mval, MPFR_RNDN);
    return Elementary<posit16>::FromMPFR(mval);
}

template <>
bool Elementary<posit16>::ComputeSpecialCase(posit16 x, posit16& res) {
    if (x.value == 0x0 || x.value >= 0x8000) {
        // If x == 0, NaR, or negative, then resutl should be NaR
        res.value = 0x8000;
        return true;
    }
    
    return false;
}

template <>
double Elementary<posit16>::RangeReduction(posit16 x, double& modifier) {
    // Extract exponent and mantissa (where mantissa is between [1, 2))
    int m;
    float fx = frexpf(x.toDouble(), &m);
    fx *= 2.0;
    m--;
    
    // Cody and Waite Transformation on input
    double dx = (double)fx;
    double codyX = (dx - 1) / (dx + 1);
    modifier = m;
    return codyX;
}

template <class T>
double Elementary<T>::RangePropagation(double y, double modifier) {
    return (y + modifier) / 3.321928094887362181708567732130177319049835205078125;
}

template <class T>
double Elementary<T>::ReverseRangePropagation(double y, double modifier) {
    return y * 3.321928094887362181708567732130177319049835205078125 - modifier;
}

template <class T>
bool Elementary<T>::FlipLbAndUb(double modifier) {
    return false;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 2000);
    
    PolySynth<posit16, Elementary<posit16>> p16log10;
    p16log10.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    p16log10.FindPolynomials({1, 3, 5, 7, 9});
    printf("\tCOMPLETED\n\n");
    
    p16log10.poly->PrintPiecewiseInfo();
    p16log10.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}

