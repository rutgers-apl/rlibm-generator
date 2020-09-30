#include "PolySynth.h"
#include "Elementary.h"

template <>
posit16 Elementary<posit16>::MpfrCalculateFunction(posit16 x) {
    mpfr_set_d(mval, x.toDouble(), MPFR_RNDN);
    mpfr_exp2(mval, mval, MPFR_RNDN);
    posit16 retVal = Elementary<posit16>::FromMPFR(mval);
    return retVal;
}

template <>
bool Elementary<posit16>::ComputeSpecialCase(posit16 x, posit16& res) {
    if (x.value > 0x8000 & x.value <= 0x8d3f) {
        // Take care of when result is minpos. exp(x) for posit should never
        // return 0, because exp(x) is always > 0 as long as x != -infinity
        res.value = 0x1;
        return true;
    } else if (x.value >= 0x72c1 && x.value < 0x8000) {
        // Take care of maxpos case.
        res.value = 0x7FFF;
        return true;
    } else if (x.value == 0x8000) {
        // Take care of NaR
        res.value = 0x8000;
        return true;
    } else if ((x.value >= 0xff64) || (x.value <= 0xdc)) {
        // The values in these range return 1.0.
        res = 1.0f;
        return true;
    }
    
    return false;
}

template <>
double Elementary<posit16>::RangeReduction(posit16 x, double& modifier) {
    double xprime = x.toDouble();
    double f = modf(xprime, &modifier);
    if (f < 0.0f) {
        f += 1.0;
        modifier -= 1.0;
    }
    return f;
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
    
    PolySynth<posit16, Elementary<posit16>> p16exp(SAMPLE, 200);
    p16exp.SetUpSamplingOption(0, 1);
    p16exp.CalcIntervals();
    
    do {
        printf("FINDING POLYNOMIALS\n");
        p16exp.ResetPolynomials();
        p16exp.FindPolynomials({0, 1, 2, 3, 4, 5, 6});
        printf("\tCOMPLETED\n\n");
    } while (p16exp.TestAndAddSamplePoints());
    
    p16exp.poly->PrintPiecewiseInfo();
    p16exp.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}
