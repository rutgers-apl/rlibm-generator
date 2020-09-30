#include "PolySynth.h"
#include "Elementary.h"

template <>
bfloat16 Elementary<bfloat16>::MpfrCalculateFunction(bfloat16 x) {
    mpfr_set_d(mval, (double)x, MPFR_RNDN);
    mpfr_exp10(mval, mval, MPFR_RNDN);
    return Elementary<bfloat16>::FromMPFR(mval);
}

template <>
bool Elementary<bfloat16>::ComputeSpecialCase(bfloat16 x, bfloat16& res) {
    if ((float)x <= -40.5f) {
        // Take care of when result is 0
        res.val = 0;
        return true;
    } else if ((float)x >= 38.75f) {
        // Take care of infinity case
        res.val = 0x7F80;
        return true;
    } else if ((x.val & 0x7FFF) > 0x7F80) {
        // Take care of NaN
        res.val = x.val;
        return true;
    } else if ((x >= -8.4686279296875e-04) && (x <= 1.68609619140625e-03)) {
        // The values in these range return 1.0.
        //printf("%.100e\n", (float)x);
        res = 1.0f;
        return true;
    }
    
    return false;
}

template <>
double Elementary<bfloat16>::RangeReduction(bfloat16 x, double& modifier) {
    double xprime = (double)x * 3.321928094887362181708567732130177319049835205078125;
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
bool Elementary<T>::Elementary<T>::FlipLbAndUb(double y) {
    return false;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 2000);
    
    PolySynth<bfloat16, Elementary<bfloat16>> bf16exp10;
    bf16exp10.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    bf16exp10.FindPolynomials({0, 1, 2, 3, 4});
    printf("\tCOMPLETED\n\n");
    
    bf16exp10.poly->PrintPiecewiseInfo();
    bf16exp10.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}
