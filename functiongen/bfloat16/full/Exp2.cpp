#include "PolySynth.h"
#include "Elementary.h"

template <>
bfloat16 Elementary<bfloat16>::MpfrCalculateFunction(bfloat16 x) {
    mpfr_set_d(mval, (double)x, MPFR_RNDN);
    mpfr_exp2(mval, mval, MPFR_RNDN);
    return Elementary<bfloat16>::FromMPFR(mval);
}

template <>
bool Elementary<bfloat16>::ComputeSpecialCase(bfloat16 x, bfloat16& res) {
    if ((float)x <= -134.0f) {
        res.val = 0;
        return true;
    } else if ((float)x >= 128) {
        // Take care of infinity case
        res.val = 0x7F80;
        return true;
    } else if ((x.val & 0x7FFF) > 0x7F80) {
        // Take care of NaN
        res.val = x.val;
        return true;
    }
    
    // Take care of when result is 1 :
    if (x >= -2.8076171875e-03 && x <= 2.8076171875e-03) {
        res = 1.0;
        return true;
    }
    
    return false;
}

template <>
double Elementary<bfloat16>::RangeReduction(bfloat16 x, double& modifier) {
    double f = modf((float)x, &modifier);
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
bool Elementary<T>::FlipLbAndUb(double y) {
    return false;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 2000);
    
    PolySynth<bfloat16, Elementary<bfloat16>> bf16exp2;
    bf16exp2.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    bf16exp2.FindPolynomials({0, 1, 2, 3, 4});
    printf("\tCOMPLETED\n\n");
    
    bf16exp2.poly->PrintPiecewiseInfo();
    bf16exp2.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}

