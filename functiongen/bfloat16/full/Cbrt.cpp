#include "PolySynth.h"
#include "Elementary.h"

template <>
bfloat16 Elementary<bfloat16>::MpfrCalculateFunction(bfloat16 x) {
    mpfr_set_d(mval, (float)x, MPFR_RNDN);
    mpfr_cbrt(mval, mval, MPFR_RNDN);
    return Elementary<bfloat16>::FromMPFR(mval);
}

template <>
bool Elementary<bfloat16>::ComputeSpecialCase(bfloat16 x, bfloat16& res) {
    if (x.val == 0x0 || x.val == 0x8000) {
        // If x == 0, then it should be -inf
        res.val = x.val;
        return true;
    } else if (x.val == 0x7f80) {
        // If x == inf, then it should be infinity
        res = x;
        return true;
    } else if (x.val > 0x7F80) {
        // If x == NaN or negative, then it should be NaN
        res.val = 0xFFFF;
        return true;
    }
    
    return false;
}

template <>
double Elementary<bfloat16>::RangeReduction(bfloat16 x, double& modifier) {
    // Extract exponent and mantissa (where mantissa is between [1, 8))
    int m;
    float fx = frexpf((float)x, &m);
    fx *= 2.0;
    m--;
    
    int leftOver = m % 3;
    if (leftOver != 0) {
        if (leftOver < 0) leftOver += 3;
        
        floatint fi;
        fi.x = 0x3f800000 + (leftOver << 23);
        fx *= fi.f;
        m -= leftOver;
    }

    
    modifier = (m / 3);
    
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
bool Elementary<T>::FlipLbAndUb(double y) {
    return false;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 2000);
    PolySynth<bfloat16, Elementary<bfloat16>> bf16cbrt;
    bf16cbrt.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    bf16cbrt.FindPolynomials({0, 1, 2, 3, 4, 5, 6});
    printf("\tCOMPLETED\n\n");
    
    bf16cbrt.poly->PrintPiecewiseInfo();
    bf16cbrt.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}
