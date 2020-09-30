#include "PolySynth.h"
#include "Elementary.h"
    
template <>
bfloat16 Elementary<bfloat16>::MpfrCalculateFunction(bfloat16 x) {
    mpfr_set_d(mval, (float)x, MPFR_RNDN);
    mpfr_log(mval, mval, MPFR_RNDN);
    return Elementary<bfloat16>::FromMPFR(mval);
}

template <>
bool Elementary<bfloat16>::ComputeSpecialCase(bfloat16 x, bfloat16& res) {
    if (x.val == 0x0 || x.val == 0x8000) {
        // If x == 0, then it should be -inf
        res.val = 0xFF80;
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
    // Extract exponent and mantissa (where mantissa is between [1, 2))
    int m;
    float fx = frexpf((float)x, &m);
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
    return (y + modifier) / 1.442695040888963387004650940070860087871551513671875;
}

template <class T>
double Elementary<T>::ReverseRangePropagation(double y, double modifier) {
    return y * 1.442695040888963387004650940070860087871551513671875 - modifier;
}

template <class T>
bool Elementary<T>::FlipLbAndUb(double y) {
    return false;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, 2000);
    
    PolySynth<bfloat16, Elementary<bfloat16>> bf16log;
    bf16log.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    bf16log.FindPolynomials({1, 3, 5, 7});
    printf("\tCOMPLETED\n\n");
    
    bf16log.poly->PrintPiecewiseInfo();
    bf16log.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}

