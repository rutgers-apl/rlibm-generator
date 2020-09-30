#include "PolySynth.h"
#include "Elementary.h"

template <>
bfloat16 Elementary<bfloat16>::MpfrCalculateFunction(bfloat16 x) {
    mpfr_const_pi(mval, MPFR_RNDN);
    mpfr_mul_d(mval, mval, (double)x, MPFR_RNDN);
    mpfr_cos(mval, mval, MPFR_RNDN);
    return Elementary<bfloat16>::FromMPFR(mval);
}

template <>
bool Elementary<bfloat16>::ComputeSpecialCase(bfloat16 x, bfloat16& res) {
    if ((x.val & 0x7FFF) > 0x7F80) {
        // Take care of NaN
        res.val = x.val;
        return true;
    } else if ((x.val & 0x7FFF) == 0x7F80) {
        res.val = 0x7FFF;
        return true;
    } else if (x >= 256.0f || x <= -256.0f) {
        res = 1.0f;
        return true;
    }
    
    return false;
}

template <>
double Elementary<bfloat16>::RangeReduction(bfloat16 x, double& modifier) {
    modifier = 1;

    if (x < 0.0f) x = -(float)x;
    
    double intPart;
    double frac = modf((double)x, &intPart);
    int iIntPart = intPart;
    
    if (iIntPart % 2 == 1) modifier *= -1;
    
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
    PolySynth<bfloat16, Elementary<bfloat16>> bf16cospi;
    bf16cospi.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    bf16cospi.FindPolynomialOnce({0});
    bf16cospi.FindPolynomialsUntil({0, 2, 4, 6}, 16127);
    bf16cospi.FindPolynomials({0});
    printf("\tCOMPLETED\n\n");
    
    bf16cospi.poly->PrintPiecewiseInfo();
    bf16cospi.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}
