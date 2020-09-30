#include "PolySynth.h"
#include "Elementary.h"

template <>
float Elementary<float>::MpfrCalculateFunction(float x) {
    // Code that computes the correctly rounded result for into x
}

template <>
bool Elementary<float>::ComputeSpecialCase(float x, float& res) {
    // x is a special value:
    // 1. Return true.
    // 2. store the correct result in res.
    // x is not a special value:
    // return false.
}

template <>
double Elementary<float>::RangeReduction(float x, double& modifier) {
    // Range reduction function.
    // return value = reduced input
    // modifier = value that is used for output compensation.
}
    
template <class T>
double Elementary<T>::RangePropagation(double yp, double modifier) {
    // Output compensation function.
}

template <class T>
double Elementary<T>::ReverseRangePropagation(double y, double modifier) {
    // Inverse of output compensation function.
}

template <class T>
bool Elementary<T>::FlipLbAndUb(double y) {
    // When computing interval, should we flip lb and ub?
    // This happens when range reduction function is a decreasing function.
    // So in essence:
    // 1. Return false if range reduction function is increasing function.
    // 2. Return true if range reduction function is decreasing function.
}

int main(int argc, char** argv) {
    // Initialize mval value with an appropriate precision. mval is used in
    // MpfrCalculateFunction
    mpfr_init2(mval, 2000);
    
    PolySynth<float, Elementary<float>> synth(SAMPLE, 5000);
    synth.SetUpSamplingOption(0.75, 1.5);
    synth.CalcIntervals();

    do {
        printf("FINDING POLYNOMIALS\n");
        fp32log2.ResetPolynomials();
        // The degree of polynomial in array form
        synth.FindPolynomials(/* EX {1, 3, 5, 7} */);
        printf("\tCOMPLETED\n\n");
    
        synth.poly->PrintPiecewiseInfo();
    } while (synth.TestAndAddSamplePoints());
    
    synth.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}
