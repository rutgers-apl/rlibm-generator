#include "PolySynth.h"
#include "Elementary.h"

template <>
posit16 Elementary<posit16>::MpfrCalculateFunction(posit16 x) {
    // Code that computes the correctly rounded result for into x
}

template <>
bool Elementary<posit16>::ComputeSpecialCase(posit16 x, posit16& res) {
    // x is a special value:
    // 1. Return true.
    // 2. store the correct result in res.
    // x is not a special value:
    // return false.
}

template <>
double Elementary<posit16>::RangeReduction(posit16 x, double& modifier) {
    // Range reduction function.
    // return value = reduced input
    // modifier = value that is used for output compensation.
}

template <class T>
double Elementary<T>::RangePropagation(double y, double modifier) {
    // Output compensation function.
}

template <class T>
double Elementary<T>::ReverseRangePropagation(double y, double modifier) {
    // Inverse of output compensation function.
}

template <class T>
bool Elementary<T>::FlipLbAndUb(double modifier) {
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
    
    PolySynth<posit16, Elementary<posit16>> p16log;
    p16log.CalcIntervals();
    
    printf("FINDING POLYNOMIALS\n");
    // The degree of polynomial in array form
    p16log.FindPolynomials(/* EX {1, 3, 5, 7} */);
    printf("\tCOMPLETED\n\n");
    
    p16log.poly->PrintPiecewiseInfo();
    p16log.PerformErrorAnalysis();
    mpfr_clear(mval);
    return 0;
}

