#ifndef GRADIENT_H
#define GRADIENT_H

#include "matrix.h"

class ConjugateGradientMethod{
private:
    Matrix returns;
    int n;

public:
    ConjugateGradientMethod(const Matrix &returns);

    Matrix build_Q() const;

    Matrix build_b(double target_return);

    Matrix build_x0();
    
    Matrix run(double target_return, double tolerance);

    Matrix train_test(int train_window, int test_window, double target_return, double tolerance);

};

#endif