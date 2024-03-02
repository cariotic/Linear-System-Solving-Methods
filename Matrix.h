#pragma once
#include <math.h>

class Matrix {
    int n;

public:
    double** matrix;

    Matrix(int n);
    void initializeMatrixA(double a1, double a2, double a3);
    Matrix* multiply(Matrix* matrix);
    double* multiplyByVec(double* vector);
    void copy(Matrix* matrix);
    int getN();
    ~Matrix();
};