#include "Matrix.h"
#include <iostream>


Matrix::Matrix(int n) {
    this->n = n;
    this->matrix = new double* [n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n];
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0.0;
        }
    }
}

void Matrix::initializeMatrixA(double a1, double a2, double a3) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                matrix[i][j] = a1;
            }
            else if (i == j + 1 || i == j - 1) {
                matrix[i][j] = a2;
            }
            else if (i == j + 2 || i == j - 2) {
                matrix[i][j] = a3;
            }
        }
    }
}

Matrix* Matrix::multiply(Matrix* matrix) {
    if (this->n != matrix->n) {
        return NULL;
    }
    Matrix* result = new Matrix(matrix->n);
    double temp = 0;
    for (int i = 0; i < result->n; i++) {
        for (int j = 0; j < result->n; j++) {
            temp = 0;
            for (int k = 0; k < result->n; k++) {
                result->matrix[i][j] += this->matrix[i][k] * matrix->matrix[k][j];
            }
        }
    }
    return result;
}

double* Matrix::multiplyByVec(double* vector) {
    double* result = new double[n];
    double tmp = 0;

    for (int i = 0; i < n; i++) {
        tmp = 0;
        for (int j = 0; j < n; j++) {
            tmp += matrix[i][j] * vector[j];
        }
        result[i] = tmp;
    }
    return result;
}

void Matrix::copy(Matrix* matrix) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            this->matrix[i][j] = matrix->matrix[i][j];
        }
    }
}

int Matrix::getN() {
    return n;
}

Matrix::~Matrix() {
    for (int i = 0; i < n; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}