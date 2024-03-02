#include "Matrix.h"
#include "Solution.h"
#include <iostream>
#include <chrono>
#include <fstream>

#define N 941

using namespace std;


double* initializeVectorB(int n) {
    double* b = new double[n];
    for (int i = 0; i < n; i++) {
        b[i] = sin(i * 9.0);
    }
    return b;
}

void calculateResidual(Matrix* A, double* b, Solution* solution) {
    solution->residual = A->multiplyByVec(solution->x);
    for (int i = 0; i < A->getN(); i++) {
        solution->residual[i] -= b[i];
    }
}

double calculateNorm(double* res, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        norm += (res[i] * res[i]);
    }
    return sqrt(norm);
}

Solution* solveJacobi(Matrix* A, double* b, double maxNorm) {
    Solution* solution = new Solution(A->getN());
    double* xPrev = new double[A->getN()];
    double tmp = 0.0;
    int iterations = 0;

    for (int i = 0; i < A->getN(); i++) {
        xPrev[i] = 1;
    }

    auto start = chrono::system_clock::now();
    while (true) {
        for (int i = 0; i < A->getN(); i++) {
            tmp = 0;
            for (int j = 0; j < A->getN(); j++) {
                if (i != j) {
                    tmp += A->matrix[i][j] * xPrev[j];
                }
            }
            solution->x[i] = (b[i] - tmp) / A->matrix[i][i];
        }
        for (int i = 0; i < A->getN(); i++) {
            xPrev[i] = solution->x[i];
        }
        iterations++;
        calculateResidual(A, b, solution);
        solution->setNorm(calculateNorm(solution->residual, A->getN()));
        solution->normValues.push_back(solution->getNorm());
        if (solution->getNorm() <= maxNorm || isinf(solution->getNorm())) {
            solution->setIterations(iterations);
            break;
        }
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsedTime = end - start;
    solution->setTime(elapsedTime.count());
    delete[] xPrev;
    return solution;
}

Solution* solveGaussSeidel(Matrix* A,  double* b, double maxNorm) {
    Solution* solution = new Solution(A->getN());
    double tmp = 0.0;
    int iterations = 0;
    
    auto start = chrono::system_clock::now();
    while (true) {
        for (int i = 0; i < A->getN(); i++) {
            tmp = 0;
            for (int j = 0; j < A->getN(); j++) {
                if (i != j) {
                    tmp += A->matrix[i][j] * solution->x[j];
                }
            }
            solution->x[i] = (b[i] - tmp) / A->matrix[i][i];
        }
        iterations++;
        calculateResidual(A, b, solution);
        solution->setNorm(calculateNorm(solution->residual, A->getN()));
        solution->normValues.push_back(solution->getNorm());
        if (solution->getNorm() <= maxNorm || isinf(solution->getNorm())) {
            solution->setIterations(iterations);
            break;
        }
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsedTime = end - start;
    solution->setTime(elapsedTime.count());
    return solution;
}

Solution* LUfactorization(Matrix* A, double* b) {
    Solution* solution = new Solution(A->getN());
    Matrix* L, * U;
    double* y;
    double tmp = 0.0;

    y = new double[A->getN()];
    L = new Matrix(A->getN());
    L->initializeMatrixA(1.0, 0.0, 0.0);
    U = new Matrix(A->getN());
    U->copy(A);

    for (int i = 0; i < A->getN(); i++) {
        y[i] = 1;
    }

    auto start = chrono::system_clock::now();
    for (int i = 0; i < A->getN() - 1; i++) {
        for (int j = i + 1; j < A->getN(); j++) {
            L->matrix[j][i] = U->matrix[j][i] / U->matrix[i][i];
            for (int k = i; k < A->getN(); k++) {
                U->matrix[j][k] = U->matrix[j][k] - (L->matrix[j][i] * U->matrix[i][k]);
            }
        }
    }

    // forward substitution for Ly = b
    for (int i = 0; i < A->getN(); i++) {
        tmp = 0.0;
        for (int j = 0; j < i; j++) {
            tmp += L->matrix[i][j] * y[j];
        }
        y[i] = (b[i] - tmp) / L->matrix[i][i];
    }

    //backward substitution for Ux = y
    for (int i = A->getN() - 1; i >= 0; i--) {
        tmp = 0.0;
        for (int j = i + 1; j < A->getN(); j++) {
            tmp += U->matrix[i][j] * solution->x[j];
        }
        solution->x[i] = (y[i] - tmp) / U->matrix[i][i];
    }

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsedTime = end - start;
    solution->setTime(elapsedTime.count());
    calculateResidual(A, b, solution);
    solution->setNorm(calculateNorm(solution->residual, A->getN()));
    delete L;
    delete U;
    delete[] y;
    return solution;
}

int main()
{
    Matrix* A;
    double* b;
    Solution* solutionJacobi, *solutionGS, *solutionLU;
    list<double>::iterator it;
    int n;

    ofstream file("data.txt");
    if (!file.is_open()) {
        return -1;
    }

    // A
    A = new Matrix(N);
    A->initializeMatrixA(13.0, -1.0, -1.0);
    b = initializeVectorB(N);

    // B
    cout << " -- Zad. B -- \n N = 941  a1 = 13  a2 = a3 = -1\n";
    solutionJacobi = solveJacobi(A, b, 1e-9);
    cout << " Metoda Jacobiego:\n";
    solutionJacobi->printSolution();

    file << solutionJacobi->getIterations() << "\n";
    for (it = solutionJacobi->normValues.begin(); it != solutionJacobi->normValues.end(); it++) {
        file << *it << "\n";
    }
    file << "\n\n";


    solutionGS = solveGaussSeidel(A, b, 1e-9);
    cout << " Metoda Gaussa-Seidla:\n";
    solutionGS->printSolution();

    file << solutionGS->getIterations() << "\n";
    for (it = solutionGS->normValues.begin(); it != solutionGS->normValues.end(); it++) {
        file << *it << "\n";
    }
    file << "\n\n";
    delete solutionJacobi;
    delete solutionGS;


    // C
    A->initializeMatrixA(3.0, -1.0, -1.0);

    cout << " -- Zad. C -- \n N = 941  a1 = 3  a2 = a3 = -1\n";
    solutionJacobi = solveJacobi(A, b, 1e-9);
    cout << " Metoda Jacobiego:\n";
    solutionJacobi->printSolution();

    file << solutionJacobi->getIterations() << "\n";
    for (it = solutionJacobi->normValues.begin(); it != solutionJacobi->normValues.end(); it++) {
        file << *it << "\n";
    }
    file << "\n\n";


    solutionGS = solveGaussSeidel(A, b, 1e-9);
    cout << " Metoda Gaussa-Seidla:\n";
    solutionGS->printSolution();

    file << solutionGS->getIterations() << "\n";
    for (it = solutionGS->normValues.begin(); it != solutionGS->normValues.end(); it++) {
        file << *it << "\n";
    }
    file << "\n\n";


    // D
    cout << " -- Zad. D --\n";
    cout << " Metoda faktoryzacji LU:\n";
    solutionLU = LUfactorization(A, b);
    solutionLU->printLUsolution();


    // E
    cout << " -- Zad. E -- \n";

    for (int i = 0; i < 6; i++) {
        if (i == 0) {
            n = 100;
        }
        else if (i == 1) {
            n = 500;
        }
        else {
            n = (i - 1) * 1000;
        }
        
        delete A;
        delete[] b;
        delete solutionJacobi;
        delete solutionGS;
        delete solutionLU;

        A = new Matrix(n);
        A->initializeMatrixA(13, -1, -1);
        b = initializeVectorB(n);
        solutionJacobi = solveJacobi(A, b, 1e-9);
        solutionGS = solveGaussSeidel(A, b, 1e-9);
        solutionLU = LUfactorization(A, b);

        cout << " N = " << n << "\n";
        cout << " Metoda Jacobiego:\n";
        solutionJacobi->printSolution();
        cout << " Metoda Gaussa-Seidla:\n";
        solutionGS->printSolution();
        cout << " Metoda faktoryzacji LU:\n";
        solutionLU->printLUsolution();

        file << n << "\n";
        file << solutionJacobi->getTime() << "\n";
        file << solutionGS->getTime() << "\n";
        file << solutionLU->getTime() << "\n\n";
    }
    
    delete A;
    delete[] b;
    delete solutionJacobi;
    delete solutionGS;
    delete solutionLU;

    file.close();
    return 0;
}

