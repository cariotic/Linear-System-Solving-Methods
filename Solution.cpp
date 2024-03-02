#include "Solution.h"

Solution::Solution(int n) {
	x = new double[n];
	residual = new double[n];
	for (int i = 0; i < n; i++) {
		x[i] = 1;
		residual[i] = 0;
	}
	norm = 0.0;
	time = 0.0;
	iterations = 0;
}

void Solution::printSolution() {
	std::cout << " Czas trwania: " << time << "s\n";
	std::cout << " Liczba iteracji: " << iterations << "\n";
	std::cout << " Norma z residuum: " << norm << "\n\n";
}

void Solution::printLUsolution() {
	std::cout << " Czas trwania: " << time << "s\n";
	std::cout << " Norma z residuum: " << norm << "\n\n";
}

int Solution::getIterations() {
	return iterations;
}

void Solution::incrementIterations() {
	iterations++;
}

void Solution::setIterations(int iterations) {
	this->iterations = iterations;
}

double Solution::getTime() {
	return time;
}

void Solution::setTime(double time) {
	this->time = time;
}

double Solution::getNorm() {
	return norm;
}

void Solution::setNorm(double norm) {
	this->norm = norm;
}

Solution::~Solution() {
	delete[] residual;
	delete[] x;
}