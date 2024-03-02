#pragma once
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <list>

class Solution {
	double norm;
	double time;
	int iterations;

public:
	double* x;
	double* residual;
	std::list<double> normValues;

	Solution(int n);
	void printSolution();
	void printLUsolution();
	int getIterations();
	void incrementIterations();
	void setIterations(int iterations);
	double getTime();
	void setTime(double time);
	double getNorm();
	void setNorm(double norm);
	~Solution();
};