// 200710SAPLA.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include "stdafx.h"
#include "pch.h"
#include "SHARE_TOOL.h"
#include "GEOMETRY_TOOL.h"
#include "CAPCA.h"
#include "CAPCA_KNN.h"
#include "CPLA.h"
#include "CCHEBYSHEV.h"
#include "CAPLA.h"
#include "APLA_ICDE07.h"
#include "MULTI_DIMENSION.h"
#include "EVALUATION.h"
#include "./lib/saxquantizer.hpp"

typedef double ValueType;


int main()
{
	int option_tree = 0;
	cout << "Input: 0 Rtree. 1 updated Rtree" << endl;
	cin >> option_tree;
	Evaluation evaluation;
	evaluation.evaluate_multi_KNN_speed(option_tree);
	system("pause");
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
