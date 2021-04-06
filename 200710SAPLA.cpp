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

int main()
{
	Evaluation evaluation;
	evaluation.evaluate_multi_KNN_speed();
	system("pause");
	return 0;
}

