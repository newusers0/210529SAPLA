#pragma once
#ifndef APLA_ICDE07_H
#define APLA_ICDE07_H
//#define DEBUG
#define TIME_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "GEOMETRY_TOOL.h"
#include "CPLA.h"
#include "CAPLA.h"

////#include <algorithm>
//#include <iterator>
//#include <queue>
////#include <string>
//#include <cassert>
//#include <iostream>
//#include <iomanip>
//#include <functional>
//#include <numeric>
//
//#include <time.h>
//#include <chrono>
//#include <ctime>
//#include <Windows.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <list>
////#include <vector>
//#include <fstream>
//#include <sstream>
//#include <cstdint>
//#include <cmath>
//#include <math.h>
//#include <complex>
//#include <random>
//#include <type_traits>
//#include <iosfwd>
//#include <limits>

using namespace std;


#define APLA_ICDE APLA_ICDE07<DataType>

TEMPLATE
class APLA_ICDE07 :virtual public GEOMETRY, virtual public TOOL, virtual public PLA_QUAL, virtual public APLA {
	//friend class CAPLA<DataType>;
public:
	struct TOOL::INPUT_ARGUMENT input_argument;
	struct TOOL::OUTPUT_ARGUMENT output_argument;//181214

	APLA_ICDE07() {};

	template<typename T>
	APLA_ICDE07(const T& const n, const T& const N);
	//191016
	//template<typename T>
	//double getSegmentMaxDifference(const vector<DataType>& const original_time_series_vector, const T& const segment_pla);
	//191016
	void getAPLA_ICDE07(typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<DataType>& const original_time_series_vector);// From paper APLA ICDE 07

	template<typename T> //PLA_QUAL::PLA_C 191121 has return value
	void getAPLA_ICDE07(typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<DataType>& const original_time_series_vector, DoublyLinkedList<T>& const F_result_vector);// From paper APLA ICDE 07 
};

#include "APLA_ICDE07.cpp"

#endif