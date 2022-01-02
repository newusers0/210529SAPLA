// 200710SAPLA.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include <windows.h>
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
#include "CDFT.h"
#include "EVALUATION.h"
//#include "./lib/dtwrecoge.h"
#include "./lib/saxquantizer.hpp"

typedef double ValueType;

struct Rect
{
	Rect() {}

	Rect(DataType a_minX, DataType a_minY, DataType a_maxX, DataType a_maxY)
	{
		min[0] = a_minX;
		min[1] = a_minY;

		max[0] = a_maxX;
		max[1] = a_maxY;
	}


	DataType min[2];
	DataType max[2];
};

struct Rect rects[] =
{
  Rect(0, 0, 2, 2), // xmin, ymin, xmax, ymax (for 2 dimensional RTree)
  Rect(5, 5, 7, 7),
  Rect(8, 5, 9, 6),
  Rect(7, 1, 9, 2),
};

int nrects = sizeof(rects) / sizeof(rects[0]);

Rect search_rect(6, 4, 10, 6); // search will find above rects that this one overlaps


bool MySearchCallback(ValueType id)
{
	cout << "Hit data rect " << id << "\n";
	return true; // keep going
}

int main()
{
	std::cout << "Hello World!\n";

	SYSTEM_INFO si;
	GetSystemInfo(&si);

	printf("The page size for this system is %u bytes.\n", si.dwPageSize);// page size is 4096 bytes
	cout << "Size of DataType: " << sizeof(DataType) << " byte" << endl;//DataType 8 byte
	cout << "Size of char: " << sizeof(char) << " byte" << endl;//1 byte
	cout << "Size of int: " << sizeof(int) << " bytes" << endl;//4 byte
	cout << "Size of float: " << sizeof(float) << " bytes" << endl;//4 byte
	cout << "Size of double: " << sizeof(double) << " bytes" << endl;//8 byte
	cout << "Size of long double: " << sizeof(long double) << " bytes" << endl;//8 byte


	/*............................................................................*/
#ifdef _DEBUG
	cout << "Test Release.\n";
	//assert(0);
#endif
	/*............................................................................*/

	/////////////////////////////////////////////////////

	
		vector<string> str_suffix_vector = { "1", "2", "3" };
		int tree_type_0 = 0;
		int tree_type_1 = 2;
		int option_homogenous_data_type = 7;
		int size_file = 117;// 45 45 27
		int file_id_begin = 0;// 0 45 90
		int size_query_time_series = 5;// Defualt 5
		int n = 1024;
		int initial_N = 12;//default 12
		int final_N = 24;
		int initial_K = 4;
		int final_K = 64;
		int point_number = 100;//default number: 50
		int max_node = 5;
		vector<int> representation_option_vector = { int(INF)};
		

		typename TOOL::template EVALUATION_ARGUMENT<int> evaluation_argument_struct_0(str_suffix_vector[0], tree_type_0, option_homogenous_data_type, size_file,
			file_id_begin, size_query_time_series, n, initial_N, final_N, initial_K, final_K, point_number, max_node, representation_option_vector);

		typename TOOL::template EVALUATION_ARGUMENT<int> evaluation_argument_struct_1(str_suffix_vector[1], tree_type_1, option_homogenous_data_type, size_file,
			file_id_begin, size_query_time_series, n, initial_N, final_N, initial_K, final_K, point_number, max_node, representation_option_vector);

		/*##############################################################################################################################*/

		Evaluation evaluation;

		evaluation.evaluate_multi_KNN_speed(evaluation_argument_struct_0);

		evaluation.evaluate_multi_KNN_speed(evaluation_argument_struct_1);

		
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
