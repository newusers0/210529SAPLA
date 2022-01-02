#pragma once
#ifndef _CPLA_CPP_
#define _CPLA_CPP_

#include "CPLA.h"

////TEMPLATE
////struct PLA_QUAL::TIME {
////	_LARGE_INTEGER time_start;
////	_LARGE_INTEGER time_over;   //finish time
////	double dqFrequency = NULL;      //timer frequency
////	double run_time = NULL;
////
////	~TIME() {
////		time_start.QuadPart = NULL;
////		time_over.QuadPart = NULL;   //finish time
////		dqFrequency = NULL;      //timer frequency
////		run_time = NULL;
////	}
////};
//
////TEMPLATE
////struct PLA_QUAL::INPUT_ARGUMENT {
////	int time_series_length = INF;//n
////	int point_dimension = INF;//N
////	int point_number = INF;
////	int remainder = INF;
////	int segment_length_first = INF;//l=n/m+1
////	int segment_length_second = INF;//l=n/m
////
////	int rtree_max_nodes = INF;
////	int K = INF;
////	string read_file_name;
////	string* write_file_name = nullptr;
////
////	//for Chebyshev
////	int degree_m = INF;//point_dimension for APCA & PLA
////	int arity_d = INF;// dimension of every point, for multi-dimensional case
////
////	//for KNN
////	double pruning_power = INF;
////	double sum_distance_euc = INF;
////
////	//for time
////	double build_rtree_time = INF;
////	double approximation_query_time = INF;
////	double knn_rest_part_time = INF;
////	double knn_total_time = INF;
////
////	//    three part time
////	double navigate_index_time = INF;// navigate time
////	double distance_lowbound_time = INF; // distance chebyshev, PLA, APCA time
////	double distance_euc_time = INF;// distance euclidean time
////
////	double result_accuracy = INF;// KNN result accuracy
////
////	//I/O cost
////	double IO_cost = NULL;
////
////	~INPUT_ARGUMENT() {
////		time_series_length = INF;//n
////		point_dimension = INF;//N
////		point_number = INF;
////		remainder = INF;
////		segment_length_first = INF;//l=n/m+1
////		segment_length_second = INF;//l=n/m
////
////		rtree_max_nodes = INF;
////		K = INF;
////
////		pruning_power = INF;
////		sum_distance_euc = INF;
////
////		build_rtree_time = INF;
////		approximation_query_time = INF;
////		knn_rest_part_time = INF;
////		knn_total_time = INF;
////
////		//    three part time
////		navigate_index_time = NULL;// navigate time
////		distance_lowbound_time = NULL; // distance chebyshev, PLA, APCA time
////		distance_euc_time = NULL;// distance euclidean time
////
////		IO_cost = NULL;
////
////		result_accuracy = NULL;// KNN result accuracy
////
////		read_file_name.clear();
////		read_file_name.shrink_to_fit();
////
////		/*write_file_name.clear();
////		write_file_name.shrink_to_fit();*/
////		write_file_name = nullptr;
////		//read_file_name = nullptr;
////
////		//for Chebyshev
////		degree_m = NULL;
////	}
////};
//
////TEMPLATE
////struct PLA_QUAL::PLA {
////	DataType* a = nullptr;
////	DataType* b = nullptr;
////	int segmentNum = NULL; //the dimension of the index(eg. MBR, imput parameter)
////
////	/*================================================*/
////	double right_endpoint = INF;//190619
////	int segment_width = INF; //190619 the number of points in every segment
////	double sum_value = INF; //191021
////	double segment_area = INF;//190619 area of
////	double segment_density = INF;//190619
////	double segment_area0 = INF;//190619 area of
////	double segment_density0 = INF;//190619
////	double parallelogram_height = INF;//190619
////	/*...............................................*/
////
////	~PLA() {
////		segmentNum = NULL;
////
////		right_endpoint = INF;
////		segment_width = INF;
////		sum_value = INF; //191021
////		segment_area = INF;
////		segment_density = INF;
////		parallelogram_height = INF;
////		/*if (a != nullptr) {
////			delete a;
////			a = nullptr;
////		}
////		if (b != nullptr) {
////			delete b;
////			b = nullptr;
////		}*/
////	}
////};
//
////**********************************************************************************************************************************
//// Method:PLA_C
//// Qualifier:No pointer, use vairable a&b to instead pointer 
//// Input:
//// Output:
//// date:191022
//// author:
////**********************************************************************************************************************************
//TEMPLATE
//struct PLA_QUAL::PLA_C {
//	TOOL::APLA_COEFFICIENT apla;
//	/*DataType a = INF;
//	DataType b = INF;*/
//	/*================================================*/
//	double right_endpoint = INF;//190619
//	double rectangle_width = INF; //190619 the number of points in every segment
//	double sum_value = INF; //191021
//	/*...............................................*/
//	~PLA_C() {
//		/*a = INF;
//		b = INF;*/
//		right_endpoint = INF;
//		rectangle_width = INF;
//		sum_value = INF; //191021
//	}
//};
//
//TEMPLATE
//struct PLA_QUAL::POINT {
//	DataType u_A_1 = NULL;
//	DataType u_A_2 = NULL;
//	DataType v_A_1 = NULL;
//	DataType v_A_2 = NULL;
//
//	DataType u_B_1 = NULL;
//	DataType u_B_2 = NULL;
//	DataType v_B_1 = 0;
//	DataType v_B_2 = 0;
//
//	DataType u_B_C_1 = NULL; //perpendicular x
//	DataType v_B_C_1 = NULL; //perpendicular x
//	DataType u_B_C_2 = NULL; //perpendicular x
//	DataType v_B_C_2 = NULL; //perpendicular y
//};
//
//TEMPLATE
//struct PLA_QUAL::PLA_NODE_PAIR {								//queue data structure.
//	double d_dist = NULL;										//dist.
//	int original_time_series_id = NULL;                  // APCA point ID.
//	PLA* p_PLA_point = nullptr;						//data point. if NULL, this is internal node.
//	RTREE::Node* p_rtree_node = nullptr;            //subNode, if NULL, this is apca point.
//};
//
//TEMPLATE
//struct PLA_QUAL::priorityIncrement {//small to big
//	virtual bool operator ()(const PLA_NODE_PAIR& a, const PLA_NODE_PAIR& b) {
//		return a.d_dist > b.d_dist;
//	}
//};
//
//TEMPLATE
//struct PLA_QUAL::ORIGINAL_TIME_SERIES_PAIR {  //for temp queue
//	double d_dist = NULL;
//	//double *p_original_time_series = nullptr;
//	int original_time_series_id = NULL;
//};
//
//TEMPLATE
//struct PLA_QUAL::priorityDistanceEUC {//big to small
//	bool operator ()(const ORIGINAL_TIME_SERIES_PAIR& a, const ORIGINAL_TIME_SERIES_PAIR& b) {
//		return a.d_dist > b.d_dist;
//	}
//};
//
//TEMPLATE
//bool PLA_QUAL::compare(ORIGINAL_TIME_SERIES_PAIR& first, ORIGINAL_TIME_SERIES_PAIR& second)
//{
//	return first.d_dist < second.d_dist;
//}
//
////TEMPLATE
////PLA_QUAL::CPLA(const int& n, const int& N) {
////	//For PLA input_argument
////	input_argument.time_series_length = n;
////	input_argument.point_dimension = N;
////
////	input_argument.remainder = int(n) % int(N);
////	double integerDividend = n - input_argument.remainder;
////	input_argument.segment_length_second = integerDividend / N;
////	input_argument.segment_length_first = input_argument.segment_length_second + 1;
////
////	assert(input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
////
////	//For TOOL input_argument
////	tool_input_argument.time_series_length = n;
////	tool_input_argument.point_dimension = N;
////
////	tool_input_argument.remainder = int(n) % int(N);
////	double tool_integerDividend = n - tool_input_argument.remainder;
////	tool_input_argument.segment_length_second = tool_integerDividend / N;
////	tool_input_argument.segment_length_first = tool_input_argument.segment_length_second + 1;
////
////	assert(tool_input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
////}
//
//TEMPLATE
//PLA_QUAL::CPLA(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const string& read_file_name, string*& const write_file_name) {
//	assert(K <= point_number);
//
//	//For PLA Argument
//	//****************PLA*********************************************************
//	input_argument.remainder = int(n) % int(N);//For PLA
//	double integerDividend = n - input_argument.remainder;
//	input_argument.segment_length_second = integerDividend / N;
//	input_argument.segment_length_first = input_argument.segment_length_second + 1;
//	assert(input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
//	//*******************************************************************************
//
//	input_argument.time_series_length = n;
//	input_argument.point_dimension = N;
//	input_argument.point_number = point_number;
//	input_argument.rtree_max_nodes = rtree_max_nodes;
//	input_argument.K = K;
//	input_argument.read_file_name = read_file_name;
//	input_argument.write_file_name = write_file_name;
//
//	input_argument.pruning_power = 0.0;
//	input_argument.sum_distance_euc = 0.0;
//	//time
//	input_argument.build_rtree_time = 0.0;
//	input_argument.approximation_query_time = 0.0;
//	input_argument.knn_rest_part_time = 0.0;
//	input_argument.knn_total_time = 0.0;
//
//	input_argument.navigate_index_time = 0;// navigate time
//	input_argument.distance_lowbound_time = 0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0;// distance euclidean time
//
//	input_argument.IO_cost = 0.0;
//
//	input_argument.result_accuracy = 0;// KNN result accuracy
//
//	/*cout << "input_argument.time_series_length: "<< input_argument.time_series_length <<endl;
//	cout << "input_argument.point_dimension: " << input_argument.point_dimension << endl;
//	cout << "input_argument.remainder: " << input_argument.remainder << endl;
//	cout << "input_argument.segment_length_first: " << input_argument.segment_length_first << endl;
//	cout << "input_argument.segment_length_second: " << input_argument.segment_length_second << endl;*/
//
//	//For TOOL Argument
//
//	//****************PLA*********************************************************
//	tool_input_argument.remainder = int(n) % int(N);//For PLA
//	integerDividend = n - tool_input_argument.remainder;
//	tool_input_argument.segment_length_second = integerDividend / N;
//	tool_input_argument.segment_length_first = tool_input_argument.segment_length_second + 1;
//	assert(tool_input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
//	//*******************************************************************************
//
//	tool_input_argument.time_series_length = n;
//	tool_input_argument.point_dimension = N;
//	tool_input_argument.point_number = point_number;
//	tool_input_argument.rtree_max_nodes = rtree_max_nodes;
//	tool_input_argument.K = K;
//	tool_input_argument.read_file_name = read_file_name;
//	tool_input_argument.write_file_name = write_file_name;
//
//	tool_input_argument.pruning_power = 0.0;
//	tool_input_argument.sum_distance_euc = 0.0;
//	//time
//	tool_input_argument.build_rtree_time = 0.0;
//	tool_input_argument.approximation_query_time = 0.0;
//	tool_input_argument.knn_rest_part_time = 0.0;
//	tool_input_argument.knn_total_time = 0.0;
//
//	tool_input_argument.navigate_index_time = 0;// navigate time
//	tool_input_argument.distance_lowbound_time = 0; // distance chebyshev, PLA, APCA time
//	tool_input_argument.distance_euc_time = 0;// distance euclidean time
//
//	tool_input_argument.IO_cost = 0.0;
//
//	tool_input_argument.result_accuracy = 0;// KNN result accuracy
//}
//
//TEMPLATE
//void PLA_QUAL::initialPLA(PLA& pla, const int& N) {
//	pla.segmentNum = N;
//	pla.a = new DataType[N];
//	pla.b = new DataType[N];
//}
//
//TEMPLATE
//void PLA_QUAL::deletePLA(PLA& pla) {
//	pla.segmentNum = INF;
//
//	delete[] pla.a;
//	pla.a = nullptr;
//
//	delete[] pla.b;
//	pla.b = nullptr;
//}
//
//TEMPLATE
//void PLA_QUAL::initialPLAArray(const INPUT_ARGUMENT& input_argument, PLA*& const pla_array) {
//	for (int i = 0; i < input_argument.point_number; i++) {
//		pla_array[i].segmentNum = input_argument.point_dimension;
//		pla_array[i].a = new DataType[input_argument.point_dimension];
//		pla_array[i].b = new DataType[input_argument.point_dimension];
//		//initialPLA(pla_array[i], input_argument.point_dimension);
//		fill_n(pla_array[i].a, input_argument.point_dimension, INF);
//		fill_n(pla_array[i].b, input_argument.point_dimension, INF);
//	}
//}
//
//TEMPLATE
//void PLA_QUAL::initialPLAArray(const typename TOOL::INPUT_ARGUMENT& input_argument, PLA*& const pla_array) {
//	pla_array = new PLA[input_argument.point_number];
//
//	for (int i = 0; i < input_argument.point_number; i++) {
//		pla_array[i].segmentNum = input_argument.point_dimension;
//		pla_array[i].a = new DataType[input_argument.point_dimension];
//		pla_array[i].b = new DataType[input_argument.point_dimension];
//		//initialPLA(pla_array[i], input_argument.point_dimension);
//		fill_n(pla_array[i].a, input_argument.point_dimension, INF);
//		fill_n(pla_array[i].b, input_argument.point_dimension, INF);
//	}
//}
//
//TEMPLATE
//void PLA_QUAL::initialPLAArray(const int& point_number, const int& point_dimension, PLA*& const pla_array) {
//	pla_array = new PLA[point_number];
//
//	for (int i = 0; i < point_number; i++) {
//		pla_array[i].segmentNum = point_dimension;
//		pla_array[i].a = new DataType[point_dimension];
//		pla_array[i].b = new DataType[point_dimension];
//		//initialPLA(pla_array[i], input_argument.point_dimension);
//		fill_n(pla_array[i].a, point_dimension, INF);
//		fill_n(pla_array[i].b, point_dimension, INF);
//	}
//}
//
//TEMPLATE
//void PLA_QUAL::deletePLAArray(const INPUT_ARGUMENT& input_argument, PLA*& const pla_array) {
//	for (int i = 0; i < input_argument.point_number; i++) {
//		delete[] pla_array[i].a;
//		pla_array[i].a = nullptr;
//
//		delete[] pla_array[i].b;
//		pla_array[i].b = nullptr;
//	}
//	delete[] pla_array;
//	pla_array = nullptr;
//}
//
//TEMPLATE
//void PLA_QUAL::deletePLAArray(const typename TOOL::INPUT_ARGUMENT& input_argument, PLA*& const pla_array) {
//	for (int i = 0; i < input_argument.point_number; i++) {
//		delete[] pla_array[i].a;
//		pla_array[i].a = nullptr;
//
//		delete[] pla_array[i].b;
//		pla_array[i].b = nullptr;
//	}
//	delete[] pla_array;
//	pla_array = nullptr;
//}
//
//TEMPLATE
//void PLA_QUAL::getPLA(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, const PLA& pla) {
//	//printf("getPLA()\n");
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.time_series_length != INF && input_argument.point_dimension != INF && input_argument.segment_length_first != INF && input_argument.segment_length_second != INF);
//	/*cout << "input_argument.time_series_length: " << input_argument.time_series_length << endl;
//	cout << "input_argument.point_dimension: " << input_argument.point_dimension << endl;
//	cout << "input_argument.remainder: " << input_argument.remainder << endl;
//	cout << "input_argument.segment_length_first: " << input_argument.segment_length_first << endl;
//	cout << "input_argument.segment_length_second: " << input_argument.segment_length_second << endl;*/
//
//	//a
//	double a_first_divisor = input_argument.segment_length_second * input_argument.segment_length_first * (input_argument.segment_length_first + 1.0); //n(n-1)(n+1)
//	double a_second_divisor = input_argument.segment_length_second * (input_argument.segment_length_second - 1.0) * input_argument.segment_length_first;  //n(n-1)(n+1)
//	double a_first_minuend = input_argument.segment_length_second / 2.0;//(n-1)/2
//	double a_second_minuend = (input_argument.segment_length_second - 1.0) / 2.0;//(n-1)/2
//	  //b
//	double b_first_divisor = input_argument.segment_length_first * (1.0 + input_argument.segment_length_first); //n(1+n)
//	double b_second_divisor = input_argument.segment_length_second * input_argument.segment_length_first; //n(1+n)
//	double b_first_minuend = 2.0 * input_argument.segment_length_first - 1.0;//2n-1
//	double b_second_minuend = 2.0 * input_argument.segment_length_second - 1.0;//2n-1
//
//	int j = 0;
//	int endOfLongSegment = 0, indexOfLongSegment = 0;
//	/*int account[1280];
//	int* pointer = account;
//	int index = 0;*/
//
//	double a_sum = NULL;
//	double b_sum = NULL;
//	double t = NULL;
//	for (int i = 1; i <= input_argument.point_dimension; i++) {
//		//cout<<"i: "<<i << endl;
//		a_sum = 0;
//		b_sum = 0;
//		t = 0;
//		if (i <= input_argument.remainder) {//first part
//			for (j = int((i - 1) * input_argument.segment_length_first); j < i * input_argument.segment_length_first; j++) {
//				//cout << "    j: "<<j << endl;
//				a_sum += (t - a_first_minuend) * original_time_series[j];
//				b_sum += (b_first_minuend - t * 3.0) * original_time_series[j];
//				/*account[index] = j;
//				index++;*/
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//y[i - 1] = sum / input_argument.segment_length_first;
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			pla.a[i - 1] = 12.0 * a_sum / a_first_divisor;
//			pla.b[i - 1] = 2.0 * b_sum / b_first_divisor;
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//			indexOfLongSegment = i;
//			endOfLongSegment = int(i * input_argument.segment_length_first);
//		}
//		else {//second part
//			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * input_argument.segment_length_second); j < endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second; j++) {
//				assert(j >= 0);
//				//cout << "    j: " << j << endl;
//				a_sum += (t - a_second_minuend) * original_time_series[j];
//				b_sum += (b_second_minuend - t * 3.0) * original_time_series[j];
//				//account[index] = j;
//				////cout << account[index] << endl;
//				//if (index!=0) assert(account[index]-account[index-1]==1);
//				//index++;
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			pla.a[i - 1] = 12.0 * a_sum / a_second_divisor;
//			pla.b[i - 1] = 2.0 * b_sum / b_second_divisor;
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//		}
//	}
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(account[i] - account[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//	/*APCA_KNN_QUAL::printArray(pointer,1280);*/
//}
//
//TEMPLATE//190501
//void PLA_QUAL::getPLAVector(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, vector<PLA>& const pla) {
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.time_series_length != INF && input_argument.point_dimension != INF && input_argument.segment_length_first != INF && input_argument.segment_length_second != INF);
//	//printf("getPLA()\n");
//	//assert(pla.size() == input_argument.point_dimension);
//
//	/*cout << "input_argument.time_series_length: " << input_argument.time_series_length << endl;
//	cout << "input_argument.point_dimension: " << input_argument.point_dimension << endl;
//	cout << "input_argument.remainder: " << input_argument.remainder << endl;
//	cout << "input_argument.segment_length_first: " << input_argument.segment_length_first << endl;
//	cout << "input_argument.segment_length_second: " << input_argument.segment_length_second << endl;*/
//
//	PLA temp_pla;
//	//temp_pla.a = new DataType;
//	//temp_pla.b = new DataType;
//
//	//a
//	double a_first_divisor = input_argument.segment_length_second * input_argument.segment_length_first * (input_argument.segment_length_first + 1.0); //n(n-1)(n+1)
//	double a_second_divisor = input_argument.segment_length_second * (input_argument.segment_length_second - 1.0) * input_argument.segment_length_first;  //n(n-1)(n+1)
//	double a_first_minuend = input_argument.segment_length_second / 2.0;//(n-1)/2
//	double a_second_minuend = (input_argument.segment_length_second - 1.0) / 2.0;//(n-1)/2
//	  //b
//	double b_first_divisor = input_argument.segment_length_first * (1.0 + input_argument.segment_length_first); //n(1+n)
//	double b_second_divisor = input_argument.segment_length_second * input_argument.segment_length_first; //n(1+n)
//	double b_first_minuend = 2.0 * input_argument.segment_length_first - 1.0;//2n-1
//	double b_second_minuend = 2.0 * input_argument.segment_length_second - 1.0;//2n-1
//
//	int j = 0;
//	int endOfLongSegment = 0, indexOfLongSegment = 0;
//	/*int account[1280];
//	int* pointer = account;
//	int index = 0;*/
//
//	double a_sum = NULL;
//	double b_sum = NULL;
//	double t = NULL;
//	for (int i = 1; i <= input_argument.point_dimension; i++) {
//		//cout<<"i: "<<i << endl;
//		a_sum = 0;
//		b_sum = 0;
//		t = 0;
//
//		if (i <= input_argument.remainder) {//first part
//			for (j = int((i - 1) * input_argument.segment_length_first); j < i * input_argument.segment_length_first; j++) {
//				//cout << "    j: "<<j << endl;
//				a_sum += (t - a_first_minuend) * original_time_series[j];
//				b_sum += (b_first_minuend - t * 3.0) * original_time_series[j];
//				/*account[index] = j;
//				index++;*/
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//y[i - 1] = sum / input_argument.segment_length_first;
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			//pla.a[i - 1] = 12.0 * a_sum / a_first_divisor;
//			//pla.b[i - 1] = 2.0 * b_sum / b_first_divisor;
//			temp_pla.a = new DataType;
//			temp_pla.b = new DataType;
//			*temp_pla.a = 12.0 * a_sum / a_first_divisor;
//			*temp_pla.b = 2.0 * b_sum / b_first_divisor;
//			pla.push_back(temp_pla);
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//			indexOfLongSegment = i;
//			endOfLongSegment = int(i * input_argument.segment_length_first);
//		}
//		else {//second part
//			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * input_argument.segment_length_second); j < endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second; j++) {
//				assert(j >= 0);
//				//cout << "    j: " << j << endl;
//				a_sum += (t - a_second_minuend) * original_time_series[j];
//				b_sum += (b_second_minuend - t * 3.0) * original_time_series[j];
//				//account[index] = j;
//				////cout << account[index] << endl;
//				//if (index!=0) assert(account[index]-account[index-1]==1);
//				//index++;
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			//pla.a[i - 1] = 12.0 * a_sum / a_second_divisor;
//			//pla.b[i - 1] = 2.0 * b_sum / b_second_divisor;
//			temp_pla.a = new DataType;
//			temp_pla.b = new DataType;
//			*temp_pla.a = 12.0 * a_sum / a_second_divisor;
//			*temp_pla.b = 2.0 * b_sum / b_second_divisor;
//			pla.push_back(temp_pla);
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//		}
//	}
//
//	/*for (auto&&i : pla) {190501
//		cout << *i.a<<" "<< *i.b << endl;
//	}*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(account[i] - account[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//	/*APCA_KNN_QUAL::printArray(pointer,1280);*/
//}
//
//TEMPLATE
//void PLA_QUAL::getPLAVectorNormal(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, vector<PLA>& const pla) {//190604
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.time_series_length != INF && input_argument.point_dimension != INF && input_argument.segment_length_first != INF && input_argument.segment_length_second != INF);
//	//printf("getPLA()\n");
//	//assert(pla.size() == input_argument.point_dimension);
//
//	/*cout << "input_argument.time_series_length: " << input_argument.time_series_length << endl;
//	cout << "input_argument.point_dimension: " << input_argument.point_dimension << endl;
//	cout << "input_argument.remainder: " << input_argument.remainder << endl;
//	cout << "input_argument.segment_length_first: " << input_argument.segment_length_first << endl;
//	cout << "input_argument.segment_length_second: " << input_argument.segment_length_second << endl;*/
//
//	PLA temp_pla;
//	//temp_pla.a = new DataType;
//	//temp_pla.b = new DataType;
//
//	//a
//	double a_first_divisor = input_argument.segment_length_second * input_argument.segment_length_first * (input_argument.segment_length_first + 1.0); //n(n-1)(n+1)
//	double a_second_divisor = input_argument.segment_length_second * (input_argument.segment_length_second - 1.0) * input_argument.segment_length_first;  //n(n-1)(n+1)
//	double a_first_minuend = input_argument.segment_length_second / 2.0;//(n-1)/2
//	double a_second_minuend = (input_argument.segment_length_second - 1.0) / 2.0;//(n-1)/2
//	  //b
//	double b_first_divisor = input_argument.segment_length_first * (1.0 + input_argument.segment_length_first); //n(1+n)
//	double b_second_divisor = input_argument.segment_length_second * input_argument.segment_length_first; //n(1+n)
//	double b_first_minuend = 2.0 * input_argument.segment_length_first - 1.0;//2n-1
//	double b_second_minuend = 2.0 * input_argument.segment_length_second - 1.0;//2n-1
//
//	int j = 0;
//	int endOfLongSegment = 0, indexOfLongSegment = 0;
//	/*int account[1280];
//	int* pointer = account;
//	int index = 0;*/
//
//	double a_sum = NULL;
//	double b_sum = NULL;
//	double t = NULL;
//	for (int i = 1; i <= input_argument.point_dimension; i++) {
//		//cout<<"i: "<<i << endl;
//		a_sum = 0;
//		b_sum = 0;
//		t = 0;
//
//		if (i <= input_argument.remainder) {//first part
//			for (j = int((i - 1) * input_argument.segment_length_first); j < i * input_argument.segment_length_first; j++) {
//				//cout << "    j: "<<j << endl;
//				a_sum += (t - (input_argument.segment_length_second / 2.0)) * original_time_series[j];
//				b_sum += ((2.0 * input_argument.segment_length_first - 1.0) - t * 3.0) * original_time_series[j];
//				/*account[index] = j;
//				index++;*/
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//y[i - 1] = sum / input_argument.segment_length_first;
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			//pla.a[i - 1] = 12.0 * a_sum / a_first_divisor;
//			//pla.b[i - 1] = 2.0 * b_sum / b_first_divisor;
//			temp_pla.a = new DataType;
//			temp_pla.b = new DataType;
//			*temp_pla.a = 12.0 * a_sum / (input_argument.segment_length_second * input_argument.segment_length_first * (input_argument.segment_length_first + 1.0));
//			*temp_pla.b = 2.0 * b_sum / (input_argument.segment_length_first * (1.0 + input_argument.segment_length_first));
//			pla.push_back(temp_pla);
//			/*===================================================================*/
//			pla[i - 1].segment_width = input_argument.segment_length_first;
//			pla[i - 1].right_endpoint = i * input_argument.segment_length_first - 1;
//			/*...................................................................*/
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//			indexOfLongSegment = i;
//			endOfLongSegment = int(i * input_argument.segment_length_first);
//		}
//		else {//second part
//			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * input_argument.segment_length_second); j < endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second; j++) {
//				assert(j >= 0);
//				//cout << "    j: " << j << endl;
//				a_sum += (t - ((input_argument.segment_length_second - 1.0) / 2.0)) * original_time_series[j];
//				b_sum += ((2.0 * input_argument.segment_length_second - 1.0) - t * 3.0) * original_time_series[j];
//				//account[index] = j;
//				////cout << account[index] << endl;
//				//if (index!=0) assert(account[index]-account[index-1]==1);
//				//index++;
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			//pla.a[i - 1] = 12.0 * a_sum / a_second_divisor;
//			//pla.b[i - 1] = 2.0 * b_sum / b_second_divisor;
//			temp_pla.a = new DataType;
//			temp_pla.b = new DataType;
//			*temp_pla.a = 12.0 * a_sum / (input_argument.segment_length_second * (input_argument.segment_length_second - 1.0) * input_argument.segment_length_first);
//			*temp_pla.b = 2.0 * b_sum / (input_argument.segment_length_second * input_argument.segment_length_first);
//			pla.push_back(temp_pla);
//			/*===================================================================*/
//			pla[i - 1].segment_width = input_argument.segment_length_second;
//			pla[i - 1].right_endpoint = endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second - 1;
//			/*...................................................................*/
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//		}
//	}
//
//	/*for (auto&&i : pla) {190501
//		cout << *i.a<<" "<< *i.b << endl;
//	}*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(account[i] - account[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//	/*APCA_KNN_QUAL::printArray(pointer,1280);*/
//}
//
////**********************************************************************************************************************************
//// Method:getPLALinkedListNormal
//// Qualifier:Use Linked List ot instead Vector
//// Input:
//// Output:
//// date:191104
//// author:
////**********************************************************************************************************************************
//TEMPLATE
//template<typename T>
//void PLA_QUAL::getPLALinkedListNormal(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, DoublyLinkedList<T>& const pla) {//191104 11:12
//	assert(input_argument.time_series_length != INF && input_argument.point_dimension != INF && input_argument.segment_length_first != INF && input_argument.segment_length_second != INF);
//	//printf("getPLA()\n");
//	//assert(pla.size() == input_argument.point_dimension);
//
//	/*cout << "input_argument.time_series_length: " << input_argument.time_series_length << endl;
//	cout << "input_argument.point_dimension: " << input_argument.point_dimension << endl;
//	cout << "input_argument.remainder: " << input_argument.remainder << endl;
//	cout << "input_argument.segment_length_first: " << input_argument.segment_length_first << endl;
//	cout << "input_argument.segment_length_second: " << input_argument.segment_length_second << endl;*/
//
//	PLA temp_pla;
//	//temp_pla.a = new DataType;
//	//temp_pla.b = new DataType;
//
//	//a
//	double a_first_divisor = input_argument.segment_length_second * input_argument.segment_length_first * (input_argument.segment_length_first + 1.0); //n(n-1)(n+1)
//	double a_second_divisor = input_argument.segment_length_second * (input_argument.segment_length_second - 1.0) * input_argument.segment_length_first;  //n(n-1)(n+1)
//	double a_first_minuend = input_argument.segment_length_second / 2.0;//(n-1)/2
//	double a_second_minuend = (input_argument.segment_length_second - 1.0) / 2.0;//(n-1)/2
//	  //b
//	double b_first_divisor = input_argument.segment_length_first * (1.0 + input_argument.segment_length_first); //n(1+n)
//	double b_second_divisor = input_argument.segment_length_second * input_argument.segment_length_first; //n(1+n)
//	double b_first_minuend = 2.0 * input_argument.segment_length_first - 1.0;//2n-1
//	double b_second_minuend = 2.0 * input_argument.segment_length_second - 1.0;//2n-1
//
//	int j = 0;
//	int endOfLongSegment = 0, indexOfLongSegment = 0;
//	/*int account[1280];
//	int* pointer = account;
//	int index = 0;*/
//
//	double a_sum = NULL;
//	double b_sum = NULL;
//	double t = NULL;
//	for (int i = 1; i <= input_argument.point_dimension; i++) {
//		//cout<<"i: "<<i << endl;
//		a_sum = 0;
//		b_sum = 0;
//		t = 0;
//
//		if (i <= input_argument.remainder) {//first part
//			for (j = int((i - 1) * input_argument.segment_length_first); j < i * input_argument.segment_length_first; j++) {
//				//cout << "    j: "<<j << endl;
//				a_sum += (t - (input_argument.segment_length_second / 2.0)) * original_time_series[j];
//				b_sum += ((2.0 * input_argument.segment_length_first - 1.0) - t * 3.0) * original_time_series[j];
//				/*account[index] = j;
//				index++;*/
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//y[i - 1] = sum / input_argument.segment_length_first;
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			//pla.a[i - 1] = 12.0 * a_sum / a_first_divisor;
//			//pla.b[i - 1] = 2.0 * b_sum / b_first_divisor;
//			temp_pla.a = new DataType;
//			temp_pla.b = new DataType;
//			*temp_pla.a = 12.0 * a_sum / (input_argument.segment_length_second * input_argument.segment_length_first * (input_argument.segment_length_first + 1.0));
//			*temp_pla.b = 2.0 * b_sum / (input_argument.segment_length_first * (1.0 + input_argument.segment_length_first));
//			pla.add(temp_pla);
//			/*===================================================================*/
//			pla[i - 1].segment_width = input_argument.segment_length_first;
//			pla[i - 1].right_endpoint = i * input_argument.segment_length_first - 1;
//			/*...................................................................*/
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//			indexOfLongSegment = i;
//			endOfLongSegment = int(i * input_argument.segment_length_first);
//		}
//		else {//second part
//			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * input_argument.segment_length_second); j < endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second; j++) {
//				assert(j >= 0);
//				//cout << "    j: " << j << endl;
//				a_sum += (t - ((input_argument.segment_length_second - 1.0) / 2.0)) * original_time_series[j];
//				b_sum += ((2.0 * input_argument.segment_length_second - 1.0) - t * 3.0) * original_time_series[j];
//				//account[index] = j;
//				////cout << account[index] << endl;
//				//if (index!=0) assert(account[index]-account[index-1]==1);
//				//index++;
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			//pla.a[i - 1] = 12.0 * a_sum / a_second_divisor;
//			//pla.b[i - 1] = 2.0 * b_sum / b_second_divisor;
//			temp_pla.a = new DataType;
//			temp_pla.b = new DataType;
//			*temp_pla.a = 12.0 * a_sum / (input_argument.segment_length_second * (input_argument.segment_length_second - 1.0) * input_argument.segment_length_first);
//			*temp_pla.b = 2.0 * b_sum / (input_argument.segment_length_second * input_argument.segment_length_first);
//			pla.add(temp_pla);
//			/*===================================================================*/
//			pla[i - 1].segment_width = input_argument.segment_length_second;
//			pla[i - 1].right_endpoint = endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second - 1;
//			/*...................................................................*/
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//		}
//	}
//
//	/*for (auto&&i : pla) {190501
//		cout << *i.a<<" "<< *i.b << endl;
//	}*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(account[i] - account[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//	/*APCA_KNN_QUAL::printArray(pointer,1280);*/
//
//}
//
////************************************
//// Method:getPLAVectorPreMemory
//// Qualifier:Pre-define memory of vector
//// Input:
//// Output:
//// date:190611
//// author:
////************************************
//TEMPLATE
//void PLA_QUAL::getPLAVectorPreMemory(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, vector<PLA>& const pla) {//190611
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.time_series_length != INF && input_argument.point_dimension != INF && input_argument.segment_length_first != INF && input_argument.segment_length_second != INF);
//	//printf("getPLA()\n");
//	//assert(pla.size() == input_argument.point_dimension);
//
//	/*cout << "input_argument.time_series_length: " << input_argument.time_series_length << endl;
//	cout << "input_argument.point_dimension: " << input_argument.point_dimension << endl;
//	cout << "input_argument.remainder: " << input_argument.remainder << endl;
//	cout << "input_argument.segment_length_first: " << input_argument.segment_length_first << endl;
//	cout << "input_argument.segment_length_second: " << input_argument.segment_length_second << endl;*/
//	/*=========================Initial PLA Vector=======================================*/
//	PLA temp_pla;
//	temp_pla.a = new DataType;
//	temp_pla.b = new DataType;
//	pla.resize(input_argument.point_dimension, PLA());
//	for (auto&& au : pla) {
//		au.a = new DataType;
//		au.b = new DataType;
//	}
//	/*................................................................................*/
//
//	//a
//	double a_first_divisor = input_argument.segment_length_second * input_argument.segment_length_first * (input_argument.segment_length_first + 1.0); //n(n-1)(n+1)
//	double a_second_divisor = input_argument.segment_length_second * (input_argument.segment_length_second - 1.0) * input_argument.segment_length_first;  //n(n-1)(n+1)
//	double a_first_minuend = input_argument.segment_length_second / 2.0;//(n-1)/2
//	double a_second_minuend = (input_argument.segment_length_second - 1.0) / 2.0;//(n-1)/2
//	  //b
//	double b_first_divisor = input_argument.segment_length_first * (1.0 + input_argument.segment_length_first); //n(1+n)
//	double b_second_divisor = input_argument.segment_length_second * input_argument.segment_length_first; //n(1+n)
//	double b_first_minuend = 2.0 * input_argument.segment_length_first - 1.0;//2n-1
//	double b_second_minuend = 2.0 * input_argument.segment_length_second - 1.0;//2n-1
//
//	int j = 0;
//	int endOfLongSegment = 0, indexOfLongSegment = 0;
//	/*int account[1280];
//	int* pointer = account;
//	int index = 0;*/
//
//	double a_sum = NULL;
//	double b_sum = NULL;
//	double t = NULL;
//	for (int i = 1; i <= input_argument.point_dimension; i++) {
//		//cout<<"i: "<<i << endl;
//		a_sum = 0;
//		b_sum = 0;
//		t = 0;
//
//		if (i <= input_argument.remainder) {//first part
//			/*===================================================================*/
//			pla[i - 1].segment_width = input_argument.segment_length_first;
//			pla[i - 1].right_endpoint = i * input_argument.segment_length_first - 1;
//			/*...................................................................*/
//			for (j = int((i - 1) * input_argument.segment_length_first); j < i * input_argument.segment_length_first; j++) {
//				//cout << "    j: "<<j << endl;
//				a_sum += (t - a_first_minuend) * original_time_series[j];
//				b_sum += (b_first_minuend - t * 3.0) * original_time_series[j];
//				/*account[index] = j;
//				index++;*/
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//y[i - 1] = sum / input_argument.segment_length_first;
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			//pla.a[i - 1] = 12.0 * a_sum / a_first_divisor;
//			//pla.b[i - 1] = 2.0 * b_sum / b_first_divisor;
//			//temp_pla.a = new DataType;
//			//temp_pla.b = new DataType;
//			*temp_pla.a = 12.0 * a_sum / a_first_divisor;
//			*temp_pla.b = 2.0 * b_sum / b_first_divisor;
//			//pla.push_back(temp_pla);
//
//			*pla[i - 1].a = *temp_pla.a;
//			*pla[i - 1].b = *temp_pla.b;
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//			indexOfLongSegment = i;
//			endOfLongSegment = int(i * input_argument.segment_length_first);
//		}
//		else {//second part
//			/*===================================================================*/
//			pla[i - 1].segment_width = input_argument.segment_length_second;
//			pla[i - 1].right_endpoint = endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second - 1;
//			/*...................................................................*/
//			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * input_argument.segment_length_second); j < endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second; j++) {
//#ifdef _DEBUG
//				assert(j >= 0);
//#endif
//				//cout << "    j: " << j << endl;
//				a_sum += (t - a_second_minuend) * original_time_series[j];
//				b_sum += (b_second_minuend - t * 3.0) * original_time_series[j];
//				//account[index] = j;
//				////cout << account[index] << endl;
//				//if (index!=0) assert(account[index]-account[index-1]==1);
//				//index++;
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			//pla.a[i - 1] = 12.0 * a_sum / a_second_divisor;
//			//pla.b[i - 1] = 2.0 * b_sum / b_second_divisor;
//			//temp_pla.a = new DataType;
//			//temp_pla.b = new DataType;
//			*temp_pla.a = 12.0 * a_sum / a_second_divisor;
//			*temp_pla.b = 2.0 * b_sum / b_second_divisor;
//			//pla.push_back(temp_pla);
//
//			*pla[i - 1].a = *temp_pla.a;
//			*pla[i - 1].b = *temp_pla.b;
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//		}
//	}
//
//	/*for (auto&&i : pla) {///190501
//		cout << *i.a<<" "<< *i.b << endl;
//	}*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(account[i] - account[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//	/*APCA_KNN_QUAL::printArray(pointer,1280);*/
//}
//
//TEMPLATE
//void PLA_QUAL::getPLA(const typename TOOL::INPUT_ARGUMENT& input_argument, const DataType* original_time_series, const PLA& pla) {
//	//printf("getPLA()\n");
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//	/*cout << "input_argument.time_series_length: " << input_argument.time_series_length << endl;
//	cout << "input_argument.point_dimension: " << input_argument.point_dimension << endl;
//	cout << "input_argument.remainder: " << input_argument.remainder << endl;
//	cout << "input_argument.segment_length_first: " << input_argument.segment_length_first << endl;
//	cout << "input_argument.segment_length_second: " << input_argument.segment_length_second << endl;*/
//
//	//a
//	double a_first_divisor = input_argument.segment_length_second * input_argument.segment_length_first * (input_argument.segment_length_first + 1.0); //n(n-1)(n+1)
//	double a_second_divisor = input_argument.segment_length_second * (input_argument.segment_length_second - 1.0) * input_argument.segment_length_first;  //n(n-1)(n+1)
//	double a_first_minuend = input_argument.segment_length_second / 2.0;//(n-1)/2
//	double a_second_minuend = (input_argument.segment_length_second - 1.0) / 2.0;//(n-1)/2
//	  //b
//	double b_first_divisor = input_argument.segment_length_first * (1.0 + input_argument.segment_length_first); //n(1+n)
//	double b_second_divisor = input_argument.segment_length_second * input_argument.segment_length_first; //n(1+n)
//	double b_first_minuend = 2.0 * input_argument.segment_length_first - 1.0;//2n-1
//	double b_second_minuend = 2.0 * input_argument.segment_length_second - 1.0;//2n-1
//
//	int j = 0;
//	int endOfLongSegment = 0, indexOfLongSegment = 0;
//	/*int account[1280];
//	int* pointer = account;
//	int index = 0;*/
//
//	double a_sum = NULL;
//	double b_sum = NULL;
//	double t = NULL;
//	for (int i = 1; i <= input_argument.point_dimension; i++) {
//		//cout<<"i: "<<i << endl;
//		a_sum = 0;
//		b_sum = 0;
//		t = 0;
//		if (i <= input_argument.remainder) {//first part
//			for (j = int((i - 1) * input_argument.segment_length_first); j < i * input_argument.segment_length_first; j++) {
//				//cout << "    j: "<<j << endl;
//				a_sum += (t - a_first_minuend) * original_time_series[j];
//				b_sum += (b_first_minuend - t * 3.0) * original_time_series[j];
//				/*account[index] = j;
//				index++;*/
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//y[i - 1] = sum / input_argument.segment_length_first;
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			pla.a[i - 1] = 12.0 * a_sum / a_first_divisor;
//			pla.b[i - 1] = 2.0 * b_sum / b_first_divisor;
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//			indexOfLongSegment = i;
//			endOfLongSegment = int(i * input_argument.segment_length_first);
//		}
//		else {//second part
//			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * input_argument.segment_length_second); j < endOfLongSegment + (i - indexOfLongSegment) * input_argument.segment_length_second; j++) {
//				assert(j >= 0);
//				//cout << "    j: " << j << endl;
//				a_sum += (t - a_second_minuend) * original_time_series[j];
//				b_sum += (b_second_minuend - t * 3.0) * original_time_series[j];
//				//account[index] = j;
//				////cout << account[index] << endl;
//				//if (index!=0) assert(account[index]-account[index-1]==1);
//				//index++;
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			pla.a[i - 1] = 12.0 * a_sum / a_second_divisor;
//			pla.b[i - 1] = 2.0 * b_sum / b_second_divisor;
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//		}
//	}
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(account[i] - account[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//	/*APCA_KNN_QUAL::printArray(pointer,1280);*/
//}
//
//TEMPLATE
//void PLA_QUAL::getPLA(const int& time_series_length, const int& segment_number, const DataType* original_time_series, const PLA& pla) {
//	//printf("getPLA()\n");
//	assert(pla.segmentNum == segment_number);
//	/*cout << "input_argument.time_series_length: " << input_argument.time_series_length << endl;
//	cout << "input_argument.point_dimension: " << input_argument.point_dimension << endl;
//	cout << "input_argument.remainder: " << input_argument.remainder << endl;
//	cout << "input_argument.segment_length_first: " << input_argument.segment_length_first << endl;
//	cout << "input_argument.segment_length_second: " << input_argument.segment_length_second << endl;*/
//
//	double remainder = int(time_series_length) % int(segment_number);//For PLA
//	double integerDividend = time_series_length - remainder;
//	double segment_length_second = integerDividend / segment_number;
//	double segment_length_first = segment_length_second + 1;
//	assert(segment_length_second > 1);//l(l-1)(l+1), so l != 1
//
//	//a
//	double a_first_divisor = segment_length_second * segment_length_first * (segment_length_first + 1.0); //n(n-1)(n+1)
//	double a_second_divisor = segment_length_second * (segment_length_second - 1.0) * segment_length_first;  //n(n-1)(n+1)
//	double a_first_minuend = segment_length_second / 2.0;//(n-1)/2
//	double a_second_minuend = (segment_length_second - 1.0) / 2.0;//(n-1)/2
//	  //b
//	double b_first_divisor = segment_length_first * (1.0 + segment_length_first); //n(1+n)
//	double b_second_divisor = segment_length_second * segment_length_first; //n(1+n)
//	double b_first_minuend = 2.0 * segment_length_first - 1.0;//2n-1
//	double b_second_minuend = 2.0 * segment_length_second - 1.0;//2n-1
//
//	int j = 0;
//	int endOfLongSegment = 0, indexOfLongSegment = 0;
//	/*int account[1280];
//	int* pointer = account;
//	int index = 0;*/
//
//	double a_sum = NULL;
//	double b_sum = NULL;
//	double t = NULL;
//	for (int i = 1; i <= pla.segmentNum; i++) {
//		//cout<<"i: "<<i << endl;
//		a_sum = 0;
//		b_sum = 0;
//		t = 0;
//		if (i <= remainder) {//first part
//			for (j = int((i - 1) * segment_length_first); j < i * segment_length_first; j++) {
//				//cout << "    j: "<<j << endl;
//				a_sum += (t - a_first_minuend) * original_time_series[j];
//				b_sum += (b_first_minuend - t * 3.0) * original_time_series[j];
//				/*account[index] = j;
//				index++;*/
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//y[i - 1] = sum / segment_length_first;
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			pla.a[i - 1] = 12.0 * a_sum / a_first_divisor;
//			pla.b[i - 1] = 2.0 * b_sum / b_first_divisor;
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//			indexOfLongSegment = i;
//			endOfLongSegment = int(i * segment_length_first);
//		}
//		else {//second part
//			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * segment_length_second); j < endOfLongSegment + (i - indexOfLongSegment) * segment_length_second; j++) {
//				assert(j >= 0);
//				//cout << "    j: " << j << endl;
//				a_sum += (t - a_second_minuend) * original_time_series[j];
//				b_sum += (b_second_minuend - t * 3.0) * original_time_series[j];
//				//cout<< b_sum <<endl;
//				//account[index] = j;
//				////cout << account[index] << endl;
//				//if (index!=0) assert(account[index]-account[index-1]==1);
//				//index++;
//				t++;
//				//cout << "x[" << j << "] = " << x[j] << " ";
//			}
//			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
//			pla.a[i - 1] = 12.0 * a_sum / a_second_divisor;
//			pla.b[i - 1] = 2.0 * b_sum / b_second_divisor;
//			//cout << endl << "v[" << i - 1 << "] = " << italicC->v[i - 1] << ", r[" << i - 1 << "] = " << italicC->r[i - 1] << endl;
//		}
//	}
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(account[i] - account[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//	/*APCA_KNN_QUAL::printArray(pointer,1280);*/
//}
//
//TEMPLATE
//double& PLA_QUAL::getPLADistance(const INPUT_ARGUMENT& input_argument, const PLA& pla_array, const PLA& pla_array_qeury, double& pla_distance) {
//
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//
//	double multipulier_two_first = input_argument.segment_length_first * input_argument.segment_length_second;//l(l-1)
//	double multipulier_one_first = (input_argument.segment_length_first * 2 - 1) * multipulier_two_first / 6;//l(l-1)(2l-1)/6
//
//	double multipulier_two_second = input_argument.segment_length_second * (input_argument.segment_length_second - 1);//l(l-1)
//	double multipulier_one_second = (input_argument.segment_length_second * 2 - 1) * multipulier_two_second / 6;//l(l-1)(2l-1)/6
//
//	double sum = 0;
//
//	for (int i = 0; i < input_argument.point_dimension; i++) {
//		//cout << endl;
//		if (i < input_argument.remainder) {//first part
//			sum += multipulier_one_first * (pla_array.a[i] - pla_array_qeury.a[i]) * (pla_array.a[i] - pla_array_qeury.a[i]) + multipulier_two_first * (pla_array.a[i] - pla_array_qeury.a[i]) * (pla_array.b[i] - pla_array_qeury.b[i]) + input_argument.segment_length_first * (pla_array.b[i] - pla_array_qeury.b[i]) * (pla_array.b[i] - pla_array_qeury.b[i]);
//			//cout << sum << endl;
//		}
//		else {//second part
//			sum += multipulier_one_second * (pla_array.a[i] - pla_array_qeury.a[i]) * (pla_array.a[i] - pla_array_qeury.a[i]) + multipulier_two_second * (pla_array.a[i] - pla_array_qeury.a[i]) * (pla_array.b[i] - pla_array_qeury.b[i]) + input_argument.segment_length_second * (pla_array.b[i] - pla_array_qeury.b[i]) * (pla_array.b[i] - pla_array_qeury.b[i]);
//			//cout << sum << endl;
//		}
//	}
//
//	pla_distance = sqrt(sum);
//
//	//cout << pla_distance << endl;
//	return pla_distance;
//}
//
//TEMPLATE
//double& PLA_QUAL::getPLADistance(const int& time_series_length, const PLA& pla_array, const PLA& pla_array_qeury, double& pla_distance) {
//	assert(pla_array_qeury.segmentNum == pla_array.segmentNum);
//
//	double remainder = int(time_series_length) % int(pla_array.segmentNum);//For PLA
//	double integerDividend = time_series_length - remainder;
//	double segment_length_second = integerDividend / pla_array.segmentNum;
//	double segment_length_first = segment_length_second + 1;
//	assert(segment_length_second > 1);//l(l-1)(l+1), so l != 1
//
//	double multipulier_two_first = segment_length_first * segment_length_second;//l(l-1)
//	double multipulier_one_first = (segment_length_first * 2 - 1) * multipulier_two_first / 6;//l(l-1)(2l-1)/6
//
//	double multipulier_two_second = segment_length_second * (segment_length_second - 1);//l(l-1)
//	double multipulier_one_second = (segment_length_second * 2 - 1) * multipulier_two_second / 6;//l(l-1)(2l-1)/6
//
//	double sum = 0;
//
//	for (int i = 0; i < pla_array.segmentNum; i++) {
//		//cout << endl;
//		if (i < remainder) {//first part
//			sum += multipulier_one_first * (pla_array.a[i] - pla_array_qeury.a[i]) * (pla_array.a[i] - pla_array_qeury.a[i]) + multipulier_two_first * (pla_array.a[i] - pla_array_qeury.a[i]) * (pla_array.b[i] - pla_array_qeury.b[i]) + segment_length_first * (pla_array.b[i] - pla_array_qeury.b[i]) * (pla_array.b[i] - pla_array_qeury.b[i]);
//			//cout << sum << endl;
//		}
//		else {//second part
//			sum += multipulier_one_second * (pla_array.a[i] - pla_array_qeury.a[i]) * (pla_array.a[i] - pla_array_qeury.a[i]) + multipulier_two_second * (pla_array.a[i] - pla_array_qeury.a[i]) * (pla_array.b[i] - pla_array_qeury.b[i]) + segment_length_second * (pla_array.b[i] - pla_array_qeury.b[i]) * (pla_array.b[i] - pla_array_qeury.b[i]);
//			//cout << sum << endl;
//		}
//	}
//
//	pla_distance = sqrt(sum);
//
//	//cout << pla_distance << endl;
//	return pla_distance;
//}
//
//TEMPLATE
//inline double PLA_QUAL::getPointDistanceSquare(const DataType& point_a_x, const DataType& point_a_y, const DataType& point_b_x, const DataType& point_b_y) {
//	return (point_a_x - point_b_x) * (point_a_x - point_b_x) + (point_a_y - point_b_y) * (point_a_y - point_b_y);
//}
//
//TEMPLATE
//double PLA_QUAL::getNearestDistance(const DataType& point_a_x, const DataType& point_a_y, const DataType& point_b_x, const DataType& point_b_y, const DataType& point_q_x, const DataType& point_q_y) {
//	//----------point on the segment--------------------
//	double a = NULL, b = NULL, c = NULL;
//	double a_sqrt = NULL, b_sqrt = NULL, c_sqrt = NULL;
//	a = getPointDistanceSquare(point_b_x, point_b_y, point_q_x, point_q_y);//a
//	a_sqrt = sqrt(a);
//	if (a_sqrt <= 0.00001)
//		return 0.0;
//	b = getPointDistanceSquare(point_a_x, point_a_y, point_q_x, point_q_y);//b
//	b_sqrt = sqrt(b);
//	if (b_sqrt <= 0.00001)
//		return 0.0;
//	c = getPointDistanceSquare(point_a_x, point_a_y, point_b_x, point_b_y);//c
//	c_sqrt = sqrt(c);
//	if (c_sqrt <= 0.00001)
//		return a_sqrt;//if point_a == point_b £¬return distance
//	//------------------------------
//
//	if (a >= b + c) {//--------p_q_x<p_a_x--------
//		//cout << "aaa" << endl;
//		return b;
//	}     //if obtuse angle, return b
//	if (b >= a + c) {//--------p_q_x>p_b_x-------
//		//cout << "bb" << endl;
//		return a;
//	}      //if obtuse angle, return a
//
//					   //-------p_a_x<p_q_x<p_b_x && p_q_y != p_a_y&&p_b_y -----
//	double l = (a_sqrt + b_sqrt + c_sqrt) / 2.0;     //circumference /2
//	double s = sqrt(l * (l - a_sqrt) * (l - b_sqrt) * (l - c_sqrt));  //Heron's formula
//
//	//cout << "cc" << endl;
//	return pow(2.0 * s / c_sqrt, 2.0); // ^2
//}
//
////ax+by+c=0,
////(m,n),
////((b*b*m-a*b*n-a*c)/(a*a+b*b),(a*a*n-a*b*m-b*c)/(a*a+b*b))
////n=0,c=0
//
//TEMPLATE
//double& PLA_QUAL::getPLAMBRSegmentDistance(const typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const int& segment_id, const PLA& pla_array_qeury, double& pla_MBR_segment_distance) {
//	//assert(MBR.segmentNum == 2 * pla_array_qeury.segmentNum);
//
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//
//	POINT point;
//	double sqrt_segment_legngth = sqrt(double(input_argument.segment_length_second));//l^0.5
//	double perpendicular_b = (double(input_argument.segment_length_second) - 1.0) / 2.0; // (l-1)/2;  -b
//	double argument_u_A = sqrt_segment_legngth * perpendicular_b;//l^0.5 * (l-1)/2
//
//	double perpendicular_a = sqrt((input_argument.segment_length_second * input_argument.segment_length_second - 1.0) / 12.0);//((l^2-1)/12)^0.5 a
//	double argument_v_B = perpendicular_a / perpendicular_b;//    ((l^2-1)/12)^0.5  /  (l-1)/2
//
//	point.u_A_1 = argument_u_A * (MBR.m_min[segment_id * 2] - pla_array_qeury.a[segment_id]);//a_min amin-qa
//	point.u_A_2 = argument_u_A * (MBR.m_max[segment_id * 2] - pla_array_qeury.a[segment_id]);//a_max amax-qa
//	point.v_A_1 = argument_v_B * point.u_A_1;
//	point.v_A_2 = argument_v_B * point.u_A_2;
//
//	point.u_B_1 = sqrt_segment_legngth * (pla_array_qeury.b[segment_id] - MBR.m_max[segment_id * 2 + 1]);//B1
//	point.u_B_2 = sqrt_segment_legngth * (pla_array_qeury.b[segment_id] - MBR.m_min[segment_id * 2 + 1]);//B2
//
//	double perpendicular_denominator = perpendicular_a * perpendicular_a + perpendicular_b * perpendicular_b;//a^2+b^2
//	point.u_B_C_1 = -perpendicular_b * -perpendicular_b * point.u_B_1 / perpendicular_denominator; //b^2*uB_1 / (a^2+b^2)
//	point.v_B_C_1 = perpendicular_a * perpendicular_b * point.u_B_1 / perpendicular_denominator; //-ab*uB_1/(a^2+b^2)
//	point.u_B_C_2 = -perpendicular_b * -perpendicular_b * point.u_B_2 / perpendicular_denominator; //b*b*uB_2 / (a^2+b^2)
//	point.v_B_C_2 = perpendicular_a * perpendicular_b * point.u_B_2 / perpendicular_denominator; //a*-b*uB_2 / (a^2+b^2)
//
//	if (point.u_A_2 < 0) {//case 1   Line segment A1A2 is completely contained in the third quadrant of the u - v space.
//		if (point.u_B_2 < point.u_A_2) {//case 1.3-1.5
//			//cout << "a" << endl;
//			return pla_MBR_segment_distance = getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_2, point.v_B_2);
//		}
//		else if (point.u_B_2 >= point.u_A_2) {//case 1.1-1.2
//			//cout << "b" << endl;
//			return pla_MBR_segment_distance = getNearestDistance(point.u_B_1, point.v_B_1, point.u_B_2, point.v_B_2, point.u_A_2, point.v_A_2);
//		}
//		else assert(0);
//	}
//	else if (point.u_A_1 <= 0 && point.u_A_2 >= 0) {//case 2: Line segment A1A2 is partially contained in the first and third quadrants of the u - v space
//		if (point.u_B_1 > 0) {
//			if (point.u_B_C_1 >= point.u_A_2) {//case 2.1; uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 > 0, uC ¡Ý uA2
//				//cout << "c" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);//|A2B1|^2 // case2.1
//			}
//			else if (point.u_B_C_1 < point.u_A_2) {//case 2.2; uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 > 0, uC < uA2
//				//cout << "d" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_B_1, point.v_B_1, point.u_B_C_1, point.v_B_C_1);//|B1C|^2 //case 2.2
//			}
//			else assert(0);
//		}
//		else if (point.u_B_1 <= 0 && point.u_B_2 >= 0) {//uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 ¡Ü 0, uB2 ¡Ý 0
//			//cout <<"dist: "<< point.u_B_1 << ", " << point.u_B_2 << endl;
//			//cout << "e" << endl;
//			return pla_MBR_segment_distance = 0;//case 2.3
//		}
//		else if (point.u_B_2 < 0) {
//			if (point.u_B_C_2 > point.u_A_1) {//case2.4: uA1 ¡Ü 0, uA2 ¡Ý 0, uB2 < 0, uC > uA1
//				//cout <<"dist: "<< point.u_B_1 << ", " << point.u_B_2 << endl;
//				//cout << "f" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_B_2, point.v_B_2, point.u_B_C_2, point.v_B_C_2);;//|B2C|^2 //case 2.4
//			}
//			else if (point.u_B_C_2 <= point.u_A_1) {//case: 2.5 uA1 ¡Ü 0, uA2 ¡Ý 0, uB2 < 0, uC ¡Ü uA1
//				//cout << "g" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_B_2, point.v_B_2);//|A1B1|^2 //case 2.5
//			}
//			else assert(0);
//		}
//		else assert(0);
//	}
//	else if (point.u_A_1 > 0) {//case 3: Line segment A1A2 is completely contained in the first quadrant of the u-v space.
//		if (point.u_B_1 > point.u_A_1) {//case 3.1-3.3
//			//cout << sqrt(getPointDistanceSquare(point.u_B_2, 0.0, point.u_B_C_2, point.v_B_C_2)) << endl;
//			return pla_MBR_segment_distance = getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);
//		}
//		else if (point.u_B_1 <= point.u_A_1) {//case 3.4-3.5
//			//cout << "i" << endl;
//			return pla_MBR_segment_distance = getNearestDistance(point.u_B_1, point.v_B_1, point.u_B_2, point.v_B_2, point.u_A_1, point.v_A_1);
//		}
//		else assert(0);
//	}
//	else assert(0);
//
//	assert(0);
//}
//
//TEMPLATE
//double& PLA_QUAL::getPLAMBRSegmentDistance(INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const int& segment_id, const PLA& pla_array_qeury, double& pla_MBR_segment_distance) {
//	//assert(MBR.segmentNum == 2 * pla_array_qeury.segmentNum);
//
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//
//	POINT point;
//	double sqrt_segment_legngth = sqrt(double(input_argument.segment_length_second));//l^0.5
//	double perpendicular_b = (double(input_argument.segment_length_second) - 1.0) / 2.0; // (l-1)/2;  -b
//	double argument_u_A = sqrt_segment_legngth * perpendicular_b;//l^0.5 * (l-1)/2
//
//	double perpendicular_a = sqrt((input_argument.segment_length_second * input_argument.segment_length_second - 1.0) / 12.0);//((l^2-1)/12)^0.5 a
//	double argument_v_B = perpendicular_a / perpendicular_b;//    ((l^2-1)/12)^0.5  /  (l-1)/2
//
//	point.u_A_1 = argument_u_A * (MBR.m_min[segment_id * 2] - pla_array_qeury.a[segment_id]);//a_min amin-qa
//	point.u_A_2 = argument_u_A * (MBR.m_max[segment_id * 2] - pla_array_qeury.a[segment_id]);//a_max amax-qa
//	point.v_A_1 = argument_v_B * point.u_A_1;
//	point.v_A_2 = argument_v_B * point.u_A_2;
//
//	point.u_B_1 = sqrt_segment_legngth * (pla_array_qeury.b[segment_id] - MBR.m_max[segment_id * 2 + 1]);//B1
//	point.u_B_2 = sqrt_segment_legngth * (pla_array_qeury.b[segment_id] - MBR.m_min[segment_id * 2 + 1]);//B2
//
//	double perpendicular_denominator = perpendicular_a * perpendicular_a + perpendicular_b * perpendicular_b;//a^2+b^2
//	point.u_B_C_1 = -perpendicular_b * -perpendicular_b * point.u_B_1 / perpendicular_denominator; //b^2*uB_1 / (a^2+b^2)
//	point.v_B_C_1 = perpendicular_a * perpendicular_b * point.u_B_1 / perpendicular_denominator; //-ab*uB_1/(a^2+b^2)
//	point.u_B_C_2 = -perpendicular_b * -perpendicular_b * point.u_B_2 / perpendicular_denominator; //b*b*uB_2 / (a^2+b^2)
//	point.v_B_C_2 = perpendicular_a * perpendicular_b * point.u_B_2 / perpendicular_denominator; //a*-b*uB_2 / (a^2+b^2)
//
//	if (point.u_A_2 < 0) {//case 1   Line segment A1A2 is completely contained in the third quadrant of the u - v space.
//		if (point.u_B_2 < point.u_A_2) {//case 1.3-1.5
//			//cout << "a" << endl;
//			return pla_MBR_segment_distance = getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_2, point.v_B_2);
//		}
//		else if (point.u_B_2 >= point.u_A_2) {//case 1.1-1.2
//			//cout << "b" << endl;
//			return pla_MBR_segment_distance = getNearestDistance(point.u_B_1, point.v_B_1, point.u_B_2, point.v_B_2, point.u_A_2, point.v_A_2);
//		}
//		else assert(0);
//	}
//	else if (point.u_A_1 <= 0 && point.u_A_2 >= 0) {//case 2: Line segment A1A2 is partially contained in the first and third quadrants of the u - v space
//		if (point.u_B_1 > 0) {
//			if (point.u_B_C_1 >= point.u_A_2) {//case 2.1; uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 > 0, uC ¡Ý uA2
//				//cout << "c" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);//|A2B1|^2 // case2.1
//			}
//			else if (point.u_B_C_1 < point.u_A_2) {//case 2.2; uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 > 0, uC < uA2
//				//cout << "d" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_B_1, point.v_B_1, point.u_B_C_1, point.v_B_C_1);//|B1C|^2 //case 2.2
//			}
//			else assert(0);
//		}
//		else if (point.u_B_1 <= 0 && point.u_B_2 >= 0) {//uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 ¡Ü 0, uB2 ¡Ý 0
//			//cout <<"dist: "<< point.u_B_1 << ", " << point.u_B_2 << endl;
//			//cout << "e" << endl;
//			return pla_MBR_segment_distance = 0;//case 2.3
//		}
//		else if (point.u_B_2 < 0) {
//			if (point.u_B_C_2 > point.u_A_1) {//case2.4: uA1 ¡Ü 0, uA2 ¡Ý 0, uB2 < 0, uC > uA1
//				//cout <<"dist: "<< point.u_B_1 << ", " << point.u_B_2 << endl;
//				//cout << "f" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_B_2, point.v_B_2, point.u_B_C_2, point.v_B_C_2);;//|B2C|^2 //case 2.4
//			}
//			else if (point.u_B_C_2 <= point.u_A_1) {//case: 2.5 uA1 ¡Ü 0, uA2 ¡Ý 0, uB2 < 0, uC ¡Ü uA1
//				//cout << "g" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_B_2, point.v_B_2);//|A1B1|^2 //case 2.5
//			}
//			else assert(0);
//		}
//		else assert(0);
//	}
//	else if (point.u_A_1 > 0) {//case 3: Line segment A1A2 is completely contained in the first quadrant of the u-v space.
//		if (point.u_B_1 > point.u_A_1) {//case 3.1-3.3
//			//cout << sqrt(getPointDistanceSquare(point.u_B_2, 0.0, point.u_B_C_2, point.v_B_C_2)) << endl;
//			return pla_MBR_segment_distance = getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);
//		}
//		else if (point.u_B_1 <= point.u_A_1) {//case 3.4-3.5
//			//cout << "i" << endl;
//			return pla_MBR_segment_distance = getNearestDistance(point.u_B_1, point.v_B_1, point.u_B_2, point.v_B_2, point.u_A_1, point.v_A_1);
//		}
//		else assert(0);
//	}
//	else assert(0);
//
//	assert(0);
//}
//
//TEMPLATE//from paper
//double& PLA_QUAL::getPLAMBRSegmentDistanceBase(const RTREE::Rect& MBR, const int& segment_id, const PLA& pla_array_qeury, double& pla_MBR_segment_distance) {
//	//assert(MBR.segmentNum == 2 * pla_array_qeury.segmentNum);
//
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//	POINT point;
//	double sqrt_segment_legngth = sqrt(double(input_argument.segment_length_second));//l^0.5
//	double perpendicular_b = (double(input_argument.segment_length_second) - 1.0) / 2.0; // (l-1)/2;  -b
//	double argument_u_A = sqrt_segment_legngth * perpendicular_b;//l^0.5 * (l-1)/2
//
//	double perpendicular_a = sqrt((input_argument.segment_length_second * input_argument.segment_length_second - 1.0) / 12.0);//((l^2-1)/12)^0.5 a
//	double argument_v_B = perpendicular_a / perpendicular_b;//    ((l^2-1)/12)^0.5  /  (l-1)/2
//
//	point.u_A_1 = argument_u_A * (MBR.m_min[segment_id * 2] - pla_array_qeury.a[segment_id]);//a_min amin-qa
//	point.u_A_2 = argument_u_A * (MBR.m_max[segment_id * 2] - pla_array_qeury.a[segment_id]);//a_max amax-qa
//	point.v_A_1 = argument_v_B * point.u_A_1;
//	point.v_A_2 = argument_v_B * point.u_A_2;
//
//	point.u_B_1 = sqrt_segment_legngth * (pla_array_qeury.b[segment_id] - MBR.m_max[segment_id * 2 + 1]);//B1
//	point.u_B_2 = sqrt_segment_legngth * (pla_array_qeury.b[segment_id] - MBR.m_min[segment_id * 2 + 1]);//B2
//
//	double perpendicular_denominator = perpendicular_a * perpendicular_a + perpendicular_b * perpendicular_b;//a^2+b^2
//	point.u_B_C_1 = -perpendicular_b * -perpendicular_b * point.u_B_1 / perpendicular_denominator; //b^2*uB_1 / (a^2+b^2)
//	point.v_B_C_1 = perpendicular_a * perpendicular_b * point.u_B_1 / perpendicular_denominator; //-ab*uB_1/(a^2+b^2)
//	point.u_B_C_2 = -perpendicular_b * -perpendicular_b * point.u_B_2 / perpendicular_denominator; //b*b*uB_2 / (a^2+b^2)
//	point.v_B_C_2 = perpendicular_a * perpendicular_b * point.u_B_2 / perpendicular_denominator; //a*-b*uB_2 / (a^2+b^2)
//
//	if (point.u_A_2 < 0) {//case 1   Line segment A1A2 is completely contained in the third quadrant of the u - v space.
//		if (point.u_B_1 >= point.u_A_2) {//case 1.1; uA2 < 0, uB1 > uA2
//			//cout << "11" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);//|A2B1|2
//		}
//		else if (point.u_B_1 <= point.u_A_2 && point.u_B_2 >= point.u_A_2) {//case 1.2; uA2 < 0, uB1 ¡Ü uA2, uB2 > uA2
//			//cout << "22" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_A_2, 0.0);//|A2C|2
//		}
//		else if (point.u_A_2 > point.u_B_2 && point.u_B_C_2 >= point.u_A_2) {//case 1.3; uA2 < 0, uB2 < uA2, uC ¡Ý uA2
//			//cout << "33" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_B_2, point.v_B_2);//|A2B2|2
//		}
//		else if (point.u_B_2 <= point.u_A_2 && point.u_B_C_2 > point.u_A_1 && point.u_B_C_2 < point.u_A_2) {//case 1.4; uA2 < 0, uB2 ¡Ü uA2, uA1 < uC < uA2
//			//cout << "44" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_B_2, point.v_B_2, point.u_B_C_2, point.v_B_C_2);//|B2C|2
//		}
//		else if (point.u_B_2 <= point.u_A_1 && point.u_B_C_2 <= point.u_A_1) {//case 1.5; uA2 < 0, uB2 < uA1, uC ¡Ü uA1
//			//cout << "55" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_B_2, point.v_B_2);//|A1B2|2
//		}
//		else {
//			cout << "segment_id=" << segment_id << endl;
//			cout << "pla_array_qeury.a[segment_id]=" << pla_array_qeury.a[segment_id] << endl;
//			cout << "pla_array_qeury.b[segment_id]=" << pla_array_qeury.b[segment_id] << endl;
//			cout << "MBR.m_min[segment_id * 2]=" << MBR.m_min[segment_id * 2] << endl;
//			cout << "MBR.m_max[segment_id * 2]=" << MBR.m_max[segment_id * 2] << endl;
//			cout << "MBR.m_min[segment_id * 2+1]=" << MBR.m_min[segment_id * 2 + 1] << endl;
//			cout << "MBR.m_max[segment_id * 2+1]=" << MBR.m_max[segment_id * 2 + 1] << endl;
//
//			cout << "sqrt_segment_legngth=" << sqrt_segment_legngth << endl;//l^0.5
//			cout << "argument_u_A=" << argument_u_A << endl;
//			cout << "argument_v_B=" << argument_v_B << endl;
//			cout << "perpendicular_b=" << perpendicular_b << endl;
//			cout << "perpendicular_a=" << perpendicular_a << endl;
//			cout << "perpendicular_denominator=" << perpendicular_denominator << endl;
//			cout << "ua1=" << point.u_A_1 << endl;
//			cout << "ua2=" << point.u_A_2 << endl;
//			cout << "va1=" << point.v_A_1 << endl;
//			cout << "va2=" << point.v_A_2 << endl;
//			cout << "ub1=" << point.u_B_1 << endl;
//			cout << "ub2=" << point.u_B_2 << endl;
//			cout << "ubc1=" << point.u_B_C_1 << endl;
//			cout << "ubc2=" << point.u_B_C_2 << endl;
//			cout << "vbc1=" << point.v_B_C_1 << endl;
//			cout << "vbc2=" << point.v_B_C_2 << endl;
//
//			cout << "bc1 = " << sqrt(getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_1, 0.0)) << endl;
//			cout << "bc2 = " << sqrt(getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_2, 0.0)) << endl;
//			cout << sqrt(getPointDistanceSquare(point.u_B_1, 0.0, point.u_B_C_1, point.v_B_C_1)) << endl;
//			cout << sqrt(getPointDistanceSquare(point.u_B_2, 0.0, point.u_B_C_2, point.v_B_C_2)) << endl;
//
//			assert(0);
//		}
//	}
//	else if (point.u_A_1 <= 0 && point.u_A_2 >= 0) {//case 2: Line segment A1A2 is partially contained in the first and third quadrants of the u - v space
//		if (point.u_B_1 > 0) {
//			if (point.u_B_C_1 >= point.u_A_2) {//case 2.1; uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 > 0, uC ¡Ý uA2
//				//cout << "66" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);//|A2B1|^2 // case2.1
//			}
//			else if (point.u_B_C_1 < point.u_A_2) {//case 2.2; uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 > 0, uC < uA2
//				//cout << "77" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_B_1, point.v_B_1, point.u_B_C_1, point.v_B_C_1);//|B1C|^2 //case 2.2
//			}
//			else assert(0);
//		}
//		else if (point.u_B_1 <= 0 && point.u_B_2 >= 0) {//uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 ¡Ü 0, uB2 ¡Ý 0
//			//cout << "88" << endl;
//			//cout <<"base: "<<point.u_B_1 << ", " << point.u_B_2 << endl;
//			return pla_MBR_segment_distance = 0;//case 2.3
//		}
//		else if (point.u_B_2 < 0) {
//			if (point.u_B_C_2 > point.u_A_1) {//case2.4: uA1 ¡Ü 0, uA2 ¡Ý 0, uB2 < 0, uC > uA1
//				//cout << "99" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_B_2, point.v_B_2, point.u_B_C_2, point.v_B_C_2);;//|B2C|^2 //case 2.4
//			}
//			else if (point.u_B_C_2 <= point.u_A_1) {//case: 2.5 uA1 ¡Ü 0, uA2 ¡Ý 0, uB2 < 0, uC ¡Ü uA1
//				//cout << "1010" << endl;
//				return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_B_2, point.v_B_2);//|A1B1|^2 //case 2.5
//			}
//			else assert(0);
//		}
//		else assert(0);
//	}
//	else if (point.u_A_1 > 0) {//case 3: Line segment A1A2 is completely contained in the first quadrant of the u-v space.
//		if (point.u_B_2 <= point.u_A_1) {//case 3.5; uA1 > 0, uB2 < uA1
//			//cout << "1111" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_B_2, point.v_B_2);//|A1B2|2
//		}
//		else if (point.u_B_2 > 0 && point.u_B_1 <= point.u_A_1 && point.u_B_2 >= point.u_A_1) {//case 3.4; uA2 < 0, uB1 ¡Ü uA2, uB2 > uA2
//			//cout << "1212" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_A_1, 0.0);//|A1C|2
//		}
//		else if (point.u_A_1 <= point.u_B_1 && point.u_B_C_1 <= point.u_A_1) {//case 3.3; uA1 > 0, uA1 ¡Ü uB1, uC <= uA1
//			//cout << "1313" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_B_1, point.v_B_1);//|A1B1|2
//		}
//		else if (point.u_B_1 > point.u_A_1 && point.u_B_C_1 < point.u_A_2 && point.u_B_C_1 > point.u_A_1) {//case 3.2; uA1 > 0, uB1 > uA1, uA1 < uC < uA2
//			//cout << "1414" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_B_1, point.v_B_1, point.u_B_C_1, point.v_B_C_1);//|B1C|2
//		}
//		else if (point.u_B_1 > point.u_A_2 && point.u_B_C_1 >= point.u_A_2) {//case 3.1; uA1 > 0, uB1 > uA2, uC >= uA2
//			//cout << "1515" << endl;
//			return pla_MBR_segment_distance = getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);//|A2B1|2
//		}
//		else assert(0);
//	}
//	else assert(0);
//
//	assert(0);
//}
//
//TEMPLATE
//double& PLA_QUAL::getPLAMBRDistance(typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const PLA& pla_array_qeury, double& pla_MBR_distance) {
//
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//	double sum = 0;
//	//double pla_MBR_distance_base = NULL;// to evaluate result distance
//	for (int i = 0; i < input_argument.point_dimension; i++) {
//		//pla_MBR_distance_base = getPLAMBRSegmentDistanceBase(MBR, i, pla_array_qeury, pla_MBR_distance_base);//to evaluate result distance
//		sum += getPLAMBRSegmentDistance(input_argument, MBR, i, pla_array_qeury, pla_MBR_distance);
//
//		//assert(pla_MBR_distance_base - pla_MBR_distance<0.0000000001&&pla_MBR_distance_base - pla_MBR_distance>-0.0000000001);//to evaluate result distance
//	}
//	pla_MBR_distance = sqrt(sum);
//	return pla_MBR_distance;
//}
//
//TEMPLATE
//double& PLA_QUAL::getPLAMBRDistance(INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const PLA& pla_array_qeury, double& pla_MBR_distance) {
//
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//	double sum = 0;
//	//double pla_MBR_distance_base = NULL;// to evaluate result distance
//	for (int i = 0; i < input_argument.point_dimension; i++) {
//		//pla_MBR_distance_base = getPLAMBRSegmentDistanceBase(MBR, i, pla_array_qeury, pla_MBR_distance_base);//to evaluate result distance
//		sum += getPLAMBRSegmentDistance(input_argument, MBR, i, pla_array_qeury, pla_MBR_distance);
//
//		//assert(pla_MBR_distance_base - pla_MBR_distance<0.0000000001&&pla_MBR_distance_base - pla_MBR_distance>-0.0000000001);//to evaluate result distance
//	}
//	pla_MBR_distance = sqrt(sum);
//	return pla_MBR_distance;
//}
//
//TEMPLATE//Wrong, for PLA paper, it translate every point at begining
//void PLA_QUAL::approximateOriginalFunctionPLA(const INPUT_ARGUMENT& input_argument, const PLA& pla, DataType*& const approximate_array) {
//	assert("This method has problem, use improved one!!!!");
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//	assert(input_argument.segment_length_second != NULL);
//	//assert(input_argument.time_series_length % input_argument.point_dimension == 0);
//
//	int endOfLongSegment = input_argument.remainder * input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	//double difference = NULL;
//	double sum = 0;
//	double point_id = NULL;
//	int count[1280];
//	int i = 0;
//	int* pointer = count;
//
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		/*for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//		point_id = input_argument.segment_length_second*segment_id + interval_id;
//		difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//		deviation_max = difference > deviation_max ? difference : deviation_max;
//		sum += difference * difference;
//		}*/
//
//		if (segment_id < input_argument.remainder) {
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
//				count[i] = point_id;
//				i++;
//				//difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximate_array[int(point_id)] = pla.a[segment_id] * point_id + pla.b[segment_id];
//				//cout << "difference: "<<difference<<endl;
//			}
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//			assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
//			//assert(indexOfLongSegment== input_argument.remainder-1);
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment) * input_argument.segment_length_second + interval_id;// time series id
//				assert(point_id >= 0);
//				//cout << "id: " << point_id << endl;
//				count[i] = point_id;
//				if (i != 0) {
//					assert(count[i] - count[i - 1] == 1);
//				}
//				i++;
//				//difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximate_array[int(point_id)] = pla.a[segment_id] * point_id + pla.b[segment_id];
//				/*if (difference>50) {
//				cout << pla.a[segment_id] * point_id + pla.b[segment_id]<<" "<< original_array[int(point_id)] << endl;
//				}*/
//				//cout << "difference: " << difference << endl;
//			}
//		}
//	}
//	//APCA_KNN_QUAL::printArray(pointer,1280);
//	/*cout << "approximation:" << endl;
//	APCA_KNN_QUAL::printArray(approximation, 1280);*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(count[i] - count[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//}
//
//TEMPLATE
//void PLA_QUAL::approximateOriginalFunctionPLA1(const INPUT_ARGUMENT& input_argument, const PLA& pla, DataType*& const approximate_array) {
//	assert(input_argument.segment_length_second != NULL);
//	//assert(input_argument.time_series_length % input_argument.point_dimension == 0);
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//	int endOfLongSegment = input_argument.remainder * input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	//double difference = NULL;
//	double sum = 0;
//	double point_id = NULL;
//	int count[1280];
//	int i = 0;
//	int* pointer = count;
//
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		/*for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//		point_id = input_argument.segment_length_second*segment_id + interval_id;
//		difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//		deviation_max = difference > deviation_max ? difference : deviation_max;
//		sum += difference * difference;
//		}*/
//		//cout << "i: " << segment_id <<endl;
//		if (segment_id < input_argument.remainder) {
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
//				//cout << "    j: " << point_id << endl;
//				count[i] = point_id;
//				i++;
//				//difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximate_array[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				//cout << "difference: "<<difference<<endl;
//			}
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//			assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
//			//assert(indexOfLongSegment== input_argument.remainder-1);
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment) * input_argument.segment_length_second + interval_id;// time series id
//				//cout << "    j: " << point_id << endl;
//				assert(point_id >= 0);
//				//cout << "id: " << point_id << endl;
//				count[i] = point_id;
//				if (i != 0) {
//					assert(count[i] - count[i - 1] == 1);
//				}
//				i++;
//				//difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximate_array[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				/*if (difference>50) {
//				cout << pla.a[segment_id] * point_id + pla.b[segment_id]<<" "<< original_array[int(point_id)] << endl;
//				}*/
//				//cout << "difference: " << difference << endl;
//			}
//		}
//	}
//	//APCA_KNN_QUAL::printArray(pointer,1280);
//	/*cout << "approximation:" << endl;
//	APCA_KNN_QUAL::printArray(approximation, 1280);*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(count[i] - count[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//}
//
//TEMPLATE
//void PLA_QUAL::approximateOriginalFunctionPLA1(const typename TOOL::INPUT_ARGUMENT& const input_argument, const PLA& pla, DataType*& const approximate_array) {//right: because PLA paper translate every point to begin.
//	assert(input_argument.segment_length_second != NULL);
//	//assert(input_argument.time_series_length % input_argument.point_dimension == 0);
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//	int endOfLongSegment = input_argument.remainder * input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	//double difference = NULL;
//	double sum = 0;
//	double point_id = NULL;
//	int count[1280];
//	int i = 0;
//	int* pointer = count;
//
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		/*for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//		point_id = input_argument.segment_length_second*segment_id + interval_id;
//		difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//		deviation_max = difference > deviation_max ? difference : deviation_max;
//		sum += difference * difference;
//		}*/
//		//cout << "i: " << segment_id <<endl;
//		if (segment_id < input_argument.remainder) {
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
//				//cout << "    j: " << point_id << endl;
//				count[i] = point_id;
//				i++;
//				//difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximate_array[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				//cout << "difference: "<<difference<<endl;
//			}
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//			assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
//			//assert(indexOfLongSegment== input_argument.remainder-1);
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment) * input_argument.segment_length_second + interval_id;// time series id
//				//cout << "    j: " << point_id << endl;
//				assert(point_id >= 0);
//				//cout << "id: " << point_id << endl;
//				count[i] = point_id;
//				if (i != 0) {
//					assert(count[i] - count[i - 1] == 1);
//				}
//				i++;
//				//difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximate_array[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				/*if (difference>50) {
//				cout << pla.a[segment_id] * point_id + pla.b[segment_id]<<" "<< original_array[int(point_id)] << endl;
//				}*/
//				//cout << "difference: " << difference << endl;
//			}
//		}
//	}
//	//APCA_KNN_QUAL::printArray(pointer,1280);
//	/*cout << "approximation:" << endl;
//	APCA_KNN_QUAL::printArray(approximation, 1280);*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(count[i] - count[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//}
//
////191206 get sum deviation
//TEMPLATE
//double PLA_QUAL::get_pla_sum_deviation(const typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<DataType>& const original_time_series, const PLA& const pla) {
//	assert(input_argument.segment_length_second != NULL);
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//
//	//assert(input_argument.time_series_length % input_argument.point_dimension == 0);
//	double deviation_sum = INF;
//	long double deviation_max=-INF;
//	int endOfLongSegment = input_argument.remainder * input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	double difference = NULL;
//	double sum = 0;
//	double point_id = NULL;
//	int count[1280];
//	int i = 0;
//	int* pointer = count;
//
//	vector<DataType> approximation(input_argument.time_series_length, INF);
//
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		/*for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//			point_id = input_argument.segment_length_second*segment_id + interval_id;
//			difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//			deviation_max = difference > deviation_max ? difference : deviation_max;
//			sum += difference * difference;
//		}*/
//
//		if (segment_id < input_argument.remainder) {
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
//				//cout << "id: " << point_id <<endl;
//				count[i] = point_id;
//				i++;
//				difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_time_series[int(point_id)]);
//				approximation[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				//cout << "difference: "<<difference<<endl;
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//			}
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//			assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
//			//assert(indexOfLongSegment== input_argument.remainder-1);
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment) * input_argument.segment_length_second + interval_id;// time series id
//				assert(point_id >= 0);
//				//cout << "id: " << point_id << endl;
//				count[i] = point_id;
//				if (i != 0) {
//					assert(count[i] - count[i - 1] == 1);
//				}
//				i++;
//				difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_time_series[int(point_id)]);
//				approximation[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				/*if (difference>50) {
//					cout << pla.a[segment_id] * point_id + pla.b[segment_id]<<" "<< original_array[int(point_id)] << endl;
//				}*/
//				//cout << "difference: " << difference << endl;
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//			}
//		}
//	}
//	//APCA_KNN_QUAL::printArray(pointer,1280);
//	/*cout << "approximation:" << endl;
//	APCA_KNN_QUAL::printArray(approximation, 1280);*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(count[i] - count[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//
//
//	deviation_sum = sqrt(sum);
//
//	double test_deviation = TOOL::getDeviation(original_time_series, approximation);
//	assert(float(deviation_sum) == float(test_deviation));
//
//	approximation.clear();
//	approximation.shrink_to_fit();
//	return deviation_sum;
//}
//
//TEMPLATE
//template<typename T, typename Y, typename U, typename T1>
//long double PLA_QUAL::get_pla_sum_deviation(const T& const input_argument, const vector<Y>& const original_time_series, const T1& const pla, U& const result_collection) {
//	assert(input_argument.segment_length_second != NULL);
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//
//	//assert(input_argument.time_series_length % input_argument.point_dimension == 0);
//	long double deviation_sum = INF;
//	long double deviation_max = -INF;
//	int endOfLongSegment = input_argument.remainder * input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	long double difference = NULL;
//	long double sum = 0;
//	double point_id = NULL;
//	int count[1280];
//	int i = 0;
//	int* pointer = count;
//
//	result_collection.sum_deviation = 0;
//	result_collection.max_deviation = 0;
//	result_collection.max_deviation_multiple_width = 0;
//
//	vector<DataType> approximation(input_argument.time_series_length, INF);
//
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		/*for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//			point_id = input_argument.segment_length_second*segment_id + interval_id;
//			difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//			deviation_max = difference > deviation_max ? difference : deviation_max;
//			sum += difference * difference;
//		}*/
//		deviation_max = -INF;
//		if (segment_id < input_argument.remainder) {
//			
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
//				//cout << "id: " << point_id <<endl;
//				count[i] = point_id;
//				i++;
//				difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_time_series[int(point_id)]);
//				approximation[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				//cout << "difference: "<<difference<<endl;
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//			}
//
//			result_collection.max_deviation += deviation_max;
//			result_collection.max_deviation_multiple_width += deviation_max * input_argument.segment_length_first;
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//			assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
//			//assert(indexOfLongSegment== input_argument.remainder-1);
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment) * input_argument.segment_length_second + interval_id;// time series id
//				assert(point_id >= 0);
//				//cout << "id: " << point_id << endl;
//				count[i] = point_id;
//				if (i != 0) {
//					assert(count[i] - count[i - 1] == 1);
//				}
//				i++;
//				difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_time_series[int(point_id)]);
//				approximation[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				/*if (difference>50) {
//					cout << pla.a[segment_id] * point_id + pla.b[segment_id]<<" "<< original_array[int(point_id)] << endl;
//				}*/
//				//cout << "difference: " << difference << endl;
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//			}
//
//			result_collection.max_deviation += deviation_max;
//			result_collection.max_deviation_multiple_width += deviation_max * input_argument.segment_length_second;
//		}
//	}
//	//APCA_KNN_QUAL::printArray(pointer,1280);
//	/*cout << "approximation:" << endl;
//	APCA_KNN_QUAL::printArray(approximation, 1280);*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(count[i] - count[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//
//
//	deviation_sum = sqrt(sum);
//	result_collection.sum_deviation = deviation_sum;
//
//	double test_deviation = TOOL::getDeviation(original_time_series, approximation);
//	assert(float(deviation_sum) == float(test_deviation));
//
//	approximation.clear();
//	approximation.shrink_to_fit();
//	return deviation_sum;
//}
//
//TEMPLATE
//double& PLA_QUAL::getReconstructionErrorPLA(const INPUT_ARGUMENT& input_argument, DataType*& const original_array, const PLA& pla, double& const deviation_sum, double& const deviation_max) {
//	assert(input_argument.segment_length_second != NULL);
//	assert(pla.segmentNum == input_argument.point_dimension);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//
//	//assert(input_argument.time_series_length % input_argument.point_dimension == 0);
//
//	int endOfLongSegment = input_argument.remainder * input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	double difference = NULL;
//	double sum = 0;
//	double point_id = NULL;
//	deviation_max = 0;
//	int count[1280];
//	int i = 0;
//	int* pointer = count;
//
//	DataType* approximation = new DataType[input_argument.time_series_length];
//
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		/*for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//			point_id = input_argument.segment_length_second*segment_id + interval_id;
//			difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//			deviation_max = difference > deviation_max ? difference : deviation_max;
//			sum += difference * difference;
//		}*/
//
//		if (segment_id < input_argument.remainder) {
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
//				//cout << "id: " << point_id <<endl;
//				count[i] = point_id;
//				i++;
//				difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximation[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				//cout << "difference: "<<difference<<endl;
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//			}
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//			assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
//			//assert(indexOfLongSegment== input_argument.remainder-1);
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment) * input_argument.segment_length_second + interval_id;// time series id
//				assert(point_id >= 0);
//				//cout << "id: " << point_id << endl;
//				count[i] = point_id;
//				if (i != 0) {
//					assert(count[i] - count[i - 1] == 1);
//				}
//				i++;
//				difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximation[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				/*if (difference>50) {
//					cout << pla.a[segment_id] * point_id + pla.b[segment_id]<<" "<< original_array[int(point_id)] << endl;
//				}*/
//				//cout << "difference: " << difference << endl;
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//			}
//		}
//	}
//	//APCA_KNN_QUAL::printArray(pointer,1280);
//	/*cout << "approximation:" << endl;
//	APCA_KNN_QUAL::printArray(approximation, 1280);*/
//
//	//for (int i = 1; i < 1280; i++) {
//	//	//cout << test_array[i] << ", ";
//	//	assert(count[i] - count[i - 1] == 1);
//	//	//assert(test_array[i+1]- test_array[i]==1);
//
//	//}
//	//cout << endl;
//
//	deviation_sum = sqrt(sum);
//	delete[] approximation;
//	return deviation_sum;
//}
//
////************************************
//// Method:getReconstructionErrorPLA
//// Qualifier: Get max and sum deviation
//// Input: Orginal time series
//// Output:max and sum deviation with reconstructed PLA
//// date:181212
//// author:
////************************************
//TEMPLATE
//double& PLA_QUAL::getReconstructionErrorPLA(const INPUT_ARGUMENT& const input_argument, double*& const original_array, double& const deviation_sum, double& const deviation_max) {// 181212 deviation_sum deviation_max
//	assert(input_argument.segment_length_second != NULL);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//	//assert(input_argument.time_series_length % input_argument.point_dimension == 0);
//	deviation_sum = 0;
//	deviation_max = 0;
//
//	int endOfLongSegment = input_argument.remainder * input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	double difference = NULL;
//	double sum = 0;
//	double point_id = NULL;
//
//	DataType* reconstructed_time_series = new DataType[input_argument.time_series_length];
//	DataType* segment_deviation = new DataType[input_argument.point_dimension];
//
//	/*=====================================get PLA=====================================*/
//	TOOL::recordStartTime(TOOL::time_record[2]);
//	PLA pla;
//	initialPLA(pla, input_argument.point_dimension);
//	getPLA(input_argument, original_array, pla);
//	output_argument.run_time = TOOL::recordFinishTime(TOOL::time_record[2]);
//	/*..................................................................................*/
//	//cout << "PLA running Time: " << output_argument.run_time << endl;// compare percentage time
//
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		segment_deviation[segment_id] = 0;
//		if (segment_id < input_argument.remainder) {
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
//				//cout << "id: " << point_id <<endl;
//				difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_array[int(point_id)]);
//				reconstructed_time_series[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				//cout << "difference: "<<difference<<endl;
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//				segment_deviation[segment_id] += difference * difference;
//			}
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//			assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment) * input_argument.segment_length_second + interval_id;// time series id
//				assert(point_id >= 0);
//				//cout << "id: " << point_id << endl;
//				difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_array[int(point_id)]);
//				reconstructed_time_series[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//				segment_deviation[segment_id] += difference * difference;
//			}
//		}
//		segment_deviation[segment_id] = sqrt(segment_deviation[segment_id]);
//	}
//
//	//print Array
//	/*cout << "PLA Reconstructed time series:" << endl;
//	TOOL::printArray(reconstructed_time_series, input_argument.time_series_length);*/
//	cout << "PLA  deviation: ";
//	TOOL::printArray(segment_deviation, input_argument.point_dimension);
//	//Write Array
//	TOOL::writeSingleResult("./200706AllAPLAEvaluation/ReconstructPLA181218", reconstructed_time_series, input_argument.time_series_length);
//
//	deviation_sum = sqrt(sum);
//
//	TOOL::getDeviation(input_argument, original_array, reconstructed_time_series, output_argument);
//	assert(deviation_sum == output_argument.sum_deviation && deviation_max == output_argument.max_deviation);
//
//	delete[] reconstructed_time_series;
//	deletePLA(pla);
//
//	TOOL::deleteArray(segment_deviation);
//
//	return deviation_sum;
//}
//
////************************************
//// Method:getReconstructionErrorPLAVector
//// Qualifier: Get max and sum deviation
//// Input: Orginal time series
//// Output:max and sum deviation with reconstructed PLA
//// date:190501
//// author:
////************************************
//TEMPLATE
//double& PLA_QUAL::getReconstructionErrorPLAVector(const INPUT_ARGUMENT& const input_argument, double*& const original_array, double& const deviation_sum, double& const deviation_max) {// 190501 Vector to instead array
//	TOOL::recordStartTime(TOOL::time_record[2]);
//#ifdef _DEBUG
//	assert(input_argument.segment_length_second != NULL);
//#endif
//
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//	//assert(input_argument.time_series_length % input_argument.point_dimension == 0);
//	//DataType* reconstructed_time_series = new DataType[input_argument.time_series_length];
//	vector<DataType> reconstructed_time_series;//190501
//	//DataType* segment_deviation = new DataType[input_argument.point_dimension];
//	vector<DataType> segment_deviation;//190501
//	deviation_sum = 0;
//	deviation_max = 0;
//	int endOfLongSegment = input_argument.remainder * input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	double difference = NULL;
//	double sum = 0;
//	double absolut_difference_sum = 0;
//	double point_id = NULL;
//
//	/*============================================*/
//	double test_parallelogram_height = 0;//190619
//	vector<double> absolute_difference_vector;//190619
//	vector<double> segment_series;
//	/*.................................................*/
//
//	/*cout << "Original array: " << endl;
//	for (int i = 0; i < input_argument.time_series_length;i++) {
//		cout << original_array[i] << ",";
//	}
//	cout << endl;*/
//
//	/*=====================================get PLA=====================================*/
//
//	//vector<PLA> pla;//190501
//	DoublyLinkedList<PLA> pla = DoublyLinkedList<PLA>();//191104 use linked list to insead vector
//	//initialPLA(pla, input_argument.point_dimension);
//	//getPLA(input_argument, original_array, pla);
//
//	//getPLAVector(input_argument, original_array, pla);//190501
//	//getPLAVectorNormal(input_argument, original_array, pla);//190604
//	getPLALinkedListNormal(input_argument, original_array, pla);//191104 use linked list to instead vector
//	//getPLAVectorPreMemory(input_argument, original_array, pla);//190611
//#ifdef _DEBUG
//	assert(pla.size() == input_argument.point_dimension);
//#endif
//	output_argument.run_time = TOOL::recordFinishTime(TOOL::time_record[2]);
//	/*..................................................................................*/
//	//cout << "PLA running Time: " << output_argument.run_time << endl;// compare percentage time
//	output_argument.sum_area = 0;
//	output_argument.sum_density = 0;
//	output_argument.sum_area0 = 0;
//	output_argument.sum_density0 = 0;
//	//cout << "Reconstrcut: " << endl;
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		//segment_deviation[segment_id] = 0;
//		segment_deviation.push_back(0.0);
//
//		/*==================================================================================*/
//		output_argument.sum_area += computeParallelogram(original_array, pla[segment_id]);//190619
//		output_argument.sum_density += 1.0 / pla[segment_id].parallelogram_height;//190619
//		/*....................................................................................*/
//
//		if (segment_id < input_argument.remainder) {
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
//				//cout << "id: " << point_id <<endl;
//				//difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_array[int(point_id)]);
//				difference = fabs(*pla[segment_id].a * interval_id + *pla[segment_id].b - original_array[int(point_id)]);//190501
//				absolute_difference_vector.push_back(*pla[segment_id].a * interval_id + *pla[segment_id].b - original_array[int(point_id)]);//190619
//				segment_series.push_back(original_array[int(point_id)]);
//				//reconstructed_time_series[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				reconstructed_time_series.push_back(*pla[segment_id].a * interval_id + *pla[segment_id].b);//190501
//				//cout << *pla[segment_id].a * interval_id + *pla[segment_id].b << ",";
//				//cout << "difference: "<<difference<<endl;
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//				absolut_difference_sum += difference;
//				segment_deviation[segment_id] += difference * difference;
//			}
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//#ifdef _DEBUG
//			assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
//#endif
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment) * input_argument.segment_length_second + interval_id;// time series id
//#ifdef _DEBUG
//				assert(point_id >= 0);
//#endif
//				//cout << "id: " << point_id << endl;
//				//difference = fabs(pla.a[segment_id] * interval_id + pla.b[segment_id] - original_array[int(point_id)]);
//				difference = fabs(*pla[segment_id].a * interval_id + *pla[segment_id].b - original_array[int(point_id)]);//190501
//				absolute_difference_vector.push_back(*pla[segment_id].a * interval_id + *pla[segment_id].b - original_array[int(point_id)]);//190619
//				segment_series.push_back(original_array[int(point_id)]);
//				//reconstructed_time_series[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				//reconstructed_time_series.push_back(pla.a[segment_id] * interval_id + pla.b[segment_id]);
//				reconstructed_time_series.push_back(*pla[segment_id].a * interval_id + *pla[segment_id].b);//190501
//				//cout << *pla[segment_id].a * interval_id + *pla[segment_id].b << ",";
//				deviation_max = difference > deviation_max ? difference : deviation_max;
//				sum += difference * difference;
//				segment_deviation[segment_id] += difference * difference;
//				absolut_difference_sum += difference;
//			}
//		}
//		auto min_max_difference = minmax_element(absolute_difference_vector.begin(), absolute_difference_vector.end());
//		auto min_max_series = minmax_element(segment_series.begin(), segment_series.end());
//		double min_max_height = *min_max_series.second - *min_max_series.first;
//		double sum_test = *min_max_difference.second - *min_max_difference.first;
//		assert(fabs(float(pla[segment_id].parallelogram_height) - float(sum_test)) < (std::numeric_limits<float>::min)());
//		output_argument.sum_density0 += 1.0 / min_max_height;
//		output_argument.sum_area0 += min_max_height * pla[segment_id].segment_width;
//		absolute_difference_vector.clear();
//		segment_series.clear();
//		segment_deviation[segment_id] = sqrt(segment_deviation[segment_id]);
//	}
//	//cout << endl;
//
//	//print Array
//	/*cout << "PLA Reconstructed time series:" << endl;
//	TOOL::printArray(reconstructed_time_series, input_argument.time_series_length);*/
//	cout << "PLA  deviation: ";
//	TOOL::printArray(segment_deviation, input_argument.point_dimension);
//	//Write Array
//	TOOL::writeSingleResult("./200706AllAPLAEvaluation/ReconstructPLA181218", reconstructed_time_series);
//
//	deviation_sum = sqrt(sum);
//
//	TOOL::getDeviation(input_argument, original_array, reconstructed_time_series, output_argument);
//	assert(deviation_sum == output_argument.sum_deviation && deviation_max == output_argument.max_deviation);
//
//	//output_argument.sum_deviation = absolut_difference_sum;
//	//delete[] reconstructed_time_series;
//	//deletePLA(pla);
//
//	/*=======Delete Memory======*/
//	/*for (auto&& i : pla) {
//		if (i.a != nullptr) {
//			delete i.a;
//			i.a = nullptr;
//		}
//
//		if (i.b != nullptr) {
//			delete i.b;
//			i.b = nullptr;
//		}
//	}*/
//
//	for (int segment_id = 0; segment_id < pla.size(); segment_id++) {
//		if (pla[segment_id].a != nullptr) {
//			delete pla[segment_id].a;
//			pla[segment_id].a = nullptr;
//		}
//
//		if (pla[segment_id].b != nullptr) {
//			delete pla[segment_id].b;
//			pla[segment_id].b = nullptr;
//		}
//	}
//	/*..................*/
//
//	//TOOL::deleteArray(segment_deviation);
//	pla.clear();
//	//pla.shrink_to_fit();
//	reconstructed_time_series.clear();
//	reconstructed_time_series.shrink_to_fit();
//	segment_deviation.clear();
//	segment_deviation.shrink_to_fit();
//
//	return deviation_sum;
//}
//
//TEMPLATE
//void PLA_QUAL::getDeviationIterationPLA(DataType*& const original_time_series, const int& const array_length, string*& const write_file_name, const int& segment_begin, const int& segment_end, const int& segment_interval) {
//	double deviation_sum = 0;
//	double deviation_max = 0;
//
//	DataType* original_normalization = new DataType[array_length];
//	fill_n(original_normalization, array_length, NULL);
//	DataType* PLA_approximation = new DataType[array_length];
//	fill_n(PLA_approximation, array_length, NULL);
//	DataType* PLA_normalization = new DataType[array_length];
//	fill_n(PLA_normalization, array_length, NULL);
//
//	//typename APCA_QUAL::APCA apca;
//	PLA_QUAL::PLA pla_variable;
//
//	APCA_KNN_QUAL::normalizeStandard(original_time_series, array_length, original_normalization);
//
//	//APCA_KNN_QUAL::printArray(original_time_series, array_length);
//	cout << "original_normalization " << endl;
//	APCA_KNN_QUAL::printArray(original_normalization, array_length);
//
//	for (int segment_number = segment_begin; segment_number <= segment_end; segment_number += segment_interval) {
//		cout << "segment_number: " << segment_number << endl;
//		//APCA_QUAL::initialAPCA(apca, segment_number);
//		initialPLA(pla_variable, segment_number);
//
//		PLA_QUAL pla(array_length, segment_number, NULL, NULL, NULL, "ssss");
//
//		//APCA_QUAL::divideRemainderPAA(original_normalization, apca, array_length, segment_number);
//		//APCA_QUAL::getAPCAPoint(original_normalization, array_length, segment_number, apca);
//		pla.getPLA(pla.input_argument, original_normalization, pla_variable);
//		pla.approximateOriginalFunctionPLA1(pla.input_argument, pla_variable, PLA_approximation);
//		/*cout << "approximation:" << endl;
//		APCA_KNN_QUAL::printArray(PLA_approximation, array_length);*/
//
//		APCA_KNN_QUAL::normalizeStandard(PLA_approximation, array_length, PLA_normalization);
//
//		APCA_KNN_QUAL::getReconstructionError(original_normalization, PLA_normalization, array_length, deviation_sum, deviation_max);
//		/*cout << "PLA.a: ";
//		APCA_KNN_QUAL::printArray(pla_variable.a, segment_number);
//		cout << "PLA.b: ";
//		APCA_KNN_QUAL::printArray(pla_variable.b, segment_number);*/
//
//		//getReconstructionErrorPLA(pla.input_argument, original_normalization, pla_variable, deviation_sum, deviation_max);
//
//		//approximateOriginalFunction(input_argument, original_normalization, chebyshev_share, PAA_approximation);
//		//APCA_KNN_QUAL::distanceAE(original_normalization, array_length, apca, deviation_sum, deviation_max);
//		//APCA_KNN_QUAL::normalizeStandard(PAA_approximation, array_length, PAA_normalization);
//		/*APCA_KNN_QUAL::getAverage(chebyshve_normalization, array_length,average);
//		APCA_KNN_QUAL::getVariance(chebyshve_normalization, array_length,variance);*/
//		/*cout << "normalized Chebyshev approximation: ";
//		*/
//		//APCA_KNN_QUAL::getReconstructionError(original_normalization, PAA_normalization, array_length, deviation_sum, deviation_max);
//
//		APCA_KNN_QUAL::writeSingleResult(input_argument, write_file_name[6], deviation_sum);
//		APCA_KNN_QUAL::writeSingleResult(input_argument, write_file_name[7], deviation_max);
//
//		//writeApproximationResult(input_argument, PAA_normalization, deviation_sum, deviation_max);
//
//		//APCA_QUAL::deleteAPCA(apca);
//		deletePLA(pla_variable);
//	}
//
//	delete[] original_normalization;
//	original_normalization = nullptr;
//	delete[] PLA_approximation;
//	PLA_approximation = nullptr;
//	delete[] PLA_normalization;
//	PLA_normalization = nullptr;
//}
//
//TEMPLATE//for build RTree //191107 notes: needs improvement
//typename PLA_QUAL::PLA& PLA_QUAL::getPLAMBR(const PLA& pla, PLA& pla_MBR) {
//#ifdef _DEBUG
//	assert(pla.segmentNum << 1 == pla_MBR.segmentNum);
//#endif
//	for (int i = 0; i < pla_MBR.segmentNum; i++) {
//		if (i & 1) {//odd
//			pla_MBR.a[i] = pla_MBR.b[i] = pla.b[(i - 1) / 2];
//#ifdef _DEBUG
//			assert(pla_MBR.a[i] != NULL);
//#endif
//		}
//		else {//even
//			pla_MBR.a[i] = pla_MBR.b[i] = pla.a[i / 2];
//		}
//		//cout <<"PLA MBR: a: " <<pla_MBR.a[i] <<", b:" << pla_MBR.b[i] << endl;
//	}
//
//	return pla_MBR;
//}
//
//TEMPLATE
//RTREE& PLA_QUAL::buidRTreeIndex(INPUT_ARGUMENT& const input_argument, RTREE& PLARTree, PLA* PLA_array_accumulate, const string& file_name) {
//	//RTree<DataType, ElementType>& buidRTreeIndex(RTree<DataType, ElementType> &APCARTree, double* g_query_time_series, const double& g_time_series_length, const double& g_index_point_number, double(&test_d_original_time_series)[ROW][COLUMN], Link *APCALinkOriginal) {
//	assert(input_argument.build_rtree_time == 0.0);
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//	printf(">>>***Build RTree Index***<<<\n");
//
//	TOOL::recordStartTime(TOOL::time_record[0]);//whole build_rtree_time
//
//	int i = NULL, j = NULL, f_insert_count = NULL;
//	DataType* original_time_series = new DataType[input_argument.time_series_length];
//
//	PLA pla_array_MBR; //Temp MBR
//	assert(input_argument.point_dimension * 2 == input_argument.point_dimension << 1);
//	initialPLA(pla_array_MBR, input_argument.point_dimension << 1);
//
//	/*cout << "Query Point : ";
//	for (i = 0; i < g_time_series_length; i++) cout << g_query_time_series[i] << ", ";
//	cout << endl;*/
//
//	string fs_row_string;//for file stream
//	string fs_row_number;//for every data in file stream
//	ifstream file_stream = ifstream(file_name);
//	assert(file_stream);
//
//	f_insert_count = 0;
//	for (i = 0; i < input_argument.point_number && (!file_stream.eof()) && file_stream.is_open() && file_stream.good(); i++) {
//		//getRandomAPCAPoint(APCALinkOriginal[i].originalLink, g_time_series_length);
//		//APCALinkOriginal[i].originalLink = test_d_original_time_series[i];
//		//getNormalArray(APCALinkOriginal[i].originalLink, g_time_series_length);
//
//		file_stream >> fs_row_string;
//		//memory_account[2] = fs_row_string.size();
//		stringstream sstr(fs_row_string);
//
//		int string_id = -1;
//		while (getline(sstr, fs_row_number, ',') && string_id < input_argument.time_series_length) {
//			if (string_id > -1) {
//				original_time_series[string_id] = stod(fs_row_number);
//				//cout << original_time_series[string_id] << ", ";
//			}
//			string_id++;
//		}
//		//cout << endl;
//
//		/*if (PAA_or_APCA == 0) {
//			getAPCAPoint(original_time_series, g_time_series_length, APCARTree.NUMDIMS / 2, APCALinkOriginal[i].APCALink);
//		}
//		else {
//			divideRemainderPAA(original_time_series, APCALinkOriginal[i].APCALink, g_time_series_length, APCARTree.NUMDIMS / 2);
//		}*/
//
//		TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);//z-score normalization
//
//		getPLA(input_argument, original_time_series, PLA_array_accumulate[i]);
//		getPLAMBR(PLA_array_accumulate[i], pla_array_MBR);
//
//		/*cout << "PLA MBR.a\n";
//		APCA_KNN_QUAL::printArray(pla_array_MBR.a, pla_array_MBR.segmentNum);*/
//
//		/*cout << "PLA MBR.b\n";
//		APCA_KNN_QUAL::printArray(pla_array_MBR.b, pla_array_MBR.segmentNum);*/
//
//		PLARTree.Insert(pla_array_MBR.a, pla_array_MBR.b, i);// a min, b max
//
//		f_insert_count++;
//	}
//
//	file_stream.close();
//	fs_row_string.clear();
//	fs_row_string.shrink_to_fit();
//	fs_row_number.clear();
//	fs_row_number.shrink_to_fit();
//
//	/*cout << "Root Node : sub node number = " << PLARTree.m_root->m_count << " Root level = : " << PLARTree.m_root->m_level << "\n\nBegin to build a RTree:\n";
//	cout << "\n RTree conclusion\n The number of RTree Data Point = : " << PLARTree.Count() << endl;*/
//
//	deletePLA(pla_array_MBR);
//	delete[] original_time_series;
//	original_time_series = nullptr;
//
//	//***** Delete Leaf node memory   ****//
//	RTREE::Iterator it;
//	for (PLARTree.GetFirst(it); !PLARTree.IsNull(it); PLARTree.GetNext(it)) {
//		it.deleteLeafNodeMemory();
//	}
//
//	TOOL::recordFinishTime(TOOL::time_record[0], input_argument.build_rtree_time);
//	//cout << "RTree build time: " << input_argument.build_rtree_time << " us" << endl;
//	TOOL::writeSingleResult(input_argument.write_file_name[2], input_argument.build_rtree_time);
//
//	return PLARTree;
//}
//
//TEMPLATE
//void PLA_QUAL::SimpleBaseKNNSearch(const INPUT_ARGUMENT& input_argument, const DataType* g_query_time_series, const int& K, const string& file_name, priority_queue<ORIGINAL_TIME_SERIES_PAIR, vector<ORIGINAL_TIME_SERIES_PAIR>, priorityDistanceEUC >& q_base_queue) {
//	printf("<///////**  Base KNN Begin  **////////>\n");
//
//	DataType* original_time_series = new DataType[int(input_argument.time_series_length)];
//	//priority_queue<DataType, vector<DataType>, greater<DataType> > q_base_queue;
//
//	ORIGINAL_TIME_SERIES_PAIR id_distance;
//
//	string fs_row_string;
//	string fs_row_number;
//	ifstream file_stream = ifstream(file_name);
//	assert(file_stream);
//
//	int i = 0, j = NULL;
//	while (!file_stream.eof() && i < input_argument.point_number) {
//		file_stream >> fs_row_string;
//		stringstream sstr(fs_row_string);
//		j = -1;
//		while (getline(sstr, fs_row_number, ',') && j < input_argument.time_series_length) {
//			if (j > -1) {
//				original_time_series[j] = stod(fs_row_number);
//			}
//			j++;
//		}
//		id_distance.original_time_series_id = i;
//		id_distance.d_dist = APCA_KNN_QUAL::distanceEUC(g_query_time_series, input_argument.time_series_length, original_time_series, input_argument.time_series_length);
//		q_base_queue.push(id_distance);
//		i++;
//	}
//
//	file_stream.close();
//	delete[] original_time_series;
//	original_time_series = nullptr;
//
//	/*cout << "Top K: \n";
//	for (i = 0; i < K; i++) {
//		cout << q_base_queue.top().d_dist << "\n";
//		q_base_queue.pop();
//	}
//	cout << endl;*/
//}
//
//TEMPLATE
//bool PLA_QUAL::PLAKNNSearch(INPUT_ARGUMENT& const  input_argument, const DataType* g_query_time_series, const RTREE& PLARTree, PLA* PLA_array_accumulate, const int& K, const string& file_name, list<ORIGINAL_TIME_SERIES_PAIR>& result) {
//	printf("<///////** PLA KNN Begin  **////////>\n");
//	//g_n_account_apca_point = 0;//pruning power
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//	assert(input_argument.sum_distance_euc == 0.0);
//	input_argument.navigate_index_time = 0.0;// navigate time
//	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0.0;// distance euclidean time
//
//	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[3]);//for total KNN time
//
//	int i = NULL, j = NULL;
//	assert(K <= input_argument.point_number);
//
//	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[1]);//for approximate query time
//	PLA PLA_query;
//	initialPLA(PLA_query, input_argument.point_dimension);
//	getPLA(input_argument, g_query_time_series, PLA_query);
//	APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[1], input_argument.approximation_query_time);//for approximate query time
//	cout << "PLA approximation query time : " << input_argument.approximation_query_time << " us" << endl;
//
//	cout << "KNN for PLA" << endl;
//	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[2]);//for rest part of KNN time
//
//	priority_queue <PLA_NODE_PAIR, vector<PLA_NODE_PAIR>, priorityIncrement > queue;
//	PLA_NODE_PAIR f_APCA_Root, f_temp_APCA_Pair;
//
//	list<ORIGINAL_TIME_SERIES_PAIR> temp;// , result;
//	ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair;
//
//	f_APCA_Root.p_rtree_node = PLARTree.m_root;
//	f_APCA_Root.d_dist = 0;
//	queue.push(f_APCA_Root);
//	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
//	//printf("<///////**    KNN Begin   **////////>\n");
//
//	int n_data_point_count = 0;
//	int n_distanceLB_count = 0;
//	int n_push_leaf_node_count = 0;
//
//	//g_n_time_getMBR_count = 0;
//	/*g_n_time_apca_point_count = 0;
//	g_n_time_count1 = 0;
//	g_n_time_count2 = 0;
//	g_n_time_count_while_first_part = 0;
//
//	g_n_time_loop_result_count = 0;
//	g_n_time_leaf_node_push_count = 0;
//	g_n_time_leaf_node_distanceLB_count = 0;
//	g_n_time_leaf_node_assignment_count = 0;
//	g_n_time_child_node_MINDIST_count = 0;
//	g_n_time_child_node_push_count = 0;*/
//
//	/*g_d_time_whole_first_run = recordFinishTime(whole_first_time_start[0], whole_first_time_over[0], whole_first_dqFreq[0], whole_first_run_time);
//	g_d_time_whole_first_run = recordFinishTime(time_record[0], whole_first_run_time);*/
//
//	double temp_navigate_time = 0;
//	TOOL::recordStartTime(TOOL::time_record[4]);//navigate index time
//
//	int n_result_count = 0;
//	while (!queue.empty()) {
//		/*while (!queue.empty()|| result.size() != K) {*/
//		//while (result.size() != K) {
//		//cout << "KNN: the " << count << " turn : \n";
//		double run_time3;
//		//recordStartTime(time_record[1]);
//		/*if (count == K) {
//		cout << "////*********EMPTY  Found " << K << " result!!!!!!********///" << endl;break;}*/
//		double run_time5;
//		//recordStartTime(time_record[2]);
//		PLA_NODE_PAIR m_temp_queue_top = queue.top();
//
//		/*if (queue.size() > memory_account[9]) {
//			memory_account[9] = queue.size();
//		}*/
//
//		/*cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator.dist: ";
//		for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
//		cout << endl;*/
//
//		//if (memory_account[10] < temp.size()) { memory_account[10] = temp.size(); }
//
//		for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
//			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
//			if (plist->d_dist <= m_temp_queue_top.d_dist) {
//				//cout << "           <= " << endl;
//				result.push_back(*plist);
//				plist = temp.erase(plist);
//			}
//			else plist++;
//			if (K == result.size()) {
//				//rest part of KNN
//				APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[2], input_argument.knn_rest_part_time);
//				cout << "rest part of KNN time : " << input_argument.knn_rest_part_time << " us" << endl;
//
//				//time of whole KNN procedure
//				APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[3], input_argument.knn_total_time);
//				//cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
//				input_argument.pruning_power = input_argument.sum_distance_euc / double(input_argument.point_number);
//				//cout << "pruning power: " << input_argument.pruning_power << endl;
//
//				input_argument.IO_cost = input_argument.sum_distance_euc;
//				//cout << "I/O cost: " << input_argument.IO_cost << endl;
//
//				TOOL::writeSingleResult(input_argument.write_file_name[0], input_argument.pruning_power);
//				TOOL::writeSingleResult(input_argument.write_file_name[1], input_argument.IO_cost);
//				TOOL::writeSingleResult(input_argument.write_file_name[3], input_argument.navigate_index_time);
//				TOOL::writeSingleResult(input_argument.write_file_name[4], input_argument.distance_lowbound_time);
//				TOOL::writeSingleResult(input_argument.write_file_name[5], input_argument.distance_euc_time);
//				TOOL::writeSingleResult(input_argument.write_file_name[6], input_argument.knn_total_time);
//
//				cout << "*************************Find result!!!!!!!!!!!!! result.size: " << result.size() << endl;
//				result.sort([](const ORIGINAL_TIME_SERIES_PAIR& first, const  ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist; });
//
//				/*cout << "Result list: \n";
//				for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//					cout << "id:" << it->original_time_series_id << ", dist: " << it->d_dist << "\n";
//				}*/
//
//				deletePLA(PLA_query);
//				priority_queue<PLA_NODE_PAIR, vector<PLA_NODE_PAIR>, priorityIncrement >().swap(queue);
//				temp.clear();
//				//result.clear();
//				list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//				//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//
//				return true;
//			}
//		}
//
//		//g_n_time_loop_result_count += recordFinishTime(time_record[2], run_time5);
//
//		//printf("\nrun_time = %f us\n", run_time0);
//		//cout << "Before Pop: Queue.size(): " << queue.size() << " Pop: top.dist:" << queue.top().d_dist << " " << queue.top().original_time_series_id << endl;
//		queue.pop();
//		//cout << "After Pop: Queue.size(): " << queue.size();
//		/*if (!queue.empty()) cout << " Pop: top.dist:" << queue.top().d_dist << " " << queue.top().original_time_series_id << endl;
//		else cout << endl;*/
//
//		//g_n_time_count_while_first_part = recordFinishTime(time_record[1], run_time3);
//		//printf("\nrun_time = %f us\n", run_time0);
//
//		if (m_temp_queue_top.p_rtree_node == nullptr) { //is PLA data point
//			TOOL::recordFinishTime(TOOL::time_record[4], temp_navigate_time);
//			input_argument.navigate_index_time += temp_navigate_time;
//
//			//cout << "    queue.top is data point\n";
//			DataType* original_time_series = new DataType[input_argument.time_series_length];
//
//			//g_n_account_apca_point++;
//			double run_time0;
//
//			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
//			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
//			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;
//			APCA_KNN_QUAL::getFileStreamByID(file_name, input_argument.time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//
//			double distance_euclidean_time = 0;
//			TOOL::recordStartTime(TOOL::time_record[6]);//time for euclidean distance
//
//			tempOriginalTimeSeriesPair.d_dist = typename APCA_KNN_QUAL::distanceEUC(g_query_time_series, input_argument.time_series_length, original_time_series, input_argument.time_series_length);
//
//			TOOL::recordFinishTime(TOOL::time_record[6], distance_euclidean_time);
//			input_argument.distance_euc_time += distance_euclidean_time;
//			input_argument.sum_distance_euc++;
//
//			delete[] original_time_series;
//			original_time_series = nullptr;
//
//			//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
//			temp.push_back(tempOriginalTimeSeriesPair);
//
//			//printf("\nrun_time = %f us\n", run_time0);
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
//			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			typename APCA_QUAL::APCA QProjection;
//			typename APCA_QUAL::initialAPCA(QProjection, PLARTree.NUMDIMS / 2);
//
//			double run_time1;
//			recordStartTime(time_record[4]);
//
//			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {
//				//g_n_account_leaf_node++;
//
//				f_temp_APCA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);
//				//f_temp_APCA_Pair.p_APCA_point = &APCALinkOriginal[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)].APCALink;//
//				f_temp_APCA_Pair.p_PLA_point = &PLA_array_accumulate[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)];
//
//				//printf("\nrun_time = %f us\n", run_time0);
//				double distance_PLA_time = 0;
//				TOOL::recordStartTime(TOOL::time_record[5]);//distance PLA & MBR time
//
//				//f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point), *f_temp_APCA_Pair.p_APCA_point);
//				//f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point, QProjection), *f_temp_APCA_Pair.p_APCA_point);
//				f_temp_APCA_Pair.d_dist = getPLADistance(input_argument, *f_temp_APCA_Pair.p_PLA_point, PLA_query, f_temp_APCA_Pair.d_dist);
//
//				TOOL::recordFinishTime(TOOL::time_record[5], distance_PLA_time);
//				input_argument.distance_lowbound_time += distance_PLA_time;
//
//				f_temp_APCA_Pair.p_rtree_node = nullptr; //attach data point label
//
//				double queue_run_time;
//				recordStartTime(time_record[7]);
//
//				//cout << "            Push apca data. Branch id: " << i << " DLB: " << f_temp_APCA_Pair.d_dist << ", time series ID: " << f_temp_APCA_Pair.original_time_series_id << endl;
//				queue.push(f_temp_APCA_Pair);
//
//				n_push_leaf_node_count++;
//				//g_n_time_leaf_node_push_count += recordFinishTime(time_record[7], queue_run_time);
//
//				//printf("\nrun_time = %f us\n", run_time0);
//				//cout << "KNN : queue.size() = " << queue.size() << endl;
//			}
//			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
//			//cout << tempOriginalTimeSeriesPair.key;
//			//g_n_time_count1 += recordFinishTime(time_record[4], run_time1);
//			//printf("\nrun_time = %f us\n", run_time0);
//			APCA_QUAL::deleteAPCA(QProjection);
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {															//is internal Node
//			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			double run_time2;
//			recordStartTime(time_record[8]);
//
//			PLA_NODE_PAIR tempApcaPair;
//			typename APCA_KNN_QUAL::REGION fs_region_G;
//			APCA_KNN_QUAL::initialREGION(fs_region_G, PLARTree.NUMDIMS / 2);
//			//cout << "regionNum = " << G.regionNum << endl;
//
//			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
//				//g_n_account_child_node++;
//
//				double mindistQR_run_time;
//				recordStartTime(time_record[9]);
//
//				double distance_PLA_time = 0;
//				TOOL::recordStartTime(TOOL::time_record[5]);//distance PLA & MBR time
//
//				//tempApcaPair.d_dist = MINDISTQR(g_query_time_series, g_time_series_length, getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//				tempApcaPair.d_dist = getPLAMBRDistance(input_argument, m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, PLA_query, tempApcaPair.d_dist);
//				/*for (int i = 0; i < fs_region_G.regionNum; i++) {
//				cout << "            G1[" << i << "] = " << fs_region_G.G1[i] << ", G2[" << i << "] = " << fs_region_G.G2[i] << ", G3[" << i << "] = " << fs_region_G.G3[i] << ", G4[" << i << "] = " << fs_region_G.G4[i] << endl;
//				}*/
//
//				TOOL::recordFinishTime(TOOL::time_record[5], distance_PLA_time);
//				input_argument.distance_lowbound_time += distance_PLA_time;
//
//				//g_n_time_child_node_MINDIST_count += recordFinishTime(time_record[9], mindistQR_run_time);
//				//printf("\nrun_time = %f us\n", run_time0);
//
//				tempApcaPair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//
//				double push_run_time;
//				recordStartTime(time_record[10]);
//				//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//				queue.push(tempApcaPair);
//				//cout << "KNN : queue.size() = " << queue.size() << endl;
//				//g_n_time_child_node_push_count += recordFinishTime(time_record[10], push_run_time);
//				//printf("\nrun_time = %f us\n", run_time0);
//			}
//
//			APCA_KNN_QUAL::deleteREGION(fs_region_G);
//			//g_n_time_count2 += recordFinishTime(time_record[8], run_time2);
//			//printf("\nrun_time = %f us\n", run_time0);
//		}
//		else {
//			cout << "WRONG!!!!!!!  WRONG!!!!!!!!!!" << endl;
//			assert(0);
//		}
//	}
//
//	cout << "??????????????????     PLA KNN failed    ??????????????????" << endl;
//	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	ofstream outfile(input_argument.write_file_name[8] + ".txt", ios::app);
//	assert(outfile.is_open());
//	outfile << "??????????????????     PLA KNN failed    ??????????????????" << endl;
//	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
//	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	outfile.close();
//	//assert(0);
//
//	deletePLA(PLA_query);
//	priority_queue<PLA_NODE_PAIR, vector<PLA_NODE_PAIR>, priorityIncrement >().swap(queue);
//	temp.clear();
//	list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//	//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//	return false;
//}
//
//TEMPLATE
//bool PLA_QUAL::PLAKNNMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const DataType* g_query_time_series, const RTREE& PLARTree, PLA* PLA_array_accumulate, const int& K) {
//	int n_push_leaf_node_count = 0;
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//#ifdef _DEBUG
//	printf("<///////** PLA KNN Begin  **////////>\n");
//	//g_n_account_apca_point = 0;//pruning power
//
//	assert(input_argument.sum_distance_euc == 0.0);
//
//	int n_data_point_count = 0;
//	int n_distanceLB_count = 0;
//
//
//	/*g_n_time_getMBR_count = 0;
//	g_n_time_apca_point_count = 0;
//	g_n_time_count1 = 0;
//	g_n_time_count2 = 0;
//	g_n_time_count_while_first_part = 0;
//
//	g_n_time_loop_result_count = 0;
//	g_n_time_leaf_node_push_count = 0;
//	g_n_time_leaf_node_distanceLB_count = 0;
//	g_n_time_leaf_node_assignment_count = 0;
//	g_n_time_child_node_MINDIST_count = 0;
//	g_n_time_child_node_push_count = 0;*/
//
//	input_argument.knn_total_time = 0.0;
//	input_argument.navigate_index_time = 0.0;// navigate time
//	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0.0;// distance euclidean time
//	//APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[3]);//for total KNN time
//
//	assert(K <= input_argument.point_number);
//#endif
//
//	/*----------------------------------Evaluation: Mathlab Bar Chart---------------------------------------------*/
//	TOOL::recordStartTime(TOOL::time_record[14]);
//	input_argument.IO_cost = 0;// measure I/O cost
//	input_argument.sum_distance_euc = 0.0;
//	/*------------------------------------------------------------------------------------------------------------*/
//	int i = NULL, j = NULL;
//	//APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[1]);//for approximate query time
//
//	PLA PLA_query;
//#ifdef _DEBUG
//	assert(input_argument.time_series_length * input_argument.arity_d == input_argument.point_multi_single_length);
//	assert(input_argument.point_dimension * input_argument.arity_d == input_argument.point_multi_single_dimension);
//#endif
//	initialPLA(PLA_query, input_argument.point_multi_single_dimension);//????????180918 , this multi ot single has problem
//	getPLA(input_argument.point_multi_single_length, input_argument.point_multi_single_dimension, g_query_time_series, PLA_query);
//
//	//APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[1], input_argument.approximation_query_time);//for approximate query time
//	//cout << "PLA approximation query time : " << input_argument.approximation_query_time << " us" << endl;
//#ifdef _DEBUG
//	cout << "KNN for PLA" << endl;
//#endif
//	//APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[2]);//for rest part of KNN time
//
//	priority_queue <PLA_NODE_PAIR, vector<PLA_NODE_PAIR>, priorityIncrement > queue;
//	PLA_NODE_PAIR f_APCA_Root, f_temp_APCA_Pair;
//
//	list<ORIGINAL_TIME_SERIES_PAIR> temp;
//	list<ORIGINAL_TIME_SERIES_PAIR> result;
//	ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair;
//
//	f_APCA_Root.p_rtree_node = PLARTree.m_root;
//	f_APCA_Root.d_dist = 0;
//	queue.push(f_APCA_Root);
//	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
//	//printf("<///////**    KNN Begin   **////////>\n");
//
//	//g_d_time_whole_first_run = recordFinishTime(whole_first_time_start[0], whole_first_time_over[0], whole_first_dqFreq[0], whole_first_run_time);
//#ifdef _DEBUG
//	//g_d_time_whole_first_run = recordFinishTime(time_record[0], whole_first_run_time);
//
//	double temp_navigate_time = 0;
//	TOOL::recordStartTime(TOOL::time_record[4]);//navigate index time
//
//	int n_result_count = 0;
//#endif
//	while (!queue.empty()) {
//		/*while (!queue.empty()|| result.size() != K) {*/
//		//while (result.size() != K) {
//		//cout << "KNN: the " << count << " turn : \n";
//#ifdef _DEBUG
//		double run_time3;
//		//recordStartTime(time_record[1]);
//		/*if (count == K) {
//		cout << "////*********EMPTY  Found " << K << " result!!!!!!********///" << endl;break;}*/
//		double run_time5;
//#endif
//		//recordStartTime(time_record[2]);
//		PLA_NODE_PAIR m_temp_queue_top = queue.top();
////#ifdef _DEBUG
////		if (queue.size() > memory_account[9])
////			memory_account[9] = queue.size();
////#endif
//		/*cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator.dist: ";
//		for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
//		cout << endl;*/
////#ifdef _DEBUG
////		if (memory_account[10] < temp.size()) { memory_account[10] = temp.size(); }
////#endif
//
//		for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
//			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
//			if (plist->d_dist <= m_temp_queue_top.d_dist) {
//				//cout << "           <= " << endl;
//				result.push_back(*plist);
//				plist = temp.erase(plist);
//			}
//			else plist++;
//			if (K == result.size()) {
//				input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//				//rest part of KNN
//				//APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[2], input_argument.knn_rest_part_time);
//				//cout << "rest part of KNN time : " << input_argument.knn_rest_part_time << " us" << endl;
//
//				//time of whole KNN procedure
//				//APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[3], input_argument.knn_total_time);
//				//cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
//				input_argument.pruning_power = input_argument.sum_distance_euc / double(input_argument.point_number);
//				//cout << "pruning power: " << input_argument.pruning_power << endl;
//#ifdef _DEBUG
//				input_argument.IO_cost = input_argument.sum_distance_euc;
//				//cout << "I/O cost: " << input_argument.IO_cost << endl;
//
//				cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
//
//				cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;
//
//				cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;
//
//				cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;
//				assert(input_argument.pruning_power != INF);
//
//				cout << "pruning power: " << input_argument.pruning_power << endl;
//
//				assert(input_argument.IO_cost == input_argument.sum_distance_euc);
//				cout << "I/O cost: " << input_argument.IO_cost << endl;
//
//				/*TOOL::writeSingleResult(input_argument,input_argument.write_file_name[0], input_argument.pruning_power);
//				TOOL::writeSingleResult(input_argument,input_argument.write_file_name[1], input_argument.IO_cost);
//				TOOL::writeSingleResult(input_argument,input_argument.write_file_name[3], input_argument.navigate_index_time);
//				TOOL::writeSingleResult(input_argument,input_argument.write_file_name[4], input_argument.distance_lowbound_time);
//				TOOL::writeSingleResult(input_argument,input_argument.write_file_name[5], input_argument.distance_euc_time);
//				TOOL::writeSingleResult(input_argument,input_argument.write_file_name[6], input_argument.knn_total_time);*/
//
//				cout << "!!!!!!!!!!!!!!!!!!!!!!!!            PLA Find result     !!!!!!!!!!!!!      result list size: " << result.size() << endl;
//				result.sort([](const ORIGINAL_TIME_SERIES_PAIR& first, const  ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist; });
//
//				/*cout << "Result list: \n";
//				for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//					cout << "id:" << it->original_time_series_id << ", dist: " << it->d_dist << "\n";
//				}*/
//#endif
//				deletePLA(PLA_query);
//				priority_queue<PLA_NODE_PAIR, vector<PLA_NODE_PAIR>, priorityIncrement >().swap(queue);
//				temp.clear();
//				//result.clear();
//				list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//				//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//
//				return true;
//			}
//		}
//
//		//g_n_time_loop_result_count += recordFinishTime(time_record[2], run_time5);
//
//		//printf("\nrun_time = %f us\n", run_time0);
//		//cout << "Before Pop: Queue.size(): " << queue.size() << " Pop: top.dist:" << queue.top().d_dist << " " << queue.top().original_time_series_id << endl;
//		queue.pop();
//		//cout << "After Pop: Queue.size(): " << queue.size();
//		/*if (!queue.empty()) cout << " Pop: top.dist:" << queue.top().d_dist << " " << queue.top().original_time_series_id << endl;
//		else cout << endl;*/
//
//		//g_n_time_count_while_first_part = recordFinishTime(time_record[1], run_time3);
//		//printf("\nrun_time = %f us\n", run_time0);
//
//		if (m_temp_queue_top.p_rtree_node == nullptr) { //is PLA data point
//#ifdef _DEBUG
//			TOOL::recordFinishTime(TOOL::time_record[4], temp_navigate_time);
//			input_argument.navigate_index_time += temp_navigate_time;
//#endif
//			//cout << "    queue.top is data point\n";
//			DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
//#ifdef _DEBUG
//			//g_n_account_apca_point++;
//			double run_time0;
//#endif
//			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
//			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
//			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;
//			//APCA_KNN_QUAL::getFileStreamByID(file_name, input_argument.time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			//TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//
//			if (input_argument.read_multiple_file_name) {
//				TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			}
//			else {
//				//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id));
//				TOOL::getFileStreamByID(TOOL::getStringByID(TOOL::file_address, input_argument.file_id), input_argument.time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			}
//			TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);
//			double distance_euclidean_time = 0;
//			TOOL::recordStartTime(TOOL::time_record[6]);//time for euclidean distance
//			tempOriginalTimeSeriesPair.d_dist = typename APCA_KNN_QUAL::distanceEUC(g_query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);
//
//			TOOL::recordFinishTime(TOOL::time_record[6], distance_euclidean_time);
//			input_argument.distance_euc_time += distance_euclidean_time;
//			input_argument.sum_distance_euc++;
//			input_argument.IO_cost++;
//
//			delete[] original_time_series;
//			original_time_series = nullptr;
//
//			//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
//			temp.push_back(tempOriginalTimeSeriesPair);
//
//			//printf("\nrun_time = %f us\n", run_time0);
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
//			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			//typename APCA_QUAL::APCA QProjection;
//			assert(PLARTree.NUMDIMS / 2 == input_argument.point_multi_single_dimension);
//			//typename APCA_QUAL::initialAPCA(QProjection, PLARTree.NUMDIMS / 2);
//
//			double run_time1;
//			recordStartTime(time_record[4]);
//
//			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {
//				//g_n_account_leaf_node++;
//
//				f_temp_APCA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);
//				//f_temp_APCA_Pair.p_APCA_point = &APCALinkOriginal[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)].APCALink;//
//				f_temp_APCA_Pair.p_PLA_point = &PLA_array_accumulate[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)];
//
//				//printf("\nrun_time = %f us\n", run_time0);
//				double distance_PLA_time = 0;
//				TOOL::recordStartTime(TOOL::time_record[5]);//distance PLA & MBR time
//
//				//f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point), *f_temp_APCA_Pair.p_APCA_point);
//				//f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point, QProjection), *f_temp_APCA_Pair.p_APCA_point);
//				//f_temp_APCA_Pair.d_dist = getPLADistance(input_argument, *f_temp_APCA_Pair.p_PLA_point, PLA_query, f_temp_APCA_Pair.d_dist);
//				f_temp_APCA_Pair.d_dist = getPLADistance(input_argument.point_multi_single_length, *f_temp_APCA_Pair.p_PLA_point, PLA_query, f_temp_APCA_Pair.d_dist);
//
//				TOOL::recordFinishTime(TOOL::time_record[5], distance_PLA_time);
//				input_argument.distance_lowbound_time += distance_PLA_time;
//
//				f_temp_APCA_Pair.p_rtree_node = nullptr; //attach data point label
//
//				double queue_run_time;
//				recordStartTime(time_record[7]);
//
//				//cout << "            Push apca data. Branch id: " << i << " DLB: " << f_temp_APCA_Pair.d_dist << ", time series ID: " << f_temp_APCA_Pair.original_time_series_id << endl;
//				queue.push(f_temp_APCA_Pair);
//
//				n_push_leaf_node_count++;
//				//g_n_time_leaf_node_push_count += recordFinishTime(time_record[7], queue_run_time);
//
//				//printf("\nrun_time = %f us\n", run_time0);
//				//cout << "KNN : queue.size() = " << queue.size() << endl;
//			}
//			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
//			//cout << tempOriginalTimeSeriesPair.key;
//
//			//g_n_time_count1 += recordFinishTime(time_record[4], run_time1);
//			//printf("\nrun_time = %f us\n", run_time0);
//			//APCA_QUAL::deleteAPCA(QProjection);
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {															//is internal Node
//			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			double run_time2;
//			recordStartTime(time_record[8]);
//
//			PLA_NODE_PAIR tempApcaPair;
//			//typename APCA_KNN_QUAL::REGION fs_region_G;
//			//APCA_KNN_QUAL::initialREGION(fs_region_G, PLARTree.NUMDIMS / 2);
//			//cout << "regionNum = " << G.regionNum << endl;
//
//			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
//				//g_n_account_child_node++;
//
//				double mindistQR_run_time;
//				recordStartTime(time_record[9]);
//
//				double distance_PLA_time = 0;
//				TOOL::recordStartTime(TOOL::time_record[5]);//distance PLA & MBR time
//
//				//tempApcaPair.d_dist = MINDISTQR(g_query_time_series, g_time_series_length, getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//				tempApcaPair.d_dist = getPLAMBRDistance(input_argument, m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, PLA_query, tempApcaPair.d_dist);
//				/*for (int i = 0; i < fs_region_G.regionNum; i++) {
//				cout << "            G1[" << i << "] = " << fs_region_G.G1[i] << ", G2[" << i << "] = " << fs_region_G.G2[i] << ", G3[" << i << "] = " << fs_region_G.G3[i] << ", G4[" << i << "] = " << fs_region_G.G4[i] << endl;
//				}*/
//
//				TOOL::recordFinishTime(TOOL::time_record[5], distance_PLA_time);
//				input_argument.distance_lowbound_time += distance_PLA_time;
//
//				//g_n_time_child_node_MINDIST_count += recordFinishTime(time_record[9], mindistQR_run_time);
//				//printf("\nrun_time = %f us\n", run_time0);
//
//				tempApcaPair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//
//				double push_run_time;
//				recordStartTime(time_record[10]);
//				//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//				queue.push(tempApcaPair);
//				//cout << "KNN : queue.size() = " << queue.size() << endl;
//				//g_n_time_child_node_push_count += recordFinishTime(time_record[10], push_run_time);
//				//printf("\nrun_time = %f us\n", run_time0);
//			}
//
//			//APCA_KNN_QUAL::deleteREGION(fs_region_G);
//			//g_n_time_count2 += recordFinishTime(time_record[8], run_time2);
//			//printf("\nrun_time = %f us\n", run_time0);
//		}
//		else {
//			cout << "WRONG!!!!!!!  WRONG!!!!!!!!!!" << endl;
//			assert(0);
//		}
//	}
//	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//	input_argument.pruning_power = 2;
//
//#ifdef _DEBUG
//	cout << "??????????????????     PLA KNN failed    ??????????????????" << endl;
//	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	assert(input_argument.pruning_power != INF);
//#endif
//
//	/*ofstream outfile(input_argument.write_file_name[8] + ".txt", ios::app);
//	assert(outfile.is_open());
//	outfile << "??????????????????     PLA KNN failed    ??????????????????" << endl;
//	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
//	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	outfile.close();*/
//	//assert(0);
//
//	deletePLA(PLA_query);
//	priority_queue<PLA_NODE_PAIR, vector<PLA_NODE_PAIR>, priorityIncrement >().swap(queue);
//	temp.clear();
//	list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//	//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//	return false;
//}
//
////************************************
//// Method:computeParallelogram
//// Qualifier: Get PLA a&b, parallelogram area and deviation point of Parallelogram
//// date:190619
//// author:
////************************************
//TEMPLATE
//double& PLA_QUAL::getAAndBByPLA(DataType*& const original_time_series, PLA& const temp_coefficient) {//190619 return a;
//#ifdef _DEBUG
//	assert(temp_coefficient.segment_width > 0 && temp_coefficient.segment_width != INF && temp_coefficient.right_endpoint > 0 && temp_coefficient.right_endpoint != INF);
//#endif
//	double temp_a;
//	double temp_b;
//
//	if (temp_coefficient.segment_width == 1) {
//		assert(0);
//		//return getAandB(original_time_series, temp_coefficient);
//	}
//
//	//printf("getPLA()\n");
//	//assert(temp_coefficient.size() == input_argument.point_dimension);
//
//	int segment_left_id = temp_coefficient.right_endpoint - temp_coefficient.segment_width + 1;
//
//	if (temp_coefficient.segment_width == 2) {
//#ifdef _DEBUG
//		assert(temp_coefficient.right_endpoint - segment_left_id == 1);
//#endif
//		//*temp_coefficient.a = (original_time_series[int(temp_coefficient.right_endpoint)] - original_time_series[segment_left_id]);
//		//*temp_coefficient.b = original_time_series[segment_left_id];
//		temp_a = (original_time_series[int(temp_coefficient.right_endpoint)] - original_time_series[segment_left_id]);
//		temp_b = original_time_series[segment_left_id];
//		assert(*temp_coefficient.a == temp_a && *temp_coefficient.b == temp_b);
//		return *temp_coefficient.a;
//	}
//
//	//coefficient of equation
//	//a
//	double a_minuend = NULL; //(l-1)/2
//	double a_divisor = NULL; //l(l-1)(l+1)
//	//b
//	double b_minuend = NULL; //2l-1
//	double b_divisor = NULL; //l(l+1)
//
//	auto a_sum = 0.0;
//	auto b_sum = 0.0;
//	double variable_id = NULL; //[0-segment_length)
//	int array_id = NULL;//time series id
//
//	a_minuend = (temp_coefficient.segment_width - 1) / 2.0;//(l-1)/2
//	a_divisor = (temp_coefficient.segment_width - 1) * (temp_coefficient.segment_width + 1) * temp_coefficient.segment_width;//l(l-1)(l+1)
//
//	b_minuend = 2.0 * temp_coefficient.segment_width - 1;//2l-1
//	b_divisor = (temp_coefficient.segment_width + 1) * temp_coefficient.segment_width;//l(l+1)
//
//#ifdef _DEBUG
//	assert(a_divisor != 0);
//	assert(b_divisor != 0);
//#endif
//
//	for (variable_id = 0, array_id = segment_left_id; array_id <= temp_coefficient.right_endpoint; array_id++, variable_id++) {
//		a_sum += (variable_id - a_minuend) * original_time_series[array_id];
//		b_sum += (b_minuend - variable_id * 3.0) * original_time_series[array_id];
//	}
//	temp_a = 12.0 * a_sum / a_divisor;
//	temp_b = 2.0 * b_sum / b_divisor;
//
//#ifdef _DEBUG
//	if (temp_coefficient.segment_width == 2) {
//		assert(float(*temp_coefficient.a) == float(temp_a));
//		assert(float(*temp_coefficient.b) == float(temp_b));
//	}
//#endif
//
//	return *temp_coefficient.a;
//}
//
////************************************
//// Method:computeParallelogram
//// Qualifier: Get PLA a&b, parallelogram area and deviation point of Parallelogram
//// date:190619
//// author:
////************************************
//TEMPLATE
//double& PLA_QUAL::computeParallelogram(DataType*& const original_time_series, PLA& const  temp_coefficient) {//190619
//#ifdef _DEBUG
//	assert(temp_coefficient.segment_width != NULL);
//#endif
//
//	double b = NULL;//y=ax+b;
//	int segment_left_id = temp_coefficient.right_endpoint - temp_coefficient.segment_width + 1;
//	auto right_endpoint_id = temp_coefficient.right_endpoint;
//
//	if (right_endpoint_id - segment_left_id == 1) {
//		assert(0);
//		/*temp_coefficient.parallelogram_height = 0;
//		temp_coefficient.deviation_point.id = temp_coefficient.max_point.id;
//		return temp_coefficient.rectangle_area = temp_coefficient.rectangle_height * temp_coefficient.segment_width;*/
//	}
//
//	deque<DataType> left_y_array;
//	deque<DataType> test_deviation_array;
//	//deque<DataType> left_y_array0;
//	//deque<DataType> right_y_array;
//
//	//GEOMETRY::POINT  up_deviation_point;
//	//GEOMETRY::POINT  down_deviation_point;
//
//	/*cout << "Original time series: ";
//	for (int i = segment_left_id; i <= right_endpoint_id; i++) {
//		cout <<i<<" "<< original_time_series[i] << "; ";
//	}
//	cout << endl;*/
//
//	//a = (temp_coefficient.max_point.value-temp_coefficient.min_point.value)/(temp_coefficient.max_point.id - temp_coefficient.min_point.id);
//	auto a = getAAndBByPLA(original_time_series, temp_coefficient);//compute PLA a
//
//	for (int array_id = segment_left_id, zero_id = 0; array_id <= right_endpoint_id; array_id++, zero_id++) {
//		//b = original_time_series[array_id] - a * array_id;
//		//cout << original_time_series[array_id] << ","<<endl;
//		double left_y0 = a * (0 - zero_id) + original_time_series[array_id];
//		//cout << left_y0 << endl;
//		//cout << a* zero_id + *temp_coefficient.b << endl;
//		double point_deviation = a * zero_id + *temp_coefficient.b - original_time_series[array_id];
//
//#ifdef _DEBUG
//		double left_y = a * (segment_left_id - array_id) + original_time_series[array_id];//y_2=a*(x_2-x_1)+y_1
//		assert(left_y == left_y0);
//		assert(0 - zero_id == segment_left_id - array_id && left_y0 - original_time_series[array_id] == left_y - original_time_series[array_id]);
//#endif
//		//left_y_array0.push_back(left_y0);
//		//cout <<"left_y: "<< left_y0 << endl;
//		left_y_array.push_back(left_y0);
//		test_deviation_array.push_back(point_deviation);
//	}
//	//cout << endl;
//	auto left_projection_point = minmax_element(left_y_array.begin(), left_y_array.end());
//	auto test_deviaion_point = minmax_element(test_deviation_array.begin(), test_deviation_array.end());
//
//	//**get parallelogram height & base and area.
//	//sort(left_y_array.begin(), left_y_array.end());
//	//auto parallelogram_base_value0 = minmax_element(left_y_array0.begin(), left_y_array0.end());
//
//	temp_coefficient.parallelogram_height = *test_deviaion_point.second - *test_deviaion_point.first;
//	double test_parallelogram_height = *left_projection_point.second - *left_projection_point.first;
//#ifdef _DEBUG
//	//assert(fabs(float(test_parallelogram_height) - float(temp_coefficient.parallelogram_height) )< (std::numeric_limits<float>::min)());
//	//assert(*test_deviaion_point.second >= 0 && *test_deviaion_point.first <= 0);
//#endif
//	//cout << "!!!!!!!: " << *left_projection_point.second << " " << *left_projection_point.first << " " << original_time_series[segment_left_id] << endl;
//
//	//assert(temp_coefficient.parallelogram_height!=NULL);
//
//	/*============================================================================================================================*/
//	//get split point For some rectangle, we have to split it to two small rectangle to reduce deviation.
//	//up_deviation_point.id = left_projection_point.second - left_y_array.begin() + segment_left_id;
//	//up_deviation_point.value = fabs(*left_projection_point.second - original_time_series[int(segment_left_id)]);//height
//	//up_deviation_point.value = fabs(*left_projection_point.second - temp_coefficient.apla.b);//height
//	//down_deviation_point.id = left_projection_point.first - left_y_array.begin() + segment_left_id;
//	//down_deviation_point.value = fabs(*left_projection_point.first - original_time_series[int(segment_left_id)]);//height
//	//down_deviation_point.value = fabs(*left_projection_point.first - temp_coefficient.apla.b);//height
//	//temp_coefficient.deviation_point.id = down_deviation_point.value > up_deviation_point.value ? down_deviation_point.id : up_deviation_point.id;
//	//auto temp_left_deviation = max(up_deviation_point.value, down_deviation_point.value);
//	//assert(temp_left_deviation == *test_deviaion_point.second);
//	//cout << "ORI parallelogram_height:" << temp_coefficient.parallelogram_height<<"deviation point id: "<< temp_coefficient.deviation_point.id <<" down_deviation_point value: "<< down_deviation_point.value <<"up_deviation_point value"<< up_deviation_point.value << endl;
//	//cout << "down_deviation_point.id: " << down_deviation_point.id << ";    up_deviation_point.id: " << up_deviation_point.id << ". down_deviation_point.value: " << down_deviation_point.value << ";    up_deviation_point.value: " << up_deviation_point.value << endl;
//	/*.................................................................................................................................*/
//
//	//assert(*parallelogram_base_value.second - *parallelogram_base_value.first == *parallelogram_base_value0.second - *parallelogram_base_value0.first);
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	//for (int array_id = left_endpoint_id; array_id <= right_endpoint_id; array_id++) {
//	//	//b = original_time_series[array_id] - a * array_id;
//	//	double right_y = a * (right_endpoint_id - array_id) + original_time_series[array_id];//y_2=a*(x_2-x_1)+y_1
//	//	right_y_array.push_back(right_y);
//	//}
//	//sort(left_y_array.begin(), left_y_array.end());
//	//auto parallelogram_base_value0 = minmax_element(right_y_array.begin(), right_y_array.end());
//	//cout << *parallelogram_base_value.second << " " << *parallelogram_base_value.first << " "<< *parallelogram_base_value.second - *parallelogram_base_value.first << endl;
//	//cout << *parallelogram_base_value0.second << " "<< *parallelogram_base_value0.first <<" " <<*parallelogram_base_value0.second - *parallelogram_base_value0.first <<endl;
//	//assert(*parallelogram_base_value.second - *parallelogram_base_value.first == *parallelogram_base_value0.second -  *parallelogram_base_value0.first);
//	//auto test = temp_coefficient.parallelogram_height * temp_coefficient.rectangle_width;
//
//	return temp_coefficient.segment_area = temp_coefficient.parallelogram_height * temp_coefficient.segment_width;
//}
//
////************************************
//// Method:getAAndBByPLASegment
//// Qualifier: Get PLA a&b
//// Input:
//// Output:
//// date:191016
//// author:
////************************************
//TEMPLATE
//template<typename T>
//void PLA_QUAL::getAAndBByPLASegment(const vector<DataType>& const original_time_series_vector, T& const temp_coefficient) {//191016 get a&b for one segment
//#ifdef _DEBUG
//	assert(temp_coefficient.rectangle_width > 0 && temp_coefficient.rectangle_width != INF && temp_coefficient.right_endpoint != INF);
//#endif
//	const int segment_left_id = temp_coefficient.right_endpoint - temp_coefficient.rectangle_width + 1;
//
//	if (temp_coefficient.rectangle_width == 1) {
//
//		//a
//		temp_coefficient.apla.a_minuend = 0;//(l-1)/2
//		temp_coefficient.apla.a_divisor = 0;//l(l-1)(l+1)
//		//b
//		temp_coefficient.apla.b_minuend = 1;//2l-1
//		temp_coefficient.apla.b_divisor = 2;//l(l+1)
//
//		temp_coefficient.apla.a = 0;
//		temp_coefficient.apla.b = original_time_series_vector[temp_coefficient.right_endpoint];
//		return;
//	}
//	if (temp_coefficient.rectangle_width == 2) {
//#ifdef _DEBUG
//		assert(temp_coefficient.right_endpoint - segment_left_id == 1);
//#endif
//		//a
//		temp_coefficient.apla.a_minuend = 0.5;//(l-1)/2
//		temp_coefficient.apla.a_divisor = 6;//l(l-1)(l+1)
//		//b
//		temp_coefficient.apla.b_minuend = 3;//2l-1
//		temp_coefficient.apla.b_divisor = 6;//l(l+1)
//
//		temp_coefficient.apla.a = (original_time_series_vector[int(temp_coefficient.right_endpoint)] - original_time_series_vector[segment_left_id]);
//		temp_coefficient.apla.b = original_time_series_vector[segment_left_id];
//#ifdef _DEBUG
//		assert(temp_coefficient.apla.a != INF && temp_coefficient.apla.b != INF);
//#endif
//		return;
//	}
//	//coefficient of equation
//
//	auto a_sum = 0.0;
//	auto b_sum = 0.0;
//	double variable_id = NULL; //[0-segment_length)
//	int array_id = NULL;//time series id
//	//a
//	temp_coefficient.apla.a_minuend = (temp_coefficient.rectangle_width - 1) / 2.0;//(l-1)/2
//	temp_coefficient.apla.a_divisor = (temp_coefficient.rectangle_width - 1) * (temp_coefficient.rectangle_width + 1) * temp_coefficient.rectangle_width;//l(l-1)(l+1)
//	//b
//	temp_coefficient.apla.b_minuend = 2.0 * temp_coefficient.rectangle_width - 1;//2l-1
//	temp_coefficient.apla.b_divisor = (temp_coefficient.rectangle_width + 1) * temp_coefficient.rectangle_width;//l(l+1)
//
//#ifdef _DEBUG
//	assert(temp_coefficient.apla.a_divisor != 0);
//	assert(temp_coefficient.apla.b_divisor != 0);
//#endif
//	//cout << "endpoint: " << temp_coefficient.right_endpoint << " length: " << temp_coefficient.rectangle_width << endl;
//	for (variable_id = 0, array_id = segment_left_id; array_id <= temp_coefficient.right_endpoint; array_id++, variable_id++) {
//		//cout << "variable_id " << variable_id << "  array_id "<< array_id <<endl;
//		a_sum += (variable_id - temp_coefficient.apla.a_minuend) * original_time_series_vector[array_id];
//		b_sum += (temp_coefficient.apla.b_minuend - variable_id * 3.0) * original_time_series_vector[array_id];
//	}
//	temp_coefficient.apla.a = 12.0 * a_sum / temp_coefficient.apla.a_divisor;
//	temp_coefficient.apla.b = 2.0 * b_sum / temp_coefficient.apla.b_divisor;
//
//#ifdef _DEBUG
//	assert(variable_id == temp_coefficient.rectangle_width && array_id == temp_coefficient.right_endpoint + 1);
//	assert(temp_coefficient.apla.a != INF && temp_coefficient.apla.b != INF);
//#endif
//}
//
////191215
////************************************
//// Method:get_a_b_by_endpoint
//// Qualifier: Get PLA a&b
//// Input:
//// Output:
//// date:191215
//// author:
////************************************
//TEMPLATE
//template<typename T>
//T& PLA_QUAL::get_a_b_by_endpoint(const vector<DataType>& const original_time_series_vector, T& const temp_coefficient) {//191215 get a&b for one segment
//#ifdef _DEBUG
//	assert(temp_coefficient.rectangle_width > 0 && temp_coefficient.rectangle_width != INF && temp_coefficient.right_endpoint != INF);
//#endif
//
//	int segment_left_id = temp_coefficient.right_endpoint - temp_coefficient.rectangle_width + 1;
//
//	if (temp_coefficient.rectangle_width == 1) {
//		temp_coefficient.apla.a = 0;
//		temp_coefficient.apla.b = original_time_series_vector[temp_coefficient.right_endpoint];
//		return temp_coefficient;
//	}
//	if (temp_coefficient.rectangle_width == 2) {
//#ifdef _DEBUG
//		assert(temp_coefficient.right_endpoint - segment_left_id == 1);
//#endif
//		//*temp_coefficient.apla.a = (original_time_series[int(temp_coefficient.right_endpoint)] - original_time_series[segment_left_id]);
//		//*temp_coefficient.apla.b = original_time_series[segment_left_id];
//		temp_coefficient.apla.a = (original_time_series_vector[int(temp_coefficient.right_endpoint)] - original_time_series_vector[segment_left_id]);
//		temp_coefficient.apla.b = original_time_series_vector[segment_left_id];
//#ifdef _DEBUG
//		assert(temp_coefficient.apla.a != INF && temp_coefficient.apla.b != INF);
//#endif
//		return temp_coefficient;
//	}
//
//	temp_coefficient.apla.a = (original_time_series_vector[int(temp_coefficient.right_endpoint)] - original_time_series_vector[segment_left_id]) / temp_coefficient.rectangle_width;
//	temp_coefficient.apla.b = original_time_series_vector[segment_left_id];
//
//#ifdef _DEBUG
//	assert(temp_coefficient.apla.a != INF && temp_coefficient.apla.b != INF);
//#endif
//	return temp_coefficient;
//}
//
////************************************
//// Method:getAAndBByPLASegment
//// Qualifier: Get PLA a&b
//// Input:a1,b1,sum1,l1 & st
//// Output:a,b, sum, l.
//// date:191021
//// author:
////************************************
//TEMPLATE
//template<typename T>
//T& PLA_QUAL::getAAndBByAccumulation(const vector<DataType>& const original_time_series_vector, const T& const last_segment, T& const temp_coefficient) {//191021 the length of segment is added one by one 
//	/*...................................................................................................*/
//#ifdef _DEBUG
//	assert(temp_coefficient.rectangle_width > 0 && temp_coefficient.rectangle_width != INF && temp_coefficient.right_endpoint > 0 && temp_coefficient.right_endpoint != INF);
//#endif
//	/*...................................................................................................*/
//
//	/*####################################################################################################################################################*/
//	int segment_left_id = temp_coefficient.right_endpoint - temp_coefficient.rectangle_width + 1;
//
//	if (temp_coefficient.rectangle_width == 1) {
//		assert(0);
//		temp_coefficient.apla.a = 0;
//		temp_coefficient.apla.b = original_time_series_vector[temp_coefficient.right_endpoint];
//		return temp_coefficient;
//	}
//
//	if (temp_coefficient.rectangle_width == 2) {
//
//		/*.......................................................................................................................................................*/
//#ifdef _DEBUG
//		assert(temp_coefficient.right_endpoint - segment_left_id == 1);
//#endif
//		/*.......................................................................................................................................................*/
//		
//		temp_coefficient.apla.a_minuend = 0.5;//(l-1)/2
//		temp_coefficient.apla.a_divisor = 6;//l(l-1)(l+1)
//		temp_coefficient.apla.b_minuend = 3;//2l-1
//		temp_coefficient.apla.b_divisor = 6;//l(l+1)
//
//		//200316 sum value
//		//temp_coefficient.sum_value = original_time_series_vector[segment_left_id] + original_time_series_vector[int(temp_coefficient.right_endpoint)];
//		temp_coefficient.apla.a = (original_time_series_vector[int(temp_coefficient.right_endpoint)] - original_time_series_vector[segment_left_id]);
//		temp_coefficient.apla.b = original_time_series_vector[segment_left_id];
//
//		/*.......................................................................................................................................................*/
//#ifdef _DEBUG
//		assert(temp_coefficient.apla.a != INF && temp_coefficient.apla.b != INF);
//#endif
//		/*.......................................................................................................................................................*/
//
//		return temp_coefficient;
//	}
//	/*####################################################################################################################################################*/
//
//	/*.......................................................................................................................................................*/
//#ifdef _DEBUG
//	assert(last_segment.rectangle_width > 1 && last_segment.rectangle_width == temp_coefficient.rectangle_width - 1 && last_segment.rectangle_width > 0 && last_segment.rectangle_width != INF && last_segment.right_endpoint > 0 && last_segment.right_endpoint != INF && last_segment.apla.a != INF && last_segment.apla.b != INF);
//#endif
//	/*.......................................................................................................................................................*/
//
//	//200316 sum value
//	//temp_coefficient.sum_value = last_segment.sum_value + original_time_series_vector[temp_coefficient.right_endpoint];// sum value
//
//	//200316
//	//double temp_length = temp_coefficient.rectangle_width + 1;//l + 1
//	//double temp_coefficient0 = last_segment.rectangle_width - 1;// l1 - 1
//	//double temp_coefficient1 = last_segment.rectangle_width / temp_length;// l1 / l+1
//	//double temp_coefficient2 = temp_length * temp_coefficient.rectangle_width; // l*(l+1)
//	//temp_coefficient.apla.a = last_segment.apla.a * (temp_coefficient0 / temp_length) + 6 * (last_segment.rectangle_width * original_time_series_vector[temp_coefficient.right_endpoint] - last_segment.sum_value) / (last_segment.rectangle_width * temp_coefficient2);
//	//temp_coefficient.apla.b = last_segment.apla.b * temp_coefficient1 + (4 * last_segment.sum_value + 2 * (1 - last_segment.rectangle_width) * original_time_series_vector[temp_coefficient.right_endpoint]) / temp_coefficient2;
//
//	/*#####################################################    200827 If width >= 3   ##############################################################################*/
//
//	temp_coefficient.apla.a_minuend = (temp_coefficient.rectangle_width - 1) / 2.0;//(l-1)/2
//	temp_coefficient.apla.a_divisor = (temp_coefficient.rectangle_width - 1) * (temp_coefficient.rectangle_width + 1) * temp_coefficient.rectangle_width;//l(l-1)(l+1)
//	temp_coefficient.apla.b_minuend = 2.0 * temp_coefficient.rectangle_width - 1;//2l-1
//	temp_coefficient.apla.b_divisor = (temp_coefficient.rectangle_width + 1) * temp_coefficient.rectangle_width;//l(l+1)
//
//	//double temp_length = temp_coefficient.rectangle_width + 1;//l + 1
//	 const long double temp_coefficient0 = last_segment.rectangle_width - 1;// l1 - 1
//	//double temp_coefficient1 = last_segment.rectangle_width / temp_length;// l1 / l+1
//	//double temp_coefficient2 = temp_length * temp_coefficient.rectangle_width; // l*(l+1)
//	//temp_coefficient.apla.a = last_segment.apla.a * (temp_coefficient0 / temp_length) + 6 * (last_segment.rectangle_width * original_time_series_vector[temp_coefficient.right_endpoint] - last_segment.sum_value) / (last_segment.rectangle_width * temp_coefficient2);
//	//temp_coefficient.apla.b = last_segment.apla.b * temp_coefficient1 + (4 * last_segment.sum_value + 2 * (1 - last_segment.rectangle_width) * original_time_series_vector[temp_coefficient.right_endpoint]) / temp_coefficient2;
//
//	//double divide_coefficients = temp_coefficient.rectangle_width * (temp_coefficient.rectangle_width + 1);
//	temp_coefficient.apla.a = ((temp_coefficient.rectangle_width - 3) * temp_coefficient0 * last_segment.apla.a + 6 * (original_time_series_vector[temp_coefficient.right_endpoint] - last_segment.apla.b)) / temp_coefficient.apla.b_divisor;
//	temp_coefficient.apla.b = (2 * temp_coefficient0 * (last_segment.apla.a * last_segment.rectangle_width - original_time_series_vector[temp_coefficient.right_endpoint]) + (temp_coefficient.rectangle_width + 4) * last_segment.rectangle_width * last_segment.apla.b) / temp_coefficient.apla.b_divisor;
//	/*###############################################################################################################################################################*/
//
//	/*.......................................................................................................................................................*/
//#ifdef _DEBUG
//	assert(temp_coefficient.apla.a != INF && temp_coefficient.apla.b != INF);
//	T test_coefficient = temp_coefficient;
//	//test_coefficient.right_endpoint = temp_coefficient.right_endpoint;
//	//test_coefficient.rectangle_width = temp_coefficient.rectangle_width;
//	PLA_QUAL::getAAndBByPLASegment(original_time_series_vector, test_coefficient);
//	//cout << *temp_coefficient.apla.a << " : " << temp_a << "      " << *temp_coefficient.apla.b << " : " << temp_b << endl;
//	assert(fabs(float(temp_coefficient.apla.a) - float(test_coefficient.apla.a)) <= MIN_D && fabs(float(temp_coefficient.apla.b) - float(test_coefficient.apla.b)) <= MIN_D);
//#endif
//	/*.......................................................................................................................................................*/
//
//	return temp_coefficient;
//}
//
////************************************
//// Method:getSegmentMaxDifference
//// Qualifier: Get Max Difference of one segment
//// Input:
//// Output:
//// date:191016
//// author:
////************************************
//TEMPLATE
//template<typename T>
//long double PLA_QUAL::getSegmentMaxDifference(const vector<DataType>& const original_time_series_vector, const T& const segment_pla) {
//#ifdef _DEBUG
//	//assert(area_vector.size() == input_argument.point_dimension);
//	assert(segment_pla.rectangle_width != INF && segment_pla.right_endpoint != INF && segment_pla.apla.a != INF && segment_pla.apla.b != INF);
//#endif
//
//	int variable_id = NULL; //[0-segment_length
//	int array_id = 0;
//	int segment_left_id = segment_pla.right_endpoint - segment_pla.rectangle_width + 1;
//	long double max_difference = -INF;
//
//#ifdef _DEBUG
//	assert(segment_left_id > -1);
//#endif
//	for (variable_id = 0, array_id = segment_left_id; array_id <= segment_pla.right_endpoint; variable_id++, array_id++) {
//		max_difference = max(max_difference, fabs(segment_pla.apla.a * variable_id + segment_pla.apla.b - original_time_series_vector[array_id]));
//	}
//#ifdef _DEBUG
//	assert(variable_id == segment_pla.rectangle_width && array_id == segment_pla.right_endpoint + 1 && max_difference != -INF);
//#endif
//	/*=================================Evaluation Result=====================================================================*/
////#ifdef _DEBUG
////	typename TOOL::getDeviation(input_argument, original_time_series, reconstruct_time_series, output_argument);
////	assert(sum_deviation == output_argument.sum_deviation);
////#endif
//	/*..............................................................................................................*/
//	//TOOL::deleteArray(reconstruct_time_series);
//	return max_difference;
//}
//
////************************************
//// Method:getSegmentSumDifference
//// Qualifier:  get segment sum difference
//// Input:
//// Output:
//// date:191018
//// author:
////************************************
//TEMPLATE
//template<typename T>
//double PLA_QUAL::getSegmentSumDifferenceSquare(const vector<DataType>& const original_time_series_vector, const T& const segment_pla, vector<DataType>& const rescontruct_time_series) {
//#ifdef _DEBUG
//	//assert(area_vector.size() == input_argument.point_dimension);
//	assert(segment_pla.apla.a != INF && segment_pla.apla.b != INF);
//	T test_segment = segment_pla;
//	PLA_QUAL::getAAndBByPLASegment(original_time_series_vector, test_segment);
//	//assert(fabs(float(segment_pla.apla.a) - float(test_segment.apla.a)) <= MIN_D && float(segment_pla.apla.b) == float(test_segment.apla.b));
//#endif
//
//	int variable_id = NULL; //[0-segment_length
//	int array_id = 0;
//	int segment_left_id = segment_pla.right_endpoint - segment_pla.rectangle_width + 1;
//	double sum_difference = 0;
//
//#ifdef _DEBUG
//	assert(segment_left_id > -1);
//#endif
//	for (variable_id = 0, array_id = segment_left_id; array_id <= segment_pla.right_endpoint; variable_id++, array_id++) {
//		sum_difference += pow(segment_pla.apla.a * variable_id + segment_pla.apla.b - original_time_series_vector[array_id], 2.0);// difference ^2
//		/*--------------------------------------------------------------------------------------------------------------------------------------*/
//
//		rescontruct_time_series.push_back(segment_pla.apla.a * variable_id + segment_pla.apla.b);
//		//cout << segment_pla.apla.a * variable_id + segment_pla.apla.b << ",";
//
//		/*--------------------------------------------------------------------------------------------------------------------------------------*/
//	}
//#ifdef _DEBUG
//	assert(variable_id == segment_pla.rectangle_width && array_id == segment_pla.right_endpoint + 1);
//#endif
//	/*=================================Evaluation Result=====================================================================*/
////#ifdef _DEBUG
////	typename TOOL::getDeviation(input_argument, original_time_series, reconstruct_time_series, output_argument);
////	assert(sum_deviation == output_argument.sum_deviation);
////#endif
//	/*===============================================================================================================*/
//	//TOOL::deleteArray(reconstruct_time_series);
//	return sum_difference;
//}
////************************************
//// Method:getSegmentSumDifference
//// Qualifier:  get segment sum difference
//// Input:
//// Output:
//// date:191018
//// author:
////************************************
//TEMPLATE
//template<typename T>
//double PLA_QUAL::getSegmentSumDifferenceSquare(const vector<DataType>& const original_time_series_vector, const T& const segment_pla) {
//#ifdef _DEBUG
//	//assert(area_vector.size() == input_argument.point_dimension);
//	assert(segment_pla.a != INF && segment_pla.b != INF);
//	T test_segment = segment_pla;
//	PLA_QUAL::getAAndBByPLASegment(original_time_series_vector, test_segment);
//	assert(fabs(float(segment_pla.a) - float(test_segment.a)) <= MIN_D && float(segment_pla.b) == float(test_segment.b));
//#endif
//
//	int variable_id = NULL; //[0-segment_length
//	int array_id = 0;
//	int segment_left_id = segment_pla.right_endpoint - segment_pla.segment_width + 1;
//	double sum_difference = 0;
//
//#ifdef _DEBUG
//	assert(segment_left_id > -1);
//#endif
//	for (variable_id = 0, array_id = segment_left_id; array_id <= segment_pla.right_endpoint; variable_id++, array_id++) {
//		sum_difference += pow(segment_pla.a * variable_id + segment_pla.b - original_time_series_vector[array_id], 2.0);// difference ^2
//		/*--------------------------------------------------------------------------------------------------------------------------------------*/
//#ifdef _DEBUG
//		//rescontruct_time_series.push_back(segment_pla.a * variable_id + segment_pla.b);
//		//cout << segment_pla.a * variable_id + segment_pla.b << ",";
//#endif
//		/*--------------------------------------------------------------------------------------------------------------------------------------*/
//	}
//#ifdef _DEBUG
//	assert(variable_id == segment_pla.segment_width && array_id == segment_pla.right_endpoint + 1);
//#endif
//	/*=================================Evaluation Result=====================================================================*/
////#ifdef _DEBUG
////	typename TOOL::getDeviation(input_argument, original_time_series, reconstruct_time_series, output_argument);
////	assert(sum_deviation == output_argument.sum_deviation);
////#endif
//	/*===============================================================================================================*/
//	//TOOL::deleteArray(reconstruct_time_series);
//	return sum_difference;
//}
//
////************************************
//// Method:getPLASumDeviation
//// Qualifier:  get segment sum difference, Except for first element. begin() + 1
//// Input:
//// Output:
//// date:191018
//// author:
////************************************
//TEMPLATE
//template<typename T>
//double PLA_QUAL::getPLASumDeviation(const vector<DataType>& const original_time_series_vector, const vector<T>& const segment_vector) {
//#ifdef _DEBUG
//	assert(!segment_vector.empty());
//	cout << "reconstruct time series: " << segment_vector.size() << endl;
//#endif
//	vector<DataType> rescontruct_time_series;
//
//	double sum_deviation = 0;
//	for_each(segment_vector.begin(), segment_vector.end(), [&](auto segment) {
//#ifdef _DEBUG
//		assert(segment.rectangle_width != INF && segment.right_endpoint != INF && segment.apla.a != INF && segment.apla.b != INF);
//#endif
//		// << *segment.apla.a << " "<< *segment.apla.b <<endl;
//		sum_deviation += getSegmentSumDifferenceSquare(original_time_series_vector, segment, rescontruct_time_series);//sum += difference^2
//	});
//#ifdef _DEBUG
//	//cout << endl;
//	double temp_deviation = sqrtl(sum_deviation);
//	double test_deviation = TOOL::getDeviation(original_time_series_vector, rescontruct_time_series);
//	assert(original_time_series_vector.size() == rescontruct_time_series.size() && float(temp_deviation) == float(test_deviation));
//#endif
//	return std::sqrtl(sum_deviation);
//}
//
////191206 Linked List get segment sum difference,Except for first element. begin() + 1
////************************************
//// Method:getPLASumDeviation
//// Qualifier:  get segment sum difference, Except for first element. begin() + 1
//// Input: Linked list
//// Output:
//// date:191206
//// author:
////************************************
//TEMPLATE
//template<typename T>
//double PLA_QUAL::getPLASumDeviation(const vector<DataType>& const original_time_series_vector, const DoublyLinkedList<T>& const segment_list) {
//#ifdef _DEBUG
//	assert(!segment_list.empty());
//	//cout << "reconstruct time series: " << segment_list.size() << endl;
//#endif
//	vector<DataType> rescontruct_time_series;
//
//	double sum_deviation = 0;
//	for_each(segment_list.begin(), segment_list.end(), [&](auto segment) {
//#ifdef _DEBUG
//		assert(segment.rectangle_width != INF && segment.right_endpoint != INF && segment.apla.a != INF && segment.apla.b != INF);
//#endif
//		// << *segment.apla.a << " "<< *segment.apla.b <<endl;
//		sum_deviation += getSegmentSumDifferenceSquare(original_time_series_vector, segment, rescontruct_time_series);//sum += difference^2
//	});
//	//#ifdef _DEBUG
//		//cout << endl;
//	double temp_deviation = sqrtl(sum_deviation);
//	double test_deviation = TOOL::getDeviation(original_time_series_vector, rescontruct_time_series);
//	//cout << temp_deviation << "       " << test_deviation << endl;
//	assert(original_time_series_vector.size() == rescontruct_time_series.size() && float(temp_deviation) == float(test_deviation));
//	//#endif
//	rescontruct_time_series.clear();
//	rescontruct_time_series.shrink_to_fit();
//	return std::sqrtl(sum_deviation);
//}
//
//TEMPLATE
//template<typename T, typename Y, typename U>
//long double PLA_QUAL::getPLASumDeviation(const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const segment_list, U& const result_collection) {
//#ifdef _DEBUG
//	assert(!segment_list.empty());
//	//cout << "reconstruct time series: " << segment_list.size() << endl;
//#endif
//	vector<DataType> rescontruct_time_series;
//
//	result_collection.max_deviation = 0;
//	result_collection.max_deviation_multiple_width = 0;
//	result_collection.sum_deviation = 0;
//	long double sum_deviation = 0;
//	long double max_deviation = 0;
//	for_each(segment_list.begin(), segment_list.end(), [&](auto segment) {
//
//#ifdef _DEBUG
//		assert(segment.rectangle_width != INF && segment.right_endpoint != INF && segment.apla.a != INF && segment.apla.b != INF);
//#endif
//		// << *segment.apla.a << " "<< *segment.apla.b <<endl;
//		sum_deviation += getSegmentSumDifferenceSquare(original_time_series_vector, segment, rescontruct_time_series);//sum += difference^2
//		max_deviation = getSegmentMaxDifference(original_time_series_vector, segment);
//		result_collection.max_deviation += max_deviation;
//		result_collection.max_deviation_multiple_width += max_deviation * segment.rectangle_width;
//	});
//	//#ifdef _DEBUG
//		//cout << endl;
//	long double temp_deviation = sqrtl(sum_deviation);
//	double test_deviation = TOOL::getDeviation(original_time_series_vector, rescontruct_time_series);
//	//cout << temp_deviation << "       " << test_deviation << endl;
//	assert(original_time_series_vector.size() == rescontruct_time_series.size() && float(temp_deviation) == float(test_deviation));
//	//#endif
//	rescontruct_time_series.clear();
//	rescontruct_time_series.shrink_to_fit();
//	result_collection.sum_deviation = std::sqrtl(sum_deviation);
//	return std::sqrtl(sum_deviation);
//}
//
////************************************
//// Method:recordStartTime
//// Qualifier: Get PLA a&b
//// Input:
//// Output:
//// date:191016
//// author:
////************************************
//TEMPLATE
//void PLA_QUAL::recordStartTime(TIME& time) {
//	LARGE_INTEGER whole_first_f;    //timer frequency
//	QueryPerformanceFrequency(&whole_first_f);
//	time.dqFrequency = (double)whole_first_f.QuadPart;
//	QueryPerformanceCounter(&time.time_start);
//}
//
//TEMPLATE
//double PLA_QUAL::recordFinishTime(TIME& time, double& whole_first_run_time) {
//	QueryPerformanceCounter(&time.time_over);    //Finish recording time  5
//	whole_first_run_time = 1000000 * (time.time_over.QuadPart - time.time_start.QuadPart) / time.dqFrequency;
//	//cout << whole_first_run_time << endl;
//
//	return whole_first_run_time;
//}
//
//TEMPLATE
//void PLA_QUAL::writeResult(const INPUT_ARGUMENT& input_argument, const string& write_file_name, list<ORIGINAL_TIME_SERIES_PAIR>& result, priority_queue<ORIGINAL_TIME_SERIES_PAIR, vector<ORIGINAL_TIME_SERIES_PAIR>, priorityDistanceEUC >& q_base_queue) {
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
//	double count = 0;
//	time_t now = time(0);// system time
//	char* dt = ctime(&now);// from now to string
//	ofstream outfile(write_file_name + ".txt", ios::app);
//	outfile << endl << dt << "file:" << input_argument.read_file_name << "; " << write_file_name << endl;
//	outfile << " n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
//	//outfile << "  count apca point = " << g_n_account_apca_point << " times, p= " << g_n_account_apca_point / input_argument.point_number << endl;
//	//outfile << "  APCA Memory = " << memory_account[0] + memory_account[1] + memory_account[2] + memory_account[3] << "  RTree Memeory = " << memory_account[4] + memory_account[5] + memory_account[6] + memory_account[7] + memory_account[8] << "  KNN Memory = " << memory_account[9] << endl;
//	for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//		if (it->d_dist == q_base_queue.top().d_dist) count++;
//		outfile << "id: " << it->original_time_series_id << " , " << q_base_queue.top().original_time_series_id << "; dist: " << it->d_dist << " , " << q_base_queue.top().d_dist << endl;
//		q_base_queue.pop();
//	}
//
//	outfile << "accuracy: " << count / double(result.size()) << endl;
//	outfile.close();
//}



#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////