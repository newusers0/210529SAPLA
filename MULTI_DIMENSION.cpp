#pragma once
#ifndef _MULTI_DIMENSION_CPP_
#define _MULTI_DIMENSION_CPP_

#include "MULTI_DIMENSION.h"
/////////MULTI///////////////////////////////////////////////////////////////////////////////////////
//
//TEMPLATE
//struct MULTI::MULTI_DIM {
//	DataType** multi_dimension;
//	DataType* projected_dimension;
//};
//
////***************************************************************
//// Method:APLA_NODE_PAIR
//// Qualifier:  for KNN search in internal node.
//// Input: 1 d_dist distance, 2 time seris( APLA point ) id. 3 pointer to APCA point. 4 pointer to RTree node.
//// Output: 
//// date:191112
//// author:
////***************************************************************
//TEMPLATE
//template<typename T>
//struct MULTI::APLA_NODE_PAIR {//queue data structure.
//	double d_dist = INF;//dist.
//	int original_time_series_id = INF;// APCA point ID.
//	DoublyLinkedList<T>* approximation_pointer = nullptr;//data point. if NULL, this is internal node.
//	RTREE::Node* p_rtree_node = nullptr;//subNode, if NULL, this is apca point.
//};
//
////***************************************************************
//// Method:RTREE_NODE_PAIR
//// Qualifier:  for KNN search in internal node.
//// Input: 1 d_dist distance, 2 time seris( APLA point ) id. 3 pointer to RTree node.
//// Output: 
//// date:191128
//// author:
////***************************************************************
//TEMPLATE
//struct MULTI::RTREE_NODE_PAIR {
//	double d_dist = INF;//dist.
//	int original_time_series_id = INF;// time series ID.
//	RTREE::Node* p_rtree_node = nullptr;//subNode, if NULL, this is apca point.
//};
//
////***************************************************************
//// Method: RTREE_NODE_PAIR_PARTITION
//// Qualifier:
//// Input: 1 d_dist distance, 2 time seris( APLA point ) id. 3 pointer to RTree node.
//// Output: 
//// date:210618
//// author:
////***************************************************************
//TEMPLATE
//template<typename T, typename Y, typename U>
//struct MULTI::RTREE_NODE_PAIR_PARTITION {
//	T d_dist = INF;//dist.
//	Y original_time_series_id = INF;// time series ID.
//	U* p_rtree_node = nullptr;//subNode, if NULL, this is apca point.
//};
//
////***************************************************************
//// Method:priorityIncrement
//// Qualifier: for priority queue, top is smallest, small to big
//// Input: 
//// Output: 
//// date:191112
//// author:
////***************************************************************
//TEMPLATE
//template<typename T>
//struct MULTI::priorityIncrement {//small to big
//	bool operator ()(const T& const a, const T& const b) {
//		return a.d_dist > b.d_dist;
//		//return a.distance > b.distance;
//	}
//};
//
//TEMPLATE
//MULTI::MULTI_DIMENSION(const int& n, const int& N, const int& d, string*& const multi_file_name, const int& representation_option) {
//	input_argument.remainder = int(n) % int(N);//For PLA
//	double integerDividend = n - input_argument.remainder;
//	input_argument.segment_length_second = integerDividend / N;
//	input_argument.segment_length_first = input_argument.segment_length_second + 1;
//	assert(input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
//
//	input_argument.time_series_length = n;
//	input_argument.point_dimension = N;//N = N;//point_dimension for APCA & PLA n=m+1
//	input_argument.degree_m = N - 1;//degree of polymonail n=m+1
//	input_argument.arity_d = d;
//	input_argument.point_multi_single_dimension = d * N;
//	input_argument.read_multiple_file_name = multi_file_name;
//	input_argument.representation_option = representation_option;
//}
//
////***************************************************************
//// Method:MULTI_DIMENSION
//// Qualifier: initial of ever coefficient
//// 1£º the N of different Approximation is different at initialization period. APLA = 4, PAA : 12 , PLA: 6, APCA 6. Cheby 12, ICDE07: 4
//// Input:
//// Output: 
//// date:191107
//// author:
////***************************************************************
//TEMPLATE
//MULTI::MULTI_DIMENSION(const int& const n, const int& const N, const int& const file_id, const int& const d, const int& const point_number, const int& const query_time_series_id, const int& const rtree_max_nodes, const int& const K, const int& const representation_option, string*& const multi_file_name, string*& const write_file_name) {
//#ifdef _DEBUG
//	assert(representation_option != INF && representation_option > 0 && point_number >= K && file_id != INF && d > 0 && N > 0 && input_argument.time_series_length % 2 == 0 && query_time_series_id > -1 && query_time_series_id < point_number&& point_number != INF);//n is even
//	cout << "============ File id: " << file_id + 1 << ", representation option: " << representation_option << ", n: " << n << ", Original N: " << N << ", K: " << K << ", point number: " << point_number << ", query time series id: " << query_time_series_id << ", Rtree max sub node number: " << rtree_max_nodes << "============================" << endl;
//#endif
//	input_argument.representation_option = representation_option;
//	/*====================================   inital N   =============================================*/
//	switch (input_argument.representation_option) {
//	case 1: //APLA
//	case 6: {//ICDE07
//		//input_argument.initial_N = initial_N / 3;
//		input_argument.point_dimension = N / 3;
//		break;
//	}
//	case 2: //PLA
//	case 3: {//APCA
//		input_argument.point_dimension = N / 2;
//		break;
//	}
//	case 4: //PAA
//	case 7: //PAALM
//	case 5: {//Cheby
//		input_argument.point_dimension = N;
//		break;
//	}
//	default:
//		assert(0);
//		break;
//	}
//#ifdef _DEBUG
//	assert(input_argument.point_dimension > 0);
//#endif
//	/*==============================================================================================*/
//
//	input_argument.time_series_length = n;
//	input_argument.degree_m = input_argument.point_dimension - 1;//degree of chebyshev polymonail n=m+1
//	input_argument.file_id = file_id;
//	input_argument.arity_d = d;//dimension of time series, >0
//	input_argument.point_multi_single_dimension = d * input_argument.point_dimension;// project high dimension to single dimension For Rtree point dimension and initializaiton
//	input_argument.point_multi_single_length = d * n;
//
//	/*===================================For PLA coefficient====================================================*/
//	input_argument.remainder = int(n) % int(input_argument.point_dimension);//For PLA
//	input_argument.segment_length_second = (n - input_argument.remainder) / input_argument.point_dimension;
//	input_argument.segment_length_first = input_argument.segment_length_second + 1;
//#ifdef _DEBUG
//	assert(input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
//#endif
//	/*==========================================================================================================*/
//	input_argument.query_time_series_id = query_time_series_id;
//	input_argument.point_number = point_number;
//	input_argument.rtree_max_nodes = rtree_max_nodes;
//	input_argument.K = K;
//	input_argument.read_multiple_file_name = multi_file_name;
//	//input_argument.read_file_name = read_file_name;
//	input_argument.write_file_name = write_file_name;
//
//	//*********Initial object****************//
//	initialMultiDimension();
//	//***************************************//
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
//	input_argument.result_accuracy = 0;
//}
//
////***************************************************************
//// Method:MULTI_DIMENSION
//// Qualifier: initial of ever coefficient
//// 1£º the N of different Approximation is different at initialization period. APLA = 4, PAA : 12 , PLA: 6, APCA 6. Cheby 12, ICDE07: 4
//// Input:
//// Output: 
//// date:191107
//// author:
////***************************************************************
//TEMPLATE
//MULTI::MULTI_DIMENSION(const typename TOOL::DATA_SOURCE& const data_source, const int& const n, const int& const N, const double& const initial_N, const int& const file_id, const bool& const change_file, const int& const d, const int& const point_number, const int& const query_time_series_id, const int& const rtree_max_nodes, const int& const K, const int& const representation_option, string*& const multi_file_name, string*& const write_file_name) {
//
//	assert(representation_option != INF && representation_option > 0 && point_number >= K && file_id != INF && file_id >= 0 && data_source.time_series_dimension == d && d > 0 && N > 0 && query_time_series_id > -1 && query_time_series_id < point_number&& point_number != INF);//n is even
//	assert(query_time_series_id < point_number&& query_time_series_id >= 0);
//	assert(data_source.data_list.size() == data_source.single_file_number && data_source.data_type != INF && !data_source.data_list.empty() && data_source.single_file_number != INF && data_source.single_point_number != INF && data_source.single_time_series_length != INF);
//	assert(data_source.single_time_series_length == n && data_source.point_number == point_number && query_time_series_id < data_source.point_number&& K < data_source.point_number&& file_id <= data_source.file_operation_number);
//
//	input_argument.change_file = change_file;
//	this->data_source = data_source;
//	input_argument.representation_option = representation_option;
//	input_argument.time_series_length = n;
//	input_argument.file_id = file_id;
//	input_argument.arity_d = d;//dimension of time series, >0, 1 or 3
//
//	/*====================================   inital N   =============================================*/
//	switch (input_argument.representation_option) {
//	case 1: //APLA
//	case 6: //ICDE07
//	case 8: {//initial200706
//		/*------200330 initial_N Evaluation-------*/
//		if (initial_N != INF)
//			input_argument.initial_N = initial_N;
//		/*----------------------------------------*/
//		input_argument.point_dimension = N / 3;
//		break;
//	}
//	case 2: //PLA
//	case 3: {//APCA
//		input_argument.point_dimension = N / 2;
//		break;
//	}
//	case 4: //PAA
//	case 7: //PAALM
//	case 5: //Cheby
//	case 9: {//SAX
//		input_argument.point_dimension = N;
//		break;
//	}
//	default:
//		assert(0);
//		break;
//	}
//#ifdef _DEBUG
//	assert(input_argument.point_dimension > 0 && this->data_source.multi_single_time_series_length == d * n);
//#endif
//	/*==============================================================================================*/
//
//
//	input_argument.point_multi_single_dimension = input_argument.point_dimension;
//	//input_argument.point_multi_single_dimension = d * input_argument.point_dimension;// project high dimension to single dimension For Rtree point dimension and initializaiton
//	input_argument.point_multi_single_length = this->data_source.multi_single_time_series_length;
//
//	input_argument.time_series_length = input_argument.point_multi_single_length;
//	input_argument.point_dimension = input_argument.point_multi_single_dimension;
//
//
//	/*===================================For Chebyshev coefficient====================================================*/
//	input_argument.degree_m = input_argument.point_dimension - 1;//degree of chebyshev polymonail n=m+1
//	/*==========================================================================================================*/
//	/*===================================For PLA coefficient====================================================*/
//	input_argument.remainder = int(input_argument.point_multi_single_length) % int(input_argument.point_dimension);//For PLA
//	input_argument.segment_length_second = (input_argument.point_multi_single_length - input_argument.remainder) / input_argument.point_dimension;
//	input_argument.segment_length_first = input_argument.segment_length_second + 1;
//#ifdef _DEBUG
//	assert(input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
//#endif
//	/*==========================================================================================================*/
//	input_argument.query_time_series_id = query_time_series_id;
//	input_argument.point_number = data_source.point_number;
//	input_argument.rtree_max_nodes = rtree_max_nodes;
//	input_argument.K = K;
//	input_argument.read_multiple_file_name = multi_file_name;
//	//input_argument.read_file_name = read_file_name;
//	input_argument.write_file_name = write_file_name;
//
//	/*==========================================Read file name================================================================*/
//	//TOOL::get_read_multi_file_address(this->data_source, input_argument.file_id);
//	switch (data_source.data_type) {
//	case 0: // 0 manipulate data set;    
//		assert(0);
//		break;
//	case 1:// 1 real data set, homogenous(single) data;
//
//		//switch (data_source.option_has_burst_data) {
//		//case 0: // 0 has burst data 21
//		//	input_argument.read_file_name = TOOL::getStringByID(TOOL::file_address, input_argument.file_id); //file has 21 files for DistLB.
//		//	
//		//	break;
//		//case 1: // 1 no burst data 17
//		//	input_argument.read_file_name = TOOL::getStringByID(TOOL::no_burst_single_file_address, input_argument.file_id);//No burst time series, has 17 files.
//		//	break;
//		//default:
//		//	assert(0);
//		//}
//		input_argument.read_file_name = data_source.read_file_address_vector.front();
//		assert(input_argument.read_file_name == data_source.read_file_address_vector.front() && d == 1);
//		break;
//	case 2:// 2 real data set, heterogeneous(mixed) data.
//		input_argument.read_file_name = data_source.read_file_address;
//		assert(input_argument.read_file_name == data_source.read_file_address_vector.front() && d == 1);
//#ifdef _DEBUG
//		assert(!data_source.data_list.empty() && data_source.single_point_number * data_source.data_list.size() == input_argument.point_number && input_argument.point_number == data_source.point_number);
//#endif
//		break;
//	case 3:// 3 multi dimenson single data.
//		assert(data_source.read_file_address_vector.size() == data_source.time_series_dimension && d > 1);
//		break;
//	case 4:// 4 multi mixed data
//		assert(0);
//		break;
//	case 5:
//		assert(0);
//		break;
//	default:
//		assert(0);
//		break;
//	}
//	/*===========================================================================================================*/
//
//	//*********Initial object****************//
//	initialMultiDimension();
//	//***************************************//
//
//	/*---------   Evaluation     -------*/
//	input_argument.representation_time = 0;
//	input_argument.build_rtree_time = 0;
//	input_argument.knn_total_time = 0;
//	input_argument.whole_run_time = 0;
//	input_argument.IO_cost = 0;
//	input_argument.pruning_power = 0;
//	input_argument.result_accuracy = 0;
//	/*----------------------------------*/
//
//	input_argument.sum_distance_euc = 0.0;
//	input_argument.build_rtree_time = 0.0;
//	input_argument.approximation_query_time = 0.0;
//	input_argument.knn_rest_part_time = 0.0;
//	input_argument.navigate_index_time = 0;// navigate time
//	input_argument.distance_lowbound_time = 0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0;// distance euclidean time
//
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//		cout << "=====================!!!!!!!!  Begin  SAPLA  !!!!!!!!=====================\n";
//		break;
//	case 2: //PLA
//		cout << "=====================!!!!!!!!  Begin  PLA   !!!!!!!!====================\n";
//		break;
//	case 3: //APCA
//		cout << "=====================!!!!!!!!  Begin  APCA  !!!!!!!!====================\n";
//		break;
//	case 4: //PAA
//		cout << "=====================!!!!!!!!  Begin  PAA   !!!!!!!!=====================\n";
//		break;
//	case 5://CHEBY
//		cout << "=====================!!!!!!!!  Begin  CHEBY  !!!!!!!!=====================\n";
//		break;
//	case 6: //ICDE07
//		cout << "=====================!!!!!!!!  Begin  ICDE07  !!!!!!!!================\n";
//		break;
//	case 7: //PAALM
//		cout << "=====================!!!!!!!!  Begin  PAALM  !!!!!!!!================\n";
//		break;
//	case 8://Initial_part
//		cout << "=====================!!!!!!!!  Begin  Initial 200706  !!!!!!!!=====================\n";
//		break;
//	case 9://SAX
//		cout << "=====================!!!!!!!!  Begin  SAX  !!!!!!!!=====================\n";
//		break;
//	default: assert(0);
//		break;
//	}
//	cout << "File id: " << input_argument.file_id + 1 << ", change file: " << input_argument.change_file << ", K: " << K << ", Original  N: " << N << ", representation option: " << input_argument.representation_option << ", n: " << input_argument.point_multi_single_length << ", point number: " << input_argument.point_number << ", query time series id: " << query_time_series_id << ", Rtree max sub node number: " << rtree_max_nodes << "============================" << endl;
//
//}
//
////***************************************************************
//// Method: ~MULTI_DIMENSION
//// Qualifier: 
//// Input:
//// Output: 
//// date:
//// author:
////***************************************************************
//TEMPLATE
//MULTI::~MULTI_DIMENSION() {
//	assert(input_argument.arity_d > 0);
//	input_argument.time_series_length = INF;
//	input_argument.degree_m = INF;//point_dimension for APCA & PLA n=m+1
//	input_argument.arity_d = INF;
//	deleteMultiDimension();
//	RTree.RemoveAll();
//
//	/*--------- 210618 RTree_partition ---------*/
//	RTree_partition_object.RemoveAll();
//	/*------------------------------------------*/
//}
//
////***************************************************************
//// Method:initialMultiDimension
//// Qualifier: 
//// Input:
//// Output: 
//// date:
//// author:
////***************************************************************
//TEMPLATE
//void MULTI::initialMultiDimension() {
//	/*............................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.representation_option != NULL && input_argument.segment_length_second > 1 && input_argument.point_multi_single_dimension > 0 && input_argument.point_number != INF);
//#endif
//	/*............................................................................................*/
//
//	/*---------   Evaluation     -------*/
//	input_argument.representation_time = 0;
//	input_argument.build_rtree_time = 0;
//	input_argument.knn_total_time = 0;
//	input_argument.whole_run_time = 0;
//	input_argument.IO_cost = 0;
//	input_argument.pruning_power = 0;
//	input_argument.result_accuracy = 0;
//	output_argument.sum_deviation = 0;
//	/*----------------------------------*/
//
//	switch (input_argument.representation_option) {//1 MSPLA  2 PLA  3 APCA   4 PAA  5 CHEBY  6  ICDE07
//	case 1://MSPLA
//	case 6: //ICDE07
//	case 8: {//initial 200706
//		/*=============================================  APLA & ICDE07 ==============================================================================*/
//		//cout << "Begin APLA & ICDE07 representation: " << endl;
//		RTree.initialRTree(input_argument.point_multi_single_dimension * 2, input_argument.rtree_max_nodes, input_argument.representation_option);
//
//		/*------------------------------------------ 210618 RTree_partition ------------------------------------------------------------------------------*/
//		RTree_partition_object.initialRTree(input_argument.point_multi_single_dimension * 2, input_argument.rtree_max_nodes, input_argument.representation_option);
//		/*------------------------------------------------------------------------------------------------------------------------------------------------*/
//
//		/*---------------------------------------Original version Rtree insertion--------------------------------------------------*/
//		APCA_QUAL::initial_rect_vector(input_argument.point_number, input_argument.point_multi_single_dimension * 2, rtree_rectangle_vector);
//		/*-------------------------------------------------------------------------------------------------------------------------*/
//		//APCA_QUAL::initialAPCAArray(input_argument.point_number, input_argument.point_multi_single_dimension, apca_array);
//		//191111 APLA vector
//		//apla_array_vector.resize(input_argument.point_number);
//		apla_array_vector.resize(input_argument.point_number, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>());
//		/*--------------------------------------------------------APCA version--------------------------------------------------*/
//		//APCA_QUAL::initialAPCAArray(input_argument.point_number, input_argument.point_multi_single_dimension * 2, apca_MBR);
//		/*-----------------------------------------------------------------------------------------------------------------------*/
//		break;
//	}
//	case 2: {//PLA
//		//cout << "Begin PLA representation: " << endl;
//		RTree.initialRTree(input_argument.point_multi_single_dimension * 2, input_argument.rtree_max_nodes, input_argument.representation_option);
//		PLA_QUAL::initialPLAArray(input_argument.point_number, input_argument.point_multi_single_dimension, pla_array);
//		break;
//	}
//	case 3://APCA
//	case 7://PAALM
//	case 9://SAX
//	case 4: {//PAA
//		//cout << "Begin PAA & APCA representation: " << endl;
//		RTree.initialRTree(input_argument.point_multi_single_dimension * 2, input_argument.rtree_max_nodes, input_argument.representation_option);
//		APCA_QUAL::initialAPCAArray(input_argument.point_number, input_argument.point_multi_single_dimension, apca_point_array);
//		//APCA_QUAL::initialAPCAArray(input_argument.point_number, input_argument.point_multi_single_dimension, apca_array);
//		//APCA_QUAL::initialAPCAArray(input_argument.point_number, input_argument.point_multi_single_dimension * 2, apca_MBR);
//		//APCA_QUAL::initial_rect_vector(input_argument.point_number, input_argument.point_multi_single_dimension * 2, rtree_rectangle_vector);
//		break;
//	}
//	case 5: {//Cheby
//		//cout << "Begin CHEBY representation: " << endl;
//		RTree.initialRTree(input_argument.point_multi_single_dimension, input_argument.rtree_max_nodes, input_argument.representation_option);
//		CHEBYSHEV_QUAL::initialCHEBYSHEVArray(input_argument.point_number, input_argument.time_series_length, input_argument.point_multi_single_dimension, chebyshev_array);
//		CHEBYSHEV_QUAL::initialCHEBYSHEV_SHARE(input_argument, chebyshev_share);
//		CHEBYSHEV_QUAL::getCHEBYSHEV_SHARE(input_argument, chebyshev_share);
//		break;
//	}
//	default: 
//		assert(0);
//		break;
//	}
//	TOOL::printInputArgument(input_argument);
//}
//
//TEMPLATE
//void MULTI::deleteMultiDimension() {
//	/*for (int dimension_id = 0; dimension_id < input_argument.arity_d; dimension_id++) {
//		delete[] multi_dim.multi_dimension[dimension_id];
//		multi_dim.multi_dimension[dimension_id] = nullptr;
//	}
//
//	delete[] multi_dim.multi_dimension;
//	multi_dim.multi_dimension = nullptr;
//
//	delete[] multi_dim.projected_dimension;
//	multi_dim.projected_dimension = nullptr;*/
//
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//	case 6: {//ICDE07
//		//cout << "Finish APLA ICDE07 representation." << endl;
//		//APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_array);
//		/*--------------------------------*/
//		APCA_QUAL::delete_rect_vector(rtree_rectangle_vector);
//		for (auto&& au : apla_array_vector) {
//			au.clear();
//		}
//		apla_array_vector.clear();
//		apla_array_vector.shrink_to_fit();
//		/*--------------------------------*/
//		APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_MBR);
//		cout << endl;
//		break;
//	}
//	case 8: {//Initial 200706
//		//cout << "Finish APLA ICDE07 representation." << endl;
//		//APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_array);
//		/*--------------------------------*/
//		APCA_QUAL::delete_rect_vector(rtree_rectangle_vector);
//		for (auto&& au : apla_array_vector) {
//			for (auto auu : au) {
//				if (auu.right_subsegment != nullptr) {
//					delete auu.right_subsegment;
//					auu.right_subsegment = nullptr;
//				}
//			}
//			au.clear();
//		}
//		apla_array_vector.clear();
//		apla_array_vector.shrink_to_fit();
//		/*--------------------------------*/
//		APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_MBR);
//		cout << endl;
//		break;
//	}
//	case 2: {
//		//cout << "Finish PLA representation." << endl;
//		PLA_QUAL::deletePLAArray(input_argument, pla_array);
//		cout << endl;
//		break;
//	}
//	case 3://APCA
//	case 7://APCA
//	case 9://SAX
//	case 4: {//PAA
//		//cout << "Finish PAA & APCA representation." << endl;
//		APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_point_array);
//		/*APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_array);
//		APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_MBR);
//		APCA_QUAL::delete_rect_vector(rtree_rectangle_vector);*/
//		break;
//	}
//	case 5: {//CHEBY
//		//cout << "Finish CHEBY representation." << endl;
//		CHEBYSHEV_QUAL::deleteCHEBYSHEV_SHARE(chebyshev_share);
//		CHEBYSHEV_QUAL::deleteCHEBYSHEVArray(input_argument.point_number, chebyshev_array);
//		cout << endl;
//		break;
//	}
//	default:
//		assert(0);
//		break;
//	}
//}
//
////************************************
//// Method:  get_cmin_cmax
//// Qualifier: Pget MBR : min id, min value, max id, max value to insert into Rtree.
//// Input:
//// Output:
//// Notice: 
//// date: 191121
//// author:
////************************************
//TEMPLATE
//template<typename T>
//void MULTI::get_cmin_cmax(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, typename APCA_QUAL::APCA& const CminParameter, typename APCA_QUAL::APCA& const CmaxParameter) {
//#ifdef _DEBUG
//	assert(time_series_length != INF && linked_list.size() == CmaxParameter.segmentNum);
//	//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//
//	for (int segment_id = 0; segment_id < linked_list.size(); segment_id++) {
//		auto& const segment = linked_list[segment_id];
//#ifdef _DEBUG
//		// assert right ednpoint & min & max point
//		APLA::assert_segment_minmax(segment);
//		assert(time_series_length != INF && segment.right_endpoint != INF && segment.rectangle_width != INF && segment.min_point.id != INF && segment.max_point.id != INF && segment.min_point.value != INF && segment.max_point.value != INF);
//		//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//		CmaxParameter.r[segment_id] = segment.right_endpoint;
//		CminParameter.r[segment_id] = segment.right_endpoint - segment.rectangle_width;//191117 differentr from APCA paper: min.r == max.r
//		CminParameter.v[segment_id] = segment.min_point.value;
//		CmaxParameter.v[segment_id] = segment.max_point.value;
//		/*--------------------For high dimension dataset--------------------*/
//		/*
//		  so far most data is one dimension time series
//		*/
//		if (dimension_id > 0) {
//			int extend_right_endpoint = time_series_length * dimension_id;
//			segment.right_endpoint += extend_right_endpoint;
//			CminParameter.r[segment_id] += extend_right_endpoint;
//			CmaxParameter.r[segment_id] = segment.right_endpoint;
//		}
//		/*-----------------------------------------------------------------*/
//	}
//}
//
//
////************************************
//// Method:  get_cmin_cmax
//// Qualifier: get MBR : min id, min value, max id, max value to insert into Rtree.
//// Input: min id < max id
//// Output:
//// Notice: 
//// date: 191121
//// author:
////************************************
//TEMPLATE
//template<typename T, typename Y>
//void MULTI::get_cmin_cmax_original(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle) {
//
//#ifdef _DEBUG
//	assert(time_series_length == linked_list.back().right_endpoint + 1);
//	//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//
//	rtree_rectangle.m_min[0] = 0; // even min id
//	rtree_rectangle.m_min[1] = linked_list.front().min_point.value;// odd first is min value
//	rtree_rectangle.m_max[0] = linked_list.front().right_endpoint; // even max id
//	rtree_rectangle.m_max[1] = linked_list.front().max_point.value;// odd second is max value
//
//	for (int segment_id = 1; segment_id < linked_list.size(); segment_id++) {
//		auto& const segment = linked_list[segment_id];
//#ifdef _DEBUG
//		// assert right ednpoint & min & max point
//		APLA::assert_segment_minmax(segment);
//		assert(time_series_length != INF && segment.right_endpoint != INF && segment.rectangle_width != INF && segment.min_point.id != INF && segment.max_point.id != INF && segment.min_point.value != INF && segment.max_point.value != INF);
//		//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//		rtree_rectangle.m_min[segment_id * 2] = segment.right_endpoint - segment.rectangle_width + 1;/// even min id
//		rtree_rectangle.m_min[segment_id * 2 + 1] = segment.min_point.value;// odd first is min value
//		rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint; //even max id
//		rtree_rectangle.m_max[segment_id * 2 + 1] = segment.max_point.value;// odd second is max value
//
//#ifdef _DEBUG
////cout << rtree_rectangle.m_min[segment_id * 2] << " " << rtree_rectangle.m_min[segment_id * 2 + 1] << " " << rtree_rectangle.m_max[segment_id * 2] << " " << rtree_rectangle.m_max[segment_id * 2 + 1] << endl;
//		assert(rtree_rectangle.m_min[segment_id * 2] <= rtree_rectangle.m_max[segment_id * 2] && rtree_rectangle.m_min[segment_id * 2 + 1] <= rtree_rectangle.m_max[segment_id * 2 + 1]);
//#endif
//		/*--------------------For high dimension dataset--------------------*/
//		/*so far most data is one dimension time series*/
//		if (dimension_id > 0) {
//			int extend_right_endpoint = time_series_length * dimension_id;
//			segment.right_endpoint += extend_right_endpoint;
//			rtree_rectangle.m_min[segment_id * 2] += extend_right_endpoint;
//			rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;
//		}
//		/*-----------------------------------------------------------------*/
//	}
//}
//
////************************************
//// Method:  get_cmin_cmax_apca
//// Qualifier: get MBR : min id, min value, max id, max value to insert into Rtree.
//// Input: min_id == max id. min value, max value to insert into Rtree.
//// Output:
//// Notice: 
//// date: 191127
//// author:
////************************************ 
//TEMPLATE
//template<typename T, typename Y>
//void MULTI::get_cmin_cmax_apca(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle) {
//
//#ifdef _DEBUG
//	assert(time_series_length == linked_list.back().right_endpoint + 1);
//	//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//
//	for (int segment_id = 0; segment_id < linked_list.size(); segment_id++) {
//		auto& const segment = linked_list[segment_id];
//#ifdef _DEBUG
//		// assert right ednpoint & min & max point
//		APLA::assert_segment_minmax(segment);
//		assert(time_series_length != INF && segment.right_endpoint != INF && segment.rectangle_width != INF && segment.min_point.id != INF && segment.max_point.id != INF && segment.min_point.value != INF && segment.max_point.value != INF && segment.min_point.value <= segment.max_point.value);
//		//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//		rtree_rectangle.m_min[segment_id * 2] = rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;/// even min id //even max id
//		rtree_rectangle.m_min[segment_id * 2 + 1] = segment.min_point.value;// odd first is min value
//		rtree_rectangle.m_max[segment_id * 2 + 1] = segment.max_point.value;// odd second is max value
//
//#ifdef _DEBUG
////cout << rtree_rectangle.m_min[segment_id * 2] << " " << rtree_rectangle.m_min[segment_id * 2 + 1] << " " << rtree_rectangle.m_max[segment_id * 2] << " " << rtree_rectangle.m_max[segment_id * 2 + 1] << endl;
//		assert(rtree_rectangle.m_min[segment_id * 2] <= rtree_rectangle.m_max[segment_id * 2] && rtree_rectangle.m_min[segment_id * 2 + 1] <= rtree_rectangle.m_max[segment_id * 2 + 1]);
//#endif
//		/*--------------------For high dimension dataset--------------------*/
//		/*so far most data is one dimension time series*/
//		if (dimension_id > 0) {
//			int extend_right_endpoint = time_series_length * dimension_id;
//			segment.right_endpoint += extend_right_endpoint;
//			rtree_rectangle.m_min[segment_id * 2] += extend_right_endpoint;
//			rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;
//		}
//		/*-----------------------------------------------------------------*/
//	}
//}
//
////************************************
//// Method:  get_cmin_cmax_apca
//// Qualifier: get MBR : min_id == max id. min value, max value to insert into Rtree.
//// Input: min_id == max id. min value, max value to insert into Rtree.
//// Output:
//// Notice: 
//// date: 191128
//// author:
////************************************ 
//TEMPLATE
//template<typename T, typename Y, typename U>
//void MULTI::get_cmin_cmax_apca(const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const linked_list, U& const rtree_rectangle) {
//
//#ifdef _DEBUG
//	assert(original_time_series_vector.size() == linked_list.back().right_endpoint + 1);
//	//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//
//	for (int segment_id = 0; segment_id < linked_list.size(); segment_id++) {
//		const auto& const segment = linked_list[segment_id];
//
//#ifdef _DEBUG
//		// assert right ednpoint & min & max point
//		//APLA::assert_segment_minmax(segment);
//		assert(segment.right_endpoint != INF && segment.rectangle_width != INF);//&& segment.min_point.id != INF && segment.max_point.id != INF && segment.min_point.value != INF && segment.max_point.value != INF && segment.min_point.value <= segment.max_point.value);
//		//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//
//		/*.......................           200314    get min&max value of segment   ....................*/
//		const int id_coefficient = segment.right_endpoint + 1;
//		const auto& const value_minmax = minmax_element(original_time_series_vector.begin() + int(id_coefficient - segment.rectangle_width), original_time_series_vector.begin() + id_coefficient);
//		/*...............................................................................................*/
//
//		rtree_rectangle.m_min[segment_id * 2] = rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;// even min id //even max id
//		rtree_rectangle.m_min[segment_id * 2 + 1] = *value_minmax.first;// odd first is min value
//		rtree_rectangle.m_max[segment_id * 2 + 1] = *value_minmax.second;// odd second is max value
//		//rtree_rectangle.m_min[segment_id * 2 + 1] = segment.min_point.value;// odd first is min value
//		//rtree_rectangle.m_max[segment_id * 2 + 1] = segment.max_point.value;// odd second is max value
//
//#ifdef _DEBUG
////cout << rtree_rectangle.m_min[segment_id * 2] << " " << rtree_rectangle.m_min[segment_id * 2 + 1] << " " << rtree_rectangle.m_max[segment_id * 2] << " " << rtree_rectangle.m_max[segment_id * 2 + 1] << endl;
//		assert(rtree_rectangle.m_min[segment_id * 2] <= rtree_rectangle.m_max[segment_id * 2] && rtree_rectangle.m_min[segment_id * 2 + 1] <= rtree_rectangle.m_max[segment_id * 2 + 1]);
//#endif
//
//	}
//}
//
////************************************
//// Method:  get_cmin_cmax_pla
//// Qualifier: get MBR : min_id == max id. min value, max value to insert into Rtree.
//// Input: min_id == max id. min value, max value to insert into Rtree.
//// Output:get MBR : min_id == max id. PLA line into Rtree.
//// Notice: 
//// date: 191129
//// author:
////************************************ 
//TEMPLATE
//template<typename T, typename Y>
//void MULTI::get_cmin_cmax_pla(const int& const time_series_length, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle) {
//#ifdef _DEBUG
//	assert(time_series_length == linked_list.back().right_endpoint + 1);
//	//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//
//	for (int segment_id = 0; segment_id < linked_list.size(); segment_id++) {
//		auto& const segment = linked_list[segment_id];
//#ifdef _DEBUG
//		// assert right ednpoint & min & max point
//		//APLA::assert_segment_minmax(segment);
//		//assert(time_series_length != INF && segment.right_endpoint != INF && segment.rectangle_width != INF && segment.min_point.id != INF && segment.max_point.id != INF && segment.min_point.value != INF && segment.max_point.value != INF && segment.min_point.value <= segment.max_point.value);
//		assert(segment.apla.a != INF && segment.apla.b != INF && segment.right_endpoint != INF && segment.rectangle_width != INF);
//		//cout << segment.min_point.id << " " << segment.min_point.value << " " << segment.max_point.id << " " << segment.max_point.value << " " << segment.right_endpoint << endl;
//#endif
//		auto result = minmax(segment.apla.b, segment.apla.a * (segment.rectangle_width - 1) + segment.apla.b);
//		rtree_rectangle.m_min[segment_id * 2] = rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;/// even min id //even max id
//		rtree_rectangle.m_min[segment_id * 2 + 1] = result.first;// odd first is min value
//		rtree_rectangle.m_max[segment_id * 2 + 1] = result.second;// odd second is max value
//
//#ifdef _DEBUG
////cout << rtree_rectangle.m_min[segment_id * 2] << " " << rtree_rectangle.m_min[segment_id * 2 + 1] << " " << rtree_rectangle.m_max[segment_id * 2] << " " << rtree_rectangle.m_max[segment_id * 2 + 1] << endl;
//		assert(result.first <= result.second);
//		assert(rtree_rectangle.m_min[segment_id * 2] <= rtree_rectangle.m_max[segment_id * 2] && rtree_rectangle.m_min[segment_id * 2 + 1] <= rtree_rectangle.m_max[segment_id * 2 + 1]);
//#endif
//
//	}
//
//}
//
//TEMPLATE
//template<typename T>
//double& MULTI::getPLAMBRSegmentDistance(const typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const int& segment_id, const T& pla_array_qeury, double& pla_MBR_segment_distance) {
//	//assert(MBR.segmentNum == 2 * pla_array_qeury.segmentNum);
//
//	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
//
//
//	typename PLA_QUAL::POINT point;
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
//			return pla_MBR_segment_distance = PLA_QUAL::getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_2, point.v_B_2);
//		}
//		else if (point.u_B_2 >= point.u_A_2) {//case 1.1-1.2
//			//cout << "b" << endl;
//			return pla_MBR_segment_distance = PLA_QUAL::getNearestDistance(point.u_B_1, point.v_B_1, point.u_B_2, point.v_B_2, point.u_A_2, point.v_A_2);
//		}
//		else assert(0);
//	}
//	else if (point.u_A_1 <= 0 && point.u_A_2 >= 0) {//case 2: Line segment A1A2 is partially contained in the first and third quadrants of the u - v space
//		if (point.u_B_1 > 0) {
//			if (point.u_B_C_1 >= point.u_A_2) {//case 2.1; uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 > 0, uC ¡Ý uA2
//				//cout << "c" << endl;
//				return pla_MBR_segment_distance = PLA_QUAL::getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);//|A2B1|^2 // case2.1
//			}
//			else if (point.u_B_C_1 < point.u_A_2) {//case 2.2; uA1 ¡Ü 0, uA2 ¡Ý 0, uB1 > 0, uC < uA2
//				//cout << "d" << endl;
//				return pla_MBR_segment_distance = PLA_QUAL::getPointDistanceSquare(point.u_B_1, point.v_B_1, point.u_B_C_1, point.v_B_C_1);//|B1C|^2 //case 2.2
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
//				return pla_MBR_segment_distance = PLA_QUAL::getPointDistanceSquare(point.u_B_2, point.v_B_2, point.u_B_C_2, point.v_B_C_2);;//|B2C|^2 //case 2.4
//			}
//			else if (point.u_B_C_2 <= point.u_A_1) {//case: 2.5 uA1 ¡Ü 0, uA2 ¡Ý 0, uB2 < 0, uC ¡Ü uA1
//				//cout << "g" << endl;
//				return pla_MBR_segment_distance = PLA_QUAL::getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_B_2, point.v_B_2);//|A1B1|^2 //case 2.5
//			}
//			else assert(0);
//		}
//		else assert(0);
//	}
//	else if (point.u_A_1 > 0) {//case 3: Line segment A1A2 is completely contained in the first quadrant of the u-v space.
//		if (point.u_B_1 > point.u_A_1) {//case 3.1-3.3
//			//cout << sqrt(getPointDistanceSquare(point.u_B_2, 0.0, point.u_B_C_2, point.v_B_C_2)) << endl;
//			return pla_MBR_segment_distance = PLA_QUAL::getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);
//		}
//		else if (point.u_B_1 <= point.u_A_1) {//case 3.4-3.5
//			//cout << "i" << endl;
//			return pla_MBR_segment_distance = PLA_QUAL::getNearestDistance(point.u_B_1, point.v_B_1, point.u_B_2, point.v_B_2, point.u_A_1, point.v_A_1);
//		}
//		else assert(0);
//	}
//	else assert(0);
//
//	assert(0);
//}
//
//TEMPLATE
//template<typename T>
//double& MULTI::getPLAMBRDistance(typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const T& pla_array_qeury, double& pla_MBR_distance) {
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
////************************************
//// Method:  projectMultiToSingle
//// Qualifier: PAA APCA ,PLA Chebyshev
//// date: 180916
//// author:
////************************************
//TEMPLATE
//void MULTI::projectMultiToSingle() {
//	assert(input_argument.point_number != NULL && input_argument.read_multiple_file_name != nullptr && input_argument.arity_d != NULL);
//	input_argument.representation_time = 0;
//	DataType* original_time_series = new DataType[input_argument.time_series_length];
//
//	string row_string;
//	string row_number;
//	ifstream file_stream;
//
//	//from multi files
//	for (int file_id = 0; file_id < input_argument.arity_d; file_id++) {
//		file_stream = ifstream(input_argument.read_multiple_file_name[file_id]);
//		assert(file_stream);
//		cout << "Successfull read file: " << input_argument.read_multiple_file_name[file_id] << endl;
//
//		for (int point_id = 0; point_id < input_argument.point_number && (!file_stream.eof()) && file_stream.is_open() && file_stream.good(); point_id++) {
//			fill_n(original_time_series, input_argument.time_series_length, NULL);
//
//			file_stream >> row_string;
//			//memory_account[2] = row_string.size();
//			stringstream sstr(row_string);
//
//			int string_id = -1;
//			while (getline(sstr, row_number, ',') && string_id < input_argument.time_series_length) {
//				if (string_id > -1) {
//					original_time_series[string_id] = stod(row_number);
//					assert(original_time_series[string_id] != NULL);
//					//multi_dimention_trajectories[array_id][string_id] = stod(fs_row_number);
//					/*cout << original_time_series[string_id] << ", ";*/
//				}
//				string_id++;
//			}
//			/*cout << endl;*/
//
//			TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);//z-score normalization
//
//			if (input_argument.representation_option == 1) {//PAA
//				typename APCA_QUAL::APCA apca;
//				typename APCA_QUAL::APCA CminParameter, CmaxParameter, MBRParameter; //Temp MBR
//
//				APCA_QUAL::initialAPCA(apca, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(CminParameter, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(CmaxParameter, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(MBRParameter, input_argument.point_dimension * 2);
//
//				TOOL::recordStartTime(TOOL::time_record[0]);
//				APCA_QUAL::divideRemainderPAA(original_time_series, apca, input_argument.time_series_length, input_argument.point_dimension);
//				input_argument.representation_time += TOOL::recordFinishTime(TOOL::time_record[0]);
//
//				APCA_QUAL::getCmin(original_time_series, apca, CminParameter);
//				APCA_QUAL::getCmax(original_time_series, apca, CmaxParameter);
//
//				if (point_id == 0 && input_argument.arity_d == 1) {
//					if (input_argument.representation_option == 1)
//						TOOL::writeArray("Original time series", original_time_series, input_argument.time_series_length);
//
//					cout << "Original time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					APCA_QUAL::approximateOriginalFunction(apca, original_time_series);
//					cout << "Approximated PAA time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					TOOL::writeArray("PAA Approximation", original_time_series, input_argument.time_series_length);
//				}
//
//				for (int array_id = 0; array_id < apca.segmentNum; array_id++) {
//					apca.r[array_id] += input_argument.time_series_length * file_id;
//					CminParameter.r[array_id] += input_argument.time_series_length * file_id;
//					CmaxParameter.r[array_id] += input_argument.time_series_length * file_id;
//				}
//
//				APCA_QUAL::getMBR(CminParameter, CmaxParameter, MBRParameter);
//
//				copy_n(apca.r, apca.segmentNum, apca_array[point_id].APCALink.r + input_argument.point_dimension * file_id);
//				copy_n(apca.v, apca.segmentNum, apca_array[point_id].APCALink.v + input_argument.point_dimension * file_id);
//				copy_n(MBRParameter.r, MBRParameter.segmentNum, apca_MBR[point_id].r + MBRParameter.segmentNum * file_id);
//				copy_n(MBRParameter.v, MBRParameter.segmentNum, apca_MBR[point_id].v + MBRParameter.segmentNum * file_id);
//
//				/*TOOL::printArray(apca_MBR[point_id].r, apca_MBR[point_id].segmentNum);
//				TOOL::printArray(apca_MBR[point_id].v, apca_MBR[point_id].segmentNum);*/
//
//				APCA_QUAL::deleteAPCA(CminParameter);
//				APCA_QUAL::deleteAPCA(CmaxParameter);
//				APCA_QUAL::deleteAPCA(MBRParameter);
//				APCA_QUAL::deleteAPCA(apca);
//			}
//			else if (input_argument.representation_option == 2) {//APCA
//				typename APCA_QUAL::APCA apca;
//
//				typename APCA_QUAL::APCA CminParameter, CmaxParameter, MBRParameter; //Temp MBR
//				APCA_QUAL::initialAPCA(apca, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(CminParameter, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(CmaxParameter, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(MBRParameter, input_argument.point_dimension * 2);
//
//				//APCA_QUAL::divideRemainderPAA(original_time_series, apca, input_argument.time_series_length, input_argument.point_dimension);
//				TOOL::recordStartTime(TOOL::time_record[0]);
//				APCA_QUAL::getAPCAPoint(original_time_series, input_argument.time_series_length, input_argument.point_dimension, apca);
//				input_argument.representation_time += TOOL::recordFinishTime(TOOL::time_record[0]);
//
//				APCA_QUAL::getCmin(original_time_series, apca, CminParameter);
//				APCA_QUAL::getCmax(original_time_series, apca, CmaxParameter);
//
//				for (int array_id = 0; array_id < apca.segmentNum; array_id++) {
//					apca.r[array_id] += input_argument.time_series_length * file_id;
//					CminParameter.r[array_id] += input_argument.time_series_length * file_id;
//					CmaxParameter.r[array_id] += input_argument.time_series_length * file_id;
//				}
//
//				APCA_QUAL::getMBR(CminParameter, CmaxParameter, MBRParameter);
//
//				if (point_id == 0 && input_argument.arity_d == 1) {
//					cout << "Original time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					APCA_QUAL::approximateOriginalFunction(apca, original_time_series);
//					cout << "Approximated APCA time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					TOOL::writeArray("APCA Approximation", original_time_series, input_argument.time_series_length);
//				}
//
//				copy_n(apca.r, apca.segmentNum, apca_array[point_id].APCALink.r + input_argument.point_dimension * file_id);
//				copy_n(apca.v, apca.segmentNum, apca_array[point_id].APCALink.v + input_argument.point_dimension * file_id);
//				copy_n(MBRParameter.r, MBRParameter.segmentNum, apca_MBR[point_id].r + MBRParameter.segmentNum * file_id);
//				copy_n(MBRParameter.v, MBRParameter.segmentNum, apca_MBR[point_id].v + MBRParameter.segmentNum * file_id);
//
//				/*TOOL::printArray(apca_MBR[point_id].r, apca_MBR[point_id].segmentNum);
//				TOOL::printArray(apca_MBR[point_id].v, apca_MBR[point_id].segmentNum);*/
//
//				APCA_QUAL::deleteAPCA(CminParameter);
//				APCA_QUAL::deleteAPCA(CmaxParameter);
//				APCA_QUAL::deleteAPCA(MBRParameter);
//				APCA_QUAL::deleteAPCA(apca);
//			}
//			else if (input_argument.representation_option == 3) {//PLA
//				typename PLA_QUAL::PLA pla;
//				PLA_QUAL::initialPLA(pla, input_argument.point_dimension);
//
//				TOOL::recordStartTime(TOOL::time_record[0]);
//				PLA_QUAL::getPLA(input_argument, original_time_series, pla);
//
//				input_argument.representation_time += TOOL::recordFinishTime(TOOL::time_record[0]);
//
//				copy_n(pla.a, pla.segmentNum, pla_array[point_id].a + input_argument.point_dimension * file_id);
//				copy_n(pla.b, pla.segmentNum, pla_array[point_id].b + input_argument.point_dimension * file_id);
//
//				if (point_id == 0 && input_argument.arity_d == 1) {
//					cout << "Original time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					PLA_QUAL::approximateOriginalFunctionPLA1(input_argument, pla, original_time_series);
//					cout << "Approximated PLA time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					TOOL::writeArray("PLA Approximation", original_time_series, input_argument.time_series_length);
//				}
//
//				PLA_QUAL::deletePLA(pla);
//			}
//			else if (input_argument.representation_option == 4) {//CHEBY
//				if (point_id == 0 && input_argument.arity_d == 1) {
//					DataType* approximated_time_series = new DataType[input_argument.time_series_length];
//
//					cout << "Original time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					CHEBYSHEV_QUAL::approximateOriginalFunction(input_argument, original_time_series, chebyshev_share, approximated_time_series);
//					cout << "Approximated CHEBY time series: " << endl;
//					TOOL::printArray(approximated_time_series, input_argument.time_series_length);
//					TOOL::writeArray("CHEBY Approximation", approximated_time_series, input_argument.time_series_length);
//
//					delete[] approximated_time_series;
//					approximated_time_series = nullptr;
//				}
//
//				typename CHEBYSHEV_QUAL::CHEBYSHEV chebyshev;
//				CHEBYSHEV_QUAL::initialCHEBYSHEV(input_argument, chebyshev);
//
//				TOOL::recordStartTime(TOOL::time_record[0]);
//				CHEBYSHEV_QUAL::getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev);
//				input_argument.representation_time += TOOL::recordFinishTime(TOOL::time_record[0]);
//
//				copy_n(chebyshev.coefficient, chebyshev.segmentNum, chebyshev_array[point_id].coefficient + input_argument.point_dimension * file_id);
//
//				CHEBYSHEV_QUAL::deleteCHEBYSHEV(chebyshev);
//			}
//			else {
//				assert(0);
//			}
//		}
//	}
//
//	cout << "representation time: " << input_argument.representation_time << " us" << endl;
//	file_stream.close();
//	row_string.clear();
//	row_string.shrink_to_fit();
//	row_number.clear();
//	row_number.shrink_to_fit();
//	/*cout << "Root Node : sub node number = " << RTree.m_root->m_count << " Root level = : " << RTree.m_root->m_level << "\n\nBegin to build a RTree:\n";
//	cout << "\n RTree conclusion\n The number of RTree Data Point = : " << RTree.Count() << endl;*/
//
//	//deletePLA(pla_array_MBR);
//	//deleteCHEBYSHEV(chebyshev_MBR);
//	//deleteMultiDimensionTrajectory(input_argument, multi_dimention_trajectories);
//
//	delete[] original_time_series;
//	original_time_series = nullptr;
//}
//
////***************************************************************
//// Method:project_multi_to_single
//// Qualifier: optimization of above projectMultiToSingle()
//// Input:
//// Output: Approximated time series
//// Note: for differnt approximate input_argument, the N should be different from multi class input_argument.  APLA = 4, PAA : 12 , PLA: 6, APCA 6. Cheby 12, ICDE07: 4
//// date:191106
//// author:
////***************************************************************
//TEMPLATE
//template<typename T, typename Y>
//void MULTI::project_multi_to_single(vector<Y>& const multi_y_projection_argument, vector<DoublyLinkedList<T>>& const multi_all_linked_list, vector<DoublyLinkedList<T>>& const multi_cluster_linked_list) {
//#ifdef _DEBUG
//	assert(input_argument.point_number != NULL && input_argument.file_id != INF && input_argument.arity_d != NULL && multi_y_projection_argument.size() == input_argument.point_number && multi_all_linked_list.size() == input_argument.point_number && multi_all_linked_list.size() == multi_cluster_linked_list.size());
//	input_argument.representation_time = 0;
//#endif
//
//	DataType* original_time_series = new DataType[input_argument.time_series_length];
//	vector<DataType> normalized_series_vector(input_argument.time_series_length, INF);
//
//	ifstream file_stream;// for string reading
//	string row_string;
//	string row_number;
//
//	//from multi files
//	for (int dimension_id = 0; dimension_id < input_argument.arity_d; dimension_id++) {//data dimension
//		/*======================================= Read file from dataset(txt file ) =================================================*/
//		if (input_argument.read_multiple_file_name) {
//			file_stream = ifstream(input_argument.read_multiple_file_name[dimension_id]);
//		}
//		else {
//			file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id));
//		}
//		assert(file_stream);
//
//		//cout << "Successfull read file: " << input_argument.read_multiple_file_name[dimension_id] << endl;
//		/*============================================================================================================================*/
//		// point_id := time series id
//		for (int point_id = 0; point_id < input_argument.point_number && (!file_stream.eof()) && file_stream.is_open() && file_stream.good(); point_id++) {
//#ifdef _DEBUG
//			//cout << "point id: " << point_id << endl;
//			input_argument.point_id = point_id;
//#endif
//			fill_n(original_time_series, input_argument.time_series_length, INF);
//
//			/*=====================================Read file from dataset(txt file)==================================================*/
//			TOOL::get_normalized_series_by_stream(input_argument, file_stream, original_time_series);
//			copy_n(original_time_series, input_argument.time_series_length, normalized_series_vector.begin());
//
//			//if (point_id == 0 && input_argument.arity_d == 1) {
//				/*if (input_argument.representation_option == 1)
//					TOOL::writeArray("Original time series", original_time_series, input_argument.time_series_length);
//
//				cout << "Original time series: " << endl;
//				TOOL::printArray(original_time_series, input_argument.time_series_length);
//				APCA_QUAL::approximateOriginalFunction(apca, original_time_series);
//				cout << "Approximated PAA time series: " << endl;
//				TOOL::printArray(original_time_series, input_argument.time_series_length);
//				TOOL::writeArray("PAA Approximation", original_time_series, input_argument.time_series_length);*/
//				/*for (auto au : normalized_series_vector) {
//					cout << au << ",";
//				}
//				cout << endl;*/
//				//}
//			/*===========================================================================================================================*/
//
//		/*=====================================================================Approximation===========================================================================================*/
//			input_argument.whole_run_time = 0;
//			TOOL::recordStartTime(TOOL::time_record[14]);
//
//			if (input_argument.representation_option == 1 || input_argument.representation_option == 2) {
//				/*================================= PAA & APCA =====================================================================*/
//#ifdef _DEBUG
//				assert(rtree_rectangle_vector.size() == input_argument.point_number);
//#endif
//				typename APCA_QUAL::APCA apca;
//				typename APCA_QUAL::APCA CminParameter, CmaxParameter, MBRParameter; //Temp MBR
//				APCA_QUAL::initialAPCA(apca, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(CminParameter, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(CmaxParameter, input_argument.point_dimension);
//				APCA_QUAL::initialAPCA(MBRParameter, input_argument.point_dimension * 2);
//
//				//APCA_QUAL::divideRemainderPAA(original_time_series, apca, input_argument.time_series_length, input_argument.point_dimension);
//#ifdef _DEBUG
//				TOOL::recordStartTime(TOOL::time_record[0]);
//#endif          
//				switch (input_argument.representation_option) {
//				case 1: //PAA
//					APCA_QUAL::divideRemainderPAA(original_time_series, apca, input_argument.time_series_length, input_argument.point_dimension);
//
//					/*----------------Evaluation: Use mean value to instead min&max value. getCmin() getCmax() -------------------------*/
//					//CminParameter.segmentNum = apca.segmentNum;
//					//CmaxParameter.segmentNum = apca.segmentNum;
//					//for (int i = 0; i < apca.segmentNum; i++) {
//					//	CminParameter.r[i] = apca.r[i];
//					//	CmaxParameter.r[i] = apca.r[i];
//					//	CminParameter.v[i] = apca.v[i];//min[begin, end);
//					//	CmaxParameter.v[i] = apca.v[i];
//					//	//cout << "r[" << i << "] = " << Cmin.r[i] << ", Cmin0[" << i << "] = " << Cmin.v[i] << endl << endl;
//					//}
//					/*-------------------------------------------------------------------------------------------------------  ---------*/
//
//					break;
//				case 2://APCA
//					APCA_QUAL::getAPCAPoint(original_time_series, input_argument.time_series_length, input_argument.point_dimension, apca);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//#ifdef _DEBUG
//				input_argument.representation_time += TOOL::recordFinishTime(TOOL::time_record[0]);
//#endif
//				APCA_QUAL::getCmin(original_time_series, apca, CminParameter);
//				APCA_QUAL::getCmax(original_time_series, apca, CmaxParameter);
//				/*---------------------------------191120 Original Rtree MBR-------------------------------------------*/
//				//RTREE::Rect rtree_rectangle;
//				//APCA_QUAL::get_minmax_original(normalized_series_vector, dimension_id, apca, rtree_rectangle_vector[point_id]);
//				/*--------------------------------------------------------------------------------------------------*/
//				/*---------------------------------191127 APCA Paper Rtree MBR-------------------------------------------*/
//				APCA_QUAL::get_minmax_apca(normalized_series_vector, dimension_id, apca, rtree_rectangle_vector[point_id]);
//				/*--------------------------------------------------------------------------------------------------*/
//#ifdef _DEBUG
//				/*for (int segment_id = 0; segment_id < apca.segmentNum; segment_id++ ) {
//					int point_end_id = segment_id << 1;
//					int value_id = (segment_id << 1) + 1;
//
//					int c_max_id = CmaxParameter.r[segment_id];
//					int r_max_id = rtree_rectangle_vector[point_id].m_max[point_end_id];
//					int c_min_value = CminParameter.v[segment_id];
//					int r_min_value = rtree_rectangle_vector[point_id].m_min[value_id];
//					int c_max_value = CmaxParameter.v[segment_id];
//					int r_max_value = rtree_rectangle_vector[point_id].m_max[value_id];
//
//					assert(rtree_rectangle_vector[point_id].m_min[point_end_id] <= r_max_id && c_max_id == r_max_id && c_min_value == r_min_value && c_max_value == r_max_value);
//				}*/
//#endif
//				for (int array_id = 0; array_id < apca.segmentNum; array_id++) {
//					apca.r[array_id] += input_argument.time_series_length * dimension_id;
//					CminParameter.r[array_id] += input_argument.time_series_length * dimension_id;
//					CmaxParameter.r[array_id] += input_argument.time_series_length * dimension_id;
//				}
//
//				APCA_QUAL::getMBR(CminParameter, CmaxParameter, MBRParameter);
//#ifdef _DEBUG
//				if (point_id == 0 && input_argument.arity_d == 1) {
//					/*	cout << "Original time series: " << endl;
//						TOOL::printArray(original_time_series, input_argument.time_series_length);
//						APCA_QUAL::approximateOriginalFunction(apca, original_time_series);
//						cout << "Approximated APCA time series: " << endl;
//						TOOL::printArray(original_time_series, input_argument.time_series_length);
//						TOOL::writeArray("APCA Approximation", original_time_series, input_argument.time_series_length);*/
//				}
//#endif
//				copy_n(apca.r, apca.segmentNum, apca_array[point_id].APCALink.r + input_argument.point_dimension * dimension_id);
//				copy_n(apca.v, apca.segmentNum, apca_array[point_id].APCALink.v + input_argument.point_dimension * dimension_id);
//				copy_n(MBRParameter.r, MBRParameter.segmentNum, apca_MBR[point_id].r + MBRParameter.segmentNum * dimension_id);
//				copy_n(MBRParameter.v, MBRParameter.segmentNum, apca_MBR[point_id].v + MBRParameter.segmentNum * dimension_id);
//
//				/*TOOL::printArray(apca_MBR[point_id].r, apca_MBR[point_id].segmentNum);
//				TOOL::printArray(apca_MBR[point_id].v, apca_MBR[point_id].segmentNum);*/
//
//				APCA_QUAL::deleteAPCA(CminParameter);
//				APCA_QUAL::deleteAPCA(CmaxParameter);
//				APCA_QUAL::deleteAPCA(MBRParameter);
//				APCA_QUAL::deleteAPCA(apca);
//			}
//			else if (input_argument.representation_option == 3) {//PLA
//				/*================================= PLA =====================================================================*/
//				typename PLA_QUAL::PLA pla;
//				PLA_QUAL::initialPLA(pla, input_argument.point_dimension);
//#ifdef _DEBUG
//				TOOL::recordStartTime(TOOL::time_record[0]);
//#endif
//				PLA_QUAL::getPLA(input_argument, original_time_series, pla);
//#ifdef _DEBUG
//				input_argument.representation_time += TOOL::recordFinishTime(TOOL::time_record[0]);
//#endif
//				copy_n(pla.a, pla.segmentNum, pla_array[point_id].a + input_argument.point_dimension * dimension_id);
//				copy_n(pla.b, pla.segmentNum, pla_array[point_id].b + input_argument.point_dimension * dimension_id);
//
//				//if (point_id == 0 && input_argument.arity_d == 1) {
//					/*cout << "Original time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					PLA_QUAL::approximateOriginalFunctionPLA1(input_argument, pla, original_time_series);
//					cout << "Approximated PLA time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					TOOL::writeArray("PLA Approximation", original_time_series, input_argument.time_series_length);*/
//					//}
//
//				PLA_QUAL::deletePLA(pla);
//			}
//			else if (input_argument.representation_option == 4) {//CHEBY
//			/*=================================  CHEBY  =====================================================================*/
//			//	if (point_id == 0 && input_argument.arity_d == 1) {
//					//DataType* approximated_time_series = new DataType[input_argument.time_series_length];
//
//				/*	cout << "Original time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					CHEBYSHEV_QUAL::approximateOriginalFunction(input_argument, original_time_series, chebyshev_share, approximated_time_series);
//					cout << "Approximated CHEBY time series: " << endl;
//					TOOL::printArray(approximated_time_series, input_argument.time_series_length);
//					TOOL::writeArray("CHEBY Approximation", approximated_time_series, input_argument.time_series_length);*/
//
//					//delete[] approximated_time_series;
//					//approximated_time_series = nullptr;
//			//	}
//
//				typename CHEBYSHEV_QUAL::CHEBYSHEV chebyshev;
//				CHEBYSHEV_QUAL::initialCHEBYSHEV(input_argument, chebyshev);
//#ifdef _DEBUG
//				TOOL::recordStartTime(TOOL::time_record[0]);
//#endif
//				CHEBYSHEV_QUAL::getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev);
//#ifdef _DEBUG
//				input_argument.representation_time += TOOL::recordFinishTime(TOOL::time_record[0]);
//#endif
//				copy_n(chebyshev.coefficient, chebyshev.segmentNum, chebyshev_array[point_id].coefficient + input_argument.point_dimension * dimension_id);
//
//				CHEBYSHEV_QUAL::deleteCHEBYSHEV(chebyshev);
//			}
//			else if (input_argument.representation_option == 5 || input_argument.representation_option == 6) {//APLA & ICDE07
//			/*======================================================== APLA & ICDE07 ==============================================================*/
//#ifdef _DEBUG
//				assert(input_argument.point_number == apla_array_vector.size() && rtree_rectangle_vector.size() == input_argument.point_number);
//				TOOL::recordStartTime(TOOL::time_record[0]);
//#endif
//				switch (input_argument.representation_option) {
//				case 5: {//APLA
//#ifdef _DEBUG
//					//cout << "APLA:" << endl;
//#endif
//					/*------------------------------------------Get APLA approximation----------------------------------------------*/
//					//DoublyLinkedList<APLA::AREA_COEFFICIENT> doubly_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();
//					//apla_array_vector[point_id] = DoublyLinkedList<APLA::AREA_COEFFICIENT>();
//					APLA::get_APLA_point(input_argument, original_time_series, multi_y_projection_argument[point_id], multi_all_linked_list[point_id], multi_cluster_linked_list[point_id], apla_array_vector[point_id]);
//					//multi_all_linked_list.clear();  
//					//multi_cluster_linked_list.clear();
//					/*--------------------------------------------------------------------------------------------------------------*/
//					break;
//				}
//				case 6: {//ICDE07
//#ifdef _DEBUG
//					//cout << "ICDE07:" << endl;
//#endif
//					//apla_array_vector[point_id].resize(input_argument.point_dimension);
//					APLA_ICDE07<DataType>::getAPLA_ICDE07(input_argument, normalized_series_vector, apla_array_vector[point_id]);
//					for (int segment_id = 0; segment_id < apla_array_vector[point_id].size(); segment_id++) {
//						APLA::getSegmentMinMaxPoint(original_time_series, apla_array_vector[point_id][segment_id]);
//					};
//					break;
//				}
//				default:
//					assert(0);
//					break;
//				}
//#ifdef _DEBUG
//				for (auto&& au : apla_array_vector[point_id]) {
//					APLA::assert_segment_minmax(au);
//				}
//				input_argument.representation_time += TOOL::recordFinishTime(TOOL::time_record[0]);
//				assert(input_argument.point_dimension == apla_array_vector[point_id].size());
//#endif
//				/*--------------------------------APCA version: Initial RTree insertion coefficient-------------------------------------------*/
//				//typename APCA_QUAL::APCA CminParameter, CmaxParameter, MBRParameter; //Temp MBR
//				//APCA_QUAL::initialAPCA(apca, input_argument.point_dimension);
//				//APCA_QUAL::initialAPCA(CminParameter, input_argument.point_dimension);
//				//APCA_QUAL::initialAPCA(CmaxParameter, input_argument.point_dimension);
//				//APCA_QUAL::initialAPCA(MBRParameter, input_argument.point_dimension * 2);
//				/*--------------------------------------------------------------------------------------------------------------*/
//				/*------------------------------------------Get Segment Min&Max point. min id < max id-------------------------------------------*/
//				//get_cmin_cmax_original(input_argument.time_series_length, dimension_id, apla_array_vector[point_id], rtree_rectangle_vector[point_id]);
//				/*-------------------------------------------------------------------------------------------------------------------------------*/
//				/*------------------------------------------Original version: Get Segment Min&Max point Cmin Cmax---------------------------------------------*/
//				//get_cmin_cmax(input_argument.time_series_length, dimension_id, apla_array_vector[point_id],  CminParameter, CmaxParameter);
//				/*--------------------------------------------------------------------------------------------------------------------------------------------*/
//				/*------------------------------------------191127 APCA version: Get Segment Min&Max point Cmin Cmax-------------------------------------------*/
//				get_cmin_cmax_apca(input_argument.time_series_length, dimension_id, apla_array_vector[point_id], rtree_rectangle_vector[point_id]);
//				/*---------------------------------------------------------------------------------------------------------------------------------------------*/
//				/*--------------------APCA version:  get MBR for Rtree min max point -------------------------------------------------------------*/
//				//APCA_QUAL::getMBR(CminParameter, CmaxParameter, MBRParameter);
//				/*-------------------------------------------------------------------------------------------------------------------------------*/
//				//copy_n(apca.r, apca.segmentNum, apca_array[point_id].APCALink.r + input_argument.point_dimension * file_id);
//				//copy_n(apca.v, apca.segmentNum, apca_array[point_id].APCALink.v + input_argument.point_dimension * file_id);.
//				/*--------------------APCA version:  get MBR for Rtree min max point -------------------------------------------------------------*/
//				//copy_n(MBRParameter.r, MBRParameter.segmentNum, apca_MBR[point_id].r + MBRParameter.segmentNum * dimension_id);
//				//copy_n(MBRParameter.v, MBRParameter.segmentNum, apca_MBR[point_id].v + MBRParameter.segmentNum * dimension_id);
//				/*-----------------------------------------------------------------------------------------------------------------------------------*/
//#ifdef _DEBUG
//				/*for (int segment_id = 0; segment_id < MBRParameter.segmentNum; segment_id++ ) {
//					cout << apca_MBR[point_id].r[segment_id] << " " << apca_MBR[point_id].r[segment_id] << endl;
//				}*/
//#endif
//				/*TOOL::printArray(apca_MBR[point_id].r, apca_MBR[point_id].segmentNum);
//				TOOL::printArray(apca_MBR[point_id].v, apca_MBR[point_id].segmentNum);*/
//
//				/*------------------------------------------Print function--------------------------------------------------------*/
//				/*if (point_id == 0 && input_argument.arity_d == 1) {
//					cout << "Original time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					APCA_QUAL::approximateOriginalFunction(apca, original_time_series);
//					cout << "Approximated APCA time series: " << endl;
//					TOOL::printArray(original_time_series, input_argument.time_series_length);
//					TOOL::writeArray("APCA Approximation", original_time_series, input_argument.time_series_length);
//				}*/
//				/*-----------------------------------------------------------------------------------------------------------------*/
//				/*-----------------------------------------Evaluate funciton---------------------------------------------------------*/
//#ifdef _DEBUG
//				assert(input_argument.point_number == apla_array_vector.size() && apla_array_vector[point_id].size() == input_argument.point_dimension);
//#endif
//				/*-----------------------------------------------------------------------------------------------------------------*/
//				/*-----------------------------------------APCA version: Clear funciton---------------------------------------------------------*/
//				//APCA_QUAL::deleteAPCA(CminParameter);
//				//APCA_QUAL::deleteAPCA(CmaxParameter);
//				//APCA_QUAL::deleteAPCA(MBRParameter);
//				/*-----------------------------------------------------------------------------------------------------------------*/
//
//			}
//			else {
//				assert(0);
//			}
//
//			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//			/*=============================================================================================================================================*/
//		}
//	}
//
//	cout << "representation time: " << input_argument.representation_time << " us" << endl;
//	file_stream.close();
//
//	/*cout << "Root Node : sub node number = " << RTree.m_root->m_count << " Root level = : " << RTree.m_root->m_level << "\n\nBegin to build a RTree:\n";
//	cout << "\n RTree conclusion\n The number of RTree Data Point = : " << RTree.Count() << endl;*/
//
//	//deletePLA(pla_array_MBR);
//	//deleteCHEBYSHEV(chebyshev_MBR);
//	//deleteMultiDimensionTrajectory(input_argument, multi_dimention_trajectories);
//
//	delete[] original_time_series;
//	original_time_series = nullptr;
//	normalized_series_vector.clear();
//	normalized_series_vector.shrink_to_fit();
//}
//
////***************************************************************
//// Method:approximation_build_rtree
//// Qualifier: optimization of above projectMultiToSingle()
//// Input:
//// Output: inserted rtree index
//// Note: 
//// date:191128
//// author:
////***************************************************************
//TEMPLATE
//template<typename T, typename Y>
//void MULTI::all_approximation_build_rtree(vector<T>& const multi_y_projection_argument, vector<DoublyLinkedList<Y>>& const multi_all_linked_list, vector<DoublyLinkedList<Y>>& const multi_cluster_linked_list) {
//
//	/*..........................................................................................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.point_number != NULL && input_argument.file_id != INF && input_argument.arity_d != NULL && multi_y_projection_argument.size() == input_argument.point_number && multi_all_linked_list.size() == input_argument.point_number && multi_all_linked_list.size() == multi_cluster_linked_list.size());
//#endif
//	/*..........................................................................................................................................................*/
//
//	/*-----------------   Time Evaluation   -------------------*/
//	input_argument.representation_time = 0;
//	input_argument.build_rtree_time = 0;
//	input_argument.whole_run_time = 0;
//	input_argument.whole_run_time_has_IO = 0;//210606
//	output_argument.sum_deviation = 0;
//	output_argument.max_deviation = 0;
//	output_argument.max_deviation_multiple_width = 0;
//	/*--------------------------------------------------------*/
//
//	DataType* original_time_series = new DataType[input_argument.time_series_length];
//	vector<DataType> normalized_series_vector(input_argument.time_series_length, INF);
//
//	//ifstream file_stream;// for string reading
//	string row_string;
//	string row_number;
//
//	//from multi files
//	//for (int dimension_id = 0; dimension_id < input_argument.arity_d; dimension_id++) {//data dimension
//
//	// point_id := time series id. The number of time series
//	for (int point_id = 0; point_id < input_argument.point_number; point_id++) {
//		/*..........................................*/
//#ifdef _DEBUG
//		//cout << "point id: " << point_id << endl;
//		input_argument.point_id = point_id;
//#endif
//		/*..........................................*/
//
//		fill_n(original_time_series, input_argument.time_series_length, INF);
//
//		/*========================================= Read dataset from file; get normalized data set ==========================================================*/
//		//already normalized and store in new file, no need skip first point
//		TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, point_id, normalized_series_vector);
//		copy_n(normalized_series_vector.begin(), normalized_series_vector.size(), original_time_series);
//
//		/*==============================================     200901 Print each time series results     ==================================================*/
//		/*if (input_argument.print_each_result == true) {
//			cout << point_id << " Normalized Time Series: \n";
//			for (auto&& au : normalized_series_vector) {
//				cout << au << ",";
//			}
//			cout << endl;
//		}*/
//		/*===============================================================================================================================================*/
//
//		/* if (data_source.data_type == 1 || data_source.data_type == 2) {
//			 TOOL::get_normalized_series_by_stream(input_argument, file_stream, original_time_series);
//			 copy_n(original_time_series, input_argument.time_series_length, normalized_series_vector.begin());
// #ifdef _DEBUG
//			 for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
//				 assert(normalized_series_vector[array_id] == original_time_series[array_id]);
//			 }
// #endif
//		 }*/
//		 /*===================================================================================================================================================*/
//
//		/*=====================================================================Approximation===========================================================================================*/
//
//	    /*-------------------------Time Evaluation---------------------------------*/
//		TOOL::recordStartTime(TOOL::time_record[14]);
//		TOOL::recordStartTime(TOOL::time_record[0]);//representation time
//		TOOL::recordStartTime(TOOL::time_record[16]);
//		/*-------------------------------------------------------------------------*/
//
//		/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  SAPLA & ICDE07  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//		if (input_argument.representation_option == 1 || input_argument.representation_option == 6 || input_argument.representation_option == 8) {//SAPLA & ICDE07 & initial200706
//		
//			/*...........................................................................................................................*/
//#ifdef _DEBUG
//			assert(input_argument.point_number == apla_array_vector.size() && rtree_rectangle_vector.size() == input_argument.point_number);
//#endif
//			/*...........................................................................................................................*/
//
//			typename TOOL::RECTANGLE rectangle_insertion(input_argument.point_dimension << 1);
//
//			switch (input_argument.representation_option) {
//			case 1: {//MSPLA
//
//				/*++++++++++++++++++++++++++++++++++++++++++            Get MSPLA approximation         +++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//				//APLA::get_APLA_point(input_argument, original_time_series, multi_y_projection_argument[point_id], multi_all_linked_list[point_id], multi_cluster_linked_list[point_id], apla_array_vector[point_id]);
//				//191209
//				//APLA::get_APLA_point_new_y_projection(input_argument, original_time_series, multi_y_projection_argument[point_id], multi_all_linked_list[point_id], multi_cluster_linked_list[point_id], apla_array_vector[point_id]);
//				//APLA::get_APLA_point_new_y_projection(input_argument, normalized_series_vector, multi_y_projection_argument[point_id], multi_all_linked_list[point_id], multi_cluster_linked_list[point_id], apla_array_vector[point_id], output_argument);
//				APLA::initial_SAPLA_200706(input_argument, normalized_series_vector, apla_array_vector[point_id], output_argument);
//				/*----------------------------------------200113 Test lower bound distance--------------------------*/
//				/*for (auto au : apla_array_vector[point_id]) {
//					APLA::get_segment_sum_deviation(normalized_series_vector, au);
//					assert(au.rec_deviation != INF);
//				}*/
//				/*---------------------------------------------------------------------------------------------------*/
//				//multi_all_linked_list.clear();  
//				//multi_cluster_linked_list.clear();
//				/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//
//				break;
//			}
//			case 6: {//ICDE07
//				/*++++++++++++++++++++++++++++++++++++++++++         Get ICDE07 approximation       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//				//apla_array_vector[point_id].resize(input_argument.point_dimension);
//				APLA_ICDE07<DataType>::getAPLA_ICDE07(input_argument, normalized_series_vector, apla_array_vector[point_id]);
//
//				/*for (int segment_id = 0; segment_id < apla_array_vector[point_id].size(); segment_id++) { 210203
//					APLA::getSegmentMinMaxPoint(normalized_series_vector, apla_array_vector[point_id][segment_id]);
//				}*/ 
//
//				/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//
//				break;
//			}
//			case 8: {//initial 200706
//				/*++++++++++++++++++++++++++++++++++++++++++         Get initial 200706 approximation       +++++++++++++++++++++++++++++++++++++++++++++++++++*/
//				APLA::initial_SAPLA_200706(input_argument, normalized_series_vector, apla_array_vector[point_id], output_argument);
//				//apla_array_vector[point_id].resize(input_argument.point_dimension);
//				/*APLA_ICDE07<DataType>::getAPLA_ICDE07(input_argument, normalized_series_vector, apla_array_vector[point_id]);
//				for (int segment_id = 0; segment_id < apla_array_vector[point_id].size(); segment_id++) {
//					APLA::getSegmentMinMaxPoint(normalized_series_vector, apla_array_vector[point_id][segment_id]);
//				}*/
//				/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//				break;
//			}
//			default:
//				assert(0);
//				break;
//			}
//
//			/*.................................................................................*/
//#ifdef _DEBUG
//			/*for (auto&& au : apla_array_vector[point_id]) { 210203
//				APLA::assert_segment_minmax(au);
//			}*/
//			assert(input_argument.point_dimension == apla_array_vector[point_id].size());
//#endif
//			/*.................................................................................*/
//
//			/*===================================================================        RTree Insertion     ========================================================================================*/
//			/*-------------------------------- APCA version: Initial RTree insertion coefficient -------------------------------------------*/
//			//typename APCA_QUAL::APCA CminParameter, CmaxParameter, MBRParameter; //Temp MBR
//			//APCA_QUAL::initialAPCA(apca, input_argument.point_dimension);
//			//APCA_QUAL::initialAPCA(CminParameter, input_argument.point_dimension);
//			//APCA_QUAL::initialAPCA(CmaxParameter, input_argument.point_dimension);
//			//APCA_QUAL::initialAPCA(MBRParameter, input_argument.point_dimension * 2);
//			/*-------------------------------------------------------------------------------------------------------------------------------*/
//			/*------------------------------------------Get Segment Min&Max point. min id < max id-------------------------------------------*/
//			//get_cmin_cmax_original(input_argument.time_series_length, dimension_id, apla_array_vector[point_id], rtree_rectangle_vector[point_id]);
//			/*-------------------------------------------------------------------------------------------------------------------------------*/
//			/*------------------------------------------Original version: Get Segment Min&Max point Cmin Cmax---------------------------------------------*/
//			//get_cmin_cmax(input_argument.time_series_length, dimension_id, apla_array_vector[point_id],  CminParameter, CmaxParameter);
//			/*---------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*------------------------------------------191127 APCA version: Get Segment Min&Max point Cmin Cmax-------------------------------------------*/
//			get_cmin_cmax_apca(normalized_series_vector, apla_array_vector[point_id], rectangle_insertion);
//			/*---------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*------------------------------------------191129 APCA version: PLA version-------------------------------------------------------------------*/
//			//get_cmin_cmax_pla(input_argument.time_series_length, apla_array_vector[point_id], rectangle_insertion);
//			/*---------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*--------------------APCA version:  get MBR for Rtree min max point ---------------------------------------------------------------*/
//			//APCA_QUAL::getMBR(CminParameter, CmaxParameter, MBRParameter);
//			/*----------------------------------------------------------------------------------------------------------------------------------*/
//			//copy_n(apca.r, apca.segmentNum, apca_array[point_id].APCALink.r + input_argument.point_dimension * file_id);
//			//copy_n(apca.v, apca.segmentNum, apca_array[point_id].APCALink.v + input_argument.point_dimension * file_id);.
//			/*--------------------APCA version:  get MBR for Rtree min max point ----------------------------------------------------------------*/
//			//copy_n(MBRParameter.r, MBRParameter.segmentNum, apca_MBR[point_id].r + MBRParameter.segmentNum * dimension_id);
//			//copy_n(MBRParameter.v, MBRParameter.segmentNum, apca_MBR[point_id].v + MBRParameter.segmentNum * dimension_id);
//			/*-----------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*.........................................................................................................*/
//#ifdef _DEBUG
//			assert(rectangle_insertion.m_min.size() == input_argument.point_dimension << 1);
//			for (int segment_id = 0; segment_id < input_argument.point_dimension << 1; segment_id++) {
//				//cout << apca_MBR[point_id].r[segment_id] << " " << apca_MBR[point_id].r[segment_id] << endl;
//				assert(rectangle_insertion.m_min[segment_id] != INF && rectangle_insertion.m_max[segment_id] != INF);
//			}
//#endif
//			/*.........................................................................................................*/
//
//			/*TOOL::printArray(apca_MBR[point_id].r, apca_MBR[point_id].segmentNum);
//			TOOL::printArray(apca_MBR[point_id].v, apca_MBR[point_id].segmentNum);*/
//			/*-------------------------Time Evaluation-----------------------------------*/
//			result_record.representation_time = TOOL::recordFinishTime(TOOL::time_record[0]);
//			input_argument.representation_time += result_record.representation_time;
//			/*---------------------------------------------------------------------------*/
//			/*------------------------------Rtree inserion-------------------------------*/
//			TOOL::recordStartTime(TOOL::time_record[1]);
//			RTree.Insert(rectangle_insertion.m_min, rectangle_insertion.m_max, point_id);
//			input_argument.build_rtree_time += TOOL::recordFinishTime(TOOL::time_record[1]);
//
//			/*------------------------------ 210618 RTree_partition -------------------------------*/
//			RTree_partition_object.Insert(rectangle_insertion.m_min, rectangle_insertion.m_max, point_id, apla_array_vector[point_id], normalized_series_vector);
//			/*------------------------------------------------------------------------------------*/
//
//			/*------------------------- Time Evaluation ---------------------------------*/
//			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//			input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);
//			/*---------------------------------------------------------------------------*/
//
//			/*----------------------------------------------------------------------------*/
//			/*=====================================================================================================================================================================*/
//
//			/*-------------------------------  Sum Deviation & Max Deviation   -------------------------------------------------*/
//			//result_record.sum_deviation = PLA_QUAL::getPLASumDeviation(normalized_series_vector, apla_array_vector[point_id]);
//			//output_argument.sum_deviation += result_record.sum_deviation;// except fisrt element
//			/*PLA_QUAL::getPLASumDeviation(normalized_series_vector, apla_array_vector[point_id], result_record);
//			output_argument.max_deviation += result_record.max_deviation;
//			output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
//			output_argument.sum_deviation += result_record.sum_deviation;*/
//
//			APLA::get_sum_deviation_no_ab(normalized_series_vector, apla_array_vector[point_id], result_record);
//			output_argument.max_deviation += result_record.max_deviation;
//			output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
//			output_argument.sum_deviation += result_record.sum_deviation;
//			/*----------------------------------------------------------------------------------------------------------------*/
//
//			/*.................................     Evaluate funciton  .......................................................*/
//#ifdef _DEBUG
//			assert(input_argument.point_number == apla_array_vector.size() && apla_array_vector[point_id].size() == input_argument.point_dimension);
//#endif
//			/*................................................................................................................*/
//
//			/*----    APCA version: Clear funciton  -------*/
//			TOOL::delete_rectangle(rectangle_insertion);
//			//APCA_QUAL::deleteAPCA(CminParameter);
//			//APCA_QUAL::deleteAPCA(CmaxParameter);
//			//APCA_QUAL::deleteAPCA(MBRParameter);
//			/*---------------------------------------------*/
//
//			//cout << "1 SAPLA Single sum deviation: " << result_record.sum_deviation << ", Accumulation sum deviation: " << output_argument.sum_deviation << ", max deviation: " << result_record.max_deviation << ", Accumulation max deviation: " << output_argument.max_deviation << endl;
//
//		}
//		/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//		else if (input_argument.representation_option == 2) {//PLA
//			/*================================= PLA =====================================================================*/
//
//			/*-------------------------------initial coefficient-----------------------------*/
//			//typename PLA_QUAL::PLA pla;
//			typename PLA_QUAL::PLA pla_MBR; //Temp MBR like APCA MBR
//			//PLA_QUAL::initialPLA(pla, input_argument.point_dimension);
//			PLA_QUAL::initialPLA(pla_MBR, input_argument.point_multi_single_dimension * 2);
//			/*-------------------------------------------------------------------------------*/
//			/*------------------------------get PLA approximation----------------------------*/
//			PLA_QUAL::getPLA(input_argument, original_time_series, pla_array[point_id]);
//			//PLA_QUAL::getPLA(input_argument, original_time_series, pla);
//			/*-------------------------------------------------------------------------------*/
//			/*------------------------------------get MBR------------------------------------*/
//			PLA_QUAL::getPLAMBR(pla_array[point_id], pla_MBR);
//			/*-------------------------------------------------------------------------------*/
//
//#ifdef _DEBUG
//			assert(input_argument.point_multi_single_dimension * 2 == input_argument.point_multi_single_dimension << 1);
//			for (int i = 0; i < pla_MBR.segmentNum; i++) {
//				//cout << pla_MBR.a[i] << ""<< pla_MBR.b[i] <<endl;
//				assert(pla_MBR.a[i] == pla_MBR.b[i] && pla_MBR.b[i] != INF);
//			}
//#endif
//			//(pla.a, pla.segmentNum, pla_array[point_id].a + input_argument.point_dimension * dimension_id);
//			//copy_n(pla.b, pla.segmentNum, pla_array[point_id].b + input_argument.point_dimension * dimension_id);
//			/*-------------------------Time Evaluation---------------------------------*/
//			result_record.representation_time = TOOL::recordFinishTime(TOOL::time_record[0]);
//			input_argument.representation_time += result_record.representation_time;
//			/*-------------------------------------------------------------------------*/
//			/*-----------------------------------Rtree inserion-------------------------------*/
//			TOOL::recordStartTime(TOOL::time_record[1]);
//			RTree.Insert(pla_MBR.a, pla_MBR.b, point_id);// a min, b max
//			input_argument.build_rtree_time += TOOL::recordFinishTime(TOOL::time_record[1]);
//			/*-------------------------Time Evaluation---------------------------------*/
//			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//			input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);
//			/*---------------------------------------------------------------------------*/
//			/*--------------------------------------------------------------------------------*/
//
//			/*--------------------------- Sum Deviation & max deviation -------------------------------------------*/
//			/*result_record.sum_deviation = PLA_QUAL::get_pla_sum_deviation(input_argument, normalized_series_vector, pla_array[point_id]);
//			output_argument.sum_deviation += result_record.sum_deviation;*/
//
//			PLA_QUAL::get_pla_sum_deviation(input_argument, normalized_series_vector, pla_array[point_id], result_record);
//			output_argument.max_deviation += result_record.max_deviation;
//			output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
//			output_argument.sum_deviation += result_record.sum_deviation;
//			/*----------------------------------------------------------------------------------------------------*/
//			//if (point_id == 0 && input_argument.arity_d == 1) {
//				/*cout << "Original time series: " << endl;
//				TOOL::printArray(original_time_series, input_argument.time_series_length);
//				PLA_QUAL::approximateOriginalFunctionPLA1(input_argument, pla, original_time_series);
//				cout << "Approximated PLA time series: " << endl;
//				TOOL::printArray(original_time_series, input_argument.time_series_length);
//				TOOL::writeArray("PLA Approximation", original_time_series, input_argument.time_series_length);*/
//				//}
//
//			//PLA_QUAL::deletePLA(pla);
//			PLA_QUAL::deletePLA(pla_MBR);
//
//			//cout << "1 PLA Single sum deviation: " << result_record.sum_deviation << ", Accumulation sum deviation: " << output_argument.sum_deviation << ", max deviation: " << result_record.max_deviation << ", Accumulation max deviation: " << output_argument.max_deviation << endl;
//
//		}
//		else if (input_argument.representation_option == 3 || input_argument.representation_option == 4 || input_argument.representation_option == 7 || input_argument.representation_option == 9) {// PAA & APCA & PAALM & SAX
//			/*================================= PAA & APCA =====================================================================*/
//			typename TOOL::RECTANGLE rectangle_insertion(input_argument.point_dimension << 1);
//
//			//typename APCA_QUAL::APCA apca;
//			//typename APCA_QUAL::APCA CminParameter, CmaxParameter, MBRParameter; //Temp MBR
//			//APCA_QUAL::initialAPCA(apca, input_argument.point_dimension);
//			//APCA_QUAL::initialAPCA(CminParameter, input_argument.point_dimension);
//			//APCA_QUAL::initialAPCA(CmaxParameter, input_argument.point_dimension);
//			///APCA_QUAL::initialAPCA(MBRParameter, input_argument.point_dimension * 2);
//
//			switch (input_argument.representation_option) {
//			case 3: {//APCA
//				APCA_QUAL::getAPCAPoint(original_time_series, input_argument.time_series_length, input_argument.point_dimension, apca_point_array[point_id]);
//				//APCA_QUAL::getAPCAPoint(original_time_series, input_argument.time_series_length, input_argument.point_dimension, apca);
//				break;
//			case 4: //PAA
//				APCA_QUAL::divideRemainderPAA(original_time_series, apca_point_array[point_id], input_argument.time_series_length, input_argument.point_dimension);
//				//APCA_QUAL::divideRemainderPAA(original_time_series, apca, input_argument.time_series_length, input_argument.point_dimension);
//
//				/*----------------Evaluation: Use mean value to instead min&max value. getCmin() getCmax() -------------------------*/
//				//CminParameter.segmentNum = apca.segmentNum;
//				//CmaxParameter.segmentNum = apca.segmentNum;
//				//for (int i = 0; i < apca.segmentNum; i++) {
//				//	CminParameter.r[i] = apca.r[i];
//				//	CmaxParameter.r[i] = apca.r[i];
//				//	CminParameter.v[i] = apca.v[i];//min[begin, end);
//				//	CmaxParameter.v[i] = apca.v[i];
//				//	//cout << "r[" << i << "] = " << Cmin.r[i] << ", Cmin0[" << i << "] = " << Cmin.v[i] << endl << endl;
//				//}
//				/*-------------------------------------------------------------------------------------------------------  ---------*/
//				break;
//			case 7:// PAALM 200108
//				APCA_QUAL::divideRemainderPAA(original_time_series, apca_point_array[point_id], input_argument.time_series_length, input_argument.point_dimension);
//				//APCA_QUAL::divideRemainderPAA(original_time_series, apca, input_argument.time_series_length, input_argument.point_dimension);
//				APCA_QUAL::get_paa_lagrangian(apca_point_array[point_id]);
//				break;
//			case 9:// SAX 200923
//				SaxQuantizer::SAX sax(input_argument.point_dimension);
//				/*APCA_QUAL::divideRemainderPAA(original_time_series, apca_point_array[point_id], input_argument.time_series_length, input_argument.point_dimension);
//				vector<double> mean_value_vector;
//				for (int i = 0; i < input_argument.point_dimension; i++) {
//					mean_value_vector.emplace_back(apca_point_array[point_id].v[i]);
//				}
//				vector<int> quantize_SAX_vector;
//				sax.quantize(mean_value_vector, &quantize_SAX_vector, false);
//				copy_n(quantize_SAX_vector.begin(), quantize_SAX_vector.size(), apca_point_array[point_id].v);*/
//				sax.get_SAX(normalized_series_vector, input_argument.point_dimension, apca_point_array[point_id]);
//				break;
//			}
//			default:
//				assert(0);
//				break;
//			}
//
//			//APCA_QUAL::getCmin(original_time_series, apca, CminParameter);
//			//APCA_QUAL::getCmax(original_time_series, apca, CmaxParameter);
//			/*---------------------------------191120 Original Rtree MBR-------------------------------------------*/
//			//RTREE::Rect rtree_rectangle;
//			//APCA_QUAL::get_minmax_original(normalized_series_vector, dimension_id, apca, rtree_rectangle_vector[point_id]);
//			/*-------------------------------------------------------------------------------------------------------*/
//			/*---------------------------------191127 APCA Paper Rtree MBR-------------------------------------------*/
//			APCA_QUAL::get_minmax_apca(normalized_series_vector, apca_point_array[point_id], rectangle_insertion);
//			//APCA_QUAL::get_minmax_apca(normalized_series_vector, dimension_id, apca, rtree_rectangle_vector[point_id]);
//			/*-------------------------------------------------------------------------------------------------------*/
//
//			/*..........................................................................................*/
//#ifdef _DEBUG
//			assert(rectangle_insertion.m_min.size() == input_argument.point_dimension << 1 && rectangle_insertion.m_min.size() == rectangle_insertion.m_max.size());
//			for (auto&& au : rectangle_insertion.m_min) {
//				assert(au != INF);
//			}
//			for (auto&& au : rectangle_insertion.m_max) {
//				assert(au != INF);
//			}
//#endif
//			/*..........................................................................................*/
//
//			//for (int array_id = 0; array_id < apca.segmentNum; array_id++) {
//				//apca.r[array_id] += input_argument.time_series_length * dimension_id;
//				//CminParameter.r[array_id] += input_argument.time_series_length * dimension_id;
//				//CmaxParameter.r[array_id] += input_argument.time_series_length * dimension_id;
//			//}
//			//APCA_QUAL::getMBR(CminParameter, CmaxParameter, MBRParameter);
//			//copy_n(apca.r, apca.segmentNum, apca_array[point_id].APCALink.r + input_argument.point_dimension * dimension_id);
//			//copy_n(apca.v, apca.segmentNum, apca_array[point_id].APCALink.v + input_argument.point_dimension * dimension_id);
//			//copy_n(MBRParameter.r, MBRParameter.segmentNum, apca_MBR[point_id].r + MBRParameter.segmentNum * dimension_id);
//			//copy_n(MBRParameter.v, MBRParameter.segmentNum, apca_MBR[point_id].v + MBRParameter.segmentNum * dimension_id);
//			/*-------------------------  Time Evaluation ---------------------------------*/
//			result_record.representation_time = TOOL::recordFinishTime(TOOL::time_record[0]);
//			input_argument.representation_time += result_record.representation_time;
//			/*----------------------------------------------------------------------------*/
//			/*--------------------------------------191128 Rtree insertion-------------------------------------------*/
//			TOOL::recordStartTime(TOOL::time_record[1]);
//			RTree.Insert(rectangle_insertion.m_min, rectangle_insertion.m_max, point_id);//191128 APCA paper MBR structure, min id == max id 
//			//RTree.Insert(MBRParameter.m_min, MBRParameter.m_max, point_id);//191128 APCA paper MBR structure, min id == max id 
//			input_argument.build_rtree_time += TOOL::recordFinishTime(TOOL::time_record[1]);
//			/*------------------------- Time Evaluation ---------------------------------*/
//			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//			input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);
//			/*---------------------------------------------------------------------------*/
//			/*-------------------------------------------------------------------------------------------------------*/
//			
//
//			/*--------------------201106  SAX has no sum deviation ---------------------------------*/
//			if (input_argument.representation_option == 9) {//SAX
//				output_argument.sum_deviation = 0;
//				output_argument.max_deviation = 0;
//				output_argument.max_deviation_multiple_width = 0;
//			}
//			else {
//				/*result_record.sum_deviation = APCA_KNN_QUAL::distanceAE(normalized_series_vector, normalized_series_vector.size(), apca_point_array[point_id]);
//				output_argument.sum_deviation += result_record.sum_deviation;*/
//
//				APCA_KNN_QUAL::distanceAE(normalized_series_vector, input_argument.time_series_length, apca_point_array[point_id], result_record);
//				output_argument.max_deviation += result_record.max_deviation;
//				output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
//				output_argument.sum_deviation += result_record.sum_deviation;
//			}
//			/*--------------------------------------------------------------------------------------*/
//			//APCA_QUAL::deleteAPCA(apca);
//			TOOL::delete_rectangle(rectangle_insertion);
//			//APCA_QUAL::deleteAPCA(CminParameter);
//			//APCA_QUAL::deleteAPCA(CmaxParameter);
//			//APCA_QUAL::deleteAPCA(MBRParameter);
//			//APCA_QUAL::deleteAPCA(apca);
//		}
//		else if (input_argument.representation_option == 5) {//CHEBY
//		/*=================================  CHEBY  =====================================================================*/
//			/*--------------------------------initial coefficient-----------------------------*/
//			assert(input_argument.point_multi_single_dimension > 0);
//			typename CHEBYSHEV_QUAL::CHEBYSHEV chebyshev_MBR;
//			CHEBYSHEV_QUAL::initialCHEBYSHEV(input_argument.time_series_length, input_argument.point_multi_single_dimension, chebyshev_MBR);
//			//
//			//typename CHEBYSHEV_QUAL::CHEBYSHEV chebyshev;
//			//CHEBYSHEV_QUAL::initialCHEBYSHEV(input_argument.time_series_length, input_argument.point_multi_single_dimension, chebyshev);
//			/*-------------------------------------------------------------------------------*/	
//			/*-----------------------------------------get chebyshev -----------------------------------------------------*/
//			CHEBYSHEV_QUAL::getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev_array[point_id]);
//			//CHEBYSHEV_QUAL::getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev);
//		    /*------------------------------------------------------------------------------------------------------------*/
//			/*----------------------------------------------get MBR ------------------------------------------------------*/
//			CHEBYSHEV_QUAL::getChebyshevMBR(chebyshev_array[point_id], chebyshev_MBR);
//			/*-------------------------------------------------------------------------------------------------------------*/
//			//
//			//				copy_n(chebyshev.coefficient, chebyshev.segmentNum, chebyshev_array[point_id].coefficient + input_argument.point_dimension * dimension_id);
//#ifdef _DEBUG
//			for (int i = 0; i < chebyshev_MBR.segmentNum; i++) {
//				//cout << chebyshev_MBR.f[i] << " "<< chebyshev_MBR.coefficient[i] <<endl;
//				assert(chebyshev_MBR.f[i] == chebyshev_MBR.coefficient[i] && chebyshev_MBR.coefficient[i] != INF);
//			}
//#endif
//			/*-------------------------   Time Evaluation   ---------------------------------*/
//			result_record.representation_time = TOOL::recordFinishTime(TOOL::time_record[0]);
//			input_argument.representation_time += result_record.representation_time;
//			/*--------------------------------------------------------------------------------*/
//				/*-----------------------------------Rtree inserion--------------------------------*/
//			TOOL::recordStartTime(TOOL::time_record[1]);
//			RTree.Insert(chebyshev_MBR.f, chebyshev_MBR.coefficient, point_id);// a min, b max
//			input_argument.build_rtree_time += TOOL::recordFinishTime(TOOL::time_record[1]);
//			/*-------------------------Time Evaluation---------------------------------*/
//			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//			input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);
//			/*---------------------------------------------------------------------------*/
//				/*---------------------------------------------------------------------------------*/
//			/*-----------------------------------------Sum Deviation----------------------------------------------*/
//			/*result_record.sum_deviation = CHEBYSHEV_QUAL::get_sum_deviation(input_argument, normalized_series_vector, chebyshev_share);
//			output_argument.sum_deviation += result_record.sum_deviation;*/
//			CHEBYSHEV_QUAL::get_sum_deviation(input_argument, normalized_series_vector, chebyshev_share, result_record);
//			output_argument.max_deviation += result_record.max_deviation;
//			output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
//			output_argument.sum_deviation += result_record.sum_deviation;
//			/*----------------------------------------------------------------------------------------------------*/
////
//				//CHEBYSHEV_QUAL::deleteCHEBYSHEV(chebyshev);
//			CHEBYSHEV_QUAL::deleteCHEBYSHEV(chebyshev_MBR);
//		}
//		else {
//			assert(0);
//		}
//		/*=============================================================================================================================================*/
//
//		/*==============================================     200901 Print each time series results     ==================================================*/
//		if ( input_argument.print_each_result == true ) {
//			cout << "!!!!!! Time Series id: " << point_id << endl;
//			//print_result_each_time_series(input_argument, normalized_series_vector, apla_array_vector[point_id], pla_array[point_id], apca_point_array[point_id], chebyshev_array[point_id], output_argument, result_record);
//			//print_result_each_time_series(input_argument, normalized_series_vector, apla_array_vector[point_id], pla_array[point_id], apca_point_array[point_id], chebyshev_array[point_id], output_argument, result_record);
//			TOOL::printInputArgument(input_argument);
//			switch (input_argument.representation_option) {
//			case 1:
//			case 6:
//			case 8: {//initial 200706
//				TOOL::print_each_segment_right_endpoint(apla_array_vector[point_id]);
//				break;
//			}
//			case 2: {//PLA
//				cout << " ======================================    PLA   ===========================================\n";
//				break;
//			}
//			case 3: //APCA
//			case 4: //PAA
//			case 7:// PAALM 200108
//			case 9: {// SAX
//				cout << "each right endpoint: ";
//				for (int segment_id = 0; segment_id < apca_point_array[point_id].segmentNum; segment_id++) {
//#ifdef _DEBUG
//					assert(apca_point_array[point_id].r[segment_id] != INF && apca_point_array[point_id].v[segment_id] != INF);
//#endif
//					cout << apca_point_array[point_id].r[segment_id] << ",";
//				}
//				cout << endl;
//				break;
//			}
//			case 5:
//				cout << " ======================================    Chebyshev   ===========================================\n";
//				break;
//			default:
//				assert(0);
//				break;
//			}	
//		}
//		/*===============================================================================================================================================*/
//
//		//cout << "2 Single sum deviation: " << result_record.sum_deviation<< ", Accumulation sum deviation: "<< output_argument.sum_deviation << ", max deviation: " << result_record.max_deviation << ", Accumulation max deviation: " << output_argument.max_deviation << endl;
//	}
//	//}
//	/*----- RTree: Delete Leaf node memory -------*/
//	RTree.deleteAllLeafNodeMemory();
//	/*--------------------------------------------*/
//	//file_stream.close();
//
//	/*cout << "Root Node : sub node number = " << RTree.m_root->m_count << " Root level = : " << RTree.m_root->m_level << "\n\nBegin to build a RTree:\n";
//	cout << "\n RTree conclusion\n The number of RTree Data Point = : " << RTree.Count() << endl;*/
//
//	//deletePLA(pla_array_MBR);
//	//deleteCHEBYSHEV(chebyshev_MBR);
//	//deleteMultiDimensionTrajectory(input_argument, multi_dimention_trajectories);
//	cout << "Sum Deviation: " << output_argument.sum_deviation << ", Max Deviation: " << output_argument.max_deviation << ", representation time: " << input_argument.representation_time << " us, build rtree time: " << input_argument.build_rtree_time << "us, whole time: " << input_argument.whole_run_time << endl;
//	delete[] original_time_series;
//	original_time_series = nullptr;
//	normalized_series_vector.clear();
//	normalized_series_vector.shrink_to_fit();
//}
//
////************************************
//// Method:buidRTreeIndex
//// Qualifier: KNN
//// date: 180916 1 MSPLA  2 PLA  3 APCA   4 PAA  5 CHEBY  6  ICDE07
//// author:
////************************************
//TEMPLATE
//RTREE& MULTI::buidRTreeIndex() {
//	//RTree<DataType, ElementType>& buidRTreeIndex(RTree<DataType, ElementType> &APCARTree, double* g_query_time_series, const double& g_time_series_length, const double& g_index_point_number, double(&test_d_original_time_series)[ROW][COLUMN], Link *APCALinkOriginal) {
//	//getCHEBYSHEV_SHARE(input_argument, chebyshev_share);
//#ifdef _DEBUG
//	assert(input_argument.build_rtree_time == 0.0);
//	printf(">>>***Build RTree Index***<<<\n");
//	//APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[0]);//whole build_rtree_time
//#endif
//	int point_id = NULL, j = NULL, f_insert_count = NULL;
//	//DataType* original_time_series = new DataType[input_argument.time_series_length];
//
//	/*cout << "Query Point : ";
//	for (i = 0; i < g_time_series_length; i++) cout << g_query_time_series[i] << ", ";
//	cout << endl;*/
//	//1 MSPLA  2 PLA  3 APCA   4 PAA  5 CHEBY  6  ICDE07
//	switch (input_argument.representation_option) {
//	case 1: {
//		for (point_id = 0; point_id < input_argument.point_number; point_id++) {
//			RTree.Insert(apca_MBR[point_id].r, apca_MBR[point_id].v, point_id);
//			f_insert_count++;
//		}
//		break;
//	}
//	case 2: {
//		for (point_id = 0; point_id < input_argument.point_number; point_id++) {
//			RTree.Insert(apca_MBR[point_id].r, apca_MBR[point_id].v, point_id);
//			f_insert_count++;
//		}
//		break;
//	}
//	case 3: {
//		typename PLA_QUAL::PLA pla_MBR; //Temp MBR like APCA MBR
//		PLA_QUAL::initialPLA(pla_MBR, input_argument.point_multi_single_dimension * 2);
//		for (point_id = 0; point_id < input_argument.point_number; point_id++) {
//			PLA_QUAL::getPLAMBR(pla_array[point_id], pla_MBR);
//			RTree.Insert(pla_MBR.a, pla_MBR.b, point_id);// a min, b max
//			f_insert_count++;
//		}
//
//		PLA_QUAL::deletePLA(pla_MBR);
//		break;
//	}
//	case 4: {
//		typename CHEBYSHEV_QUAL::CHEBYSHEV chebyshev_MBR;
//		CHEBYSHEV_QUAL::initialCHEBYSHEV(input_argument.time_series_length, input_argument.point_multi_single_dimension, chebyshev_MBR);
//		f_insert_count = 0;
//		for (point_id = 0; point_id < input_argument.point_number; point_id++) {
//			CHEBYSHEV_QUAL::getChebyshevMBR(chebyshev_array[point_id], chebyshev_MBR);
//			RTree.Insert(chebyshev_MBR.f, chebyshev_MBR.coefficient, point_id);// a min, b max
//			f_insert_count++;
//		}
//		CHEBYSHEV_QUAL::deleteCHEBYSHEV(chebyshev_MBR);
//		break;
//	}
//	default: assert(0); break;
//	}
//
//	//***print RTree***//
//	//RTree.printRTree();
//	/*------------------------------------- RTree: Delete Leaf node memory   ------------------------------------------------------------*/
//	RTree.deleteAllLeafNodeMemory();
//	/*-----------------------------------------------------------------------------------------------------------------------------------*/
//#ifdef _DEBUG
//	//APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[0], input_argument.build_rtree_time);
//	cout << "RTree build time: " << input_argument.build_rtree_time << " us" << endl;
//	//TOOL::writeSingleResult(input_argument.write_file_name[2], input_argument.build_rtree_time);
//#endif
//	/*================================================   KNN  =================================================================================*/
//	DataType* query_time_series = new DataType[input_argument.point_multi_single_length];
//	TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, 360, query_time_series);
//	TOOL::normalizeStandard(input_argument.time_series_length, query_time_series);//z-score normalization
//
//	switch (input_argument.representation_option) {
//	case 1:
//	case 2: {
//		APCA_KNN_QUAL::APCAKNNMulti(input_argument, query_time_series, input_argument.point_number, input_argument.time_series_length, input_argument.arity_d, RTree, apca_array, input_argument.K, input_argument.read_multiple_file_name);
//		//Base KNN
//		//priority_queue<TOOL::ID_DIST, vector<TOOL::ID_DIST>, TOOL::priorityDistanceEUC > base_queue;
//		//TOOL::SimpleBaseKNNSearch(input_argument, query_time_series, base_queue);
//		break;
//	}
//	case 3: {
//		PLA_QUAL::PLAKNNMulti(input_argument, query_time_series, RTree, pla_array, input_argument.K);
//		break;
//	}
//	case 4: {
//		CHEBYSHEV_QUAL::KNNCHEBYMulti(input_argument, query_time_series, chebyshev_share, chebyshev_array, RTree);
//		break;
//	}
//	default: assert(0); break;
//	}
//	/*==========================================================================================================================================*/
//	delete[] query_time_series;
//	query_time_series = nullptr;
//	return RTree;
//}
//
////************************************
//// Method:buid_rtree_knn
//// Qualifier: build Rtree & KNN. Add APLA & ICDE07
//// Input:
//// Output:
//// Notice:
//// date: 191108
//// author:
////************************************
//TEMPLATE
//RTREE& MULTI::buid_rtree_knn() {
//	//RTree<DataType, ElementType>& buidRTreeIndex(RTree<DataType, ElementType> &APCARTree, double* g_query_time_series, const double& g_time_series_length, const double& g_index_point_number, double(&test_d_original_time_series)[ROW][COLUMN], Link *APCALinkOriginal) {
//	//getCHEBYSHEV_SHARE(input_argument, chebyshev_share);
//#ifdef _DEBUG
//	assert(input_argument.build_rtree_time == 0.0);
//	printf(">>>***Build RTree Index***<<<\n");
//	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[0]);//whole build_rtree_time
//#endif
//	int point_id = NULL, j = NULL, f_insert_count = NULL;
//	//DataType* original_time_series = new DataType[input_argument.time_series_length];
//
//	/*cout << "Query Point : ";
//	for (i = 0; i < g_time_series_length; i++) cout << g_query_time_series[i] << ", ";
//	cout << endl;*/
//
//	/*================================Build Rtree ==================================================================*/
//	TOOL::recordStartTime(TOOL::time_record[14]);
//	switch (input_argument.representation_option) {
//	case 1: {
//		/*-----------------------------PAA-----------------------------------------*/
////#ifdef _DEBUG
////		cout << "PAA Insert: \n";
////		assert(input_argument.point_number == rtree_rectangle_vector.size());
////#endif
////		for (point_id = 0; point_id < input_argument.point_number; point_id++) {
////#ifdef _DEBUG
////			for (int i = 0; i < apca_MBR[point_id].segmentNum; i++) {
////				//cout << apca_MBR[point_id].r[i] << ", " << apca_MBR[point_id].v[i] << endl;
////				if (i & 1) {//odd value
////					assert(apca_MBR[point_id].r[i] <= apca_MBR[point_id].v[i]);
////					assert(rtree_rectangle_vector[point_id].m_min[i] <= rtree_rectangle_vector[point_id].m_max[i]);
////				}
////				else {//even id
////					assert(apca_MBR[point_id].r[i] == apca_MBR[point_id].v[i]);
////					assert(rtree_rectangle_vector[point_id].m_min[i] <= rtree_rectangle_vector[point_id].m_max[i]);
////				}
////			}
////#endif
////			/*-------------------------------------------old version----------------------------------------------------------*/
////			RTree.Insert(apca_MBR[point_id].r, apca_MBR[point_id].v, point_id); //old, APCA paper MBR ,min id == max id
////			/*----------------------------------------------------------------------------------------------------------------*/
////			/*-------------------------------------------new version----------------------------------------------------------*/
////			//RTree.Insert(rtree_rectangle_vector[point_id].m_min, rtree_rectangle_vector[point_id].m_max, point_id);//191119 new origianl MBR structure, min id < max id 
////			/*----------------------------------------------------------------------------------------------------------------*/
////				f_insert_count++;
////		}
//		//break;
//	}
//	case 2: {
//		/*-----------------------------PAA & APCA-----------------------------------------*/
//#ifdef _DEBUG
//		cout << "PAA & APCA Insert: \n";
//		assert(input_argument.point_number == rtree_rectangle_vector.size());
//#endif
//		for (point_id = 0; point_id < input_argument.point_number; point_id++) {
//#ifdef _DEBUG
//			for (int i = 0; i < apca_MBR[point_id].segmentNum; i++) {
//				//cout << apca_MBR[point_id].r[i] << ", " << apca_MBR[point_id].v[i] << endl;
//				if (i & 1) {//odd value
//					assert(apca_MBR[point_id].r[i] <= apca_MBR[point_id].v[i]);
//					assert(rtree_rectangle_vector[point_id].m_min[i] <= rtree_rectangle_vector[point_id].m_max[i]);
//				}
//				else {      //even id
//					assert(apca_MBR[point_id].r[i] == apca_MBR[point_id].v[i]);
//					assert(rtree_rectangle_vector[point_id].m_min[i] <= rtree_rectangle_vector[point_id].m_max[i]);
//				}
//			}
//#endif
//			/*-------------------------------------------old version----------------------------------------------------------*/
//			//RTree.Insert(apca_MBR[point_id].r, apca_MBR[point_id].v, point_id);
//			/*----------------------------------------------------------------------------------------------------------------*/
//			/*-------------------------------------------new version----------------------------------------------------------*/
//			RTree.Insert(rtree_rectangle_vector[point_id].m_min, rtree_rectangle_vector[point_id].m_max, point_id);//191119 new origianl MBR structure, min id < max id 
//			/*----------------------------------------------------------------------------------------------------------------*/
//			f_insert_count++;
//		}
//		break;
//	}
//	case 3: {
//		/*-----------------------------PLA-----------------------------------------*/
//#ifdef _DEBUG
//		cout << "PLA Insert: \n";
//#endif
//		typename PLA_QUAL::PLA pla_MBR; //Temp MBR like APCA MBR
//		PLA_QUAL::initialPLA(pla_MBR, input_argument.point_multi_single_dimension * 2);
//		for (point_id = 0; point_id < input_argument.point_number; point_id++) {
//			PLA_QUAL::getPLAMBR(pla_array[point_id], pla_MBR);
//
//#ifdef _DEBUG
//			for (int i = 0; i < pla_MBR.segmentNum; i++) {
//				//cout << pla_MBR.a[i] << ""<< pla_MBR.b[i] <<endl;
//				assert(pla_MBR.a[i] == pla_MBR.b[i]);
//			}
//#endif
//			RTree.Insert(pla_MBR.a, pla_MBR.b, point_id);// a min, b max
//			f_insert_count++;
//		}
//		PLA_QUAL::deletePLA(pla_MBR);
//		break;
//	}
//	case 4: {
//		/*-----------------------------Chebyshev-----------------------------------------*/
//#ifdef _DEBUG
//		cout << "Chebyshev Insert: \n";
//#endif
//		typename CHEBYSHEV_QUAL::CHEBYSHEV chebyshev_MBR;
//		CHEBYSHEV_QUAL::initialCHEBYSHEV(input_argument.point_multi_single_dimension, chebyshev_MBR);
//		f_insert_count = 0;
//		for (point_id = 0; point_id < input_argument.point_number; point_id++) {
//			CHEBYSHEV_QUAL::getChebyshevMBR(chebyshev_array[point_id], chebyshev_MBR);
//#ifdef _DEBUG
//			for (int i = 0; i < chebyshev_MBR.segmentNum; i++) {
//				//cout << chebyshev_MBR.f[i] << ""<< chebyshev_MBR.coefficient[i] <<endl;
//				assert(chebyshev_MBR.f[i] == chebyshev_MBR.coefficient[i]);
//			}
//#endif
//			RTree.Insert(chebyshev_MBR.f, chebyshev_MBR.coefficient, point_id);// a min, b max
//			f_insert_count++;
//		}
//		CHEBYSHEV_QUAL::deleteCHEBYSHEV(chebyshev_MBR);
//		break;
//	}
//	case 5:
//	case 6: {
//		/*-----------------------------APLA & ICDE07-----------------------------------------*/
//#ifdef _DEBUG
//		cout << "APLA & ICDE07 Insert: \n";
//#endif
//		for (point_id = 0; point_id < input_argument.point_number; point_id++) {// time series id 
//			/*---------------------------------APCA version------------------------------------------------*/
//			//RTree.Insert(apca_MBR[point_id].r, apca_MBR[point_id].v, point_id);//insert MBR, the same as APCA
//			/*---------------------------------------------------------------------------------------------*/
//			RTree.Insert(rtree_rectangle_vector[point_id].m_min, rtree_rectangle_vector[point_id].m_max, point_id);
//			f_insert_count++;
//		}
//		break;
//	}
//	default: assert(0);
//		break;
//	}
//	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//	/*==============================================================================================================================*/
//	//***print RTree***//
//	//RTree.printRTree();
//	/*------------------------------------- RTree: Delete Leaf node memory   ------------------------------------------------------------*/
//	RTree.deleteAllLeafNodeMemory();
//	/*-----------------------------------------------------------------------------------------------------------------------------------*/
//	APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[0], input_argument.build_rtree_time);
//	cout << "RTree build time: " << input_argument.build_rtree_time << " us" << endl;
//	//TOOL::writeSingleResult(input_argument.write_file_name[2], input_argument.build_rtree_time);
//
//	/*================================================   KNN  =================================================================================*/
//
//	/*--------------------------------------------------get query time series ---------------------------------------------------*/
//	DataType* query_time_series = new DataType[input_argument.point_multi_single_length];
//
//	if (input_argument.read_multiple_file_name) {
//		TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, input_argument.query_time_series_id, query_time_series);
//	}
//	else {
//		//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id)); input_argument.query_time_series_id
//		TOOL::getFileStreamByID(TOOL::getStringByID(TOOL::file_address, input_argument.file_id), input_argument.time_series_length, input_argument.query_time_series_id, query_time_series);
//	}
//
//
//	TOOL::normalizeStandard(input_argument.time_series_length, query_time_series);//z-score normalization
//	/*---------------------------------------------------------------------------------------------------------------------------*/
//
//	switch (input_argument.representation_option) {
//	case 1://PAA
//	case 7://PAALM
//	case 2: {//APCA
//		APCA_KNN_QUAL::APCAKNNMulti(input_argument, query_time_series, input_argument.point_number, input_argument.time_series_length, input_argument.arity_d, RTree, apca_array, input_argument.K, input_argument.read_multiple_file_name);
//		//Base KNN
//		//priority_queue<TOOL::ID_DIST, vector<TOOL::ID_DIST>, TOOL::priorityDistanceEUC > base_queue;
//		//TOOL::SimpleBaseKNNSearch(input_argument, query_time_series, base_queue);
//		break;
//	}
//	case 3: {//PLA
//		PLA_QUAL::PLAKNNMulti(input_argument, query_time_series, RTree, pla_array, input_argument.K);
//		break;
//	}
//	case 4: {//Chebyshev
//		CHEBYSHEV_QUAL::KNNCHEBYMulti(input_argument, query_time_series, chebyshev_share, chebyshev_array, RTree);
//		break;
//	}
//	case 5:
//	case 6: {//APLA & ICDE07
//		apla_knn_multi(input_argument, query_time_series, input_argument.point_number, input_argument.time_series_length, input_argument.arity_d, RTree, apla_array_vector, input_argument.K, input_argument.read_multiple_file_name);
//		break;
//	}
//	default: assert(0);
//		break;
//	}
//	/*==========================================================================================================================================*/
//	delete[] query_time_series;
//	query_time_series = nullptr;
//	return RTree;
//}
//
////200906 compute KNN accuracy
////191128 Optimization KNN. Add APLA & ICDE07 
////************************************
//// Method:compute_knn_accuracy
//// Qualifier:
//// Input:
//// Output:
//// Notice: Accuracy: Best case is 1, the worst case is 0, (Here is 0.1/K)
//// date: 200906
//// author:
////************************************
//TEMPLATE
//template<typename T>
//T MULTI::compute_knn_accuracy(const int& const K, const multiset<pair<T, int>>& const squential_scan_result_set, const multiset<pair<T, int>>& const knn_result_set){
//	assert(squential_scan_result_set.size() > K && knn_result_set.size() >= 0 && K > 0 && K != INF);
//	T knn_accuracy = 0;
//
//	set<int> sequential_time_series_id_set;
//	int result_id = 0;
//
//	/*########################  get K sequential time series id #######################################*/
//	for (multiset<pair<double, int>>::iterator squential_it = squential_scan_result_set.begin(); result_id < input_argument.K && squential_it != squential_scan_result_set.end(); squential_it++, result_id++) {
//		sequential_time_series_id_set.emplace(squential_it->second);
//	}
//	assert(sequential_time_series_id_set.size() == K);
//	/*#################################################################################################*/
//
//	/*######################## find if result set id in sequential result #############################################################*/
//	if (knn_result_set.size() <= K) {
//
//		for (multiset<pair<double, int>>::iterator knn_it = knn_result_set.begin(); knn_it != knn_result_set.end(); knn_it++) {
//			if (sequential_time_series_id_set.find(knn_it->second) != sequential_time_series_id_set.end()) {
//				knn_accuracy++;
//			}
//		}
//
//		knn_accuracy /= double(K);// The best case is 1
//	}
//	else {//ChebyShev result may set > K
//
//		result_id = 0;
//
//		for (multiset<pair<double, int>>::iterator knn_it = knn_result_set.begin(); result_id < K; knn_it++, result_id++) {
//			if (sequential_time_series_id_set.find(knn_it->second) != sequential_time_series_id_set.end()) {
//				knn_accuracy++;
//			}
//		}
//
//		knn_accuracy /= knn_result_set.size();// The best case is 1
//	}
//	
//	/*################################################################################################################################*/
//	
//	/*########################################         when accuracy = 0      ########################################################*/
//	if (knn_accuracy == 0) {
//		knn_accuracy = 0.01 / double(K);
//	}
//	/*################################################################################################################################*/
//
//	assert(knn_accuracy > 0 && knn_accuracy <= 1);
//
//	return knn_accuracy;
//}
//
////191128 Optimization KNN. Add APLA & ICDE07 
////************************************
//// Method:all_knn
//// Qualifier: KNN for 1 MSPLA  2 PLA  3 APCA   4 PAA  5 CHEBY  6  ICDE07
//// Input:
//// Output:
//// Notice: Use same index structure
//// date: 191128
//// author:
////************************************
//TEMPLATE
//template<typename T, typename Y, typename U>
//RTREE& MULTI::all_knn(typename TOOL::DATA_SOURCE& const data_source, const vector<T>& const query_time_series_vector, T*& query_time_series, const vector<Y>& const multi_y_projection_argument, const std::multiset<pair<U, int>>& const squential_scan_result_set) {
//
//	/*.........................................................................................................................................*/
//#ifdef _DEBUG
//	assert(query_time_series_vector.size() == input_argument.time_series_length);
//	for (int array_id = 0; array_id < query_time_series_vector.size(); array_id++) {
//		assert(query_time_series_vector[array_id] == query_time_series[array_id] && query_time_series_vector[array_id] != INF && query_time_series[array_id] != INF);
//	}
//#endif
//	/*........................................................................................................................................*/
//
//	/*-------------------------Time Evaluation---------------------------------*/
//	input_argument.knn_total_time = 0;
//	input_argument.result_accuracy = 0;
//	/*--------------------------------------------------------------------------*/
//	/*================================================  Approximation KNN  =================================================================================*/
//
//	/*--------------------------------------------------  get normalized query time series  ---------------------------------------------------*/
//	//assert(input_argument.file_id != INF);
//	//DataType* query_time_series = new DataType[input_argument.point_multi_single_length];
//	//if (input_argument.read_multiple_file_name) {
//	//	TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, input_argument.query_time_series_id, query_time_series);
//	//}
//	//else {
//	//	//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id)); input_argument.query_time_series_id
//	//	TOOL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, input_argument.query_time_series_id, query_time_series);
//	//}
//	//TOOL::normalizeStandard(input_argument.time_series_length, query_time_series);//z-score normalization
//	//vector<DataType> query_time_series_vector(input_argument.time_series_length, INF);
//	//copy_n(query_time_series, input_argument.time_series_length, query_time_series_vector.begin());
//	/*---------------------------------------------------------------------------------------------------------------------------------------*/
//
//	std::multiset<pair<U, int>> result_set;
//	switch (input_argument.representation_option) {
//	case 1:
//	case 6: //APLA & ICDE07
//	case 8: {//initial 200706
//		//result_set = all_knn_multi(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, apla_array_vector);
//		result_set = all_knn_multi_210511(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, apla_array_vector);
//		/*--------------------------210618-------------------------------------*/
//		result_set = all_knn_multi_210618(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree_partition_object, apla_array_vector);
//		/*---------------------------------------------------------------------*/ 
//		//result_set = all_knn_multi_210512(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, apla_array_vector);
//		break;
//	}
//	case 2: {//PLA
//		 //PLA_QUAL::PLAKNNMulti(input_argument, query_time_series, RTree, pla_array, input_argument.K);
//		//result_set = all_knn_multi(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, pla_array);
//		result_set = all_knn_multi_210511(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, pla_array);
//		break;
//	}
//	case 3: //APCA
//	case 4://PAA 
//	case 7: //PAALM
//	case 9: {//SAX
//		//result_set = all_knn_multi(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, apca_array);
//		result_set = all_knn_multi_210511(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, apca_array);
//		//Base KNN
//		//priority_queue<TOOL::ID_DIST, vector<TOOL::ID_DIST>, TOOL::priorityDistanceEUC > base_queue;
//		break;
//	}
//	case 5: {//Chebyshev
//		result_set = CHEBYSHEV_QUAL::KNNCHEBYMulti(input_argument, data_source, query_time_series, chebyshev_share, chebyshev_array, RTree);
//		break;
//	}
//	default: assert(0);
//		break;
//	}
//	/*==========================================================================================================================================*/
//
//	/*################################     get prune power & accuracy    #################################*/
//	//assert(squential_scan_result_set.size() > input_argument.K && result_set.size() >= 0 && input_argument.K > 0 && input_argument.K != INF);
//	//input_argument.result_accuracy = 0;
//
//	//int result_id = 0;
//	//for (multiset<pair<double, int>>::iterator squential_it = squential_scan_result_set.begin(), knn_it = result_set.begin(); result_id < input_argument.K && knn_it != result_set.end() && squential_it != squential_scan_result_set.end(); squential_it++, knn_it++, result_id++) {
//	//	if (float(squential_it->first) == float(knn_it->first)) {
//	//		input_argument.result_accuracy++;
//	//	}
//	//}
//	//input_argument.result_accuracy /= input_argument.K;// The best case is 1, 
//	//assert(input_argument.result_accuracy <= 1);
//
//	/*+++++++++++++++++++++++++++++++++++ Comput new  prune power +++++++++++++++++++++++++++++++++++++++++++++++*/
//	if (input_argument.representation_option != 5) {// Not Chebyshev, Chebyshev prune power may > 1
//		assert(input_argument.pruning_power > 0);//&& input_argument.pruning_power <= 1);//210603
//	}
//
//	input_argument.result_accuracy = compute_knn_accuracy(input_argument.K, squential_scan_result_set, result_set);
//	assert(input_argument.result_accuracy > 0 && input_argument.result_accuracy <= 1);
//	input_argument.prune_power_combine = input_argument.pruning_power / input_argument.result_accuracy;
//	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//
//	//switch (input_argument.representation_option) {
//	//case 1: {// SAPLA
//	//	//assert(input_argument.result_accuracy > 0);
//	//	if (input_argument.result_accuracy == 0) { input_argument.pruning_power = 20; }
//	//	break;
//	//}
//	//case 2: {// 2 PLA
//	//	if (input_argument.result_accuracy == 0) { input_argument.pruning_power = 20; }
//	//	break;
//	//}
//	//case 3: {// 3 APCA
//	//	if (input_argument.result_accuracy == 0) { input_argument.pruning_power = 20; }
//	//	break;
//	//}
//	//case 4:// 4 PAA
//	//	if (input_argument.result_accuracy == 0) { input_argument.pruning_power = 20; }
//	//	break;
//	//case 5:// CHEBY
//	//	if (input_argument.result_accuracy == 0) { input_argument.pruning_power = 20; }
//	//	break;
//	//case 6:// ICDE07
//	//	if (input_argument.result_accuracy == 0) { input_argument.pruning_power = 20; }
//	//	break;
//	//case 7:// PAALM
//	//	if (input_argument.result_accuracy == 0) { input_argument.pruning_power = 20; }
//	//	break;
//	//default:
//	//	assert(0);
//	//	break;
//	//}
//	/*###########################################################################################################*/
//
//	/*-----------------------------------------Evaluation------------------------------------*/
//	output_argument.run_time = input_argument.whole_run_time;
//	cout << "Accuracy: " << input_argument.result_accuracy << ", KNN time: " << input_argument.knn_total_time<< ", KNN time has IO: " << input_argument.knn_total_time_has_IO << "us, whole time: " << input_argument.whole_run_time << endl;
//	/*---------------------------------------------------------------------------------------*/
//	/*delete[] query_time_series;
//	query_time_series = nullptr;*/
//	result_set.clear();
//
//	return RTree;
//}
//
////************************************
//// Method:apla_knn_multi
//// Qualifier: KNN for APLA & ICDE07
//// Input:
//// Output:
//// Notice: Basic KNN, most copy from APCA KNN 
//// date: 191111
//// author:
////************************************
//TEMPLATE
//template<typename T>
//bool MULTI::apla_knn_multi(typename TOOL::INPUT_ARGUMENT& const input_argument, const DataType* g_query_time_series, const DataType& g_index_point_number, const DataType& g_time_series_length, const int& arity_d, const RTREE& apcaRTree, std::vector<DoublyLinkedList<T>>& apla_array_vector, const int& K, string*& const multi_file_name) {
//	cout << "APLA & ICDE07 KNN Search function: " << endl;
//#ifdef _DEBUG
//	assert(input_argument.K <= input_argument.point_number && input_argument.time_series_length == g_time_series_length && input_argument.point_number == g_index_point_number && input_argument.arity_d == arity_d && input_argument.K == K);
//#endif
//	/*-------------------------------------change pointer array to vector array--------------------------------*/
//	int multi_to_single_series_length = input_argument.time_series_length * input_argument.arity_d;
//	vector<DataType> query_time_series_vector;
//	query_time_series_vector.resize(multi_to_single_series_length, INF);
//	std::copy_n(g_query_time_series, multi_to_single_series_length, query_time_series_vector.begin());
//	/*---------------------------------------------------------------------------------------------------------*/
//
//	/*=======================================Evaluation: plot & bar chart=======================================================*/
//#ifdef _DEBUG
//	//g_n_account_apca_point = 0;
//	input_argument.knn_total_time = 0.0;
//	input_argument.navigate_index_time = 0.0;// navigate time
//	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0.0;// distance euclidean time
//#endif
//	input_argument.IO_cost = 0;// measure I/O cost
//	input_argument.sum_distance_euc = 0;
//	/*===========================================================================================================================*/
//	/*========================================Variable construct=================================================================*/
//	int i = NULL, j = NULL;
//	priority_queue< MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>, vector< MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>>, MULTI::priorityIncrement<MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT> >> queue;// <Rtree node, distance>
//	APLA_NODE_PAIR<APLA::AREA_COEFFICIENT> f_APLA_Root, f_temp_APLA_Pair;// distance, APCA point(original time series) id, APCA point pointer, Rtree sub node pointer
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> temp; // <dist, original time series id>
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> result;// <dist, original time series id>
//	typename APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair; // <dist, original time series id>
//	/*===========================================================================================================================*/
//	TOOL::recordStartTime(TOOL::time_record[14]);
//#ifdef _DEBUG
//	TOOL::recordStartTime(TOOL::time_record[13]);//for total KNN time
//#endif
//	/*========================================== Variable Initial ===============================================================*/
//	f_APLA_Root.p_rtree_node = apcaRTree.m_root;
//	f_APLA_Root.d_dist = 0;
//	queue.push(f_APLA_Root);
//	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
//	//printf("<///////**    KNN Begin   **////////>\n");
//	/*===========================================================================================================================*/
//#ifdef _DEBUG
//	TOOL::recordStartTime(TOOL::time_record[4]);//navigate index time
//#endif
//	/*========================================================================================Begin K-NN=======================================================================================================================================*/
//	while (!queue.empty()) {
//		APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>  m_temp_queue_top = queue.top();
//
//		//cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator: ";
//		//for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
//		//cout << endl;
//
//		/*===========================================================================Scan temp to get ruslt=========================================================================================================*/
//		for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
//			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
//			if (plist->d_dist <= m_temp_queue_top.d_dist) {
//				//cout << "           <= " << endl;
//				result.push_back(*plist);
//				plist = temp.erase(plist);
//			}
//			else
//				plist++;
//
//			if (K == result.size()) {
//				input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//				input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
//#ifdef _DEBUG
//				input_argument.knn_total_time = TOOL::recordFinishTime(TOOL::time_record[13]);
//				assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF);
//				assert(input_argument.IO_cost == input_argument.sum_distance_euc);
//				/*-------------------------------------------- Print Result ---------------------------------*/
//				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! APLA Find result !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    result list size: " << result.size() << endl;
//
//				cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
//
//				cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;
//
//				cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;
//
//				cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;
//
//
//				cout << "pruning power: " << input_argument.pruning_power << endl;
//
//				cout << "I/O cost: " << input_argument.IO_cost << endl;
//
//				result.sort([](const APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& first, const  APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist; });//small to big
//
//				for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//					cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//				}
//#endif
//				/*--------------------------------------------------------------------------------------------*/
//				/*--------------------------------------     Clear memory    -----------------------------------------------*/
//				f_temp_APLA_Pair.approximation_pointer = nullptr;
//				f_APLA_Root.approximation_pointer = nullptr;
//				//priority_queue< MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>, vector< MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>>, MULTI::priorityIncrement<MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>>>().swap(queue);
//
//				temp.clear();
//				result.clear();
//				//list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//				//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//				/*---------------------------------------------------------------------------------------------------------*/
//				return true;
//			}
//		}
//		/*==========================================================================================================================================================================*/
//		/*======================================================================Pop top node ine queue==============================================================================*/
//		queue.pop();
//		/*==========================================================================================================================================================================*/
//		/*====================================================================== Check APLA point ==================================================================================*/
//		if (m_temp_queue_top.p_rtree_node == nullptr) { //is APCA data point
//			input_argument.navigate_index_time += TOOL::recordFinishTime(TOOL::time_record[4]);
//
//			//cout << "    queue.top is data point\n";
//			DataType* original_time_series = new DataType[g_time_series_length * arity_d];
//
//			/*------------------------test time -----------------------------*/
//#ifdef _DEBUG
//			double run_time0;
//			//recordStartTime(whole_first_time_start[3], whole_first_dqFreq[3]);
//			TOOL::recordStartTime(TOOL::time_record[3]);
//#endif
//			/*---------------------------------------------------------------*/
//
//			/*------------------------------------------------get original time series from database (file)----------------------------------------------*/
//			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
//			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
//			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;//get top queue time sereis id.
//			//getFileStreamByID(file_name, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			// get original time sires from database(txt file)
//
//			if (input_argument.read_multiple_file_name) {
//				TOOL::getMultiFoldToSingleByID(multi_file_name, arity_d, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			}
//			else {
//				//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id));
//				TOOL::getFileStreamByID(TOOL::getStringByID(TOOL::file_address, input_argument.file_id), input_argument.time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			}
//			TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);
//			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*===============Evaluation==============*/
//			input_argument.sum_distance_euc++;
//			input_argument.IO_cost++;
//#ifdef _DEBUG
//			//g_n_account_apca_point++;
//			TOOL::recordStartTime(TOOL::time_record[6]);//time for euclidean distance
//#endif
//			/*=======================================*/
//
//			/*------------------------------------------------compute Euclidean distance from query times and candidate original time series-----------------------------------------*/
//			tempOriginalTimeSeriesPair.d_dist = APCA_KNN_QUAL::distanceEUC(g_query_time_series, g_time_series_length * arity_d, original_time_series, g_time_series_length * arity_d);
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*===============Evaluation==============*/
//#ifdef _DEBUG
//			input_argument.distance_euc_time += TOOL::recordFinishTime(TOOL::time_record[6]);
//#endif
//			/*=======================================*/
//
//			/*------------------------------------------------temp: insert time seris dist into temp---------------------------------------------------------------------------------*/
//			//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
//			temp.push_back(tempOriginalTimeSeriesPair);
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*-------------------------------------------------Clear memory----------------------------------------------------------------------------------------------------------*/
//			delete[] original_time_series;
//			original_time_series = nullptr;
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//		/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
//			/*====================================================================== Check Leaf Node, insert APCA point into queue===============================================================================*/
//			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS / 2 == input_argument.point_dimension);
//#endif
//			//typename APCA_QUAL::APCA QProjection;
//			//APCA_QUAL::initialAPCA(QProjection, input_argument.point_dimension);
//			vector<APLA::AREA_COEFFICIENT> query_apla_projection_vector(input_argument.point_dimension);
//
//			/*-----------------------------------------------------------------------queue: insert leaf node. compute distance LB--------------------------------------------------------------------------------------------*/
//			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {// scan sub node in Rtree
//#ifdef _DEBUG
//				//g_n_account_leaf_node++;
//#endif
//				/*-----------------------------------------get ACPA(origitnal time series) point & id----------------------------------------------------*/
//				f_temp_APLA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);// get APCA(original time series) id
//				//cout << tempAPCAPair.id_originalTimeSeries << endl;
//				//f_temp_APLA_Pair.approximation_pointer = &apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)];// pointer to APCA point
//				/*----------------------------------------------------------------------------------------------------------------------------------------*/
//				//cout << "KNN : data point ID = " << tempAPCAPair.id_originalTimeSeries << endl;
//				//cout << tempAPCAPair.APCAValue << endl;
//#ifdef _DEBUG
//				TOOL::recordStartTime(TOOL::time_record[5]);
//#endif
//				//f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point), *f_temp_APCA_Pair.p_APCA_point);
//				//f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceLB(APCA_KNN_QUAL::QAPCAProjection(g_query_time_series, g_time_series_length * arity_d, *f_temp_APLA_Pair.approximation_pointer, QProjection), *f_temp_APLA_Pair.approximation_pointer);
//				f_temp_APLA_Pair.d_dist = APLA::get_distance_LB(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_apla_projection_vector));
//#ifdef _DEBUG
//				input_argument.distance_lowbound_time += TOOL::recordFinishTime(TOOL::time_record[5]);
//#endif
//				//printf("\nrun_time = %f us\n", run_time0);
//
//				//cout << tempAPCAPair.key;
//				f_temp_APLA_Pair.p_rtree_node = nullptr;
//				//cout << tempAPCAPair.key;
//#ifdef _DEBUG   
//				DoublyLinkedList<T> test_distance_LB = DoublyLinkedList<T>();
//				test_distance_LB.copy(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//				assert(f_temp_APLA_Pair.d_dist == APLA::get_distance_LB_by_series_apla(query_time_series_vector, test_distance_LB));
//				double queue_run_time;
//				TOOL::recordStartTime(TOOL::time_record[7]);
//#endif
//				//cout << "            Push apca data. Branch id: " << i << " DLB: " << f_temp_APCA_Pair.d_dist << ", time series ID: " << f_temp_APCA_Pair.original_time_series_id << endl;
//				//f_temp_APLA_Pair.approximation_pointer = nullptr;
//				queue.push(f_temp_APLA_Pair);
//#ifdef _DEBUG
//				//g_n_time_leaf_node_push_count += TOOL::recordFinishTime(TOOL::time_record[7], queue_run_time);
//#endif
//				//printf("\nrun_time = %f us\n", run_time0);
//				//cout << "KNN : queue.size() = " << queue.size() << endl;
//				//		//cout << tempAPCAPair.key;
//			}
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
//			//	//cout << tempOriginalTimeSeriesPair.key;
//			/*------------------------------------------------------------------------Clear memory-------------------------------------------------------------------------------*/
//			//APCA_QUAL::deleteAPCA(QProjection);
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS / 2 == input_argument.point_dimension && input_argument.point_multi_single_length == g_time_series_length * arity_d);
//#endif//is internal Node
//			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//#ifdef _DEBUG	
//			double run_time2;
//			TOOL::recordStartTime(TOOL::time_record[8]);
//#endif
//			APLA_NODE_PAIR<APLA::AREA_COEFFICIENT> tempApcaPair;
//			typename APCA_KNN_QUAL::REGION fs_region_G;
//			typename APCA_KNN_QUAL::initialREGION(fs_region_G, input_argument.point_dimension);
//			////cout << "regionNum = " << G.regionNum << endl;
//
//			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
//				/*------------------------------------Evaluate-------------------------------------------*/
//#ifdef _DEBUG	
//				//g_n_account_child_node++;
//				double mindistQR_run_time;
//				TOOL::recordStartTime(TOOL::time_record[9]);
//				TOOL::recordStartTime(TOOL::time_record[5]);//distance PLA & MBR time
//#endif
//				/*---------------------------------------------------------------------------------------*/
//				// APCA paper, min id == max id
//				tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//				//original min id < max id
//				//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//#ifdef _DEBUG	
//				input_argument.distance_lowbound_time += TOOL::recordFinishTime(TOOL::time_record[5]);
//				/*for (int i = 0; i < fs_region_G.regionNum; i++) {
//				cout << "            G1[" << i << "] = " << fs_region_G.G1[i] << ", G2[" << i << "] = " << fs_region_G.G2[i] << ", G3[" << i << "] = " << fs_region_G.G3[i] << ", G4[" << i << "] = " << fs_region_G.G4[i] << endl;
//				}*/
//				//g_n_time_child_node_MINDIST_count += TOOL::recordFinishTime(TOOL::time_record[9], mindistQR_run_time);
//#endif
//				//printf("\nrun_time = %f us\n", run_time0);
//
//				tempApcaPair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//
//#ifdef _DEBUG
//				double push_run_time;
//				TOOL::recordStartTime(TOOL::time_record[10]);
//#endif
//				//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//				queue.push(tempApcaPair);
//				//cout << "KNN : queue.size() = " << queue.size() << endl;
//#ifdef _DEBUG
//				//g_n_time_child_node_push_count += TOOL::recordFinishTime(TOOL::time_record[10], push_run_time);
//#endif
//				//printf("\nrun_time = %f us\n", run_time0);
//				tempApcaPair.approximation_pointer = nullptr;
//			}
//
//			APCA_KNN_QUAL::deleteREGION(fs_region_G);
//#ifdef _DEBUG
//			//g_n_time_count2 += TOOL::recordFinishTime(TOOL::time_record[8], run_time2);
//#endif
//			//printf("\nrun_time = %f us\n", run_time0);
//		}
//		else {
//			assert(0);
//		}
//	}
//	/*=====================================================================================================================================================================================================================================================================================================================*/
//	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//	input_argument.knn_total_time = TOOL::recordFinishTime(TOOL::time_record[13]);
//	input_argument.pruning_power = 2;
//
//#ifdef _DEBUG
//	assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF);
//	cout << "??????????????????     APLA KNN failed    ??????????????????" << endl;
//	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	/*ofstream outfile(input_argument.write_file_name[16] + ".txt", ios::app);
//	assert(outfile.is_open());
//	outfile << "??????????????????     APLA KNN failed    ??????????????????   " << PAA_or_APCA << endl;
//	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
//	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	outfile.close();*/
//	TOOL::printInputArgument(input_argument);
//#endif
//
//	//assert(0);
//	priority_queue<MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>, vector< MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>>, MULTI::priorityIncrement<MULTI::APLA_NODE_PAIR<APLA::AREA_COEFFICIENT>>>().swap(queue);
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//	//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//	f_temp_APLA_Pair.approximation_pointer = nullptr;
//	f_APLA_Root.approximation_pointer = nullptr;
//	return false;
//}
//
////************************************
//// Method:all_knn_multi
//// Qualifier: KNN structure for PAA, APCA, APLA, ICDE07
//// Input:
//// Output:
//// Notice: Basic KNN, most copy from APCA KNN 
//// date: 191128
//// author:
////************************************
//TEMPLATE
//template<typename T, typename Y>
//std::multiset<pair<double, int>> MULTI::all_knn_multi(typename TOOL::DATA_SOURCE& const data_source, typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<Y>& const multi_y_projection_argument, DataType*& g_query_time_series, const RTREE& const apcaRTree, T& const approximation_array) {
//#ifdef _DEBUG
//	assert(input_argument.file_id != INF && input_argument.K <= input_argument.point_number && input_argument.time_series_length != INF && input_argument.point_number != INF && input_argument.arity_d != INF && input_argument.K != INF);
//	assert(input_argument.time_series_length == input_argument.point_multi_single_length);
//	assert(input_argument.point_dimension == input_argument.point_multi_single_dimension);
//
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//		cout << "MSPLA(SAPLA KNN: \n";
//		break;
//	case 2://PLA
//		cout << "PLA KNN: \n";
//		break;
//	case 3://APCA
//		cout << "APCA KNN: \n";
//		break;
//	case 4://PAA
//		cout << "PAA KNN: \n";
//		break;
//	case 5://Chebyshev
//		cout << "CHEBY KNN: \n";
//		break;
//	case 6://ICDE07
//		cout << "ICDE07 KNN: \n";
//		break;
//	case 7://PAALM
//		cout << "PAALM KNN: \n";
//		break;
//	case 8://initial 200706
//		cout << "Initial 200706 KNN: \n";
//		break;
//	case 9://SAX
//		cout << "SAX KNN: \n";
//		break;
//	default:
//		assert(0);
//		break;
//	}
//#endif
//	/*-------------------------      Evaluation  ---------------------------------*/
//	input_argument.knn_total_time = 0.0;
//	input_argument.IO_cost = 0;// measure I/O cost
//	input_argument.pruning_power = 0;
//	/*-------------------------------------------------------------------------------*/
//
//	/*-------------------------------------change pointer array to vector array--------------------------------*/
//	//int multi_to_single_series_length = input_argument.time_series_length * input_argument.arity_d;
//	int multi_to_single_series_length = input_argument.time_series_length;
//	vector<DataType> query_time_series_vector;
//	query_time_series_vector.resize(multi_to_single_series_length, INF);
//	std::copy_n(g_query_time_series, multi_to_single_series_length, query_time_series_vector.begin());
//
//	vector<DataType> reconstruct_query_time_series_vector;
//	reconstruct_query_time_series_vector.resize(multi_to_single_series_length, INF);
//	/*---------------------------------------------------------------------------------------------------------*/
//	/*=======================================Evaluation: plot & bar chart=======================================================*/
//#ifdef _DEBUG
//	//g_n_account_apca_point = 0;
//	input_argument.navigate_index_time = 0.0;// navigate time
//	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0.0;// distance euclidean time
//#endif
//	/*===========================================================================================================================*/
//
//	/**/
//	typename APCA_QUAL::APCA query_APCA;
//	SaxQuantizer::SAX sax(input_argument.point_dimension);
//	//APCA_QUAL::getAPCAPoint(g_query_time_series, input_argument.time_series_length, input_argument.point_dimension, query_APCA);
//	typename PLA_QUAL::PLA PLA_query;
//	/**/
//
//	/*========================================!!!!!!!!!!!   MSPLA & ICDE07 Variable construct=================================================================*/
//	
//	DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>();
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//	//APLA_ICDE07<DataType>::getAPLA_ICDE07(input_argument, normalized_series_vector, apla_array_vector[point_id]);
//	case 6: //ICDE07
//	case 8: {//Initial 200706
//		//typename TOOL::Y_PROJECTION_ARGUMENT y_projection_argument(input_argument.initial_N);
//		//DoublyLinkedList<APLA::AREA_COEFFICIENT> all_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//		//DoublyLinkedList<APLA::AREA_COEFFICIENT> cluster_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//		query_linked_list.copy(apla_array_vector[input_argument.query_time_series_id]);
//		//APLA::get_APLA_point(input_argument, g_query_time_series, y_projection_argument, all_linked_list, cluster_linked_list, query_linked_list);//191129
//
//		APLA::getAPLAReconstructSeries(query_linked_list, reconstruct_query_time_series_vector);
//
//		/*......................................*/
//#ifdef _DEBUG
//		assert(query_linked_list.size() == apla_array_vector[input_argument.query_time_series_id].size());
//		for (int segment_id = 0; segment_id < query_linked_list.size(); segment_id++) {
//			const auto& const query_segment = query_linked_list[segment_id];
//			const auto& const array_segment = apla_array_vector[input_argument.query_time_series_id][segment_id];
//			//cout <<"!!!!!!!!!!!!" <<query_segment.right_endpoint << ", "<<array_segment.right_endpoint << endl;
//			assert(query_segment.right_endpoint == array_segment.right_endpoint);
//			assert(query_segment.right_endpoint != INF && query_segment.right_endpoint >= 0 && query_segment.right_endpoint < input_argument.point_multi_single_length);
//		}
//#endif
//		/*......................................*/
//
//		break;
//	}
//	case 2: {//PLA
//		PLA_QUAL::initialPLA(PLA_query, input_argument.point_multi_single_dimension);//????????180918 , this multi ot single has problem
//		PLA_QUAL::getPLA(input_argument.point_multi_single_length, input_argument.point_multi_single_dimension, g_query_time_series, PLA_query);
//		break;
//	}
//	case 3://APCA
//	case 4://PAA
//		break;
//	case 5://CHEBY
//		assert(0);
//		break;
//	case 7://PAALM
//		break;
//	case 9: {//SAX
//		APCA_QUAL::initialAPCA(query_APCA, input_argument.point_dimension);
//		sax.get_SAX(query_time_series_vector, input_argument.point_dimension, query_APCA);
//		break;
//	}
//	default:
//		assert(0);
//		break;
//	}
//	/*=======================================================================================================================================================*/
//	/*========================================!!!!!!!!!!!   APLA & ICDE07 Variable construct=================================================================*/
//	int i = NULL, j = NULL;
//	priority_queue<RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>> queue;// <Rtree node, distance>
//	RTREE_NODE_PAIR f_APLA_Root, f_temp_APLA_Pair;// distance, APCA point(original time series) id, APCA point pointer, Rtree sub node pointer
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> temp; // <dist, original time series id>
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> result;// <dist, original time series id>
//	std::multiset<pair<double, int>> result_set;//191204 result
//	typename APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair; // <dist, original time series id>
//	RTREE_NODE_PAIR m_temp_queue_top;//
//	vector<DataType> original_time_series_vector;
//	/*==========================================================================================================================================*/
//	/*-------------------------Evaluation Time---------------------------------*/
//	TOOL::recordStartTime(TOOL::time_record[14]); //whole time
//	TOOL::recordStartTime(TOOL::time_record[2]);//knn time
//	/*--------------------------------------------------------------------------------*/
//	/*========================================== Variable Initial ===============================================================*/
//	f_APLA_Root.p_rtree_node = apcaRTree.m_root;
//	f_APLA_Root.d_dist = 0;
//	queue.push(f_APLA_Root);
//	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
//	//printf("<///////**    KNN Begin   **////////>\n");
//	/*===========================================================================================================================*/
//
//	/*========================================================================================Begin K-NN=======================================================================================================================================*/
//	while (!queue.empty()) {
//
//		m_temp_queue_top = queue.top();
//
//		//cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator: ";
//		//for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
//		//cout << endl;
//
//		/*===========================================================================Scan temp to get ruslt=========================================================================================================*/
//		for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
//			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
//			if (plist->d_dist <= m_temp_queue_top.d_dist) {
//				//cout << "           <= " << endl;
//				result.push_back(*plist);
//				plist = temp.erase(plist);
//			}
//			else
//				plist++;
//
//			if (input_argument.K == result.size()) {
//				/*-------------------------Evaluation ---------------------------------*/
//				input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
//				input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//				assert(input_argument.IO_cost != 0 && input_argument.point_number != INF);
//				input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
//				assert(input_argument.pruning_power > 0 && input_argument.pruning_power < 1);
//				/*---------------------------------------------------------------------*/
//#ifdef _DEBUG
//				assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF);
//				/*-------------------------------------------- Print Result ---------------------------------*/
//				cout << "!!!!!!!!!!!!! KNN Find result !!!!!!!!!!!!!!!!!!!!!!!!!!!!   result list size: " << result.size() << endl;
//
//				cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
//
//				cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;
//
//				cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;
//
//				cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;
//
//				cout << "pruning power: " << input_argument.pruning_power << endl;
//
//				cout << "I/O cost: " << input_argument.IO_cost << endl;
//
//				result.sort([](const APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& first, const  APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist; });//small to big
//#endif
//				for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//					//cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//					result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
//				}
//				assert(result_set.size() == input_argument.K);
//				//cout << endl;
//				/*--------------------------------------------------------------------------------------------*/
//				/*--------------------------------------     Clear memory    -----------------------------------------------*/
//				//f_temp_APLA_Pair.approximation_pointer = nullptr;
//				//f_APLA_Root.approximation_pointer = nullptr;
//				priority_queue< RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>>().swap(queue);
//
//				temp.clear();
//				result.clear();
//				list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//				//list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//				/*---------------------------------------------------------------------------------------------------------*/
//				return result_set;
//			}
//		}
//		/*==========================================================================================================================================================================*/
//		/*======================================================================Pop top node ine queue==============================================================================*/
//		queue.pop();
//		/*==========================================================================================================================================================================*/
//		/*====================================================================== Check APLA point ==================================================================================*/
//		if (m_temp_queue_top.p_rtree_node == nullptr) { //is Approximation data point
//
//			//cout << "    queue.top is data point\n";
//			//DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
//			//vector<DataType> original_time_series_vector;
//
//			/*------------------------------------------------get original time series from database (file)----------------------------------------------*/
//			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
//			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
//			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;//get top queue time sereis id.
//			//getFileStreamByID(file_name, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			// get original time sires from database(txt file)
//
//#ifdef _DEBUG
//			assert(input_argument.time_series_length == data_source.multi_single_time_series_length && input_argument.read_file_name == data_source.read_file_address);
//#endif
//
//			TOOL::read_normalized_multi_time_series(data_source, m_temp_queue_top.original_time_series_id, original_time_series_vector);
//			//			//already normalized
//			//			TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series_vector);
//			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*===============Evaluation==============*/
//			input_argument.IO_cost++;
//			/*=======================================*/
//
//			/*------------------------------------------------compute Euclidean distance from query times and candidate original time series-----------------------------------------*/
//			tempOriginalTimeSeriesPair.d_dist = TOOL::distanceEUC(query_time_series_vector, original_time_series_vector);
////#ifdef _DEBUG
////			copy_n(original_time_series_vector.begin(), original_time_series_vector.size(), original_time_series);
////			double test_distance = APCA_KNN_QUAL::distanceEUC(g_query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);
////			assert(test_distance == tempOriginalTimeSeriesPair.d_dist);
////#endif
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*------------------------------------------------temp: insert time seris dist into temp---------------------------------------------------------------------------------*/
//			//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
//			temp.push_back(tempOriginalTimeSeriesPair);
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*-------------------------------------------------Clear memory----------------------------------------------------------------------------------------------------------*/
//			//original_time_series_vector.clear();
//			//original_time_series_vector.shrink_to_fit();
//			//delete[] original_time_series;
//			//original_time_series = nullptr;
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//		/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
//			/*====================================================================== Check Leaf Node, insert APCA point into queue===============================================================================*/
//			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//
//			/*................................................*/
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS / 2 == input_argument.point_dimension && m_temp_queue_top.p_rtree_node->m_count > 0);
//#endif
//			/*................................................*/
//
//			typename APCA_QUAL::APCA QProjection;
//			APCA_QUAL::initialAPCA(QProjection, input_argument.point_dimension);
//			vector<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_apla_projection_vector(input_argument.point_dimension);
//
//			/*-----------------------------------------------------------------------queue: insert leaf node. compute distance LB--------------------------------------------------------------------------------------------*/
//			f_temp_APLA_Pair.p_rtree_node = nullptr;
//			//cout << "!!!!!!!!!!!!!!!!: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//
//			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {// scan sub node in Rtree
//
//				/*................................................*/
//#ifdef _DEBUG
//				//cout << "!!!: " << m_temp_queue_top.p_rtree_node->m_branch[i].m_data << endl;
//				//g_n_account_leaf_node++;
//#endif
//				/*................................................*/
//
//				/*-----------------------------------------get ACPA(origitnal time series) point & id----------------------------------------------------*/
//				f_temp_APLA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);// get APCA(original time series) id
//				//cout << tempAPCAPair.id_originalTimeSeries << endl;
//				//f_temp_APLA_Pair.approximation_pointer = &apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)];// pointer to APCA point
//				/*----------------------------------------------------------------------------------------------------------------------------------------*/
//				//cout << "KNN : data point ID = " << tempAPCAPair.id_originalTimeSeries << endl;
//				//cout << tempAPCAPair.APCAValue << endl;
//
//				switch (input_argument.representation_option) {
//				case 1://MSPLA
//				case 6: //ICDE07
//				case 8: {//Initial 200706
//					//if (multi_y_projection_argument[f_temp_APLA_Pair.original_time_series_id].is_y_projection) {
//						//f_temp_APLA_Pair.d_dist = APLA::get_distance_apca_LB(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_apla_projection_vector));
//					//}
//					//else{
//					
//					/*~~~~~~~~~~~~~~~~~~~ 210602  ~~~~~~~~~~~~~~~~~~~*/
//					/*-------210511 PAA APCA version has projection of query time series -----------*/
//					f_temp_APLA_Pair.d_dist = APLA::get_distance_LB(apla_array_vector[f_temp_APLA_Pair.original_time_series_id], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[f_temp_APLA_Pair.original_time_series_id], query_apla_projection_vector));
//					/*------------------------------------------------------------------------------*/
//
//					/*-------210511 PAA APCA version has projection of query time series -----------*/
//					/*TOOL::read_normalized_multi_time_series(data_source, f_temp_APLA_Pair.original_time_series_id, original_time_series_vector);
//					long double pruning_power_distance = 0;
//					f_temp_APLA_Pair.d_dist = APLA::get_distance_SAPLA(query_time_series_vector, original_time_series_vector, apla_array_vector[input_argument.query_time_series_id], apla_array_vector[f_temp_APLA_Pair.original_time_series_id], pruning_power_distance);
//					input_argument.IO_cost += pruning_power_distance;*/
//					/*------------------------------------------------------------------------------*/
//					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//					//200109 AE
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_AE(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//PLA version
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB_pla(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//double test_reconstruct_dist = f_temp_APLA_Pair.d_dist;
//					//191206 APLA distance LB speed
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB_pla_speed(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//f_temp_APLA_Pair.d_dist = APLA::get_apla_endpoint_distance(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_lower_bound(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
////#ifdef _DEBUG
//					/*==================*/
//					//assert(float(test_reconstruct_dist) == float(f_temp_APLA_Pair.d_dist));
//					//vector<DataType> test_normalized_time_series_vector;
//					//vector<DataType> reconstruct_time_series_vector(query_time_series_vector.size(), INF);
//					//APLA::getAPLAReconstructSeries(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], reconstruct_time_series_vector);
//
//				//	TOOL::read_normalized_multi_time_series(data_source, f_temp_APLA_Pair.original_time_series_id, test_normalized_time_series_vector);
//					//double 
//					//double  test_sum_deviation = PLA_QUAL::getPLASumDeviation(test_normalized_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//double  test_query_deviaiton = PLA_QUAL::getPLASumDeviation(query_time_series_vector, query_linked_list);
//					//cout <<"AAAAAAAAAAAAAA:::" << f_temp_APLA_Pair.d_dist << endl;
//				//	double test_euc = TOOL::distanceEUC(query_time_series_vector, test_normalized_time_series_vector);
//					//assert(f_temp_APLA_Pair.d_dist <= test_euc);
//					//if (f_temp_APLA_Pair.original_time_series_id != 19){
//						//cout << test_sum_deviation <<"    <    "<< f_temp_APLA_Pair.d_dist << " <  "<< test_euc <<endl;
//						//assert(test_sum_deviation <= f_temp_APLA_Pair.d_dist);
//						//if (test_sum_deviation < f_temp_APLA_Pair.d_dist) {
//						//if (test_sum_deviation < f_temp_APLA_Pair.d_dist* 2) {
//							//f_temp_APLA_Pair.d_dist = fabs(f_temp_APLA_Pair.d_dist - test_sum_deviation - test_query_deviaiton);
//						//if(test_sum_deviation  > 10)
//						//f_temp_APLA_Pair.d_dist /= (100);
//						//else {
//							//f_temp_APLA_Pair.d_dist -= test_sum_deviation;
//						//}
//						//}
//
//						/*	for (auto&& au : test_normalized_time_series_vector) {
//								cout << au << ",";
//							}
//							cout << "========================================================================="<<endl;
//							for (auto&& au : reconstruct_time_series_vector) {
//								cout << au << ",";
//							}
//							cout << endl;*/
//							//assert(0);
//						//}
//
//					   // assert(f_temp_APLA_Pair.d_dist <= test_euc);
//						//if (f_temp_APLA_Pair.d_dist > test_euc) {
//						//	cout << test_sum_deviation << "    <    " << f_temp_APLA_Pair.d_dist << " <  " << test_euc << endl;
//						//	data_source.bigger_account++;
//						//		for (auto&& au : test_normalized_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << "\n========================================================================="<<endl;
//						//	for (auto&& au : reconstruct_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;
//
//						//	/*for (auto&& au : query_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;
//						//	cout << "##########################################\n";
//						//	for (auto&& au : test_normalized_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;*/
//						//	assert(0);
//						//}
//
//					//}
//					//test_normalized_time_series_vector.clear();
//					//test_normalized_time_series_vector.shrink_to_fit();
//////#endif
//					//}
//					break;
//				}
//				case 2://PLA
//					f_temp_APLA_Pair.d_dist = PLA_QUAL::getPLADistance(input_argument.point_multi_single_length, pla_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], PLA_query, f_temp_APLA_Pair.d_dist);
//					break;
//				case 3://APCA
//				case 4://PAA
//				case 7: {//PAALM
//					//f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceAE(query_time_series_vector, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceLB(APCA_KNN_QUAL::QAPCAProjection(g_query_time_series, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], QProjection), apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					/*============================================================================================================*/
//					/*vector<double> query_reconstruct_time_series;
//					vector<double> reconstruct_time_series;
//
//					APCA_KNN_QUAL::get_APCA_reconstruction(query_APCA, query_reconstruct_time_series);
//					APCA_KNN_QUAL::get_APCA_reconstruction(apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], reconstruct_time_series);
//					f_temp_APLA_Pair.d_dist = TOOL::distanceEUC(query_reconstruct_time_series, reconstruct_time_series);
//
//					query_reconstruct_time_series.clear();
//					query_reconstruct_time_series.shrink_to_fit();
//					reconstruct_time_series.clear();
//					reconstruct_time_series.shrink_to_fit();*/
//					/*================================================================================================================================*/
//					break;
//				}
//				case 9: {//SAX
//					f_temp_APLA_Pair.d_dist = sax.distance_LB_SAX(query_APCA, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@      " << f_temp_APLA_Pair.d_dist << endl;
//					break;
//				}
//				case 5:
//					assert(0);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//				queue.push(f_temp_APLA_Pair);
//			}
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
//			//	//cout << tempOriginalTimeSeriesPair.key;
//			/*------------------------------------------------------------------------Clear memory-------------------------------------------------------------------------------*/
//			APCA_QUAL::deleteAPCA(QProjection);
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS >> 1 == input_argument.point_dimension && input_argument.point_multi_single_length == input_argument.time_series_length);
//#endif
//			
//			//is internal Node
//			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			//APCA_KNN_QUAL::APCA_NODE_PAIR tempApcaPair;
//			RTREE_NODE_PAIR temp_apla_pair;
//			typename APCA_KNN_QUAL::REGION fs_region_G;
//			typename APCA_KNN_QUAL::initialREGION(fs_region_G, input_argument.point_dimension);
//			////cout << "regionNum = " << G.regionNum << endl;
//			//cout << "================: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
//				//cout << "===: " << m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_data << endl;
//				switch (input_argument.representation_option) {
//				case 1://MSPLA
//				case 6://ICDE07
//				case 8://Initial 200706
//					// APCA paper, min id == max id
//					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(reconstruct_query_time_series_vector, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//
//
//					//original min id < max id
//					//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//					//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//					//queue.push(temp_apla_pair);
//					//temp_apla_pair.approximation_pointer = nullptr;
//					break;
//				case 2://PLA
//					temp_apla_pair.d_dist = PLA_QUAL::getPLAMBRDistance(input_argument, m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, PLA_query, temp_apla_pair.d_dist);
//					break;
//				case 3://APCA
//				case 4://PAA
//				case 7://PAALM
//				case 9://SAX
//					// APCA paper, min id == max id
//					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//original min id < max id
//					//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//					//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//					//queue.push(temp_apla_pair);
//					//tempApcaPair.approximation_pointer = nullptr;
//					break;
//				case 5://CHEBY
//					assert(0);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//				temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//				queue.push(temp_apla_pair);
//			}
//			APCA_KNN_QUAL::deleteREGION(fs_region_G);
//		}
//		else {
//			assert(0);
//		}
//	}
//	/*=====================================================================================================================================================================================================================================================================================================================*/
//
//	/*========================Evaluation==========================================*/
//	input_argument.knn_total_time = TOOL::recordFinishTime(TOOL::time_record[2]);
//	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//	//assert(input_argument.pruning_power > 0);
//	input_argument.pruning_power = 1;//200929
//	/*============================================================================*/
//
//	/*................................................................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF);
//	cout << "??????????????????     APLA KNN failed    ??????????????????" << endl;
//	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	/*ofstream outfile(input_argument.write_file_name[16] + ".txt", ios::app);
//	assert(outfile.is_open());
//	outfile << "??????????????????     APLA KNN failed    ??????????????????   " << PAA_or_APCA << endl;
//	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
//	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	outfile.close();*/
//	TOOL::printInputArgument(input_argument);
//#endif
//	/*................................................................................................................................*/
//
//	//assert(0);
//	priority_queue<RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>>().swap(queue);
//	temp.clear();
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//
//	//f_temp_APLA_Pair.approximation_pointer = nullptr;
//	//f_APLA_Root.approximation_pointer = nullptr;
//
//	for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//		//cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//		result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
//	}
//	assert(result_set.size() <= input_argument.K);
//	result.clear();
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//
//	return result_set;
//}
//
//// Improve Index
////************************************
//// Method:all_knn_multi_210511
//// Qualifier: KNN structure for PAA, APCA, APLA, ICDE07
//// Input:
//// Output:
//// Notice: Evaluate SAPLA distance. Add input_argument.IO_cost++; in leaf node computation. change KNN time computation
//// date: 210511
//// author:
////************************************
//TEMPLATE
//template<typename T, typename Y>
//std::multiset<pair<double, int>> MULTI::all_knn_multi_210511(typename TOOL::DATA_SOURCE& const data_source, typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<Y>& const multi_y_projection_argument, DataType*& g_query_time_series, const RTREE& const apcaRTree, T& const approximation_array) {
//	/*..........................................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.file_id != INF && input_argument.K <= input_argument.point_number && input_argument.time_series_length != INF && input_argument.point_number != INF && input_argument.arity_d != INF && input_argument.K != INF);
//	assert(input_argument.time_series_length == input_argument.point_multi_single_length);
//	assert(input_argument.point_dimension == input_argument.point_multi_single_dimension);
//
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//		cout << "MSPLA(SAPLA KNN: \n";
//		break;
//	case 2://PLA
//		cout << "PLA KNN: \n";
//		break;
//	case 3://APCA
//		cout << "APCA KNN: \n";
//		break;
//	case 4://PAA
//		cout << "PAA KNN: \n";
//		break;
//	case 5://Chebyshev
//		cout << "CHEBY KNN: \n";
//		break;
//	case 6://ICDE07
//		cout << "ICDE07 KNN: \n";
//		break;
//	case 7://PAALM
//		cout << "PAALM KNN: \n";
//		break;
//	case 8://initial 200706
//		cout << "Initial 200706 KNN: \n";
//		break;
//	case 9://SAX
//		cout << "SAX KNN: \n";
//		break;
//	default:
//		assert(0);
//		break;
//	}
//#endif
//	/*..........................................................................................................*/
//
//	/*-------------------------      Evaluation  ---------------------------------*/
//	input_argument.knn_total_time = 0.0;
//	input_argument.knn_total_time_has_IO = 0.0;//210606
//	input_argument.IO_cost = 0;// measure I/O cost
//	input_argument.pruning_power = 0;
//	/*-------------------------------------------------------------------------------*/
//
//	/*-------------------------------------change pointer array to vector array--------------------------------*/
//	//int multi_to_single_series_length = input_argument.time_series_length * input_argument.arity_d;
//	int multi_to_single_series_length = input_argument.time_series_length;
//	vector<DataType> query_time_series_vector;
//	query_time_series_vector.resize(multi_to_single_series_length, INF);
//	std::copy_n(g_query_time_series, multi_to_single_series_length, query_time_series_vector.begin());
//
//	vector<DataType> reconstruct_query_time_series_vector;
//	reconstruct_query_time_series_vector.resize(multi_to_single_series_length, INF);
//	/*---------------------------------------------------------------------------------------------------------*/
//
//
//	/*=======================================Evaluation: plot & bar chart=======================================================*/
//#ifdef _DEBUG
//	//g_n_account_apca_point = 0;
//	input_argument.navigate_index_time = 0.0;// navigate time
//	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0.0;// distance euclidean time
//#endif
//	/*===========================================================================================================================*/
//
//	/**/
//	typename APCA_QUAL::APCA query_APCA;
//	SaxQuantizer::SAX sax(input_argument.point_dimension);
//	//APCA_QUAL::getAPCAPoint(g_query_time_series, input_argument.time_series_length, input_argument.point_dimension, query_APCA);
//	typename PLA_QUAL::PLA PLA_query;
//	/**/
//
//	/*========================================!!!!!!!!!!!   MSPLA & ICDE07 Variable construct=================================================================*/
//
//	DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>();
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//	//APLA_ICDE07<DataType>::getAPLA_ICDE07(input_argument, normalized_series_vector, apla_array_vector[point_id]);
//	case 6: //ICDE07
//	case 8: {//Initial 200706
//		//typename TOOL::Y_PROJECTION_ARGUMENT y_projection_argument(input_argument.initial_N);
//		//DoublyLinkedList<APLA::AREA_COEFFICIENT> all_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//		//DoublyLinkedList<APLA::AREA_COEFFICIENT> cluster_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//		query_linked_list.copy(apla_array_vector[input_argument.query_time_series_id]);
//		//APLA::get_APLA_point(input_argument, g_query_time_series, y_projection_argument, all_linked_list, cluster_linked_list, query_linked_list);//191129
//
//		APLA::getAPLAReconstructSeries(query_linked_list, reconstruct_query_time_series_vector);
//
//		/*......................................*/
//#ifdef _DEBUG
//		assert(query_linked_list.size() == apla_array_vector[input_argument.query_time_series_id].size());
//		for (int segment_id = 0; segment_id < query_linked_list.size(); segment_id++) {
//			const auto& const query_segment = query_linked_list[segment_id];
//			const auto& const array_segment = apla_array_vector[input_argument.query_time_series_id][segment_id];
//			//cout <<"!!!!!!!!!!!!" <<query_segment.right_endpoint << ", "<<array_segment.right_endpoint << endl;
//			assert(query_segment.right_endpoint == array_segment.right_endpoint);
//			assert(query_segment.right_endpoint != INF && query_segment.right_endpoint >= 0 && query_segment.right_endpoint < input_argument.point_multi_single_length);
//		}
//#endif
//		/*......................................*/
//
//		break;
//	}
//	case 2: {//PLA
//		PLA_QUAL::initialPLA(PLA_query, input_argument.point_multi_single_dimension);//????????180918 , this multi ot single has problem
//		PLA_QUAL::getPLA(input_argument.point_multi_single_length, input_argument.point_multi_single_dimension, g_query_time_series, PLA_query);
//		break;
//	}
//	case 3://APCA
//	case 4://PAA
//		break;
//	case 5://CHEBY
//		assert(0);
//		break;
//	case 7://PAALM
//		break;
//	case 9: {//SAX
//		APCA_QUAL::initialAPCA(query_APCA, input_argument.point_dimension);
//		sax.get_SAX(query_time_series_vector, input_argument.point_dimension, query_APCA);
//		break;
//	}
//	default:
//		assert(0);
//		break;
//	}
//	/*=======================================================================================================================================================*/
//	/*========================================!!!!!!!!!!!   APLA & ICDE07 Variable construct=================================================================*/
//	int i = NULL, j = NULL;
//	priority_queue<RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>> queue;// <Rtree node, distance>
//	RTREE_NODE_PAIR f_APLA_Root, f_temp_APLA_Pair;// distance, APCA point(original time series) id, APCA point pointer, Rtree sub node pointer
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> temp; // <dist, original time series id>
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> result;// <dist, original time series id>
//	std::multiset<pair<double, int>> result_set;//191204 result
//	typename APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair; // <dist, original time series id>
//	RTREE_NODE_PAIR m_temp_queue_top;//
//	vector<DataType> original_time_series_vector;
//	/*==========================================================================================================================================*/
//	/*-------------------------Evaluation Time---------------------------------*/
//	TOOL::recordStartTime(TOOL::time_record[14]); //whole time
//	TOOL::recordStartTime(TOOL::time_record[16]); //whole time has IO
//	TOOL::recordStartTime(TOOL::time_record[2]);//knn time
//	TOOL::recordStartTime(TOOL::time_record[15]);//knn has IO time 210606
//	/*--------------------------------------------------------------------------------*/
//	/*========================================== Variable Initial ===============================================================*/
//	f_APLA_Root.p_rtree_node = apcaRTree.m_root;
//	f_APLA_Root.d_dist = 0;
//	queue.push(f_APLA_Root);
//	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
//	//printf("<///////**    KNN Begin   **////////>\n");
//	/*===========================================================================================================================*/
//
//	/*========================================================================================Begin K-NN=======================================================================================================================================*/
//	while (!queue.empty()) {
//
//		m_temp_queue_top = queue.top();
//
//		//cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator: ";
//		//for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
//		//cout << endl;
//
//		/*===========================================================================Scan temp to get ruslt=========================================================================================================*/
//		for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
//			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
//			if (plist->d_dist <= m_temp_queue_top.d_dist) {
//				//cout << "           <= " << endl;
//				result.push_back(*plist);
//				plist = temp.erase(plist);
//			}
//			else
//				plist++;
//
//			if (input_argument.K == result.size()) {
//				/*-------------------------Evaluation ---------------------------------*/
//				input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
//				input_argument.knn_total_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[15]);//210606
//				input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//				input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);//210606
//				assert(input_argument.IO_cost != 0 && input_argument.point_number != INF);
//				input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
//				assert(input_argument.pruning_power > 0);//&& input_argument.pruning_power < 1);//210603
//				/*---------------------------------------------------------------------*/
//
//				/*..........................................................................................................*/
//#ifdef _DEBUG
//				assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF && input_argument.knn_total_time_has_IO != INF);
//				/*-------------------------------------------- Print Result ---------------------------------*/
//				cout << "!!!!!!!!!!!!! KNN Find result !!!!!!!!!!!!!!!!!!!!!!!!!!!!   result list size: " << result.size() << endl;
//
//				cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
//
//				cout << "Total KNN time has IO : " << input_argument.knn_total_time_has_IO << " us" << endl;
//
//				cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;
//
//				cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;
//
//				cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;
//
//				cout << "pruning power: " << input_argument.pruning_power << endl;
//
//				cout << "I/O cost: " << input_argument.IO_cost << endl;
//
//				result.sort([](const APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& first, const  APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist;});//small to big
//#endif
//				/*..........................................................................................................*/
//
//				for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//					//cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//					result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
//				}
//				assert(result_set.size() == input_argument.K);
//				//cout << endl;
//				/*--------------------------------------------------------------------------------------------*/
//				/*--------------------------------------     Clear memory    -----------------------------------------------*/
//				//f_temp_APLA_Pair.approximation_pointer = nullptr;
//				//f_APLA_Root.approximation_pointer = nullptr;
//				priority_queue< RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>>().swap(queue);
//
//				temp.clear();
//				result.clear();
//				list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//				//list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//				/*---------------------------------------------------------------------------------------------------------*/
//				return result_set;
//			}
//		}
//		/*==========================================================================================================================================================================*/
//		/*======================================================================Pop top node ine queue==============================================================================*/
//		queue.pop();
//		/*==========================================================================================================================================================================*/
//		/*====================================================================== Check APLA point ==================================================================================*/
//		if (m_temp_queue_top.p_rtree_node == nullptr) { //is Approximation data point
//
//			//cout << "    queue.top is data point\n";
//			//DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
//			//vector<DataType> original_time_series_vector;
//
//			/*------------------------------------------------get original time series from database (file)----------------------------------------------*/
//			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
//			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
//			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;//get top queue time sereis id.
//			//getFileStreamByID(file_name, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			// get original time sires from database(txt file)
//
//			/*..........................................................................................................*/
//#ifdef _DEBUG
//			assert(input_argument.time_series_length == data_source.multi_single_time_series_length && input_argument.read_file_name == data_source.read_file_address);
//#endif
//			/*..........................................................................................................*/
//			input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);//210603
//			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//			TOOL::read_normalized_multi_time_series(data_source, m_temp_queue_top.original_time_series_id, original_time_series_vector);
//			TOOL::recordStartTime(TOOL::time_record[14]); //whole time
//			TOOL::recordStartTime(TOOL::time_record[2]);//knn time
//			//			//already normalized
//			//			TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series_vector);
//			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*===============Evaluation==============*/
//			input_argument.IO_cost++;
//			/*=======================================*/
//
//			/*------------------------------------------------compute Euclidean distance from query times and candidate original time series-----------------------------------------*/
//			tempOriginalTimeSeriesPair.d_dist = TOOL::distanceEUC(query_time_series_vector, original_time_series_vector);
//			//#ifdef _DEBUG
//			//			copy_n(original_time_series_vector.begin(), original_time_series_vector.size(), original_time_series);
//			//			double test_distance = APCA_KNN_QUAL::distanceEUC(g_query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);
//			//			assert(test_distance == tempOriginalTimeSeriesPair.d_dist);
//			//#endif
//						/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//						/*------------------------------------------------temp: insert time seris dist into temp---------------------------------------------------------------------------------*/
//						//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
//			temp.push_back(tempOriginalTimeSeriesPair);
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*-------------------------------------------------Clear memory----------------------------------------------------------------------------------------------------------*/
//			//original_time_series_vector.clear();
//			//original_time_series_vector.shrink_to_fit();
//			//delete[] original_time_series;
//			//original_time_series = nullptr;
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//		/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
//			/*====================================================================== Check Leaf Node, insert APCA point into queue===============================================================================*/
//			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//
//			/*................................................*/
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS / 2 == input_argument.point_dimension && m_temp_queue_top.p_rtree_node->m_count > 0);
//#endif
//			/*................................................*/
//
//			typename APCA_QUAL::APCA QProjection;
//			APCA_QUAL::initialAPCA(QProjection, input_argument.point_dimension);
//			vector<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_apla_projection_vector(input_argument.point_dimension);
//
//			/*-----------------------------------------------------------------------queue: insert leaf node. compute distance LB--------------------------------------------------------------------------------------------*/
//			f_temp_APLA_Pair.p_rtree_node = nullptr;
//			//cout << "!!!!!!!!!!!!!!!!: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//
//			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {// scan sub node in Rtree
//
//				/*................................................*/
//#ifdef _DEBUG
//				//cout << "!!!: " << m_temp_queue_top.p_rtree_node->m_branch[i].m_data << endl;
//				//g_n_account_leaf_node++;
//#endif
//				/*................................................*/
//
//				/*-----------------------------------------get ACPA(origitnal time series) point & id----------------------------------------------------*/
//				f_temp_APLA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);// get APCA(original time series) id
//				//cout << tempAPCAPair.id_originalTimeSeries << endl;
//				//f_temp_APLA_Pair.approximation_pointer = &apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)];// pointer to APCA point
//				/*----------------------------------------------------------------------------------------------------------------------------------------*/
//				//cout << "KNN : data point ID = " << tempAPCAPair.id_originalTimeSeries << endl;
//				//cout << tempAPCAPair.APCAValue << endl;
//
//				switch (input_argument.representation_option) {
//				case 1://MSPLA
//				case 6: //ICDE07
//				case 8: {//Initial 200706
//					//if (multi_y_projection_argument[f_temp_APLA_Pair.original_time_series_id].is_y_projection) {
//						//f_temp_APLA_Pair.d_dist = APLA::get_distance_apca_LB(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_apla_projection_vector));
//					//}
//					//else{
//
//					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       210511    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//					/*-------210511 PAA APCA version has projection of query time series -----------*/
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB(apla_array_vector[f_temp_APLA_Pair.original_time_series_id], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[f_temp_APLA_Pair.original_time_series_id], query_apla_projection_vector));
//					/*------------------------------------------------------------------------------*/
//
//					/*-------210511 PAA APCA version has projection of query time series -----------*/
//					input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);//210603
//					input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//					TOOL::read_normalized_multi_time_series(data_source, f_temp_APLA_Pair.original_time_series_id, original_time_series_vector);
//					TOOL::recordStartTime(TOOL::time_record[14]); //whole time
//					TOOL::recordStartTime(TOOL::time_record[2]);//knn time
//
//					long double pruning_power_distance = 0;//210603
//					f_temp_APLA_Pair.d_dist = APLA::get_distance_SAPLA(query_time_series_vector, original_time_series_vector, apla_array_vector[input_argument.query_time_series_id], apla_array_vector[f_temp_APLA_Pair.original_time_series_id], pruning_power_distance);
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_SAPLA(query_time_series_vector, f_temp_APLA_Pair.original_time_series_id, data_source.read_file_address_vector[0], apla_array_vector[input_argument.query_time_series_id], apla_array_vector[f_temp_APLA_Pair.original_time_series_id], pruning_power_distance);
//					/*===============Evaluation==============*/
//					input_argument.IO_cost += pruning_power_distance;
//					/*=======================================*/
//					
//					/*------------------------------------------------------------------------------*/
//					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//
//					//200109 AE
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_AE(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//PLA version
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB_pla(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//double test_reconstruct_dist = f_temp_APLA_Pair.d_dist;
//					//191206 APLA distance LB speed
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB_pla_speed(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//f_temp_APLA_Pair.d_dist = APLA::get_apla_endpoint_distance(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_lower_bound(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
////#ifdef _DEBUG
//					/*==================*/
//					//assert(float(test_reconstruct_dist) == float(f_temp_APLA_Pair.d_dist));
//					//vector<DataType> test_normalized_time_series_vector;
//					//vector<DataType> reconstruct_time_series_vector(query_time_series_vector.size(), INF);
//					//APLA::getAPLAReconstructSeries(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], reconstruct_time_series_vector);
//
//				//	TOOL::read_normalized_multi_time_series(data_source, f_temp_APLA_Pair.original_time_series_id, test_normalized_time_series_vector);
//					//double 
//					//double  test_sum_deviation = PLA_QUAL::getPLASumDeviation(test_normalized_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//double  test_query_deviaiton = PLA_QUAL::getPLASumDeviation(query_time_series_vector, query_linked_list);
//					//cout <<"AAAAAAAAAAAAAA:::" << f_temp_APLA_Pair.d_dist << endl;
//				//	double test_euc = TOOL::distanceEUC(query_time_series_vector, test_normalized_time_series_vector);
//					//assert(f_temp_APLA_Pair.d_dist <= test_euc);
//					//if (f_temp_APLA_Pair.original_time_series_id != 19){
//						//cout << test_sum_deviation <<"    <    "<< f_temp_APLA_Pair.d_dist << " <  "<< test_euc <<endl;
//						//assert(test_sum_deviation <= f_temp_APLA_Pair.d_dist);
//						//if (test_sum_deviation < f_temp_APLA_Pair.d_dist) {
//						//if (test_sum_deviation < f_temp_APLA_Pair.d_dist* 2) {
//							//f_temp_APLA_Pair.d_dist = fabs(f_temp_APLA_Pair.d_dist - test_sum_deviation - test_query_deviaiton);
//						//if(test_sum_deviation  > 10)
//						//f_temp_APLA_Pair.d_dist /= (100);
//						//else {
//							//f_temp_APLA_Pair.d_dist -= test_sum_deviation;
//						//}
//						//}
//
//						/*	for (auto&& au : test_normalized_time_series_vector) {
//								cout << au << ",";
//							}
//							cout << "========================================================================="<<endl;
//							for (auto&& au : reconstruct_time_series_vector) {
//								cout << au << ",";
//							}
//							cout << endl;*/
//							//assert(0);
//						//}
//
//					   // assert(f_temp_APLA_Pair.d_dist <= test_euc);
//						//if (f_temp_APLA_Pair.d_dist > test_euc) {
//						//	cout << test_sum_deviation << "    <    " << f_temp_APLA_Pair.d_dist << " <  " << test_euc << endl;
//						//	data_source.bigger_account++;
//						//		for (auto&& au : test_normalized_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << "\n========================================================================="<<endl;
//						//	for (auto&& au : reconstruct_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;
//
//						//	/*for (auto&& au : query_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;
//						//	cout << "##########################################\n";
//						//	for (auto&& au : test_normalized_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;*/
//						//	assert(0);
//						//}
//
//					//}
//					//test_normalized_time_series_vector.clear();
//					//test_normalized_time_series_vector.shrink_to_fit();
//////#endif
//					//}
//					break;
//				}
//				case 2://PLA
//					f_temp_APLA_Pair.d_dist = PLA_QUAL::getPLADistance(input_argument.point_multi_single_length, pla_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], PLA_query, f_temp_APLA_Pair.d_dist);
//					break;
//				case 3://APCA
//				case 4://PAA
//				case 7: {//PAALM
//					//f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceAE(query_time_series_vector, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceLB(APCA_KNN_QUAL::QAPCAProjection(g_query_time_series, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], QProjection), apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					/*============================================================================================================*/
//					/*vector<double> query_reconstruct_time_series;
//					vector<double> reconstruct_time_series;
//
//					APCA_KNN_QUAL::get_APCA_reconstruction(query_APCA, query_reconstruct_time_series);
//					APCA_KNN_QUAL::get_APCA_reconstruction(apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], reconstruct_time_series);
//					f_temp_APLA_Pair.d_dist = TOOL::distanceEUC(query_reconstruct_time_series, reconstruct_time_series);
//
//					query_reconstruct_time_series.clear();
//					query_reconstruct_time_series.shrink_to_fit();
//					reconstruct_time_series.clear();
//					reconstruct_time_series.shrink_to_fit();*/
//					/*================================================================================================================================*/
//					break;
//				}
//				case 9: {//SAX
//					f_temp_APLA_Pair.d_dist = sax.distance_LB_SAX(query_APCA, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@      " << f_temp_APLA_Pair.d_dist << endl;
//					break;
//				}
//				case 5:
//					assert(0);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//				queue.push(f_temp_APLA_Pair);
//			}
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
//			//	//cout << tempOriginalTimeSeriesPair.key;
//			/*------------------------------------------------------------------------Clear memory-------------------------------------------------------------------------------*/
//			APCA_QUAL::deleteAPCA(QProjection);
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS >> 1 == input_argument.point_dimension && input_argument.point_multi_single_length == input_argument.time_series_length);
//#endif
//
//			//is internal Node
//			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			//APCA_KNN_QUAL::APCA_NODE_PAIR tempApcaPair;
//			RTREE_NODE_PAIR temp_apla_pair;
//			typename APCA_KNN_QUAL::REGION fs_region_G;
//			typename APCA_KNN_QUAL::initialREGION(fs_region_G, input_argument.point_dimension);
//			////cout << "regionNum = " << G.regionNum << endl;
//			//cout << "================: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
//				//cout << "===: " << m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_data << endl;
//				switch (input_argument.representation_option) {
//				case 1://MSPLA
//				case 6://ICDE07
//				case 8://Initial 200706
//					// APCA paper, min id == max id
//					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(reconstruct_query_time_series_vector, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//
//
//					//original min id < max id
//					//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//					//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//					//queue.push(temp_apla_pair);
//					//temp_apla_pair.approximation_pointer = nullptr;
//					break;
//				case 2://PLA
//					temp_apla_pair.d_dist = PLA_QUAL::getPLAMBRDistance(input_argument, m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, PLA_query, temp_apla_pair.d_dist);
//					break;
//				case 3://APCA
//				case 4://PAA
//				case 7://PAALM
//				case 9://SAX
//					// APCA paper, min id == max id
//					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//original min id < max id
//					//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//					//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//					//queue.push(temp_apla_pair);
//					//tempApcaPair.approximation_pointer = nullptr;
//					break;
//				case 5://CHEBY
//					assert(0);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//				temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//				queue.push(temp_apla_pair);
//			}
//			APCA_KNN_QUAL::deleteREGION(fs_region_G);
//		}
//		else {
//			assert(0);
//		}
//	}
//	/*=====================================================================================================================================================================================================================================================================================================================*/
//
//	/*========================Evaluation==========================================*/
//	input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
//	input_argument.knn_total_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[15]);
//	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//	input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);
//	//assert(input_argument.pruning_power > 0);
//	input_argument.pruning_power = 1;//200929
//	/*============================================================================*/
//
//	/*................................................................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF);
//	cout << "??????????????????     APLA KNN failed    ??????????????????" << endl;
//	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	/*ofstream outfile(input_argument.write_file_name[16] + ".txt", ios::app);
//	assert(outfile.is_open());
//	outfile << "??????????????????     APLA KNN failed    ??????????????????   " << PAA_or_APCA << endl;
//	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
//	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	outfile.close();*/
//	TOOL::printInputArgument(input_argument);
//#endif
//	/*................................................................................................................................*/
//
//	//assert(0);
//	priority_queue<RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>>().swap(queue);
//	temp.clear();
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//
//	//f_temp_APLA_Pair.approximation_pointer = nullptr;
//	//f_APLA_Root.approximation_pointer = nullptr;
//
//	for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//		//cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//		result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
//	}
//	assert(result_set.size() <= input_argument.K);
//	result.clear();
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//
//	return result_set;
//}
//
//
////210618 RTree Partition. Evaluate SAPLA distance. Add input_argument.IO_cost++; in leaf node computation
//// Improve Index
////************************************
//// Method:all_knn_multi_210618
//// Qualifier: KNN structure for PAA, APCA, APLA, ICDE07
//// Input:
//// Output:
//// Notice: Evaluate SAPLA distance. Add input_argument.IO_cost++; in leaf node computation. change KNN time computation
//// date: 210618
//// author:
////************************************
//TEMPLATE
//template<typename T, typename Y, typename U>
//std::multiset<pair<double, int>> MULTI::all_knn_multi_210618(typename TOOL::DATA_SOURCE& const data_source, typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<Y>& const multi_y_projection_argument, DataType*& g_query_time_series, const U& const apcaRTree, T& const approximation_array) {
//	/*..........................................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.file_id != INF && input_argument.K <= input_argument.point_number && input_argument.time_series_length != INF && input_argument.point_number != INF && input_argument.arity_d != INF && input_argument.K != INF);
//	assert(input_argument.time_series_length == input_argument.point_multi_single_length);
//	assert(input_argument.point_dimension == input_argument.point_multi_single_dimension);
//
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//		cout << "MSPLA(SAPLA KNN: \n";
//		break;
//	case 2://PLA
//		cout << "PLA KNN: \n";
//		break;
//	case 3://APCA
//		cout << "APCA KNN: \n";
//		break;
//	case 4://PAA
//		cout << "PAA KNN: \n";
//		break;
//	case 5://Chebyshev
//		cout << "CHEBY KNN: \n";
//		break;
//	case 6://ICDE07
//		cout << "ICDE07 KNN: \n";
//		break;
//	case 7://PAALM
//		cout << "PAALM KNN: \n";
//		break;
//	case 8://initial 200706
//		cout << "Initial 200706 KNN: \n";
//		break;
//	case 9://SAX
//		cout << "SAX KNN: \n";
//		break;
//	default:
//		assert(0);
//		break;
//	}
//#endif
//	/*..........................................................................................................*/
//
//	/*-------------------------      Evaluation  ---------------------------------*/
//	input_argument.knn_total_time = 0.0;
//	input_argument.knn_total_time_has_IO = 0.0;//210606
//	input_argument.IO_cost = 0;// measure I/O cost
//	input_argument.pruning_power = 0;
//	/*-------------------------------------------------------------------------------*/
//
//	/*-------------------------------------change pointer array to vector array--------------------------------*/
//	//int multi_to_single_series_length = input_argument.time_series_length * input_argument.arity_d;
//	int multi_to_single_series_length = input_argument.time_series_length;
//	vector<DataType> query_time_series_vector;
//	query_time_series_vector.resize(multi_to_single_series_length, INF);
//	std::copy_n(g_query_time_series, multi_to_single_series_length, query_time_series_vector.begin());
//
//	vector<DataType> reconstruct_query_time_series_vector;
//	reconstruct_query_time_series_vector.resize(multi_to_single_series_length, INF);
//	/*---------------------------------------------------------------------------------------------------------*/
//
//
//	/*=======================================Evaluation: plot & bar chart=======================================================*/
//#ifdef _DEBUG
//	//g_n_account_apca_point = 0;
//	input_argument.navigate_index_time = 0.0;// navigate time
//	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0.0;// distance euclidean time
//#endif
//	/*===========================================================================================================================*/
//
//	/**/
//	typename APCA_QUAL::APCA query_APCA;
//	SaxQuantizer::SAX sax(input_argument.point_dimension);
//	//APCA_QUAL::getAPCAPoint(g_query_time_series, input_argument.time_series_length, input_argument.point_dimension, query_APCA);
//	typename PLA_QUAL::PLA PLA_query;
//	/**/
//
//	/*========================================!!!!!!!!!!!   MSPLA & ICDE07 Variable construct=================================================================*/
//
//	DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>();
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//	//APLA_ICDE07<DataType>::getAPLA_ICDE07(input_argument, normalized_series_vector, apla_array_vector[point_id]);
//	case 6: //ICDE07
//	case 8: {//Initial 200706
//		//typename TOOL::Y_PROJECTION_ARGUMENT y_projection_argument(input_argument.initial_N);
//		//DoublyLinkedList<APLA::AREA_COEFFICIENT> all_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//		//DoublyLinkedList<APLA::AREA_COEFFICIENT> cluster_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//		query_linked_list.copy(apla_array_vector[input_argument.query_time_series_id]);
//		//APLA::get_APLA_point(input_argument, g_query_time_series, y_projection_argument, all_linked_list, cluster_linked_list, query_linked_list);//191129
//
//		APLA::getAPLAReconstructSeries(query_linked_list, reconstruct_query_time_series_vector);
//
//		/*......................................*/
//#ifdef _DEBUG
//		assert(query_linked_list.size() == apla_array_vector[input_argument.query_time_series_id].size());
//		for (int segment_id = 0; segment_id < query_linked_list.size(); segment_id++) {
//			const auto& const query_segment = query_linked_list[segment_id];
//			const auto& const array_segment = apla_array_vector[input_argument.query_time_series_id][segment_id];
//			//cout <<"!!!!!!!!!!!!" <<query_segment.right_endpoint << ", "<<array_segment.right_endpoint << endl;
//			assert(query_segment.right_endpoint == array_segment.right_endpoint);
//			assert(query_segment.right_endpoint != INF && query_segment.right_endpoint >= 0 && query_segment.right_endpoint < input_argument.point_multi_single_length);
//		}
//#endif
//		/*......................................*/
//
//		break;
//	}
//	case 2: {//PLA
//		PLA_QUAL::initialPLA(PLA_query, input_argument.point_multi_single_dimension);//????????180918 , this multi ot single has problem
//		PLA_QUAL::getPLA(input_argument.point_multi_single_length, input_argument.point_multi_single_dimension, g_query_time_series, PLA_query);
//		break;
//	}
//	case 3://APCA
//	case 4://PAA
//		break;
//	case 5://CHEBY
//		assert(0);
//		break;
//	case 7://PAALM
//		break;
//	case 9: {//SAX
//		APCA_QUAL::initialAPCA(query_APCA, input_argument.point_dimension);
//		sax.get_SAX(query_time_series_vector, input_argument.point_dimension, query_APCA);
//		break;
//	}
//	default:
//		assert(0);
//		break;
//	}
//	/*=======================================================================================================================================================*/
//	/*========================================!!!!!!!!!!!   APLA & ICDE07 Variable construct=================================================================*/
//	int i = NULL, j = NULL;
//	priority_queue<RTREE_NODE_PAIR_PARTITION<double, int, U::Node>, vector<MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node>>, MULTI::priorityIncrement<MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node>>> queue;// <Rtree node, distance>
//	MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node> f_APLA_Root, f_temp_APLA_Pair;// distance, APCA point(original time series) id, APCA point pointer, Rtree sub node pointer
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> temp; // <dist, original time series id>
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> result;// <dist, original time series id>
//	std::multiset<pair<double, int>> result_set;//191204 result
//	typename APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair; // <dist, original time series id>
//	MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node> m_temp_queue_top;//
//	vector<DataType> original_time_series_vector;
//	/*==========================================================================================================================================*/
//	/*-------------------------Evaluation Time---------------------------------*/
//	TOOL::recordStartTime(TOOL::time_record[14]); //whole time
//	TOOL::recordStartTime(TOOL::time_record[16]); //whole time has IO
//	TOOL::recordStartTime(TOOL::time_record[2]);//knn time
//	TOOL::recordStartTime(TOOL::time_record[15]);//knn has IO time 210606
//	/*--------------------------------------------------------------------------------*/
//	/*========================================== Variable Initial ===============================================================*/
//	f_APLA_Root.p_rtree_node = apcaRTree.m_root;
//	f_APLA_Root.d_dist = 0;
//	queue.push(f_APLA_Root);
//	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
//	//printf("<///////**    KNN Begin   **////////>\n");
//	/*===========================================================================================================================*/
//
//	/*========================================================================================Begin K-NN=======================================================================================================================================*/
//	while (!queue.empty()) {
//
//		m_temp_queue_top = queue.top();
//
//		//cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator: ";
//		//for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
//		//cout << endl;
//
//		/*===========================================================================Scan temp to get ruslt=========================================================================================================*/
//		for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
//			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
//			if (plist->d_dist <= m_temp_queue_top.d_dist) {
//				//cout << "           <= " << endl;
//				result.push_back(*plist);
//				plist = temp.erase(plist);
//			}
//			else
//				plist++;
//
//			if (input_argument.K == result.size()) {
//				/*-------------------------Evaluation ---------------------------------*/
//				input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
//				input_argument.knn_total_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[15]);//210606
//				input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//				input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);//210606
//				assert(input_argument.IO_cost != 0 && input_argument.point_number != INF);
//				input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
//				assert(input_argument.pruning_power > 0);//&& input_argument.pruning_power < 1);//210603
//				/*---------------------------------------------------------------------*/
//
//				/*..........................................................................................................*/
//#ifdef _DEBUG
//				assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF && input_argument.knn_total_time_has_IO != INF);
//				/*-------------------------------------------- Print Result ---------------------------------*/
//				cout << "!!!!!!!!!!!!! KNN Find result !!!!!!!!!!!!!!!!!!!!!!!!!!!!   result list size: " << result.size() << endl;
//
//				cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
//
//				cout << "Total KNN time has IO : " << input_argument.knn_total_time_has_IO << " us" << endl;
//
//				cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;
//
//				cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;
//
//				cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;
//
//				cout << "pruning power: " << input_argument.pruning_power << endl;
//
//				cout << "I/O cost: " << input_argument.IO_cost << endl;
//
//				result.sort([](const APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& first, const  APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist; });//small to big
//#endif
//				/*..........................................................................................................*/
//
//				for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//					//cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//					result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
//				}
//				assert(result_set.size() == input_argument.K);
//				//cout << endl;
//				/*--------------------------------------------------------------------------------------------*/
//				/*--------------------------------------     Clear memory    -----------------------------------------------*/
//				//f_temp_APLA_Pair.approximation_pointer = nullptr;
//				//f_APLA_Root.approximation_pointer = nullptr;
//				priority_queue< MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node>, vector<MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node>>, MULTI::priorityIncrement<MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node>>>().swap(queue);
//
//				temp.clear();
//				result.clear();
//				list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//				//list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//				/*---------------------------------------------------------------------------------------------------------*/
//				return result_set;
//			}
//		}
//		/*==========================================================================================================================================================================*/
//		/*======================================================================Pop top node ine queue==============================================================================*/
//		queue.pop();
//		/*==========================================================================================================================================================================*/
//		/*====================================================================== Check APLA point ==================================================================================*/
//		if (m_temp_queue_top.p_rtree_node == nullptr) { //is Approximation data point
//
//			//cout << "    queue.top is data point\n";
//			//DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
//			//vector<DataType> original_time_series_vector;
//
//			/*------------------------------------------------get original time series from database (file)----------------------------------------------*/
//			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
//			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
//			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;//get top queue time sereis id.
//			//getFileStreamByID(file_name, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			// get original time sires from database(txt file)
//
//			/*..........................................................................................................*/
//#ifdef _DEBUG
//			assert(input_argument.time_series_length == data_source.multi_single_time_series_length && input_argument.read_file_name == data_source.read_file_address);
//#endif
//			/*..........................................................................................................*/
//			input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);//210603
//			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//			TOOL::read_normalized_multi_time_series(data_source, m_temp_queue_top.original_time_series_id, original_time_series_vector);
//			TOOL::recordStartTime(TOOL::time_record[14]); //whole time
//			TOOL::recordStartTime(TOOL::time_record[2]);//knn time
//			//			//already normalized
//			//			TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series_vector);
//			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*===============Evaluation==============*/
//			input_argument.IO_cost++;
//			/*=======================================*/
//
//			/*------------------------------------------------compute Euclidean distance from query times and candidate original time series-----------------------------------------*/
//			tempOriginalTimeSeriesPair.d_dist = TOOL::distanceEUC(query_time_series_vector, original_time_series_vector);
//			//#ifdef _DEBUG
//			//			copy_n(original_time_series_vector.begin(), original_time_series_vector.size(), original_time_series);
//			//			double test_distance = APCA_KNN_QUAL::distanceEUC(g_query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);
//			//			assert(test_distance == tempOriginalTimeSeriesPair.d_dist);
//			//#endif
//						/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//						/*------------------------------------------------temp: insert time seris dist into temp---------------------------------------------------------------------------------*/
//						//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
//			temp.push_back(tempOriginalTimeSeriesPair);
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*-------------------------------------------------Clear memory----------------------------------------------------------------------------------------------------------*/
//			//original_time_series_vector.clear();
//			//original_time_series_vector.shrink_to_fit();
//			//delete[] original_time_series;
//			//original_time_series = nullptr;
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//		/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
//			/*====================================================================== Check Leaf Node, insert APCA point into queue===============================================================================*/
//			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//
//			/*................................................*/
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS / 2 == input_argument.point_dimension && m_temp_queue_top.p_rtree_node->m_count > 0);
//#endif
//			/*................................................*/
//
//			typename APCA_QUAL::APCA QProjection;
//			APCA_QUAL::initialAPCA(QProjection, input_argument.point_dimension);
//			vector<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_apla_projection_vector(input_argument.point_dimension);
//
//			/*-----------------------------------------------------------------------queue: insert leaf node. compute distance LB--------------------------------------------------------------------------------------------*/
//			f_temp_APLA_Pair.p_rtree_node = nullptr;
//			//cout << "!!!!!!!!!!!!!!!!: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//
//			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {// scan sub node in Rtree
//
//				/*................................................*/
//#ifdef _DEBUG
//				//cout << "!!!: " << m_temp_queue_top.p_rtree_node->m_branch[i].m_data << endl;
//				//g_n_account_leaf_node++;
//#endif
//				/*................................................*/
//
//				/*-----------------------------------------get ACPA(origitnal time series) point & id----------------------------------------------------*/
//				f_temp_APLA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);// get APCA(original time series) id
//				//cout << tempAPCAPair.id_originalTimeSeries << endl;
//				//f_temp_APLA_Pair.approximation_pointer = &apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)];// pointer to APCA point
//				/*----------------------------------------------------------------------------------------------------------------------------------------*/
//				//cout << "KNN : data point ID = " << tempAPCAPair.id_originalTimeSeries << endl;
//				//cout << tempAPCAPair.APCAValue << endl;
//
//				switch (input_argument.representation_option) {
//				case 1://MSPLA
//				case 6: //ICDE07
//				case 8: {//Initial 200706
//					//if (multi_y_projection_argument[f_temp_APLA_Pair.original_time_series_id].is_y_projection) {
//						//f_temp_APLA_Pair.d_dist = APLA::get_distance_apca_LB(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_apla_projection_vector));
//					//}
//					//else{
//
//					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       210511    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//					/*-------210511 PAA APCA version has projection of query time series -----------*/
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB(apla_array_vector[f_temp_APLA_Pair.original_time_series_id], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[f_temp_APLA_Pair.original_time_series_id], query_apla_projection_vector));
//					/*------------------------------------------------------------------------------*/
//
//					/*-------210511 PAA APCA version has projection of query time series -----------*/
//					input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);//210603
//					input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//					TOOL::read_normalized_multi_time_series(data_source, f_temp_APLA_Pair.original_time_series_id, original_time_series_vector);
//					TOOL::recordStartTime(TOOL::time_record[14]); //whole time
//					TOOL::recordStartTime(TOOL::time_record[2]);//knn time
//
//					long double pruning_power_distance = 0;//210603
//					f_temp_APLA_Pair.d_dist = APLA::get_distance_SAPLA(query_time_series_vector, original_time_series_vector, apla_array_vector[input_argument.query_time_series_id], apla_array_vector[f_temp_APLA_Pair.original_time_series_id], pruning_power_distance);
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_SAPLA(query_time_series_vector, f_temp_APLA_Pair.original_time_series_id, data_source.read_file_address_vector[0], apla_array_vector[input_argument.query_time_series_id], apla_array_vector[f_temp_APLA_Pair.original_time_series_id], pruning_power_distance);
//					/*===============Evaluation==============*/
//					input_argument.IO_cost += pruning_power_distance;
//					/*=======================================*/
//
//					/*------------------------------------------------------------------------------*/
//					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//
//					//200109 AE
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_AE(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//PLA version
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB_pla(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//double test_reconstruct_dist = f_temp_APLA_Pair.d_dist;
//					//191206 APLA distance LB speed
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB_pla_speed(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//f_temp_APLA_Pair.d_dist = APLA::get_apla_endpoint_distance(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_lower_bound(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
////#ifdef _DEBUG
//					/*==================*/
//					//assert(float(test_reconstruct_dist) == float(f_temp_APLA_Pair.d_dist));
//					//vector<DataType> test_normalized_time_series_vector;
//					//vector<DataType> reconstruct_time_series_vector(query_time_series_vector.size(), INF);
//					//APLA::getAPLAReconstructSeries(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], reconstruct_time_series_vector);
//
//				//	TOOL::read_normalized_multi_time_series(data_source, f_temp_APLA_Pair.original_time_series_id, test_normalized_time_series_vector);
//					//double 
//					//double  test_sum_deviation = PLA_QUAL::getPLASumDeviation(test_normalized_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//double  test_query_deviaiton = PLA_QUAL::getPLASumDeviation(query_time_series_vector, query_linked_list);
//					//cout <<"AAAAAAAAAAAAAA:::" << f_temp_APLA_Pair.d_dist << endl;
//				//	double test_euc = TOOL::distanceEUC(query_time_series_vector, test_normalized_time_series_vector);
//					//assert(f_temp_APLA_Pair.d_dist <= test_euc);
//					//if (f_temp_APLA_Pair.original_time_series_id != 19){
//						//cout << test_sum_deviation <<"    <    "<< f_temp_APLA_Pair.d_dist << " <  "<< test_euc <<endl;
//						//assert(test_sum_deviation <= f_temp_APLA_Pair.d_dist);
//						//if (test_sum_deviation < f_temp_APLA_Pair.d_dist) {
//						//if (test_sum_deviation < f_temp_APLA_Pair.d_dist* 2) {
//							//f_temp_APLA_Pair.d_dist = fabs(f_temp_APLA_Pair.d_dist - test_sum_deviation - test_query_deviaiton);
//						//if(test_sum_deviation  > 10)
//						//f_temp_APLA_Pair.d_dist /= (100);
//						//else {
//							//f_temp_APLA_Pair.d_dist -= test_sum_deviation;
//						//}
//						//}
//
//						/*	for (auto&& au : test_normalized_time_series_vector) {
//								cout << au << ",";
//							}
//							cout << "========================================================================="<<endl;
//							for (auto&& au : reconstruct_time_series_vector) {
//								cout << au << ",";
//							}
//							cout << endl;*/
//							//assert(0);
//						//}
//
//					   // assert(f_temp_APLA_Pair.d_dist <= test_euc);
//						//if (f_temp_APLA_Pair.d_dist > test_euc) {
//						//	cout << test_sum_deviation << "    <    " << f_temp_APLA_Pair.d_dist << " <  " << test_euc << endl;
//						//	data_source.bigger_account++;
//						//		for (auto&& au : test_normalized_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << "\n========================================================================="<<endl;
//						//	for (auto&& au : reconstruct_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;
//
//						//	/*for (auto&& au : query_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;
//						//	cout << "##########################################\n";
//						//	for (auto&& au : test_normalized_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;*/
//						//	assert(0);
//						//}
//
//					//}
//					//test_normalized_time_series_vector.clear();
//					//test_normalized_time_series_vector.shrink_to_fit();
//////#endif
//					//}
//					break;
//				}
//				case 2://PLA
//					f_temp_APLA_Pair.d_dist = PLA_QUAL::getPLADistance(input_argument.point_multi_single_length, pla_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], PLA_query, f_temp_APLA_Pair.d_dist);
//					break;
//				case 3://APCA
//				case 4://PAA
//				case 7: {//PAALM
//					//f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceAE(query_time_series_vector, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceLB(APCA_KNN_QUAL::QAPCAProjection(g_query_time_series, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], QProjection), apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					/*============================================================================================================*/
//					/*vector<double> query_reconstruct_time_series;
//					vector<double> reconstruct_time_series;
//
//					APCA_KNN_QUAL::get_APCA_reconstruction(query_APCA, query_reconstruct_time_series);
//					APCA_KNN_QUAL::get_APCA_reconstruction(apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], reconstruct_time_series);
//					f_temp_APLA_Pair.d_dist = TOOL::distanceEUC(query_reconstruct_time_series, reconstruct_time_series);
//
//					query_reconstruct_time_series.clear();
//					query_reconstruct_time_series.shrink_to_fit();
//					reconstruct_time_series.clear();
//					reconstruct_time_series.shrink_to_fit();*/
//					/*================================================================================================================================*/
//					break;
//				}
//				case 9: {//SAX
//					f_temp_APLA_Pair.d_dist = sax.distance_LB_SAX(query_APCA, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@      " << f_temp_APLA_Pair.d_dist << endl;
//					break;
//				}
//				case 5:
//					assert(0);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//				queue.push(f_temp_APLA_Pair);
//			}
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
//			//	//cout << tempOriginalTimeSeriesPair.key;
//			/*------------------------------------------------------------------------Clear memory-------------------------------------------------------------------------------*/
//			APCA_QUAL::deleteAPCA(QProjection);
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS >> 1 == input_argument.point_dimension && input_argument.point_multi_single_length == input_argument.time_series_length);
//#endif
//
//			//is internal Node
//			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			//APCA_KNN_QUAL::APCA_NODE_PAIR tempApcaPair;
//			MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node> temp_apla_pair;
//			typename APCA_KNN_QUAL::REGION fs_region_G;
//			typename APCA_KNN_QUAL::initialREGION(fs_region_G, input_argument.point_dimension);
//			////cout << "regionNum = " << G.regionNum << endl;
//			//cout << "================: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
//				//cout << "===: " << m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_data << endl;
//				switch (input_argument.representation_option) {
//				case 1://MSPLA
//				case 6://ICDE07
//				case 8://Initial 200706
//					// APCA paper, min id == max id
//					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(reconstruct_query_time_series_vector, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//
//
//					//original min id < max id
//					//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//					//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//					//queue.push(temp_apla_pair);
//					//temp_apla_pair.approximation_pointer = nullptr;
//					break;
//				case 2://PLA
//					/*-----------------*/
//					assert(0);
//					//temp_apla_pair.d_dist = PLA_QUAL::getPLAMBRDistance(input_argument, m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, PLA_query, temp_apla_pair.d_dist);
//					/*-----------------*/
//					break;
//				case 3://APCA
//				case 4://PAA
//				case 7://PAALM
//				case 9://SAX
//					// APCA paper, min id == max id
//					/*-----------------*/
//					assert(0);
//					//temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					/*-----------------*/
//					//original min id < max id
//					//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//					//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//					//queue.push(temp_apla_pair);
//					//tempApcaPair.approximation_pointer = nullptr;
//					break;
//				case 5://CHEBY
//					assert(0);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//				temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//				queue.push(temp_apla_pair);
//			}
//			APCA_KNN_QUAL::deleteREGION(fs_region_G);
//		}
//		else {
//			assert(0);
//		}
//	}
//	/*=====================================================================================================================================================================================================================================================================================================================*/
//
//	/*========================Evaluation==========================================*/
//	input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
//	input_argument.knn_total_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[15]);
//	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//	input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);
//	//assert(input_argument.pruning_power > 0);
//	input_argument.pruning_power = 1;//200929
//	/*============================================================================*/
//
//	/*................................................................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF);
//	cout << "??????????????????     APLA KNN failed    ??????????????????" << endl;
//	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	/*ofstream outfile(input_argument.write_file_name[16] + ".txt", ios::app);
//	assert(outfile.is_open());
//	outfile << "??????????????????     APLA KNN failed    ??????????????????   " << PAA_or_APCA << endl;
//	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
//	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	outfile.close();*/
//	TOOL::printInputArgument(input_argument);
//#endif
//	/*................................................................................................................................*/
//
//	//assert(0);
//	priority_queue<MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node>, vector<MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node>>, MULTI::priorityIncrement<MULTI::RTREE_NODE_PAIR_PARTITION<double, int, U::Node>>>().swap(queue);
//	temp.clear();
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//
//	//f_temp_APLA_Pair.approximation_pointer = nullptr;
//	//f_APLA_Root.approximation_pointer = nullptr;
//
//	for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//		//cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//		result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
//	}
//	assert(result_set.size() <= input_argument.K);
//	result.clear();
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//
//	return result_set;
//}
//
//
//// Improve Index
////************************************
//// Method:all_knn_multi_210512
//// Qualifier: KNN structure for PAA, APCA, APLA, ICDE07
//// Input:
//// Output:
//// Notice: Evaluate SAPLA distance. Add  1 input_argument.IO_cost++; 2 Euclidean distance; in leaf node computation 
//// date: 210512
//// author:
////************************************
//TEMPLATE
//template<typename T, typename Y>
//std::multiset<pair<double, int>> MULTI::all_knn_multi_210512(typename TOOL::DATA_SOURCE& const data_source, typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<Y>& const multi_y_projection_argument, DataType*& g_query_time_series, const RTREE& const apcaRTree, T& const approximation_array) {
//	/*..........................................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.file_id != INF && input_argument.K <= input_argument.point_number && input_argument.time_series_length != INF && input_argument.point_number != INF && input_argument.arity_d != INF && input_argument.K != INF);
//	assert(input_argument.time_series_length == input_argument.point_multi_single_length);
//	assert(input_argument.point_dimension == input_argument.point_multi_single_dimension);
//
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//		cout << "MSPLA(SAPLA KNN: \n";
//		break;
//	case 2://PLA
//		cout << "PLA KNN: \n";
//		break;
//	case 3://APCA
//		cout << "APCA KNN: \n";
//		break;
//	case 4://PAA
//		cout << "PAA KNN: \n";
//		break;
//	case 5://Chebyshev
//		cout << "CHEBY KNN: \n";
//		break;
//	case 6://ICDE07
//		cout << "ICDE07 KNN: \n";
//		break;
//	case 7://PAALM
//		cout << "PAALM KNN: \n";
//		break;
//	case 8://initial 200706
//		cout << "Initial 200706 KNN: \n";
//		break;
//	case 9://SAX
//		cout << "SAX KNN: \n";
//		break;
//	default:
//		assert(0);
//		break;
//	}
//#endif
//	/*..........................................................................................................*/
//
//	/*-------------------------      Evaluation  ---------------------------------*/
//	input_argument.knn_total_time = 0.0;
//	input_argument.IO_cost = 0;// measure I/O cost
//	input_argument.pruning_power = 0;
//	/*-------------------------------------------------------------------------------*/
//
//	/*-------------------------------------change pointer array to vector array--------------------------------*/
//	//int multi_to_single_series_length = input_argument.time_series_length * input_argument.arity_d;
//	int multi_to_single_series_length = input_argument.time_series_length;
//	vector<DataType> query_time_series_vector;
//	query_time_series_vector.resize(multi_to_single_series_length, INF);
//	std::copy_n(g_query_time_series, multi_to_single_series_length, query_time_series_vector.begin());
//
//	vector<DataType> reconstruct_query_time_series_vector;
//	reconstruct_query_time_series_vector.resize(multi_to_single_series_length, INF);
//	/*---------------------------------------------------------------------------------------------------------*/
//
//
//	/*=======================================Evaluation: plot & bar chart=======================================================*/
//#ifdef _DEBUG
//	//g_n_account_apca_point = 0;
//	input_argument.navigate_index_time = 0.0;// navigate time
//	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
//	input_argument.distance_euc_time = 0.0;// distance euclidean time
//#endif
//	/*===========================================================================================================================*/
//
//	/**/
//	typename APCA_QUAL::APCA query_APCA;
//	SaxQuantizer::SAX sax(input_argument.point_dimension);
//	//APCA_QUAL::getAPCAPoint(g_query_time_series, input_argument.time_series_length, input_argument.point_dimension, query_APCA);
//	typename PLA_QUAL::PLA PLA_query;
//	/**/
//
//	/*========================================!!!!!!!!!!!   MSPLA & ICDE07 Variable construct=================================================================*/
//
//	DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>();
//	switch (input_argument.representation_option) {
//	case 1://MSPLA
//	//APLA_ICDE07<DataType>::getAPLA_ICDE07(input_argument, normalized_series_vector, apla_array_vector[point_id]);
//	case 6: //ICDE07
//	case 8: {//Initial 200706
//		//typename TOOL::Y_PROJECTION_ARGUMENT y_projection_argument(input_argument.initial_N);
//		//DoublyLinkedList<APLA::AREA_COEFFICIENT> all_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//		//DoublyLinkedList<APLA::AREA_COEFFICIENT> cluster_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//		query_linked_list.copy(apla_array_vector[input_argument.query_time_series_id]);
//		//APLA::get_APLA_point(input_argument, g_query_time_series, y_projection_argument, all_linked_list, cluster_linked_list, query_linked_list);//191129
//
//		APLA::getAPLAReconstructSeries(query_linked_list, reconstruct_query_time_series_vector);
//
//		/*......................................*/
//#ifdef _DEBUG
//		assert(query_linked_list.size() == apla_array_vector[input_argument.query_time_series_id].size());
//		for (int segment_id = 0; segment_id < query_linked_list.size(); segment_id++) {
//			const auto& const query_segment = query_linked_list[segment_id];
//			const auto& const array_segment = apla_array_vector[input_argument.query_time_series_id][segment_id];
//			//cout <<"!!!!!!!!!!!!" <<query_segment.right_endpoint << ", "<<array_segment.right_endpoint << endl;
//			assert(query_segment.right_endpoint == array_segment.right_endpoint);
//			assert(query_segment.right_endpoint != INF && query_segment.right_endpoint >= 0 && query_segment.right_endpoint < input_argument.point_multi_single_length);
//		}
//#endif
//		/*......................................*/
//
//		break;
//	}
//	case 2: {//PLA
//		PLA_QUAL::initialPLA(PLA_query, input_argument.point_multi_single_dimension);//????????180918 , this multi ot single has problem
//		PLA_QUAL::getPLA(input_argument.point_multi_single_length, input_argument.point_multi_single_dimension, g_query_time_series, PLA_query);
//		break;
//	}
//	case 3://APCA
//	case 4://PAA
//		break;
//	case 5://CHEBY
//		assert(0);
//		break;
//	case 7://PAALM
//		break;
//	case 9: {//SAX
//		APCA_QUAL::initialAPCA(query_APCA, input_argument.point_dimension);
//		sax.get_SAX(query_time_series_vector, input_argument.point_dimension, query_APCA);
//		break;
//	}
//	default:
//		assert(0);
//		break;
//	}
//	/*=======================================================================================================================================================*/
//	/*========================================!!!!!!!!!!!   APLA & ICDE07 Variable construct=================================================================*/
//	int i = NULL, j = NULL;
//	priority_queue<RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>> queue;// <Rtree node, distance>
//	RTREE_NODE_PAIR f_APLA_Root, f_temp_APLA_Pair;// distance, APCA point(original time series) id, APCA point pointer, Rtree sub node pointer
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> temp; // <dist, original time series id>
//	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> result;// <dist, original time series id>
//	std::multiset<pair<double, int>> result_set;//191204 result
//	/*-------------------*/
//	std::multiset<pair<double, int>> temp_list_set; // <dist, original time series id>
//	std::multiset<pair<double, int>> result_list_set;// <dist, original time series id>
//	/*-----------------*/
//	typename APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair; // <dist, original time series id>
//	RTREE_NODE_PAIR m_temp_queue_top;//
//	vector<DataType> original_time_series_vector;
//	/*==========================================================================================================================================*/
//	/*-------------------------Evaluation Time---------------------------------*/
//	TOOL::recordStartTime(TOOL::time_record[14]); //whole time
//	TOOL::recordStartTime(TOOL::time_record[2]);//knn time
//	/*--------------------------------------------------------------------------------*/
//	/*========================================== Variable Initial ===============================================================*/
//	f_APLA_Root.p_rtree_node = apcaRTree.m_root;
//	f_APLA_Root.d_dist = 0;
//	queue.push(f_APLA_Root);
//	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
//	//printf("<///////**    KNN Begin   **////////>\n");
//	/*===========================================================================================================================*/
//
//	/*========================================================================================Begin K-NN=======================================================================================================================================*/
//	while (!queue.empty()) {
//
//		m_temp_queue_top = queue.top();
//
//		//cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator: ";
//		//for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
//		//cout << endl;
//
//		/*===========================================================================Scan temp to get ruslt=========================================================================================================*/
//		for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
//			//cout << " Temp       Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
//			if (plist->d_dist <= m_temp_queue_top.d_dist) {
//				//cout << "           <= " << endl;
//				result.push_back(*plist);
//				plist = temp.erase(plist);
//			}
//			else
//				plist++;
//
//			if (input_argument.K == result.size()) {
//				/*-------------------------Evaluation ---------------------------------*/
//				input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
//				input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//				assert(input_argument.IO_cost != 0 && input_argument.point_number != INF);
//				input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
//				assert(input_argument.pruning_power > 0 && input_argument.pruning_power <= 1);
//				/*---------------------------------------------------------------------*/
//
//				/*..........................................................................................................*/
//#ifdef _DEBUG
//				assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF);
//				/*-------------------------------------------- Print Result ---------------------------------*/
//				cout << "!!!!!!!!!!!!! KNN Find result !!!!!!!!!!!!!!!!!!!!!!!!!!!!   result list size: " << result.size() << endl;
//
//				cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
//
//				cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;
//
//				cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;
//
//				cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;
//
//				cout << "pruning power: " << input_argument.pruning_power << endl;
//
//				cout << "I/O cost: " << input_argument.IO_cost << endl;
//
//				result.sort([](const APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& first, const  APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist; });//small to big
//#endif
//				/*..........................................................................................................*/
//
//				for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//					//cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//					result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
//				}
//				assert(result_set.size() == input_argument.K);
//				//cout << endl;
//				/*--------------------------------------------------------------------------------------------*/
//				/*--------------------------------------     Clear memory    -----------------------------------------------*/
//				//f_temp_APLA_Pair.approximation_pointer = nullptr;
//				//f_APLA_Root.approximation_pointer = nullptr;
//				priority_queue< RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>>().swap(queue);
//
//				temp.clear();
//				result.clear();
//				list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//				//list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//				/*---------------------------------------------------------------------------------------------------------*/
//				return result_set;
//			}
//		}
//		/*==========================================================================================================================================================================*/
//		/*======================================================================Pop top node ine queue==============================================================================*/
//		queue.pop();
//		/*==========================================================================================================================================================================*/
//		/*====================================================================== Check APLA point ==================================================================================*/
//		if (m_temp_queue_top.p_rtree_node == nullptr) { //is Approximation data point
//
//			//cout << "    queue.top is data point\n";
//			//DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
//			//vector<DataType> original_time_series_vector;
//
//			/*------------------------------------------------get original time series from database (file)----------------------------------------------*/
//			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
//			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
//			/*--------------Improve 210512----------------*/
//			//tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;//get top queue time sereis id.
//			/*-------------------------------------------*/
//			//getFileStreamByID(file_name, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
//			// get original time sires from database(txt file)
//
//			/*........................................Improve 210512 ....................................................*/
////#ifdef _DEBUG
////			assert(input_argument.time_series_length == data_source.multi_single_time_series_length && input_argument.read_file_name == data_source.read_file_address);
////#endif
//			/*..........................................................................................................*/
//			/*--------------Improve 210512----------------*/
//			//TOOL::read_normalized_multi_time_series(data_source, m_temp_queue_top.original_time_series_id, original_time_series_vector);
//			/*--------------------------------------------*/
//			//already normalized
//			//			TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series_vector);
//			/*--------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*===============Evaluation==============*/
//			/*--------------Improve 210512----------------*/
//			//input_argument.IO_cost++;
//			/*--------------------------------------------*/
//			/*=======================================*/
//
//			/*------------------------------------------------compute Euclidean distance from query times and candidate original time series-----------------------------------------*/
//			/*--------------Improve 210512----------------*/
//			//tempOriginalTimeSeriesPair.d_dist = TOOL::distanceEUC(query_time_series_vector, original_time_series_vector);
//			/*--------------------------------------------*/
//			//#ifdef _DEBUG
//			//			copy_n(original_time_series_vector.begin(), original_time_series_vector.size(), original_time_series);
//			//			double test_distance = APCA_KNN_QUAL::distanceEUC(g_query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);
//			//			assert(test_distance == tempOriginalTimeSeriesPair.d_dist);
//			//#endif
//						/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//						/*------------------------------------------------temp: insert time seris dist into temp---------------------------------------------------------------------------------*/
//						//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
//			/*--------------Improve 210512----------------*/
//			//temp.push_back(tempOriginalTimeSeriesPair);
//			/*--------------------------------------------*/
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//
//			/*-------------------------------------------------Clear memory----------------------------------------------------------------------------------------------------------*/
//			//original_time_series_vector.clear();
//			//original_time_series_vector.shrink_to_fit();
//			//delete[] original_time_series;
//			//original_time_series = nullptr;
//			/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//		/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
//			/*====================================================================== Check Leaf Node, insert APCA point into queue===============================================================================*/
//			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//
//			/*................................................*/
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS / 2 == input_argument.point_dimension && m_temp_queue_top.p_rtree_node->m_count > 0);
//#endif
//			/*................................................*/
//
//			typename APCA_QUAL::APCA QProjection;
//			APCA_QUAL::initialAPCA(QProjection, input_argument.point_dimension);
//			vector<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_apla_projection_vector(input_argument.point_dimension);
//
//			/*-----------------------------------------------------------------------queue: insert leaf node. compute distance LB--------------------------------------------------------------------------------------------*/
//			f_temp_APLA_Pair.p_rtree_node = nullptr;
//			//cout << "!!!!!!!!!!!!!!!!: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//
//			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {// scan sub node in Rtree
//
//				/*................................................*/
//#ifdef _DEBUG
//				//cout << "!!!: " << m_temp_queue_top.p_rtree_node->m_branch[i].m_data << endl;
//				//g_n_account_leaf_node++;
//#endif
//				/*................................................*/
//
//				/*-----------------------------------------get ACPA(origitnal time series) point & id----------------------------------------------------*/
//				f_temp_APLA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);// get APCA(original time series) id
//				//cout << tempAPCAPair.id_originalTimeSeries << endl;
//				//f_temp_APLA_Pair.approximation_pointer = &apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)];// pointer to APCA point
//				/*----------------------------------------------------------------------------------------------------------------------------------------*/
//				//cout << "KNN : data point ID = " << tempAPCAPair.id_originalTimeSeries << endl;
//				//cout << tempAPCAPair.APCAValue << endl;
//
//				switch (input_argument.representation_option) {
//				case 1://MSPLA
//				case 6: //ICDE07
//				case 8: {//Initial 200706
//					//if (multi_y_projection_argument[f_temp_APLA_Pair.original_time_series_id].is_y_projection) {
//						//f_temp_APLA_Pair.d_dist = APLA::get_distance_apca_LB(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_apla_projection_vector));
//					//}
//					//else{
//
//					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       210511    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//					/*-------210511 PAA APCA version has projection of query time series -----------*/
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB(apla_array_vector[f_temp_APLA_Pair.original_time_series_id], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[f_temp_APLA_Pair.original_time_series_id], query_apla_projection_vector));
//					/*------------------------------------------------------------------------------*/
//
//					/*-------210511 PAA APCA version has projection of query time series -----------*/
//					long double pruning_power_distance = 0;
//					TOOL::read_normalized_multi_time_series(data_source, f_temp_APLA_Pair.original_time_series_id, original_time_series_vector);
//					f_temp_APLA_Pair.d_dist = APLA::get_distance_SAPLA(query_time_series_vector, original_time_series_vector, apla_array_vector[input_argument.query_time_series_id], apla_array_vector[f_temp_APLA_Pair.original_time_series_id], pruning_power_distance);
//
//					/*===============Evaluation==============*/
//					input_argument.IO_cost++;
//					/*=======================================*/
//					tempOriginalTimeSeriesPair.original_time_series_id = f_temp_APLA_Pair.original_time_series_id;//get top queue time sereis id.
//					tempOriginalTimeSeriesPair.d_dist = TOOL::distanceEUC(query_time_series_vector, original_time_series_vector);
//					temp.push_back(tempOriginalTimeSeriesPair);
//					temp_list_set.emplace(make_pair(tempOriginalTimeSeriesPair.d_dist, tempOriginalTimeSeriesPair.original_time_series_id));
//					/*------------------------------------------------------------------------------*/
//					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//
//					//200109 AE
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_AE(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//PLA version
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB_pla(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//double test_reconstruct_dist = f_temp_APLA_Pair.d_dist;
//					//191206 APLA distance LB speed
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_LB_pla_speed(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//f_temp_APLA_Pair.d_dist = APLA::get_apla_endpoint_distance(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
//					//f_temp_APLA_Pair.d_dist = APLA::get_distance_lower_bound(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_linked_list);
////#ifdef _DEBUG
//					/*==================*/
//					//assert(float(test_reconstruct_dist) == float(f_temp_APLA_Pair.d_dist));
//					//vector<DataType> test_normalized_time_series_vector;
//					//vector<DataType> reconstruct_time_series_vector(query_time_series_vector.size(), INF);
//					//APLA::getAPLAReconstructSeries(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], reconstruct_time_series_vector);
//
//				//	TOOL::read_normalized_multi_time_series(data_source, f_temp_APLA_Pair.original_time_series_id, test_normalized_time_series_vector);
//					//double 
//					//double  test_sum_deviation = PLA_QUAL::getPLASumDeviation(test_normalized_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//double  test_query_deviaiton = PLA_QUAL::getPLASumDeviation(query_time_series_vector, query_linked_list);
//					//cout <<"AAAAAAAAAAAAAA:::" << f_temp_APLA_Pair.d_dist << endl;
//				//	double test_euc = TOOL::distanceEUC(query_time_series_vector, test_normalized_time_series_vector);
//					//assert(f_temp_APLA_Pair.d_dist <= test_euc);
//					//if (f_temp_APLA_Pair.original_time_series_id != 19){
//						//cout << test_sum_deviation <<"    <    "<< f_temp_APLA_Pair.d_dist << " <  "<< test_euc <<endl;
//						//assert(test_sum_deviation <= f_temp_APLA_Pair.d_dist);
//						//if (test_sum_deviation < f_temp_APLA_Pair.d_dist) {
//						//if (test_sum_deviation < f_temp_APLA_Pair.d_dist* 2) {
//							//f_temp_APLA_Pair.d_dist = fabs(f_temp_APLA_Pair.d_dist - test_sum_deviation - test_query_deviaiton);
//						//if(test_sum_deviation  > 10)
//						//f_temp_APLA_Pair.d_dist /= (100);
//						//else {
//							//f_temp_APLA_Pair.d_dist -= test_sum_deviation;
//						//}
//						//}
//
//						/*	for (auto&& au : test_normalized_time_series_vector) {
//								cout << au << ",";
//							}
//							cout << "========================================================================="<<endl;
//							for (auto&& au : reconstruct_time_series_vector) {
//								cout << au << ",";
//							}
//							cout << endl;*/
//							//assert(0);
//						//}
//
//					   // assert(f_temp_APLA_Pair.d_dist <= test_euc);
//						//if (f_temp_APLA_Pair.d_dist > test_euc) {
//						//	cout << test_sum_deviation << "    <    " << f_temp_APLA_Pair.d_dist << " <  " << test_euc << endl;
//						//	data_source.bigger_account++;
//						//		for (auto&& au : test_normalized_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << "\n========================================================================="<<endl;
//						//	for (auto&& au : reconstruct_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;
//
//						//	/*for (auto&& au : query_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;
//						//	cout << "##########################################\n";
//						//	for (auto&& au : test_normalized_time_series_vector) {
//						//		cout << au << ",";
//						//	}
//						//	cout << endl;*/
//						//	assert(0);
//						//}
//
//					//}
//					//test_normalized_time_series_vector.clear();
//					//test_normalized_time_series_vector.shrink_to_fit();
//////#endif
//					//}
//					break;
//				}
//				case 2://PLA
//					f_temp_APLA_Pair.d_dist = PLA_QUAL::getPLADistance(input_argument.point_multi_single_length, pla_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], PLA_query, f_temp_APLA_Pair.d_dist);
//					break;
//				case 3://APCA
//				case 4://PAA
//				case 7: {//PAALM
//					//f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceAE(query_time_series_vector, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceLB(APCA_KNN_QUAL::QAPCAProjection(g_query_time_series, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], QProjection), apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					/*============================================================================================================*/
//					/*vector<double> query_reconstruct_time_series;
//					vector<double> reconstruct_time_series;
//
//					APCA_KNN_QUAL::get_APCA_reconstruction(query_APCA, query_reconstruct_time_series);
//					APCA_KNN_QUAL::get_APCA_reconstruction(apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], reconstruct_time_series);
//					f_temp_APLA_Pair.d_dist = TOOL::distanceEUC(query_reconstruct_time_series, reconstruct_time_series);
//
//					query_reconstruct_time_series.clear();
//					query_reconstruct_time_series.shrink_to_fit();
//					reconstruct_time_series.clear();
//					reconstruct_time_series.shrink_to_fit();*/
//					/*================================================================================================================================*/
//					break;
//				}
//				case 9: {//SAX
//					f_temp_APLA_Pair.d_dist = sax.distance_LB_SAX(query_APCA, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
//					//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@      " << f_temp_APLA_Pair.d_dist << endl;
//					break;
//				}
//				case 5:
//					assert(0);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//				queue.push(f_temp_APLA_Pair);
//			}
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
//			//	//cout << tempOriginalTimeSeriesPair.key;
//			/*------------------------------------------------------------------------Clear memory-------------------------------------------------------------------------------*/
//			APCA_QUAL::deleteAPCA(QProjection);
//			/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//			/*============================================================================================================================================================================*/
//		}
//		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {
//#ifdef _DEBUG
//			assert(apcaRTree.NUMDIMS >> 1 == input_argument.point_dimension && input_argument.point_multi_single_length == input_argument.time_series_length);
//#endif
//
//			//is internal Node
//			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			//APCA_KNN_QUAL::APCA_NODE_PAIR tempApcaPair;
//			RTREE_NODE_PAIR temp_apla_pair;
//			typename APCA_KNN_QUAL::REGION fs_region_G;
//			typename APCA_KNN_QUAL::initialREGION(fs_region_G, input_argument.point_dimension);
//			////cout << "regionNum = " << G.regionNum << endl;
//			//cout << "================: " << m_temp_queue_top.p_rtree_node->m_count << endl;
//			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
//				//cout << "===: " << m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_data << endl;
//				switch (input_argument.representation_option) {
//				case 1://MSPLA
//				case 6://ICDE07
//				case 8://Initial 200706
//					// APCA paper, min id == max id
//					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(reconstruct_query_time_series_vector, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//
//
//					//original min id < max id
//					//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//					//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//					//queue.push(temp_apla_pair);
//					//temp_apla_pair.approximation_pointer = nullptr;
//					break;
//				case 2://PLA
//					temp_apla_pair.d_dist = PLA_QUAL::getPLAMBRDistance(input_argument, m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, PLA_query, temp_apla_pair.d_dist);
//					break;
//				case 3://APCA
//				case 4://PAA
//				case 7://PAALM
//				case 9://SAX
//					// APCA paper, min id == max id
//					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//original min id < max id
//					//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
//					//temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//					//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
//					//queue.push(temp_apla_pair);
//					//tempApcaPair.approximation_pointer = nullptr;
//					break;
//				case 5://CHEBY
//					assert(0);
//					break;
//				default:
//					assert(0);
//					break;
//				}
//				temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
//				queue.push(temp_apla_pair);
//			}
//			APCA_KNN_QUAL::deleteREGION(fs_region_G);
//		}
//		else {
//			assert(0);
//		}
//	}
//	/*=====================================================================================================================================================================================================================================================================================================================*/
//
//	/*========================Evaluation==========================================*/
//	input_argument.knn_total_time = TOOL::recordFinishTime(TOOL::time_record[2]);
//	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
//	//assert(input_argument.pruning_power > 0);
//	input_argument.pruning_power = 1;//200929
//	/*============================================================================*/
//
//	/*................................................................................................................................*/
//#ifdef _DEBUG
//	assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF);
//	cout << "??????????????????     APLA KNN failed    ??????????????????" << endl;
//	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	/*ofstream outfile(input_argument.write_file_name[16] + ".txt", ios::app);
//	assert(outfile.is_open());
//	outfile << "??????????????????     APLA KNN failed    ??????????????????   " << PAA_or_APCA << endl;
//	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
//	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
//	outfile.close();*/
//	TOOL::printInputArgument(input_argument);
//#endif
//	/*................................................................................................................................*/
//
//	//assert(0);
//	priority_queue<RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>>().swap(queue);
//	temp.clear();
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
//
//	//f_temp_APLA_Pair.approximation_pointer = nullptr;
//	//f_APLA_Root.approximation_pointer = nullptr;
//
//	for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//		//cout << it->d_dist << ", " << it->original_time_series_id << "; ";
//		result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
//	}
//	assert(result_set.size() <= input_argument.K);
//	result.clear();
//	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
//
//	return result_set;
//}
//
////************************************
//// Method:print_result_each_time_series
//// Qualifier:Print sum deviation, right endpoint, right endpoint of each time series 
//// Input:
//// Output:
//// Notice:
//// date: 200901 
//// author:
////************************************
//TEMPLATE
//template<typename T, typename T1, typename T2, typename T3, typename Y, typename U, typename U1, typename U2>
//void MULTI::print_result_each_time_series(const U& const input_argument, const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const doubly_linked_list, const T1& const pla, const T2& const apca_paa, const T3& const chebyshev, const U1& const output_argument, const U2& const result_record) {
//	
//	TOOL::printInputArgument(input_argument);
//
//	if (input_argument.representation_option == 1 || input_argument.representation_option == 6 || input_argument.representation_option == 8) {//APLA & ICDE07 & initial200706
//
//		switch (input_argument.representation_option) {
//		case 1: {//MSPLA
//			cout << " ====================================    SAPLA   ==========================================\n";
//			break;
//		}
//		case 6: {//ICDE07
//			cout << " ====================================    ICDE07   =========================================\n";
//			break;
//		}
//		case 8: {//initial 200706
//			cout << " ================================    Initial200706   =======================================\n";
//			break;
//		}
//		default:
//			assert(0);
//			break;
//		}
//		
//		TOOL::print_each_segment_right_endpoint(doubly_linked_list);
//	}
//	else if (input_argument.representation_option == 2) {//PLA
//		cout << " ======================================    PLA   ===========================================\n";
//	}
//	else if (input_argument.representation_option == 3 || input_argument.representation_option == 4 || input_argument.representation_option == 7 || input_argument.representation_option == 9) {// PAA & APCA & PAALM & SAX
//		
//		switch (input_argument.representation_option) {
//		case 3: {//APCA
//			cout << " ======================================    APCA   ===========================================\n";
//			break;
//		case 4: //PAA
//			cout << " ======================================    PAA    ===========================================\n";
//			break;
//		case 7:// PAALM 200108
//			cout << " =====================================    PAALM   ===========================================\n";
//			break;
//		case 9:// SAX
//			cout << " =====================================     SAX    ===========================================\n";
//			break;
//		}
//		default:
//			assert(0);
//			break;
//		}
//
//		cout << "each right endpoint: ";
//		for (int segment_id = 0; segment_id < apca_paa.segmentNum; segment_id++) {
//
//#ifdef _DEBUG
//			assert(apca_paa.r[segment_id] != INF && apca_paa.v[segment_id] != INF);
//#endif
//
//			cout << apca_paa.r[segment_id] << ",";
//		}
//		cout << endl;
//
//	}
//	else if (input_argument.representation_option == 5) {//CHEBY
//		cout << " ======================================    Chebyshev   ===========================================\n";
//	}
//	else {
//		assert(0);
//	}
//
//	cout << "Representation Time: " << result_record.representation_time << ", Sum Deviation: "<< result_record.sum_deviation << endl;
//
//}
//
////200930 Optimization KNN. Add APLA & ICDE07
////************************************
//// Method:print_result_each_time_series
//// Qualifier:Print sum deviation, right endpoint, right endpoint of each time series 
//// Input:
//// Output:
//// Notice:
//// date: 200930
//// author:
////************************************
//TEMPLATE
//template<typename T, typename Y, typename U, typename U1, typename U2>
//void MULTI::print_result_each_time_series(const U& const input_argument, const vector<T>& const original_time_series_vector, const Y& const approximation, const U1& const output_argument, const U2& const result_record) {
//	TOOL::printInputArgument(input_argument);
//
//	if (input_argument.representation_option == 1 || input_argument.representation_option == 6 || input_argument.representation_option == 8) {//APLA & ICDE07 & initial200706
//
//		switch (input_argument.representation_option) {
//		case 1: {//MSPLA
//			cout << " ====================================    SAPLA   ==========================================\n";
//			break;
//		}
//		case 6: {//ICDE07
//			cout << " ====================================    ICDE07   =========================================\n";
//			break;
//		}
//		case 8: {//initial 200706
//			cout << " ================================    Initial200706   =======================================\n";
//			break;
//		}
//		default:
//			assert(0);
//			break;
//		}
//
//		TOOL::print_each_segment_right_endpoint(approximation);
//	}
//	else if (input_argument.representation_option == 2) {//PLA
//		cout << " ======================================    PLA   ===========================================\n";
//	}
//	else if (input_argument.representation_option == 3 || input_argument.representation_option == 4 || input_argument.representation_option == 7 || input_argument.representation_option == 9) {// PAA & APCA & PAALM & SAX
//
//		switch (input_argument.representation_option) {
//		case 3: {//APCA
//			cout << " ======================================    APCA   ===========================================\n";
//			break;
//		case 4: //PAA
//			cout << " ======================================    PAA    ===========================================\n";
//			break;
//		case 7:// PAALM 200108
//			cout << " =====================================    PAALM   ===========================================\n";
//			break;
//		case 9:// SAX
//			cout << " =====================================     SAX    ===========================================\n";
//			break;
//		}
//		default:
//			assert(0);
//			break;
//		}
//
//		cout << "each right endpoint: ";
//		for (int segment_id = 0; segment_id < approximation.segmentNum; segment_id++) {
//
//#ifdef _DEBUG
//			assert(approximation.r[segment_id] != INF && approximation.v[segment_id] != INF);
//#endif
//
//			cout << approximation.r[segment_id] << ",";
//		}
//		cout << endl;
//
//	}
//	else if (input_argument.representation_option == 5) {//CHEBY
//		cout << " ======================================    Chebyshev   ===========================================\n";
//	}
//	else {
//		assert(0);
//	}
//
//	cout << "Representation Time: " << result_record.representation_time << ", Sum Deviation: " << result_record.sum_deviation << endl;
//}

#endif
