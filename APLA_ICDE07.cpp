#pragma once
#ifndef _APLA_ICDE07_CPP_
#define _APLA_ICDE07_CPP_

#include "APLA_ICDE07.h"

//************************************
// Method:getAPLA_ICDE07
// Qualifier:Use Linked List ot instead Vector. For every segment, use minmax point as endpoint. Then for every segment use same subsegemnt number to merge.
// Input:
// Output:
// date:191016
// author:
//************************************
TEMPLATE
template<typename T>
APLA_ICDE::APLA_ICDE07(const T& const n, const T& const N) {//
	input_argument.time_series_length = n;
	input_argument.point_dimension = N;
}

//************************************
// Method:getSegmentMaxDifference
// Qualifier:get absolute difference for one segment
// Input:
// Output:
// date:191016
// author:
//************************************
//TEMPLATE
//template<typename T>
//double APLA_ICDE::getSegmentMaxDifference(const vector<DataType>& const original_time_series_vector, const T& const segment_pla) {
//#ifdef _DEBUG
//	//assert(area_vector.size() == input_argument.point_dimension);
//	//assert(segment_pla.rectangle_width == segment_pla.right_endpoint + 1 && *segment_pla.apla.a != INF && *segment_pla.apla.b != INF);
//#endif
//	
//	int variable_id = NULL; //[0-segment_length
//	int array_id = 0;
//	int segment_left_id = segment_pla.right_endpoint - segment_pla.rectangle_width + 1;
//	double max_difference = -1;
//
//#ifdef _DEBUG
//	assert(segment_left_id > -1);
//#endif
//	for (variable_id = 0; variable_id < segment_pla.right_endpoint; variable_id++) {
//		max_difference = max(max_difference, fabs(*segment_pla.apla.a * variable_id + *segment_pla.apla.b - original_time_series_vector[variable_id]));
//	}
//
//	/*=================================Evaluation Result=====================================================================*/
////#ifdef _DEBUG
////	typename TOOL::getDeviation(input_argument, original_time_series, reconstruct_time_series, output_argument);
////	assert(sum_deviation == output_argument.sum_deviation);
////#endif
//	/*..............................................................................................................*/
//	//TOOL::deleteArray(reconstruct_time_series);
//
//	return max_difference;
//}

//**********************************************************************************************************************************
// Method:getAPLA_ICDE07
// Qualifier:Use Linked List ot instead Vector. For every segment, use minmax point as endpoint. Then for every segment use same subsegemnt number to merge.
// Input:
// Output:
// date:191016
// author:
//**********************************************************************************************************************************
TEMPLATE
void APLA_ICDE::getAPLA_ICDE07(typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<DataType>& const original_time_series_vector) {// From paper APLA ICDE 07 
	TOOL::recordStartTime(TOOL::time_record[49]);

	std::vector<std::vector<PLA_QUAL::PLA_C>> L_PLA_factor_vector(input_argument.time_series_length, std::vector<PLA_QUAL::PLA_C>(input_argument.point_dimension + 1));//L[0] is for intitial a&B. L Fit a line PLA L through original time series
	//PLA_by_index_vector.resize();
	std::vector<std::vector<double>> max_error_vector(input_argument.time_series_length, std::vector<double>(input_argument.point_dimension + 1, INF));// the max error of L, 0,1,2,...,right_endpoint. reconstruct_series - origitnal time seires
	std::vector<std::vector<double>> A_best_right_endpoint_vector(input_argument.time_series_length, std::vector<double>(input_argument.point_dimension + 1, 0.0));//A keeps track best result. is the index of the first point of the t_th segment of the best t-segment approximation to	points 0,...,w
	std::vector<PLA_QUAL::PLA_C> F_result_vector(input_argument.point_dimension);//F 

	double c_max_error = -INF;//c
	int w = INF;//w
	int t = INF;// t segment id

	typename PLA_QUAL::PLA_C temp_pla;
	//typename PLA_QUAL::PLA_C last_pla;

	/*=================================================max deviation for every PLA approximation 0,...,w ======================================================*/
	temp_pla.right_endpoint = 0;
	temp_pla.rectangle_width = 1;
	for (w = 1; w < input_argument.time_series_length; w++) {
		temp_pla.right_endpoint++;
		temp_pla.rectangle_width++;
		//PLA_QUAL::getAAndBByPLASegment(original_time_series_vector, temp_pla);
		L_PLA_factor_vector[w][0] = PLA_QUAL::getAAndBByAccumulation(original_time_series_vector, L_PLA_factor_vector[w - 1][0], temp_pla);
		//PLA_by_index_vector.emplace_back(std::vector<PLA_QUAL::PLA>());
		//PLA_by_index_vector[right_endpoint-1].emplace_back(temp_pla);
		//L_PLA_factor_vector[w][0] = temp_pla;
		//max_error_vector.emplace_back(std::vector<double>());
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		//if (w + 1 < input_argument.time_series_length && original_time_series_vector[w + 1] == original_time_series_vector[w]) {
			//max_error_vector[w - 1].emplace_back(INF);
			//max_error_vector[w][0] = -INF;
		//}
		//else {
			//max_error_vector[w - 1].emplace_back(getSegmentMaxDifference(original_time_series_vector, temp_pla));
		max_error_vector[w][0] = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
		//}
		//max_error_vector[w][0] = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//A_right_endpoint_vector.emplace_back(std::vector<double>());
		//A_right_endpoint_vector[w - 1].emplace_back(0);
		A_best_right_endpoint_vector[w][0] = 0;
	}
	
	for (t = 1; t < input_argument.point_dimension; t++) {// segment id
		//cout << "t = " << t <<"--"<< input_argument.point_dimension << endl;
		for (w = t + 1; w < input_argument.time_series_length; w++) {// time series id
			//cout << "w = " << w <<"--"<< input_argument.time_series_length << endl;
			for (int sub_id = t; sub_id < w; sub_id++) {//t to w-1
				/*--------------------------------Maximum max_deviaiton during [sub_id - w], while sub_id++ ---------------------------------------------------*/
				//cout << "sub id = "<< sub_id << " -- "<< w << endl;
				//cout << "----------------------- 1 Get best index-------------------------------\n";
				temp_pla.right_endpoint = w;
				//if (t == 1)
					//temp_pla.rectangle_width = temp_pla.right_endpoint + 1;
				//else
				temp_pla.rectangle_width = temp_pla.right_endpoint - sub_id;
				PLA_QUAL::getAAndBByPLASegment(original_time_series_vector, temp_pla);
				//cout << "--------------------2 compute coefficient--------------------------\n";
				//cout << "L[" << sub_id << "," << temp_pla.right_endpoint << "], segment width: " << temp_pla.rectangle_width << ", a: " << temp_pla.apla.a << ", b: " << temp_pla.apla.b << endl;
				double temp_sub_max_error = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
				c_max_error = max(max_error_vector[sub_id][t - 1], temp_sub_max_error);//maximum max
				//cout << "--------------------3 compare  error--------------------------\n";
				//cout <<"error: [" << sub_id<<","<<t-1<<"]: "<<max_error_vector[sub_id][t - 1] << " VS segment error: " << temp_sub_max_error<< ",   C: " << c_max_error << endl;
				/*---------------------------------------------------------------------------------------------------------------------------------------------*/
				/*if (w < input_argument.time_series_length && w + 1 < input_argument.time_series_length && original_time_series_vector[w + 1] == original_time_series_vector[w]) {
					max_error_vector[w][0] = -INF;
				}
				if (sub_id > 0 && original_time_series_vector[sub_id - 1] == original_time_series_vector[sub_id]) {
					max_error_vector[w][0] = -INF;
				}*/
				//cout << "------------------------4 change best index----------------------\n";
				//cout << "sub id: " << sub_id << " ==  t: " << t << endl;
				//cout << "C:" << c_max_error << " == error[" << w << "," << t << "]: " << max_error_vector[w][t] << endl;
				if (sub_id == t || c_max_error < max_error_vector[w][t]) {

					/*	if (sub_id == t) {
							cout << "sub id: " << sub_id << " ==  t: " << t << endl;
						}
						else{
							cout << " C:"<< c_max_error<<" == error["<<w<<","<<t<<"]?"<< max_error_vector[w][t] <<endl;
						}*/
					L_PLA_factor_vector[w][t] = temp_pla;//L
					max_error_vector[w][t] = c_max_error;// 
					A_best_right_endpoint_vector[w][t] = sub_id;//A

				}
			}
		}
	}
	/*=======================================================================================================*/
	t--;
	c_max_error = max_error_vector[input_argument.time_series_length - 1][t];
	w = input_argument.time_series_length - 1;
	while (t >= 0) {
		//cout<<"L_PLA: " << L_PLA_factor_vector[w][t].apla.a << "   " << L_PLA_factor_vector[w][t].apla.b << endl;
		F_result_vector[t] = L_PLA_factor_vector[w][t];
		w = A_best_right_endpoint_vector[w][t];
		t--;
	}
	output_argument.run_time = TOOL::recordFinishTime(TOOL::time_record[49]);
	//F_result_vector.pop_back();
	/*=======================================================================================================*/
	//F_result_vector[1].rectangle_width = F_result_vector[1].right_endpoint + 1;

	/*=======================================================================================================*/
	
	/*-------------------------------       Print     -------------------------------------------------------------------*/
#ifdef _DEBUG
	/*cout << "L a&b right_endpoint: \n";
	for_each(L_PLA_factor_vector.begin(), L_PLA_factor_vector.end(), [](auto&& const au) {
		for (auto&& auu: au) {
			cout << auu.right_endpoint << " " << auu.rectangle_width << " " << auu.apla.a << " " << auu.apla.b << ",          ";
		}
		cout << endl;
	});

	cout << "Error : \n";
	for_each(max_error_vector.begin(), max_error_vector.end(), [](auto&& const au) {
		for (auto&& auu : au) {
			cout << auu << "         ";
		}
		cout << endl;
	});

	cout << "Best id : \n";
	for_each(A_best_right_endpoint_vector.begin(), A_best_right_endpoint_vector.end(), [](auto&& const au) {
		for (auto&& auu : au) {
			cout << auu << " ";
		}
		cout << endl;
	});*/

	cout << endl;
	cout << "ICDE 07 a&b :  ";
	for_each(F_result_vector.begin(), F_result_vector.end(), [](auto&& const au) {
		cout << au.apla.a << " " << au.apla.b << ", ";
	});
	cout << endl;
#endif

#ifdef _DEBUG
	/*---------------------------------------------Evaluation---------------------------------------------------------------*/
	assert(c_max_error != INF && c_max_error != -INF && F_result_vector.size() == input_argument.point_dimension && F_result_vector.back().right_endpoint == input_argument.time_series_length - 1 && F_result_vector.front().right_endpoint + 1 == F_result_vector.front().rectangle_width);
	double sum_segment_width = F_result_vector.front().rectangle_width;

	for_each(F_result_vector.begin() + 1, F_result_vector.end(), [&](auto&& const au) {
		//cout << std::prev(&au, 1)->right_endpoint << " ";
		assert(au.right_endpoint != INF && au.rectangle_width == au.right_endpoint - std::prev(&au, 1)->right_endpoint);
		sum_segment_width += au.rectangle_width;
	});
	assert(sum_segment_width == input_argument.time_series_length);
	/*-------------------------------------------------------------------------------------------------------------------------*/
#endif
	/*-----------------------------------------Output Result-----------------------------------------------------*/
	cout << "ICDE07 right end point: ";
	for_each(F_result_vector.begin(), F_result_vector.end(), [](auto&& const au) {
		cout << au.right_endpoint + 1 << " ";
	});
	cout << endl;
	cout << "ICDE07 segment width: ";
	for_each(F_result_vector.begin(), F_result_vector.end(), [](auto&& const au) {
		//assert(au.right_endpoint != INF && au.rectangle_width == au.right_endpoint - std::prev(&au, 1)->right_endpoint);
		cout << au.rectangle_width << " ";
	});
	cout << endl;
	/*-------------------------------------------------------------------------------------------------------------*/
	output_argument.sum_deviation = PLA_QUAL::getPLASumDeviation(original_time_series_vector, F_result_vector);// except fisrt element
	/*----------------------------------------------------------------------------------------------------------------------*/
	L_PLA_factor_vector.clear();
	L_PLA_factor_vector.shrink_to_fit();
	max_error_vector.clear();
	max_error_vector.shrink_to_fit();
	A_best_right_endpoint_vector.clear();
	A_best_right_endpoint_vector.shrink_to_fit();
	/*=======================================================================================================*/
}


//**********************************************************************************************************************************
// Method:getAPLA_ICDE07
// Qualifier:Use Linked List ot instead Vector. For every segment, use minmax point as endpoint. Then for every segment use same subsegemnt number to merge.
// Input://PLA_QUAL::PLA_C 191121 has return value
// Output:
// date:191016
// author:
//**********************************************************************************************************************************
TEMPLATE
template<typename T> //PLA_QUAL::PLA_C 191121
void APLA_ICDE::getAPLA_ICDE07(typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<DataType>& const original_time_series_vector, DoublyLinkedList<T>& const F_result_vector) {// From paper APLA ICDE 07 
	TOOL::recordStartTime(TOOL::time_record[49]);

	std::vector<std::vector<T>> L_PLA_factor_vector(input_argument.time_series_length, std::vector<T>(input_argument.point_dimension + 1));//L[0] is for intitial a&B. L Fit a line PLA L through original time series
	/*-------------------------------------------    201021 Speed up ICDE07   ---------------------------------------------------------*/
	vector< vector<T> > store_segment_factor_vector(input_argument.time_series_length, std::vector<T>(input_argument.time_series_length));//201021 store already computed a&b, do not need compute again
	/*---------------------------------------------------------------------------------------------------------------------------------*/
    //PLA_by_index_vector.resize();
	std::vector<std::vector<double>> max_error_vector(input_argument.time_series_length, std::vector<double>(input_argument.point_dimension + 1, INF));// the max error of L, 0,1,2,...,right_endpoint. reconstruct_series - origitnal time seires
	std::vector<std::vector<double>> A_best_right_endpoint_vector(input_argument.time_series_length, std::vector<double>(input_argument.point_dimension + 1, 0));//A keeps track best result. is the index of the first point of the t_th segment of the best t-segment approximation to	points 0,...,w
	//std::vector<PLA_QUAL::PLA_C> F_result_vector(input_argument.point_dimension);//F 
	F_result_vector.resize(input_argument.point_dimension);//F 

	long double temp_sub_max_error;
	double c_max_error = -INF;//c
	int w = INF;//w
	int t = INF;// t segment id

	T temp_pla;
	//typename PLA_QUAL::PLA_C last_pla;

	/*=================================================    max deviation for every PLA approximation 0,...,w     ======================================================*/
	temp_pla.right_endpoint = 0;
	temp_pla.rectangle_width = 1;
	for (w = 1; w < input_argument.time_series_length; w++) {
		temp_pla.right_endpoint++;
		temp_pla.rectangle_width++;
		/*........................get a & b of segment.......................*/
		//PLA_QUAL::getAAndBByPLASegment(original_time_series_vector, temp_pla);

		L_PLA_factor_vector[w][0] = PLA_QUAL::getAAndBByAccumulation(original_time_series_vector, L_PLA_factor_vector[w - 1][0], temp_pla);

		/*-------------------------------------------    201021 Speed up ICDE07   ---------------------------------------------------------*/
		store_segment_factor_vector[w][0] = L_PLA_factor_vector[w][0];//201021 store already computed a&b, do not need compute again
		/*---------------------------------------------------------------------------------------------------------------------------------*/

		//200827
		//L_PLA_factor_vector[w][0] = APLA::get_ab_segment_by_accumulation(original_time_series_vector, L_PLA_factor_vector[w - 1][0], temp_pla);

		/*...................................................................*/
		//L_PLA_factor_vector[w][0] = PLA_QUAL::get_a_b_by_endpoint(original_time_series_vector, temp_pla);//191215
		/*...................................................................*/

		max_error_vector[w][0] = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);

		//A_right_endpoint_vector.emplace_back(std::vector<double>());
		//A_right_endpoint_vector[w - 1].emplace_back(0);
		A_best_right_endpoint_vector[w][0] = 0;
	}
	/*=============================================================================================================================================================*/
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	for (t = 1; t < input_argument.point_dimension; t++) {// segment id
		//cout << "t = " << t <<"--"<< input_argument.point_dimension << endl;
		for (w = t + 1; w < input_argument.time_series_length; w++) {// time series id. w = 2 : n
			//cout << "w = " << w <<"--"<< input_argument.time_series_length << endl;
			for (int sub_id = t; sub_id < w; sub_id++) {//t = 1 to w-1
				/*+++++++++++++++++++++++++++++++++++++++++++++    Maximum max_deviaiton during [sub_id - w], while sub_id++    ++++++++++++++++++++++++++++++++++++++++++*/

				/*================================    201021 get a&b by accumulation   ===========================================*/
				/*================================================================================================================*/

				/*:::::::::::::::::::::::::::::::::::::::   get initial segment  ::::::::::::::::::::::*/
				//cout << "sub id = "<< sub_id << " -- "<< w << endl;
				//cout << "----------------------- 1 Get best index-------------------------------\n";
				temp_pla.right_endpoint = w;
				//if (t == 1)
					//temp_pla.rectangle_width = temp_pla.right_endpoint + 1;
				//else
				temp_pla.rectangle_width = temp_pla.right_endpoint - sub_id;// cannot + 1
				/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

				/*===========================================     Compute a&b, max deviation of one segment. Speed up     =================================================*/
				if (store_segment_factor_vector[w][sub_id].right_endpoint == temp_pla.right_endpoint && store_segment_factor_vector[w][sub_id].rectangle_width == temp_pla.rectangle_width) {
					temp_pla = store_segment_factor_vector[w][sub_id];

					/*:::::::::::::::::::::::      get max deviation of segment      ::::::::::::::::::::::*/
					temp_sub_max_error = store_segment_factor_vector[w][sub_id].area_difference;
					/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
				}
				//accumulate right
				else if (store_segment_factor_vector[w - 1][sub_id].right_endpoint + 1 == temp_pla.right_endpoint && store_segment_factor_vector[w - 1][sub_id].rectangle_width + 1 == temp_pla.rectangle_width) {
#ifdef _DEBUG
					APLA::assert_segment_a_b(original_time_series_vector, store_segment_factor_vector[w - 1][sub_id]);
#endif
					store_segment_factor_vector[w][sub_id] = PLA_QUAL::getAAndBByAccumulation(original_time_series_vector, store_segment_factor_vector[w - 1][sub_id], temp_pla);
#ifdef _DEBUG
					APLA::assert_segment_a_b(original_time_series_vector, temp_pla);
#endif
					/*:::::::::::::::::::::::      get max deviation of segment      ::::::::::::::::::::::*/
					temp_sub_max_error = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
					store_segment_factor_vector[w][sub_id].area_difference = temp_sub_max_error;
					/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
				}
				//accumulate left
				else if (store_segment_factor_vector[w][sub_id + 1].right_endpoint == temp_pla.right_endpoint && store_segment_factor_vector[w][sub_id + 1].rectangle_width + 1 == temp_pla.rectangle_width) {
#ifdef _DEBUG
					APLA::assert_segment_a_b(original_time_series_vector, store_segment_factor_vector[w][sub_id + 1]);
#endif
					store_segment_factor_vector[w][sub_id] = APLA::get_ab_segment_by_accumulation_left(original_time_series_vector[sub_id + 1], store_segment_factor_vector[w][sub_id + 1], temp_pla);
#ifdef _DEBUG
					APLA::assert_segment_a_b(original_time_series_vector, temp_pla);
#endif
					/*:::::::::::::::::::::::      get max deviation of segment      ::::::::::::::::::::::*/
					temp_sub_max_error = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
					store_segment_factor_vector[w][sub_id].area_difference = temp_sub_max_error;
					/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
				}
				//decrement right
				else if (w + 1 < original_time_series_vector.size() && store_segment_factor_vector[w + 1][sub_id].right_endpoint - 1 == temp_pla.right_endpoint && store_segment_factor_vector[w + 1][sub_id].rectangle_width - 1 == temp_pla.rectangle_width) {
#ifdef _DEBUG
					assert(w == temp_pla.right_endpoint && w - temp_pla.rectangle_width == sub_id);
					APLA::assert_segment_a_b(original_time_series_vector, store_segment_factor_vector[w + 1][sub_id]);
#endif
					store_segment_factor_vector[w][sub_id] = APLA::get_ab_segment_by_decrement(original_time_series_vector, store_segment_factor_vector[w + 1][sub_id], temp_pla);
#ifdef _DEBUG
					APLA::assert_segment_a_b(original_time_series_vector, temp_pla);
#endif
					/*:::::::::::::::::::::::      get max deviation of segment      ::::::::::::::::::::::*/
					temp_sub_max_error = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
					store_segment_factor_vector[w][sub_id].area_difference = temp_sub_max_error;
					/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
				}
				//decrement left
				else if (store_segment_factor_vector[w][sub_id - 1].right_endpoint == temp_pla.right_endpoint && store_segment_factor_vector[w][sub_id - 1].rectangle_width - 1 == temp_pla.rectangle_width) {
#ifdef _DEBUG
					APLA::assert_segment_a_b(original_time_series_vector, store_segment_factor_vector[w][sub_id - 1]);
#endif
					store_segment_factor_vector[w][sub_id] = APLA::get_ab_segment_by_decrement_left(original_time_series_vector, store_segment_factor_vector[w][sub_id - 1], temp_pla);
#ifdef _DEBUG
					APLA::assert_segment_a_b(original_time_series_vector, temp_pla);
#endif
					/*:::::::::::::::::::::::      get max deviation of segment      ::::::::::::::::::::::*/
					temp_sub_max_error = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
					store_segment_factor_vector[w][sub_id].area_difference = temp_sub_max_error;
					/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
				}
				else {
					/*:::::::::::::::::::::::::        get a & b of segment       :::::::::::::::::::::::::*/
					PLA_QUAL::getAAndBByPLASegment(original_time_series_vector, temp_pla);
					//PLA_QUAL::get_a_b_by_endpoint(original_time_series_vector, temp_pla);//191215
					store_segment_factor_vector[w][sub_id] = temp_pla;
					/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

					/*:::::::::::::::::::::::      get max deviation of segment      ::::::::::::::::::::::*/
					temp_sub_max_error = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
					store_segment_factor_vector[w][sub_id].area_difference = temp_sub_max_error;
					/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
				}
				/*=======================================================================================================================================================*/

				/*:::::::::::::::::::::::::::::         Max Difference error       :::::::::::::::::::::::::::::*/
				//temp_sub_max_error = PLA_QUAL::getSegmentMaxDifference(original_time_series_vector, temp_pla);
				c_max_error = max(max_error_vector[sub_id][t - 1], temp_sub_max_error);//maximum max
				/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

				/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

				/*+++++++++++++++++++++++++++++++++++++++++++++++++++++       Update max area, L'    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
				/*if (w < input_argument.time_series_length && w + 1 < input_argument.time_series_length && original_time_series_vector[w + 1] == original_time_series_vector[w]) {
					max_error_vector[w][0] = -INF;
				}
				if (sub_id > 0 && original_time_series_vector[sub_id - 1] == original_time_series_vector[sub_id]) {
					max_error_vector[w][0] = -INF;
				}*/
				//cout << "----------------------------------          4 change best index         -----------------------------------------\n";
				//cout << "sub id: " << sub_id << " ==  t: " << t << endl;
				//cout << "C:" << c_max_error << " == error[" << w << "," << t << "]: " << max_error_vector[w][t] << endl;
				if (sub_id == t || c_max_error < max_error_vector[w][t]) {

					/*if (sub_id == t) {
						cout << "sub id: " << sub_id << " ==  t: " << t << endl;
					}
					else{
						cout << " C:"<< c_max_error<<" == error["<<w<<","<<t<<"]?"<< max_error_vector[w][t] <<endl;
					}*/
					L_PLA_factor_vector[w][t] = temp_pla;//L
					max_error_vector[w][t] = c_max_error;// 
					A_best_right_endpoint_vector[w][t] = sub_id;//A

					//cout<< "L["<<w<<","<< t<<"]   =  " << "L'[" << sub_id << "," << temp_pla.right_endpoint << "] segment width: " << temp_pla.rectangle_width << " a: " << temp_pla.apla.apla.a << " b: " << temp_pla.apla.apla.b << endl;
					//cout << "error[" << w << "," << t << "] : " << c_max_error << endl;
					//cout << "A[" << w << "," << t << "] : " << sub_id << endl;
				}
				
				/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
			}
		}
	}
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	t--;

	/*=======================================================================================================*/
	c_max_error = max_error_vector[input_argument.time_series_length - 1][t];
	w = input_argument.time_series_length - 1;
	while (t >= 0) {
		//cout<<"L_PLA: " << L_PLA_factor_vector[w][t].apla.apla.a << "   " << L_PLA_factor_vector[w][t].apla.apla.b << endl;
		F_result_vector[t] = L_PLA_factor_vector[w][t];
		w = A_best_right_endpoint_vector[w][t];
		t--;
	}
	/*=======================================================================================================*/

	output_argument.run_time = TOOL::recordFinishTime(TOOL::time_record[49]);

	/*----------------------------------------------------------------------------------------------------------------------*/
	L_PLA_factor_vector.clear();
	L_PLA_factor_vector.shrink_to_fit();
	max_error_vector.clear();
	max_error_vector.shrink_to_fit();
	A_best_right_endpoint_vector.clear();
	A_best_right_endpoint_vector.shrink_to_fit();
	/*=======================================================================================================*/
}

#endif
