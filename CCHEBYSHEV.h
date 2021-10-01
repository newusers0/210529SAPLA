#pragma once
#ifndef CCHEBYSHEV_H
#define CCHEBYSHEV_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "CPLA.h"
#include "CAPCA_KNN.h"
#include "./lib/chebyshev_polynomial.hpp"

TEMPLATE
class CCHEBYSHEV : public PLA_QUAL, virtual public TOOL {
public:

	//struct APCA_KNN_QUAL::TIME;
	struct APCA_KNN_QUAL::TIME time_record[20];
	struct TOOL::INPUT_ARGUMENT input_argument;
	struct CHEBYSHEV_SHARE;
	CHEBYSHEV_SHARE chebyshev_share;

	struct CHEBYSHEV;
	struct CHEBYSHEV_VECTOR;//multi-dimensional trajectories
	struct ID_DIST;//time series ID & distances
	struct priorityDecrement; //big to small
	struct priorityIncrement; //small to big
	/*struct POINT;
	struct PLA_NODE_PAIR;

	struct ORIGINAL_TIME_SERIES_PAIR;
	struct priorityDistanceEUC;*/

	CCHEBYSHEV() {};
	CCHEBYSHEV(const int& n, const int& N);
	CCHEBYSHEV(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const int& arity_d, const string& read_file_name, string*& const write_file_name);
	~CCHEBYSHEV();

	void initialMultiDimensionTrajectory(const typename TOOL::INPUT_ARGUMENT& input_argument, const DataType**& multi_dimention_trajectories);
	void deleteMultiDimensionTrajectory(const typename TOOL::INPUT_ARGUMENT& input_argument, const DataType**& multi_dimention_trajectories);

	void initialChebyshevVectorArray(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_VECTOR*& const  chebyshev_vector_array);
	void deleteChebyshevVectorArray(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_VECTOR*& const  chebyshev_vector_array);

	void initialCHEBYSHEV_SHARE(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share);
	void initialCHEBYSHEV(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV& const chebyshev);
	void initialCHEBYSHEV(const int& const time_series_length, const int& point_dimension, CHEBYSHEV& const chebyshev);
	void initialCHEBYSHEVArray(const int& array_number, CHEBYSHEV*& const chebyshev_array);
	void initialCHEBYSHEVArray(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV*& const chebyshev_array);
	void initialCHEBYSHEVArray(const int& trajectory_number, const int& const time_series_length, const int& point_dimension, CHEBYSHEV*& const chebyshev_array);

	void deleteCHEBYSHEV_SHARE(CHEBYSHEV_SHARE& const chebyshev_share);
	void deleteCHEBYSHEV(CHEBYSHEV& const chebyshev);
	void deleteCHEBYSHEVArray(const int& array_number, CHEBYSHEV*& const chebyshev_array);

	//normalize ti [-1,1]
	bool normalizeOneOne(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share);
	//interval Ii
	void getInterval(const typename TOOL::INPUT_ARGUMENT& const input_argument, CHEBYSHEV_SHARE& const chebyshev_share);
	//root_tj
	void getRoots(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share);
	// Every root in which Interval; Binary Search
	void getRootsIndex(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share);
	// Every root in which Interval; normal search
	void getRootsIndex1(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share, double*& const index_evaluation);
	// get internal coefficient of f(x), to speed up program
	void getFunctionCoefficient(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share);
	//get chebyshev_share
	void getCHEBYSHEV_SHARE(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share);
	//interval function for chebyshev coefficient, f(t)=g(t)/sqrt(w(t)*|Ii|)
	void getFunction(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& chebyshev_share, CHEBYSHEV& const chebyshev);

	//interval function for approximation, difference is f(t)=g(t)
	void getFunctionApproximation(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& chebyshev_share, CHEBYSHEV& const chebyshev);

	//specific coefficient ci for KNN
	void getSpecificCoefficient(const typename TOOL::INPUT_ARGUMENT& input_argument, const int& coefficient_index, const CHEBYSHEV_SHARE& chebyshev_share, CHEBYSHEV& const chebyshev);
	//All coefficient c1 - cN for KNN
	void getAllCoefficient(const typename TOOL::INPUT_ARGUMENT& input_argument, const CHEBYSHEV_SHARE& chebyshev_share, CHEBYSHEV& const chebyshev);

	//specific coefficient ci for approximation
	void getSpecificCoefficientApproximation(const typename TOOL::INPUT_ARGUMENT& input_argument, const int& coefficient_index, const CHEBYSHEV_SHARE& chebyshev_share, DataType*& const original_time_series, CHEBYSHEV& const chebyshev);
	//All coefficient c1 - cN for approximation
	void getAllCoefficientApproximation(const typename TOOL::INPUT_ARGUMENT& input_argument, const CHEBYSHEV_SHARE& chebyshev_share, DataType*& const original_time_series, CHEBYSHEV& const chebyshev);

	//get chebyshev f(t)=g(t)/sqrt(w(t)*|Ii|)
	void getCHEBYSHEV(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, CHEBYSHEV& const chebyshev);

	//multi project to single
	void getMultiCHEBYSHEV(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, CHEBYSHEV& const chebyshev);

	//get chebyshev approximation. difference is f(t)=g(t)
	void getChebyshevApproximation(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, CHEBYSHEV& const chebyshev);

	// Loop to compute maximum deviation and sum deviation
	void getDeviationIterationChebyshev(typename TOOL::INPUT_ARGUMENT& const input_argument, string*& const write_file_name, DataType*& const original_time_series, CHEBYSHEV_SHARE& const chebyshev_share, const int& degree_m_begin, const int& degree_m_end, const int& interval);
	// Loop to measure pruning power, execution time, I/O cost
	void getPrunIterationChebyshev(typename TOOL::INPUT_ARGUMENT& const input_argument, string*& const write_file_name, DataType*& const original_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, const int& degree_m_begin, const int& degree_m_end, const int& interval);

	//get MBR
	void getChebyshevMBR(const CHEBYSHEV& const chebyshev, CHEBYSHEV& const chebyshev_MBR);

	DataType& getSpecificApproximation(typename TOOL::INPUT_ARGUMENT& input_argument, const CHEBYSHEV& const chebyshev, CHEBYSHEV_SHARE& const chebyshev_share, const int& normalized_array_id, DataType& const specific_approximation);

	void approximateOriginalFunction(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, CHEBYSHEV_SHARE& const chebyshev_share, DataType*& const approximated_time_series);
	//191206 reconstruct time series, Sum deviation
	double get_sum_deviation(typename TOOL::INPUT_ARGUMENT& input_argument, const vector<DataType>& const original_time_series_vector, CHEBYSHEV_SHARE& const chebyshev_share);

	template<typename T, typename Y, typename U>
	long double get_sum_deviation(T& input_argument, const vector<Y>& const original_time_series_vector, CHEBYSHEV_SHARE& const chebyshev_share, U& const result_collection);

	double& getChebyshevDistance(const typename TOOL::INPUT_ARGUMENT& input_argument, const CHEBYSHEV& chebyshev_array_a, const CHEBYSHEV& chebyshev_array_b, double& const chebyshev_distance);

	double& getChebyshevDistance(const int& segment_number, const CHEBYSHEV& chebyshev_array_a, const CHEBYSHEV& chebyshev_array_b, double& const chebyshev_distance);

	RTREE& buidRTreeIndex(typename TOOL::INPUT_ARGUMENT& const  input_argument, RTREE& ChebyshevRTree, CHEBYSHEV*& const chebyshev_array);

	template<typename T>
	void rangeSearch(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, T& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue);

	template<typename T>
	void rangeSearchMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, T& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue);

	template<typename T>
	void rangeSearchMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const typename TOOL::DATA_SOURCE& const data_source, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, T& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue);

	void KNNSearch(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree, priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >& queue_result);

	std::multiset<pair<double, int>> KNNCHEBYMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree);

	//191223 add data_source
	template<typename T>
	std::multiset<pair<double, int>> KNNCHEBYMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const typename TOOL::DATA_SOURCE& const data_source, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, T& const RTree);

	void writeApproximationResult(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const result_time_series, const double& deviation_sum, const double& deviation_max);

	template<typename T>
	void writeSingleResult(const string& write_file_name, T& result);

	void getSeveralRecronstructionError();

	//friend class APCA_KNN_QUAL;
	//friend class CAPCA_KNN<DataType>;

	/*bool compare(ORIGINAL_TIME_SERIES_PAIR& first, ORIGINAL_TIME_SERIES_PAIR& second);

	INPUT_ARGUMENT input_argument;
	CPLA(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const string& read_file_name);*/

	/*void initialPLAArray(const INPUT_ARGUMENT& input_argument, PLA*& const pla_array);
	void deletePLAArrary(const INPUT_ARGUMENT& input_argument, PLA*& const pla_array);

	void getPLA(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, const PLA& pla);
	PLA& getPLAMBR(const PLA& pla, PLA& pla_MBR);
	double& getPLADistance(const INPUT_ARGUMENT& input_argument, const PLA& pla_array, const PLA& pla_array_qeury, double& pla_distance);
	inline double getPointDistanceSquare(const DataType& point_a_x, const DataType& point_a_y, const DataType& point_b_x, const DataType& point_b_y);
	double getNearestDistance(const DataType& point_a_x, const DataType& point_a_y, const DataType& point_b_x, const DataType& point_b_y, const DataType& point_q_x, const DataType& point_q_y);
	double& getPLAMBRSegmentDistance(const RTREE::Rect &MBR, const int& segment_id, const PLA& pla_array_qeury, double& pla_MBR_segment_distance);
	double& getPLAMBRSegmentDistanceBase(const RTREE::Rect &MBR, const int& segment_id, const PLA& pla_array_qeury, double& pla_MBR_segment_distance);
	double& getPLAMBRDistance(const RTREE::Rect &MBR, const PLA& pla_array_qeury, double& pla_MBR_distance);

	RTREE& buidRTreeIndex(const INPUT_ARGUMENT& input_argument, RTREE& PLARTree, PLA* PLA_array_accumulate, const string& file_name);

	bool PLAKNNSearch(const INPUT_ARGUMENT& input_argument, const DataType* g_query_time_series, const RTREE& PLARTree, PLA* PLA_array_accumulate, const int &K, const string& file_name, list<ORIGINAL_TIME_SERIES_PAIR>& result);
	void SimpleBaseKNNSearch(const INPUT_ARGUMENT& input_argument, const DataType *g_query_time_series, const int &K, const string& file_name, priority_queue<ORIGINAL_TIME_SERIES_PAIR, vector<ORIGINAL_TIME_SERIES_PAIR>, priorityDistanceEUC >& q_base_queue);

	void recordStartTime(TIME& time);
	double recordFinishTime(TIME& time, double& whole_first_run_time);
	void writeResult(const INPUT_ARGUMENT& input_argument, const string& write_file_name, list<ORIGINAL_TIME_SERIES_PAIR>& result, priority_queue<ORIGINAL_TIME_SERIES_PAIR, vector<ORIGINAL_TIME_SERIES_PAIR>, priorityDistanceEUC >& q_base_queue);*/

private:
};

//TEMPLATE
//struct CHEBYSHEV_QUAL::CHEBYSHEV_SHARE {
//	int segmentNum = NULL; //the dimension of the index(eg. MBR, imput parameter)
//
//	double* t = nullptr;// -1 <= ti <= 1
//	double* interval = nullptr;// Ii interval
//	double* root_tj = nullptr;//  roots of PN(t)
//	int* root_tj_index = nullptr;// every root belongs to which interval
//	double* function_coefficient = nullptr; // internal coefficient of f(x), sqrt(w(t)*|Ii|)
//	//double* f = nullptr;// interval function f(t)
//	double* polynomials = nullptr;// Chebyshev polynomials Pm(t)
//	//double* coefficient = nullptr;// coefficient ci
//
//	~CHEBYSHEV_SHARE() {
//		segmentNum = NULL;
//		if (t != nullptr) {
//			delete[] t;
//			t = nullptr;
//		}
//
//		if (interval != nullptr) {
//			delete[] interval;
//			interval = nullptr;
//		}
//
//		if (root_tj != nullptr) {
//			delete[] root_tj;
//			root_tj = nullptr;
//		}
//
//		if (root_tj_index != nullptr) {
//			delete[] root_tj_index;
//			root_tj_index = nullptr;
//		}
//
//		if (function_coefficient != nullptr) {
//			delete[] function_coefficient;
//			function_coefficient = nullptr;
//		}
//
//		/*if (f != nullptr) {
//			delete[] f;
//			f = nullptr;
//		}*/
//
//		if (polynomials != nullptr) {
//			delete[] polynomials;
//			polynomials = nullptr;
//		}
//
//		/*if (coefficient != nullptr) {
//			delete[] coefficient;
//			coefficient = nullptr;
//		}*/
//	}
//};

//TEMPLATE
//struct CHEBYSHEV_QUAL::CHEBYSHEV {
//	int segmentNum = NULL; //the dimension of the index(eg. MBR, imput parameter)
//
//	//double* t = nullptr;// -1 <= ti <= 1
//	//double* interval = nullptr;// Ii interval
//	//double* root_tj = nullptr;//  roots of PN(t)
//	//double* root_tj_index = nullptr;// every root belongs to which interval
//	//double* function_coefficient = nullptr; // internal coefficient of f(x), sqrt(w(t)*|Ii|)
//	double* f = nullptr;// interval function f(t)
//	//double* polynomials = nullptr;// Chebyshev polynomials Pm(t)
//	double* coefficient = nullptr;// coefficient ci length is N+1
//
//	~CHEBYSHEV() {
//		segmentNum = NULL;
//		/*if (t != nullptr) {
//			delete[] t;
//			t = nullptr;
//		}
//
//		if (interval != nullptr) {
//			delete[] interval;
//			interval = nullptr;
//		}
//
//		if (root_tj != nullptr) {
//			delete[] root_tj;
//			root_tj = nullptr;
//		}
//
//		if (root_tj_index != nullptr) {
//			delete[] root_tj_index;
//			root_tj_index = nullptr;
//		}
//
//		if (function_coefficient != nullptr) {
//			delete[] function_coefficient;
//			function_coefficient = nullptr;
//		}*/
//
//		if (f != nullptr) {
//			delete[] f;
//			f = nullptr;
//		}
//
//		/*if (polynomials != nullptr) {
//			delete[] polynomials;
//			polynomials = nullptr;
//		}*/
//
//		if (coefficient != nullptr) {
//			delete[] coefficient;
//			coefficient = nullptr;
//		}
//	}
//};

//TEMPLATE
//CHEBYSHEV_QUAL::CCHEBYSHEV(const int& n, const int& N) {
//	input_argument.time_series_length = n;
//	input_argument.point_dimension = N;//number of coefficient, for APCA & PLA n=m+1
//	input_argument.degree_m = N - 1;//degree of polymonail n=m+1
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
//
//	initialCHEBYSHEV_SHARE(input_argument, chebyshev_share);
//	getCHEBYSHEV_SHARE(input_argument, chebyshev_share);
//}

//TEMPLATE
//CHEBYSHEV_QUAL::~CCHEBYSHEV() {
//	deleteCHEBYSHEV_SHARE(chebyshev_share);
//}

//TEMPLATE
//void CHEBYSHEV_QUAL::initialCHEBYSHEV_SHARE(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& chebyshev_share) {
//	chebyshev_share.segmentNum = input_argument.point_dimension;
//
//	chebyshev_share.t = new double[input_argument.time_series_length];
//	fill_n(chebyshev_share.t, input_argument.time_series_length, NULL);
//
//	chebyshev_share.interval = new double[input_argument.time_series_length << 1];
//	fill_n(chebyshev_share.interval, input_argument.time_series_length << 1, NULL);
//
//	chebyshev_share.root_tj = new double[input_argument.time_series_length];
//	fill_n(chebyshev_share.root_tj, input_argument.time_series_length, NULL);
//
//	chebyshev_share.root_tj_index = new int[input_argument.time_series_length];
//	fill_n(chebyshev_share.root_tj_index, input_argument.time_series_length, NULL);
//
//	chebyshev_share.function_coefficient = new double[input_argument.time_series_length];
//	fill_n(chebyshev_share.function_coefficient, input_argument.time_series_length, NULL);
//
//	/*chebyshev_share.f = new double[input_argument.time_series_length];
//	fill_n(chebyshev_share.f, input_argument.time_series_length, NULL);*/
//
//	chebyshev_share.polynomials = new double[input_argument.time_series_length];
//	fill_n(chebyshev_share.polynomials, input_argument.time_series_length, NULL);
//
//	/*chebyshev_share.coefficient = new double[input_argument.time_series_length];
//	fill_n(chebyshev_share.coefficient, input_argument.time_series_length, NULL);*/
//}

//TEMPLATE
//void CHEBYSHEV_QUAL::getCHEBYSHEV_SHARE(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
//	assert(input_argument.point_dimension == chebyshev_share.segmentNum);
//	assert(input_argument.time_series_length > 1);
//	assert(chebyshev_share.t != nullptr);
//
//	//t(i) [-1,1]
//	normalizeOneOne(input_argument, chebyshev_share);
//
//	//I(i)
//	getInterval(input_argument, chebyshev_share);
//
//	//root_tj
//	getRoots(input_argument, chebyshev_share);
//
//	//root_index
//	getRootsIndex(input_argument, chebyshev_share);
//
//	////root_index1
//	//getRootsIndex1(input_argument, chebyshev_share);
//
//	//f(t) coefficient
//	getFunctionCoefficient(input_argument, chebyshev_share);
//}
//
//TEMPLATE//[-1,1]
//bool CHEBYSHEV_QUAL::normalizeOneOne(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
//	/*double* max_value = max_element(original_time_series, original_time_series + input_argument.time_series_length);
//	double* min_value = min_element(original_time_series, original_time_series + input_argument.time_series_length);*/
//
//	double max_value = input_argument.time_series_length - 1;
//	double min_value = 0;
//
//	/*double max_value = input_argument.time_series_length;
//	double min_value = 1;*/
//
//	double max_min_coefficient = 2.0 / (max_value - min_value);
//
//	for (int t = 0; t < input_argument.time_series_length; t++) {
//		chebyshev_share.t[t] = max_min_coefficient * (t - min_value) - 1.0;
//	}
//
//	/*cout << "ti: "<<endl;
//	APCA_KNN_QUAL::printArray(chebyshev_share.t, input_argument.time_series_length);*/
//
//	return true;
//}

//TEMPLATE//Ii
//void CHEBYSHEV_QUAL::getInterval(const typename TOOL::INPUT_ARGUMENT& const input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
//	assert(chebyshev_share.t[input_argument.time_series_length - 1] != NULL);
//
//	//double* interval_test = new double[input_argument.time_series_length];
//
//	chebyshev_share.interval[0] = -1.0;
//	chebyshev_share.interval[1] = (chebyshev_share.t[0] + chebyshev_share.t[1]) / 2.0;
//
//	/*interval_test[0] = -1;
//
//	for (int inerval_id = 0; inerval_id < input_argument.time_series_length; inerval_id++) {
//		interval_test[inerval_id]= (chebyshev_share.t[inerval_id] + chebyshev_share.t[i]) / 2.0
//	}*/
//
//	for (int i = 1; i <= input_argument.time_series_length - 2; i++) {
//		chebyshev_share.interval[i * 2] = (chebyshev_share.t[i - 1] + chebyshev_share.t[i]) / 2.0;
//		chebyshev_share.interval[i * 2 + 1] = (chebyshev_share.t[i] + chebyshev_share.t[i + 1]) / 2.0;
//	}
//
//	chebyshev_share.interval[input_argument.time_series_length * 2 - 2] = (chebyshev_share.t[input_argument.time_series_length - 2] + chebyshev_share.t[input_argument.time_series_length - 1]) / 2.0;
//	chebyshev_share.interval[input_argument.time_series_length * 2 - 1] = 1.0;
//
//	/*cout << "Ii: \n";
//	APCA_KNN_QUAL::printArray(chebyshev_share.interval, input_argument.time_series_length<<1);*/
//	//delete[] interval_test;
//}
//
//TEMPLATE//root tj
//void CHEBYSHEV_QUAL::getRoots(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
//	double temp_coefficient = boost::math::constants::pi<double>() / double(input_argument.time_series_length);
//
//	for (int root_index = 1; root_index <= input_argument.time_series_length; root_index++) {
//		chebyshev_share.root_tj[root_index - 1] = cos((double(root_index) - 0.5) * temp_coefficient);
//	}
//
//	/*double* root_test = new double[input_argument.time_series_length];
//
//	for (int root_index = 0; root_index < input_argument.time_series_length; root_index++) {
//		root_test[root_index] = cos((double(root_index) - 0.5)*temp_coefficient);
//	}*/
//
//	/*cout << "root_test: \n";
//	APCA_KNN_QUAL::printArray(root_test, input_argument.time_series_length);*/
//
//	/*cout<<"root_tj: \n";
//	APCA_KNN_QUAL::printArray(chebyshev_share.root_tj, input_argument.time_series_length);*/
//
//	/*delete[] root_test;*/
//}
//
//TEMPLATE//every root in which interval; Binary search
//void CHEBYSHEV_QUAL::getRootsIndex(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
//	assert(chebyshev_share.root_tj[input_argument.time_series_length - 1] != NULL);
//	assert(chebyshev_share.root_tj_index != nullptr);
//
//	double* index_evaluation = new double[input_argument.time_series_length];
//	getRootsIndex1(input_argument, chebyshev_share, index_evaluation);
//
//	int left_index = NULL;
//	int right_index = NULL;
//	int middle_index = NULL;
//
//	for (int i = 0; i < input_argument.time_series_length; i++) {
//		if (chebyshev_share.root_tj[i] == 1) {
//			chebyshev_share.root_tj_index[i] = input_argument.time_series_length - 1;
//			continue;
//		}
//
//		left_index = 0;
//		right_index = input_argument.time_series_length - 1;
//		middle_index = NULL;
//		while (left_index <= right_index) {
//			middle_index = left_index + (right_index - left_index) / 2;
//			assert(middle_index == (left_index + right_index) >> 1);//test
//
//			if (chebyshev_share.root_tj[i] >= chebyshev_share.interval[middle_index << 1] && chebyshev_share.root_tj[i] < chebyshev_share.interval[middle_index * 2 + 1]) {
//				//assert(middle_index >> 1 == middle_index / 2);
//				chebyshev_share.root_tj_index[i] = middle_index;
//				assert(middle_index < input_argument.time_series_length&& middle_index >= 0);
//				break;
//			}
//			else if (chebyshev_share.root_tj[i] < chebyshev_share.interval[middle_index << 1]) {
//				right_index = middle_index - 1;
//			}
//			else if (chebyshev_share.root_tj[i] >= chebyshev_share.interval[middle_index * 2 + 1]) {
//				left_index = middle_index + 1;
//			}
//			else {
//				assert(0);
//			}
//		}
//	}
//
//	//assert
//	for (int interval_id = 0; interval_id < input_argument.time_series_length; interval_id++) {
//		assert(chebyshev_share.root_tj_index[interval_id] == index_evaluation[interval_id]);
//	}
//
//	/*cout << "root_index: \n";
//	APCA_KNN_QUAL::printArray(chebyshev_share.root_tj_index, input_argument.time_series_length);*/
//
//	delete[] index_evaluation;
//	index_evaluation = nullptr;
//}
//
//TEMPLATE//every root in which interval normal search
//void CHEBYSHEV_QUAL::getRootsIndex1(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share, double*& const index_evaluation) {
//	assert(chebyshev_share.root_tj[input_argument.time_series_length - 1] != NULL);
//	assert(chebyshev_share.root_tj_index != nullptr);
//
//	int left_index = NULL;
//	int right_index = NULL;
//	int middle_index = NULL;
//
//	for (int root_id = 0; root_id < input_argument.time_series_length; root_id++) {
//		for (int interval_id = 0; interval_id < input_argument.time_series_length; interval_id++) {
//			if (chebyshev_share.root_tj[root_id] >= chebyshev_share.interval[interval_id * 2] && chebyshev_share.root_tj[root_id] < chebyshev_share.interval[interval_id * 2 + 1]) {
//				index_evaluation[root_id] = interval_id;
//			}
//		}
//	}
//
//	/*cout << "root_index1: \n";
//	APCA_KNN_QUAL::printArray(chebyshev_share.root_tj_index, input_argument.time_series_length);*/
//}
//
//TEMPLATE// get internal coefficient of f(x), to speed up program
//void CHEBYSHEV_QUAL::getFunctionCoefficient(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
//	assert(chebyshev_share.t[input_argument.time_series_length - 1] != NULL);
//	assert(chebyshev_share.interval[input_argument.time_series_length * 2 - 1] != NULL);
//	assert(chebyshev_share.root_tj[input_argument.time_series_length - 1] != NULL);
//	double w_t = NULL;//w(t) Chebyshev weight function
//	double I_length = NULL;//|Ii|
//	int subinterval_index = NULL;
//
//	for (int i = 0; i < input_argument.time_series_length; i++) {
//		w_t = 1.0 / sqrt(1.0 - chebyshev_share.root_tj[i] * chebyshev_share.root_tj[i]);
//		//cout << "w(t): " << w_t << endl;
//		//cout << chebyshev_share.interval[chebyshev_share.root_tj_index[i] * 2 + 1] << " "<< chebyshev_share.interval[chebyshev_share.root_tj_index[i] << 1] <<endl;
//		I_length = chebyshev_share.interval[chebyshev_share.root_tj_index[i] * 2 + 1] - chebyshev_share.interval[chebyshev_share.root_tj_index[i] << 1];
//		//cout << "|Ii|: " << I_length << endl;
//		chebyshev_share.function_coefficient[i] = sqrt(w_t * I_length);
//	}
//
//	/*cout << "f(tj)_coefficient:\n";
//	APCA_KNN_QUAL::printArray(chebyshev_share.function_coefficient, input_argument.time_series_length);*/
//}

//TEMPLATE
//void CHEBYSHEV_QUAL::deleteCHEBYSHEV_SHARE(CHEBYSHEV_SHARE& chebyshev_share) {
//	chebyshev_share.segmentNum = NULL;
//
//	delete[] chebyshev_share.t;
//	chebyshev_share.t = nullptr;
//
//	delete[] chebyshev_share.interval;
//	chebyshev_share.interval = nullptr;
//
//	if (chebyshev_share.root_tj != nullptr) {
//		delete[] chebyshev_share.root_tj;
//		chebyshev_share.root_tj = nullptr;
//	}
//
//	if (chebyshev_share.root_tj_index != nullptr) {
//		delete[] chebyshev_share.root_tj_index;
//		chebyshev_share.root_tj_index = nullptr;
//	}
//
//	if (chebyshev_share.function_coefficient != nullptr) {
//		delete[] chebyshev_share.function_coefficient;
//		chebyshev_share.function_coefficient = nullptr;
//	}
//
//	/*delete[] chebyshev_share.f;
//	chebyshev_share.f = nullptr;*/
//
//	delete[] chebyshev_share.polynomials;
//	chebyshev_share.polynomials = nullptr;
//
//	/*if (chebyshev_share.coefficient != nullptr) {
//		delete[] chebyshev_share.coefficient;
//		chebyshev_share.coefficient = nullptr;
//	}*/
//}


//#include "CCHEBYSHEV.cpp"

//template class CCHEBYSHEV<DataType>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPLATE
struct CHEBYSHEV_QUAL::CHEBYSHEV_SHARE {
	int segmentNum = NULL; //the dimension of the index(eg. MBR, imput parameter)

	double* t = nullptr;// -1 <= ti <= 1
	double* interval = nullptr;// Ii interval
	double* root_tj = nullptr;//  roots of PN(t)
	int* root_tj_index = nullptr;// every root belongs to which interval
	double* function_coefficient = nullptr; // internal coefficient of f(x), sqrt(w(t)*|Ii|)
	//double* f = nullptr;// interval function f(t)
	double* polynomials = nullptr;// Chebyshev polynomials Pm(t)
	//double* coefficient = nullptr;// coefficient ci

	~CHEBYSHEV_SHARE() {
		segmentNum = NULL;
		if (t != nullptr) {
			delete[] t;
			t = nullptr;
		}

		if (interval != nullptr) {
			delete[] interval;
			interval = nullptr;
		}

		if (root_tj != nullptr) {
			delete[] root_tj;
			root_tj = nullptr;
		}

		if (root_tj_index != nullptr) {
			delete[] root_tj_index;
			root_tj_index = nullptr;
		}

		if (function_coefficient != nullptr) {
			delete[] function_coefficient;
			function_coefficient = nullptr;
		}

		/*if (f != nullptr) {
			delete[] f;
			f = nullptr;
		}*/

		if (polynomials != nullptr) {
			delete[] polynomials;
			polynomials = nullptr;
		}

		/*if (coefficient != nullptr) {
			delete[] coefficient;
			coefficient = nullptr;
		}*/
	}
};

TEMPLATE
struct CHEBYSHEV_QUAL::CHEBYSHEV {
	int segmentNum = NULL; //the dimension of the index(eg. MBR, imput parameter)

	//double* t = nullptr;// -1 <= ti <= 1
	//double* interval = nullptr;// Ii interval
	//double* root_tj = nullptr;//  roots of PN(t)
	//double* root_tj_index = nullptr;// every root belongs to which interval
	//double* function_coefficient = nullptr; // internal coefficient of f(x), sqrt(w(t)*|Ii|)
	double* f = nullptr;// interval function f(t)
	//double* polynomials = nullptr;// Chebyshev polynomials Pm(t)
	double* coefficient = nullptr;// coefficient ci length is N+1

	~CHEBYSHEV() {
		segmentNum = NULL;
		/*if (t != nullptr) {
			delete[] t;
			t = nullptr;
		}

		if (interval != nullptr) {
			delete[] interval;
			interval = nullptr;
		}

		if (root_tj != nullptr) {
			delete[] root_tj;
			root_tj = nullptr;
		}

		if (root_tj_index != nullptr) {
			delete[] root_tj_index;
			root_tj_index = nullptr;
		}

		if (function_coefficient != nullptr) {
			delete[] function_coefficient;
			function_coefficient = nullptr;
		}*/

		if (f != nullptr) {
			delete[] f;
			f = nullptr;
		}

		/*if (polynomials != nullptr) {
			delete[] polynomials;
			polynomials = nullptr;
		}*/

		if (coefficient != nullptr) {
			delete[] coefficient;
			coefficient = nullptr;
		}
	}
};

TEMPLATE
struct CHEBYSHEV_QUAL::CHEBYSHEV_VECTOR {
	int arity_d = NULL; // airty_d of coefficient for every Chebyshev.

	CHEBYSHEV_QUAL::CHEBYSHEV* multi_dimension_chebyshev = nullptr;//d dimension

	~CHEBYSHEV_VECTOR() {
		assert(arity_d != NULL);
		if (multi_dimension_chebyshev != nullptr) {
			for (int arity_id = 0; arity_id < arity_d; arity_id++) {
				if (multi_dimension_chebyshev[arity_id].f != nullptr) {
					delete[] multi_dimension_chebyshev[arity_id].f;
					multi_dimension_chebyshev[arity_id].f = nullptr;
				}

				if (multi_dimension_chebyshev[arity_id].coefficient != nullptr) {
					delete[] multi_dimension_chebyshev[arity_id].coefficient;
					multi_dimension_chebyshev[arity_id].coefficient = nullptr;
				}
			}

			delete[] multi_dimension_chebyshev;
			multi_dimension_chebyshev = nullptr;
		}

		arity_d = NULL;
	}
};

TEMPLATE
struct CHEBYSHEV_QUAL::ID_DIST {
	double d_dist = NULL;					//dist.
	int trajectory_id = NULL;      // APCA point ID.
};

TEMPLATE
struct CHEBYSHEV_QUAL::priorityDecrement {//big to small
	bool operator ()(const ID_DIST& a, const ID_DIST& b) {
		return a.d_dist < b.d_dist;
	}
};

TEMPLATE
struct CHEBYSHEV_QUAL::priorityIncrement {//small to big
	virtual bool operator ()(const ID_DIST& a, const ID_DIST& b) {
		return a.d_dist > b.d_dist;
	}
};

TEMPLATE
CHEBYSHEV_QUAL::CCHEBYSHEV(const int& n, const int& N) {
	input_argument.time_series_length = n;
	input_argument.point_dimension = N;//number of coefficient, for APCA & PLA n=m+1
	input_argument.degree_m = N - 1;//degree of polymonail n=m+1
	input_argument.pruning_power = 0.0;
	input_argument.sum_distance_euc = 0.0;
	//time
	input_argument.build_rtree_time = 0.0;
	input_argument.approximation_query_time = 0.0;
	input_argument.knn_rest_part_time = 0.0;
	input_argument.knn_total_time = 0.0;

	input_argument.navigate_index_time = 0;// navigate time
	input_argument.distance_lowbound_time = 0; // distance chebyshev, PLA, APCA time
	input_argument.distance_euc_time = 0;// distance euclidean time

	input_argument.IO_cost = 0.0;
	input_argument.result_accuracy = 0;

	initialCHEBYSHEV_SHARE(input_argument, chebyshev_share);
	getCHEBYSHEV_SHARE(input_argument, chebyshev_share);
}

TEMPLATE
CHEBYSHEV_QUAL::CCHEBYSHEV(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const int& arity_d, const string& read_file_name, string*& const write_file_name) {
	assert(arity_d > 0 && N > 0);
	assert(point_number >= K);

	input_argument.time_series_length = n;
	input_argument.point_dimension = N;//number of coefficient, for APCA & PLA n=m+1
	input_argument.degree_m = N - 1;//degree of polymonail n=m+1
	input_argument.point_number = point_number;
	input_argument.rtree_max_nodes = rtree_max_nodes;
	input_argument.K = K;
	input_argument.arity_d = arity_d;
	input_argument.read_file_name = read_file_name;
	input_argument.write_file_name = write_file_name;

	input_argument.pruning_power = 0.0;
	input_argument.sum_distance_euc = 0.0;
	//time
	input_argument.build_rtree_time = 0.0;
	input_argument.approximation_query_time = 0.0;
	input_argument.knn_rest_part_time = 0.0;
	input_argument.knn_total_time = 0.0;

	input_argument.navigate_index_time = 0;// navigate time
	input_argument.distance_lowbound_time = 0; // distance chebyshev, PLA, APCA time
	input_argument.distance_euc_time = 0;// distance euclidean time

	input_argument.IO_cost = 0.0;
	input_argument.result_accuracy = 0;

	initialCHEBYSHEV_SHARE(input_argument, chebyshev_share);
	getCHEBYSHEV_SHARE(input_argument, chebyshev_share);
}

TEMPLATE
CHEBYSHEV_QUAL::~CCHEBYSHEV() {
	deleteCHEBYSHEV_SHARE(chebyshev_share);
}

TEMPLATE
void CHEBYSHEV_QUAL::initialMultiDimensionTrajectory(const typename TOOL::INPUT_ARGUMENT& input_argument, const DataType**& multi_dimention_trajectories) {
	assert(0);
	/*multi_dimention_trajectories = new DataType * [input_argument.airty_d];

	for (int array_id = 0; array_id < input_argument.airty_d; array_id++) {
		multi_dimention_trajectories[array_id] = new DataType[input_argument.time_series_length];
		fill_n(multi_dimention_trajectories[array_id], input_argument.time_series_length, NULL);
	}*/
}

TEMPLATE
void CHEBYSHEV_QUAL::deleteMultiDimensionTrajectory(const typename TOOL::INPUT_ARGUMENT& input_argument, const DataType**& multi_dimention_trajectories) {
	for (int array_id = 0; array_id < input_argument.arity_d; array_id++) {
		if (multi_dimention_trajectories[array_id] != nullptr) {
			delete[] multi_dimention_trajectories[array_id];
			multi_dimention_trajectories[array_id] = nullptr;
		}
	}
	delete[] multi_dimention_trajectories;
	multi_dimention_trajectories = nullptr;
}

TEMPLATE
void CHEBYSHEV_QUAL::initialChebyshevVectorArray(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_VECTOR*& const  chebyshev_vector_array) {
	for (int data_point_id = 0; data_point_id < input_argument.point_number; data_point_id++) {
		chebyshev_vector_array[data_point_id].arity_d = input_argument.arity_d;
		initialCHEBYSHEVArray(input_argument.arity_d, chebyshev_vector_array[data_point_id].multi_dimension_chebyshev);

		/*for (int array_id = 0; array_id < input_argument.arity_d; array_id++) {
			chebyshev_vector_array[data_point_id].multi_dimension_chebyshev[array_id].f = new DataType[input_argument.time_series_length];
			chebyshev_vector_array[data_point_id].multi_dimension_chebyshev[array_id].coefficient = new DataType[input_argument.degree_m+1];
		}*/
	}
}

TEMPLATE
void CHEBYSHEV_QUAL::deleteChebyshevVectorArray(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_VECTOR*& const chebyshev_vector_array) {
	for (int data_point_id = 0; data_point_id < input_argument.point_number; data_point_id++) {
		for (int arity_id = 0; arity_id < input_argument.arity_d; arity_id++) {
			/*delete[] chebyshev_vector_array[data_point_id].multi_dimension_chebyshev[arity_id].f ;
			delete[] chebyshev_vector_array[data_point_id].multi_dimension_chebyshev[arity_id].coefficient;*/

			deleteCHEBYSHEV(chebyshev_vector_array[data_point_id].multi_dimension_chebyshev[arity_id]);
		}

		if (chebyshev_vector_array[data_point_id].multi_dimension_chebyshev != nullptr) {
			delete[] chebyshev_vector_array[data_point_id].multi_dimension_chebyshev;
			chebyshev_vector_array[data_point_id].multi_dimension_chebyshev = nullptr;
		}
	}

	if (chebyshev_vector_array != nullptr) {
		delete[] chebyshev_vector_array;
		chebyshev_vector_array = nullptr;
	}
}

TEMPLATE
void CHEBYSHEV_QUAL::initialCHEBYSHEV_SHARE(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& chebyshev_share) {
	chebyshev_share.segmentNum = input_argument.point_dimension;

	chebyshev_share.t = new double[input_argument.time_series_length];
	fill_n(chebyshev_share.t, input_argument.time_series_length, NULL);

	chebyshev_share.interval = new double[input_argument.time_series_length << 1];
	fill_n(chebyshev_share.interval, input_argument.time_series_length << 1, NULL);

	chebyshev_share.root_tj = new double[input_argument.time_series_length];
	fill_n(chebyshev_share.root_tj, input_argument.time_series_length, NULL);

	chebyshev_share.root_tj_index = new int[input_argument.time_series_length];
	fill_n(chebyshev_share.root_tj_index, input_argument.time_series_length, NULL);

	chebyshev_share.function_coefficient = new double[input_argument.time_series_length];
	fill_n(chebyshev_share.function_coefficient, input_argument.time_series_length, NULL);

	/*chebyshev_share.f = new double[input_argument.time_series_length];
	fill_n(chebyshev_share.f, input_argument.time_series_length, NULL);*/

	chebyshev_share.polynomials = new double[input_argument.time_series_length];
	fill_n(chebyshev_share.polynomials, input_argument.time_series_length, NULL);

	/*chebyshev_share.coefficient = new double[input_argument.time_series_length];
	fill_n(chebyshev_share.coefficient, input_argument.time_series_length, NULL);*/
}

TEMPLATE
void CHEBYSHEV_QUAL::initialCHEBYSHEV(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV& const chebyshev) {
#ifdef _DEBUG
	assert(input_argument.time_series_length != INF && input_argument.degree_m != INF && input_argument.point_dimension != INF);
#endif

	chebyshev.segmentNum = input_argument.degree_m + 1;

	int f_length = input_argument.degree_m + 1 > input_argument.time_series_length ? input_argument.degree_m + 1 : input_argument.time_series_length; //if m+1>n

	chebyshev.f = new DataType[f_length];
	fill_n(chebyshev.f, f_length, INF);

	chebyshev.coefficient = new DataType[input_argument.degree_m + 1];//n=m+1; the number of Chebyshev coefficients
	fill_n(chebyshev.coefficient, input_argument.degree_m + 1, INF);
}

TEMPLATE
void CHEBYSHEV_QUAL::initialCHEBYSHEV(const int& const time_series_length, const int& point_dimension, CHEBYSHEV& const chebyshev) {
#ifdef _DEBUG
	assert(time_series_length != INF);
#endif
	chebyshev.segmentNum = point_dimension;

	int f_length = chebyshev.segmentNum > time_series_length ? chebyshev.segmentNum : time_series_length; //if m+1>trajectory length

	chebyshev.f = new DataType[f_length];
	fill_n(chebyshev.f, f_length, INF);

	chebyshev.coefficient = new DataType[chebyshev.segmentNum];//n=m+1; the number of Chebyshev coefficients
	fill_n(chebyshev.coefficient, chebyshev.segmentNum, INF);
}

TEMPLATE
void CHEBYSHEV_QUAL::initialCHEBYSHEVArray(const int& array_number, CHEBYSHEV*& const chebyshev_array) {
	chebyshev_array = new CHEBYSHEV[array_number];

	for (int array_id = 0; array_id < array_number; array_id++) {
		//chebyshev_array[i].segmentNum = input_argument.point_dimension;
		//chebyshev_array[i].f = new double[input_argument.time_series_length];
		//chebyshev_array[i].coefficient = new double[input_argument.degree_m+1];//n=m+1; the number of Chebyshev coefficients

		initialCHEBYSHEV(input_argument, chebyshev_array[array_id]);
	}
}

TEMPLATE
void CHEBYSHEV_QUAL::initialCHEBYSHEVArray(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV*& const chebyshev_array) {
	chebyshev_array = new CHEBYSHEV[input_argument.arity_d];

	for (int array_id = 0; array_id < input_argument.arity_d; array_id++) {
		//chebyshev_array[i].segmentNum = input_argument.point_dimension;
		//chebyshev_array[i].f = new double[input_argument.time_series_length];
		//chebyshev_array[i].coefficient = new double[input_argument.degree_m+1];//n=m+1; the number of Chebyshev coefficients

		initialCHEBYSHEV(input_argument, chebyshev_array[array_id]);
	}
}

TEMPLATE
void CHEBYSHEV_QUAL::initialCHEBYSHEVArray(const int& trajectory_number, const int& const time_series_length, const int& point_dimension, CHEBYSHEV*& const chebyshev_array) {
	assert(point_dimension);
	chebyshev_array = new CHEBYSHEV[trajectory_number];

	for (int array_id = 0; array_id < trajectory_number; array_id++) {
		//chebyshev_array[i].segmentNum = input_argument.point_dimension;
		//chebyshev_array[i].f = new double[input_argument.time_series_length];
		//chebyshev_array[i].coefficient = new double[input_argument.degree_m+1];//n=m+1; the number of Chebyshev coefficients

		initialCHEBYSHEV(time_series_length, point_dimension, chebyshev_array[array_id]);
	}
}

TEMPLATE
void CHEBYSHEV_QUAL::deleteCHEBYSHEV_SHARE(CHEBYSHEV_SHARE& chebyshev_share) {
	chebyshev_share.segmentNum = NULL;

	delete[] chebyshev_share.t;
	chebyshev_share.t = nullptr;

	delete[] chebyshev_share.interval;
	chebyshev_share.interval = nullptr;

	if (chebyshev_share.root_tj != nullptr) {
		delete[] chebyshev_share.root_tj;
		chebyshev_share.root_tj = nullptr;
	}

	if (chebyshev_share.root_tj_index != nullptr) {
		delete[] chebyshev_share.root_tj_index;
		chebyshev_share.root_tj_index = nullptr;
	}

	if (chebyshev_share.function_coefficient != nullptr) {
		delete[] chebyshev_share.function_coefficient;
		chebyshev_share.function_coefficient = nullptr;
	}

	/*delete[] chebyshev_share.f;
	chebyshev_share.f = nullptr;*/

	delete[] chebyshev_share.polynomials;
	chebyshev_share.polynomials = nullptr;

	/*if (chebyshev_share.coefficient != nullptr) {
		delete[] chebyshev_share.coefficient;
		chebyshev_share.coefficient = nullptr;
	}*/
}

TEMPLATE
void CHEBYSHEV_QUAL::deleteCHEBYSHEV(CHEBYSHEV& const chebyshev) {
	chebyshev.segmentNum = INF;

	if (chebyshev.f != nullptr) {
		delete[] chebyshev.f;
		chebyshev.f = nullptr;
	}

	if (chebyshev.coefficient != nullptr) {
		delete[] chebyshev.coefficient;
		chebyshev.coefficient = nullptr;
	}
}

TEMPLATE
void CHEBYSHEV_QUAL::deleteCHEBYSHEVArray(const int& array_number, CHEBYSHEV*& const chebyshev_array) {
	for (int i = 0; i < array_number; i++) {
		deleteCHEBYSHEV(chebyshev_array[i]);

		delete[] chebyshev_array[i].f;
		chebyshev_array[i].f = nullptr;

		if (chebyshev_array[i].coefficient != nullptr) {
			delete[] chebyshev_array[i].coefficient;
			chebyshev_array[i].coefficient = nullptr;
		}
	}
	delete[] chebyshev_array;
	chebyshev_array = nullptr;
}

TEMPLATE//[-1,1]
bool CHEBYSHEV_QUAL::normalizeOneOne(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
	/*double* max_value = max_element(original_time_series, original_time_series + input_argument.time_series_length);
	double* min_value = min_element(original_time_series, original_time_series + input_argument.time_series_length);*/

	double max_value = input_argument.time_series_length - 1;
	double min_value = 0;

	/*double max_value = input_argument.time_series_length;
	double min_value = 1;*/

	double max_min_coefficient = 2.0 / (max_value - min_value);

	for (int t = 0; t < input_argument.time_series_length; t++) {
		chebyshev_share.t[t] = max_min_coefficient * (t - min_value) - 1.0;
	}

	/*cout << "ti: "<<endl;
	APCA_KNN_QUAL::printArray(chebyshev_share.t, input_argument.time_series_length);*/

	return true;
}

TEMPLATE//Ii
void CHEBYSHEV_QUAL::getInterval(const typename TOOL::INPUT_ARGUMENT& const input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
	assert(chebyshev_share.t[input_argument.time_series_length - 1] != NULL);

	//double* interval_test = new double[input_argument.time_series_length];

	chebyshev_share.interval[0] = -1.0;
	chebyshev_share.interval[1] = (chebyshev_share.t[0] + chebyshev_share.t[1]) / 2.0;

	/*interval_test[0] = -1;

	for (int inerval_id = 0; inerval_id < input_argument.time_series_length; inerval_id++) {
		interval_test[inerval_id]= (chebyshev_share.t[inerval_id] + chebyshev_share.t[i]) / 2.0
	}*/

	for (int i = 1; i <= input_argument.time_series_length - 2; i++) {
		chebyshev_share.interval[i * 2] = (chebyshev_share.t[i - 1] + chebyshev_share.t[i]) / 2.0;
		chebyshev_share.interval[i * 2 + 1] = (chebyshev_share.t[i] + chebyshev_share.t[i + 1]) / 2.0;
	}

	chebyshev_share.interval[input_argument.time_series_length * 2 - 2] = (chebyshev_share.t[input_argument.time_series_length - 2] + chebyshev_share.t[input_argument.time_series_length - 1]) / 2.0;
	chebyshev_share.interval[input_argument.time_series_length * 2 - 1] = 1.0;

	/*cout << "Ii: \n";
	APCA_KNN_QUAL::printArray(chebyshev_share.interval, input_argument.time_series_length<<1);*/
	//delete[] interval_test;
}

TEMPLATE//root tj
void CHEBYSHEV_QUAL::getRoots(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
	double temp_coefficient = boost::math::constants::pi<double>() / double(input_argument.time_series_length);

	for (int root_index = 1; root_index <= input_argument.time_series_length; root_index++) {
		chebyshev_share.root_tj[root_index - 1] = cos((double(root_index) - 0.5) * temp_coefficient);
	}

	/*double* root_test = new double[input_argument.time_series_length];

	for (int root_index = 0; root_index < input_argument.time_series_length; root_index++) {
		root_test[root_index] = cos((double(root_index) - 0.5)*temp_coefficient);
	}*/

	/*cout << "root_test: \n";
	APCA_KNN_QUAL::printArray(root_test, input_argument.time_series_length);*/

	/*cout<<"root_tj: \n";
	APCA_KNN_QUAL::printArray(chebyshev_share.root_tj, input_argument.time_series_length);*/

	/*delete[] root_test;*/
}

TEMPLATE//every root in which interval; Binary search
void CHEBYSHEV_QUAL::getRootsIndex(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
	assert(chebyshev_share.root_tj[input_argument.time_series_length - 1] != NULL);
	assert(chebyshev_share.root_tj_index != nullptr);

	double* index_evaluation = new double[input_argument.time_series_length];
	getRootsIndex1(input_argument, chebyshev_share, index_evaluation);

	int left_index = NULL;
	int right_index = NULL;
	int middle_index = NULL;

	for (int i = 0; i < input_argument.time_series_length; i++) {
		if (chebyshev_share.root_tj[i] == 1) {
			chebyshev_share.root_tj_index[i] = input_argument.time_series_length - 1;
			continue;
		}

		left_index = 0;
		right_index = input_argument.time_series_length - 1;
		middle_index = NULL;
		while (left_index <= right_index) {
			middle_index = left_index + (right_index - left_index) / 2;
			assert(middle_index == (left_index + right_index) >> 1);//test

			if (chebyshev_share.root_tj[i] >= chebyshev_share.interval[middle_index << 1] && chebyshev_share.root_tj[i] < chebyshev_share.interval[middle_index * 2 + 1]) {
				//assert(middle_index >> 1 == middle_index / 2);
				chebyshev_share.root_tj_index[i] = middle_index;
				assert(middle_index < input_argument.time_series_length&& middle_index >= 0);
				break;
			}
			else if (chebyshev_share.root_tj[i] < chebyshev_share.interval[middle_index << 1]) {
				right_index = middle_index - 1;
			}
			else if (chebyshev_share.root_tj[i] >= chebyshev_share.interval[middle_index * 2 + 1]) {
				left_index = middle_index + 1;
			}
			else {
				assert(0);
			}
		}
	}

	//assert
	for (int interval_id = 0; interval_id < input_argument.time_series_length; interval_id++) {
		assert(chebyshev_share.root_tj_index[interval_id] == index_evaluation[interval_id]);
	}

	/*cout << "root_index: \n";
	APCA_KNN_QUAL::printArray(chebyshev_share.root_tj_index, input_argument.time_series_length);*/

	delete[] index_evaluation;
	index_evaluation = nullptr;
}

TEMPLATE//every root in which interval normal search
void CHEBYSHEV_QUAL::getRootsIndex1(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share, double*& const index_evaluation) {
	assert(chebyshev_share.root_tj[input_argument.time_series_length - 1] != NULL);
	assert(chebyshev_share.root_tj_index != nullptr);

	int left_index = NULL;
	int right_index = NULL;
	int middle_index = NULL;

	for (int root_id = 0; root_id < input_argument.time_series_length; root_id++) {
		for (int interval_id = 0; interval_id < input_argument.time_series_length; interval_id++) {
			if (chebyshev_share.root_tj[root_id] >= chebyshev_share.interval[interval_id * 2] && chebyshev_share.root_tj[root_id] < chebyshev_share.interval[interval_id * 2 + 1]) {
				index_evaluation[root_id] = interval_id;
			}
		}
	}

	/*cout << "root_index1: \n";
	APCA_KNN_QUAL::printArray(chebyshev_share.root_tj_index, input_argument.time_series_length);*/
}

TEMPLATE// get internal coefficient of f(x), to speed up program
void CHEBYSHEV_QUAL::getFunctionCoefficient(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
	assert(chebyshev_share.t[input_argument.time_series_length - 1] != NULL);
	assert(chebyshev_share.interval[input_argument.time_series_length * 2 - 1] != NULL);
	assert(chebyshev_share.root_tj[input_argument.time_series_length - 1] != NULL);
	double w_t = NULL;//w(t) Chebyshev weight function
	double I_length = NULL;//|Ii|
	int subinterval_index = NULL;

	for (int i = 0; i < input_argument.time_series_length; i++) {
		w_t = 1.0 / sqrt(1.0 - chebyshev_share.root_tj[i] * chebyshev_share.root_tj[i]);
		//cout << "w(t): " << w_t << endl;
		//cout << chebyshev_share.interval[chebyshev_share.root_tj_index[i] * 2 + 1] << " "<< chebyshev_share.interval[chebyshev_share.root_tj_index[i] << 1] <<endl;
		I_length = chebyshev_share.interval[chebyshev_share.root_tj_index[i] * 2 + 1] - chebyshev_share.interval[chebyshev_share.root_tj_index[i] << 1];
		//cout << "|Ii|: " << I_length << endl;
		chebyshev_share.function_coefficient[i] = sqrt(w_t * I_length);
	}

	/*cout << "f(tj)_coefficient:\n";
	APCA_KNN_QUAL::printArray(chebyshev_share.function_coefficient, input_argument.time_series_length);*/
}

TEMPLATE
void CHEBYSHEV_QUAL::getCHEBYSHEV_SHARE(const typename TOOL::INPUT_ARGUMENT& input_argument, CHEBYSHEV_SHARE& const chebyshev_share) {
	assert(input_argument.point_dimension == chebyshev_share.segmentNum);
	assert(input_argument.time_series_length > 1);
	assert(chebyshev_share.t != nullptr);

	//t(i) [-1,1]
	normalizeOneOne(input_argument, chebyshev_share);

	//I(i)
	getInterval(input_argument, chebyshev_share);

	//root_tj
	getRoots(input_argument, chebyshev_share);

	//root_index
	getRootsIndex(input_argument, chebyshev_share);

	////root_index1
	//getRootsIndex1(input_argument, chebyshev_share);

	//f(t) coefficient
	getFunctionCoefficient(input_argument, chebyshev_share);
}

TEMPLATE//f(t0)-f(tN-1)interval function for chebyshev coefficient, f(t)=g(t)/sqrt(w(t)*|Ii|)
void CHEBYSHEV_QUAL::getFunction(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& chebyshev_share, CHEBYSHEV& const chebyshev) {
	assert(chebyshev_share.t[input_argument.time_series_length - 1] != NULL);
	assert(chebyshev_share.interval[input_argument.time_series_length * 2 - 1] != NULL);
	assert(chebyshev_share.root_tj[input_argument.time_series_length - 1] != NULL);
	assert(chebyshev_share.function_coefficient[input_argument.time_series_length - 1] != NULL);

	//for (int point_number_index = 0; point_number_index < input_argument.point_number; point_number_index++) {
	for (int time_series_index = 0; time_series_index < input_argument.time_series_length; time_series_index++) {//0<= index <N
		//cout << original_time_series[chebyshev_share.root_tj_index[time_series_index]] << endl;
		chebyshev.f[time_series_index] = original_time_series[chebyshev_share.root_tj_index[time_series_index]] / chebyshev_share.function_coefficient[time_series_index];
		//chebyshev.f[time_series_index] = original_time_series[chebyshev_share.root_tj_index[time_series_index]];
	}
	//}
	/*cout << "f(tj):\n";
	APCA_KNN_QUAL::printArray(chebyshev.f, input_argument.time_series_length);*/
}

TEMPLATE//interval function for approximation, difference is f(t)=g(t)
void CHEBYSHEV_QUAL::getFunctionApproximation(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& chebyshev_share, CHEBYSHEV& const chebyshev) {
	assert(chebyshev_share.t[input_argument.time_series_length - 1] != NULL);
	assert(chebyshev_share.interval[input_argument.time_series_length * 2 - 1] != NULL);
	assert(chebyshev_share.root_tj[input_argument.time_series_length - 1] != NULL);
	assert(chebyshev_share.function_coefficient[input_argument.time_series_length - 1] != NULL);

	//for (int point_number_index = 0; point_number_index < input_argument.point_number; point_number_index++) {
	for (int time_series_index = 0; time_series_index < input_argument.time_series_length; time_series_index++) {//0<= index <N
	//cout << original_time_series[chebyshev_share.root_tj_index[time_series_index]] << endl;
	//chebyshev.f[time_series_index] = original_time_series[chebyshev_share.root_tj_index[time_series_index]] / chebyshev_share.function_coefficient[time_series_index];
		chebyshev.f[time_series_index] = original_time_series[chebyshev_share.root_tj_index[time_series_index]];
	}
	//}
	/*cout << "f(tj):\n";
	APCA_KNN_QUAL::printArray(chebyshev.f, input_argument.time_series_length);*/
}

TEMPLATE
void CHEBYSHEV_QUAL::getSpecificCoefficient(const typename TOOL::INPUT_ARGUMENT& input_argument, const int& coefficient_index, const CHEBYSHEV_SHARE& chebyshev_share, CHEBYSHEV& const chebyshev) {
	double temp_coefficient = 2.0 / double(input_argument.time_series_length);
	double sum = 0;

	//cout << "Pi: ";
	for (int time_series_index = 0; time_series_index < input_argument.time_series_length; time_series_index++) {
		sum += chebyshev.f[time_series_index] * boost::math::chebyshev_t(coefficient_index, chebyshev_share.root_tj[time_series_index]);
		//cout << sum << endl;
		//cout << boost::math::chebyshev_t(coefficient_index, chebyshev_share.root_tj[time_series_index]) << ", ";
	}
	//cout << endl <<"sum: " <<sum << endl;
	//chebyshev.coefficient[coefficient_index] = coefficient_index > 0 ? temp_coefficient * sum : temp_coefficient * sum / 2.0;

	chebyshev.coefficient[coefficient_index] = temp_coefficient * sum;

	if (coefficient_index == 0) {
		chebyshev.coefficient[coefficient_index] = chebyshev.coefficient[coefficient_index] / 2.0;
	}

	//cout << "c"<< coefficient_index<<": "<< chebyshev[point_number_index].coefficient[coefficient_index] <<endl;
}

TEMPLATE//1 Notice degree_m or N  2 number of ceofficient is N+1
void CHEBYSHEV_QUAL::getAllCoefficient(const typename TOOL::INPUT_ARGUMENT& input_argument, const CHEBYSHEV_SHARE& chebyshev_share, CHEBYSHEV& const chebyshev) {
	assert(chebyshev.f[input_argument.time_series_length - 1] != NULL);

	for (int coefficient_index = 0; coefficient_index <= input_argument.degree_m; coefficient_index++) {
		getSpecificCoefficient(input_argument, coefficient_index, chebyshev_share, chebyshev);
	}

	/*cout << "ci: ";
	APCA_KNN_QUAL::printArray(chebyshev.coefficient, input_argument.degree_m +1);*/
}

TEMPLATE//specific coefficient ci for approximation
void CHEBYSHEV_QUAL::getSpecificCoefficientApproximation(const typename TOOL::INPUT_ARGUMENT& input_argument, const int& coefficient_index, const CHEBYSHEV_SHARE& chebyshev_share, DataType*& const original_time_series, CHEBYSHEV& const chebyshev) {
	double temp_coefficient = 2.0 / double(input_argument.time_series_length);
	double sum = 0;
	//double weight_function = NULL;
	//double I_length = NULL;
	//cout << "Pi: ";
	for (int time_series_index = 0; time_series_index < input_argument.time_series_length; time_series_index++) {
		//weight_function= 1.0 / sqrt(1.0 - chebyshev_share.root_tj[coefficient_index] * chebyshev_share.root_tj[coefficient_index]);
		//I_length = chebyshev_share.interval[chebyshev_share.root_tj_index[time_series_index] * 2 + 1] - chebyshev_share.interval[chebyshev_share.root_tj_index[time_series_index] << 1];
		//cout<< I_length<<endl;
		//sum += original_time_series[chebyshev_share.root_tj_index[time_series_index]] * boost::math::chebyshev_t(coefficient_index, chebyshev_share.root_tj[time_series_index]);// *I_length;
		sum += original_time_series[chebyshev_share.root_tj_index[time_series_index]] * boost::math::chebyshev_t(coefficient_index, chebyshev_share.root_tj[time_series_index]);
		//sum += chebyshev.f[time_series_index] * boost::math::chebyshev_t(coefficient_index, chebyshev_share.root_tj[time_series_index]);

		//cout << boost::math::chebyshev_t(coefficient_index, chebyshev_share.root_tj[time_series_index]) << ", ";
	}
	//cout << endl <<"sum: " <<sum << endl;
	//chebyshev.coefficient[coefficient_index] = coefficient_index > 0 ? temp_coefficient * sum : temp_coefficient * sum / 2.0;

	chebyshev.coefficient[coefficient_index] = temp_coefficient * sum;

	if (coefficient_index == 0) {
		chebyshev.coefficient[coefficient_index] = chebyshev.coefficient[coefficient_index] / 2.0;
	}

	//cout << "c"<< coefficient_index<<": "<< chebyshev[point_number_index].coefficient[coefficient_index] <<endl;
}

TEMPLATE//All coefficient c1 - cN for approximation
void CHEBYSHEV_QUAL::getAllCoefficientApproximation(const typename TOOL::INPUT_ARGUMENT& input_argument, const CHEBYSHEV_SHARE& chebyshev_share, DataType*& const original_time_series, CHEBYSHEV& const chebyshev) {
	assert(chebyshev.f[input_argument.time_series_length - 1] == INF || chebyshev.f[input_argument.time_series_length - 1] == NULL);
	//cout<<"w(t): \n"<<endl;
	for (int coefficient_index = 0; coefficient_index <= input_argument.degree_m; coefficient_index++) {
		getSpecificCoefficientApproximation(input_argument, coefficient_index, chebyshev_share, original_time_series, chebyshev);
	}
	//cout << endl;
}

TEMPLATE//get chebyshev
void CHEBYSHEV_QUAL::getCHEBYSHEV(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, CHEBYSHEV& const chebyshev) {
	getFunction(input_argument, original_time_series, chebyshev_share, chebyshev);
	getAllCoefficient(input_argument, chebyshev_share, chebyshev);
}

TEMPLATE//get chebyshev
void CHEBYSHEV_QUAL::getMultiCHEBYSHEV(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, CHEBYSHEV& const chebyshev) {
	assert(0);
	DataType* single_trajectory = new DataType[input_argument.time_series_length];

	CHEBYSHEV temp_chebyshev;
	initialCHEBYSHEV(input_argument, temp_chebyshev);

	//TOOL::printArray(original_time_series, input_argument.point_multi_single_length);

	for (int arity_id = 0; arity_id < input_argument.arity_d; arity_id++) {
		copy(original_time_series + input_argument.time_series_length * arity_id, original_time_series + input_argument.time_series_length * arity_id + input_argument.time_series_length, single_trajectory);

		getFunction(input_argument, single_trajectory, chebyshev_share, temp_chebyshev);
		getAllCoefficient(input_argument, chebyshev_share, temp_chebyshev);
		//TOOL::printArray(temp_chebyshev.coefficient, input_argument.point_dimension);

		copy_n(temp_chebyshev.coefficient, temp_chebyshev.segmentNum, chebyshev.coefficient + input_argument.point_dimension * arity_id);
	}

	//TOOL::printArray(chebyshev.coefficient, input_argument.point_multi_single_dimension);

	delete[] single_trajectory;
	single_trajectory = nullptr;

	deleteCHEBYSHEV(temp_chebyshev);
}

TEMPLATE//get chebyshev approximation. difference is f(t)=g(t)
void CHEBYSHEV_QUAL::getChebyshevApproximation(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, CHEBYSHEV& const chebyshev) {
	//getFunctionApproximation(input_argument, original_time_series, chebyshev_share, chebyshev);
	getAllCoefficientApproximation(input_argument, chebyshev_share, original_time_series, chebyshev);
}

TEMPLATE//get chebyshev interation degree_m increase
void CHEBYSHEV_QUAL::getDeviationIterationChebyshev(typename  TOOL::INPUT_ARGUMENT& const input_argument, string*& const write_file_name, DataType*& const original_time_series, CHEBYSHEV_SHARE& const chebyshev_share, const int& degree_m_begin, const int& degree_m_end, const int& interval) {
	double deviation_sum = 0;
	double deviation_max = 0;

	/*double average = NULL;
	double variance = NULL;*/

	DataType* original_normalization = new DataType[input_argument.time_series_length];
	fill_n(original_normalization, input_argument.time_series_length, NULL);
	DataType* chebyshev_approximation = new DataType[input_argument.time_series_length];
	fill_n(chebyshev_approximation, input_argument.time_series_length, NULL);

	cout << "Original time series: " << endl;
	APCA_KNN_QUAL::printArray(original_time_series, input_argument.time_series_length);

	APCA_KNN_QUAL::normalizeStandard(original_time_series, input_argument.time_series_length, original_normalization);

	cout << "Normalized Original time series: " << endl;
	APCA_KNN_QUAL::printArray(original_normalization, input_argument.time_series_length);

	for (input_argument.degree_m = degree_m_begin; input_argument.degree_m <= degree_m_end; input_argument.degree_m += interval) {
		cout << "degree_m: " << input_argument.degree_m << endl;

		approximateOriginalFunction(input_argument, original_normalization, chebyshev_share, chebyshev_approximation);
		//approximateOriginalFunction(input_argument, original_time_series, chebyshev_share, chebyshev_approximation);

		cout << "Chebyshev approximation: " << endl;
		APCA_KNN_QUAL::printArray(chebyshev_approximation, input_argument.time_series_length);

		/*cout << "normalized Chebyshev approximation: " << endl;
		APCA_KNN_QUAL::printArray(chebyshev_normalization, input_argument.time_series_length);*/

		APCA_KNN_QUAL::getReconstructionError(original_normalization, chebyshev_approximation, input_argument.time_series_length, deviation_sum, deviation_max);
		//APCA_KNN_QUAL::getReconstructionError(original_time_series, chebyshev_approximation, input_argument.time_series_length, deviation_sum, deviation_max);

		TOOL::writeSingleResult(input_argument, write_file_name[0], deviation_sum);
		TOOL::writeSingleResult(input_argument, write_file_name[1], deviation_max);

		//writeApproximationResult(input_argument, chebyshev_normalization, deviation_sum,deviation_max);
	}

	delete[] original_normalization;
	original_normalization = nullptr;
	delete[] chebyshev_approximation;
	chebyshev_approximation = nullptr;
}

TEMPLATE//Loop to get pruning pwoer, execution time, I/O cost
void CHEBYSHEV_QUAL::getPrunIterationChebyshev(typename TOOL::INPUT_ARGUMENT& const input_argument, string*& const write_file_name, DataType*& const original_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, const int& degree_m_begin, const int& degree_m_end, const int& interval) {
	//for (K = 5; K < 21; K += 5) {
	//	for (n = 256; n < 1025; n <<= 1) {
	//		for (N = 16; N < 65; N <<= 1) {
	//			//DataType original[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
	//			//DataType query_time_series[] = {10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
	//			PLA_QUAL pla(n, N, point_number, max_node, K, read_file_name);
	//			DataType* query_time_series = new DataType[n];
	//			CAPCA_KNN<DataType> knn;
	//			knn.getFileStreamByID("CinC_ECG_torso_TRAIN", n, 3, query_time_series);
	//			priority_queue<PLA_QUAL::ORIGINAL_TIME_SERIES_PAIR, vector<PLA_QUAL::ORIGINAL_TIME_SERIES_PAIR>, PLA_QUAL::priorityDistanceEUC > base_queue;
	//			list<PLA_QUAL::ORIGINAL_TIME_SERIES_PAIR> result;

	//			CPLA<DataType>::PLA* pla_array_accumulate = new CPLA<DataType>::PLA[pla.input_argument.point_number];
	//			RTREE apcaRTree(pla.input_argument.point_dimension * 2, pla.input_argument.rtree_max_nodes);
	//			pla.initialPLAArray(pla.input_argument, pla_array_accumulate);
	//
	//			pla.buidRTreeIndex(pla.input_argument, apcaRTree, pla_array_accumulate, read_file_name);
	//			pla.PLAKNNSearch(pla.input_argument, query_time_series, apcaRTree, pla_array_accumulate, K, read_file_name, result);
	//
	//			apcaRTree.RemoveAll();
	//			pla.deletePLAArrary(pla.input_argument, pla_array_accumulate);

	//			pla.SimpleBaseKNNSearch(pla.input_argument, query_time_series, K, read_file_name, base_queue);

	//			delete[] query_time_series;
	//			query_time_series = nullptr;

	//			pla.writeResult(pla.input_argument, write_file_name, result, base_queue);

	//			priority_queue<PLA_QUAL::ORIGINAL_TIME_SERIES_PAIR, vector<PLA_QUAL::ORIGINAL_TIME_SERIES_PAIR>, PLA_QUAL::priorityDistanceEUC >().swap(base_queue);
	//			result.clear();

	//		}
	//	}
	//}
}

TEMPLATE//get chebyshev //191107 notes: needs improvement, maybe directly insert into Rtree
void CHEBYSHEV_QUAL::getChebyshevMBR(const CHEBYSHEV& const chebyshev, CHEBYSHEV& const chebyshev_MBR) {
	//copy_n(chebyshev.coefficient, chebyshev_MBR.segmentNum, chebyshev_MBR.f);          //min
	//copy_n(chebyshev.coefficient, chebyshev_MBR.segmentNum, chebyshev_MBR.coefficient);//max

	for (int i = 0; i < chebyshev_MBR.segmentNum; i++) {
		chebyshev_MBR.f[i] = chebyshev_MBR.coefficient[i] = chebyshev.coefficient[i];
	}

	/*cout << "mbr:" << endl;
	TOOL::printArray(chebyshev_MBR.coefficient, input_argument.degree_m + 1);
	TOOL::printArray(chebyshev_MBR.f, input_argument.degree_m + 1);*/
}

TEMPLATE//approximate original function
DataType& CHEBYSHEV_QUAL::getSpecificApproximation(typename TOOL::INPUT_ARGUMENT& input_argument, const CHEBYSHEV& const chebyshev, CHEBYSHEV_SHARE& const chebyshev_share, const int& normalized_array_id, DataType& const specific_approximation) {
	specific_approximation = 0;
	for (int degree_id = 0; degree_id <= input_argument.degree_m; degree_id++) {
		specific_approximation += chebyshev.coefficient[degree_id] * boost::math::chebyshev_t(degree_id, chebyshev_share.t[normalized_array_id]);
		//cout << chebyshev.coefficient[degree_id] << endl;
		//specific_approximation += chebyshev.coefficient[degree_id] * boost::math::chebyshev_t(degree_id, chebyshev_share.root_tj[normalized_array_id]);
		//specific_approximation += chebyshev.coefficient[degree_id] * boost::math::chebyshev_t(degree_id, normalized_array_id);
#ifdef _DEBUG
		assert(t_polynomial_value(degree_id, chebyshev_share.root_tj[normalized_array_id]) == boost::math::chebyshev_t(degree_id, chebyshev_share.root_tj[normalized_array_id]));
#endif
	}
	return specific_approximation;
}

TEMPLATE//approximate original function
void CHEBYSHEV_QUAL::approximateOriginalFunction(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, CHEBYSHEV_SHARE& const chebyshev_share, DataType*& const approximated_time_series) {
	CHEBYSHEV chebyshev;
	initialCHEBYSHEV(input_argument, chebyshev);
	double specific_approximation = NULL;

	//getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev);
	getChebyshevApproximation(input_argument, original_time_series, chebyshev_share, chebyshev);

	//cout << "coefficient: " << endl;
	//TOOL::printArray(chebyshev.coefficient, input_argument.degree_m + 1);

	for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		//specific_approximation = getSpecificApproximation(chebyshev, array_id, specific_approximation);

		approximated_time_series[array_id] = getSpecificApproximation(input_argument, chebyshev, chebyshev_share, array_id, specific_approximation);
	}

	/*cout << "approximated function: " << endl;
	APCA_KNN_QUAL::printArray(approximated_time_series, input_argument.time_series_length);*/

	deleteCHEBYSHEV(chebyshev);
}

//191206
TEMPLATE
double CHEBYSHEV_QUAL::get_sum_deviation(typename TOOL::INPUT_ARGUMENT& input_argument, const vector<DataType>& const original_time_series_vector, CHEBYSHEV_SHARE& const chebyshev_share) {
	assert(input_argument.time_series_length != INF);
	double deviation_sum = INF;
	CHEBYSHEV chebyshev;
	initialCHEBYSHEV(input_argument, chebyshev);
	double specific_approximation = NULL;

	DataType* original_time_series = new DataType[input_argument.time_series_length];
	vector<DataType> approximated_time_series(input_argument.time_series_length, INF);

	copy_n(original_time_series_vector.begin(), original_time_series_vector.size(), original_time_series);

	//getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev);
	getChebyshevApproximation(input_argument, original_time_series, chebyshev_share, chebyshev);

	/*cout << "coefficient: " << endl;
	APCA_KNN_QUAL::printArray(chebyshev.coefficient, input_argument.degree_m + 1);*/

	for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		//specific_approximation = getSpecificApproximation(chebyshev, array_id, specific_approximation);
		assert(original_time_series_vector[array_id] == original_time_series[array_id]);
		approximated_time_series[array_id] = getSpecificApproximation(input_argument, chebyshev, chebyshev_share, array_id, specific_approximation);
	}

	/*cout << "approximated function: " << endl;
	APCA_KNN_QUAL::printArray(approximated_time_series, input_argument.time_series_length);*/

	deviation_sum = TOOL::getDeviation(original_time_series_vector, approximated_time_series);
	assert(deviation_sum != INF);

	delete[] original_time_series;
	original_time_series = nullptr;

	approximated_time_series.clear();
	approximated_time_series.shrink_to_fit();

	deleteCHEBYSHEV(chebyshev);

	return deviation_sum;
}

TEMPLATE
template<typename T, typename Y, typename U>
long double CHEBYSHEV_QUAL::get_sum_deviation(T& input_argument, const vector<Y>& const original_time_series_vector, CHEBYSHEV_SHARE& const chebyshev_share, U& const result_collection) {
	assert(input_argument.time_series_length != INF);
	double deviation_sum = INF;
	CHEBYSHEV chebyshev;
	initialCHEBYSHEV(input_argument, chebyshev);
	double specific_approximation = NULL;

	DataType* original_time_series = new DataType[input_argument.time_series_length];
	vector<DataType> approximated_time_series(input_argument.time_series_length, INF);

	copy_n(original_time_series_vector.begin(), original_time_series_vector.size(), original_time_series);

	//getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev);
	getChebyshevApproximation(input_argument, original_time_series, chebyshev_share, chebyshev);

	for (int array_id = 0; array_id < original_time_series_vector.size(); array_id++) {
		//specific_approximation = getSpecificApproximation(chebyshev, array_id, specific_approximation);
		assert(original_time_series_vector[array_id] == original_time_series[array_id]);
		approximated_time_series[array_id] = getSpecificApproximation(input_argument, chebyshev, chebyshev_share, array_id, specific_approximation);
	}

	/*------ Print Reconstruct Time Series --------*/
	//cout << "approximated function: " << endl;
	//TOOL::print_vector(approximated_time_series);
	/*---------------------------------------------*/

	deviation_sum = TOOL::getDeviation(original_time_series_vector, approximated_time_series);
	assert(deviation_sum != INF && input_argument.point_dimension == chebyshev.segmentNum);

	//
	result_collection.sum_deviation = 0;
	result_collection.max_deviation = 0;
	result_collection.max_deviation_multiple_width = 0;

	long double sum = 0;
	long double difference = NULL;
	long double deviation_max = -INF;
	int point_id = 0;

	/*===================================For CHEBY coefficient====================================================*/
	input_argument.remainder = int(input_argument.time_series_length) % int(input_argument.point_dimension);//For PLA
	input_argument.segment_length_second = (input_argument.time_series_length - input_argument.remainder) / input_argument.point_dimension;
	input_argument.segment_length_first = input_argument.segment_length_second + 1;
	/*==========================================================================================================*/

	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
		deviation_max = -INF;
		if (segment_id < input_argument.remainder) {

			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
				point_id = input_argument.segment_length_first * segment_id + interval_id;// time series id
				//cout << "id: " << point_id <<endl;
				difference = fabs(approximated_time_series[point_id] - original_time_series_vector[point_id]);
				//cout << "difference: "<<difference<<endl;
				deviation_max = max(deviation_max, difference);
				sum += difference * difference;
			}

			result_collection.max_deviation += deviation_max;
			result_collection.max_deviation_multiple_width += deviation_max * input_argument.segment_length_first;
			//indexOfLongSegment = segment_id+1;
			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
		}
		else {
			//assert(input_argument.remainder * input_argument.segment_length_first == endOfLongSegment);
			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
				point_id = input_argument.remainder + input_argument.segment_length_second * segment_id + interval_id;
				//cout << "id: " << point_id << endl;
				difference = fabs(approximated_time_series[point_id] - original_time_series_vector[point_id]);
				//cout << "difference: " << difference << endl;
				deviation_max = max(deviation_max, difference);
				sum += difference * difference;
			}

			result_collection.max_deviation += deviation_max;
			result_collection.max_deviation_multiple_width += deviation_max * input_argument.segment_length_second;
		}
	}

	result_collection.max_deviation_av = result_collection.max_deviation / double(input_argument.point_dimension);
	//result_collection.max_deviation = result_collection.max_deviation_av;

	deviation_sum = sqrt(sum);
	result_collection.sum_deviation = deviation_sum;
	double test_deviation = TOOL::getDeviation(original_time_series_vector, approximated_time_series);
	assert(float(deviation_sum) == float(test_deviation));
	//

	//
	//result_collection.sum_deviation = 0;
	//result_collection.max_deviation = 0;
	//result_collection.max_deviation_multiple_width = 0;

	////int point_id = 0;
	//double segment_len = input_argument.time_series_length / (input_argument.point_dimension-1);
	////long double difference = 0;
	//for (int segment_id = 0; segment_id < input_argument.point_dimension - 1; segment_id++) {
	//	long double deviation_max = 0;
	//	for (int interval_id = 0; interval_id < segment_len; interval_id++) {
	//		point_id = segment_len * segment_id + interval_id;// time series id
	//		difference = fabs(approximated_time_series[point_id] - original_time_series[point_id]);
	//		deviation_max = max(deviation_max, difference);
	//	}
	//	result_collection.max_deviation += deviation_max;
	//	result_collection.max_deviation_multiple_width += deviation_max * segment_len;
	//}
	//result_collection.max_deviation_av = result_collection.max_deviation / double(input_argument.point_dimension);
	//result_collection.sum_deviation = deviation_sum;
	//

	delete[] original_time_series;
	original_time_series = nullptr;

	approximated_time_series.clear();
	approximated_time_series.shrink_to_fit();

	deleteCHEBYSHEV(chebyshev);

	return deviation_sum;
}

TEMPLATE // 0<=i<=m
double& CHEBYSHEV_QUAL::getChebyshevDistance(const typename TOOL::INPUT_ARGUMENT& input_argument, const CHEBYSHEV& chebyshev_array_a, const CHEBYSHEV& chebyshev_array_b, double& const chebyshev_distance) {
	double sum = 0;
	double difference_value = NULL;

	for (int coefficient_index = 0; coefficient_index <= input_argument.degree_m; coefficient_index++) {
		difference_value = chebyshev_array_a.coefficient[coefficient_index] - chebyshev_array_b.coefficient[coefficient_index];
		//cout << "difference_value: " << difference_value << endl;
		sum += difference_value * difference_value;
	}

	chebyshev_distance = sqrt(0.5 * boost::math::constants::pi<double>() * sum);
	/*cout << "chebyshev_distance: " << chebyshev_distance << endl;*/
	return chebyshev_distance;
}

TEMPLATE // 0<=i<=m
double& CHEBYSHEV_QUAL::getChebyshevDistance(const int& segment_number, const CHEBYSHEV& chebyshev_array_a, const CHEBYSHEV& chebyshev_array_b, double& const chebyshev_distance) {
	double sum = 0;
	double difference_value = NULL;

	for (int coefficient_index = 0; coefficient_index < segment_number; coefficient_index++) {
		difference_value = chebyshev_array_a.coefficient[coefficient_index] - chebyshev_array_b.coefficient[coefficient_index];
		//cout << "difference_value: " << difference_value << endl;
		sum += difference_value * difference_value;
	}

	chebyshev_distance = sqrt(0.5 * boost::math::constants::pi<double>() * sum);
	/*cout << "chebyshev_distance: " << chebyshev_distance << endl;*/
	return chebyshev_distance;
}

TEMPLATE
RTREE& CHEBYSHEV_QUAL::buidRTreeIndex(typename  TOOL::INPUT_ARGUMENT& const  input_argument, RTREE& RTree, CHEBYSHEV*& const chebyshev_array) {
	//RTree<DataType, ElementType>& buidRTreeIndex(RTree<DataType, ElementType> &APCARTree, double* g_query_time_series, const double& g_time_series_length, const double& g_index_point_number, double(&test_d_original_time_series)[ROW][COLUMN], Link *APCALinkOriginal) {
	//getCHEBYSHEV_SHARE(input_argument, chebyshev_share);
	assert(input_argument.build_rtree_time == 0.0);
	printf(">>>***Build RTree Index***<<<\n");

	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[0]);//whole build_rtree_time

	int point_id = NULL, j = NULL, f_insert_count = NULL;
	DataType* original_time_series = new DataType[input_argument.time_series_length];

	//PLA pla_array_MBR; //Temp MBR
	//initialPLA(pla_array_MBR, input_argument.point_dimension * 2);

	CHEBYSHEV chebyshev_MBR;
	initialCHEBYSHEV(input_argument, chebyshev_MBR);
	/*cout << "Query Point : ";
	for (i = 0; i < g_time_series_length; i++) cout << g_query_time_series[i] << ", ";
	cout << endl;*/

	//DataType** multi_dimention_trajectories=nullptr;
	//initialMultiDimensionTrajectory(input_argument, multi_dimention_trajectories);

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(input_argument.read_file_name);
	assert(file_stream);

	f_insert_count = 0;
	for (point_id = 0; point_id < input_argument.point_number && (!file_stream.eof()) && file_stream.is_open() && file_stream.good(); point_id++) {
		//getRandomAPCAPoint(APCALinkOriginal[i].originalLink, g_time_series_length);
		//APCALinkOriginal[i].originalLink = test_d_original_time_series[i];
		//getNormalArray(APCALinkOriginal[i].originalLink, g_time_series_length);

		//for (int array_id = 0; array_id < input_argument.airty_d; array_id++) {
		fill_n(original_time_series, input_argument.time_series_length, NULL);

		file_stream >> fs_row_string;
		//memory_account[2] = fs_row_string.size();
		stringstream sstr(fs_row_string);

		int string_id = -1;
		while (getline(sstr, fs_row_number, ',') && string_id < input_argument.time_series_length) {
			if (string_id > -1) {
				original_time_series[string_id] = stod(fs_row_number);
				assert(original_time_series[string_id] != NULL);
				//multi_dimention_trajectories[array_id][string_id] = stod(fs_row_number);
				/*cout << original_time_series[string_id] << ", ";*/
			}
			string_id++;
		}
		/*cout << endl;*/

		/*if (PAA_or_APCA == 0) {
		getAPCAPoint(original_time_series, g_time_series_length, APCARTree.NUMDIMS / 2, APCALinkOriginal[i].APCALink);
		}
		else {
		divideRemainderPAA(original_time_series, APCALinkOriginal[i].APCALink, g_time_series_length, APCARTree.NUMDIMS / 2);
		}*/

		/*getPLA(input_argument, original_time_series, PLA_array_accumulate[point_id]);
		getPLAMBR(PLA_array_accumulate[point_id], pla_array_MBR);*/

		//getCHEBYSHEV(input_argument, multi_dimention_trajectories[array_id], chebyshev_share, chebyshev_array);

		TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);//z-score normalization

		getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev_array[point_id]);
		getChebyshevMBR(chebyshev_array[point_id], chebyshev_MBR);
		/*for (int k = 0; k < pla_array_MBR.segmentNum; k++) {
		cout << "(" << pla_array_MBR.a[k] << "," << pla_array_MBR.b[k] << ") ";
		}
		cout << endl;*/
		//}

		RTree.Insert(chebyshev_MBR.f, chebyshev_MBR.coefficient, point_id);// a min, b max

		f_insert_count++;
	}

	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();

	/*cout << "Root Node : sub node number = " << RTree.m_root->m_count << " Root level = : " << RTree.m_root->m_level << "\n\nBegin to build a RTree:\n";
	cout << "\n RTree conclusion\n The number of RTree Data Point = : " << RTree.Count() << endl;*/

	//deletePLA(pla_array_MBR);
	deleteCHEBYSHEV(chebyshev_MBR);
	//deleteMultiDimensionTrajectory(input_argument, multi_dimention_trajectories);

	delete[] original_time_series;
	original_time_series = nullptr;

	//***** Delete Leaf node memory   ****//
	RTREE::Iterator it;
	for (RTree.GetFirst(it); !RTree.IsNull(it); RTree.GetNext(it)) {
		it.deleteLeafNodeMemory();
	}

	APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[0], input_argument.build_rtree_time);
	//cout << "RTree build time: " << input_argument.build_rtree_time<<" us" << endl;
	TOOL::writeSingleResult(input_argument, input_argument.write_file_name[2], input_argument.build_rtree_time);

	return RTree;
}

TEMPLATE//Range Search
template<typename T>
void CHEBYSHEV_QUAL::rangeSearch(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, T& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue) {
	//cout << "rangeSearch()" << endl;
	int leaf_id = NULL;
	double chebyshev_distance = NULL;
	DataType* original_time_series = new DataType[input_argument.time_series_length];
	ID_DIST id_dist;
	list<ID_DIST> cheby_candidate;
	CHEBYSHEV chebyshev_query;
	initialCHEBYSHEV(input_argument, chebyshev_query);
	getCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);
	/*priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement > queue;*/

	double temp_navigate_time = 0;
	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[4]);//time for navigate index nodes

	RTREE::Iterator it;
	for (RTree.GetFirst(it); !RTree.IsNull(it); RTree.GetNext(it)) {
		leaf_id = RTree.GetAt(it);

		/*cout << "coefficient: ";
		APCA_KNN_QUAL::printArray(chebyshev_array[leaf_id].coefficient, input_argument.degree_m+1);
		cout << "id: " << leaf_id <<" ";*/
		double distance_chebyshev_time = 0;
		APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[5]);//time for chebyshev distance
		if (getChebyshevDistance(input_argument, chebyshev_query, chebyshev_array[leaf_id], chebyshev_distance) < radius) {
			/*cout << "Loop: " << endl;*/
			APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[5], distance_chebyshev_time);
			input_argument.distance_lowbound_time += distance_chebyshev_time;

			APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[4], temp_navigate_time);
			input_argument.navigate_index_time += temp_navigate_time;

			id_dist.trajectory_id = leaf_id;
			APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, leaf_id, original_time_series);
			//id_dist.d_dist = chebyshev_distance;

			double distance_euclidean_time = 0;
			APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[6]);//time for euclidean distance

			id_dist.d_dist = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.time_series_length, original_time_series, input_argument.time_series_length);

			APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[6], distance_euclidean_time);
			input_argument.distance_euc_time += distance_euclidean_time;
			input_argument.sum_distance_euc++;
			cheby_candidate.push_front(id_dist);
			if (id_dist.d_dist < radius) {
				//cout << "   range id: " << id_dist.trajectory_id << ", distance euc: " << id_dist.d_dist << endl;
				queue.push(id_dist);
				cheby_candidate.pop_front();
			}
		}
		//it.deleteLeafNodeMemory();
	}
	//cout << "Range search queue.size: " << queue.size() << endl;
	cheby_candidate.sort([](const ID_DIST& first, const  ID_DIST& second) {return first.d_dist < second.d_dist; });//small to big

	//*************180625xuer improvement: in case queue.size< K.*********************************
	if (queue.size() < input_argument.K) {
		int difference = input_argument.K - queue.size();
		while (difference > 0) {
			queue.push(cheby_candidate.front());
			cheby_candidate.pop_front();
			difference--;
		}
	}
	//*******************************************
	assert(queue.size() >= input_argument.K);
	/*while (!queue.empty()) {
		cout << "id: " << queue.top().trajectory_id << ", dist: " << queue.top().d_dist << endl;
		queue.pop();
	}*/

	delete[] original_time_series;
	original_time_series = nullptr;
	//priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >().swap(queue);
	deleteCHEBYSHEV(chebyshev_query);
}

TEMPLATE//Range Search
template<typename T>
void CHEBYSHEV_QUAL::rangeSearchMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, T& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue) {
	//cout << "rangeSearch()" << endl;
	int leaf_id = NULL;
	double chebyshev_distance = NULL;
	DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
	ID_DIST id_dist;
	list<ID_DIST> cheby_candidate;
	CHEBYSHEV chebyshev_query;

	initialCHEBYSHEV(input_argument.point_multi_single_length, input_argument.point_multi_single_dimension, chebyshev_query);
	//getCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);
	getMultiCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);

	/*priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement > queue;*/
	double temp_navigate_time = 0;


	RTREE::Iterator it;
	for (RTree.GetFirst(it); !RTree.IsNull(it); RTree.GetNext(it)) {
		leaf_id = RTree.GetAt(it);

		/*cout << "coefficient: ";
		APCA_KNN_QUAL::printArray(chebyshev_array[leaf_id].coefficient, input_argument.degree_m+1);
		cout << "id: " << leaf_id <<" ";*/


		double cheby_dist = getChebyshevDistance(input_argument.point_multi_single_dimension, chebyshev_query, chebyshev_array[leaf_id], chebyshev_distance);
		//cout << "cheby dist: " << cheby_dist << ",    radius: " << radius << endl;
		if (cheby_dist < radius) {
			/*cout << "Loop: " << endl;*/


			id_dist.trajectory_id = leaf_id;
			//APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, leaf_id, original_time_series);
			//TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, leaf_id, original_time_series);
			//id_dist.d_dist = chebyshev_distance;

			if (input_argument.read_multiple_file_name) {
				TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, leaf_id, original_time_series);
			}
			else {
				//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id));
				TOOL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, leaf_id, original_time_series);
			}
			TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);

			id_dist.d_dist = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);
			/*=========================Evaluation=================================*/
			input_argument.IO_cost++;
			/*=====================================================================*/

			cheby_candidate.push_front(id_dist);

			if (id_dist.d_dist <= radius) {
				//cout << "   range id: " << id_dist.trajectory_id << ", distance euc: " << id_dist.d_dist << endl;
				queue.push(id_dist);
				cheby_candidate.pop_front();
			}
		}
		//it.deleteLeafNodeMemory();
	}
	//cout << "Range search queue.size: " << queue.size() << endl;
	cheby_candidate.sort([](const ID_DIST& first, const  ID_DIST& second) {return first.d_dist < second.d_dist; });//small to big

	//*************180625xuer improvement: in case queue.size< K.*********************************
	if (queue.size() < input_argument.K && !cheby_candidate.empty()) {
		int difference = input_argument.K - queue.size();
		while (difference > 0) {
			queue.push(cheby_candidate.front());
			cheby_candidate.pop_front();
			difference--;
		}
	}
	//********************************************************************************************
	//assert(queue.size() >= input_argument.K);
	/*while (!queue.empty()) {
		cout << "id: " << queue.top().trajectory_id << ", dist: " << queue.top().d_dist << endl;
		queue.pop();
	}*/

	delete[] original_time_series;
	original_time_series = nullptr;
	//priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >().swap(queue);
	deleteCHEBYSHEV(chebyshev_query);
}

//***************************************************************
// Method:rangeSearchMulti
// Qualifier:  //for priority queue, small to big, top is smallest
// Input:add data_source
// Output: 
// date:210603
// author:
//***************************************************************
TEMPLATE
template<typename T>
void CHEBYSHEV_QUAL::rangeSearchMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const typename TOOL::DATA_SOURCE& const data_source, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, T& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue) {
	//cout << "rangeSearch()" << endl;
	int leaf_id = NULL;
	double chebyshev_distance = NULL;
	DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
	vector<DataType> original_time_series_vector(input_argument.point_multi_single_length, INF);
	vector<DataType> query_time_series_vector(query_time_series, query_time_series + input_argument.point_multi_single_length);

#ifdef _DEBUG
	for (int array_id = 0; array_id < input_argument.point_multi_single_length; array_id++) {
		assert(query_time_series_vector[array_id] == query_time_series[array_id]);
	}
	assert(RTree.m_root->m_count > 0);
#endif

	ID_DIST id_dist;
	list<ID_DIST> cheby_candidate;
	CHEBYSHEV chebyshev_query;

	initialCHEBYSHEV(input_argument.point_multi_single_length, input_argument.point_multi_single_dimension, chebyshev_query);
	getCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);//191224
	//getMultiCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);

	/*priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement > queue;*/
	double temp_navigate_time = 0;


	typename T::Iterator it;
	for (RTree.GetFirst(it); !RTree.IsNull(it); RTree.GetNext(it)) {
		leaf_id = RTree.GetAt(it);

		/*cout << "coefficient: ";
		APCA_KNN_QUAL::printArray(chebyshev_array[leaf_id].coefficient, input_argument.degree_m+1);
		cout << "id: " << leaf_id <<" ";*/


		double cheby_dist = getChebyshevDistance(input_argument.point_multi_single_dimension, chebyshev_query, chebyshev_array[leaf_id], chebyshev_distance);
		//cout << "cheby dist: " << cheby_dist << ",    radius: " << radius << endl;
		if (cheby_dist < radius) {
			/*cout << "Loop: " << endl;*/


			id_dist.trajectory_id = leaf_id;
			//APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, leaf_id, original_time_series);
			//TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, leaf_id, original_time_series);
			//id_dist.d_dist = chebyshev_distance;


			/*-----------------------------------------------------Original Time Series---------------------------------------------------------------------------------------*/
			//already normalized
			input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[3]);
			TOOL::read_normalized_multi_time_series(data_source, leaf_id, original_time_series_vector);
			TOOL::recordStartTime(TOOL::time_record[3]);// whole time
			TOOL::recordStartTime(TOOL::time_record[2]);// KNN total time
#ifdef _DEBUG
			copy_n(original_time_series_vector.begin(), original_time_series_vector.size(), original_time_series);
#endif
			/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

			//if (input_argument.read_multiple_file_name) {
			//	TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, leaf_id, original_time_series);
			//}
			//else {
			//	//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id));
			//	TOOL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, leaf_id, original_time_series);
			//}
			//TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);

			id_dist.d_dist = TOOL::distanceEUC(query_time_series_vector, original_time_series_vector);

#ifdef _DEBUG
			double test_EU_dist = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);
			assert(id_dist.d_dist == test_EU_dist);
#endif

			/*=========================Evaluation=================================*/
			input_argument.IO_cost++;
			/*=====================================================================*/

			cheby_candidate.push_front(id_dist);

			if (id_dist.d_dist <= radius) {
				//cout << "   range id: " << id_dist.trajectory_id << ", distance euc: " << id_dist.d_dist << endl;
				queue.push(id_dist);
				cheby_candidate.pop_front();
			}
		}
		//it.deleteLeafNodeMemory();
	}
	//cout << "Range search queue.size: " << queue.size() << endl;
	cheby_candidate.sort([](const ID_DIST& first, const  ID_DIST& second) {return first.d_dist < second.d_dist; });//small to big

	//*************180625 improvement: in case queue.size< K.*********************************
	if (queue.size() < input_argument.K && !cheby_candidate.empty()) {
		int difference = input_argument.K - queue.size();
		while (difference > 0) {
			queue.push(cheby_candidate.front());
			cheby_candidate.pop_front();
			difference--;
		}
	}
	//********************************************************************************************
	//assert(queue.size() >= input_argument.K);
	/*while (!queue.empty()) {
		cout << "id: " << queue.top().trajectory_id << ", dist: " << queue.top().d_dist << endl;
		queue.pop();
	}*/

	delete[] original_time_series;
	original_time_series = nullptr;
	original_time_series_vector.clear();
	original_time_series_vector.shrink_to_fit();
	query_time_series_vector.clear();
	query_time_series_vector.shrink_to_fit();
	//priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >().swap(queue);
	deleteCHEBYSHEV(chebyshev_query);
}


//***************************************************************
// Method:KNNSearch
// Qualifier:  
// Input:
// Output: 
// date:
// author:
//***************************************************************
TEMPLATE//KNN Search
void CHEBYSHEV_QUAL::KNNSearch(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree, priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >& queue_result) {
	input_argument.sum_distance_euc = 0.0;
	input_argument.navigate_index_time = 0.0;// navigate time
	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
	input_argument.distance_euc_time = 0.0;// distance euclidean time
	cout << "KNN Chebyshev Search: " << endl;

	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[3]);//for total KNN time

	assert(input_argument.K <= input_argument.point_number);
	int leaf_id = NULL;
	double distance_euc = NULL;
	double distance_max = 0.0;
	double chebyshev_distance = NULL;
	ID_DIST id_dist;

	//priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement > queue_test;
	priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement > queue;//big to small
	//priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement > queue_result;//small to big
	DataType* original_time_series = new DataType[input_argument.time_series_length];

	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[1]);//for approximate query time
	CHEBYSHEV chebyshev_query;
	initialCHEBYSHEV(input_argument, chebyshev_query);
	getCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);
	APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[1], input_argument.approximation_query_time);
	//cout << "approximation query time : " << input_argument.approximation_query_time << " us" << endl;

	//cout << "KNN for distance_chebyshev" << endl;
	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[2]);//for rest part of KNN time

	double temp_navigate_time = 0;
	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[4]);//time for navigate index nodes

	RTREE::Iterator it;//KNN for distance_chebyshev
	for (RTree.GetFirst(it); !RTree.IsNull(it); RTree.GetNext(it)) {
		leaf_id = RTree.GetAt(it);
		id_dist.trajectory_id = leaf_id;

		double distance_chebyshev_time = 0;
		APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[5]);//time for chebyshve distance

		getChebyshevDistance(input_argument, chebyshev_query, chebyshev_array[leaf_id], id_dist.d_dist);
		//cout << "leaf id: " << id_dist.trajectory_id << ", chebyshev dist: " << id_dist.d_dist << endl;

		APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[5], distance_chebyshev_time);
		input_argument.distance_lowbound_time += distance_chebyshev_time;

		if (queue.size() < input_argument.K) {
			queue.push(id_dist);
			continue;
		}

		if (queue.top().d_dist > id_dist.d_dist) {
			queue.pop();
			queue.push(id_dist);
		}
	}

	APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[4], temp_navigate_time);
	input_argument.navigate_index_time += temp_navigate_time;
	//cout << "navigate time : " << input_argument.navigate_index_time << " us" << endl;

	/*for (int array_id = 0; array_id < input_argument.point_number; array_id++) {
		id_dist.trajectory_id = array_id;
		getChebyshevDistance(input_argument, chebyshev_query, chebyshev_array[array_id], id_dist.d_dist);
		cout << "chebyshev id: " << id_dist.trajectory_id << ", chebyshev dist: " << id_dist.d_dist << endl;
		if (queue_test.size() < input_argument.K) {
			queue_test.push(id_dist);
			continue;
		}
		if (queue_test.top().d_dist > id_dist.d_dist) {
			queue_test.pop();
			queue_test.push(id_dist);
		}
	}

	while (!queue.empty()) {
		cout << "id: " << queue.top().trajectory_id << ", dist: " << queue.top().d_dist << endl;
		cout << "id: " << queue_test.top().trajectory_id << ", dist: " << queue_test.top().d_dist << endl;
		assert(queue.top().trajectory_id== queue_test.top().trajectory_id);
		queue.pop();
		queue_test.pop();
	}*/

	assert(queue.size() == input_argument.K);
	//cout << "get max_distance: \n";
	while (!queue.empty()) {//according distance_EUC to get max distance in queue.
		APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);

		double distance_euclidean_time = 0;
		APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[6]);//time for euclidean distance

		distance_euc = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.time_series_length, original_time_series, input_argument.time_series_length);

		APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[6], distance_euclidean_time);
		input_argument.distance_euc_time += distance_euclidean_time;
		input_argument.sum_distance_euc++;
		//cout <<"distance euc: " <<distance_euc << endl;
		//cout<<queue.top().d_dist<<endl;
		//cout << "id: " << queue.top().trajectory_id << ", dist: " << queue.top().d_dist << endl;
		distance_max = distance_max < distance_euc ? distance_euc : distance_max;
		queue.pop();
	}
	assert(queue.empty());
	//cout << "maximum euclidean distance: " << distance_max << endl;

	//************rangeSearch
	rangeSearch(input_argument, query_time_series, chebyshev_share, chebyshev_array, RTree, distance_max, queue);
	//assert(queue.size() >= input_argument.K);

	//KNN for distance_EUC
	//cout << "KNN for distance_EUC" << endl;
	while (!queue.empty()) {
		id_dist.trajectory_id = queue.top().trajectory_id;
		id_dist.d_dist = queue.top().d_dist;
		/*APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);
		id_dist.d_dist = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.time_series_length, original_time_series, input_argument.time_series_length);*/
		queue.pop();
		queue_result.push(id_dist);
	}

	APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[2], input_argument.knn_rest_part_time);
	//cout << "rest part of KNN time : " << input_argument.knn_rest_part_time << " us" << endl;

	APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[3], input_argument.knn_total_time);
	cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;

	input_argument.pruning_power = input_argument.sum_distance_euc / double(input_argument.point_number);
	cout << "pruning power: " << input_argument.pruning_power << endl;

	input_argument.IO_cost = input_argument.sum_distance_euc;
	cout << "I/O cost: " << input_argument.IO_cost << endl;

	/*TOOL::writeSingleResult(input_argument,input_argument.write_file_name[0], input_argument.pruning_power);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[1], input_argument.IO_cost);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[3], input_argument.navigate_index_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[4], input_argument.distance_lowbound_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[5], input_argument.distance_euc_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[6], input_argument.knn_total_time);*/

	cout << "show result:     queue_result: " << queue_result.size() << endl;
	/*int count = 0;
	while (!queue_result.empty() && count < input_argument.K) {
		cout << "id: " << queue_result.top().trajectory_id << ", dist: " << queue_result.top().d_dist << endl;
		queue_result.pop();
		count++;
	}*/

	//191118  when result number < K, KNN search failed
	if (queue_result.size() < input_argument.K) {
		cout << "???????????????????????????????????? Chebyshev KNN failed ???????????????????????????????????" << endl;
	}

	TOOL::printInputArgument(input_argument);
	delete[] original_time_series;
	original_time_series = nullptr;
	priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >().swap(queue);
	//priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >().swap(queue_result);
	deleteCHEBYSHEV(chebyshev_query);
}

TEMPLATE//KNN Multi to Single Search
std::multiset<pair<double, int>> CHEBYSHEV_QUAL::KNNCHEBYMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree) {
	assert(input_argument.time_series_length != INF && input_argument.time_series_length != INF && input_argument.K <= input_argument.point_number);
#ifdef _DEBUG
	input_argument.navigate_index_time = 0.0;// navigate time
	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
	input_argument.distance_euc_time = 0.0;// distance euclidean time
	cout << "KNN Chebyshev Search: " << endl;
#endif

	/*-------------------------Evaluation ---------------------------------*/
	input_argument.knn_total_time = 0.0;// KNN total time
	input_argument.IO_cost = 0.0;
	TOOL::recordStartTime(TOOL::time_record[3]);// whole time
	TOOL::recordStartTime(TOOL::time_record[2]);// KNN total time
	/*--------------------------------------------------------------------*/

	int leaf_id = NULL;
	double distance_euc = NULL;
	double distance_max = 0.0;
	double chebyshev_distance = NULL;
	ID_DIST id_dist;

	//priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement > queue_test;
	priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement > queue;//big to small
	priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement > queue_result;//small to big
	std::multiset<pair<double, int>> result;//191204 result
	DataType* original_time_series = new DataType[input_argument.point_multi_single_length];

	CHEBYSHEV chebyshev_query;
	initialCHEBYSHEV(input_argument.time_series_length, input_argument.point_multi_single_dimension, chebyshev_query);
	//getCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);
	getMultiCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);

	RTREE::Iterator it;//KNN for distance_chebyshev
	for (RTree.GetFirst(it); !RTree.IsNull(it); RTree.GetNext(it)) {
		leaf_id = RTree.GetAt(it);
		id_dist.trajectory_id = leaf_id;

		//getChebyshevDistance(input_argument, chebyshev_query, chebyshev_array[leaf_id], id_dist.d_dist);
		getChebyshevDistance(input_argument.point_multi_single_dimension, chebyshev_query, chebyshev_array[leaf_id], id_dist.d_dist);
		//cout << "leaf id: " << id_dist.trajectory_id << ", chebyshev dist: " << id_dist.d_dist << endl;

		if (queue.size() < input_argument.K) {
			queue.push(id_dist);
			continue;
		}

		if (queue.top().d_dist > id_dist.d_dist) {
			queue.pop();
			queue.push(id_dist);
		}
	}

	//cout << "navigate time : " << input_argument.navigate_index_time << " us" << endl;

	/*for (int array_id = 0; array_id < input_argument.point_number; array_id++) {
		id_dist.trajectory_id = array_id;
		getChebyshevDistance(input_argument, chebyshev_query, chebyshev_array[array_id], id_dist.d_dist);
		cout << "chebyshev id: " << id_dist.trajectory_id << ", chebyshev dist: " << id_dist.d_dist << endl;
		if (queue_test.size() < input_argument.K) {
			queue_test.push(id_dist);
			continue;
		}
		if (queue_test.top().d_dist > id_dist.d_dist) {
			queue_test.pop();
			queue_test.push(id_dist);
		}
	}

	while (!queue.empty()) {
		cout << "id: " << queue.top().trajectory_id << ", dist: " << queue.top().d_dist << endl;
		cout << "id: " << queue_test.top().trajectory_id << ", dist: " << queue_test.top().d_dist << endl;
		assert(queue.top().trajectory_id== queue_test.top().trajectory_id);
		queue.pop();
		queue_test.pop();
	}*/
#ifdef _DEBUG
	assert(queue.size() == input_argument.K);
#endif
	//cout << "get max_distance: \n";
	while (!queue.empty()) {//according distance_EUC to get max distance in queue.
		//APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);

		if (input_argument.read_multiple_file_name) {
			TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);
		}
		else {
			//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id));
			TOOL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);
		}
		TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);

		distance_euc = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);

		/*-----------------Evaluation---------------------------------*/
		input_argument.IO_cost++;
		/*----------------------------------------------------------*/
		distance_max = distance_max < distance_euc ? distance_euc : distance_max;
		queue.pop();
	}
	assert(queue.empty());
	//cout << "maximum euclidean distance: " << distance_max << endl;

	//========================================Range  Search===============================================================
	rangeSearchMulti(input_argument, query_time_series, chebyshev_share, chebyshev_array, RTree, distance_max, queue);
	/*====================================================================================================================*/
	//assert(queue.size() >= input_argument.K);

	//KNN for distance_EUC
	//cout << "KNN for distance_EUC" << endl;
	while (!queue.empty()) {
		id_dist.trajectory_id = queue.top().trajectory_id;
		id_dist.d_dist = queue.top().d_dist;
		/*APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);
		id_dist.d_dist = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.time_series_length, original_time_series, input_argument.time_series_length);*/
		result.emplace(make_pair(id_dist.d_dist, id_dist.trajectory_id));
		queue.pop();
		queue_result.push(id_dist);
	}
	/*-------------------------Evaluation ---------------------------------*/
	input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[3]);
	assert(input_argument.IO_cost != INF && input_argument.IO_cost != 0 && input_argument.point_number != INF);
	input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
	/*---------------------------------------------------------------------*/

	if (input_argument.pruning_power > 1) {
		input_argument.pruning_power = 10;
	}

#ifdef _DEBUG

	cout << "============================    Chebyshev KNN ruslt   ======================================\n";

	cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;

	cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;

	cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;

	cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;

	assert(input_argument.pruning_power != INF);

	cout << "pruning power: " << input_argument.pruning_power << endl;


	cout << "I/O cost: " << input_argument.IO_cost << endl;

	/*TOOL::writeSingleResult(input_argument,input_argument.write_file_name[0], input_argument.pruning_power);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[1], input_argument.IO_cost);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[3], input_argument.navigate_index_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[4], input_argument.distance_lowbound_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[5], input_argument.distance_euc_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[6], input_argument.knn_total_time);*/

	cout << "show result:     queue_result: " << queue_result.size() << endl;
	int count = 0;
	while (!queue_result.empty() && count < input_argument.K) {
		//cout << "id: " << queue_result.top().trajectory_id << ", dist: " << queue_result.top().d_dist << endl;
		queue_result.pop();
		count++;
	}
	//TOOL::printInputArgument(input_argument);
#endif

	delete[] original_time_series;
	original_time_series = nullptr;
	priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >().swap(queue);
	//priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >().swap(queue_result);
	//deleteCHEBYSHEV(chebyshev_query);

	return result;
}

//***************************************************************
// Method:KNNCHEBYMulti
// Qualifier:  //for priority queue, small to big, top is smallest
// Input:add data_source
// Output: 
// date:210603
// author:
//***************************************************************
TEMPLATE
template<typename T>
std::multiset<pair<double, int>> CHEBYSHEV_QUAL::KNNCHEBYMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const typename TOOL::DATA_SOURCE& const data_source, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, T& const RTree) {
	assert(input_argument.time_series_length != INF && input_argument.time_series_length != INF && input_argument.K <= input_argument.point_number);
#ifdef _DEBUG
	input_argument.navigate_index_time = 0.0;// navigate time
	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
	input_argument.distance_euc_time = 0.0;// distance euclidean time
	assert(RTree.m_root->m_level > 0);
	cout << "KNN Chebyshev Search: " << endl;
#endif

	/*-------------------------Evaluation ---------------------------------*/
	input_argument.knn_total_time = 0.0;// KNN total time
	input_argument.knn_total_time_has_IO = 0.0;// KNN total time 210606
	input_argument.IO_cost = 0.0;
	TOOL::recordStartTime(TOOL::time_record[3]);// whole time
	TOOL::recordStartTime(TOOL::time_record[16]);// whole time has IO 210606
	TOOL::recordStartTime(TOOL::time_record[2]);// KNN total time
	TOOL::recordStartTime(TOOL::time_record[4]);// KNN total time has IO 210606
	/*--------------------------------------------------------------------*/

	int leaf_id = NULL;
	double distance_euc = NULL;
	double distance_max = 0.0;
	double chebyshev_distance = NULL;
	ID_DIST id_dist;

	//priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement > queue_test;
	priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement > queue;//big to small
	priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement > queue_result;//small to big
	std::multiset<pair<double, int>> result;//191204 result
	DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
	vector<DataType> original_time_series_vector;
	vector<DataType> query_time_series_vector(query_time_series, query_time_series + input_argument.point_multi_single_length);

#ifdef _DEBUG
	for (int array_id = 0; array_id < input_argument.point_multi_single_length; array_id++) {
		assert(query_time_series_vector[array_id] == query_time_series[array_id]);
	}
#endif

	CHEBYSHEV chebyshev_query;
	initialCHEBYSHEV(input_argument.time_series_length, input_argument.point_multi_single_dimension, chebyshev_query);
	getCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);//191224
	//getMultiCHEBYSHEV(input_argument, query_time_series, chebyshev_share, chebyshev_query);

	typename T::Iterator it;//KNN for distance_chebyshev
	for (RTree.GetFirst(it); !RTree.IsNull(it); RTree.GetNext(it)) {
		leaf_id = RTree.GetAt(it);
		id_dist.trajectory_id = leaf_id;

		//getChebyshevDistance(input_argument, chebyshev_query, chebyshev_array[leaf_id], id_dist.d_dist);
		getChebyshevDistance(input_argument.point_multi_single_dimension, chebyshev_query, chebyshev_array[leaf_id], id_dist.d_dist);
		//cout << "leaf id: " << id_dist.trajectory_id << ", chebyshev dist: " << id_dist.d_dist << endl;

		if (queue.size() < input_argument.K) {
			queue.push(id_dist);
			continue;
		}

		if (queue.top().d_dist > id_dist.d_dist) {
			queue.pop();
			queue.push(id_dist);
		}
	}

	//cout << "navigate time : " << input_argument.navigate_index_time << " us" << endl;

	/*for (int array_id = 0; array_id < input_argument.point_number; array_id++) {
		id_dist.trajectory_id = array_id;
		getChebyshevDistance(input_argument, chebyshev_query, chebyshev_array[array_id], id_dist.d_dist);
		cout << "chebyshev id: " << id_dist.trajectory_id << ", chebyshev dist: " << id_dist.d_dist << endl;
		if (queue_test.size() < input_argument.K) {
			queue_test.push(id_dist);
			continue;
		}
		if (queue_test.top().d_dist > id_dist.d_dist) {
			queue_test.pop();
			queue_test.push(id_dist);
		}
	}

	while (!queue.empty()) {
		cout << "id: " << queue.top().trajectory_id << ", dist: " << queue.top().d_dist << endl;
		cout << "id: " << queue_test.top().trajectory_id << ", dist: " << queue_test.top().d_dist << endl;
		assert(queue.top().trajectory_id== queue_test.top().trajectory_id);
		queue.pop();
		queue_test.pop();
	}*/
#ifdef _DEBUG
	assert(queue.size() == input_argument.K);
#endif
	//cout << "get max_distance: \n";
	while (!queue.empty()) {//according distance_EUC to get max distance in queue.
		//APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);

		/*-----------------------------------------------------Original Time Series---------------------------------------------------------------------------------------*/
		input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
		input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[3]);
		//already normalized
		TOOL::read_normalized_multi_time_series(data_source, queue.top().trajectory_id, original_time_series_vector);
		TOOL::recordStartTime(TOOL::time_record[3]);// whole time
		TOOL::recordStartTime(TOOL::time_record[2]);// KNN total time

#ifdef _DEBUG
		assert(original_time_series_vector.size() == input_argument.time_series_length);
		copy_n(original_time_series_vector.begin(), original_time_series_vector.size(), original_time_series);
#endif
		/*--------------------------------------------------------------------------------------------------------------------------------------------*/
		//if (input_argument.read_multiple_file_name) {
		//	TOOL::getMultiFoldToSingleByID(input_argument.read_multiple_file_name, input_argument.arity_d, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);
		//}
		//else {
		//	//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id));
		//	TOOL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);
		//}
		//TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);

		distance_euc = TOOL::distanceEUC(query_time_series_vector, original_time_series_vector);

#ifdef _DEBUG
		double test_distance_euc = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.point_multi_single_length, original_time_series, input_argument.point_multi_single_length);
		assert(distance_euc == test_distance_euc);
#endif

		/*-----------------Evaluation---------------------------------*/
		input_argument.IO_cost++;
		/*----------------------------------------------------------*/
		distance_max = distance_max < distance_euc ? distance_euc : distance_max;
		queue.pop();
	}
	assert(queue.empty());
	//cout << "maximum euclidean distance: " << distance_max << endl;

	//========================================Range  Search===============================================================
	rangeSearchMulti(input_argument, data_source, query_time_series, chebyshev_share, chebyshev_array, RTree, distance_max, queue);
	/*====================================================================================================================*/
	//assert(queue.size() >= input_argument.K);

	//KNN for distance_EUC
	//cout << "KNN for distance_EUC" << endl;
	while (!queue.empty()) {
		id_dist.trajectory_id = queue.top().trajectory_id;
		id_dist.d_dist = queue.top().d_dist;
		/*APCA_KNN_QUAL::getFileStreamByID(input_argument.read_file_name, input_argument.time_series_length, queue.top().trajectory_id, original_time_series);
		id_dist.d_dist = typename APCA_KNN_QUAL::distanceEUC(query_time_series, input_argument.time_series_length, original_time_series, input_argument.time_series_length);*/
		result.emplace(make_pair(id_dist.d_dist, id_dist.trajectory_id));
		queue.pop();
		queue_result.push(id_dist);
	}
	/*-------------------------Evaluation ---------------------------------*/
	input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
	input_argument.knn_total_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[4]);
	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[3]);
	input_argument.whole_run_time_has_IO += TOOL::recordFinishTime(TOOL::time_record[16]);
	assert(input_argument.IO_cost != INF && input_argument.IO_cost != 0 && input_argument.point_number != INF);
	input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
	/*---------------------------------------------------------------------*/

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 200906 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	if (input_argument.pruning_power > 1) {
		//input_argument.pruning_power = 30;
		//assert(0);
		//cout <<"!!!!!!!!!!!!!!!!!" <<input_argument.pruning_power << endl;
	}
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

#ifdef _DEBUG

	cout << "============================    Chebyshev KNN ruslt   ======================================\n";

	cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;

	cout << "Total KNN time has IO : " << input_argument.knn_total_time_has_IO << " us" << endl;

	cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;

	cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;

	cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;

	assert(input_argument.pruning_power != INF);

	cout << "pruning power: " << input_argument.pruning_power << endl;


	cout << "I/O cost: " << input_argument.IO_cost << endl;

	/*TOOL::writeSingleResult(input_argument,input_argument.write_file_name[0], input_argument.pruning_power);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[1], input_argument.IO_cost);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[3], input_argument.navigate_index_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[4], input_argument.distance_lowbound_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[5], input_argument.distance_euc_time);
	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[6], input_argument.knn_total_time);*/

	cout << "show result:     queue_result: " << queue_result.size() << endl;
	int count = 0;
	while (!queue_result.empty() && count < input_argument.K) {
		//cout << "id: " << queue_result.top().trajectory_id << ", dist: " << queue_result.top().d_dist << endl;
		queue_result.pop();
		count++;
	}
	//TOOL::printInputArgument(input_argument);
#endif

	delete[] original_time_series;
	original_time_series = nullptr;
	original_time_series_vector.clear();
	original_time_series_vector.shrink_to_fit();
	query_time_series_vector.clear();
	query_time_series_vector.shrink_to_fit();
	priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >().swap(queue);
	//priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >().swap(queue_result);
	//deleteCHEBYSHEV(chebyshev_query);
	return result;
}

TEMPLATE//Write approximation result to text file
void CHEBYSHEV_QUAL::writeApproximationResult(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const result_time_series, const double& deviation_sum, const double& deviation_max) {
	double count = 0;
	time_t now = time(0);// system time
	char dt[26];
	ctime_s(dt, sizeof dt, &now);// from now to string
	ofstream outfile(*input_argument.write_file_name + ".txt", ios::app);
	outfile << endl << dt << "read file:" << input_argument.read_file_name << "; write file: " << input_argument.write_file_name << endl;
	outfile << " time series length = " << input_argument.time_series_length << ", degree_m = " << input_argument.degree_m << ", deviation_sum = " << deviation_sum << ", deviation_max = " << deviation_max << endl;// << ", point number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;

	outfile << " result time series: " << endl;
	for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		outfile << result_time_series[array_id] << ", ";
	}

	outfile << endl;
	/*outfile << "  count apca point = " << g_n_account_apca_point << " times, p= " << g_n_account_apca_point / input_argument.point_number << endl;
	outfile << "  APCA Memory = " << memory_account[0] + memory_account[1] + memory_account[2] + memory_account[3] << "  RTree Memeory = " << memory_account[4] + memory_account[5] + memory_account[6] + memory_account[7] + memory_account[8] << "  KNN Memory = " << memory_account[9] << endl;
	for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
		if (it->d_dist == q_base_queue.top().d_dist) count++;
		outfile << "id: " << it->original_time_series_id << " , " << q_base_queue.top().original_time_series_id << "; dist: " << it->d_dist << " , " << q_base_queue.top().d_dist << endl;
		q_base_queue.pop();
	}*/

	outfile.close();
}

TEMPLATE
template<typename T>
void CHEBYSHEV_QUAL::writeSingleResult(const string& write_file_name, T& result) {
	ofstream outfile(write_file_name + ".txt", ios::app);
	assert(outfile.is_open());
	outfile << result << endl;
	outfile.close();
}

TEMPLATE
void CHEBYSHEV_QUAL::getSeveralRecronstructionError() {
}
#endif

