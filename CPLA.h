#pragma once
#ifndef CPLA_H
#define CPLA_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "CAPCA.h"
#include "CAPCA_KNN.h"

TEMPLATE
class CPLA : public APCA_QUAL, public APCA_KNN_QUAL, virtual public TOOL
{
public:

	struct TIME;
	TIME time_record[20];
	struct INPUT_ARGUMENT;
	struct PLA;
	struct PLA_C;//191022 no pointer, a&b are variable, with rightn endpoint, widht, sum value
	struct POINT;
	struct PLA_NODE_PAIR;
	struct priorityIncrement;
	struct ORIGINAL_TIME_SERIES_PAIR;
	struct priorityDistanceEUC;

	bool compare(ORIGINAL_TIME_SERIES_PAIR& first, ORIGINAL_TIME_SERIES_PAIR& second);

	INPUT_ARGUMENT input_argument;
	struct TOOL::INPUT_ARGUMENT tool_input_argument;
	struct TOOL::OUTPUT_ARGUMENT output_argument;//181214

	CPLA() {};
	CPLA(const int& n, const int& N);
	CPLA(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const string& read_file_name, string*& const write_file_name);
	void initialPLA(PLA& pla, const int& N);
	void deletePLA(PLA& pla);
	void initialPLAArray(const INPUT_ARGUMENT& input_argument, PLA*& const pla_array);
	void initialPLAArray(const typename TOOL::INPUT_ARGUMENT& input_argument, PLA*& const pla_array);
	void initialPLAArray(const int& point_number, const int& point_dimension, PLA*& const pla_array);
	void deletePLAArray(const INPUT_ARGUMENT& input_argument, PLA*& const pla_array);
	void deletePLAArray(const typename TOOL::INPUT_ARGUMENT& input_argument, PLA*& const pla_array);

	void getPLA(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, const PLA& pla);
	void getPLAVector(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, vector<PLA>& const pla);
	void getPLAVectorNormal(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, vector<PLA>& const pla);//190604
	template<typename T>
	void getPLALinkedListNormal(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, DoublyLinkedList<T>& const pla);//191104 11:12
	void getPLAVectorPreMemory(const INPUT_ARGUMENT& input_argument, const DataType* original_time_series, vector<PLA>& const pla);//190611
	void getPLA(const typename TOOL::INPUT_ARGUMENT& input_argument, const DataType* original_time_series, const PLA& pla);
	void getPLA(const int& time_series_length, const int& segment_number, const DataType* original_time_series, const PLA& pla);
	PLA& getPLAMBR(const PLA& pla, PLA& pla_MBR);
	double& getPLADistance(const INPUT_ARGUMENT& input_argument, const PLA& pla_array, const PLA& pla_array_qeury, double& pla_distance);
	double& getPLADistance(const int& time_series_length, const PLA& pla_array, const PLA& pla_array_qeury, double& pla_distance);
	inline double getPointDistanceSquare(const DataType& point_a_x, const DataType& point_a_y, const DataType& point_b_x, const DataType& point_b_y);
	double getNearestDistance(const DataType& point_a_x, const DataType& point_a_y, const DataType& point_b_x, const DataType& point_b_y, const DataType& point_q_x, const DataType& point_q_y);
	double& getPLAMBRSegmentDistance(const typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const int& segment_id, const PLA& pla_array_qeury, double& pla_MBR_segment_distance);
	double& getPLAMBRSegmentDistance(INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const int& segment_id, const PLA& pla_array_qeury, double& pla_MBR_segment_distance);
	double& getPLAMBRSegmentDistanceBase(const RTREE::Rect& MBR, const int& segment_id, const PLA& pla_array_qeury, double& pla_MBR_segment_distance);
	double& getPLAMBRDistance(typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const PLA& pla_array_qeury, double& pla_MBR_distance);
	double& getPLAMBRDistance(INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const PLA& pla_array_qeury, double& pla_MBR_distance);
	void approximateOriginalFunctionPLA(const INPUT_ARGUMENT& input_argument, const PLA& pla, DataType*& const approximate_array);//wrong
	void approximateOriginalFunctionPLA1(const INPUT_ARGUMENT& input_argument, const PLA& pla, DataType*& const approximate_array);//right: because PLA paper translate every point to begin.
	void approximateOriginalFunctionPLA1(const typename TOOL::INPUT_ARGUMENT& const input_argument, const PLA& pla, DataType*& const approximate_array);//right: because PLA paper translate every point to begin.
	//191206 get sum deviation
	double get_pla_sum_deviation(const typename  TOOL::INPUT_ARGUMENT& const input_argument, const vector<DataType>& const original_time_series, const PLA& const pla);
	//210327
	template<typename T, typename Y, typename U, typename T1>
	long double get_pla_sum_deviation(const T& const input_argument, const vector<Y>& const original_time_series, const T1& const pla, U& const result_collection);
	double& getReconstructionErrorPLA(const INPUT_ARGUMENT& input_argument, DataType*& const original_array, const PLA& pla, double& const deviation_sum, double& const deviation_max);// deviation_sum deviation_max
	double& getReconstructionErrorPLA(const INPUT_ARGUMENT& const input_argument, double*& const original_array, double& const deviation_sum, double& const deviation_max);// 181212 deviation_sum deviation_max
	double& getReconstructionErrorPLAVector(const INPUT_ARGUMENT& const input_argument, double*& const original_array, double& const deviation_sum, double& const deviation_max);// 190501 Vector to instead array
	void getDeviationIterationPLA(DataType*& const original_time_series, const int& const array_length, string*& const write_file_name, const int& segment_begin, const int& segment_end, const int& segment_interval);

	/*====================================================================Build RTree=======================================================================================*/
	RTREE& buidRTreeIndex(INPUT_ARGUMENT& const  input_argument, RTREE& PLARTree, PLA* PLA_array_accumulate, const string& file_name);
	/*=============================================================================================================================================================================*/
	/*====================================================================KNN====================================================================================================*/
	bool PLAKNNSearch(INPUT_ARGUMENT& const input_argument, const DataType* g_query_time_series, const RTREE& PLARTree, PLA* PLA_array_accumulate, const int& K, const string& file_name, list<ORIGINAL_TIME_SERIES_PAIR>& result);
	// 180919
	bool PLAKNNMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const DataType* g_query_time_series, const RTREE& PLARTree, PLA* PLA_array_accumulate, const int& K);
	void SimpleBaseKNNSearch(const INPUT_ARGUMENT& input_argument, const DataType* g_query_time_series, const int& K, const string& file_name, priority_queue<ORIGINAL_TIME_SERIES_PAIR, vector<ORIGINAL_TIME_SERIES_PAIR>, priorityDistanceEUC >& q_base_queue);
	/*=================================================================================================================================================================*/


	/*====================================================================Get a & b of one segment====================================================================================================*/
	//190619
	double& getAAndBByPLA(DataType*& const original_time_series, PLA& const temp_coefficient);//190619 return a;
	double& computeParallelogram(DataType*& const original_time_series, PLA& const temp_coefficient);//190619
	//191016
	template<typename T>
	void getAAndBByPLASegment(const vector<DataType>& const original_time_series_vector, T& const temp_coefficient);//191016 get a&b for one segment

	//191215
	template<typename T>
	T& get_a_b_by_endpoint(const vector<DataType>& const original_time_series_vector, T& const temp_coefficient);//191215 get a&b for one segment

	//191021
	template<typename T>
	T& getAAndBByAccumulation(const vector<DataType>& const original_time_series_vector, const T& const last_segment, T& const temp_coefficient);//191021 the length of segment is added one by one 
	/*=================================================================================================================================================================*/

   //191016 get segment difference with original time series
	template<typename T>
	long double getSegmentMaxDifference(const vector<DataType>& const original_time_series_vector, const T& const segment_pla);

	//191018 get segment sum difference
	template<typename T>
	double getSegmentSumDifferenceSquare(const vector<DataType>& const original_time_series_vector, const T& const segment_pla, vector<DataType>& const rescontruct_time_series);

	//191018 get segment sum difference
	template<typename T>
	double getSegmentSumDifferenceSquare(const vector<DataType>& const original_time_series_vector, const T& const segment_pla);

	//191018 get segment sum difference,Except for first element. begin() + 1
	template<typename T>
	double getPLASumDeviation(const vector<DataType>& const original_time_series_vector, const vector<T>& const segment_vector);

	//191206 Linked List get segment sum difference,Except for first element. begin() + 1
	template<typename T>
	double getPLASumDeviation(const vector<DataType>& const original_time_series_vector, const DoublyLinkedList<T>& const segment_list);

	//210327 Linked List get segment sum difference,Except for first element. begin() + 1
	template<typename T, typename Y, typename U>
	long double getPLASumDeviation(const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const segment_list, U& const result_collection);

	/*----------------------Time-----------------------------------------*/
	void recordStartTime(TIME& time);
	double recordFinishTime(TIME& time, double& whole_first_run_time);
	/*------------------------------------------------------------*/
	void writeResult(const INPUT_ARGUMENT& input_argument, const string& write_file_name, list<ORIGINAL_TIME_SERIES_PAIR>& result, priority_queue<ORIGINAL_TIME_SERIES_PAIR, vector<ORIGINAL_TIME_SERIES_PAIR>, priorityDistanceEUC >& q_base_queue);
	//friend class APCA_KNN_QUAL;
private:
};

TEMPLATE
struct PLA_QUAL::TIME {
	_LARGE_INTEGER time_start;
	_LARGE_INTEGER time_over;   //finish time
	double dqFrequency = NULL;      //timer frequency
	double run_time = NULL;

	~TIME() {
		time_start.QuadPart = NULL;
		time_over.QuadPart = NULL;   //finish time
		dqFrequency = NULL;      //timer frequency
		run_time = NULL;
	}
};

TEMPLATE
struct PLA_QUAL::INPUT_ARGUMENT {
	int time_series_length = INF;//n
	int point_dimension = INF;//N
	int point_number = INF;
	int remainder = INF;
	int segment_length_first = INF;//l=n/m+1
	int segment_length_second = INF;//l=n/m

	int rtree_max_nodes = INF;
	int K = INF;
	string read_file_name;
	string* write_file_name = nullptr;

	//for Chebyshev
	int degree_m = INF;//point_dimension for APCA & PLA
	int arity_d = INF;// dimension of every point, for multi-dimensional case

	//for KNN
	double pruning_power = INF;
	double sum_distance_euc = INF;

	//for time
	double build_rtree_time = INF;
	double approximation_query_time = INF;
	double knn_rest_part_time = INF;
	double knn_total_time = INF;

	//    three part time
	double navigate_index_time = INF;// navigate time
	double distance_lowbound_time = INF; // distance chebyshev, PLA, APCA time
	double distance_euc_time = INF;// distance euclidean time

	double result_accuracy = INF;// KNN result accuracy

	//I/O cost
	double IO_cost = NULL;

	~INPUT_ARGUMENT() {
		time_series_length = INF;//n
		point_dimension = INF;//N
		point_number = INF;
		remainder = INF;
		segment_length_first = INF;//l=n/m+1
		segment_length_second = INF;//l=n/m

		rtree_max_nodes = INF;
		K = INF;

		pruning_power = INF;
		sum_distance_euc = INF;

		build_rtree_time = INF;
		approximation_query_time = INF;
		knn_rest_part_time = INF;
		knn_total_time = INF;

		//    three part time
		navigate_index_time = NULL;// navigate time
		distance_lowbound_time = NULL; // distance chebyshev, PLA, APCA time
		distance_euc_time = NULL;// distance euclidean time

		IO_cost = NULL;

		result_accuracy = NULL;// KNN result accuracy

		read_file_name.clear();
		read_file_name.shrink_to_fit();

		/*write_file_name.clear();
		write_file_name.shrink_to_fit();*/
		write_file_name = nullptr;
		//read_file_name = nullptr;

		//for Chebyshev
		degree_m = NULL;
	}
};

TEMPLATE
struct PLA_QUAL::PLA {
	DataType* a = nullptr;
	DataType* b = nullptr;
	int segmentNum = NULL; //the dimension of the index(eg. MBR, imput parameter)

	/*================================================*/
	double right_endpoint = INF;//190619
	int segment_width = INF; //190619 the number of points in every segment
	double sum_value = INF; //191021
	double segment_area = INF;//190619 area of
	double segment_density = INF;//190619
	double segment_area0 = INF;//190619 area of
	double segment_density0 = INF;//190619
	double parallelogram_height = INF;//190619
	/*...............................................*/

	~PLA() {
		segmentNum = NULL;

		right_endpoint = INF;
		segment_width = INF;
		sum_value = INF; //191021
		segment_area = INF;
		segment_density = INF;
		parallelogram_height = INF;
	}
};

TEMPLATE
PLA_QUAL::CPLA(const int& n, const int& N) {
	//For PLA input_argument
	input_argument.time_series_length = n;
	input_argument.point_dimension = N;

	input_argument.remainder = int(n) % int(N);
	double integerDividend = n - input_argument.remainder;
	input_argument.segment_length_second = integerDividend / N;
	input_argument.segment_length_first = input_argument.segment_length_second + 1;

	assert(input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1

	//For TOOL input_argument
	tool_input_argument.time_series_length = n;
	tool_input_argument.point_dimension = N;

	tool_input_argument.remainder = int(n) % int(N);
	double tool_integerDividend = n - tool_input_argument.remainder;
	tool_input_argument.segment_length_second = tool_integerDividend / N;
	tool_input_argument.segment_length_first = tool_input_argument.segment_length_second + 1;

	assert(tool_input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
}

#include "CPLA.cpp"

#endif


