#pragma once
#ifndef SHARE_TOOL_H
#define SHARE_TOOL_H

#include "pch.h"
#include "GEOMETRY_TOOL.h"

//#pragma warning(disable : 4996)
TEMPLATE
class SHARE_TOOL {
public:
	const string file_address_24 = "./191202DataSet/SingleDataAddress200103/AllFiles24/DataAddress181218.txt";
	const string single_file_number_address_24 = "./191202DataSet/SingleDataAddress200103/AllFiles24/DataNumberLine200103.txt";
	const string single_file_length_address_24 = "./191202DataSet/SingleDataAddress200103/AllFiles24/DataLength200103.txt";
	const string file_address_UCR2018_512_41 = "./191202DataSet/SingleDataAddress200103/Single512UCR2018/DataAddressUCR201014.txt";
	const string single_file_number_address_UCR2018_512_41 = "./191202DataSet/SingleDataAddress200103/Single512UCR2018/DataNumberLineUCR201014.txt";
	const string single_file_length_address_UCR2018_512_41 = "./191202DataSet/SingleDataAddress200103/Single512UCR2018/DataLengthUCR201014.txt";
	const string file_address_UCR2018_1024_21 = "./191202DataSet/SingleDataAddress200103/Single1024UCR2018/DataAddress1024UCR2018201016.txt";
	const string single_file_number_address_UCR2018_1024_21 = "./191202DataSet/SingleDataAddress200103/Single1024UCR2018/DataNumberLine1024UCR2018201016.txt";
	const string single_file_length_address_UCR2018_1024_21 = "./191202DataSet/SingleDataAddress200103/Single1024UCR2018/DataLength1024UCR2018201016.txt";
	const string file_address = "./191202DataSet/SingleDataAddress200103/DataAddress200116.txt";
	const string single_file_number_address = "./191202DataSet/SingleDataAddress200103/DataNumberLine200116.txt";
	const string single_file_length_address = "./191202DataSet/SingleDataAddress200103/DataLength200116.txt";
	const string no_burst_single_file_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataAddress200130.txt";
	const string no_burst_single_file_number_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataNumberLine200130.txt";
	const string no_burst_single_file_length_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataLength200130.txt"; 
	struct DATA_SOURCE;
	struct INPUT_ARGUMENT;
	struct OUTPUT_ARGUMENT;
	struct RESULT_RECORD;
	struct Y_PROJECTION_ARGUMENT;
	struct EVALUATION_COEFFCIENTS_SPLIT_ID;
	struct EVALUATION_BOUND;
	struct TIME;
	TIME time_record[50];
	struct ID_DIST;
	struct priorityDistanceEUC;

	struct APLA_COEFFICIENT;//191121 For PLA Coefficient same with APLA class

	struct RECTANGLE;//191127 for Rtree insersion. get min&max id value of rectangle

	template<typename T>//for NULL
	void initialArray(T*& const test_array, const int& array_length);

	template<typename T>//for new memory
	void newArray(T*& const test_array, const int& array_length);

	template<typename T>//for NULL
	void deleteArray(T*& const test_array);
	//191128 clear memory of rectangle
	template<typename T>
	void delete_rectangle(T& const rectangle);

	//miss first point. file_stream >> fs_row_string;
	void getFileStreamByID(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, DataType*& const original_time_series);

	//191115 vector. miss first point. file_stream >> fs_row_string;
	template<typename T>
	void getFileStreamByID(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, vector<T>& const original_time_series);

	//Not miss first point. getline(file_stream, fs_row_string);
	void getFileStreamByID0(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, DataType*& const original_time_series);

	//Not miss first point.Vector 190702
	template<typename T>
	vector<T>& getFileStreamByID0Vector(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, vector<T>& const result_vector);

	//Not miss first point. No length .Vector 200103 
	template<typename T>
	vector<T>& getFileStreamByID0Vector(const string& file_name, const int& const original_time_series_id, vector<T>& const result_vector);

	//181016 get segment from specific row with specific length
	void getFileStreamSegment(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, const int& const column_id, DataType*& const original_time_series);

	void getMultiFoldToSingleByID(string*& const multi_file_name, const int& arity_d, const double& single_series_length, const int& const original_time_series_id, DataType*& const original_time_series);

	//191220 Original time series to instead vector
	template<typename T>
	vector<T>& getMultiFoldToSingleByID(string*& const multi_file_name, const int& arity_d, const double& single_series_length, const int& const original_time_series_id, vector<T>& const original_time_series_vector);

	//191223 vector file addresss vector 
	template<typename T>
	vector<T>& getMultiFoldToSingleByID(const vector<string>& const multi_file_name, const int& arity_d, const double& single_series_length, const int& const original_time_series_id, vector<T>& const original_time_series_vector);

	//191220 Porject multi dataset to long single dataset
	void project_multi_data_to_single_data(string*& const multi_file_name, const int& const arity_d, const double& const single_series_number, const double& const single_series_length, const string& const write_file_name);

	//191229 project long sigle dataset to multi dimension dataset
	void project_single_data_to_multi_data_single(const string& const single_file_name, const int& const arity_d, const double& const single_series_number, const double& const single_series_length, const string& const write_multi_file_address);

	//191229
	void project_single_data_to_multi_data_batch(const string& const single_file_address_file, const int& const file_number, const int& const arity_d, const double& const single_series_number, const double& const single_series_length, const string& const write_multi_file_address_file);

	//miss first point
	void getFileStreamByRow(const string& file_name, const double& g_time_series_length, const int& const row_number, DataType** original_time_series);
	//miss first point
	template<typename T>
	void getFileStreamByRow(const string& const file_name, const int& const g_time_series_length, const int& const row_number, T& const original_time_series);

	// return g_time_series_length + 1.
	template<typename T>
	void getFileStreamByRowNoMiss(const string& const file_name, const int& const g_time_series_length, const int& const row_number, T& const original_time_series);

	//miss first point
	//void getFileStreamByRow(const string& const file_name, const int& const g_time_series_length, const int& const row_number, vector<vector<double>>& const original_time_series);
	//get stream from multiple files
	template<typename T>
	void get_write_multi_files(const string& const file_name, const vector<T>& const data_id_list, const double& const g_time_series_length, const int& const row_number, const string& const write_file_name);

	void getTXTStreamSpace(const string& const file_name, const int& const array_length, DataType*& const original_time_series);

	string& getStringByID(const string& const file_name, const int& const row_num, string& const file_address);//181218
	//191115 same as above
	string getStringByID(const string& const file_name, const int& const row_num);

	//191222 get all string and read it into vector
	vector<string>& get_all_string_by_row(const string& const file_name, vector<string>& const string_vector);

	//191115 get normalized original time series from txt file list according to file id
	template<typename T>
	void get_normalized_series_by_stream(const TOOL::INPUT_ARGUMENT& const input_argument, ifstream& const file_stream, T*& const original_time_series);
	//191115 get normalized original time series from txt file list according to file id
	template<typename T>
	const vector<T>& const get_normalized_series_by_address_list(TOOL::INPUT_ARGUMENT& const input_argument, const int& const file_id, vector<T>& const normalized_series_vector) const;
	//191203 For KNN, get & write mixed dataset, convert file_total_number to 1
	void initial_mixed_dataset(const string& const write_file_name, int& const file_total_number);

	// 191201 covnert vector list to string to write into txt file
	template<typename T>
	string convert_vector_to_string(const vector<T>& const vector_value);

	void writeInputArgument(INPUT_ARGUMENT& input_argument, const string& write_file_name);//181207

	void writeKNNResult(TOOL::INPUT_ARGUMENT& input_argument);

	template<typename T>
	void writeArray(const string& write_file_name, T*& const Array, const int& array_length);

	//191126  cover content
	template<typename T>
	void writeSingleResult(const string& write_file_name, T& result);

	//191125 // with system time date. cover content
	template<typename T>
	void writeSingleResultTimeStamp(const string& write_file_name, T& result);

	//with system time date. not cover content
	template<typename T>
	void writeSingleResultNoCoverTimeStamp(const string& write_file_name, T& result);

	template<typename T, typename Y>
	void writeSingleResult(Y& input_argument, const string& write_file_name, T& result);//181207 Add input_argument

	template<typename T, typename Y>
	void writeSingleResult(const string& const write_file_name, T& array_data, const Y& const array_len);//181218 Write array to file

	template<typename T>
	void writeSingleResult(const string& const write_file_name, vector<T>& array_data);//190501 Write vector to file

	//automatilly + "txt"
	template<typename T>
	void writeResultNoCover(const string& const write_file_name, vector<T>& array_data);//190624 Write vector to file, no coverage

	//No suffix like "txt"
	template<typename T>
	void writeResultNoCoverNoSuffix(const string& const write_file_name, vector<T>& array_data);//190624 Write vector to file, no coverage

	template<typename T>
	void writeDataTXTCover(const string& const write_file_name, T& data);//190625 Write T to file, coverage

	void clearTXTFile(const string& const clear_file_name);//190625

	inline void recordStartTime(TIME& time);

	template<typename T>
	T& recordFinishTime(TIME& time, T& whole_first_run_time);

	inline long double recordFinishTime(TIME& time);// 180923

	template<typename T>//For sum deviation
	double& distanceEUC(const T& const input_argument, DataType*& const Q, DataType*& const C, double& const distance_euc);

	//191129 Euclidean distance.
	template<typename T>//For sum deviation
	double distanceEUC(const vector<T>& const time_series_vector1, const vector<T>& const time_series_vector2);

	template<typename T, typename Y>//For sum deviation 181214
	OUTPUT_ARGUMENT& getDeviation(const T& const input_argument, Y*& const original_time_series, Y*& const reconstruct_time_series, OUTPUT_ARGUMENT& const output_argument);

	template<typename T, typename Y>//For sum deviation 190501
	OUTPUT_ARGUMENT& getDeviation(const T& const input_argument, Y*& const original_time_series, vector<Y>& const reconstruct_time_series, OUTPUT_ARGUMENT& const output_argument);

	template<typename T, typename Y>//For sum deviation 200212
	OUTPUT_ARGUMENT& getDeviation(const T& const input_argument, const  vector<Y>& const original_time_series_vector, vector<Y>& const reconstruct_time_series_vector, OUTPUT_ARGUMENT& const output_argument);

	template<typename T, typename Y, typename U, typename T1>//For sum deviation & max deviation 210327
	long double getDeviation(const T& const input_argument, const vector<Y>& const original_time_series_vector, vector<Y>& const reconstruct_time_series_vector, const U& const segment_width, T1& const result_collection);

	template<typename Y>//For sum deviation 190619
	OUTPUT_ARGUMENT& getDeviation(Y*& const original_time_series, vector<Y>& const reconstruct_time_series, const int& const time_series_length, OUTPUT_ARGUMENT& const output_argument);

	template<typename T>//For sum deviation 191025
	double getDeviation(const vector<T>& const original_time_series, const vector<T>& const reconstruct_time_series);

	template<typename T, typename Y>
	void SimpleBaseKNNSearch(const T& const input_argument, DataType*& const g_query_time_series, priority_queue<Y, vector<Y>, priorityDistanceEUC >& q_base_queue);

	//191203 vector instead pointer. Normalized time series
	template<typename T, typename Y>
	multiset<pair<double, int>>& SimpleBaseKNNSearch(const T& const data_source, const vector<Y>& const query_time_series, multiset<pair<double, int>>& const knn_result_set);

	//191223 For multi dimension dataset. Vector instead pointer. Normalized time series
	template<typename T, typename Y, typename U>
	multiset<pair<U, int>>& SimpleBaseKNNSearchMulti(const T& const data_source, const vector<Y>& const query_time_series_vector, multiset<pair<U, int>>& const knn_result_set);

	template<typename T> //Get random point for the array 181113
	void getRandomPoint(T*& const original_array, const int const& array_length, const int const& scale);
	//191119 [0, max_value - 1]
	int get_random_max(const int& const max_value);
	//200108
	template<typename T>
	void get_random_vector(const int& const vector_size, const int& const scale, vector<T>& const time_series_vector);
	/*---------------------------------------------------------------------------------------------------------------------------------*/

	template<typename T> //Get mean vlaue between two value 181113
	T& getMeanValue(const T& const  vlaue_A, const T& const value_B, T& const  mean_value);

	template<typename T>
	void normalizeA_B(const DataType& left_endpoint, const DataType right_endpoint, T*& const original_array, const int& array_length, T*& const normalized_array);

	template<typename T>
	double& getAverage(T*& const original_array, const int& const array_length, double& average);
	//191115 vector
	template<typename T>
	double& getAverage(vector<T>& const original_array, const int& const array_length, double& average);

	template<typename T>
	double getAverage(T*& const original_array);//190307

	template<typename T>
	double& getVariance(T*& const original_array, const int& const array_length, double& variance);
	//191115 vector
	template<typename T>
	double& getVariance(vector<T>& const original_array, const int& const array_length, double& variance);


	template<typename T>
	void normalizeStandard(T*& const original_array, const int& const array_length, T*& const normalized_array);

	template<typename T>
	void normalizeStandard(const int& const array_length, T*& const normalized_array);
	//191115
	template<typename T>
	vector<T>& normalizeStandard(vector<T>& const normalized_array);

	template<typename T>
	void transfer_us_s(vector<T>& const time_vector);

	int transferSetToInt(const set<int>& const random_endpoint_set);//190623

	void combinationUtil(int arr[], int n, int r, int index, int data[], int i);//190623 /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	void printCombination(int arr[], int n, int r);//190623 // The main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()

	void combinationUtilJump(int arr[], int n, int r, int index, int data[], int i);//190624 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	void printCombinationJump(int arr[], int n, int r);//190624 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()

	void combinationUtilJumpVector(const vector<int>& const arr, int index, vector<int> data, int i, set<pair<double, vector<int>>>& const endpoint_collection);//190624 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	void printCombinationJumpVector(const vector<int>& const arr, const int& const r, set<pair<double, vector<int>>>& const endpoint_collection);//190624 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()

	long long NChooseK(int n, int k);//190624 Combination
	inline unsigned long long n_choose_k(const unsigned long long& n, const unsigned long long& k);//190624

	/*----------------------       Faled Cannot deal with large amount      --------------*/
	unsigned long long fact(int n);//190624
	unsigned long long nCr(int n, int r);//190624
	/*------------------------------------------------------------------------------------*/

	/*-----------------------     Cnk    -----------------------------------*/
	unsigned long long gcd(unsigned long long x, unsigned long long y);//190624
	unsigned long long N_choose_K(unsigned long long n, unsigned long long k);//190624
	/*----------------------------------------------------------------------*/

	void printOSPageSize();//for Windows system page_size

	template<typename T>//print argument like: n, N, K, max_node
	void printInputArgument(const T& const input_argument);

	template<typename T>
	void printArray(T*& const test_array, const int& array_length);

	template<typename T>
	void printArray(vector<T>& const test_array, const int& array_length);//190501
	//200108
	template<typename T>
	void print_vector(vector<T>& const time_series_vector);//190501

	//210122 print map
	template<typename T>
	void print_map(const T& const map_result);

	//200224 Print split id coeffifents
	// local sum deviation, local shift, local time. global sum deviation global time
	template<typename T>
	inline void print_split_coefficients(const vector<T>& const local_total_split_id_sum_deviation, const vector<T>& const local_total_split_id_shift, const vector<T>& const local_total_split_id_time, const vector<T>& const global_total_approximation_sum_deviation, const vector<T>& const global_total_approximation_time);
	//200228 Print split id coeffifents
	//local sum deviation, local shift, local time. global sum deviation global time, global prune power
	template<typename T>
	inline void print_split_coefficients(const vector<T>& const local_total_split_id_sum_deviation, const vector<T>& const local_total_split_id_shift, const vector<T>& const local_total_split_id_time, const vector<T>& const global_total_approximation_sum_deviation, const vector<T>& const global_total_approximation_time, const vector<T>& const global_total_knn_prune_power);

	//200326 initial_N coefficients
	template<typename T>
	inline void print_initial_N_coefficients(const vector<T>& const  approximation_initial_N_vector, const vector<T>& const total_initial_N_prune_power_vector, const vector<T>& const  total_initial_N_sum_deviation_vector, const vector<T>& const total_initial_N_run_time_vector, const vector<T>& const total_initial_N_approximation_time_vector, const vector<T>& const total_initial_N_knn_time_vector);

	//200331 initial_N coefficients sort by N
	template<typename T, typename Y>
	inline void print_initial_N_sort_N_coefficients(const Y& const N_size, const vector<T>& const approximation_initial_N_vector, const vector<T>& const initial_N_by_N_prune_power_vector, const vector<T>& const initial_N_by_N_sum_deviation_vector, const vector<T>& const initial_N_by_N_run_time_vector, const vector<T>& const initial_N_by_N_approximation_time_vector, const vector<T>& const initial_N_by_N_knn_time_vector);

	//200710 print right endpoint of each segment in linked list or vector
	template<typename T>
	void print_each_segment_right_endpoint(const T& const array_data_structure);

	//200817 print right endpoint of each segment in linked list or vector with order
	template<typename T>
	void print_each_segment_right_endpoint_with_order(const T& const array_data_structure);

	//200710 print width of each segment in linked list or vector
	template<typename T>
	void print_each_segment_width(const T& const array_data_structure);

	//200817 print width of each segment in linked list or vector
	template<typename T>
	void print_each_segment_width_with_order(const T& const array_data_structure);

	//200710 print a&b of each segment in linked list or vector
	template<typename T>
	void print_each_segment_ab(const T& const array_data_structure);

	//200817 print a&b of each segment in linked list or vector
	template<typename T>
	void print_each_segment_ab_with_order(const T& const array_data_structure);

	//200814 print density of each segment in linked list or vector
	template<typename T>
	void print_each_segment_density(const T& const array_data_structure);

	//200929 print area difference of each segment in linked list or vector
	template<typename T>
	void print_each_segment_area_difference(const T& const array_data_structure);

	//200817 print density of each segment in linked list or vector
	template<typename T>
	void print_each_segment_density_with_order(const T& const array_data_structure);

	//200712 print several coefficients of each segment in linked list or vector
	template<typename T>
	void print_each_segment_coefficient(const T& const array_data_structure);

	//200817 print several coefficients of each segment with order in linked list or vector
	template<typename T>
	void print_each_segment_coefficient_with_order(const T& const array_data_structure);

	//200831 Print sting in vector
	template<typename T>
	void print_string_vector(const vector<T>& const string_vector);

	//200908 Print multimap second part iteration
	template<typename T>
	void print_multimap_second(const T& const multi_map) const;
	/*========================================================================================================================================*/

	//191203 initial Data Source structure.
	void initial_data_source(DATA_SOURCE& const data_source);

	vector<string>& get_read_multi_file_address(DATA_SOURCE& const data_source, const int& const file_id);

	//191223 read and normalized multi or single dimenson time series from file by time series id
	template<typename T>
	vector<T>& read_normalized_multi_time_series(const DATA_SOURCE& const data_source, const int& const time_series_id, vector<T>& const normalized_time_series_vector);
	/*====================================================================================================================================*/
};

#include "SHARE_TOOL.cpp"
#endif


