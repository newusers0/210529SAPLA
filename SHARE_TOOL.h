#pragma once
#ifndef SHARE_TOOL_H
#define SHARE_TOOL_H

#include "pch.h"
#include "GEOMETRY_TOOL.h"

//#pragma warning(disable : 4996)
TEMPLATE
class SHARE_TOOL {
public:
	/*#################################################      All Homogenous Data set Address      ##########################################################*/

	/*::::::::::::::::::::::::::::::::::::::::::::::::::::                    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	const string file_address_24 = "./191202DataSet/SingleDataAddress200103/AllFiles24/DataAddress181218.txt"; // Has all 24 files
	const string single_file_number_address_24 = "./191202DataSet/SingleDataAddress200103/AllFiles24/DataNumberLine200103.txt";// Has all 24 files
	const string single_file_length_address_24 = "./191202DataSet/SingleDataAddress200103/AllFiles24/DataLength200103.txt";// Has all 24 files
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	/*::::::::::::::::::::::::::::::::::::::::::::::::::::         :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	const string file_address_UCR2018_512_41 = "./191202DataSet/SingleDataAddress200103/Single512UCR2018/DataAddressUCR201014.txt"; // Has all UCR2018 41 files
	const string single_file_number_address_UCR2018_512_41 = "./191202DataSet/SingleDataAddress200103/Single512UCR2018/DataNumberLineUCR201014.txt";// Has all UCR2018 41 files
	const string single_file_length_address_UCR2018_512_41 = "./191202DataSet/SingleDataAddress200103/Single512UCR2018/DataLengthUCR201014.txt";// Has all UCR2018 41 files
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	/*::::::::::::::::::::::::::::::::::::::::::::::::::::      :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	const string file_address_UCR2018_1024_21 = "./191202DataSet/SingleDataAddress200103/Single1024UCR2018/DataAddress1024UCR2018201016.txt"; // 201019 Has all UCR2018 21 files
	const string single_file_number_address_UCR2018_1024_21 = "./191202DataSet/SingleDataAddress200103/Single1024UCR2018/DataNumberLine1024UCR2018201016.txt";// 201019 Has all UCR2018 21 files
	const string single_file_length_address_UCR2018_1024_21 = "./191202DataSet/SingleDataAddress200103/Single1024UCR2018/DataLength1024UCR2018201016.txt";// 201019 Has all UCR2018 21 files
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	/*::::::::::::::::::::::::::::::::::::::::::::::::::::           ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	const string file_address_UCR2018_1024_15 = "./191202DataSet/SingleDataAddress200103/Single15File1024UCR2018/DataAddress15File1024UCR2018201016.txt"; // 201019 In UCR2018 21 files length >= 1024, now choose 15files except 5 10 11 12 18
	const string single_file_number_address_UCR2018_1024_15 = "./191202DataSet/SingleDataAddress200103/Single15File1024UCR2018/DataNumberLine15File1024UCR2018201016.txt";// 201019 In UCR2018 21 files length >= 1024, now choose 15files except 5 10 11 12 18
	const string single_file_length_address_UCR2018_1024_15 = "./191202DataSet/SingleDataAddress200103/Single15File1024UCR2018/DataLength15File1024UCR2018201016.txt";// 201019  In UCR2018 21 files length >= 1024, now choose 15files except 5 10 11 12 18
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	//const string file_address = "./191202DataSet/SingleDataAddress200103/DataAddress200115.txt"; // No file 6 and 22 . Total 22 files
	/*..................................................200330..............................................................*/
	// Has burst time series  Total 21 files No file 6 , 17, 20 and 22
	const string file_address = "./191202DataSet/SingleDataAddress200103/DataAddress200116.txt";
	// No burst time series: size = 17 . No file 6, 12, 13, 15, 16, 17 and 22. 
	//const string file_address = "./SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataAddress200130.txt";
	/*.....................................................................................................................*/
	//const string file_address = "./191202DataSet/SingleDataAddress200103/ForTrueDataSize200126/DataAddress200126.txt"; // No file 6, 17, 20 and 22. Total 20 files
	/*..................................................200330..............................................................*/
	const string single_file_number_address = "./191202DataSet/SingleDataAddress200103/DataNumberLine200116.txt";// No file 6, 17and 22. Total 21 files
	/*......................................................................................................................*/
	//const string single_file_number_address = "./191202DataSet/SingleDataAddress200103/ForTrueDataSize200126/DataNumberLine200126.txt";// No file 6 , 17, 20 and 22. Total 20 files
	/*..................................................200330..............................................................*/
	const string single_file_length_address = "./191202DataSet/SingleDataAddress200103/DataLength200116.txt";// No file 6, 17and 22. Total 21 files
	/*......................................................................................................................*/
	//const string single_file_length_address = "./191202DataSet/SingleDataAddress200103/ForTrueDataSize200126/DataLength200126.txt";// No file 6 , 17, 20 and 22. Total 20 files

	/*####################################################################################################################################################################################*/

	/*---------------------------------------------- No burst time series Homogenous Data set Address. size = 17---------------------------------------------------------------*/
	/*..................................................        200330       ..............................................................*/
	const string no_burst_single_file_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataAddress200130.txt";
	const string no_burst_single_file_number_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataNumberLine200130.txt";// No burst time series : No file 6, 12, 13, 15, 16, 17 and 22. 
	const string no_burst_single_file_length_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataLength200130.txt";// No burst time series: No file 6, 12, 13, 15, 16, 17 and 22. 
	//200409 
	//const string file_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataAddress200130.txt";
	//const string single_file_number_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataNumberLine200130.txt";
	//const string single_file_length_address = "./191202DataSet/SingleDataAddress200103/NoBurstTimeSeries200130/NoBurstDataLength200130.txt";
	/*.....................................................................................................................................*/
	/*---------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/*----------------------------------------------Three Dimension Data set Address-------------------------------------------------------------*/
	//const string three_d_file_address = "./191202DataSet/3DDataAddress200103/3DDataAddress191229.txt";//191222
	//const string three_d_file_number_address = "./191202DataSet/3DDataAddress200103/3DDataNumber200103.txt";//200103
	//const string three_d_file_length_address = "./191202DataSet/3DDataAddress200103/3DDataLength200103.txt";//200103

	const string three_d_file_address = "./191202DataSet/3DDataAddress201031/3DDataAddress201031.txt";//201031
	const string three_d_file_number_address = "./191202DataSet/3DDataAddress201031/3DDataNumber201031.txt";//201031
	const string three_d_file_length_address = "./191202DataSet/3DDataAddress201031/3DDataLength201031.txt";//201031
	/*-------------------------------------------------------------------------------------------------------------------------------------------*/

	/*----------------------------- 191229 Project long single time series to multi time series---------------------*/
	const string single_to_multi_file_address = "./191202DataSet/SingleToMulti191229/SingleDataAddress191229.txt";//191229 Split long single time series to multi dimension time series
	const string write_single_to_multi_file_address = "./191202DataSet/SingleToMulti191229/WriteMultiDataAddress191229.txt";//191229
	/*--------------------------------------------------------------------------------------------------------------*/

	//191201 manipute data, real dataset, single data set, mixed data set
	struct DATA_SOURCE;

	struct INPUT_ARGUMENT;

	struct OUTPUT_ARGUMENT;//181214

	//200902 include sum deviation, approximation time.
	struct RESULT_RECORD;

	struct Y_PROJECTION_ARGUMENT;//191030 For burst time series, Y_Projection_Method

	//200225 evaluate coeffciets: local sum deviaiton, shift from best split method, time; Global sum deviation & time. for 5 split mthods min density, binary, intersection point, middle, best split id
	struct EVALUATION_COEFFCIENTS_SPLIT_ID;

	//************************************
	// Stuct:EVALUATION_BOUND
	// Qualifier: upper bound : max difference of time series or PLA.
	// date:210121
	// author:
	//************************************
	struct EVALUATION_BOUND;

	struct TIME;
	TIME time_record[100];

	struct ID_DIST;

	struct priorityDistanceEUC;

	struct APLA_COEFFICIENT;//191121 For PLA Coefficient same with APLA class

	struct RECTANGLE;//191127 for Rtree insersion. get min&max id value of rectangle

	template<typename T, typename Y, typename U>// 210906 first segment width, second segment width, remainder.
	struct SEGMENT_WIDTH_STRUCTRE;

	template<typename T>//for NULL
	void initialArray(T*& const test_array, const int& array_length);

	template<typename T>//for new memory
	void newArray(T*& const test_array, const int& array_length);

	template<typename T>//for NULL
	void deleteArray(T*& const test_array);
	//191128 clear memory of rectangle
	template<typename T>
	void delete_rectangle(T& const rectangle);

	/*===================================================Read data from txt file========================================================================================================*/
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

	//210603 vector file addresss vector 
	template<typename T>
	vector<T>& getMultiFoldToSingleByID(const vector<string>& const multi_file_name, const int& arity_d, const double& single_series_length, const int& const original_time_series_id, vector<T>& const original_time_series_vector);

	//210603
	template<typename T>
	std::fstream& GotoLine(std::fstream& file, const T& const line_num);

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

	//191203 Not miss first point, normalized
	// return g_time_series_length + 1.
	template<typename T>
	void getFileStreamByRowNoMiss(const string& const file_name, const int& const g_time_series_length, const int& const row_number, T& const original_time_series);

	//210603 read datas from one, normalize, and write
	template<typename T, typename Y>
	void get_normalize_write_one_file(Y& const data_source, const T& const data_id);

	//210603 read datas from one, normalize, and write
	template<typename T, typename Y, typename T1, typename  T2, typename  T3>
	void get_normalize_write_one_file(const Y& const file_name, const T& const data_id, const T1& const g_time_series_length, const T2& const row_number, const T3& const write_file_name);

	//miss first point normalized 210603
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

	//210824 get string and split by delimiter
	template<typename T, typename Y, typename U>
	vector<U>& get_string_by_delimiter(const T& const file_name, const Y& const size_string, vector<U>& const string_vector);

	//191115 get normalized original time series from txt file list according to file id
	template<typename T>
	void get_normalized_series_by_stream(const TOOL::INPUT_ARGUMENT& const input_argument, ifstream& const file_stream, T*& const original_time_series);
	//191115 get normalized original time series from txt file list according to file id
	template<typename T>
	const vector<T>& const get_normalized_series_by_address_list(TOOL::INPUT_ARGUMENT& const input_argument, const int& const file_id, vector<T>& const normalized_series_vector) const;
	//191203 For KNN, get & write mixed dataset, convert file_total_number to 1
	void initial_mixed_dataset(const string& const write_file_name, int& const file_total_number);
	/*=====================================================================================================================================================================================*/

	/*===============================convert vector to string========================*/
	// 191201 covnert vector list to string to write into txt file
	template<typename T>
	string convert_vector_to_string(const vector<T>& const vector_value);
	/*===============================================================================*/

	/*=======================================           Write Result        =======================================================================*/
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

	//210828 // with system time date. cover content
	template<typename T, typename Y>
	inline void write_single_result_time_stamp_no_txt(const T& write_file_name, Y& result);

	//with system time date. not cover content
	template<typename T>
	void writeSingleResultNoCoverTimeStamp(const string& write_file_name, T& result);

	template<typename T, typename Y>
	void writeSingleResult(Y& input_argument, const string& write_file_name, T& result);//181207 Add input_argument

	template<typename T, typename Y>
	void writeSingleResult(const string& const write_file_name, T& array_data, const Y& const array_len);//181218 Write array to file

	template<typename T>
	void writeSingleResult(const string& const write_file_name, vector<T>& array_data);//190501 Write vector to file

	template<typename T, typename Y>
	inline void write_single_result_no_txt(const Y& const write_file_name, const T& const result);//210828

	template<typename T, typename Y>
	inline void write_single_result_no_txt(const Y& const write_file_name, const vector<T>& const array_data);//210828

	//automatilly + "txt"
	template<typename T>
	void writeResultNoCover(const string& const write_file_name, vector<T>& array_data);//190624 Write vector to file, no coverage

	//210603 no + "txt"
	template<typename T, typename Y>
	void writeResultNoCover_no_txt(const Y& const write_file_name, vector<T>& array_data);//190624 Write vector to file, no coverage

	//No suffix like "txt"
	template<typename T>
	void writeResultNoCoverNoSuffix(const string& const write_file_name, vector<T>& array_data);//190624 Write vector to file, no coverage

	template<typename T>
	void writeDataTXTCover(const string& const write_file_name, T& data);//190625 Write T to file, coverage
	/*==============================================================================================================================================*/

	/*=======================================        Clear content of txt file      ================================================================*/
	inline void clearTXTFile(const string& const clear_file_name);//190625

	template<typename T>
	inline void clearTXTFile_no_txt(const T& const clear_file_name);//210603
	/*==============================================================================================================================================*/

	inline void recordStartTime(TIME& time);

	template<typename T>
	T& recordFinishTime(TIME& time, T& whole_first_run_time);

	inline long double recordFinishTime(TIME& time);// 180923

	template<typename T>//For sum deviation
	double& distanceEUC(const T& const input_argument, DataType*& const Q, DataType*& const C, double& const distance_euc);

	//191129 Euclidean distance.
	template<typename T>//For sum deviation
	long double distanceEUC(const vector<T>& const time_series_vector1, const vector<T>& const time_series_vector2);

	/*========================================================= Get Deviation ==============================================================================================*/
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
	/*==============================================================================================================================================================*/

	/*==========================================================K-NN Base Function=================================================================================================*/
	template<typename T, typename Y>
	void SimpleBaseKNNSearch(const T& const input_argument, DataType*& const g_query_time_series, priority_queue<Y, vector<Y>, priorityDistanceEUC >& q_base_queue);

	//191203 vector instead pointer. Normalized time series
	template<typename T, typename Y>
	multiset<pair<double, int>>& SimpleBaseKNNSearch(const T& const data_source, const vector<Y>& const query_time_series, multiset<pair<double, int>>& const knn_result_set);

	//191223 For multi dimension dataset. Vector instead pointer. Normalized time series
	template<typename T, typename Y, typename U>
	multiset<pair<U, int>>& SimpleBaseKNNSearchMulti(const T& const data_source, const vector<Y>& const query_time_series_vector, multiset<pair<U, int>>& const knn_result_set);
	/*==============================================================================================================================================================================*/

	/*------------------------------------          get Random Value     --------------------------------------------------------------*/
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

	/*======================================    Normalize  ========================================================*/

	template<typename T>
	void normalizeStandard(T*& const original_array, const int& const array_length, T*& const normalized_array);

	template<typename T>
	void normalizeStandard(const int& const array_length, T*& const normalized_array);
	//191115
	template<typename T>
	vector<T>& normalizeStandard(vector<T>& const normalized_array);
	/*=============================================================================================================*/

	/*---------------------------------------  210906  --------------------------------------------------*/
	template<typename T, typename Y, typename U>
	inline void get_equal_segment_width(const T& const n, const Y& const N, U& const equal_segment_width_struct);
	/*---------------------------------------------------------------------------------------------------*/

	/*============================ 210303 Transfer time us to s ====================================================*/
	template<typename T>
	void transfer_us_s(vector<T>& const time_vector);
	/*=============================================================================================================*/

	int transferSetToInt(const set<int>& const random_endpoint_set);//190623

	/* arr[] ---> Input Array
	n ---> Size of input array
	r ---> Size of a combination to be printed
	index ---> Current index in data[]
	data[] ---> Temporary array to store current combination
	 i ---> index of current element in arr[] */
	void combinationUtil(int arr[], int n, int r, int index, int data[], int i);//190623 /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	void printCombination(int arr[], int n, int r);//190623 // The main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()

	void combinationUtilJump(int arr[], int n, int r, int index, int data[], int i);//190624 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	void printCombinationJump(int arr[], int n, int r);//190624 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()

	void combinationUtilJumpVector(const vector<int>& const arr, int index, vector<int> data, int i, set<pair<double, vector<int>>>& const endpoint_collection);//190624 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	void printCombinationJumpVector(const vector<int>& const arr, const int& const r, set<pair<double, vector<int>>>& const endpoint_collection);//190624 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()

	/*====================================Combination===================================================================================*/
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
	/*==================================================================================================================================*/

	/*======================================================Print Result======================================================================*/
	void printOSPageSize();//for Windows system page_size

	template<typename T>//print argument like: n, N, K, max_node
	inline void printInputArgument(const T& const input_argument);

	template<typename T>//210830
	inline void print_input_argument_KNN_result(const T& const input_argument);

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
	/*..................................................  200326 Print initial_N evaluation  .......................................................*/
	template<typename T>
	inline void print_initial_N_coefficients(const vector<T>& const  approximation_initial_N_vector, const vector<T>& const total_initial_N_prune_power_vector, const vector<T>& const  total_initial_N_sum_deviation_vector, const vector<T>& const total_initial_N_run_time_vector, const vector<T>& const total_initial_N_approximation_time_vector, const vector<T>& const total_initial_N_knn_time_vector);
	/*...............................................................................................................................................*/

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

	/**====================================================Data Source structure===========================================================*/
	//191203 initial Data Source structure.
	void initial_data_source(DATA_SOURCE& const data_source);

	vector<string>& get_read_multi_file_address(DATA_SOURCE& const data_source, const int& const file_id);

	//210603 read and normalized multi or single dimenson time series from file by time series id
	template<typename T>
	vector<T>& read_normalized_multi_time_series(const DATA_SOURCE& const data_source, const int& const time_series_id, vector<T>& const normalized_time_series_vector);
	/*====================================================================================================================================*/

	template<typename T, typename Y>
	bool assert_same_vector(const T& const vector_0, const Y& const vector_1);
};

//TEMPLATE
//struct TOOL::INPUT_ARGUMENT {
//	int time_series_length = INF;//n
//	int point_multi_single_length = INF;//n*arity_d
//	int point_dimension = INF;//N
//	double initial_N = INF;
//	int point_multi_single_dimension = INF;//N*arity_d
//	double point_number = INF;
//
//	//191124 test for KNN, record current time series id
//	int point_id = INF;
//
//	int remainder = INF;
//	int segment_length_first = INF;//l=n/m+1
//	int segment_length_second = INF;//l=n/m
//
//	int query_time_series_id = INF;//191119 for KNN
//	int rtree_max_nodes = INF;
//	int K = INF;
//	string read_file_name;
//	string* write_file_name = nullptr;
//	string* read_multiple_file_name = nullptr;
//
//	//for Chebyshev
//	int degree_m = INF;//point_dimension for APCA & PLA
//	//for multi dimension time series
//	int arity_d = INF;// dimension of every data, for multi-dimensional case
//	int file_id = INF;// file id in file address list
//
//	/*----------------  191004 For Burst Time Series  --------*/
//	double burst_frquent_value = INF;
//	/*--------------------------------------------------------*/
//
//	//for KNN
//	double pruning_power = INF;
//	double sum_distance_euc = INF;
//
//	//for time
//	double approximation_query_time = NULL;
//	double knn_rest_part_time = NULL;
//
//	/*....191204 KNN time......*/
//	double representation_time = INF;
//	double build_rtree_time = INF;
//	double knn_total_time = INF;
//	double whole_run_time = INF;
//	/*.........................*/
//
//	//three part time
//	double navigate_index_time = NULL;// navigate time
//	double distance_lowbound_time = NULL; // distance chebyshev, PLA, APCA time
//	double distance_euc_time = NULL;// distance euclidean time
//
//	double result_accuracy = NULL;// KNN result accuracy
//
//	//I/O cost
//	double IO_cost = NULL;
//
//	/*-------------Rtree volumme & area calculaiton option--------------*/
//	/*
//	1 for PAA, 2 for APCA, 3 for PLA, 4 for Chebyshev, 5 APLA, 6 ICDE07
//	*/
//	int representation_option = NULL;
//	/*------------------------------------------------------------------*/
//
//	/*.......................200225 Split ID Option......................*/
//	// 1 Min Density, 2 Binary, 3 Intersection Point, 4 Middle Point, 5 Best id
//	double option_split_method = INF;
//	/*..................................................................*/
//
//	/*=============================Y-Projection======================================*/
//	//190730 distinguish use y proejction method
//	bool change_file = true;//191204 for y projection, get all segments and cluster segemnts
//	bool is_y_projection = false;
//	//int flat_segment_number = INF; //191209 the number of falt segments
//	//double total_flat_segment_length = INF;//191209 the length of all flat semgent
//	//191028 Use to check Y projection value
//	//int initial_N = INF;
//	//191028 // min max value for flat segment
//	//double flat_segment_min = INF;
//	//double flat_segment_max = INF; //threshhold
//	/*===============================================================================*/
//
//	~INPUT_ARGUMENT() {
//		time_series_length = INF;//n
//		point_multi_single_length = INF;//n*arity_d
//		point_dimension = INF;//N
//		initial_N = INF;
//		point_multi_single_dimension = INF;//N*arity_d
//		point_number = INF;
//		remainder = INF;
//		segment_length_first = INF;//l=n/m+1
//		segment_length_second = INF;//l=n/m
//
//		rtree_max_nodes = INF;
//		K = INF;
//
//		pruning_power = NULL;
//		sum_distance_euc = NULL;
//
//		approximation_query_time = NULL;
//		knn_rest_part_time = NULL;
//
//		navigate_index_time = NULL;// navigate time
//		distance_lowbound_time = NULL; // distance chebyshev, PLA, APCA time
//		distance_euc_time = NULL;// distance euclidean time
//		/*....191204 KNN time.......*/
//		representation_time = INF;
//		build_rtree_time = INF;
//		knn_total_time = INF;
//		whole_run_time = INF;
//		/*.........................*/
//
//		IO_cost = NULL;
//
//		result_accuracy = NULL;// KNN result accuracy
//
//		read_file_name.clear();
//		read_file_name.shrink_to_fit();
//
//		/*write_file_name.clear();
//		write_file_name.shrink_to_fit();*/
//		if (write_file_name != nullptr) {
//			//delete[] write_file_name;
//			write_file_name = nullptr;
//		}
//
//		if (read_multiple_file_name != nullptr) {
//			//delete[] read_multiple_file_name;
//			read_multiple_file_name = nullptr;
//		}
//		//read_file_name = nullptr;
//
//		/*-----------------------for KNN ------------------------*/
//		arity_d = INF;// dimension of every data, for multi-dimensional case
//		file_id = INF;// file id in file address list
//		/*-------------------------------------------------------*/
//
//		/*.......................Split ID Option............................*/
//		// 1 Min Density, 2 Binary, 3 Intersection Point, 4 Middle Point, 5 Best id
//		option_split_method = INF;
//		/*..................................................................*/
//
//		//for Chebyshev
//		degree_m = INF;
//		representation_option = INF;
//		/*=========================== Y-Projection =====================================*/
//		//190730 distinguish use y proejction method
//		change_file = true;
//		is_y_projection = false;
//		//flat_segment_number = INF; //191209 the number of falt segments
//		//total_flat_segment_length = INF;//191209 the length of all flat semgent
//		//initial_N = INF;
//		//flat_segment_min = INF;
//		//flat_segment_max = INF; //threshhold
//		/*==============================================================================*/
//	}
//};

//************************************
// Method:OUTPUT_ARGUMENT
// Qualifier:
// Input:
// Output:
// date:181214
// author:
//************************************
//TEMPLATE
//struct TOOL::OUTPUT_ARGUMENT {//181214
//	double max_deviation = NULL;
//	double sum_deviation = NULL;
//	long double run_time = INF;
//	int density_max_id = INF;
//	int density_min_id = INF;
//	double same_deviation_id_count = NULL;//190515
//	double diff_deviation_id_count = NULL;//190515
//
//	double sum_area = NULL;//190619
//	double sum_density = NULL;//190619
//	double sum_area0 = NULL;//190619
//	double sum_density0 = NULL;//190619
//	//double density_min_value = INF;
//
//	typename GEOMETRY::POINT min_density_segment;
//
//	EVALUATION_COEFFCIENTS_SPLIT_ID coefficents_split_id;
//
//	~OUTPUT_ARGUMENT() {
//		max_deviation = NULL;
//		sum_deviation = NULL;
//		run_time = INF;
//		density_max_id = INF;
//		density_min_id = INF;
//		//density_min_value = INF;
//		sum_area = NULL;//190619
//		sum_density = NULL;//190619
//		sum_area0 = NULL;//190619
//		sum_density0 = NULL;//190619
//		min_density_segment.~POINT();
//	}
//};

//************************************
// Method:Y_PROJECTION_ARGUMENT
// Qualifier:For burst time series, Y_Projection_Method
// Input:
// Output:
// date:191030
// author:
//************************************
//TEMPLATE
//struct TOOL::Y_PROJECTION_ARGUMENT {//191030 
//	int initial_N = INF;// Begin N for each file
//	//bool change_file = true; // 191204 check for KNN whter change file name
//	bool is_y_projection = false;// is burst time series
//	double flat_segment_min = INF;
//	double flat_segment_max = INF; //threshhold
//	double flat_segment_number = INF;//191211
//	double total_flat_segment_length = INF;//191211
//	map<double, int> whole_difference_map; //cluster the number of y_projection value
//	map<int, double> number_y_value_map;//191208 count the number of every y - projection value
//	/*
//	parameter:
//	   initial_number: cannot be NULL;
//	*/
//	Y_PROJECTION_ARGUMENT(const int& const initial_number) : initial_N(initial_number), is_y_projection(false), flat_segment_min(INF), flat_segment_max(INF) {}
//
//	Y_PROJECTION_ARGUMENT() : initial_N(INF), is_y_projection(false), flat_segment_min(INF), flat_segment_max(INF) {}
//
//	//191204 for KNN
//	Y_PROJECTION_ARGUMENT(const bool& const change_file) {
//		initial_N = INF;// Begin N for each file
//		//this->change_file = change_file; // 191115 check for KNN whter change file name
//		is_y_projection = false;// is burst time series
//		flat_segment_min = INF;
//		flat_segment_max = INF; //threshhold
//		flat_segment_number = INF;//191211
//		total_flat_segment_length = INF;//191211
//	}
//
//	//Y_PROJECTION_ARGUMENT(const int& const initial_number, const bool& const change_file) : initial_N(initial_number), change_file(change_file), is_y_projection(false), flat_segment_min(INF), flat_segment_max(INF) {}
//
//	~Y_PROJECTION_ARGUMENT() {
//		//change_file = true;
//		is_y_projection = false;// is burst time series
//		initial_N = INF;// Begin N for each file
//		flat_segment_min = INF;
//		flat_segment_max = INF; //threshhold
//		flat_segment_number = INF;//191211
//		total_flat_segment_length = INF;//191211
//
//		if (!whole_difference_map.empty()) {
//			whole_difference_map.clear();
//		}
//		number_y_value_map.clear();
//	}
//};

//************************************
// Method:EVALUATION_COEFFCIENTS_SPLIT_ID
// Qualifier: evaluate coeffciets: local sum deviaiton, shift from best split method, time; Global sum deviation & time. for 5 split mthods min density, binary, intersection point, middle, best split id
// Input:
// Output:
// date:200225
// author:
//************************************
//TEMPLATE
//struct TOOL::EVALUATION_COEFFCIENTS_SPLIT_ID {
//	/*.......... Accuracy & Time Split point selection method...........*/
//	//200219
//	//200225 sum deviation
//	double local_split_sum_deviation = 0;
//	double split_sum_deviation_min_density = INF;
//	double split_sum_deviation_binary = INF;
//	double split_sum_deviation_direct_intersection = INF;
//	double split_sum_deviation_middle = INF;//200224 middle point
//	double split_sum_deviation_best = INF;
//
//	//accuracy
//	double local_split_shift_abs = 0;
//	double split_accuracy_min_density = INF;
//	double split_accuracy_binary = INF;
//	double split_accuracy_direct_intersection = INF;
//	double split_accuracy_middle = INF;//200224
//	double split_id_best = INF;
//
//	//time
//	double local_split_time = 0;
//	double split_time_min_density = INF;
//	double split_time_binary = INF;
//	double split_time_direct_intersection = INF;
//	double split_time_middle = INF;//200224
//	double split_time_best = INF;
//	/*..................................................................*/
//
//	EVALUATION_COEFFCIENTS_SPLIT_ID() {
//		/*.......... Accuracy & Time Split point selection method...........*/
//		//200219
//		local_split_sum_deviation = INF;
//		split_sum_deviation_min_density = INF;
//		split_sum_deviation_binary = INF;
//		split_sum_deviation_direct_intersection = INF;
//		split_sum_deviation_middle = INF;
//		split_sum_deviation_best = INF;
//
//		local_split_shift_abs = INF;
//		split_accuracy_min_density = INF;
//		split_accuracy_binary = INF;
//		split_accuracy_direct_intersection = INF;
//		split_accuracy_middle = INF;
//		split_id_best = INF;
//
//		local_split_time = INF;
//		split_time_min_density = INF;
//		split_time_binary = INF;
//		split_time_direct_intersection = INF;
//		split_time_middle = INF;
//		split_time_best = INF;
//		/*..................................................................*/
//	}
//
//
//	~EVALUATION_COEFFCIENTS_SPLIT_ID() {
//		/*.......... Accuracy & Time Split point selection method...........*/
//		//200219
//		local_split_sum_deviation = INF;
//		split_sum_deviation_min_density = INF;
//		split_sum_deviation_binary = INF;
//		split_sum_deviation_direct_intersection = INF;
//		split_sum_deviation_middle = INF;
//		split_sum_deviation_best = INF;
//
//		local_split_shift_abs = INF;
//		split_accuracy_min_density = INF;
//		split_accuracy_binary = INF;
//		split_accuracy_direct_intersection = INF;
//		split_accuracy_middle = INF;
//		split_id_best = INF;
//
//		local_split_time = INF;
//		split_time_min_density = INF;
//		split_time_binary = INF;
//		split_time_direct_intersection = INF;
//		split_time_middle = INF;
//		split_time_best = INF;
//		/*..................................................................*/
//	}
//
//};

//TEMPLATE
//struct TOOL::TIME {
//	_LARGE_INTEGER time_start;
//	_LARGE_INTEGER time_over;   //finish time
//	double dqFrequency = NULL;      //timer frequency
//	double run_time = NULL;
//
//	~TIME() {
//		time_start.QuadPart = NULL;
//		time_over.QuadPart = NULL;   //finish time
//		dqFrequency = NULL;      //timer frequency
//		run_time = NULL;
//	}
//};

//#include "SHARE_TOOL.cpp"

//template class SHARE_TOOL<DataType>;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/******************************************************************************************************************************************/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//************************************
// Method: DATA_SOURCE
// Qualifier: KNN for APLA & ICDE07
// Input:manipute data, real dataset, single data set, mixed data set
// Output:
// date:191201
// author:
//************************************
TEMPLATE
struct TOOL::DATA_SOURCE {
	int data_type = INF; //0 manipulate data set; 1 real data set, homogenous(single) data; 2 real data set, heterogeneous(mixed) data. 3 multi single data. 4 multi mixed data

	/***    0 has burst data 21 / 1 no burst data 17  ***/
	int option_has_burst_data = INF;
	/****************************************************/

	int time_series_dimension = INF;//191223 the dimenson of time series
	vector<size_t> data_list;//id of single data address
	int single_file_number = INF;// number of file data
	int single_point_number = INF;// number of time series in single file
	double single_time_series_length = INF;// length of time series length
	string file_address_all;// all file address 210603
	const string* file_name_pointer = nullptr;//data file name
	vector<string> file_address_vector;//191222 vector of string file name 

	int bigger_account = 0;

	/*------------------------mixed data-------------------*/
	int file_operation_number = INF;// for mixed dataset, and multi to single dataset: single file number / dataset dimenson
	int point_number = INF;
	/*-----------------------------------------------------*/
	string read_file_address;

	/*--------------muilti dimension data------------------*/
	vector<string> read_file_address_vector;//191223 
	double multi_single_time_series_length = INF;//191223
	/*-----------------------------------------------------*/

	vector<int> time_series_length_vector;
	vector<int> time_series_number_vector;

	DATA_SOURCE() {};
	DATA_SOURCE(const int& const data_type, const int& const time_series_dimension, const vector<size_t>& const data_list, const int& const single_file_number, const int& const single_point_number, const int& const single_time_series_length, const string& const data_file_name) : data_type(data_type), time_series_dimension(time_series_dimension), data_list(data_list), single_file_number(single_file_number), single_point_number(single_point_number), single_time_series_length(single_time_series_length), file_name_pointer(&data_file_name) {}

	~DATA_SOURCE() {
		data_type = INF; //0 manipulate data set; 1 real data set, homogenous(single) data; 2 real data set, heterogeneous(mixed) data. 3 multi single data. 4 multi mixed data
		time_series_dimension = INF;
		data_list.clear();//id of single data address
		data_list.shrink_to_fit();
		single_file_number = INF;// number of file data
		single_point_number = INF;// number of time series in single file
		single_time_series_length = INF;// length of single time series length
		file_name_pointer = nullptr;//mixed data file name
		file_address_vector.clear();
		file_address_vector.shrink_to_fit();
		file_operation_number = INF;
		point_number = INF;
		read_file_address.clear();
		read_file_address.shrink_to_fit();
		read_file_address_vector.clear();
		read_file_address_vector.shrink_to_fit();
		multi_single_time_series_length = INF;//191223

		time_series_length_vector.clear();//200103
		time_series_length_vector.shrink_to_fit();//200103
		time_series_number_vector.clear();//200103
		time_series_number_vector.shrink_to_fit();//200103
	}
};

TEMPLATE
struct TOOL::INPUT_ARGUMENT {

	int time_series_length = INF;//n
	int point_multi_single_length = INF;//n*arity_d
	int point_dimension = INF;//N
	double initial_N = INF;
	int point_multi_single_dimension = INF;//N*arity_d
	double point_number = INF;

	//191124 test for KNN, record current time series id
	int point_id = INF;

	int remainder = INF;
	int segment_length_first = INF;//l=n/m+1
	int segment_length_second = INF;//l=n/m

	int query_time_series_id = INF;//191119 for KNN
	int rtree_max_nodes = INF;
	int K = INF;
	string read_file_name;
	string* write_file_name = nullptr;
	string* read_multiple_file_name = nullptr;

	//for Chebyshev
	int degree_m = INF;//point_dimension for APCA & PLA
	//for multi dimension time series
	int arity_d = INF;// dimension of every data, for multi-dimensional case
	int file_id = INF;// file id in file address list

	/*--  191004 For Burst Time Series --*/
	double burst_frquent_value = INF;
	/*-----------------------------------*/

	/*--------     for KNN      ---------*/
	long double prune_power_combine = INF;// prune power combine 201221
	long double pruning_power = INF;
	double sum_distance_euc = INF;
	/*-----------------------------------*/

	//for time
	double approximation_query_time = NULL;
	double knn_rest_part_time = NULL;

	/*------   191204 KNN time   ------*/
	double representation_time = INF;
	double build_rtree_time = INF;
	double knn_total_time = INF;
	double knn_total_time_has_IO = INF;//210606
	double whole_run_time = INF;
	double whole_run_time_has_IO = INF;//210606
	/*---------------------------------*/

	//three part time
	double navigate_index_time = NULL;// navigate time
	double distance_lowbound_time = NULL; // distance chebyshev, PLA, APCA time
	double distance_euc_time = NULL;// distance euclidean time

	double result_accuracy = NULL;// KNN result accuracy

	//I/O cost
	long double IO_cost = NULL;

	/*-------------Rtree volumme & area calculaiton option--------------*/
	/*1 for PAA, 2 for APCA, 3 for PLA, 4 for Chebyshev, 5 APLA, 6 ICDE07*/
	int representation_option = NULL;
	/*------------------------------------------------------------------*/

	/*------------------------ 210618 Tree  ----------------------------*/
	int option_tree = INF;
	/*------------------------------------------------------------------*/

	/*-------------          200225 Split ID Option        -------------*/
	// 1 Min Density, 2 Binary, 3 Intersection Point, 4 Middle Point, 5 Best id
	double option_split_method = INF;
	/*------------------------------------------------------------------*/

	/*=============================Y-Projection======================================*/
	//190730 distinguish use y proejction method
	bool change_file = true;//191204 for y projection, get all segments and cluster segemnts
	bool is_y_projection = false;
	//int flat_segment_number = INF; //191209 the number of falt segments
	//double total_flat_segment_length = INF;//191209 the length of all flat semgent
	//191028 Use to check Y projection value
	//int initial_N = INF;
	//191028 // min max value for flat segment
	//double flat_segment_min = INF;
	//double flat_segment_max = INF; //threshhold
	/*===============================================================================*/

	/*----------      200901 option print each time series result      --------------*/
	bool print_each_result = false;
	/*-------------------------------------------------------------------------------*/

	bool is_MS_effect = false;
	bool is_optimization_effect = false;
	bool is_same_best_case = false;

	/*--201028 Count probability of max deviation point is minmax left right points--*/
	long double number_point_max_deviation_true = INF;
	long double number_point_max_deviation_false = INF;

	long double number_not_smaller_than_sum_deviation = INF;
	long double number_smaller_than_sum_deviation = INF;

	long double number_not_smaller_than_sum_deviation_pow = INF;
	long double number_smaller_than_sum_deviation_pow = INF;
	/*-------------------------------------------------------------------------------*/

	INPUT_ARGUMENT() {
	}

	template<typename T, typename Y>
	INPUT_ARGUMENT(const T& const n, const Y& const N) {
	
		time_series_length = n;
		point_dimension = N;
		/*===================================For Chebyshev coefficient=============================================*/
		degree_m = point_dimension - 1;//degree of chebyshev polymonail n=m+1
		/*==========================================================================================================*/
		/*===================================For PLA coefficient====================================================*/
		remainder = int(time_series_length) % int(point_dimension);//For PLA
		segment_length_second = (time_series_length - remainder) / point_dimension;
		segment_length_first = segment_length_second + 1;
	    /*==========================================================================================================*/
	}

	~INPUT_ARGUMENT() {

		time_series_length = INF;// n
		point_multi_single_length = INF;// n * arity_d
		point_dimension = INF;// N
		initial_N = INF;
		point_multi_single_dimension = INF;// N * arity_d
		point_number = INF;
		remainder = INF;
		segment_length_first = INF;// l = n / m+1
		segment_length_second = INF;// l = n / m

		rtree_max_nodes = INF;
		K = INF;

		/*-----   for KNN  -------*/
		prune_power_combine = INF;// prune power combine 201221
		pruning_power = NULL;
		sum_distance_euc = NULL;
		/*------------------------*/

		approximation_query_time = NULL;
		knn_rest_part_time = NULL;

		navigate_index_time = NULL;// navigate time
		distance_lowbound_time = NULL; // distance chebyshev, PLA, APCA time
		distance_euc_time = NULL;// distance euclidean time

		/*....191204 KNN time.......*/
		representation_time = INF;
		option_tree = INF;
		build_rtree_time = INF;
		knn_total_time = INF;
		knn_total_time_has_IO = INF;
		whole_run_time = INF;
		whole_run_time_has_IO = INF;//210606
		/*.........................*/

		IO_cost = NULL;

		result_accuracy = NULL;// KNN result accuracy

		read_file_name.clear();
		read_file_name.shrink_to_fit();

		/*write_file_name.clear();
		write_file_name.shrink_to_fit();*/
		if (write_file_name != nullptr) {
			//delete[] write_file_name;
			write_file_name = nullptr;
		}

		if (read_multiple_file_name != nullptr) {
			//delete[] read_multiple_file_name;
			read_multiple_file_name = nullptr;
		}
		//read_file_name = nullptr;

		/*-----------------------for KNN ------------------------*/
		arity_d = INF;// dimension of every data, for multi-dimensional case
		file_id = INF;// file id in file address list
		/*-------------------------------------------------------*/

		/*.......................Split ID Option............................*/
		// 1 Min Density, 2 Binary, 3 Intersection Point, 4 Middle Point, 5 Best id
		option_split_method = INF;
		/*..................................................................*/

		//for Chebyshev
		degree_m = INF;
		representation_option = INF;
		/*=========================== Y-Projection =====================================*/
		//190730 distinguish use y proejction method
		change_file = true;
		is_y_projection = false;
		//flat_segment_number = INF; //191209 the number of falt segments
		//total_flat_segment_length = INF;//191209 the length of all flat semgent
		//initial_N = INF;
		//flat_segment_min = INF;
		//flat_segment_max = INF; //threshhold
		/*==============================================================================*/

		/*..............200901 option print each time series result.....................*/
		print_each_result = false;
		/*...............................................................................*/

		/*--201028 Count probability of max deviation point is minmax left right points--*/
		number_point_max_deviation_true = INF;
		number_point_max_deviation_false = INF;

		number_not_smaller_than_sum_deviation = INF;
		number_smaller_than_sum_deviation = INF;

		number_not_smaller_than_sum_deviation_pow = INF;
		number_smaller_than_sum_deviation_pow = INF;
		/*-------------------------------------------------------------------------------*/
	}
};
//
////************************************
//// Method:OUTPUT_ARGUMENT
//// Qualifier:
//// Input:
//// Output:
//// date:181214
//// author:
////************************************
TEMPLATE
struct TOOL::OUTPUT_ARGUMENT {//181214
	double max_deviation = NULL;
	long double max_deviation_av = NULL;
	double sum_deviation = NULL;
	long double run_time = INF;
	int density_max_id = INF;
	int density_min_id = INF;
	double same_deviation_id_count = NULL;//190515
	double diff_deviation_id_count = NULL;//190515

	double sum_area = NULL;//190619
	double sum_density = NULL;//190619
	double sum_area0 = NULL;//190619
	double sum_density0 = NULL;//190619
	//double density_min_value = INF;

	typename GEOMETRY::POINT min_density_segment;

	EVALUATION_COEFFCIENTS_SPLIT_ID coefficents_split_id;

	EVALUATION_BOUND evaluation_bound;

	long double max_deviation_multiple_width = INF;

	long double sum_deviation_icde07 = INF;
	long double max_deviation_icde07 = INF;
	long double max_deviation_av_icde07 = INF;

	~OUTPUT_ARGUMENT() {
		max_deviation = NULL;
		max_deviation_av = NULL;
		sum_deviation = NULL;
		run_time = INF;
		density_max_id = INF;
		density_min_id = INF;
		//density_min_value = INF;
		sum_area = NULL;//190619
		sum_density = NULL;//190619
		sum_area0 = NULL;//190619
		sum_density0 = NULL;//190619
		min_density_segment.~POINT();

		max_deviation_multiple_width = INF;

		sum_deviation_icde07 = INF;
		max_deviation_icde07 = INF;
		max_deviation_av_icde07 = INF;
	}
};

//************************************
// Method:RESULT_RECORD
// Qualifier: include sum deviation, approximation time.
// Input:
// Output:
// date:200902 
// author:
//************************************
TEMPLATE
struct TOOL::RESULT_RECORD {
	long double max_deviation = INF;
	long double max_deviation_av = INF;//210912
	long double max_deviation_multiple_width = INF;

	long double sum_deviation = INF;
	long double representation_time = INF;

	~RESULT_RECORD() {
		max_deviation = INF;
		max_deviation_av = INF;
		max_deviation_multiple_width = INF;
		sum_deviation = INF;
		representation_time = INF;
	}
};

//************************************
// Method:Y_PROJECTION_ARGUMENT
// Qualifier:For burst time series, Y_Projection_Method
// Input:
// Output:
// date:191030
// author:
//************************************
TEMPLATE
struct TOOL::Y_PROJECTION_ARGUMENT {//191030 
	int initial_N = INF;// Begin N for each file
	//bool change_file = true; // 191204 check for KNN whter change file name
	bool is_y_projection = false;// is burst time series
	double flat_segment_min = INF;
	double flat_segment_max = INF; //threshhold
	double flat_segment_number = INF;//191211
	double total_flat_segment_length = INF;//191211
	map<double, int> whole_difference_map; //cluster the number of y_projection value
	map<int, double> number_y_value_map;//191208 count the number of every y - projection value
	/*
	parameter:
	   initial_number: cannot be NULL;
	*/
	Y_PROJECTION_ARGUMENT(const int& const initial_number) : initial_N(initial_number), is_y_projection(false), flat_segment_min(INF), flat_segment_max(INF) {}

	Y_PROJECTION_ARGUMENT() : initial_N(INF), is_y_projection(false), flat_segment_min(INF), flat_segment_max(INF) {}

	//191204 for KNN
	Y_PROJECTION_ARGUMENT(const bool& const change_file) {
		initial_N = INF;// Begin N for each file
		//this->change_file = change_file; // 191115 check for KNN whter change file name
		is_y_projection = false;// is burst time series
		flat_segment_min = INF;
		flat_segment_max = INF; //threshhold
		flat_segment_number = INF;//191211
		total_flat_segment_length = INF;//191211
	}

	//Y_PROJECTION_ARGUMENT(const int& const initial_number, const bool& const change_file) : initial_N(initial_number), change_file(change_file), is_y_projection(false), flat_segment_min(INF), flat_segment_max(INF) {}

	~Y_PROJECTION_ARGUMENT() {
		//change_file = true;
		is_y_projection = false;// is burst time series
		initial_N = INF;// Begin N for each file
		flat_segment_min = INF;
		flat_segment_max = INF; //threshhold
		flat_segment_number = INF;//191211
		total_flat_segment_length = INF;//191211

		if (!whole_difference_map.empty()) {
			whole_difference_map.clear();
		}
		number_y_value_map.clear();
	}
};
//
////************************************
//// Method:EVALUATION_COEFFCIENTS_SPLIT_ID
//// Qualifier: evaluate coeffciets: local sum deviaiton, shift from best split method, time; Global sum deviation & time. for 5 split mthods min density, binary, intersection point, middle, best split id
//// Input:
//// Output:
//// date:200225
//// author:
////************************************
TEMPLATE
struct TOOL::EVALUATION_COEFFCIENTS_SPLIT_ID {
	/*.......... Accuracy & Time Split point selection method...........*/
	//200219
	//200225 sum deviation
	double local_split_sum_deviation = 0;
	double split_sum_deviation_min_density = INF;
	double split_sum_deviation_binary = INF;
	double split_sum_deviation_direct_intersection = INF;
	double split_sum_deviation_middle = INF;//200224 middle point
	double split_sum_deviation_best = INF;

	//accuracy
	double local_split_shift_abs = 0;
	double split_accuracy_min_density = INF;
	double split_accuracy_binary = INF;
	double split_accuracy_direct_intersection = INF;
	double split_accuracy_middle = INF;//200224
	double split_id_best = INF;

	//time
	double local_split_time = 0;
	double split_time_min_density = INF;
	double split_time_binary = INF;
	double split_time_direct_intersection = INF;
	double split_time_middle = INF;//200224
	double split_time_best = INF;
	/*..................................................................*/

	EVALUATION_COEFFCIENTS_SPLIT_ID() {
		/*.......... Accuracy & Time Split point selection method...........*/
		//200219
		local_split_sum_deviation = INF;
		split_sum_deviation_min_density = INF;
		split_sum_deviation_binary = INF;
		split_sum_deviation_direct_intersection = INF;
		split_sum_deviation_middle = INF;
		split_sum_deviation_best = INF;

		local_split_shift_abs = INF;
		split_accuracy_min_density = INF;
		split_accuracy_binary = INF;
		split_accuracy_direct_intersection = INF;
		split_accuracy_middle = INF;
		split_id_best = INF;

		local_split_time = INF;
		split_time_min_density = INF;
		split_time_binary = INF;
		split_time_direct_intersection = INF;
		split_time_middle = INF;
		split_time_best = INF;
		/*..................................................................*/
	}


	~EVALUATION_COEFFCIENTS_SPLIT_ID() {
		/*.......... Accuracy & Time Split point selection method...........*/
		//200219
		local_split_sum_deviation = INF;
		split_sum_deviation_min_density = INF;
		split_sum_deviation_binary = INF;
		split_sum_deviation_direct_intersection = INF;
		split_sum_deviation_middle = INF;
		split_sum_deviation_best = INF;

		local_split_shift_abs = INF;
		split_accuracy_min_density = INF;
		split_accuracy_binary = INF;
		split_accuracy_direct_intersection = INF;
		split_accuracy_middle = INF;
		split_id_best = INF;

		local_split_time = INF;
		split_time_min_density = INF;
		split_time_binary = INF;
		split_time_direct_intersection = INF;
		split_time_middle = INF;
		split_time_best = INF;
		/*..................................................................*/
	}

};

//************************************
// Stuct:EVALUATION_BOUND
// Qualifier: upper bound : max difference of time series or PLA.
// date:210121
// author:
//************************************
TEMPLATE
struct TOOL::EVALUATION_BOUND {
	map<size_t, size_t> map_bound_increment;
	map<size_t, size_t> map_bound_merge;
	map<size_t, size_t> map_bound_split;
	map<size_t, size_t> map_bound_move;

	long double difference_id_max_deviation_vs_height_diff = NULL;

	~EVALUATION_BOUND() {
		map_bound_increment.clear();
		map_bound_merge.clear();
		map_bound_split.clear();
		map_bound_move.clear();

		difference_id_max_deviation_vs_height_diff = NULL;
	}
};

TEMPLATE
struct TOOL::TIME {
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
struct TOOL::ID_DIST {
	double dist = NULL;			//distance.
	int id = NULL;      // trajectory ID.
};

TEMPLATE
struct TOOL::priorityDistanceEUC {//small to big
	bool operator ()(const ID_DIST& a, const ID_DIST& b) {
		return a.dist > b.dist;
	}
};

//************************************
// Stuct:APLA_COEFFICIENT
// Qualifier:
// date:191112
// author:
//************************************
TEMPLATE
struct TOOL::APLA_COEFFICIENT {//181212 For PLA Coefficient
	long double a = INF;//ax+b
	long double b = INF;//ax+b
	long double right_endpoint = INF;
	int segmentNum = INF; //the dimension of the index(eg. MBR, imput parameter)

	long double a_minuend = INF;//(l-1)/2
	long double a_divisor = INF;//l(l-1)(l+1)
	long double b_minuend = INF;//2l-1
	long double b_divisor = INF;//l(l+1)

	APLA_COEFFICIENT() {
		a = INF;//ax+b
		b = INF;//ax+b
		right_endpoint = INF;
		segmentNum = INF; //the dimension of the index(eg. MBR, imput parameter)
		//because initial segment length is 2, l=2
		a_minuend = 0.5;//(l-1)/2
		a_divisor = 6;//l(l-1)(l+1)
		b_minuend = 3;//2l-1
		b_divisor = 6;//l(l+1)
	}

	~APLA_COEFFICIENT() {
		segmentNum = INF;
		a = INF;
		b = INF;
		right_endpoint = INF;

		a_minuend = INF;//(l-1)/2
		a_divisor = INF;//l(l-1)(l+1)
		b_minuend = INF;//2l-1
		b_divisor = INF;//l(l+1)
	}
};

//************************************
// Stuct:RECTANGLE
// Qualifier:for Rtree insersion. get min&max id value of rectangle
// date:191127
// author:
//************************************
TEMPLATE
struct TOOL::RECTANGLE {//191127 for Rtree insersion. get min&max id value of rectangle
	vector<DataType> m_min;
	vector<DataType> m_max;

	RECTANGLE(const size_t& const vector_size) {
		m_min.resize(vector_size, INF);
		m_max.resize(vector_size, INF);
	}

	~RECTANGLE() {
		m_min.clear();
		m_min.shrink_to_fit();
		m_max.clear();
		m_max.shrink_to_fit();
	}
};

TEMPLATE
template<typename T, typename Y, typename U>// 210906 first segment width, second segment width, remainder.
struct TOOL::SEGMENT_WIDTH_STRUCTRE {
	T  remainder = INF;
	Y segment_length_first = INF;
	U segment_length_second = INF;

	SEGMENT_WIDTH_STRUCTRE() {
		remainder = INF;
		segment_length_first = INF;
		segment_length_second = INF;
	}

	template<typename T1, typename T2>
	SEGMENT_WIDTH_STRUCTRE(const T1& const n, const T2& const N){
		remainder = int(n) % int(N);//For PLA
		segment_length_second = (n - remainder) / N;
		segment_length_first = segment_length_second + 1;
	}

	~SEGMENT_WIDTH_STRUCTRE() {
		remainder = INF;
		segment_length_first = INF;
		segment_length_second = INF;
	}

};

/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

TEMPLATE
template<typename T>//for NULL
void TOOL::initialArray(T*& const test_array, const int& array_length) {
	fill_n(test_array, array_length, NULL);
}

TEMPLATE
template<typename T>//for new memory 181112
void TOOL::newArray(T*& const test_array, const int& array_length) {
	test_array = new T[array_length];
	fill_n(test_array, array_length, NULL);
}

TEMPLATE
template<typename T>//for delete  181112
void TOOL::deleteArray(T*& const test_array) {
	if (test_array != nullptr) {
		delete[] test_array;
		test_array = nullptr;
	}
}

//191128 clear memory of rectangle
TEMPLATE
template<typename T>
void TOOL::delete_rectangle(T& const rectangle) {
	rectangle.m_min.clear();
	rectangle.m_min.shrink_to_fit();
	rectangle.m_max.clear();
	rectangle.m_max.shrink_to_fit();
}

//miss first point. Do not read first point
TEMPLATE
void TOOL::getFileStreamByID(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, DataType*& const original_time_series) {

#ifdef _DEBUG
	assert(original_time_series_id != INF);
#endif

	string fs_row_string;
	string fs_row_number;

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& get split delimiter &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	ifstream file_stream_temp = ifstream(file_name);
	assert(file_stream_temp);
	char delimiter;
	while (file_stream_temp.get(delimiter) && delimiter != '\t' && delimiter != ',') {
	}
	file_stream_temp.close();
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	int i = 0;
	while (!file_stream.eof() && file_stream.is_open())
	{
		file_stream.good();
		file_stream.fail();

		/*=====================================*/
		//file_stream >> fs_row_string;
		getline(file_stream, fs_row_string, '\n');
		/*=====================================*/

		if (i == original_time_series_id) {
			stringstream sstr(fs_row_string);
			int string_id = -1;// miss first point
			while (getline(sstr, fs_row_number, delimiter) && string_id < g_time_series_length) {
				//while ((getline(sstr, fs_row_number, ',') || getline(sstr, fs_row_number, '\t')) && string_id < g_time_series_length) {
				if (string_id > -1) {
					original_time_series[string_id] = stod(fs_row_number);
					assert(original_time_series[string_id] != INF);
				}
				string_id++;
			}
			break;
		}

		fs_row_string.clear();
		i++;
	}

	assert(original_time_series[int(g_time_series_length) - 1] != INF);

	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();
}

//191115 vector. miss first point. file_stream >> fs_row_string;
TEMPLATE
template<typename T>
void TOOL::getFileStreamByID(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, vector<T>& const original_time_series) {

#ifdef _DEBUG
	assert(original_time_series.size() == g_time_series_length && original_time_series_id != INF);
#endif

	string line_string;
	string point_value_string;


	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& get split delimiter &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	ifstream file_stream_temp = ifstream(file_name);
	assert(file_stream_temp);
	char delimiter;
	while (file_stream_temp.get(delimiter) && delimiter != '\t' && delimiter != ',') {
	}
	file_stream_temp.close();
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	int id_line = 0;
	while (!file_stream.eof() && file_stream.is_open())
	{
		file_stream.good();
		file_stream.fail();

		/*=========================================*/
		//file_stream >> line_string;
		getline(file_stream, line_string, '\n');
		/*=========================================*/

		if (id_line == original_time_series_id) {
			stringstream sstr(line_string);
			int string_id = -1;
			while (getline(sstr, point_value_string, delimiter) && string_id < g_time_series_length) {
				//while ((getline(sstr, fs_row_number, ',') || getline(sstr, fs_row_number,'\t')) && string_id < g_time_series_length) {
				if (string_id > -1) {
					original_time_series[string_id] = stod(point_value_string);
					assert(original_time_series[string_id] != INF);
				}
				string_id++;
			}
			break;
		}

		line_string.clear();
		id_line++;
	}

	assert(original_time_series[int(g_time_series_length) - 1] != INF);

#ifdef _DEBUG
	for (auto&& au : original_time_series) {
		assert(au != INF);
	}
#endif

	file_stream.close();
	line_string.clear();
	line_string.shrink_to_fit();
	point_value_string.clear();
	point_value_string.shrink_to_fit();
}

//do not miss first point.
TEMPLATE
void TOOL::getFileStreamByID0(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, DataType*& const original_time_series) {
	assert(0);

#ifdef _DEBUG
	assert(original_time_series_id != INF);
#endif

	string fs_row_string;
	string fs_row_number;

	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	int i = 0;
	while (!file_stream.eof() && file_stream.is_open())
	{
		file_stream.good();
		file_stream.fail();

		//file_stream >> fs_row_string;
		getline(file_stream, fs_row_string);//different from above methods
		if (i == original_time_series_id) {
			stringstream sstr(fs_row_string);
			int string_id = 0;
			while (getline(sstr, fs_row_number, ',') && string_id < g_time_series_length) {
				original_time_series[string_id] = stod(fs_row_number);
				assert(original_time_series[string_id] != NULL);
				string_id++;
			}
			break;
		}

		fs_row_string.clear();
		i++;
	}
	assert(original_time_series[int(g_time_series_length) - 1] != NULL);
	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();
}

//Not miss first point.Vector 190702
TEMPLATE
template<typename T>
vector<T>& TOOL::getFileStreamByID0Vector(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, vector<T>& const result_vector) {
	assert(0);

	result_vector.clear();
	result_vector.shrink_to_fit();
#ifdef _DEBUG
	assert(original_time_series_id != INF && result_vector.empty());
#endif

	string fs_row_string;
	string fs_row_number;

	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	int i = 0;
	while (!file_stream.eof() && file_stream.is_open())
	{
		file_stream.good();
		file_stream.fail();

		//file_stream >> fs_row_string;
		getline(file_stream, fs_row_string);//different from above methods
		if (i == original_time_series_id) {
			stringstream sstr(fs_row_string);
			int string_id = 0;
			while (getline(sstr, fs_row_number, ',') && string_id < g_time_series_length) {
				result_vector.emplace_back(stod(fs_row_number));
				assert(result_vector[string_id] != INF);
				string_id++;
			}
			break;
		}

		fs_row_string.clear();
		i++;
	}
	assert(result_vector[int(g_time_series_length) - 1] != INF);
	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();

	return result_vector;
}


//***************************************************************
// Method:getFileStreamByID0Vector
// Qualifier:Not miss first point. No length .Vector 200103 
// Input:
// Output:
// date:210603
// author:
//***************************************************************
TEMPLATE//Not miss first point. No length .Vector 200103 
template<typename T>
vector<T>& TOOL::getFileStreamByID0Vector(const string& file_name, const int& const original_time_series_id, vector<T>& const result_vector) {

#ifdef _DEBUG
	assert(original_time_series_id != INF);
#endif

	result_vector.clear();
	result_vector.shrink_to_fit();

	string fs_row_string;
	string fs_row_number;

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	ifstream file_stream_temp = ifstream(file_name);
	char delimiter = '\t';
	while (file_stream_temp.get(delimiter) && delimiter != '\t' && delimiter != ',') {
	}
	file_stream_temp.close();
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	int id_line = 0;
	while (!file_stream.eof() && file_stream.is_open())
	{
		file_stream.good();
		file_stream.fail();

		//getline(file_stream, fs_row_string);//different from above methods
		getline(file_stream, fs_row_string, '\n');

		if (id_line == original_time_series_id) {
			stringstream sstr(fs_row_string);
			//int string_id = 0;
			while (getline(sstr, fs_row_number, delimiter)) {
				//while (getline(sstr, fs_row_number, ',')) {
				result_vector.emplace_back(stod(fs_row_number));
				//assert(result_vector[string_id] != INF);
				//string_id++;
			}
			break;
		}

		fs_row_string.clear();
		id_line++;
	}

	//#ifdef _DEBUG
	for (auto&& au : result_vector) {
		assert(au != INF);
	}
	//#endif

	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();

	return result_vector;
}

TEMPLATE// 181016 get segment from specific row with specific length
void TOOL::getFileStreamSegment(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, const int& const column_id, DataType*& const original_time_series) {
	assert(0);
#ifdef _DEBUG
	assert(column_id < g_time_series_length&& column_id != INF && g_time_series_length != INF && original_time_series_id != INF);
#endif

	string fs_row_string;
	string fs_row_number;

	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	int i = 0;
	int difference = 0;
	while (!file_stream.eof() && file_stream.is_open())
	{
		file_stream.good();
		file_stream.fail();

		//file_stream >> fs_row_string;
		getline(file_stream, fs_row_string);//different from above methods
		if (i == original_time_series_id) {
			stringstream sstr(fs_row_string);
			int string_id = 0;
			difference = 0;
			while (getline(sstr, fs_row_number, ',') && difference < g_time_series_length) {
				difference = string_id - column_id;
				if (difference >= 0 && difference < g_time_series_length) {
					original_time_series[difference] = stod(fs_row_number);
					assert(original_time_series[difference] != NULL);
				}
				string_id++;
			}
			break;
		}

		fs_row_string.clear();
		i++;
	}
	assert(original_time_series[int(g_time_series_length) - 1] != NULL);
	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();
}

TEMPLATE //For multiple dimension trajectories
void TOOL::getMultiFoldToSingleByID(string*& const multi_file_name, const int& arity_d, const double& single_series_length, const int& const original_time_series_id, DataType*& const original_time_series) {
	assert(0);
#ifdef _DEBUG
	assert(arity_d > 0 && single_series_length != INF && original_time_series_id != INF);
#endif

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream;
	fill_n(original_time_series, int(single_series_length * arity_d), INF);

	assert(file_stream);
	for (int file_id = 0; file_id < arity_d; file_id++) {
		int i = 0;
		file_stream = ifstream(multi_file_name[file_id]);
		while (!file_stream.eof() && file_stream.is_open())
		{
			file_stream.good();
			file_stream.fail();

			//file_stream >> fs_row_string;
			getline(file_stream, fs_row_string);//different from above methods
			if (i == original_time_series_id) {
				stringstream sstr(fs_row_string);
				int string_id = -1;
				while (getline(sstr, fs_row_number, ',') && string_id < single_series_length) {
					if (string_id > -1) {
						original_time_series[string_id + int(single_series_length * file_id)] = stod(fs_row_number);
						//assert(original_time_series[string_id] != NULL);
					}
					string_id++;
				}
				break;
			}

			fs_row_string.clear();
			i++;
		}
	}
	//assert(original_time_series[int(single_series_length) - 1] != NULL);
	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();
}


//191220 Original time series to instead vector, miss first data
TEMPLATE
template<typename T>
vector<T>& TOOL::getMultiFoldToSingleByID(string*& const multi_file_name, const int& arity_d, const double& single_series_length, const int& const original_time_series_id, vector<T>& const original_time_series_vector) {
	assert(0);
#ifdef _DEBUG
	assert(original_time_series_vector.empty() && arity_d > 0 && arity_d != INF && single_series_length != INF && single_series_length > 0 && original_time_series_id != INF && original_time_series_id > -1);
#endif

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream;
	//fill_n(original_time_series, int(single_series_length * arity_d), NULL);

	for (int file_id = 0; file_id < arity_d; file_id++) {
		int i = 0;
		file_stream = ifstream(multi_file_name[file_id]);
		assert(file_stream);
		while (!file_stream.eof() && file_stream.is_open())
		{
			file_stream.good();
			file_stream.fail();

			//file_stream >> fs_row_string;
			getline(file_stream, fs_row_string);//different from above methods
			if (i == original_time_series_id) {
				stringstream sstr(fs_row_string);
				int string_id = -1;
				while (getline(sstr, fs_row_number, ',') && string_id < single_series_length) {
					if (string_id > -1) {
						//original_time_series[string_id + int(single_series_length * file_id)] = stod(fs_row_number);
						original_time_series_vector.emplace_back(stod(fs_row_number));
						//assert(original_time_series[string_id] != NULL);
					}
					string_id++;
				}
				break;
			}

			fs_row_string.clear();
			i++;
		}
	}
	//assert(original_time_series[int(single_series_length) - 1] != NULL);
	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();
#ifdef _DEBUG
	assert(original_time_series_vector.size() == arity_d * single_series_length);
#endif
	return original_time_series_vector;
}

//************************************
// Method:getMultiFoldToSingleByID
// Qualifier:vector file addresss vector. //No miss first point value. Already normalized in advance, no normalize here
// Input:
// Output:
// date:210603
// author:
//************************************
TEMPLATE
template<typename T>
vector<T>& TOOL::getMultiFoldToSingleByID(const vector<string>& const multi_file_name, const int& arity_d, const double& single_series_length, const int& const original_time_series_id, vector<T>& const original_time_series_vector) {

	original_time_series_vector.clear();
	original_time_series_vector.shrink_to_fit();

#ifdef _DEBUG
	assert(multi_file_name.size() == arity_d && original_time_series_vector.empty() && arity_d > 0 && arity_d != INF && single_series_length != INF && single_series_length > 0 && original_time_series_id != INF && original_time_series_id > -1);
#endif

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream;
	//fill_n(original_time_series, int(single_series_length * arity_d), NULL);

	bool skip_first_data = false;

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& get split delimiter &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	ifstream file_stream_temp = ifstream(multi_file_name[0]);
	assert(file_stream_temp);
	char delimiter;
	while (file_stream_temp.get(delimiter) && delimiter != '\t' && delimiter != ',') {
	}
	file_stream_temp.close();
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	/*======================1 get original time series =====================*/
	for (int file_id = 0; file_id < arity_d; file_id++) {//dimension
		int i = 0;

		file_stream = ifstream(multi_file_name[file_id]);
		assert(file_stream);

		while (!file_stream.eof() && file_stream.is_open()) {// line string

			file_stream.good();
			file_stream.fail();

			//file_stream >> fs_row_string;
			//getline(file_stream, fs_row_string);//different from above methods
			getline(file_stream, fs_row_string, '\n');//201018 different from above methods

			if (i == original_time_series_id) {
				stringstream sstr(fs_row_string);

				int string_id = single_series_length + 1;// miss first point value

				if (skip_first_data) string_id = -1;
				else string_id = 0;

				while (getline(sstr, fs_row_number, delimiter) && string_id < single_series_length) {
					if (string_id > -1) {
						//original_time_series[string_id + int(single_series_length * file_id)] = stod(fs_row_number);
						original_time_series_vector.emplace_back(stod(fs_row_number));
						//assert(original_time_series[string_id] != NULL);
					}
					string_id++;
				}
				break;
			}

			fs_row_string.clear();
			i++;
		}

	}
	/*==============================================================*/


	//assert(original_time_series[int(single_series_length) - 1] != NULL);
	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();

#ifdef _DEBUG
	assert(original_time_series_vector.size() == arity_d * single_series_length);
	/*if (original_time_series_id == 5) {
		cout << "multi dimenson time series: ";
		for (auto&& au : original_time_series_vector) {
			cout << au << ",";
		}
		cout << endl;
	}*/
#endif

	for (auto&& au : original_time_series_vector) {
		assert(au != INF);
	}

	/*======================2 normalized =====================*/
	//TOOL::normalizeStandard(original_time_series_vector);
	/*cout << "noranlized multi dimenson time series: ";
	if (original_time_series_id == 5) {
		for (auto&& au : original_time_series_vector) {
			cout << au << ",";
		}
		cout << endl;
	}*/
	/*=====================================================*/

	return original_time_series_vector;
}

//210603
//TEMPLATE
//template<typename T>
//std::fstream& TOOL::GotoLine(std::fstream& file, const T& const line_num) {
//	file.seekg(std::ios::beg);
//	for (int i = 0; i < line_num - 1; ++i) {
//		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//	}
//	return file;
//}

//191220 Porject multi dataset to long single dataset
TEMPLATE
void TOOL::project_multi_data_to_single_data(string*& const multi_file_name, const int& const arity_d, const double& const single_series_number, const double& const single_series_length, const string& const write_file_name) {

#ifdef _DEBUG
	assert(multi_file_name != nullptr && arity_d != INF && arity_d > 0 && single_series_length != INF && single_series_length > 0 && single_series_length != INF && single_series_number != INF && single_series_number > 0);
#endif
	const int multi_time_series_length = arity_d * single_series_length;

	//DataType* original_time_series = new DataType[multi_time_series_length];
	vector<DataType> original_time_series_vector;
	for (int original_time_series_id = 0; original_time_series_id < single_series_number; original_time_series_id++) {
		assert(original_time_series_vector.empty());
		getMultiFoldToSingleByID(multi_file_name, arity_d, single_series_length, original_time_series_id, original_time_series_vector);
		assert(original_time_series_vector.size() == multi_time_series_length);
		writeResultNoCover(write_file_name, original_time_series_vector);
		original_time_series_vector.clear();
		original_time_series_vector.shrink_to_fit();
	}


	//delete[] original_time_series;
	//original_time_series = nullptr;
	original_time_series_vector.clear();
	original_time_series_vector.shrink_to_fit();
}


//191229 project long sigle dataset to multi dimension dataset
//************************************
// Method:project_multi_data_to_single_data
// Qualifier:vector file addresss vector
// Input:
// Output:
// date:191229
// author:
//************************************
TEMPLATE
void TOOL::project_single_data_to_multi_data_single(const string& const single_file_name, const int& const arity_d, const double& const single_series_number, const double& const single_series_length, const string& const write_multi_file_address) {
#ifdef _DEBUG
	assert(!single_file_name.empty() && arity_d != INF && arity_d > 0 && single_series_length != INF && single_series_length > 0 && single_series_number != INF && single_series_number > 0 && single_series_number > 0);
	assert(int(single_series_length) % int(arity_d) == 0);
#endif
	string write_name_coefficient[] = { "_X","_Y","_Z" };

	int multi_time_series_length = single_series_length / arity_d;
	vector<double> single_time_series_vector(single_series_length, INF);

	vector<double> multi_time_series_vector(multi_time_series_length + 1, INF);

	for (int point_id = 0; point_id < single_series_number; point_id++) {
		getFileStreamByID(single_file_name, single_series_length, point_id, single_time_series_vector);

		for (int dimension_id = 0; dimension_id < arity_d; dimension_id++) {
			int initial_id = dimension_id * multi_time_series_length;
			//int end_id = initial_id + multi_time_series_length;
			copy_n(single_time_series_vector.begin() + initial_id, multi_time_series_length, multi_time_series_vector.begin() + 1);
			string write_name = write_multi_file_address + write_name_coefficient[dimension_id];
			//writeResultNoCover(write_name, multi_time_series_vector);
			writeResultNoCoverNoSuffix(write_name, multi_time_series_vector);// no suffix "txt"
			fill_n(multi_time_series_vector.begin(), multi_time_series_vector.size(), INF);
		}
	}

	multi_time_series_vector.clear();
	multi_time_series_vector.shrink_to_fit();
	single_time_series_vector.clear();
	single_time_series_vector.shrink_to_fit();
}

//191229
TEMPLATE
void TOOL::project_single_data_to_multi_data_batch(const string& const single_file_address_file, const int& const file_number, const int& const arity_d, const double& const single_series_number, const double& const single_series_length, const string& const write_multi_file_address_file) {
#ifdef _DEBUG
	assert(!single_file_address_file.empty() && file_number != INF && file_number > 0 && arity_d != INF && arity_d > 0 && single_series_length != INF && single_series_length > 0 && single_series_number != INF && single_series_number > 0 && single_series_number > 0);
	assert(int(single_series_length) % int(arity_d) == 0);
#endif
	vector<string> read_file_address_vector;
	vector<string> write_file_address_vector;
	get_all_string_by_row(single_file_address_file, read_file_address_vector);
	get_all_string_by_row(write_multi_file_address_file, write_file_address_vector);

	for (int file_id = 0; file_id < file_number; file_id++) {
		project_single_data_to_multi_data_single(read_file_address_vector[file_id], arity_d, single_series_number, single_series_length, write_file_address_vector[file_id]);
	}
}

TEMPLATE
void TOOL::getFileStreamByRow(const string& file_name, const double& g_time_series_length, const int& const row_number, DataType** original_time_series) {
	assert(0);
#ifdef _DEBUG
	assert(g_time_series_length != INF && row_number != INF);
#endif

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(file_name);
	assert(file_stream);
	int row_id = 0;
	while (!file_stream.eof() && file_stream.is_open() && row_id < row_number) {
		file_stream.good();
		file_stream.fail();
		file_stream >> fs_row_string;

		stringstream sstr(fs_row_string);
		int string_id = -1;
		while (getline(sstr, fs_row_number, ',') && string_id < g_time_series_length) {
			if (string_id > -1) {
				original_time_series[row_id][string_id] = stod(fs_row_number);
				assert(original_time_series[row_id][string_id] != NULL);
			}
			string_id++;
		}

		fs_row_string.clear();
		row_id++;
	}
	assert(original_time_series[row_number - 1][int(g_time_series_length) - 1] != NULL);
	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();
}

//210603 miss first point arrray, normalizae
TEMPLATE
template<typename T>
void TOOL::getFileStreamByRow(const string& const file_name, const int& const g_time_series_length, const int& const row_number, T& const original_time_series) {
	//assert(0);
#ifdef _DEBUG
	assert(g_time_series_length != INF && row_number != INF);
#endif
	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& get split delimiter &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	ifstream file_stream_temp = ifstream(file_name);
	assert(file_stream_temp);
	char delimiter;
	while (file_stream_temp.get(delimiter) && delimiter != '\t' && delimiter != ',') {
	}
	file_stream_temp.close();
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	int row_id = 0;
	while (!file_stream.eof() && file_stream.is_open() && row_id < row_number) {
#ifdef _DEBUG
		assert(original_time_series[row_id].size() == g_time_series_length);
#endif
		file_stream.good();
		file_stream.fail();
		//file_stream >> fs_row_string;
		getline(file_stream, fs_row_string, '\n');

		stringstream sstr(fs_row_string);
		int string_id = -1;
		while (getline(sstr, fs_row_number, delimiter) && string_id < g_time_series_length) {
			if (string_id > -1) {
				original_time_series[row_id][string_id] = stod(fs_row_number);
				assert(original_time_series[row_id][string_id] != INF);
			}
			string_id++;
		}

		fs_row_string.clear();
		/*==============210603 normalized =====================*/
		TOOL::normalizeStandard(original_time_series[row_id]);
		/*=====================================================*/
		row_id++;
	}
	assert(original_time_series[row_number - 1][int(g_time_series_length) - 1] != INF);

#ifdef _DEBUG
	for (int row_id = 0; row_id < row_number; row_id++) {
		for (int array_id = 0; array_id < g_time_series_length; array_id++) {
			assert(original_time_series[row_id][array_id] != INF);
		}
	}
#endif

	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();
}

//191203 Not miss first point, normalized
// return g_time_series_length + 1.
TEMPLATE
template<typename T>
void TOOL::getFileStreamByRowNoMiss(const string& const file_name, const int& const g_time_series_length, const int& const row_number, T& const original_time_series) {

#ifdef _DEBUG
	assert(g_time_series_length != INF && row_number != INF);
#endif
	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& get split delimiter &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	ifstream file_stream_temp = ifstream(file_name);
	assert(file_stream_temp);
	char delimiter;
	while (file_stream_temp.get(delimiter) && delimiter != '\t' && delimiter != ',') {
	}
	file_stream_temp.close();
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	int row_id = 0;
	while (!file_stream.eof() && file_stream.is_open() && row_id < row_number) {
#ifdef _DEBUG
		assert(original_time_series[row_id].size() == g_time_series_length + 1);
#endif
		file_stream.good();
		file_stream.fail();
		//file_stream >> fs_row_string;
		getline(file_stream, fs_row_string, '\n');//201018 different from above methods

		stringstream sstr(fs_row_string);
		//int string_id = -1;
		int string_id = 0;
		while (getline(sstr, fs_row_number, delimiter) && string_id <= g_time_series_length) {// < to <=
			//if (string_id > -1) {
			original_time_series[row_id][string_id] = stod(fs_row_number);
			assert(original_time_series[row_id][string_id] != INF);
			//}

			string_id++;
		}

		fs_row_string.clear();
		/*======================210603 normalized =====================*/
		TOOL::normalizeStandard(original_time_series[row_id]);
		/*=====================================================*/
		row_id++;
	}
	assert(original_time_series[row_number - 1][int(g_time_series_length) - 1] != INF);

#ifdef _DEBUG
	for (int row_id = 0; row_id < row_number; row_id++) {
		for (int array_id = 0; array_id <= g_time_series_length; array_id++) {
			assert(original_time_series[row_id][array_id] != INF);
		}
	}
#endif

	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();
}

//210603 read datas from one, normalize, and write
//************************************
// Stuct:get_normalize_write_one_file_d
// Qualifier:get stream from one specific file and write to single file normalized, miss first point 210603
// date:210603
// Input: 1 list of file
//        2 id of file
//        3 number of row
//        4 time series length
//        5 write file name
// Output: mixed data file
// author:
//************************************
TEMPLATE
template<typename T, typename Y>
void TOOL::get_normalize_write_one_file(Y& const data_source, const T& const data_id) {
	TOOL::clearTXTFile_no_txt(data_source.read_file_address);//210603
	get_normalize_write_one_file(data_source.file_address_all, data_id, data_source.single_time_series_length, data_source.single_point_number, data_source.read_file_address);
}

//210603 read datas from one, normalize, and write
//************************************
// Stuct:get_normalize_write_one_file_d
// Qualifier:get stream from one specific file and write to single; file normalized; miss first point 210603
// date:210603
// Input: 1 list of file
//        2 id of file
//        3 number of row
//        4 time series length
//        5 write file name
// Output: mixed data file
// author:
//************************************
TEMPLATE
template<typename T, typename Y, typename T1, typename  T2, typename  T3>
void TOOL::get_normalize_write_one_file(const Y& const file_name, const T& const data_id, const T1& const g_time_series_length, const T2& const row_number, const T3& const write_file_name) {
	//DataType original_time_series[row_number][g_time_series_length];
	vector<vector<long double>> original_time_series(row_number, vector<long double>(g_time_series_length, INF));//210603 miss first point

	//not miss first point, so length is length +1, normalized 210603
	//getFileStreamByRowNoMiss(getStringByID(file_name, data_id), g_time_series_length, row_number, original_time_series);
	getFileStreamByRow(getStringByID(file_name, data_id), g_time_series_length, row_number, original_time_series);

	for (int row_id = 0; row_id < row_number; row_id++) {
		writeResultNoCover_no_txt(write_file_name, original_time_series[row_id]);
	}
}

//************************************
// Stuct:get_write_multi_files
// Qualifier:get stream from several files and write to single file normalized 210603
// date:191201
// Input: 1 list of file
//        2 number of file
//        3 number of row
//        4 time series length
//        5 write file name
// Output: mixed data file
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::get_write_multi_files(const string& const file_name, const vector<T>& const data_id_list, const double& const g_time_series_length, const int& const row_number, const string& const write_file_name) {
	//DataType original_time_series[row_number][g_time_series_length];
	vector<vector<long double>> original_time_series(row_number, vector<long double>(g_time_series_length + 1, INF));

	for (int file_id = 0; file_id < data_id_list.size(); file_id++) {
		//not miss first point, so length is length +1, normalized 210603
		getFileStreamByRowNoMiss(getStringByID(file_name, data_id_list[file_id]), g_time_series_length, row_number, original_time_series);
		for (int row_id = 0; row_id < row_number; row_id++) {
			writeResultNoCover_no_txt(write_file_name, original_time_series[row_id]);
		}
	}
}

TEMPLATE
void TOOL::getTXTStreamSpace(const string& const file_name, const int& const array_length, DataType*& const original_time_series) {
	ifstream infile;

	infile.open(file_name);
	assert(infile.is_open());

	int array_id = 0;
	double tmp = NULL;
	while (!infile.eof() && array_id < array_length) {
		infile >> tmp;
		original_time_series[array_id] = tmp;
		assert(original_time_series[array_id] != NULL);
		array_id++;
	}

	infile.close();
}

//************************************
// Method:getStringByID
// Qualifier:
// date:181207
// author:
//************************************
TEMPLATE
string& TOOL::getStringByID(const string& const file_name, const int& const row_num, string& const file_address) {//181218
	string item_name;
	ifstream nameFile;
	nameFile.open(file_name);
	assert(nameFile);
	string line;
	auto row_id = 0;
	while (std::getline(nameFile, file_address) && row_id < row_num)
	{
		row_id++;
	}
	return file_address;
}

//************************************
// Method:getStringByID
// Qualifier:
// date:191115
// author:
//************************************
TEMPLATE
string TOOL::getStringByID(const string& const file_name, const int& const row_num) {
	string file_address;
	ifstream nameFile;
	nameFile.open(file_name);

	assert(nameFile);

	auto row_id = 0;
	while (std::getline(nameFile, file_address) && row_id < row_num)
	{
		row_id++;
	}
	return file_address;
}

//************************************
// Method:getStringByID
// Qualifier:
// date:191115
// author:
//************************************
TEMPLATE
//191222 get all string and read it into vector
vector<string>& TOOL::get_all_string_by_row(const string& const file_name, vector<string>& const string_vector) {
	assert(string_vector.empty());
	string row_string;

	ifstream nameFile;
	nameFile.open(file_name);
	assert(nameFile);

	while (std::getline(nameFile, row_string)) {
		string_vector.emplace_back(row_string);
	}

	row_string.clear();
	row_string.shrink_to_fit();

	return string_vector;
}

//210824 get string and split by delimiter
TEMPLATE
template<typename T, typename Y, typename U>
vector<U>& TOOL::get_string_by_delimiter(const T& const file_name, const Y& const size_string, vector<U>& const string_vector) {

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream;

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& get split delimiter &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	ifstream file_stream_temp = ifstream(file_name);
	assert(file_stream_temp);
	char delimiter;
	while (file_stream_temp.get(delimiter) && delimiter != '\t' && delimiter != ',') {
	}
	file_stream_temp.close();
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	string_vector.clear();
	string_vector.shrink_to_fit();

	file_stream = ifstream(file_name);
	assert(file_stream);

	while (!file_stream.eof() && file_stream.is_open()) {// line string

		file_stream.good();
		file_stream.fail();

		//file_stream >> fs_row_string;
		//getline(file_stream, fs_row_string);//different from above methods
		getline(file_stream, fs_row_string, '\n');//201018 different from above methods
		stringstream sstr(fs_row_string);

		int string_id = 0;// miss first point value
		while (getline(sstr, fs_row_number, delimiter) && string_id < size_string) {
			
			string_vector.emplace_back(fs_row_number);
			string_id++;
		}
		break;
		fs_row_string.clear();
	}

	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();

	return string_vector;
}

//************************************
// Method:get_normalized_series_by_address_list
// Qualifier: get normalized original time series from ifstream
// Notice: missed first point
// Input: ifstream , point number, time seris length.
// Output: vector : normalized time sereis
// date:191115
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::get_normalized_series_by_stream(const TOOL::INPUT_ARGUMENT& const input_argument, ifstream& const file_stream, T*& const original_time_series) {
	assert(0);
	string row_string;
	string row_number;
	assert(input_argument.time_series_length != INF);
	file_stream >> row_string;
	//memory_account[2] = row_string.size();
	stringstream sstr(row_string);
	int array_id = -1;// jump the first data in file
	while (getline(sstr, row_number, ',') && array_id < input_argument.time_series_length) {
		if (array_id > -1) {
			original_time_series[array_id] = stod(row_number);
			assert(original_time_series[array_id] != INF);
			//multi_dimention_trajectories[array_id][string_id] = stod(fs_row_number);
			//cout << original_time_series[array_id] << ", ";
		}
		array_id++;
	}

#ifdef _DEBUG
	for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		//cout << original_time_series[array_id] << ",";
		assert(original_time_series[array_id] != INF);
	}
#endif

	normalizeStandard(input_argument.time_series_length, original_time_series);//z-score normalization

#ifdef _DEBUG
	for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		//cout << original_time_series[array_id] << ",";
		assert(original_time_series[array_id] != INF);
	}
#endif

	row_string.clear();
	row_string.shrink_to_fit();
	row_number.clear();
	row_number.shrink_to_fit();
}

//************************************
// Method:get_normalized_series_by_address_list
// Qualifier: get normalized original time series from txt file list according to file id
// Input: file id, point number, time seris length.
// Output: vector : normalized time sereis
// date:191115
// author:
//************************************
TEMPLATE
template<typename T>
const vector<T>& const TOOL::get_normalized_series_by_address_list(TOOL::INPUT_ARGUMENT& const input_argument, const int& const file_id, vector<T>& const normalized_series_vector) const {
#ifdef _DEBUG
	assert(input_argument.time_series_length != INF && input_argument.point_number != INF);
#endif
	normalized_series_vector.resize(input_argument.time_series_length, INF);
	string file_name;
	getStringByID(file_address, file_id, file_name);
	getFileStreamByID(file_name, input_argument.time_series_length, 10, normalized_series_vector);
	normalizeStandard(normalized_series_vector);

	return normalized_series_vector;
}

//************************************
// Method:vonvert_vector_to_string
// Qualifier: Write!! covnert vector list to string to write into txt file
// Input: file id, point number, time seris length.
// Output: string
// date:191201
// author:
//************************************
TEMPLATE
template<typename T>
string TOOL::convert_vector_to_string(const vector<T>& const vector_value) {
	assert(!vector_value.empty());

	std::ostringstream vector_to_string;
	// Convert all but the last element to avoid a trailing "," 
	std::copy(vector_value.begin(), vector_value.end() - 1, std::ostream_iterator<T>(vector_to_string, ","));
	// Now add the last element with no delimiter 
	vector_to_string << vector_value.back();

	return vector_to_string.str();
}

//************************************
// Method:writeInputArgument
// Qualifier:
// date:181207
// author:
//************************************
TEMPLATE
void TOOL::writeInputArgument(INPUT_ARGUMENT& input_argument, const string& write_file_name) {
	time_t now = time(0);// system time
	char dt[10];
	ctime_s(dt, sizeof dt, &now);
	//char* dt = ctime(&now);// from now to string
	ofstream outfile(write_file_name + ".txt", ios::app);
	assert(outfile.is_open());

	outfile << endl << dt << "read file:" << input_argument.read_file_name << "; write file: " << write_file_name << endl;;
	/*outfile << dimension_reduction_name[PAA_or_APCA] << " n = " << mg_file_time_series_length << ", N = " << mg_APCA_point_dimension << ", number = " << mg_d_index_point_number << ", K = " << K << ", MAXNODES = " << mg_max_nodes << endl;*/
	outfile << " n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension <<  input_argument.point_number << ", K = " << input_argument.K << endl;
	/*outfile << "  count apca point = " << g_n_account_apca_point << " times, p= " << g_n_account_apca_point / double(input_argument.point_number) << endl;
	outfile << "  APCA Memory = " << memory_account[0] + memory_account[1] + memory_account[2] + memory_account[3] << "  RTree Memeory = " << memory_account[4] + memory_account[5] + memory_account[6] + memory_account[7] + memory_account[8] << "  KNN Memory = " << memory_account[9] << endl;*/

	outfile.close();
}

TEMPLATE
void TOOL::writeKNNResult(TOOL::INPUT_ARGUMENT& input_argument) {
	time_t now = time(0);// system time
	char dt[26];
	ctime_s(dt, sizeof dt, &now);// from now to string
	//char* dt = ctime(&now);// from now to string
	//string file_name = *input_argument.write_file_name + ".txt";
	ofstream outfile(*input_argument.write_file_name + ".txt", ios::app);
	assert(outfile.is_open());

	outfile << endl << dt;
	/*outfile << dimension_reduction_name[PAA_or_APCA] << " n = " << mg_file_time_series_length << ", N = " << mg_APCA_point_dimension << ", number = " << mg_d_index_point_number << ", K = " << K << ", MAXNODES = " << mg_max_nodes << endl;*/
	outfile << " n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
	//outfile << "  count apca point = " << g_n_account_apca_point << " times, p= " << g_n_account_apca_point / double(input_argument.point_number) << endl;
	//outfile << "  APCA Memory = " << memory_account[0] + memory_account[1] + memory_account[2] + memory_account[3] << "  RTree Memeory = " << memory_account[4] + memory_account[5] + memory_account[6] + memory_account[7] + memory_account[8] << "  KNN Memory = " << memory_account[9] << endl;
	outfile.close();
}

TEMPLATE
template<typename T>
void TOOL::writeArray(const string& write_file_name, T*& const Array, const int& array_length) {
	//auto file_name = write_file_name + ".txt";
	ofstream outfile(write_file_name + ".txt", ios::out | ios::trunc);
	assert(outfile.is_open());
	for (int row_id = 0; row_id < array_length; row_id++)
		outfile << Array[row_id] << endl;
	outfile.close();
}

//191126 cover content without time stamp
TEMPLATE
template<typename T>
void TOOL::writeSingleResult(const string& write_file_name, T& result) {
	//time_t now = time(0);// system time
	//char* dt = ctime(&now);// from now to string
	ofstream outfile(write_file_name + ".txt", ios::trunc);
	assert(outfile.is_open());
	outfile << result << endl;
	outfile.close();
}

//191125 cover content with time stamp
TEMPLATE
template<typename T>
void TOOL::writeSingleResultTimeStamp(const string& write_file_name, T& result) {
	time_t now = time(0);// system time
	char dt[26];
	ctime_s(dt, sizeof dt, &now);// from now to string
	ofstream outfile(write_file_name + ".txt", ios::trunc);
	assert(outfile.is_open());

	outfile << dt << "  " << result << endl;
	outfile.close();
}

//210828 // with system time date. cover content
TEMPLATE
template<typename T, typename Y>
inline void TOOL::write_single_result_time_stamp_no_txt(const T& write_file_name, Y& result) {
	time_t now = time(0);// system time
	char dt[26];
	ctime_s(dt, sizeof dt, &now);// from now to string
	ofstream outfile(write_file_name, ios::trunc);
	assert(outfile.is_open());

	outfile << dt << "  " << result << endl;
	outfile.close();
}

//191126 Not cover content with time stamp
TEMPLATE
template<typename T>
void TOOL::writeSingleResultNoCoverTimeStamp(const string& write_file_name, T& result) {
	time_t now = time(0);// system time
	char* dt = ctime(&now);// from now to string
	ofstream outfile(write_file_name + ".txt", ios::app);
	assert(outfile.is_open());

	outfile << dt << "  " << result << endl;
	outfile.close();
}

TEMPLATE
template<typename T, typename Y>
void TOOL::writeSingleResult(Y& input_argument, const string& write_file_name, T& result) {//181207 Add input_argument
	writeInputArgument(input_argument, write_file_name);

	ofstream outfile(write_file_name + ".txt", ios::app);
	assert(outfile.is_open());
	outfile << result << endl;
	outfile.close();
}

//************************************
// Method:writeSingleResult
// Qualifier:
// Input: array
// Output: array in file
// date:181218
// author:
//************************************
TEMPLATE
template<typename T, typename Y>
void TOOL::writeSingleResult(const string& const write_file_name, T& array_data, const Y& const array_len) {//181218 Write array to file
		// Write data to a file.
	std::ofstream   outfile(write_file_name + ".txt");
	assert(outfile);

	//time_t now = time(0);// system time
	//char* dt = ctime(&now);// from now to string
	////ofstream outfile(write_file_name + ".txt", ios::app);
	//assert(outfile.is_open());
	//outfile << endl << dt;

	//std::copy(array, array + array_len,std::ostream_iterator<double>(outfile, " \n"));
	std::copy(array_data, array_data + array_len, std::ostream_iterator<double>(outfile, ","));
	//std::copy(array.begin(), array.end(), std::ostream_iterator<T>(outfile, ","));
	outfile.close();
}

//************************************
// Method:writeSingleResult
// Qualifier: Clear & Cover
// Input: vector
// Output: vector in file
// date:190501
// author:
//************************************
TEMPLATE
template<typename T>// will cover original content
void TOOL::writeSingleResult(const string& const write_file_name, vector<T>& array_data) {//190501 Write vector to file, Cover
		// Write data to a file.
	std::ofstream   outfile(write_file_name + ".txt");
	assert(outfile);
	std::copy(array_data.begin(), array_data.end(), std::ostream_iterator<T>(outfile, ","));
	outfile.close();
}

TEMPLATE
template<typename T, typename Y>
inline void TOOL::write_single_result_no_txt(const Y& const write_file_name, const T& const  result) {//210828
	ofstream outfile(write_file_name, ios::app);
	assert(outfile.is_open());
	outfile << result << endl;
	outfile.close();
}

TEMPLATE
template<typename T, typename Y>
inline void TOOL::write_single_result_no_txt(const Y& const write_file_name, const vector<T>& const array_data) {//210828
	std::ofstream   outfile(write_file_name);
	assert(outfile);
	std::copy(array_data.begin(), array_data.end(), std::ostream_iterator<T>(outfile, ","));
	outfile.close();
}

//************************************
// Method:writeResultNoCover
// Qualifier:No coverage. no need + "txt"
// Input: vector
// Output: vector in file, line by line
// date:190624
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::writeResultNoCover(const string& const write_file_name, vector<T>& array_data) {//190624 Write vector to file, no coverage

	//for (auto&& au : array_data) {
		//assert(au != INF);
	//}

	std::ofstream outfile(write_file_name + ".txt", ios_base::app);
	assert(outfile);
	std::copy(array_data.begin(), array_data.end(), std::ostream_iterator<T>(outfile, ","));
	outfile << endl;
	outfile.close();
}

//210603 no + "txt"
TEMPLATE
template<typename T, typename Y>
void TOOL::writeResultNoCover_no_txt(const Y& const write_file_name, vector<T>& array_data) {//Write vector to file, no coverage

	std::ofstream outfile(write_file_name, ios_base::app);
	assert(outfile);
	std::copy(array_data.begin(), array_data.end(), std::ostream_iterator<T>(outfile, ","));
	outfile << endl;
	outfile.close();
}

//No suffix like "txt"
//************************************
// Method:writeResultNoCoverNoSuffix
// Qualifier:No coverage. no need + "txt"
// Input: vector
// Output: vector in file, line by line
// date:191229
// author:
//************************************
TEMPLATE
template<typename T>
void  TOOL::writeResultNoCoverNoSuffix(const string& const write_file_name, vector<T>& array_data) {//190624 Write vector to file, no coverage

	for (auto&& au : array_data) {
		//assert(au != INF);
	}

	std::ofstream outfile(write_file_name, ios_base::app);
	assert(outfile);
	std::copy(array_data.begin(), array_data.end(), std::ostream_iterator<T>(outfile, ","));
	outfile << endl;
	outfile.close();
}

//************************************
// Method:writeDataTXTCover
// Qualifier: coverage
// Input:
// Output:
// date:190625
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::writeDataTXTCover(const string& const write_file_name, T& data) {//190625 Write T to file, coverage
	std::ofstream   outfile(write_file_name + ".txt");
	assert(outfile);
	outfile << data << endl;
	outfile.close();
}

//************************************
// Method:clearTXTFile
// Qualifier: clear content of txt file
// Input:
// Output:
// date:190625
// author:
//************************************
TEMPLATE
inline void TOOL::clearTXTFile(const string& const clear_file_name) {//190625
	std::ofstream  clear_file;
	clear_file.open(clear_file_name + ".txt", std::ofstream::out | std::ofstream::trunc);
	assert(clear_file);
	clear_file.close();
}

//210603
TEMPLATE
template<typename T>
inline void TOOL::clearTXTFile_no_txt(const T& const clear_file_name) {
	std::ofstream  clear_file;
	clear_file.open(clear_file_name, std::ofstream::out | std::ofstream::trunc);
	assert(clear_file);
	clear_file.close();
}

//initial time == 0
TEMPLATE
inline void TOOL::recordStartTime(TIME& time) {
	LARGE_INTEGER whole_first_f;    //timer frequency
	QueryPerformanceFrequency(&whole_first_f);
	time.dqFrequency = (double)whole_first_f.QuadPart;
	time.time_over.QuadPart = 0;
	time.time_start.QuadPart = 0;

	QueryPerformanceCounter(&time.time_start);
}

TEMPLATE
template<typename T>
T& TOOL::recordFinishTime(TIME& time, T& whole_first_run_time) {
	QueryPerformanceCounter(&time.time_over);    //Finish recording time  5
	whole_first_run_time = 1000000 * (time.time_over.QuadPart - time.time_start.QuadPart) / time.dqFrequency;
	//cout << whole_first_run_time << endl;
	time.~TIME();
	return whole_first_run_time;
}

TEMPLATE
inline long double TOOL::recordFinishTime(TIME& time) {
	QueryPerformanceCounter(&time.time_over);    //Finish recording time
	return 1000000.0 * (time.time_over.QuadPart - time.time_start.QuadPart) / time.dqFrequency;
}

TEMPLATE
template<typename T>
double& TOOL::distanceEUC(const T& const input_argument, DataType*& const Q, DataType*& const C, double& const distance_euc) {
	//cout << "\n\n distanceEUC()" << endl;
	assert(input_argument.time_series_length != NULL);
	int array_id = 0;
	double sum = 0;
	double difference = 0;
	for (array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		//cout << Q[i]<<" " << C[i]<<" " << Q[i] - C[i] << endl;
		difference = Q[array_id] - C[array_id];
		sum += difference * difference;
	}
	distance_euc = sqrt(sum);
	//cout << "distanceEUC  : " << distance << endl;
	return distance_euc;
}

//191129 Euclidean distance.
TEMPLATE
template<typename T>//For sum deviation
long double TOOL::distanceEUC(const vector<T>& const time_series_vector1, const vector<T>& const time_series_vector2) {
#ifdef _DEBUG
	assert(time_series_vector1.size() == time_series_vector2.size());//n is even
#endif
	long double euclidean_distance = INF;
	long double sum = 0;
	long double difference = INF;
	for (int segment_id = 0; segment_id < time_series_vector1.size(); segment_id++) {
		assert(time_series_vector1[segment_id] != INF && time_series_vector2[segment_id] != INF);
		difference = time_series_vector1[segment_id] - time_series_vector2[segment_id];
		sum += difference * difference;
		assert(difference != INF && sum != INF);
	}
	return euclidean_distance = sqrt(sum);
}

TEMPLATE
template<typename T, typename Y>//For sum deviation 181214
typename TOOL::OUTPUT_ARGUMENT& TOOL::getDeviation(const T& const input_argument, Y*& const original_time_series, Y*& const reconstruct_time_series, OUTPUT_ARGUMENT& const output_argument) {
	double sum = 0;
	double difference = NULL;
	output_argument.max_deviation = 0;

	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
	for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		difference = fabs(original_time_series[array_id] - reconstruct_time_series[array_id]);
		//cout <<"original time series: "<< original_time_series[array_id] <<" reconstruct: "<< reconstruct_time_series[array_id] << "  difffffff: " << difference << endl;
		output_argument.max_deviation = max(output_argument.max_deviation, difference);
		//deviation_sum += difference;
		sum += difference * difference;
	}

	output_argument.sum_deviation = sqrt(sum);

	assert(output_argument.sum_deviation >= 0 && output_argument.max_deviation >= 0);

	return output_argument;
}

TEMPLATE
template<typename T, typename Y>//For sum deviation 190501
typename TOOL::OUTPUT_ARGUMENT& TOOL::getDeviation(const T& const input_argument, Y*& const original_time_series, vector<Y>& const reconstruct_time_series, OUTPUT_ARGUMENT& const output_argument) {
	double sum = 0;
	double difference = NULL;
	output_argument.max_deviation = 0;
	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
	for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		difference = fabs(original_time_series[array_id] - reconstruct_time_series[array_id]);
		//cout <<"original time series: "<< original_time_series[array_id] <<" reconstruct: "<< reconstruct_time_series[array_id] << "  difffffff: " << difference << endl;
		output_argument.max_deviation = max(output_argument.max_deviation, difference);
		//deviation_sum += difference;
		sum += difference * difference;
	}

	output_argument.sum_deviation = sqrt(sum);

	assert(output_argument.sum_deviation > 0 && output_argument.max_deviation > 0);

	return output_argument;
}

//************************************
// Method:getDeviation
// Qualifier: For sum deviation
// Input:
// Output:
// date:200212
// author:
//************************************
TEMPLATE
template<typename T, typename Y>//For sum deviation 200212
typename TOOL::OUTPUT_ARGUMENT& TOOL::getDeviation(const T& const input_argument, const  vector<Y>& const original_time_series_vector, vector<Y>& const reconstruct_time_series_vector, OUTPUT_ARGUMENT& const output_argument) {
	double sum = 0;
	double difference = NULL;
	output_argument.max_deviation = 0;
	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF && input_argument.remainder != INF);
	for (int array_id = 0; array_id < input_argument.time_series_length; array_id++) {
		difference = fabs(original_time_series_vector[array_id] - reconstruct_time_series_vector[array_id]);
		//cout <<"original time series: "<< original_time_series[array_id] <<" reconstruct: "<< reconstruct_time_series[array_id] << "  difffffff: " << difference << endl;
		output_argument.max_deviation = max(output_argument.max_deviation, difference);
		//deviation_sum += difference;
		sum += difference * difference;
	}

	output_argument.sum_deviation = sqrt(sum);

	assert(output_argument.sum_deviation >= 0 && output_argument.max_deviation >= 0);

	return output_argument;
}

TEMPLATE
template<typename T, typename Y, typename U, typename T1>//For sum deviation & max deviation 210327
long double TOOL::getDeviation(const T& const input_argument, const vector<Y>& const original_time_series_vector, vector<Y>& const reconstruct_time_series_vector, const U& const segment_width, T1& const result_collection) {
	result_collection.sum_deviation = 0;
	result_collection.max_deviation = 0;
	result_collection.max_deviation_multiple_width = 0;

	int point_id = 0;
	long double difference = INF;
	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
		long double deviation_max = -INF;
		for (int interval_id = 0; interval_id < segment_width; interval_id++) {
			point_id = segment_width * segment_id + interval_id;// time series id
			difference = fabs(reconstruct_time_series_vector[point_id] - original_time_series_vector[point_id]);
			deviation_max = max(deviation_max, difference);
		}

		result_collection.max_deviation += deviation_max;
		result_collection.max_deviation_multiple_width += deviation_max * segment_width;
	}
}

TEMPLATE
template<typename Y>//For sum deviation 190501
typename TOOL::OUTPUT_ARGUMENT& TOOL::getDeviation(Y*& const original_time_series, vector<Y>& const reconstruct_time_series, const int& const time_series_length, OUTPUT_ARGUMENT& const output_argument) {
	double sum = 0;
	double difference = NULL;
	output_argument.max_deviation = 0;

	for (int array_id = 0; array_id < time_series_length; array_id++) {
		difference = fabs(original_time_series[array_id] - reconstruct_time_series[array_id]);
		//cout <<"original time series: "<< original_time_series[array_id] <<" reconstruct: "<< reconstruct_time_series[array_id] << "  difffffff: " << difference << endl;
		output_argument.max_deviation = max(output_argument.max_deviation, difference);
		//deviation_sum += difference;
		sum += difference * difference;
	}

	output_argument.sum_deviation = sqrt(sum);

	assert(output_argument.sum_deviation > 0 && output_argument.max_deviation > 0);

	return output_argument;
}

//**********************************************************************************************************************************
// Method:getDeviation
// Qualifier:get sum deviaiton of two vector
// Input:
// Output:
// date:191025
// author:
//**********************************************************************************************************************************
TEMPLATE
template<typename T>//For sum deviation 191025
double TOOL::getDeviation(const vector<T>& const original_time_series, const vector<T>& const reconstruct_time_series) {
	assert(original_time_series.size() == reconstruct_time_series.size());
	double sum_deviation = 0;
	for (int array_id = 0; array_id < original_time_series.size(); array_id++) {
		assert(original_time_series[array_id] != INF && reconstruct_time_series[array_id] != INF);
		sum_deviation += pow(original_time_series[array_id] - reconstruct_time_series[array_id], 2.0);
	}

	return sqrtl(sum_deviation);
}

TEMPLATE
template<typename T, typename Y>
void TOOL::SimpleBaseKNNSearch(const T& const input_argument, DataType*& const g_query_time_series, priority_queue<Y, vector<Y>, priorityDistanceEUC >& q_base_queue) {
	assert(0);
	cout << "===================================  Base squencial scan KNN  ===========================================\n";

	DataType* original_time_series = new DataType[int(input_argument.time_series_length)];
	//priority_queue<DataType, vector<DataType>, greater<DataType> > q_base_queue;
	double distance_euc = NULL;
	Y id_distance;

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(input_argument.read_file_name);
	assert(file_stream);

	int i = 0, j = NULL;
	while (!file_stream.eof() && i < input_argument.point_number) {
		file_stream >> fs_row_string;
		stringstream sstr(fs_row_string);
		j = -1;
		while (getline(sstr, fs_row_number, ',') && j < input_argument.time_series_length) {
			if (j > -1) {
				original_time_series[j] = stod(fs_row_number);
			}
			j++;
		}
		id_distance.id = i;
		id_distance.dist = distanceEUC(input_argument, g_query_time_series, original_time_series, distance_euc);
		q_base_queue.push(id_distance);
		i++;
	}

	file_stream.close();
	delete[] original_time_series;
	original_time_series = nullptr;


	/*cout << "Top K: \n";
	for (i = 0; i < K; i++) {
	cout << q_base_queue.top().d_dist << "\n";
	q_base_queue.pop();
	}
	cout << endl;*/
}

//**********************************************************************************************************************************
// Method:SimpleBaseKNNSearch
// Qualifier: get KNN by squencial scan  vector instead pointer
// Input: Normalized time series and query time series
// Output:Y is data struct
// date:191203
// author:
//**********************************************************************************************************************************
TEMPLATE
template<typename T, typename Y>
multiset<pair<double, int>>& TOOL::SimpleBaseKNNSearch(const T& const data_source, const vector<Y>& const query_time_series, multiset<pair<double, int>>& const knn_result_set) {
	assert(0);
	assert(query_time_series.size() == data_source.single_time_series_length && data_source.point_number != INF && data_source.single_time_series_length != INF);
	cout << "===================================  Base squencial scan KNN  ===========================================\n";
	cout << "File name: " << data_source.read_file_address << ", n: " << data_source.single_time_series_length << ", point number: " << data_source.point_number << endl;
	//priority_queue<Y, vector<Y>, priorityDistanceEUC> q_base_queue;
	//set<pair<size_t, double>> knn_result_set;
	vector<double> original_time_series_vector(int(data_source.single_time_series_length), INF);
	//priority_queue<DataType, vector<DataType>, greater<DataType> > q_base_queue;
	double distance_euc = INF;
	//Y id_distance;

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(data_source.read_file_address);
	assert(file_stream);

	int i = 0, j = NULL;
	while (!file_stream.eof() && i < data_source.point_number) {
		file_stream >> fs_row_string;
		stringstream sstr(fs_row_string);
		j = -1;
		while (getline(sstr, fs_row_number, ',') && j < data_source.single_time_series_length) {
			if (j > -1) {
				original_time_series_vector[j] = stod(fs_row_number);
				assert(original_time_series_vector[j] != INF);
			}
			j++;
		}
		TOOL::normalizeStandard(original_time_series_vector);

		distance_euc = distanceEUC(query_time_series, original_time_series_vector);
		knn_result_set.emplace(make_pair(distance_euc, i));
		assert(i != INF && distance_euc != INF && i >= 0 && distance_euc >= 0);

		//q_base_queue.push(id_distance);
		i++;
	}

	file_stream.close();
	original_time_series_vector.clear();
	original_time_series_vector.shrink_to_fit();

	return knn_result_set;
}


//191223 For multi dimension dataset. Vector instead pointer. Normalized time series
//**********************************************************************************************************************************
// Method:SimpleBaseKNNSearchMulti
// Qualifier: For multi dimension dataset get KNN by squencial scan  vector instead pointer
// Input: Normalized time series and query time series
// Output:Y is data struct
// date:191223
// author:
//**********************************************************************************************************************************
TEMPLATE
template<typename T, typename Y, typename U>
multiset<pair<U, int>>& TOOL::SimpleBaseKNNSearchMulti(const T& const data_source, const vector<Y>& const query_time_series_vector, multiset<pair<U, int>>& const knn_result_set) {
	assert(data_source.point_number != INF && data_source.point_number > 0 && data_source.time_series_dimension * data_source.single_time_series_length == query_time_series_vector.size() && data_source.multi_single_time_series_length == query_time_series_vector.size() && knn_result_set.empty());
	assert(data_source.read_file_address_vector.size() == data_source.time_series_dimension);
	vector<Y> normalized_time_series_vector;
	long double distance_euc = INF;
	for (int point_id = 0; point_id < data_source.point_number; point_id++) {
		/*=========================================get normalized data set==========================================================*/
		//already normalized
		TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, point_id, normalized_time_series_vector);
		assert(normalized_time_series_vector.size() == data_source.single_time_series_length * data_source.time_series_dimension && data_source.multi_single_time_series_length == normalized_time_series_vector.size());
		/*====================================================================================================================*/
		distance_euc = distanceEUC(query_time_series_vector, normalized_time_series_vector);
		knn_result_set.emplace(make_pair(distance_euc, point_id));
		assert(point_id != INF && distance_euc != INF && point_id >= 0 && distance_euc >= 0);
	}
	assert(!knn_result_set.empty() && knn_result_set.size() == data_source.point_number);
	normalized_time_series_vector.clear();
	normalized_time_series_vector.shrink_to_fit();
	return knn_result_set;
}


//**********************************************************************************************************************************
// Method:getRandomPoint
// Qualifier:
// Input:
// Output:
// date:181113
// author:
//**********************************************************************************************************************************
TEMPLATE
template<typename T> //Get random point for the array 181113
void TOOL::getRandomPoint(T*& const original_array, const int const& array_length, const int const& scale) {
	//printf("createRandomAPCAPoint()\n");

	srand((unsigned)time(NULL) * array_length);//seed

	for (int i = 0; i < array_length; i++) {
		//original_array[i] = rand() % 1001;
		original_array[i] = rand() % scale;
	}

	//TOOL::normalizeStandard(array_length, original_array);
	TOOL::printArray(original_array, array_length);
	//return getNormalArray(originalArray, pointNumber);
}

//**********************************************************************************************************************************
// Method:get_random_max
// Qualifier:get random value from [0, max_value - 1]
// Input: max value
// Output: randomw value
// date:191119
// author:
//**********************************************************************************************************************************
TEMPLATE
int TOOL::get_random_max(const int& const max_value) {
	srand(time(NULL));//seed
	return rand() % max_value;
}

//**********************************************************************************************************************************
// Method:get_random_max
// Qualifier:get random value from [0, max_value - 1]
// Input: length ,scale
// Output: randomw value
// date:200108
// author:
//**********************************************************************************************************************************
//200108
TEMPLATE
template<typename T>
void TOOL::get_random_vector(const int& const vector_size, const int& const scale, vector<T>& const time_series_vector) {
	srand((unsigned)time(NULL) * vector_size);//seed

	for (int i = 0; i < vector_size; i++) {
		//original_array[i] = rand() % 1001;
		time_series_vector.emplace_back(rand() % scale);
	}

	//TOOL::normalizeStandard(array_length, original_array);
	//print_vector(time_series_vector);
	//return getNormalArray(originalArray, pointNumber);
}

TEMPLATE
template<typename T> //Get mean vlaue between two value 181113
T& TOOL::getMeanValue(const T& const  vlaue_A, const T& const value_B, T& const  mean_value) {
	mean_value = (vlaue_A + value_B) / 2.0;

	return mean_value;
}

TEMPLATE
template<typename T>
void TOOL::normalizeA_B(const DataType& left_endpoint, const DataType right_endpoint, T*& const original_array, const int& array_length, T*& const normalized_array) {
	assert(left_endpoint < right_endpoint);
	T max_value = *max_element(original_array, original_array + array_length);
	T min_value = *min_element(original_array, original_array + array_length);

	cout << max_value << " " << min_value << endl;

	double  function_coefficient = (right_endpoint - left_endpoint) / (max_value - min_value);//(right-left)/(max-min)

	for (int array_id = 0; array_id < array_length; array_id++) {
		normalized_array[array_id] = function_coefficient * (original_array[array_id] - min_value) + T(left_endpoint);
	}

	/*printArray(normalized_array, array_length);*/
}

TEMPLATE
template<typename T>
double& TOOL::getAverage(T*& const original_array, const int& const array_length, double& average) {
	double sum = 0;

	for (int array_id = 0; array_id < array_length; array_id++) {
#ifdef _DEBUG
		assert(original_array[array_id] != INF);
#endif
		sum += original_array[array_id];
	}

	average = sum / double(array_length);
	//cout << "average: " << average << endl;
#ifdef _DEBUG
	assert(average != INF);
#endif
	return average;
}

//**********************************************************************************************************************************
// Method:getDeviation
// Qualifier: vector get sum deviaiton of two vector
// Input:
// Output:
// date:191115
// author:
//**********************************************************************************************************************************
TEMPLATE
template<typename T>
double& TOOL::getAverage(vector<T>& const original_array, const int& const array_length, double& average) {
#ifdef _DEBUG
	assert(original_array.size() == array_length);//n is even
	double test_sum = 0;
	double test_average = INF;
	for (int array_id = 0; array_id < original_array.size(); array_id++) {
		test_sum += original_array[array_id];
	}
	test_average = test_sum / double(original_array.size());
#endif

	average = std::accumulate(original_array.begin(), original_array.end(), 0.0) / double(original_array.size());

	//cout << "average: " << average << endl;
#ifdef _DEBUG
	assert(test_average == average);//n is even
#endif
	return average;
}

TEMPLATE
template<typename T>
double TOOL::getAverage(T*& const original_array) {//190307
	return std::accumulate(original_array.begin(), original_array.end(), 0.0) / original_array.size();
}

//**********************************************************************************************************************************
// Method:getVariance
// Qualifier: array get variance
// Input:
// Output:
// date:191115
// author:
//**********************************************************************************************************************************
TEMPLATE
template<typename T>
double& TOOL::getVariance(T*& const original_array, const int& const array_length, double& variance) {
	double sum = 0;
	double average = NULL;
	getAverage(original_array, array_length, average);

	for (int array_id = 0; array_id < array_length; array_id++) {
		sum += (original_array[array_id] - average) * (original_array[array_id] - average);
	}

	variance = sqrt(sum / double(array_length - 1));
	assert(variance != 0);
	//cout << "variance: " << variance << endl;
	return variance;
}

//**********************************************************************************************************************************
// Method:getVariance
// Qualifier: vector get variance
// Input:
// Output:
// date:191115
// author:
//**********************************************************************************************************************************
TEMPLATE
template<typename T>
double& TOOL::getVariance(vector<T>& const original_array, const int& const array_length, double& variance) {
#ifdef _DEBUG
	assert(original_array.size() == array_length);
#endif
	double sum = 0;
	double average = NULL;
	getAverage(original_array, original_array.size(), average);

	for (int array_id = 0; array_id < original_array.size(); array_id++) {
		sum += (original_array[array_id] - average) * (original_array[array_id] - average);
	}

	variance = sqrt(sum / double(original_array.size() - 1));

#ifdef _DEBUG
	assert(variance != 0 && variance != INF);
#endif
	//cout << "variance: " << variance << endl;
	return variance;
}

TEMPLATE//z-score normalization
template<typename T>
void TOOL::normalizeStandard(T*& const original_array, const int& const array_length, T*& const normalized_array) {
	double variance = NULL;
	double average = NULL;

	getAverage(original_array, array_length, average);

	getVariance(original_array, array_length, variance);

	for (int array_id = 0; array_id < array_length; array_id++) {
		normalized_array[array_id] = (original_array[array_id] - average) / variance;
	}

	/*printArray(normalized_array, array_length);*/
}

TEMPLATE//z-score normalization
template<typename T>
void TOOL::normalizeStandard(const int& const array_length, T*& const normalized_array) {

	T* copy_array = new T[array_length];
	double variance = NULL;
	double average = NULL;

#ifdef _DEBUG
	for (int id = 0; id < array_length; id++) {
		assert(normalized_array[id] != INF);
	}
#endif

	copy_n(normalized_array, array_length, copy_array);

	getAverage(copy_array, array_length, average);

	/*for (int id = 0; id < array_length; id++) {
		cout << copy_array[id] << ",";
	}*/

	getVariance(copy_array, array_length, variance);

	fill_n(normalized_array, array_length, NULL);

	for (int array_id = 0; array_id < array_length; array_id++) {
		normalized_array[array_id] = (copy_array[array_id] - average) / variance;

#ifdef _DEBUG
		assert(normalized_array[array_id] != INF);
#endif
	}

	delete[] copy_array;
	copy_array = nullptr;
}


//************************************
// Method:normalizeStandard
// Qualifier:
// Input:
// Output:
// date:191115
// author:
//************************************
TEMPLATE
template<typename T>
vector<T>& TOOL::normalizeStandard(vector<T>& const normalized_array) {

#ifdef _DEBUG
	assert(!normalized_array.empty());
#endif

	vector<T> copy_array(normalized_array.size(), INF);
	double variance = INF;
	double average = INF;

	copy_n(normalized_array.begin(), normalized_array.size(), copy_array.begin());

	getAverage(copy_array, normalized_array.size(), average);

	getVariance(copy_array, normalized_array.size(), variance);

	//fill_n(normalized_array, normalized_array.size(), INF);
#ifdef _DEBUG
	assert(average != INF && variance != INF);
#endif
	for (int array_id = 0; array_id < normalized_array.size(); array_id++) {
		normalized_array[array_id] = (copy_array[array_id] - average) / variance;
#ifdef _DEBUG
		assert(normalized_array[array_id] != INF && copy_array[array_id] != INF);
#endif
	}

	copy_array.clear();
	copy_array.shrink_to_fit();

	return normalized_array;
}

/*--  210906 --*/
TEMPLATE
template<typename T, typename Y, typename U>
inline void TOOL::get_equal_segment_width(const T& const n, const Y& const N, U& const equal_segment_width_struct) {
	equal_segment_width_struct.remainder = int(n) % int(N);//For PLA
	equal_segment_width_struct.segment_length_second = (n - equal_segment_width_struct.remainder) / N;
	equal_segment_width_struct.segment_length_first = equal_segment_width_struct.segment_length_second + 1;
}

//************************************
// Method:transfer_us_s
// Qualifier: transfer time us to s
// Input:us
// Output:minutes
// date:201104
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::transfer_us_s(vector<T>& const time_vector) {
	for (auto&& au : time_vector) {
		au /= 1000000;
	}
}

//************************************
// Method:transferSetToInt
// Qualifier:
// Input:
// Output:
// date:190623
// author:
//************************************
TEMPLATE
int TOOL::transferSetToInt(const set<int>& const random_endpoint_set) {//190623
	assert(!random_endpoint_set.empty());
	int integer_number = 0;
	for (auto&& au : random_endpoint_set) {
		integer_number += au;
		integer_number *= 10;
	}

	return integer_number;
}

//************************************
// Method:combinationUtil
// Qualifier:/* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */
// Input:
// Output:
// date:190623
// author:
//************************************
TEMPLATE
void TOOL::combinationUtil(int arr[], int n, int r, int index, int data[], int i) {//190623 /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */
	 // Current cobination is ready, print it
	if (index == r)
	{
		for (int j = 0; j < r; j++)
			cout << data[j] << " ";
		cout << endl;
		return;
	}

	// When no more elements are there to put in data[]
	if (i >= n)
		return;

	// current is included, put next at next location
	data[index] = arr[i];
	combinationUtil(arr, n, r, index + 1, data, i + 1);

	// current is excluded, replace it with next (Note that
	// i+1 is passed, but index is not changed)
	combinationUtil(arr, n, r, index, data, i + 1);
}

//************************************
// Method:printCombination
// Qualifier:The main function that prints all combinations of size r  in arr[] of size n. This function mainly uses combinationUtil()
// Input:
// Output:
// date:190623
// author:
//************************************
TEMPLATE
void TOOL::printCombination(int arr[], int n, int r) {//190623 The main function that prints all combinations of size r  in arr[] of size n. This function mainly uses combinationUtil()
	// A temporary array to store
	// all combination one by one
	int* data = new int[r];

	// Print all combination using
	// temprary array 'data[]'
	combinationUtil(arr, n, r, 0, data, 0);
}

//************************************
// Method:combinationUtilJump
// Qualifier:/* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */
// Input:No adjacent combination
// Output:
// date:190624
// author:
//************************************
TEMPLATE
void TOOL::combinationUtilJump(int arr[], int n, int r, int index, int data[], int i) {//190624 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */
	 // Current cobination is ready, print it
	if (index == r)
	{
		for (int j = 0; j < r; j++)
			cout << data[j] << " ";
		cout << endl;
		return;
	}

	// When no more elements are there to put in data[]
	if (i >= n)
		return;

	// current is included, put next at next location
	data[index] = arr[i];
	combinationUtilJump(arr, n, r, index + 1, data, i + 2);

	// current is excluded, replace it with next (Note that
	// i+1 is passed, but index is not changed)
	combinationUtilJump(arr, n, r, index, data, i + 1);
}

//************************************
// Method:printCombinationJump
// Qualifier:The main function that prints all combinations of size r  in arr[] of size n. This function mainly uses combinationUtil()
// Input:No adjacent combination
// Output:
// date:190624
// author:
//************************************
TEMPLATE
void TOOL::printCombinationJump(int arr[], int n, int r) {//190624 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()
	// A temporary array to store
	// all combination one by one
	int* data = new int[r];

	// Print all combination using
	// temprary array 'data[]'
	combinationUtilJump(arr, n, r, 0, data, 0);
}

//************************************
// Method:combinationUtilJumpVector
// Qualifier:The main function that prints all combinations of size r  in arr[] of size n. This function mainly uses combinationUtil()
// Input:No adjacent combination
// Output:
// date:190624
// author:
//************************************
TEMPLATE
void TOOL::combinationUtilJumpVector(const vector<int>& const arr, int index, vector<int> data, int i, set<pair<double, vector<int>>>& const endpoint_collection) {//190624 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */
	assert(!data.empty());
	// Current cobination is ready, print it
	if (index == data.size()) {
		endpoint_collection.insert(make_pair(data[index - 1], data));
		/*vector< APLA::AREA_COEFFICIENT>  area_vector;
		APLA::getRightEndpointAndWidth(data, area_vector);*/

		for (int j = 0; j < data.size(); j++)
			cout << data[j] << " ";
		cout << endl;

		return;
	}

	// When no more elements are there to put in data[]
	if (i >= arr.size())
		return;

	// current is included, put next at next location
	data[index] = arr[i];
	combinationUtilJumpVector(arr, index + 1, data, i + 2, endpoint_collection);

	// current is excluded, replace it with next (Note that
	// i+1 is passed, but index is not changed)
	combinationUtilJumpVector(arr, index, data, i + 1, endpoint_collection);
}

//************************************
// Method:printCombinationJumpVector
// Qualifier:The main function that prints all combinations of size r  in arr[] of size n. This function mainly uses combinationUtil()
// Input:No adjacent combination
// Output:
// date:190624
// author:
//************************************
TEMPLATE
void TOOL::printCombinationJumpVector(const vector<int>& const arr, const int& const r, set<pair<double, vector<int>>>& const endpoint_collection) {//190624 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()
	assert(!arr.empty() && r > 0 && arr.size() > r);

	// A temporary array to store
	// all combination one by one
	vector<int>  data;
	data.resize(r);

	// Print all combination using
	// temprary array 'data[]'
	combinationUtilJumpVector(arr, 0, data, 0, endpoint_collection);
}

//************************************
// Method:NChooseK
// Qualifier: Combination
// Input:
// Output:
// date:190624
// author:
//************************************
TEMPLATE
long long TOOL::NChooseK(int n, int k) {//190624 Combination
	if (k == 0) return 1;

	return (n * NChooseK(n - 1, k - 1)) / k;
}

//************************************
// Method:NChooseK
// Qualifier: Combination
// Input:
// Output:
// date:190624
// author:
//************************************
TEMPLATE
inline unsigned long long TOOL::n_choose_k(const unsigned long long& n, const unsigned long long& k) {//190624
	if (n < k) return 0;
	if (0 == n) return 0;
	if (0 == k) return 1;
	if (n == k) return 1;
	if (1 == k) return n;
	typedef unsigned long long value_type;
	value_type* table = new value_type[static_cast<std::size_t>(n * n)];
	std::fill_n(table, n * n, 0);
	class n_choose_k_impl
	{
	public:

		n_choose_k_impl(value_type* table, const value_type& dimension)
			: table_(table),
			dimension_(dimension)
		{}

		inline value_type& lookup(const value_type& n, const value_type& k)
		{
			return table_[dimension_ * n + k];
		}

		inline value_type compute(const value_type& n, const value_type& k)
		{
			if ((0 == k) || (k == n))
				return 1;
			value_type v1 = lookup(n - 1, k - 1);
			if (0 == v1)
				v1 = lookup(n - 1, k - 1) = compute(n - 1, k - 1);
			value_type v2 = lookup(n - 1, k);
			if (0 == v2)
				v2 = lookup(n - 1, k) = compute(n - 1, k);
			return v1 + v2;
		}

		value_type* table_;
		value_type dimension_;
	};
	value_type result = n_choose_k_impl(table, n).compute(n, k);
	delete[] table;
	return result;
}

//************************************
// Method:fact
// Qualifier: Combination
// Input:
// Output:// Returns factorial of n
// date:190624
// author:
//************************************
TEMPLATE
unsigned long long TOOL::fact(int n) {
	unsigned long long res = 1;
	for (int i = 2; i <= n; i++) {
		res *= i;
		assert(res > 0);
		//cout << res << endl;
	}

	return res;
}

//************************************
// Method:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!nCr Cannot deal with larg number
// Qualifier: Combination
// Input:
// Output:
// date:190624
// author:
//************************************
TEMPLATE
unsigned long long TOOL::nCr(int n, int r) {
	return fact(n) / (fact(r) * fact(n - r));
}

//************************************
// Method:nCr
// Qualifier: Combination
// Input:
// Output:
// date:190624
// author:
//************************************
TEMPLATE
unsigned long long TOOL::gcd(unsigned long long x, unsigned long long y) {//190624
	while (y != 0) {
		unsigned long long t = x % y;
		x = y;
		y = t;
	}
	return x;
}

//************************************
// Method:nCr
// Qualifier: Combination
// Input:
// Output:
// date:190624
// author:
//************************************
TEMPLATE
unsigned long long TOOL::N_choose_K(unsigned long long n, unsigned long long k) {//190624
	if (k > n)
		throw std::invalid_argument("invalid argument in choose");
	unsigned long long r = 1;
	for (unsigned long long d = 1; d <= k; ++d, --n)
	{
		unsigned long long g = gcd(r, d);
		r /= g;
		unsigned long long t = n / (d / g);
		if (r > (std::numeric_limits<unsigned long long>::max)() / t)
			throw std::overflow_error("overflow in choose");
		r *= t;
	}
	return r;
}

TEMPLATE
void TOOL::printOSPageSize() {
	SYSTEM_INFO si;
	GetSystemInfo(&si);

	printf("The page size for this system is %u bytes.\n", si.dwPageSize);//byte
	//cout << "boost: " << boost::interprocess::mapped_region::get_page_size() << endl;
}

TEMPLATE
template<typename T>
inline void TOOL::printInputArgument(const T& const input_argument) {
	cout << "Tree: " << input_argument.option_tree << " K : " << input_argument.K << " n : " << input_argument.time_series_length << " N : " << input_argument.point_dimension << endl;// << " read file name : " << input_argument.read_file_name << endl;
	//cout << "prunning power: " << input_argument.pruning_power << ", I/O cost: " << input_argument.IO_cost << endl;
}

TEMPLATE
template<typename T>
inline void TOOL::print_input_argument_KNN_result(const T& const input_argument) {

	assert(input_argument.pruning_power != INF && input_argument.knn_total_time != INF && input_argument.knn_total_time_has_IO != INF);

	cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;

	cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;

	cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;

	cout << "pruning power: " << input_argument.pruning_power << endl;
}

TEMPLATE
template<typename T>
void TOOL::printArray(T*& const test_array, const int& array_length) {
	assert(test_array != nullptr);

	for (int i = 0; i < array_length; i++) {
		cout << i << ": " << test_array[i] << ", ";
		//assert(test_array[i+1]- test_array[i]==1);
	}
	cout << endl;
}

TEMPLATE
template<typename T>
void TOOL::printArray(vector<T>& const test_array, const int& array_length) {//190501
	//assert(test_array != nullptr);

	for (int i = 0; i < array_length; i++) {
		cout << i << ": " << test_array[i] << ", ";
		//assert(test_array[i+1]- test_array[i]==1);
	}
	cout << endl;
}

//200108
TEMPLATE
template<typename T>
void TOOL::print_vector(vector<T>& const time_series_vector) {//190501
	cout << endl;
	for (auto&& au : time_series_vector) {
		cout << au << ",";
	}
	cout << endl;
}

//210122
TEMPLATE
template<typename T>
void TOOL::print_map(const T& const map_result) {
	cout << endl << "{  ";
	for (auto&& au : map_result) {
		cout << " ( " << au.first << " , " << au.second << " ) ; ";
	}
	cout << "   }" << endl;
}

//************************************
// Method:print_split_coefficients
// Qualifier:Print split id coeffifents
// date://200224 
// author:
//************************************
TEMPLATE
template<typename T>
inline void TOOL::print_split_coefficients(const vector<T>& const local_total_split_id_sum_deviation, const vector<T>& const local_total_split_id_shift, const vector<T>& const local_total_split_id_time, const vector<T>& const global_total_approximation_sum_deviation, const vector<T>& const global_total_approximation_time) {
	assert(local_total_split_id_sum_deviation.size() == local_total_split_id_shift.size() && local_total_split_id_time.size() == global_total_approximation_sum_deviation.size() && global_total_approximation_sum_deviation.size() == global_total_approximation_time.size());

	
	cout << endl;
	/*.....................................................*/
}

//200228 Print split id coeffifents
//************************************
// Method:print_split_coefficients
// Qualifier: print //local sum deviation, local shift, local time. global sum deviation global time, global prune power
// Input:
// Output:
// date:200228
// author:
//************************************
TEMPLATE
template<typename T>
inline void TOOL::print_split_coefficients(const vector<T>& const local_total_split_id_sum_deviation, const vector<T>& const local_total_split_id_shift, const vector<T>& const local_total_split_id_time, const vector<T>& const global_total_approximation_sum_deviation, const vector<T>& const global_total_approximation_time, const vector<T>& const global_total_knn_prune_power) {
	assert(local_total_split_id_sum_deviation.size() == local_total_split_id_shift.size() && local_total_split_id_time.size() == global_total_approximation_sum_deviation.size() && global_total_approximation_sum_deviation.size() == global_total_approximation_time.size() && global_total_approximation_time.size() == global_total_knn_prune_power.size());

	

}

//200326 initial_N coefficients
/*..................................................200326 Print initial_N evaluation.......................................................*/
//************************************
// Method:print_initial_N_coefficients
// Qualifier: Print initial N coefficients.
// Input:
// Output:
// date:191203
// author:
//************************************
TEMPLATE
template<typename T>
inline void TOOL::print_initial_N_coefficients(const vector<T>& const  approximation_initial_N_vector, const vector<T>& const  total_initial_N_prune_power_vector, const vector<T>& const  total_initial_N_sum_deviation_vector, const vector<T>& const total_initial_N_run_time_vector, const vector<T>& const total_initial_N_approximation_time_vector, const vector<T>& const total_initial_N_knn_time_vector) {
	assert(total_initial_N_prune_power_vector.size() == total_initial_N_sum_deviation_vector.size() && total_initial_N_sum_deviation_vector.size() == total_initial_N_run_time_vector.size() && total_initial_N_run_time_vector.size() == total_initial_N_approximation_time_vector.size() && total_initial_N_approximation_time_vector.size() == total_initial_N_knn_time_vector.size());

	
	cout << endl;
	/*.....................................................*/
}

/*..............................................................................................................................................*/

//200331 initial_N coefficients sort by N
//************************************
// Method:print_initial_N_sort_N_coefficients
// Qualifier: Print initial N coefficients sort by N.
// Input:
// Output:
// date:200331
// author:
//************************************
TEMPLATE
template<typename T, typename Y>
inline void TOOL::print_initial_N_sort_N_coefficients(const Y& const N_size, const vector<T>& const approximation_initial_N_vector, const vector<T>& const initial_N_by_N_prune_power_vector, const vector<T>& const initial_N_by_N_sum_deviation_vector, const vector<T>& const initial_N_by_N_run_time_vector, const vector<T>& const initial_N_by_N_approximation_time_vector, const vector<T>& const initial_N_by_N_knn_time_vector) {
	assert(N_size != INF && initial_N_by_N_prune_power_vector.size() == initial_N_by_N_sum_deviation_vector.size() && initial_N_by_N_sum_deviation_vector.size() == initial_N_by_N_run_time_vector.size() && initial_N_by_N_run_time_vector.size() == initial_N_by_N_approximation_time_vector.size() && initial_N_by_N_approximation_time_vector.size() == initial_N_by_N_knn_time_vector.size());

	/*.............200219 split id evaluation..............*/
	int N_id = 1;

	cout << "  <- Initial_N segment number ->\n";
	for (auto&& au : approximation_initial_N_vector) {
		cout << au << " | ";
		if (N_id % N_size == 0) cout << endl;
		N_id++;
	}
	cout << endl;
	N_id = 1;
	cout << "  <- Initial_N sort by N prune power ->\n";
	cout << "0 | 1 | 3 | 4 | 5\n";
	for (auto&& au : initial_N_by_N_prune_power_vector) {
		cout << au << ",  ";
		if (N_id % N_size == 0) cout << endl;
		N_id++;
	}
	cout << endl;
	N_id = 1;
	cout << "  <- Initial_N sort by N  sum deviation ->\n";
	cout << "0 | 1 | 3 | 4 | 5\n";
	for (auto&& au : initial_N_by_N_sum_deviation_vector) {
		cout << au << ",  ";
		if (N_id % N_size == 0) cout << endl;
		N_id++;
	}
	cout << endl;
	N_id = 1;
	cout << "  <- Initial_N sort by N  whole run time ->\n";
	cout << "0 | 1 | 3 | 4 | 5\n";
	for (auto&& au : initial_N_by_N_run_time_vector) {
		cout << au << ",  ";
		if (N_id % N_size == 0) cout << endl;
		N_id++;
	}
	cout << endl;
	N_id = 1;
	cout << "  <- Initial_N sort by N  approximation run time ->\n";
	cout << "0 | 1 | 3 | 4 | 5\n";
	for (auto&& au : initial_N_by_N_approximation_time_vector) {
		cout << au << ",  ";
		if (N_id % N_size == 0) cout << endl;
		N_id++;
	}
	cout << endl;
	N_id = 1;
	cout << "  <- Initial_N sort by N  knn run time ->\n";
	cout << "0 | 1 | 3 | 4 | 5\n";
	for (auto&& au : initial_N_by_N_knn_time_vector) {
		cout << au << ",  ";
		if (N_id % N_size == 0) cout << endl;
		N_id++;
	}
	cout << endl;
	/*.....................................................*/
}

//************************************
// Method:print_each_segment_right_endpoint
// Qualifier: print right endpoint of each segment in linked list or vector
// Input:
// Output:
// date:200710 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_right_endpoint(const T& const array_data_structure) {
	cout << "Each right endpoint: ";
	for (auto&& au : array_data_structure) {
		cout << au.right_endpoint << ", ";
	}
	cout << endl;
}

//************************************
// Method:print_each_segment_right_endpoint_with_order
// Qualifier: print right endpoint of each segment in linked list or vector with order
// Input:
// Output:
// date: 200817  
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_right_endpoint_with_order(const T& const array_data_structure) {
	cout << "each right endpoint: ";
	for (int i = 0; i < array_data_structure.size(); i++) {
		cout << i + 1 << ": " << array_data_structure[i].right_endpoint << " , ";
	}
	cout << endl;
}

//************************************
// Method:print_each_segment_segment_width
// Qualifier: print width of each segment in linked list or vector
// Input:
// Output:
// date:200710 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_width(const T& const array_data_structure) {
	cout << "each segment width: ";
	for (auto&& au : array_data_structure) {
		cout << au.rectangle_width << ", ";
	}
	cout << endl;
}

//************************************
// Method:print_each_segment_segment_width
// Qualifier: print width of each segment in linked list or vector
// Input:
// Output:
// date:200817 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_width_with_order(const T& const array_data_structure) {
	cout << "each segment width: ";
	for (int i = 0; i < array_data_structure.size(); i++) {
		cout << i + 1 << ": " << array_data_structure[i].rectangle_width << " , ";
	}
	cout << endl;
}

//************************************
// Method:print_each_segment_segment_ab
// Qualifier: print a&b of each segment in linked list or vector
// Input:
// Output:
// date:200710 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_ab(const T& const array_data_structure) {
	cout << "each segment a&b: ";
	for (auto&& au : array_data_structure) {
		cout << au.apla.a << ", " << au.apla.b << ";  ";
	}
	cout << endl;
}

//
//************************************
// Method:print_each_segment_ab_with_order
// Qualifier: print a&b of each segment in linked list or vector
// Input:
// Output:
// date:200817
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_ab_with_order(const T& const array_data_structure) {
	cout << "each segment a&b: ";
	for (int i = 0; i < array_data_structure.size(); i++) {
		cout << i + 1 << ": " << array_data_structure[i].apla.a << "," << array_data_structure[i].apla.b << ";  ";
	}
	cout << endl;
}

//************************************
// Method:print_each_segment_density
// Qualifier: print density of each segment in linked list or vector
// Input:
// Output:
// date:200814 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_density(const T& const array_data_structure) {
	cout << "Each segment density: ";
	for (auto&& au : array_data_structure) {
		if (au.right_subsegment)
			cout << 1.0 / au.right_subsegment->segment_density << ";  ";
	}
	cout << endl;
}

//200929 print area difference of each segment in linked list or vector
//************************************
// Method:print_each_segment_density
// Qualifier: print density of each segment in linked list or vector
// Input:
// Output:
// date:200814 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_area_difference(const T& const array_data_structure) {
	cout << "Each segment area difference: ";
	for (auto&& au : array_data_structure) {
		cout << au.area_difference << ";  ";
	}
	cout << endl;
}

//200929 print area difference of each segment in linked list or vector
//************************************
// Method:print_each_segment_sum_deviation
// Qualifier: print density of each segment in linked list or vector
// Input:
// Output:
// date:201007 
// author:
//************************************
//TEMPLATE
//template<typename T>
//void TOOL::print_each_segment_sum_deviation(const T& const array_data_structure) {
//	cout << "Each segment sum deviation: ";
//	for (auto&& au : array_data_structure) {
//		cout << au. << ";  ";
//	}
//	cout << endl;
//}

//************************************
// Method:print_each_segment_density_with_order
// Qualifier: print density of each segment in linked list or vector
// Input:
// Output:
// date:200817
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_density_with_order(const T& const array_data_structure) {
	cout << "Each segment density: ";
	for (int i = 0; i < array_data_structure.size(); i++) {
		cout << i + 1 << ": " << array_data_structure[i].segment_density << " , ";
	}
	cout << endl;
}

//************************************
// Method:print_each_segment_segment_ab
// Qualifier: print several coefficients of each segment in linked list or vector
// Input:
// Output:
// date:200712 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_coefficient(const T& const array_data_structure) {
	// right endpoint
	TOOL::print_each_segment_right_endpoint(array_data_structure);
	// width
	TOOL::print_each_segment_width(array_data_structure);
	// a&b
	TOOL::print_each_segment_ab(array_data_structure);
	// segment density
	TOOL::print_each_segment_density(array_data_structure);
	// segment density
	TOOL::print_each_segment_area_difference(array_data_structure);

	// segment sum deviaiton
	//TOOL::print_each_segment_sum_deviation(array_data_structure);
}

//************************************
// Method:print_each_segment_coefficient_with_order
// Qualifier: print several coefficients of each segment with order in linked list or vector
// Input:
// Output:
// date:200817 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_each_segment_coefficient_with_order(const T& const array_data_structure) {
	cout << "      *** Print coefficients with order: ***\n";
	// right endpoint
	TOOL::print_each_segment_right_endpoint_with_order(array_data_structure);
	// width
	TOOL::print_each_segment_width_with_order(array_data_structure);
	// a&b
	TOOL::print_each_segment_ab_with_order(array_data_structure);
	// segment density
	TOOL::print_each_segment_density_with_order(array_data_structure);
}


//************************************
// Method:print_string_vector
// Qualifier: print string title in vector
// Input:
// Output:
// date:2008131 
// author:
//************************************
TEMPLATE
template<typename T>
void TOOL::print_string_vector(const vector<T>& const string_vector) {
	for (auto&& au : string_vector) {
		cout << au << "   ";
	}
	cout << endl;
}

//200908 Print multimap second part iteration
TEMPLATE
template<typename T>
void TOOL::print_multimap_second(const T& const multi_map) const {

	cout << "Print multi-MAP second: ";

	for (auto&& au : multi_map) {
		cout << 1 / au.first << " , ";
	}

	cout << endl;
}

/*=============================================================================================================================================================*/

//************************************
// Method:initial_data_source
// Qualifier:initial Data Source structure.
// Input:
// Output:
// date:191203
// author:
//************************************
TEMPLATE
void TOOL::initial_data_source(DATA_SOURCE& const data_source) {

	/*.........................................................................................................................................................................................*/
	assert(data_source.read_file_address_vector.empty() && data_source.data_list.size() == data_source.single_file_number && data_source.data_type != INF && !data_source.data_list.empty() && data_source.single_file_number != INF && data_source.single_point_number != INF && data_source.single_time_series_length != INF);
	/*.........................................................................................................................................................................................*/

	switch (data_source.data_type) {
	case 1: {//single dataset

		//assert(data_source.file_name_pointer == nullptr);210603
		data_source.read_file_address = *data_source.file_name_pointer;//210603 Write normalized homogerous data to one specific data addresss
		TOOL::clearTXTFile_no_txt(data_source.read_file_address);//210603
		data_source.file_operation_number = data_source.single_file_number;
		data_source.point_number = data_source.single_point_number;
		data_source.multi_single_time_series_length = data_source.single_time_series_length;

		/* ||||||||||||||||||||||    Data sets Address, 0 has burst data 21, 1 no burst data 17, 2 all 24 datasets   ||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
		//input_argument.read_file_name = TOOL::getStringByID(TOOL::file_address, input_argument.file_id); Get all address.
		switch (data_source.option_has_burst_data) {
		case 0:// 0 has burst data 21
			/*........................................... Has burst time series .............................................................*/
			assert(data_source.single_file_number == 21);
			// Has burst time series: size = 21
			data_source.file_address_all = TOOL::file_address;//210603
			get_all_string_by_row(TOOL::file_address, data_source.file_address_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_number_address, 0, data_source.time_series_number_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_length_address, 0, data_source.time_series_length_vector);
			/*...............................................................................................................................*/
			break;
		case 1:// 1 no burst data 17
			/*........................................... No burst time series...............................................................*/
			assert(data_source.single_file_number == 17);
			// No burst time series : size = 17.No file 6, 17, 20 and 22
			data_source.file_address_all = TOOL::no_burst_single_file_address;//210603
			get_all_string_by_row(TOOL::no_burst_single_file_address, data_source.file_address_vector);
			TOOL::getFileStreamByID0Vector(TOOL::no_burst_single_file_number_address, 0, data_source.time_series_number_vector);
			TOOL::getFileStreamByID0Vector(TOOL::no_burst_single_file_length_address, 0, data_source.time_series_length_vector);
			/*...............................................................................................................................*/
		case 2:// Has all time series: size = 24
			/*........................................... Has all time series: size = 24 ........................................................*/
			assert(data_source.single_file_number == 24);
			data_source.file_address_all = TOOL::file_address_24;//210603
			get_all_string_by_row(TOOL::file_address_24, data_source.file_address_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_number_address_24, 0, data_source.time_series_number_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_length_address_24, 0, data_source.time_series_length_vector);
			/*...................................................................................................................................*/
			break;
		case 3:// 201016 UCR2018 512 41 datasets;
			/*........................................... Has all time series: size = 41 ........................................................*/
			assert(data_source.single_file_number == 41);
			data_source.file_address_all = TOOL::file_address_UCR2018_512_41;//210603
			get_all_string_by_row(TOOL::file_address_UCR2018_512_41, data_source.file_address_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_number_address_UCR2018_512_41, 0, data_source.time_series_number_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_length_address_UCR2018_512_41, 0, data_source.time_series_length_vector);
			/*...................................................................................................................................*/
			break;
		case 4:// 201019 UCR2018 1024 21 datasets;
			/*=================================    Has all time series: size = 21, length = 1024   ===============================================*/
			assert(data_source.single_file_number == 21);
			data_source.file_address_all = TOOL::file_address_UCR2018_1024_21;//210603
			get_all_string_by_row(TOOL::file_address_UCR2018_1024_21, data_source.file_address_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_number_address_UCR2018_1024_21, 0, data_source.time_series_number_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_length_address_UCR2018_1024_21, 0, data_source.time_series_length_vector);
			/*====================================================================================================================================*/
			break;
		case 5:// 201020 UCR2018 choose specific datasets; get all 41 files
			/*==============================    get all time series: size = Adaptive, length = Adaptive   =========================================*/
			data_source.file_address_all = TOOL::file_address_UCR2018_512_41;//210603
			get_all_string_by_row(TOOL::file_address_UCR2018_512_41, data_source.file_address_vector);// get all fiel address
			TOOL::getFileStreamByID0Vector(TOOL::single_file_number_address_UCR2018_512_41, 0, data_source.time_series_number_vector);// time series number
			TOOL::getFileStreamByID0Vector(TOOL::single_file_length_address_UCR2018_512_41, 0, data_source.time_series_length_vector);// time series length
			/*=====================================================================================================================================*/
			break;
		case 6:// 210305 just for each step in SAPLA. UCR2018 21 files; get all 41 files
			/*==============================    get all time series: size = Adaptive, length = Adaptive   =========================================*/
			data_source.file_address_all = TOOL::file_address_UCR2018_512_41;//210603
			get_all_string_by_row(TOOL::file_address_UCR2018_512_41, data_source.file_address_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_number_address_UCR2018_512_41, 0, data_source.time_series_number_vector);
			TOOL::getFileStreamByID0Vector(TOOL::single_file_length_address_UCR2018_512_41, 0, data_source.time_series_length_vector);
			/*=====================================================================================================================================*/
			break;
		default:
			assert(0);
		}

		/*........................................................................................................................................*/
		//assert(data_source.file_address_vector.size() == data_source.file_operation_number && data_source.time_series_number_vector.size() == data_source.file_operation_number);
		assert(data_source.time_series_number_vector.size() == data_source.time_series_length_vector.size());
		/*........................................................................................................................................*/
		//file_name_vector
		/*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

		break;
	}
	case 2: {//mixed dataset
		//assert(mixed_data_file_name != nullptr);
		data_source.file_operation_number = 1;
		data_source.point_number = data_source.single_file_number * data_source.single_point_number;
		data_source.multi_single_time_series_length = data_source.single_time_series_length;
		data_source.read_file_address = *data_source.file_name_pointer;// Write homogerous data to heterogerous data addresss
		data_source.time_series_number_vector.emplace_back(data_source.point_number);
		data_source.time_series_length_vector.emplace_back(data_source.multi_single_time_series_length);
		data_source.file_address_all = TOOL::file_address_UCR2018_1024_21;//210603
		/*-------------------------------------------Rewrite mixed data file--------------------------------------------------------*/
		TOOL::clearTXTFile_no_txt(data_source.read_file_address);//191201 The Write destination of heterogerous data address
		//skip first element// file_address_UCR2018_1024_15
		TOOL::get_write_multi_files(TOOL::file_address_UCR2018_1024_21, data_source.data_list, data_source.single_time_series_length, data_source.single_point_number, data_source.read_file_address);//191201
		//data_source.read_file_address = data_source.read_file_address + ".txt";//210603
		assert(data_source.file_address_vector.empty() && !data_source.read_file_address.empty());
		data_source.file_address_vector.emplace_back(data_source.read_file_address);
		/*---------------------------------------------------------------------------------------------------------------------------*/
		break;
	}
	case 3: {//multi dimension data set
		assert(data_source.multi_single_time_series_length == INF);
		assert(data_source.time_series_dimension != INF && data_source.single_file_number % data_source.time_series_dimension == 0 && data_source.single_time_series_length != INF && data_source.file_name_pointer == nullptr);
		data_source.file_operation_number = data_source.single_file_number / data_source.time_series_dimension;
		data_source.point_number = data_source.single_point_number;
		//data_source.read_file_address = *data_source.file_name_pointer;
		data_source.multi_single_time_series_length = data_source.single_time_series_length * data_source.time_series_dimension;
		/*-------------------------------------------get all address of three dimension data set--------------------------------------------------------*/
		data_source.file_address_all = TOOL::three_d_file_address;//210603
		get_all_string_by_row(TOOL::three_d_file_address, data_source.file_address_vector);
		TOOL::getFileStreamByID0Vector(TOOL::three_d_file_number_address, 0, data_source.time_series_number_vector);
		TOOL::getFileStreamByID0Vector(TOOL::three_d_file_length_address, 0, data_source.time_series_length_vector);
		assert(data_source.file_address_vector.size() == data_source.single_file_number && data_source.file_address_vector.size() % data_source.time_series_dimension == 0 && data_source.time_series_number_vector.size() == data_source.time_series_length_vector.size());
		/*----------------------------------------------------------------------------------------------------------------------------------------------*/
		break;
	}
	default:
		assert(0);
		break;
	}

}

//************************************
// Method:get_read_multi_file_address
// Qualifier:get multi dimension file address,  address size == data dimension
// Input:
// Output:
// date:191223
// author:
//************************************
TEMPLATE
vector<string>& TOOL::get_read_multi_file_address(DATA_SOURCE& const data_source, const int& const file_id) {
	assert(data_source.time_series_dimension != INF && file_id != INF && data_source.read_file_address_vector.empty() && data_source.data_type > 0);

	switch (data_source.data_type) {
	case 0: // 0 manipulate data set;    
		assert(0);
		break;
	case 1: {// 1 real data set, homogenous(single) data;
		assert(data_source.file_address_vector.size() >= data_source.single_file_number);

		/*----------- 210603 no normalize ----------*/
		//data_source.read_file_address = data_source.file_address_vector[file_id];//210603
		//data_source.read_file_address_vector.emplace_back(data_source.file_address_vector[file_id]); 210603
		////assert(TOOL::getStringByID(TOOL::no_burst_single_file_address, file_id) == data_source.read_file_address && data_source.multi_single_time_series_length == data_source.single_time_series_length);
		//assert((TOOL::getStringByID(TOOL::file_address, file_id) == data_source.read_file_address || TOOL::getStringByID(TOOL::no_burst_single_file_address, file_id) == data_source.read_file_address || TOOL::getStringByID(TOOL::file_address_24, file_id) == data_source.read_file_address || TOOL::getStringByID(TOOL::file_address_UCR2018_512_41, file_id) == data_source.read_file_address) && data_source.multi_single_time_series_length == data_source.single_time_series_length);
		/*------------------------------------*/

		/*----------- 210603 normalize ----------*/
		data_source.read_file_address_vector.emplace_back(data_source.read_file_address);//210603
		/*------------------------------------*/
		break;
	}
	case 2: {// 2 real data set, heterogeneous(mixed) data.
		assert(data_source.file_name_pointer != nullptr && !data_source.read_file_address.empty() && data_source.multi_single_time_series_length == data_source.single_time_series_length);
		//data_source.read_file_address = *data_source.file_name_pointer;
		data_source.read_file_address_vector.emplace_back(data_source.read_file_address);
		break;
	}
	case 3: {// 3 multi dimension homogenous(single) data.
		assert(data_source.file_operation_number * data_source.time_series_dimension == data_source.file_address_vector.size() && data_source.file_address_vector.size() == data_source.single_file_number);
		int initial_file_id = file_id * data_source.time_series_dimension;
		for (int dimension_id = 0; dimension_id < data_source.time_series_dimension; dimension_id++) {
			data_source.read_file_address_vector.emplace_back(data_source.file_address_vector[initial_file_id + dimension_id]);
		}
		assert(data_source.multi_single_time_series_length == data_source.time_series_dimension * data_source.single_time_series_length && data_source.read_file_address_vector.size() == data_source.time_series_dimension);
		break;
	}
	case 4:// 4 multi mixed data
		assert(0);
		break;
	case 5://APCA
		assert(0);
		break;
	default:
		assert(0);
		break;
	}

	assert(data_source.read_file_address_vector.size() == data_source.time_series_dimension);
	return data_source.read_file_address_vector;
}

//191223 read and normalized multi or single dimenson time series from file by time series id
//************************************
// Method:read_normalized_multi_time_series
// Qualifier:get multi dimension file address,  address size == data dimension
// Input:
// Output:
// date:210603
// author:
//************************************
TEMPLATE
template<typename T>
vector<T>& TOOL::read_normalized_multi_time_series(const DATA_SOURCE& const data_source, const int& const time_series_id, vector<T>& const normalized_time_series_vector) {
	normalized_time_series_vector.clear();
	normalized_time_series_vector.shrink_to_fit();

	/*...............................................................................................................................................*/
#ifdef _DEBUG
	assert(time_series_id != INF && normalized_time_series_vector.empty() && data_source.time_series_dimension != INF && data_source.single_time_series_length != INF && !data_source.read_file_address_vector.empty());
#endif
	/*...............................................................................................................................................*/

	/*###################################################################################################*/
	//already normalized //get data from new file. already normalized in advanced, no miss first point.
	TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, time_series_id, normalized_time_series_vector);
	/*###################################################################################################*/

	/*............................................................................................................*/
#ifdef _DEBUG
	bool miss_first_piont = false;
	DataType* original_time_series = new DataType[data_source.multi_single_time_series_length];
	vector<DataType> original_time_series_vector_evaluation(normalized_time_series_vector.size(), INF);

	assert(normalized_time_series_vector.size() == data_source.single_time_series_length * data_source.time_series_dimension && data_source.multi_single_time_series_length == normalized_time_series_vector.size());
	switch (data_source.data_type) {
	case 0: // 0 manipulate data set;    
		assert(0);
		break;
	case 1: {// 1 real data set, homogenous(single) data;
		//data_source.read_file_address = data_source.file_address_vector[file_id];
		//assert(TOOL::getStringByID(TOOL::file_address, file_id) == data_source.read_file_address && data_source.multi_single_time_series_length == data_source.single_time_series_length);
		if (miss_first_piont)
			TOOL::getFileStreamByID(data_source.read_file_address, data_source.multi_single_time_series_length, time_series_id, original_time_series);
		else
			TOOL::getFileStreamByID0Vector(data_source.read_file_address, time_series_id, original_time_series_vector_evaluation);

		//TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);
		break;
	}
	case 2: {// 2 real data set, heterogeneous(mixed) data.
		//assert(!data_source.read_file_address.empty() && data_source.multi_single_time_series_length == data_source.single_time_series_length);
		if (miss_first_piont)
			TOOL::getFileStreamByID(data_source.read_file_address, data_source.multi_single_time_series_length, time_series_id, original_time_series);
		else
			TOOL::getFileStreamByID0Vector(data_source.read_file_address, time_series_id, original_time_series_vector_evaluation);
		//TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);
		break;
	}
	case 3: {// 3 multi dimension homogenous(single) data.
		assert(data_source.read_file_address_vector.size() == data_source.time_series_dimension && data_source.file_operation_number * data_source.time_series_dimension == data_source.file_address_vector.size() && data_source.time_series_dimension != INF);
		normalized_time_series_vector.clear();
		normalized_time_series_vector.shrink_to_fit();
		//already normalized in this function

		TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, time_series_id, normalized_time_series_vector);

		assert(normalized_time_series_vector.size() == data_source.multi_single_time_series_length && data_source.multi_single_time_series_length == data_source.time_series_dimension * data_source.single_time_series_length);
		copy_n(normalized_time_series_vector.begin(), normalized_time_series_vector.size(), original_time_series);
		break;
	}
	case 4:// 4 multi mixed data
		assert(0);
		break;
	case 5://APCA
		assert(0);
		break;
	default:
		assert(0);
		break;
	}

	if (miss_first_piont)
		TOOL::normalizeStandard(data_source.multi_single_time_series_length, original_time_series);
	else
		TOOL::normalizeStandard(original_time_series_vector_evaluation);
	/*if (m_temp_queue_top.original_time_series_id == 5) {
		cout << "normalized original time series: \n";
		for (int array_id = 0; array_id < data_source.multi_single_time_series_length; array_id++) {
			cout << original_time_series[array_id] << ",";
		}
		cout << endl;
	}*/
	for (int array_id = 0; array_id < normalized_time_series_vector.size(); array_id++) {
		//cout << "================================\n";
		//cout << original_time_series_vector[array_id] << endl;
		//cout << original_time_series[array_id] << endl;
		if (miss_first_piont)
			assert(float(normalized_time_series_vector[array_id]) == float(original_time_series[array_id]));
		else
			assert(fabs(float(normalized_time_series_vector[array_id]) - float(original_time_series_vector_evaluation[array_id])) < 0.01);
	}

	delete[] original_time_series;
	original_time_series = nullptr;
	original_time_series_vector_evaluation.clear();
	original_time_series_vector_evaluation.shrink_to_fit();
#endif
	/*....................................................................................................................*/

	return normalized_time_series_vector;
}


//210823 two vectors have same value and size
TEMPLATE
template<typename T, typename Y>
bool TOOL::assert_same_vector(const T& const vector_0, const Y& const vector_1) {
	assert(vector_0.size() == vector_1.size());
	for (int id_point = 0; id_point < vector_0.size(); id_point++) {
		assert(vector_0[id_point] == vector_1[id_point]);
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**************************************************************************************************************************************************************/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif


