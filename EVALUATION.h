#pragma once

#ifndef EVALUATION_H
#define EVALUATION_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "APLA_ICDE07.h"
#include "CAPLA.h"
#include "MULTI_DIMENSION.h"

TEMPLATE
class EVALUATION : virtual public TOOL, virtual public APLA, virtual public CHEBYSHEV_QUAL{
public:

	//210822 N, K, name
	template<typename T>
	struct RESULT_COEFFICIENT;

	struct ADDRESS_STRING;//210822 Address for writing

	template<typename T>
	struct METHOD_RESULT;//210820 for each method

	template<typename T>
	struct METHOD_FILE_RESULT;//210820

	template<typename T>
	struct TOTAL_RESULT;//210820 for each bar

	template<typename T, typename Y>
	inline void initial_address(const T& const suffix, Y& const Address_String);//210823 for each bar

	template<typename T, typename Y>
	inline void initial_address_no_txt(const T& const suffix, Y& const Address_String);//210823 for each bar

	template<typename T, typename Y, typename U>
    void assign_pointer_vector(const T& const size_assign, const Y& const value_assign, U& const pointer_vector);//210827 for pointer vector

	//211219 for clearing pointer vector
	template<typename T>
	void clear_pointer_vector(T& const pointer_vector);

	//191106 evaluate ICDE07
	void evaluate_ICDE07();

	//191106 evaluate approxiamtion of APLA PLA, APCA, PAA & ICDE07_APLA
	void evaluate_approximation();

	//191106 evaluate KNN search for multi & single time series of APLA PLA, APCA, PAA & ICDE07_APLA
	void evaluate_multi_KNN();

	//200101 change the order of loop to speed up algorithm
	template<typename T>
	void evaluate_multi_KNN_speed(const T& const evaluation_argument_struct);

	template<typename T, typename Y, typename U>
	inline void print_input_output(const T& const input_argument, const Y& const output_argument, const U& const linear_scan_struct);
};


//#include "EVALUATION.cpp"

//210822 Address for writing
TEMPLATE
struct Evaluation::ADDRESS_STRING {

	vector<string*> address_string_coefficient_vector;

	vector<string*> address_method_pointer_vector;

	vector<string*> address_file_method_pointer_vector;

	vector<string*> address_total_pointer_vector;

	/*-------------------------------------------------             String Name           --------------------------------------------------------*/
	// 1K_total_number , 2N_total_number , 3file_total_number, 4method_total_number, 5point number, 6query id, 7max id, 8 single time series length ,9 data type, 10 evaluation type, 11 split methods number, 12 option homogenous data type
	string  str_1 = "./191120KNNMatlab/BarChartCoefficient"; //
	string  str_25 = "./191120KNNMatlab/NameRepresentationMethod";//200831 the name of representation method, string vector.
	string  str_99 = "./191120KNNMatlab/IDQuery";//211217 the ID of Query time series vector.
	string  str_100 = "./191120KNNMatlab/IDRepresentationMethod";//210902 the ID of representation method, number vector.
	string  str_101 = "./191120KNNMatlab/OptionBarchartTree";//210902 the option of Tree for representation method.
	string  str_26 = "./191120KNNMatlab/IDFileChoosed";//201020 Choose File id: file_id_choose_vector
	string  str_4 = "./191120KNNMatlab/KCoefficient";//191122 K records. eg. K = 2, 4 , 8
	string  str_11 = "./191120KNNMatlab/NCoefficient";//191128 N records. eg. N = 18, 36, 72, 144, 288
	string  str_5 = "./191120KNNMatlab/ReadMe";//191122
	string  str_10 = "./191120KNNMatlab/StringSuffix";//191126 prefix of file.

	string  str_6 = "./191120KNNMatlab/MethodPrunePower";//191123 size is total method number. prune power of every method for matlab barchart
	string  str_35 = "./191120KNNMatlab/MethodPrunePowerCombine";//201221  prune power combine. prune power / accuracy
	string  str_41 = "./191120KNNMatlab/MethodMaxDeviation";//210327
	string  str_42 = "./191120KNNMatlab/MethodMaxWidthDeviation";//210327
	string  str_20 = "./191120KNNMatlab/MethodSumDeviation";//191206
	string  str_16 = "./191120KNNMatlab/MethodAccuracy";//191204 Accuracy of prune power
	string  str_7 = "./191120KNNMatlab/MethodKNNRunTime";//191123 size is total method number. total run time of every method for matlab barchart
	string  str_44 = "./191120KNNMatlab/MethodKNNRunTimeHasIO";//210606 size is total method number. total run time of every method for matlab barchart
	string  str_12 = "./191120KNNMatlab/MethodApproximationTime";//191123 size is total method number. total run time of every method for matlab barchart
	string  str_50 = "./191120KNNMatlab/MethodBuildTreeTime";//build tree time
	string  str_56 = "./191120KNNMatlab/MethodIngestDataTime";//211211 ingest data time
	string  str_13 = "./191120KNNMatlab/MethodOnlyKNNTime";//191123 size is total method number. total run time of every method for matlab barchart
	string  str_45 = "./191120KNNMatlab/MethodOnlyKNNTimeHasIO";//210606 size is total method number. total run time of every method for matlab barchart
	string  str_53 = "./191120KNNMatlab/MethodIOCost";//211213 ID cost
	string  str_58 = "./191120KNNMatlab/MethodKNNCPUCost";//211213 KNN CPU
	string  str_61 = "./191120KNNMatlab/MethodInternalNodeSize";//211215 Internal node size
	string  str_64 = "./191120KNNMatlab/MethodLeafNodeSize";//211215 Leaf node size
	string  str_67 = "./191120KNNMatlab/MethodTotalNodeSize";//211215 Total node size
	string  str_70 = "./191120KNNMatlab/MethodIndexHeight";//211215 Tree height

	string  str_8 = "./191120KNNMatlab/FileMethodPrunePower";//191125 size is total file number. prune power of every file for matlab barchart
	string  str_36 = "./191120KNNMatlab/FileMethodPrunePowerCombine";//201221  prune power combine. prune power / accuracy
	string  str_37 = "./191120KNNMatlab/FileMethodMaxDeviation";//210327 Max deviaiton 
	string  str_38 = "./191120KNNMatlab/FileMethodMaxWidthDeviation";////210327 Max deviaiton  * width
	string  str_21 = "./191120KNNMatlab/FileMethodSumDeviation";//191206
	string  str_17 = "./191120KNNMatlab/FileMethodAccuracy";//191204 File Accuracy of prune power
	string  str_9 = "./191120KNNMatlab/FileMethodKNNRunTime";//191125 size is total file number. total run time of every file for matlab barchart
	string  str_46 = "./191120KNNMatlab/FileMethodKNNRunTimeHasIO";//210606 size is total file number. total run time of every file for matlab barchart
	string  str_14 = "./191120KNNMatlab/FileMethodApproximationTime";//191125 size is total file number. total run time of every file for matlab barchart
	string  str_51 = "./191120KNNMatlab/FileMethodBuildTreeTime";
	string  str_57 = "./191120KNNMatlab/FileMethodIngestDataTime";// ingest time
	string  str_15 = "./191120KNNMatlab/FileMethodOnlyKNNTime";//191125 size is total file number. total run time of every file for matlab barchart
	string  str_47 = "./191120KNNMatlab/FileMethodOnlyKNNTimeHasIO";//210606 size is total file number. total run time of every file for matlab barchart
	string  str_54 = "./191120KNNMatlab/FileMethodIOCost";//211212 IO
	string  str_59 = "./191120KNNMatlab/FileMethodKNNCPUCost";//211213 KNN CPU
	string  str_62 = "./191120KNNMatlab/FileMethodInternalNodeSize";//211215 Internal node size
	string  str_65 = "./191120KNNMatlab/FileMethodLeafNodeSize";//211215 Leaf node size
	string  str_68 = "./191120KNNMatlab/FileMethodTotalNodeSize";//211215 Total node size
	string  str_71 = "./191120KNNMatlab/FileMethodIndexHeight";//211215 Tree height

	string  str_2 = "./191120KNNMatlab/TotalPrunePower";//191122 prune power
	string  str_34 = "./191120KNNMatlab/TotalPrunePowerCombine";//201221 prune power combine. prune power / accuracy
	string  str_39 = "./191120KNNMatlab/TotalMaxDeviation";//210327 max deviation
	string  str_40 = "./191120KNNMatlab/TotalMaxWidthDeviation";//191204 sum deviation
	string  str_19 = "./191120KNNMatlab/TotalSumDeviation";//191204 sum deviation
	string  str_18 = "./191120KNNMatlab/TotalAccuracy";//191205 Accuracy
	string  str_3 = "./191120KNNMatlab/TotalKNNRunTime";//191122 whole run time
	string  str_48 = "./191120KNNMatlab/TotalKNNRunTimeHasIO";//210606 whole run time has IO
	string  str_22 = "./191120KNNMatlab/TotalApproximationRunTime";//200108 approximation run time
	string  str_49 = "./191120KNNMatlab/TotalBuildTreeTime";
	string  str_55 = "./191120KNNMatlab/TotalIngestDataTime";// 211210 ingest data time
	string  str_23 = "./191120KNNMatlab/TotalIndexRunTime";//200108 index(KNN time) run time
	string  str_43 = "./191120KNNMatlab/TotalIndexRunTimeHasIO";//210606 whole knn run time has IO
	string  str_52 = "./191120KNNMatlab/TotalIOCost";//211212 IO
	string  str_60 = "./191120KNNMatlab/TotalKNNCPUCost";//211213 KNN CPU
	string  str_63 = "./191120KNNMatlab/TotalInternalNodeSize";//211215 Internal node size
	string  str_66 = "./191120KNNMatlab/TotalLeafNodeSize";//211215 Leaf node size
	string  str_69 = "./191120KNNMatlab/TotalTotalNodeSize";//211215 Total node size
	string  str_72 = "./191120KNNMatlab/TotalMethodIndexHeight";//211215 Tree height

   /*----------------------------------------------------------------------------------------------------------------------------------------------*/

	ADDRESS_STRING() {
		/*-------------------------------------------------             String Name           --------------------------------------------------------*/
	// 1K_total_number , 2N_total_number , 3file_total_number, 4method_total_number, 5point number, 6query id, 7max id, 8 single time series length ,9 data type, 10 evaluation type, 11 split methods number, 12 option homogenous data type
		  str_1 = "./191120KNNMatlab/BarChartCoefficient"; //
		  str_25 = "./191120KNNMatlab/NameRepresentationMethod";//200831 the name of representation method,  vector.
		  str_99 = "./191120KNNMatlab/IDQuery";// Query time seirs id vector
		  str_100 = "./191120KNNMatlab/IDRepresentationMethod";//210902 the ID of representation method, number vector.
		  str_101 = "./191120KNNMatlab/OptionBarchartTree";//210902 the option of Tree for representation method.
		  str_26 = "./191120KNNMatlab/IDFileChoosed";//201020 Choose File id: file_id_choose_vector
		  str_4 = "./191120KNNMatlab/KCoefficient";//191122 K records. eg. K = 2, 4 , 8
		  str_11 = "./191120KNNMatlab/NCoefficient";//191128 N records. eg. N = 18, 36, 72, 144, 288
		  str_5 = "./191120KNNMatlab/ReadMe";//191122
		  str_10 = "./191120KNNMatlab/StringSuffix";//191126 prefix of file.

		  str_6 = "./191120KNNMatlab/MethodPrunePower";//191123 size is total method number. prune power of every method for matlab barchart
		  str_35 = "./191120KNNMatlab/MethodPrunePowerCombine";//201221  prune power combine. prune power / accuracy
		  str_41 = "./191120KNNMatlab/MethodMaxDeviation";//210327
		  str_42 = "./191120KNNMatlab/MethodMaxWidthDeviation";//210327
		  str_20 = "./191120KNNMatlab/MethodSumDeviation";//191206
		  str_16 = "./191120KNNMatlab/MethodAccuracy";//191204 Accuracy of prune power
		  str_7 = "./191120KNNMatlab/MethodKNNRunTime";//191123 size is total method number. total run time of every method for matlab barchart
		  str_44 = "./191120KNNMatlab/MethodKNNRunTimeHasIO";//210606 size is total method number. total run time of every method for matlab barchart
		  str_12 = "./191120KNNMatlab/MethodApproximationTime";//191123 size is total method number. total run time of every method for matlab barchart
		  str_50 = "./191120KNNMatlab/MethodBuildTreeTime";
		  str_56 = "./191120KNNMatlab/MethodIngestDataTime";//211211 ingest data time
		  str_13 = "./191120KNNMatlab/MethodOnlyKNNTime";//191123 size is total method number. total run time of every method for matlab barchart
		  str_45 = "./191120KNNMatlab/MethodOnlyKNNTimeHasIO";//210606 size is total method number. total run time of every method for matlab barchart
		  str_53 = "./191120KNNMatlab/MethodIOCost";//211213 ID cost
		  str_58 = "./191120KNNMatlab/MethodKNNCPUCost";//211213 KNN CPU
		  str_61 = "./191120KNNMatlab/MethodInternalNodeSize";//211215 Internal node size
		  str_64 = "./191120KNNMatlab/MethodLeafNodeSize";//211215 Leaf node size
		  str_67 = "./191120KNNMatlab/MethodTotalNodeSize";//211215 Total node size
		  str_70 = "./191120KNNMatlab/MethodIndexHeight";//211215 Tree height

		  str_8 = "./191120KNNMatlab/FileMethodPrunePower";//191125 size is total file number. prune power of every file for matlab barchart
		  str_36 = "./191120KNNMatlab/FileMethodPrunePowerCombine";//201221  prune power combine. prune power / accuracy
		  str_37 = "./191120KNNMatlab/FileMethodMaxDeviation";//210327 Max deviaiton 
		  str_38 = "./191120KNNMatlab/FileMethodMaxWidthDeviation";////210327 Max deviaiton  * width
		  str_21 = "./191120KNNMatlab/FileMethodSumDeviation";//191206
		  str_17 = "./191120KNNMatlab/FileMethodAccuracy";//191204 File Accuracy of prune power
		  str_9 = "./191120KNNMatlab/FileMethodKNNRunTime";//191125 size is total file number. total run time of every file for matlab barchart
		  str_46 = "./191120KNNMatlab/FileMethodKNNRunTimeHasIO";//210606 size is total file number. total run time of every file for matlab barchart
		  str_14 = "./191120KNNMatlab/FileMethodApproximationTime";//191125 size is total file number. total run time of every file for matlab barchart
		  str_51 = "./191120KNNMatlab/FileMethodBuildTreeTime";
		  str_57 = "./191120KNNMatlab/FileMethodIngestDataTime";// ingest time
		  str_15 = "./191120KNNMatlab/FileMethodOnlyKNNTime";//191125 size is total file number. total run time of every file for matlab barchart
		  str_47 = "./191120KNNMatlab/FileMethodOnlyKNNTimeHasIO";//210606 size is total file number. total run time of every file for matlab barchart
		  str_54 = "./191120KNNMatlab/FileMethodIOCost";//211212 IO
		  str_59 = "./191120KNNMatlab/FileMethodKNNCPUCost";//211213 KNN CPU
		  str_62 = "./191120KNNMatlab/FileMethodInternalNodeSize";//211215 Internal node size
		  str_65 = "./191120KNNMatlab/FileMethodLeafNodeSize";//211215 Leaf node size
		  str_68 = "./191120KNNMatlab/FileMethodTotalNodeSize";//211215 Total node size
		  str_71 = "./191120KNNMatlab/FileMethodIndexHeight";//211215 Tree height

		  str_2 = "./191120KNNMatlab/TotalPrunePower";//191122 prune power
		  str_34 = "./191120KNNMatlab/TotalPrunePowerCombine";//201221 prune power combine. prune power / accuracy
		  str_39 = "./191120KNNMatlab/TotalMaxDeviation";//210327 max deviation
		  str_40 = "./191120KNNMatlab/TotalMaxWidthDeviation";//191204 sum deviation
		  str_19 = "./191120KNNMatlab/TotalSumDeviation";//191204 sum deviation
		  str_18 = "./191120KNNMatlab/TotalAccuracy";//191205 Accuracy
		  str_3 = "./191120KNNMatlab/TotalKNNRunTime";//191122 whole run time
		  str_48 = "./191120KNNMatlab/TotalKNNRunTimeHasIO";//210606 whole run time has IO
		  str_22 = "./191120KNNMatlab/TotalApproximationRunTime";//200108 approximation run time
		  str_49 = "./191120KNNMatlab/TotalBuildTreeTime";
		  str_55 = "./191120KNNMatlab/TotalIngestDataTime";
		  str_23 = "./191120KNNMatlab/TotalIndexRunTime";//200108 index(KNN time) run time
		  str_43 = "./191120KNNMatlab/TotalIndexRunTimeHasIO";//210606 whole knn run time has IO
		  str_52 = "./191120KNNMatlab/TotalIOCost";//211212 IO
		  str_60 = "./191120KNNMatlab/TotalKNNCPUCost";//211213 KNN CPU
		  str_63 = "./191120KNNMatlab/TotalInternalNodeSize";//211215 Internal node size
		  str_66 = "./191120KNNMatlab/TotalLeafNodeSize";//211215 Leaf node size
		  str_69 = "./191120KNNMatlab/TotalTotalNodeSize";//211215 Total node size
		  str_72 = "./191120KNNMatlab/TotalMethodIndexHeight";//211215 Tree height
	   /*----------------------------------------------------------------------------------------------------------------------------------------------*/
	}

	template<typename T>
	ADDRESS_STRING(const T& const suffix) {
		const string s = suffix + ".txt";

		/*-------------------------------------------------             String Name           --------------------------------------------------------*/
		str_1 += s; //191122 chart coefficient
		str_25 += s;//200831 the name of representation method, string vector.
		str_99 += s;
		str_100 += s;//210901 210902 the ID of representation method, number vector.
		str_101 += s;
		str_26 += s;//201020 Choose File id: file_id_choose_vector
		str_4 += s;//191122 K records. eg. K = 2, 4 , 8
		str_11 += s;//191128 N records. eg. N = 18, 36, 72, 144, 288
		str_5 += suffix;//191122 Read me
		str_10 += ".txt";//191126 prefix of file.

		str_6 += s;//191123 size is total method number. prune power of every method for matlab barchart
		str_35 += s;//201221  prune power combine. prune power / accuracy
		str_41 += s;//210327
		str_42 += s;//210327
		str_20 += s;//191206
		str_16 += s;//191204 Accuracy of prune power
		str_7 += s;//191123 size is total method number. total run time of every method for matlab barchart
		str_44 += s;//210606 size is total method number. total run time of every method for matlab barchart
		str_12 += s;//191123 size is total method number. total run time of every method for matlab barchart
		str_50 += s;
		str_56 += s;//211211 ingest data time
		str_13 += s;//191123 size is total method number. total run time of every method for matlab barchart
		str_45 += s;//210606 size is total method number. total run time of every method for matlab barchart
		str_53 += s;
		str_58 += s;
		str_61 += s;//211215 Internal node size
		str_64 += s;//211215 Leaf node size
		str_67 += s;//211215 Total node size
		str_70 += s;//211215 Tree height

		str_8 += s;//191125 size is total file number. prune power of every file for matlab barchart
		str_36 += s;//201221  prune power combine. prune power / accuracy
		str_37 += s;//210327 Max deviaiton 
		str_38 += s;////210327 Max deviaiton  * width
		str_21 += s;//191206
		str_17 += s;//191204 File Accuracy of prune power
		str_9 += s;//191125 size is total file number. total run time of every file for matlab barchart
		str_46 += s;//210606 size is total file number. total run time of every file for matlab barchart
		str_14 += s;//191125 size is total file number. total run time of every file for matlab barchart
		str_51 += s;
		str_57 += s;
		str_15 += s;//191125 size is total file number. total run time of every file for matlab barchart
		str_47 += s;//210606 size is total file number. total run time of every file for matlab barchart
		str_54 += s;
		str_59 += s;
		str_62 += s;//211215 Internal node size
		str_65 += s;//211215 Leaf node size
		str_68 += s;//211215 Total node size
		str_71 += s;//211215 Tree height

		str_2 += s;//191122 prune power
		str_34 += s;//201221 prune power combine. prune power / accuracy
		str_39 += s;//210327 max deviation
		str_40 += s;//191204 sum deviation
		str_19 += s;//191204 sum deviation
		str_18 += s;//191205 Accuracy
		str_3 += s;//191122 whole run time
		str_48 += s;//210606 whole run time has IO
		str_22 += s;//200108 approximation run time
		str_49 += s;
		str_55 += s;
		str_23 += s;//200108 index(KNN time) run time
		str_43 += s;//210606 whole knn run time has IO
		str_52 += s;
		str_60 += s;
		str_63 += s;//211215 Internal node size
		str_66 += s;//211215 Leaf node size
		str_69 += s;//211215 Total node size
		str_72 += s;//211215 Tree height

		/*----------------------------------------------------------------------------------------------------------------------------------------------*/
	}
};

//210820 for each method
TEMPLATE
template<typename T>
struct Evaluation::RESULT_COEFFICIENT {
	/*********************************           Experiment Argument            ****************************************/
	string str_suffix;
	// 1 K_total_number , 2 N_total_number , 3 file_total_number, 4 method_total_number, 5 point number, 6 query size, 7 max id, 
	// 8 single time series length ,9 data type, 10 evaluation type, 11 split methods number, 12 option homogenous data type
	vector<T> barchart_coefficient_vector;
	vector<T> file_id_choose_vector;//201020 store specific file id
	vector<T> data_list_heterogeneous;//210603 For mixed dataset id
	vector<T> N_coefficient_vector;//191128  N coefficient for Matlab Barchart
	vector<T> K_coefficient_vector;//191120 K coefficient for Matlab Barchart
	vector<T> representation_option_vector;// id of representation old or new
	vector<string> name_representation_vector_string_vector;
	vector<T> type_index_vector;//211213 barchart color
	vector<T> id_query_series_vector;//211217 id of query
	/*************************************************************************************************************/

	RESULT_COEFFICIENT() {
	}

	~RESULT_COEFFICIENT() {
		str_suffix.clear();
		str_suffix.shrink_to_fit();
		barchart_coefficient_vector.clear();
		barchart_coefficient_vector.shrink_to_fit();
		file_id_choose_vector.clear();//201020 store specific file id
		file_id_choose_vector.shrink_to_fit();//201020 store specific file id
		data_list_heterogeneous.clear();//210603 mixed dataset id
		data_list_heterogeneous.shrink_to_fit();//210603 mixed dataset id
		N_coefficient_vector.clear();//191128  N coefficient for Matlab Barchart
		N_coefficient_vector.shrink_to_fit();//191128  N coefficient for Matlab Barchart
		K_coefficient_vector.clear();//191120 K coefficient for Matlab Barchart
		K_coefficient_vector.shrink_to_fit();//191120 K coefficient for Matlab Barchart
		representation_option_vector.clear();
		representation_option_vector.shrink_to_fit();
		name_representation_vector_string_vector.clear();
		name_representation_vector_string_vector.shrink_to_fit();
		type_index_vector.clear();
		type_index_vector.shrink_to_fit();
		id_query_series_vector.clear();
		id_query_series_vector.shrink_to_fit();
	}
};

//210820 for each method
TEMPLATE
template<typename T>
struct Evaluation::METHOD_RESULT {

	vector<vector<T>*> method_result_pointer_vector;

	/*********************************           For Each Mehtod            ****************************************/
	vector<T> method_sum_prune_power_combine_vector;//201221 old p / accuracy
	vector<T> method_sum_prune_power_vector;//191123 totl prune power of every method for matlab barchart
	vector<T> method_sum_sum_deviation_vector;//191206 sum deviation vector
	vector<T> method_sum_max_deviation_vector;//191206 sum deviation vector
	vector<T> method_sum_max_width_deviation_vector;//191206 sum deviation vector
	vector<T> method_sum_accuracy_vector;//191204 accuracy of prune power
	vector<T> method_run_time_vector;//191123 totl run time of every method for matlab barchart
	vector<T> method_run_time_has_IO_vector;//210606 totl run time of every method for matlab barchart
	vector<T> method_approximation_time_vector;//191204 For approximation time of every method
	vector<T> method_build_tree_time_vector;
	vector<T> method_ingest_data_time_vector;
	vector<T> method_knn_time_vector;//191204 For knn time of every method
	vector<T> method_knn_time_has_IO_vector;//210606 For knn time of every method
	vector<T> method_IO_vector;
	vector<T> method_knn_cpu_time_vector;//58
	vector<T> method_internal_node_size_vector;//61
	vector<T> method_leaf_node_size_vector;//64
	vector<T> method_total_node_size_vector;//67
	vector<T> method_index_height_vector;//70
	/*************************************************************************************************************/

	inline void point_to_vector() {
		method_result_pointer_vector = {
			&method_sum_max_deviation_vector, &method_approximation_time_vector, &method_sum_prune_power_vector, &method_sum_accuracy_vector,
			&method_build_tree_time_vector, &method_knn_time_vector, &method_knn_time_has_IO_vector, &method_run_time_vector,
			&method_run_time_has_IO_vector,&method_ingest_data_time_vector, &method_IO_vector, &method_knn_cpu_time_vector,
			&method_internal_node_size_vector, &method_leaf_node_size_vector, &method_total_node_size_vector, &method_index_height_vector
		};
	}

	template<typename T, typename Y, typename U>
	void assign_pointer_vector(const T& const size_assign, const Y& const value_assign, U& const pointer_vector) {//210827 for pointer vector
		for (int id_pointer = 0; id_pointer < pointer_vector.size(); id_pointer++) {
			(*pointer_vector[id_pointer]).assign(size_assign, value_assign);
		}
	}

	METHOD_RESULT() {
		point_to_vector();
	}
	
	template<typename Y>
	METHOD_RESULT(const Y& const method_total_number) {

		 point_to_vector();
		 assign_pointer_vector(method_total_number, 0, method_result_pointer_vector);

		 method_sum_prune_power_combine_vector.assign(method_total_number, 0);//201221 old p / accuracy
		 method_sum_prune_power_vector.assign(method_total_number, 0);//191123 totl prune power of every method for matlab barchart
		 method_sum_sum_deviation_vector.assign(method_total_number, 0);//191206 sum deviation vector
		 method_sum_max_deviation_vector.assign(method_total_number, 0);//191206 sum deviation vector
		 method_sum_max_width_deviation_vector.assign(method_total_number, 0);//191206 sum deviation vector
		 method_sum_accuracy_vector.assign(method_total_number, 0);//191204 accuracy of prune power
		 method_run_time_vector.assign(method_total_number, 0);//191123 totl run time of every method for matlab barchart
		 method_run_time_has_IO_vector.assign(method_total_number, 0);//210606 totl run time of every method for matlab barchart
		 method_approximation_time_vector.assign(method_total_number, 0);//191204 For approximation time of every method
		 method_build_tree_time_vector.assign(method_total_number, 0);
		 method_ingest_data_time_vector.assign(method_total_number, 0);
		 method_knn_time_vector.assign(method_total_number, 0);//191204 For knn time of every method
		 method_knn_time_has_IO_vector.assign(method_total_number, 0);//210606 For knn time of every method
		 method_IO_vector.assign(method_total_number, 0);
		 method_knn_cpu_time_vector.assign(method_total_number, 0);//58
	}

	~METHOD_RESULT() {
		method_result_pointer_vector.clear();
		method_result_pointer_vector.shrink_to_fit();

		method_sum_prune_power_vector.clear();//191123 totl prune power of every method for matlab barchart
		method_sum_prune_power_vector.shrink_to_fit();
		method_sum_prune_power_combine_vector.clear();//201221
		method_sum_prune_power_combine_vector.shrink_to_fit();
		method_sum_sum_deviation_vector.clear();
		method_sum_sum_deviation_vector.shrink_to_fit();
		method_sum_accuracy_vector.clear();//191204 accuracy of prune power
		method_sum_accuracy_vector.shrink_to_fit();
		method_run_time_vector.clear();//191123 totl run time of every method for matlab barchart
		method_run_time_vector.shrink_to_fit();
		method_run_time_has_IO_vector.clear();//210606 totl run time of every method for matlab barchart
		method_run_time_has_IO_vector.shrink_to_fit();
		method_approximation_time_vector.clear();//191204 For approximation time of every method
		method_approximation_time_vector.shrink_to_fit();
		method_build_tree_time_vector.clear();
		method_build_tree_time_vector.shrink_to_fit();
		method_ingest_data_time_vector.clear();
		method_ingest_data_time_vector.shrink_to_fit();
		method_knn_time_vector.clear();//191204 For knn time of every method
		method_knn_time_vector.shrink_to_fit();
		method_knn_time_has_IO_vector.clear();//210606For knn time of every method
		method_knn_time_has_IO_vector.shrink_to_fit();
		method_IO_vector.clear();
		method_IO_vector.shrink_to_fit();
		method_knn_cpu_time_vector.clear();//58
		method_knn_cpu_time_vector.shrink_to_fit();
	}
};

//210820 Result for each file
TEMPLATE
template<typename T>
struct Evaluation::METHOD_FILE_RESULT {

	vector<vector<T>*> method_file_result_pointer_vector;

	/**********************************      For Each method, Each File      *************************************/
	vector<T> method_file_prune_power_combine_vector;// 201221 old p / accuracy
	vector<T> method_file_prune_power_vector;//191125 5*23 sum prune power for every file [mehtod number, file number]
	vector<T> method_file_sum_deviation_vector;//191206
	vector<T> method_file_max_deviation_vector;//191206
	vector<T> method_file_max_width_deviation_vector;//191206
	vector<T> method_file_accuracy_vector;//191204 accuracy of prune power
	vector<T> method_file_run_time_vector;//191125 5*23 sum run time for every file [mehtod number, file number]
	vector<T> method_file_run_time_has_IO_vector;//210606 5*23 sum run time for every file [mehtod number, file number]
	vector<T> method_file_approximation_time_vector;//191204 5*23 sum approximation time for every file [mehtod number, file number]
	vector<T> method_file_build_tree_time_vector;
	vector<T> method_file_ingest_data_time_vector;
	vector<T> method_file_knn_time_vector;//191204 5*23 sum knn time for every file [mehtod number, file number]
	vector<T> method_file_knn_time_has_IO_vector;//210606 5*23 sum knn time for every file [mehtod number, file number]
	vector<T> method_file_IO_vector;
	vector<T> method_file_knn_cpu_time_vector;//59
	vector<T> method_file_internal_node_size_vector;//62
	vector<T> method_file_leaf_node_size_vector;//65
	vector<T> method_file_total_node_size_vector;//68
	vector<T> method_file_index_height_vector;//71
	/*************************************************************************************************************/

	inline void point_to_vector() {

		method_file_result_pointer_vector = {
			&method_file_max_deviation_vector, &method_file_approximation_time_vector, &method_file_prune_power_vector, &method_file_accuracy_vector,
			&method_file_build_tree_time_vector, &method_file_knn_time_vector, &method_file_knn_time_has_IO_vector, &method_file_run_time_vector,
			&method_file_run_time_has_IO_vector, &method_file_ingest_data_time_vector, &method_file_IO_vector, &method_file_knn_cpu_time_vector,
			&method_file_internal_node_size_vector, &method_file_leaf_node_size_vector, &method_file_total_node_size_vector, &method_file_index_height_vector
		};
	}

	template<typename T, typename Y, typename U>
	void assign_pointer_vector(const T& const size_assign, const Y& const value_assign, U& const pointer_vector) {//210827 for pointer vector
		for (int id_pointer = 0; id_pointer < pointer_vector.size(); id_pointer++) {
			(*pointer_vector[id_pointer]).assign(size_assign, value_assign);
		}
	}

	METHOD_FILE_RESULT() {
		point_to_vector();
		/**********************************      For Each method, Each File      *************************************/
		/*************************************************************************************************************/
	}

	template<typename T, typename Y>
	METHOD_FILE_RESULT(const T& const method_total_number, const Y& const file_total_number) {
		point_to_vector();
		assign_pointer_vector(method_total_number * file_total_number, 0, method_file_result_pointer_vector);

		/**********************************      For Each method, Each File      *************************************/
		method_file_prune_power_combine_vector.assign(method_total_number * file_total_number, 0.0);// 201221 old p / accuracy
		method_file_prune_power_vector.assign(method_total_number * file_total_number, 0.0);//191125 5*23 sum prune power for every file [mehtod number, file number]
		method_file_sum_deviation_vector.assign(method_total_number * file_total_number, 0.0);//191206
		method_file_max_deviation_vector.assign(method_total_number * file_total_number, 0.0);//191206
		method_file_max_width_deviation_vector.assign(method_total_number * file_total_number, 0.0);//191206
		method_file_accuracy_vector.assign(method_total_number * file_total_number, 0.0);//191204 accuracy of prune power
		method_file_run_time_vector.assign(method_total_number * file_total_number, 0.0);//191125 5*23 sum run time for every file [mehtod number, file number]
		method_file_run_time_has_IO_vector.assign(method_total_number * file_total_number, 0.0);//210606 5*23 sum run time for every file [mehtod number, file number]
		method_file_approximation_time_vector.assign(method_total_number * file_total_number, 0.0);//191204 5*23 sum approximation time for every file [mehtod number, file number]
		method_file_build_tree_time_vector.assign(method_total_number * file_total_number, 0.0);
		method_file_ingest_data_time_vector.assign(method_total_number * file_total_number, 0.0);
		method_file_knn_time_vector.assign(method_total_number * file_total_number, 0.0);//191204 5*23 sum knn time for every file [mehtod number, file number]
		method_file_knn_time_has_IO_vector.assign(method_total_number * file_total_number, 0.0);//210606 5*23 sum knn time for every file [mehtod number, file number]
		method_file_IO_vector.assign(method_total_number * file_total_number, 0.0);
		method_file_knn_cpu_time_vector.assign(method_total_number * file_total_number, 0.0);//59
		/*************************************************************************************************************/
	}

	template<typename T>
	METHOD_FILE_RESULT(const T& const method_file_total_number) {
		point_to_vector();
		assign_pointer_vector(method_file_total_number, 0, method_file_result_pointer_vector);

		/**********************************      For Each method, Each File      *************************************/
		method_file_prune_power_combine_vector.assign(method_file_total_number, 0.0);// 201221 old p / accuracy
		method_file_prune_power_vector.assign(method_file_total_number, 0.0);//191125 5*23 sum prune power for every file [mehtod number, file number]
		method_file_sum_deviation_vector.assign(method_file_total_number, 0.0);//191206
		method_file_max_deviation_vector.assign(method_file_total_number, 0.0);//191206
		method_file_max_width_deviation_vector.assign(method_file_total_number, 0.0);//191206
		method_file_accuracy_vector.assign(method_file_total_number, 0.0);//191204 accuracy of prune power
		method_file_run_time_vector.assign(method_file_total_number, 0.0);//191125 5*23 sum run time for every file [mehtod number, file number]
		method_file_run_time_has_IO_vector.assign(method_file_total_number, 0.0);//210606 5*23 sum run time for every file [mehtod number, file number]
		method_file_approximation_time_vector.assign(method_file_total_number, 0.0);//191204 5*23 sum approximation time for every file [mehtod number, file number]
		method_file_build_tree_time_vector.assign(method_file_total_number, 0.0);
		method_file_ingest_data_time_vector.assign(method_file_total_number, 0.0);
		method_file_knn_time_vector.assign(method_file_total_number, 0.0);//191204 5*23 sum knn time for every file [mehtod number, file number]
		method_file_knn_time_has_IO_vector.assign(method_file_total_number, 0.0);//210606 5*23 sum knn time for every file [mehtod number, file number]
		method_file_IO_vector.assign(method_file_total_number, 0.0);
		method_file_knn_cpu_time_vector.assign(method_file_total_number, 0.0);//59
		/*************************************************************************************************************/
	}

	~METHOD_FILE_RESULT() {
		method_file_prune_power_vector.clear();//191125 5*23 sum prune power for every file [mehtod number, file number]
		method_file_prune_power_vector.shrink_to_fit();
		method_file_prune_power_combine_vector.clear();//201221
		method_file_prune_power_combine_vector.shrink_to_fit();
		method_file_sum_deviation_vector.clear();
		method_file_sum_deviation_vector.shrink_to_fit();
		method_file_accuracy_vector.clear();//191204 accuracy of prune power
		method_file_accuracy_vector.shrink_to_fit();
		method_file_run_time_vector.clear();//191125 5*23 sum run time for every file [mehtod number, file number]
		method_file_run_time_vector.shrink_to_fit();
		method_file_run_time_has_IO_vector.clear();//210606 5*23 sum run time for every file [mehtod number, file number]
		method_file_run_time_has_IO_vector.shrink_to_fit();
		method_file_approximation_time_vector.clear();//191204 5*23 sum approximation time for every file [mehtod number, file number]
		method_file_approximation_time_vector.shrink_to_fit();
		method_file_build_tree_time_vector.clear();
		method_file_build_tree_time_vector.shrink_to_fit();
		method_file_ingest_data_time_vector.clear();
		method_file_ingest_data_time_vector.shrink_to_fit();
		method_file_knn_time_vector.clear();//191204 5*23 sum knn time for every file [mehtod number, file number]
		method_file_knn_time_vector.shrink_to_fit();
		method_file_knn_time_has_IO_vector.clear();//210606 5*23 sum knn time for every file [mehtod number, file number]
		method_file_knn_time_has_IO_vector.shrink_to_fit();
		method_file_IO_vector.clear();
		method_file_IO_vector.shrink_to_fit();
		method_file_knn_cpu_time_vector.clear();
		method_file_knn_cpu_time_vector.shrink_to_fit();
	}
};

//210820 for each bar
TEMPLATE
template<typename T>
struct Evaluation::TOTAL_RESULT {

	vector<vector<T>*> total_result_pointer_vector;

	/*******************  Each Detail results   **********************/
	vector<T> total_prune_power_combine_vector;//201221 old p / accuracy
	vector<T> total_prune_power_vector;//191119
	vector<T> total_sum_deviation_vector;//191206
	vector<T> total_max_deviation_vector;//191206
	vector<T> total_max_width_deviation_vector;//191206
	vector<T> total_accuracy_vector;//191204
	vector<T> total_run_time_vector;//191119
	vector<T> total_run_time_has_IO_vector;//210606
	vector<T> total_approximation_time_vector;//200108
	vector<T> total_build_tree_time_vector;
	vector<T> total_ingest_data_time_vector;
	vector<T> total_knn_time_vector;//200108
	vector<T> total_knn_time_has_IO_vector;//210606
	vector<T> total_IO_vector;
	vector<T> total_knn_cpu_time_vector;//60 211213
	vector<T> total_internal_node_size_vector;//63
	vector<T> total_leaf_node_size_vector;//66
	vector<T> total_total_node_size_vector;//69
	vector<T> total_index_height_vector;//72
	/***************************************************************/

	inline void point_to_vector() {
		total_result_pointer_vector = {
			&total_max_deviation_vector, &total_approximation_time_vector, &total_prune_power_vector, &total_accuracy_vector,
			&total_build_tree_time_vector, &total_knn_time_vector, &total_knn_time_has_IO_vector, &total_run_time_vector,
			&total_run_time_has_IO_vector, &total_ingest_data_time_vector, &total_IO_vector, &total_knn_cpu_time_vector,
			&total_internal_node_size_vector, &total_leaf_node_size_vector, &total_total_node_size_vector, &total_index_height_vector
		};
	}

	template<typename T, typename Y, typename U>
	void assign_pointer_vector(const T& const size_assign, const Y& const value_assign, U& const pointer_vector) {//210827 for pointer vector
		for (int id_pointer = 0; id_pointer < pointer_vector.size(); id_pointer++) {
			(*pointer_vector[id_pointer]).assign(size_assign, value_assign);
		}
	}

	TOTAL_RESULT() {
		point_to_vector();
	}

	template<typename T>//210825
	TOTAL_RESULT(const T& const total_number_in_vector) {
		point_to_vector();
		assign_pointer_vector(total_number_in_vector, 0, total_result_pointer_vector);

		total_prune_power_vector.assign(total_number_in_vector, 0.0);//191119
		total_prune_power_combine_vector.assign(total_number_in_vector, 0.0);//201221
		total_sum_deviation_vector.assign(total_number_in_vector, 0.0);//191206
		total_max_deviation_vector.assign(total_number_in_vector, 0.0);//191206
		total_max_width_deviation_vector.assign(total_number_in_vector, 0.0);//191206
		total_accuracy_vector.assign(total_number_in_vector, 0.0);//191204
		total_run_time_vector.assign(total_number_in_vector, 0.0);//191119
		total_run_time_has_IO_vector.assign(total_number_in_vector, 0.0);//210606
		total_approximation_time_vector.assign(total_number_in_vector, 0.0);
		total_build_tree_time_vector.assign(total_number_in_vector, 0.0);
		total_ingest_data_time_vector.assign(total_number_in_vector, 0.0);
		total_knn_time_vector.assign(total_number_in_vector, 0.0);
		total_knn_time_has_IO_vector.assign(total_number_in_vector, 0.0);//210606
		total_IO_vector.assign(total_number_in_vector, 0.0);//211213
		total_knn_cpu_time_vector.assign(total_number_in_vector, 0.0);//60 211213
	}

	~TOTAL_RESULT() {
		total_prune_power_vector.clear();//191119
		total_prune_power_vector.shrink_to_fit();
		total_prune_power_combine_vector.clear();//201221
		total_prune_power_combine_vector.shrink_to_fit();
		total_sum_deviation_vector.clear();//191206
		total_sum_deviation_vector.shrink_to_fit();
		total_max_deviation_vector.clear();
		total_max_deviation_vector.shrink_to_fit();
		total_max_width_deviation_vector.clear();
		total_max_width_deviation_vector.shrink_to_fit();
		total_accuracy_vector.clear();//191204
		total_accuracy_vector.shrink_to_fit();
		total_run_time_vector.clear();//191119
		total_run_time_vector.shrink_to_fit();
		total_run_time_has_IO_vector.clear();//210606
		total_run_time_has_IO_vector.shrink_to_fit();
		total_approximation_time_vector.clear();
		total_approximation_time_vector.shrink_to_fit();
		total_build_tree_time_vector.clear();
		total_build_tree_time_vector.shrink_to_fit();
		total_ingest_data_time_vector.clear();
		total_ingest_data_time_vector.shrink_to_fit();
		total_knn_time_vector.clear();
		total_knn_time_vector.shrink_to_fit();
		total_knn_time_has_IO_vector.clear();//210606
		total_knn_time_has_IO_vector.shrink_to_fit();
		total_IO_vector.clear();
		total_IO_vector.shrink_to_fit();
		total_knn_cpu_time_vector.clear();//60 211213
		total_knn_cpu_time_vector.shrink_to_fit();//60 211213
	}
};

TEMPLATE
template<typename T, typename Y>
inline void Evaluation::initial_address(const T& const suffix, Y& const Address_String) {//210823 Add suffix and "txt?

	const string s = suffix + ".txt";

	/*-------------------------------------------------             String Name           --------------------------------------------------------*/
	Address_String.str_1 += s; //191122 chart coefficient
	Address_String.str_25  += s;//200831 the name of representation method, string vector.
	Address_String.str_99 += s;
	Address_String.str_100 += s;
	Address_String.str_101 += s;
	Address_String.str_26  += s;//201020 Choose File id: file_id_choose_vector
	Address_String.str_4  += s;//191122 K records. eg. K = 2, 4 , 8
	Address_String.str_11  += s;//191128 N records. eg. N = 18, 36, 72, 144, 288
	Address_String.str_5 += suffix;//191122
	Address_String.str_10  += ".txt";//191126 prefix of file.

	Address_String.str_6  += s;//191123 size is total method number. prune power of every method for matlab barchart
	Address_String.str_35 += s;//201221  prune power combine. prune power / accuracy
	Address_String.str_41 += s;//210327
	Address_String.str_42 += s;//210327
	Address_String.str_20 += s;//191206
	Address_String.str_16 += s;//191204 Accuracy of prune power
	Address_String.str_7  += s;//191123 size is total method number. total run time of every method for matlab barchart
	Address_String.str_44 += s;//210606 size is total method number. total run time of every method for matlab barchart
	Address_String.str_12 += s;//191123 size is total method number. total run time of every method for matlab barchart
	Address_String.str_50 += s;
	Address_String.str_56 += s;//211211 ingest data time
	Address_String.str_13 += s;//191123 size is total method number. total run time of every method for matlab barchart
	Address_String.str_45 += s;//210606 size is total method number. total run time of every method for matlab barchart
	Address_String.str_53 += s;//211213 ID cost
	Address_String.str_58 += s;//211213 KNN CPU
	Address_String.str_61 += s;//211215 Internal node size
	Address_String.str_64 += s;//211215 Leaf node size
	Address_String.str_67 += s;//211215 Total node size
	Address_String.str_70 += s;//211215 Tree height

	Address_String.str_8   +=s;//191125 size is total file number. prune power of every file for matlab barchart
	Address_String.str_36   +=s;//201221  prune power combine. prune power / accuracy
	Address_String.str_37  +=s;//210327 Max deviaiton 
	Address_String.str_38   +=s;////210327 Max deviaiton  * width
	Address_String.str_21  +=s;//191206
	Address_String.str_17  +=s;//191204 File Accuracy of prune power
	Address_String.str_9  +=s;//191125 size is total file number. total run time of every file for matlab barchart
	Address_String.str_46  +=s;//210606 size is total file number. total run time of every file for matlab barchart
	Address_String.str_14  +=s;//191125 size is total file number. total run time of every file for matlab barchart
	Address_String.str_51 += s;
	Address_String.str_15  +=s;//191125 size is total file number. total run time of every file for matlab barchart
	Address_String.str_47  +=s;//210606 size is total file number. total run time of every file for matlab barchart
	Address_String.str_57 += s;
	Address_String.str_54 += s;
	Address_String.str_59 += s;
	Address_String.str_62 += s;//211215 Internal node size
	Address_String.str_65 += s;//211215 Leaf node size
	Address_String.str_68 += s;//211215 Total node size
	Address_String.str_71 += s;//211215 Tree height


	Address_String.str_2  +=s;//191122 prune power
	Address_String.str_34  +=s;//201221 prune power combine. prune power / accuracy
	Address_String.str_39  +=s;//210327 max deviation
	Address_String.str_40  +=s;//191204 sum deviation
	Address_String.str_19  +=s;//191204 sum deviation
	Address_String.str_18  +=s;//191205 Accuracy
	Address_String.str_3  +=s;//191122 whole run time
	Address_String.str_48  +=s;//210606 whole run time has IO
	Address_String.str_22  +=s;//200108 approximation run time
	Address_String.str_49 += s;// index build
	Address_String.str_55 += s;// data ingest
	Address_String.str_23  +=s;//200108 index(KNN time) run time
	Address_String.str_43  +=s;//210606 whole knn run time has IO
	Address_String.str_52 += s;
	Address_String.str_60 += s;
	Address_String.str_63 += s;//211215 Internal node size
	Address_String.str_66 += s;//211215 Leaf node size
	Address_String.str_69 += s;//211215 Total node size
	Address_String.str_72 += s;//211215 Tree height

	/*----------------------------------------------------------------------------------------------------------------------------------------------*/
}

TEMPLATE
template<typename T, typename Y, typename U>
void Evaluation::assign_pointer_vector(const T& const size_assign, const Y& const value_assign, U& const pointer_vector) {//210827 for pointer vector
	for (int id_pointer = 0; id_pointer < pointer_vector.size(); id_pointer++) {
		(*pointer_vector[id_pointer]).assign(size_assign, value_assign);
	}
}

//211219 for clearing pointer vector
TEMPLATE
template<typename T>
void Evaluation::clear_pointer_vector(T& const pointer_vector) {
	for (int id_pointer = 0; id_pointer < pointer_vector.size(); id_pointer++) {
		(*pointer_vector[id_pointer]).clear();
	}
}

TEMPLATE
template<typename T, typename Y>
inline void Evaluation::initial_address_no_txt(const T& const suffix, Y& const Address_String) {//210823 for each bar
	const string s = suffix;

	/*-------------------------------------------------             String Name           --------------------------------------------------------*/
	Address_String.str_1 += s; //191122 chart coefficient
	Address_String.str_25 += s;//200831 the name of representation method, string vector.
	Address_String.str_99 += s;
	Address_String.str_100 += s;
	Address_String.str_101 += s;
	Address_String.str_26 += s;//201020 Choose File id: file_id_choose_vector
	Address_String.str_4 += s;//191122 K records. eg. K = 2, 4 , 8
	Address_String.str_11 += s;//191128 N records. eg. N = 18, 36, 72, 144, 288
	Address_String.str_5 += s;//191122
	Address_String.str_10;//191126 prefix of file.

	Address_String.str_6 += s;//191123 size is total method number. prune power of every method for matlab barchart
	Address_String.str_35 += s;//201221  prune power combine. prune power / accuracy
	Address_String.str_41 += s;//210327
	Address_String.str_42 += s;//210327
	Address_String.str_20 += s;//191206
	Address_String.str_16 += s;//191204 Accuracy of prune power
	Address_String.str_7 += s;//191123 size is total method number. total run time of every method for matlab barchart
	Address_String.str_44 += s;//210606 size is total method number. total run time of every method for matlab barchart
	Address_String.str_12 += s;//191123 size is total method number. total run time of every method for matlab barchart
	Address_String.str_50 += s;
	Address_String.str_56 += s;//211211 ingest data time
	Address_String.str_13 += s;//191123 size is total method number. total run time of every method for matlab barchart
	Address_String.str_45 += s;//210606 size is total method number. total run time of every method for matlab barchart
	Address_String.str_53 += s;//211213 ID cost
	Address_String.str_58 += s;//211213 KNN CPU
	Address_String.str_61 += s;//211215 Internal node size
	Address_String.str_64 += s;//211215 Leaf node size
	Address_String.str_67 += s;//211215 Total node size
	Address_String.str_70 += s;//211215 Tree height


	Address_String.str_8 += s;//191125 size is total file number. prune power of every file for matlab barchart
	Address_String.str_36 += s;//201221  prune power combine. prune power / accuracy
	Address_String.str_37 += s;//210327 Max deviaiton 
	Address_String.str_38 += s;////210327 Max deviaiton  * width
	Address_String.str_21 += s;//191206
	Address_String.str_17 += s;//191204 File Accuracy of prune power
	Address_String.str_9 += s;//191125 size is total file number. total run time of every file for matlab barchart
	Address_String.str_46 += s;//210606 size is total file number. total run time of every file for matlab barchart
	Address_String.str_14 += s;//191125 size is total file number. total run time of every file for matlab barchart
	Address_String.str_51 += s;
	Address_String.str_15 += s;//191125 size is total file number. total run time of every file for matlab barchart
	Address_String.str_47 += s;//210606 size is total file number. total run time of every file for matlab barchart
	Address_String.str_57 += s;
	Address_String.str_54 += s;
	Address_String.str_59 += s;
	Address_String.str_62 += s;//211215 Internal node size
	Address_String.str_65 += s;//211215 Leaf node size
	Address_String.str_68 += s;//211215 Total node size
	Address_String.str_71 += s;//211215 Tree height

	Address_String.str_2 += s;//191122 prune power
	Address_String.str_34 += s;//201221 prune power combine. prune power / accuracy
	Address_String.str_39 += s;//210327 max deviation
	Address_String.str_40 += s;//191204 sum deviation
	Address_String.str_19 += s;//191204 sum deviation
	Address_String.str_18 += s;//191205 Accuracy
	Address_String.str_3 += s;//191122 whole run time
	Address_String.str_48 += s;//210606 whole run time has IO
	Address_String.str_22 += s;//200108 approximation run time
	Address_String.str_49 += s;
	Address_String.str_55 += s;
	Address_String.str_23 += s;//200108 index(KNN time) run time
	Address_String.str_43 += s;//210606 whole knn run time has IO
	Address_String.str_52 += s;
	Address_String.str_60 += s;
	Address_String.str_63 += s;//211215 Internal node size
	Address_String.str_66 += s;//211215 Leaf node size
	Address_String.str_69 += s;//211215 Total node size
	Address_String.str_72 += s;//211215 Tree height

	/*----------------------------------------------------------------------------------------------------------------------------------------------*/

}


//***************************************************************
// Method:evaluate_ICDE07
// Qualifier: test ICDE07_APLA approximation algorithm
// Input:
// Output: 
// date:191106
// author:
//***************************************************************
TEMPLATE
void Evaluation::evaluate_ICDE07() {
	vector<DataType> original_time_series_vector = { 1,1,1,1,1,1,1,1,10 };

	APLA_ICDE icde07(int(original_time_series_vector.size()), 3);
	icde07.getAPLA_ICDE07(icde07.input_argument, original_time_series_vector);
}


//***************************************************************
// Method:evaluate_approximation
// Qualifier: Evaluate Approxiamtion & Index algorithm
// Input:
// Output: sum deviaiton, max deviation, time , I/O cost and so on.
// date:191106
// author:
//***************************************************************
TEMPLATE
void Evaluation::evaluate_approximation() {
	string file_address = "DataAddress181218.txt";
	string read_file_name = "";
	TOOL ttool;
	APCA_QUAL paa5;
	APCA_QUAL apca5;
	typename APCA_QUAL::APCA paa_array;
	APCA_KNN_QUAL knn;
	APCA_KNN_QUAL knn_PAA;//190619

	int n5 = 512;
	int N5 = INF;
	double total_same_deviation_id_count = 0;
	double total_diff_deviation_id_count = 0;
	//int* ID = new int[20];

	int N_total_number = 4;
	int file_total_number = 23;
	int method_total_number = 5;//APLA PLA APCA PAA ICDE07
	int total_size = N_total_number * file_total_number * method_total_number;

	double deviation_sum;
	double deviation_max;
	vector<double> total_sum_deviation;//190514
	total_sum_deviation.resize(total_size);
	vector<double> total_max_deviation;//190514
	total_max_deviation.resize(total_size);
	vector<double> total_run_time;//190514
	total_run_time.resize(total_size);
	vector<double> total_sum_density;//190619
	total_sum_density.resize(total_size);
	vector<double> total_sum_area;//190619
	total_sum_area.resize(total_size);

	int begin_id = 0;
	/*---------------------190924 Accuracy Difference & Time Difference-----------------------------------*/
	double single_file_accuracy_diffrerence_PLA = 0;
	double single_file_accuracy_diffrerence_APCA = 0;
	double single_file_accuracy_diffrerence_PAA = 0;
	double single_file_accuracy_diffrerence_ICDE07 = 0;
	vector<double> array_PLA_accuracy_difference;
	array_PLA_accuracy_difference.resize(file_total_number, 0);
	vector<double> array_APCA_accuracy_difference;
	array_APCA_accuracy_difference.resize(file_total_number, 0);
	vector<double> array_PAA_accuracy_difference;
	array_PAA_accuracy_difference.resize(file_total_number, 0);
	vector<double> array_ICDE07_accuracy_difference;
	array_ICDE07_accuracy_difference.resize(file_total_number, 0);

	double single_file_time_diffrerence_PLA = 0;
	double single_file_time_diffrerence_APCA = 0;
	double single_file_time_diffrerence_PAA = 0;
	double single_file_time_diffrerence_ICDE07 = 0;
	vector<double> array_PLA_time_difference;
	array_PLA_time_difference.resize(file_total_number, 0);
	vector<double> array_APCA_time_difference;
	array_APCA_time_difference.resize(file_total_number, 0);
	vector<double> array_PAA_time_difference;
	array_PAA_time_difference.resize(file_total_number, 0);
	vector<double> array_ICDE07_time_difference;
	array_ICDE07_time_difference.resize(file_total_number, 0);

	vector<double> array_difference;//190925  Sum Accuracy & Time Difference for PLA - APLA,APCA - APLA,PAA - APLA
	array_difference.resize((method_total_number - 1) * 2, 0);
	/*----------------------------------------------------------------------------------------------------*/
	//DataType test_array[] = { 1,4,5,0,3,6,1,4,1,3,0,4 };
	//DataType test_array[] = {6,7,7,9,5,4,2,1,0,3,6,5};
	//DataType test_array[] = { 5,6,5,4,1,0,3,0,2,2,9,0};
	//DataType* original_time_series = test_array;

	ttool.clearTXTFile("AllSumDeviation");
	ttool.clearTXTFile("./200706AllAPLAEvaluation/AllAPLARightEndpoint");
	ttool.clearTXTFile("AllAPLAInitialSplitID");//190725
	ttool.clearTXTFile("./200706AllAPLAEvaluation/AllAPLAInitialSegmentRightEndpoint");//191004
	ttool.clearTXTFile("./200706AllAPLAEvaluation/APLAReconstructSeries");
	ttool.clearTXTFile("AllAPLAIncreamentArea");
	ttool.clearTXTFile("AllAPLASegmentDensity");
	ttool.clearTXTFile("AllAPLASegmentDeviation");
	ttool.clearTXTFile("./200706AllAPLAEvaluation/AllAPLASumDeviation");
	ttool.clearTXTFile("./200706AllAPLAEvaluation/AllAPLAInitialEndpoint");
	ttool.clearTXTFile("./200706AllAPLAEvaluation/TotalAccuracyDifference");//190924
	ttool.clearTXTFile("./200706AllAPLAEvaluation/TotalTimeDifference");//190925

	//int file_id = 1;
	for (int file_id = 0; file_id < file_total_number; file_id++) {//file 0-24 file_total_number = 24
		single_file_accuracy_diffrerence_PLA = 0;
		single_file_accuracy_diffrerence_APCA = 0;
		single_file_accuracy_diffrerence_PAA = 0;
		single_file_accuracy_diffrerence_ICDE07 = 0;
		single_file_time_diffrerence_PLA = 0;
		single_file_time_diffrerence_APCA = 0;
		single_file_time_diffrerence_PAA = 0;
		single_file_time_diffrerence_ICDE07 = 0;

		ttool.getStringByID(file_address, file_id, read_file_name);
		DataType* original_time_series = new DataType[n5];
		ttool.getFileStreamByID(read_file_name, n5, 10, original_time_series);
		ttool.normalizeStandard(n5, original_time_series);
		cout << "Original time series: " << endl;
		for (int i = 0; i < n5; i++) {
			//cout <<i<<" "<< original_time_series[i] << ",";
			cout << original_time_series[i] << ",";
		}
		cout << endl;

		int initial_N = 12;//6
		/*=========================================191028 Check if Y-Projection Method======================================================================*/
		//map<double, int> whole_difference_map;
		//bool is_y_projection = false;
		typename TOOL::Y_PROJECTION_ARGUMENT y_projection_argument(initial_N / 3);
		DoublyLinkedList<APLA::AREA_COEFFICIENT> all_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
		DoublyLinkedList<APLA::AREA_COEFFICIENT> cluster_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
		/*==================================================================================================================================================*/
		for (N5 = initial_N; N5 <= 96; N5 *= 2) {//N  5-80; 6 96 - 192;384 2
			int APLA_N = N5 / 3;
			int PLA_N = N5 / 2;
			int APCA_N = N5 / 2;
			APLA apla(n5, APLA_N), apla0(n5, APLA_N);
			PLA_QUAL pla(n5, PLA_N);
			APLA_ICDE icde07(n5, APLA_N);//191016 ICDE 07

			ttool.writeSingleResult("NormalizeTImeSeries181218", original_time_series, n5);

			cout << "====================================================================================================" << endl;
			cout << "File id: " << file_id + 1 << " PAA N: " << N5 << " APLA N: " << apla.input_argument.point_dimension << " PLA N: " << PLA_N << " APCA N: " << APCA_N << endl;

			//ttool.getRandomPoint(original_time_series, 50,51);

			//apla.getAPLADifference(apla.input_argument, original_time_series, ID);

			/*===================================================APLA===============================================*/
			//apla.getAPLAArea4LogLoop(apla.input_argument, original_time_series);
			//apla.getAPLAArea5LowerBound(apla.input_argument, original_time_series);
			//apla.getAPLAArea5LowerBoundSpeed(apla.input_argument, original_time_series);
			//apla.getAPLAArea3(apla.input_argument, original_time_series);//baseline function 190514
			//apla.getAPLAArea5LowerBoundSpeedNoErase(apla.input_argument, original_time_series);//190417
			//apla.getAPLAAreaPLAImprove(apla.input_argument, original_time_series);//190605
			//apla.getAPLAByMinMax(apla.input_argument, original_time_series);//190617 Vector
			apla.getAPLAByMinMaxLinkedList(apla.input_argument, original_time_series, y_projection_argument, all_linked_list, cluster_linked_list, apla.output_argument);//190911 Linked List
			assert(apla.input_argument.point_dimension == N5 / 3);//190611 /3
			//apla.getAPLAArea5LowerBoundSpeedNoEraseArray(apla.input_argument, original_time_series);//190418 failed
			//apla.getAPLADeviation(apla.input_argument, original_time_series);
			/*==========================================================================================================*/

			/*=======================================================PLA=================================================*/
			//PLA ReconstructAPLAMinId181221
			//pla.getReconstructionErrorPLA(pla.input_argument, original_time_series, deviation_sum, deviation_max);
			//191104 Use linked list to instead vector
			pla.getReconstructionErrorPLAVector(pla.input_argument, original_time_series, deviation_sum, deviation_max);//190610
			/*===========================================================================================================*/

			/*===============================================APCA========================================*/
			ttool.recordStartTime(ttool.time_record[3]);
			//apca5.initialAPCA(apca, APCA_N);
			//apca5.getAPCAPoint(original_time_series, n5, APCA_N, apca);

			vector<APCA_QUAL::APCA> apca_argument_vector;//190501
			DoublyLinkedList<APCA_QUAL::APCA> apca_argument = DoublyLinkedList<APCA_QUAL::APCA>();//191104 use linked list to instead vector
			//apca5.getAPCAPointVector(original_time_series, n5, APCA_N, apca_argument);//190501
			apca5.getAPCAPointLinkedList(original_time_series, n5, APCA_N, apca_argument);//191104 use linked list to instead vector
			//apca5.getAPCAPointVectorPreMemory(original_time_series, n5, APCA_N, apca_argument);//190611
			apca5.output_argument.run_time = ttool.recordFinishTime(ttool.time_record[3]);
			/*sume deviation*/
			//knn.distanceAE(original_time_series, n5, apca, deviation_sum, deviation_max);
			for (int segment_id = 0; segment_id < apca_argument.size(); segment_id++) {
				apca_argument_vector.emplace_back(apca_argument[segment_id]);
			}
			knn.distanceAEVector(original_time_series, n5, apca_argument_vector, deviation_sum, deviation_max);//190501

			apca5.output_argument.sum_deviation = deviation_sum;
			apca5.output_argument.max_deviation = deviation_max;
			//paa5.deleteAPCA(apca_argument);
			apca_argument.clear();
			//apca_argument.shrink_to_fit();
			/*==========================================================================================*/

			/*=====================================          PAA          ============================================*/
			ttool.recordStartTime(ttool.time_record[4]);
			//paa5.initialAPCA(paa_array, N5);
			//paa5.divideRemainderPAA(original_time_series, paa_array,n5, N5);
			//paa5.getAPCAPoint(original_time_series, n5, N5, paa_array);
			vector<APCA_QUAL::APCA> paa_argument_vector;//190501
			DoublyLinkedList<APCA_QUAL::APCA> paa_argument = DoublyLinkedList<APCA_QUAL::APCA>();//191104 use linked list to instead vector
			//paa5.getPAAPointVector(original_time_series, n5, N5, paa_argument);//190501
			paa5.getPAAPointLinkedList(original_time_series, n5, N5, paa_argument);//191104 use linked list to instead vector
			//paa5.getPAAPointVectorPreMemory(original_time_series, n5, N5, paa_argument);//190611
			paa5.output_argument.run_time = ttool.recordFinishTime(ttool.time_record[4]);
			//knn_PAA.distanceAE(original_time_series, n5, paa_array, deviation_sum, deviation_max);
			for (int segment_id = 0; segment_id < paa_argument.size(); segment_id++) {
				paa_argument_vector.emplace_back(paa_argument[segment_id]);
			}
			knn_PAA.distanceAEVector(original_time_series, n5, paa_argument_vector, deviation_sum, deviation_max);//190501
			paa5.output_argument.sum_deviation = deviation_sum;
			paa5.output_argument.max_deviation = deviation_max;
			//cout <<"PAA deviation: "<< deviation_sum << " " << deviation_max << endl;
			//paa5.deleteAPCA(paa_argument);
			//paa_argument.clear();
			//paa_argument.shrink_to_fit();
			/*========================================================================================================*/

			/*============================================        APLA ICDE07        ============================================*/
			vector<DataType> original_time_series_vector;
			original_time_series_vector.resize(icde07.input_argument.time_series_length);
			std::copy_n(original_time_series, icde07.input_argument.time_series_length, original_time_series_vector.begin());
			//ttool.recordStartTime(ttool.time_record[5]);
			icde07.getAPLA_ICDE07(icde07.input_argument, original_time_series_vector);
			//icde07.output_argument.sum_deviation = apla.output_argument.sum_deviation;
			//icde07.output_argument.run_time = ttool.recordFinishTime(ttool.time_record[4]);
			/*=====================================================================================================================*/
			//ttool.printArray(apca.r, N5);

			/*cout << "normalized time series : ";
			ttool.printArray(original_time_series, n5);*/

			/*===============================Print deviation & time==================================================================================*/
			apla.output_argument.sum_density = double(n5) / apla.output_argument.sum_area;
			pla.output_argument.sum_density = double(n5) / pla.output_argument.sum_area;
			knn.output_argument.sum_density = double(n5) / knn.output_argument.sum_area;
			knn_PAA.output_argument.sum_density = double(n5) / knn_PAA.output_argument.sum_area;
			cout << "-----------------------------------------------------------------------------------" << endl;
			cout << "APLA sum deviation: " << apla.output_argument.sum_deviation << " max diviaiton: " << apla.output_argument.max_deviation << " sum density: " << apla.output_argument.sum_density << " sum_area: " << apla.output_argument.sum_area << endl;
			cout << "PLA sum deviation: " << pla.output_argument.sum_deviation << " max diviaiton: " << pla.output_argument.max_deviation << " sum density: " << pla.output_argument.sum_density << " sum_area: " << pla.output_argument.sum_area << endl;
			cout << "APCA sum deviation: " << apca5.output_argument.sum_deviation << " max diviaiton: " << apca5.output_argument.max_deviation << " sum density: " << knn.output_argument.sum_density << " sum_area: " << knn.output_argument.sum_area << endl;
			//cout << "APLA0 deviation: " << apla0.output_argument.sum_deviation << " " << apla0.output_argument.max_deviation << endl;
			cout << "PAA sum deviation: " << paa5.output_argument.sum_deviation << " max diviaiton: " << paa5.output_argument.max_deviation << " sum density: " << knn_PAA.output_argument.sum_density << " sum_area: " << knn_PAA.output_argument.sum_area << endl;
			cout << "ICDE07 sum deviation: " << icde07.output_argument.sum_deviation << endl;

			cout << "APLA time: " << apla.output_argument.run_time << endl;
			cout << "PLA time: " << pla.output_argument.run_time << endl;
			cout << "APCA time: " << apca5.output_argument.run_time << endl;
			//cout << "APLA0 time: " << apla0.output_argument.run_time << endl;
			cout << "PAA time: " << paa5.output_argument.run_time << endl;
			cout << "ICDE07 time: " << icde07.output_argument.run_time << endl;
			cout << "APLA deviation proportion: " << apla.output_argument.same_deviation_id_count << " " << apla.output_argument.diff_deviation_id_count << endl;
			total_same_deviation_id_count += apla.output_argument.same_deviation_id_count;
			total_diff_deviation_id_count += apla.output_argument.diff_deviation_id_count;
			cout << "......................................................................................" << endl;
			/*========================================================================================================================================*/

			/*===============================     Write deviation & time   ========================================*/
			//vector<long double> run_time = { paa5.output_argument.run_time, pla.output_argument.run_time, apla.output_argument.run_time, apca5.output_argument.run_time};
			vector<long double> run_time = { apla.output_argument.run_time, pla.output_argument.run_time, apca5.output_argument.run_time, paa5.output_argument.run_time, icde07.output_argument.run_time };
			vector<double> array_sum_deviation = { apla.output_argument.sum_deviation, pla.output_argument.sum_deviation, apca5.output_argument.sum_deviation, paa5.output_argument.sum_deviation, icde07.output_argument.sum_deviation };
			vector<double> array_max_deviation = { apla.output_argument.max_deviation, pla.output_argument.max_deviation, apca5.output_argument.max_deviation, paa5.output_argument.max_deviation };
			vector<double> array_sum_density = { apla.output_argument.sum_density, pla.output_argument.sum_density, knn.output_argument.sum_density, knn_PAA.output_argument.sum_density };//190619
			vector<double> array_sum_area = { apla.output_argument.sum_area, pla.output_argument.sum_area, knn.output_argument.sum_area, knn_PAA.output_argument.sum_area };//190619

			ttool.writeResultNoCover("AllSumDeviation", array_sum_deviation);//190625
			ttool.writeResultNoCover("./200706AllAPLAEvaluation/AllAPLASumDeviation", array_sum_deviation);//190701

			time_t now_time = time(0);// system time
			char dt[26];
			ctime_s(dt, sizeof dt, &now_time);// from now to string
			ofstream outfile("Result_compare.txt", ios::app);
			assert(outfile.is_open());
			outfile << endl << dt << "file id: " << file_id << " N: " << N5 << endl;
			//outfile << "    PAA       PLA      APLA     APCA " << endl;
			//std::copy_n(run_time.begin(), run_time.size(), std::ostream_iterator<long double>(outfile, ",   "));
			//outfile << endl;
			outfile << "    APLA      PLA      APCA      PAA     ICDE07" << endl;
			std::copy_n(array_sum_deviation.begin(), array_sum_deviation.size(), std::ostream_iterator<double>(outfile, ",   "));
			outfile << endl;
			std::copy_n(array_max_deviation.begin(), array_max_deviation.size(), std::ostream_iterator<double>(outfile, ",   "));
			outfile << endl;
			std::copy_n(run_time.begin(), run_time.size(), std::ostream_iterator<long double>(outfile, ",   "));
			outfile << endl;
			outfile << endl;
			std::copy_n(array_sum_density.begin(), array_sum_density.size(), std::ostream_iterator<double>(outfile, ",   "));
			outfile << endl;
			outfile << endl;
			std::copy_n(array_sum_area.begin(), array_sum_area.size(), std::ostream_iterator<double>(outfile, ",   "));
			outfile << endl;
			outfile.close();
			/*=======================================================================================================*/

			/*=================================Compute Accuracy&Time Difference: PAA,APCA,PLA - APLA=========================================*/
			//190924
			single_file_accuracy_diffrerence_PLA += array_sum_deviation[1] - array_sum_deviation[0];
			single_file_accuracy_diffrerence_APCA += array_sum_deviation[2] - array_sum_deviation[0];
			single_file_accuracy_diffrerence_PAA += array_sum_deviation[3] - array_sum_deviation[0];
			single_file_accuracy_diffrerence_ICDE07 += array_sum_deviation[4] - array_sum_deviation[0];

			single_file_time_diffrerence_PLA += run_time[1] - run_time[0];
			single_file_time_diffrerence_APCA += run_time[2] - run_time[0];
			single_file_time_diffrerence_PAA += run_time[3] - run_time[0];
			single_file_time_diffrerence_ICDE07 += run_time[4] - run_time[0];
			/*===============================================================================================================================*/

			/*===================================== Result Vector ============================================*/
			std::copy_n(array_sum_deviation.begin(), array_sum_deviation.size(), total_sum_deviation.begin() + begin_id);
			std::copy_n(array_max_deviation.begin(), array_max_deviation.size(), total_max_deviation.begin() + begin_id);
			std::copy_n(run_time.begin(), run_time.size(), total_run_time.begin() + begin_id);
			std::copy_n(array_sum_density.begin(), array_sum_density.size(), total_sum_density.begin() + begin_id);//190619
			std::copy_n(array_sum_area.begin(), array_sum_area.size(), total_sum_area.begin() + begin_id);//190619
			begin_id += method_total_number;
			/*.................................................................................................*/
			//ttool.deleteArray(ID);
		}
		ttool.deleteArray(original_time_series);
		all_linked_list.clear();
		cluster_linked_list.clear();
		/*================================= 190924 Compute Accuracy&Time Difference: PAA,APCA,PLA - APLA=========================================*/
		//190924
		array_PLA_accuracy_difference[file_id] = single_file_accuracy_diffrerence_PLA;
		array_APCA_accuracy_difference[file_id] = single_file_accuracy_diffrerence_APCA;
		array_PAA_accuracy_difference[file_id] = single_file_accuracy_diffrerence_PAA;
		array_ICDE07_accuracy_difference[file_id] = single_file_accuracy_diffrerence_ICDE07;

		array_PLA_time_difference[file_id] = single_file_time_diffrerence_PLA;
		array_APCA_time_difference[file_id] = single_file_time_diffrerence_APCA;
		array_PAA_time_difference[file_id] = single_file_time_diffrerence_PAA;
		array_ICDE07_time_difference[file_id] = single_file_time_diffrerence_ICDE07;

		array_difference[0] += single_file_accuracy_diffrerence_PLA;
		array_difference[1] += single_file_accuracy_diffrerence_APCA;
		array_difference[2] += single_file_accuracy_diffrerence_PAA;
		array_difference[3] += single_file_accuracy_diffrerence_ICDE07;
		array_difference[4] += single_file_time_diffrerence_PLA;
		array_difference[5] += single_file_time_diffrerence_APCA;
		array_difference[6] += single_file_time_diffrerence_PAA;
		array_difference[7] += single_file_time_diffrerence_ICDE07;
		/*=======================================================================================================================================*/
		/*============================================190924 Write Accuracy & Time Difference============================================================*/
		vector<double> single_compare_sum_deviation = { single_file_accuracy_diffrerence_PLA, single_file_accuracy_diffrerence_APCA, single_file_accuracy_diffrerence_PAA , single_file_accuracy_diffrerence_ICDE07 };
		ttool.writeResultNoCover("./200706AllAPLAEvaluation/TotalAccuracyDifference", single_compare_sum_deviation);//190924

		vector<double> single_compare_time = { single_file_time_diffrerence_PLA, single_file_time_diffrerence_APCA, single_file_time_diffrerence_PAA, single_file_time_diffrerence_ICDE07 };
		ttool.writeResultNoCover("./200706AllAPLAEvaluation/TotalTimeDifference", single_compare_time);//190925
		/*===============================================================================================================================================*/

		/*=====================================================Y-Projection Method=======================================================================*/
		all_linked_list.clear();
		cluster_linked_list.clear();
		/*===============================================================================================================================================*/
	}

	//cout <<"Deviation Proportion: " <<total_same_deviation_id_count << " "<< total_diff_deviation_id_count <<endl;
	/*===============================Write N deviation & time & All PLA,APCA,PAA accuracy diff========================================*/
	ttool.writeSingleResult("TotalSumDeviation", total_sum_deviation);
	ttool.writeSingleResult("TotalMaxDeviation", total_max_deviation);
	ttool.writeSingleResult("TotalRunTime", total_run_time);
	ttool.writeSingleResult("TotalSumDensity", total_sum_density);//190619
	ttool.writeSingleResult("TotalSumArea", total_sum_area);//190619

	ttool.writeSingleResult("TotalPLAAccuracyDifference", array_PLA_accuracy_difference);//190924
	ttool.writeSingleResult("TotalAPCAAccuracyDifference", array_APCA_accuracy_difference);//190924
	ttool.writeSingleResult("TotalPAAAccuracyDifference", array_PAA_accuracy_difference);//190924
	ttool.writeSingleResult("TotalICDE07AccuracyDifference", array_ICDE07_accuracy_difference);//191016
	ttool.writeSingleResult("TotalPLATimeDifference", array_PLA_time_difference);//190925
	ttool.writeSingleResult("TotalAPCATimeDifference", array_APCA_time_difference);//190925
	ttool.writeSingleResult("TotalPAATimeDifference", array_PAA_time_difference);//190925
	ttool.writeSingleResult("TotalICDE07TimeDifference", array_ICDE07_time_difference);//191016

	ttool.writeSingleResult("SumAccuracyTimeDifference", array_difference);//190925 Sum Accuracy Time Difference for PLA, APCA, PAA
	for_each_n(array_difference.begin(), array_difference.size(), [&](auto&& au) {au /= double(file_total_number); });
	ttool.writeSingleResult("AverageAccuracyTimeDifference", array_difference);//190925 Average Accuracy Time Difference for PLA, APCA, PAA
	/*.................................................................................................*/
	system("pause");
}

//***************************************************************
// Method:evaluate_multi_KNN
// Qualifier:  evaluate KNN search for multi & single time series of APLA PLA, APCA, PAA & ICDE07_APLA
// Input:
// Output: sum deviaiton, max deviation, time , I/O cost and so on.
// date:191106
// author:
//***************************************************************
TEMPLATE
void Evaluation::evaluate_multi_KNN() {
//	assert(0);
//	//string file_address = "DataAddress181218.txt";
//	string read_file_name = "CinC_ECG_torso_TEST";
//
//	//bool change_file = false;
//	int query_time_series_id = INF;
//	int initial_file_od = INF;//191124
//	int file_id = INF; // the id of file(data set / file)
//	int K = INF;// for KNN, the number of nearest neighbor for query time series 
//	int representation_option = INF;// 1 MSPLA, 2 PLA, 3APCA, 4 PAA, 5 Chebyshev, 6 ICDE07, 7 PLALM
//	int N = INF;//N segment number. degree_m n=m+1. APLA = 4, PAA : 12 , PLA: 6, APCA 6. Cheby 12, ICDE07: 4
//	int initial_N = INF;// for y projection method in APLA, cannot be NULL;
//	int initial_K = INF;//191204 for Yprojction
//	int data_final_dimension = INF;// 1 or 3. Final dimension of data, is single dimension or multi dimension data
//	int arity_d0 = INF;//dimension of time series, most are one dimension
//
//	int n0 = 200;// time series length single is 512, multi dimenson is 300
//
//	/*=========================Data Type========================================*/
//	const size_t& const data_type = 3;//1 single dataset, 2 mixed dataset, 3 multi dimension data set
//	/*==========================================================================*/
//
//	string dataset_attribute[] = { "","Single dataset","Mixed dataset","Multi dimension data set" };
//	/*------------------------------------Mixed DataSet--------------------------------------------------*/
//	string mixed_data_file_name = "./191202DataSet/04mixed_data_set";
//	/*---------------------------------------------------------------------------------------------------*/
//	/*------------------------------------Multi Dimension Dataset----------------------------------------*/
//	string multi_file_name1[] = { "./191202DataSet/MultiDateSet191220/Cricket_X_TEST","./191202DataSet/MultiDateSet191220/Cricket_Y_TEST","./191202DataSet/MultiDateSet191220/Cricket_Z_TEST" };
//	string multi_file_name2[] = { "./191202DataSet/MultiDateSet191220/uWaveGestureLibrary_X_TEST","./191202DataSet/MultiDateSet191220/uWaveGestureLibrary_Y_TEST","./191202DataSet/MultiDateSet191220/uWaveGestureLibrary_Z_TEST" };
//	string multi_dimension_file_name[] = {
//		"Cricket_X_TEST","Cricket_Y_TEST","Cricket_Z_TEST",
//		"uWaveGestureLibrary_X_TEST", "uWaveGestureLibrary_Y_TEST", "uWaveGestureLibrary_Z_TEST",
//	};
//	string* file_address_pointer = nullptr;
//	//string multi_to_single_file_name = "./191202DataSet/MultiToSingle191220/uWaveGestureLibrary";
//	//191220 project multi data set to single data set
//	//TOOL::project_multi_data_to_single_data(multi_file_pointer, 3, 3582, 315, multi_to_single_file_name);
//	//TOOL::project_single_data_to_multi_data_batch(TOOL::single_to_multi_file_address, 13, 3, 40, 600, TOOL::write_single_to_multi_file_address);
//	/*---------------------------------------------------------------------------------------------------*/
//
//	switch (data_type) {
//	case 0: // 0 manipulate data set;    
//		assert(0);
//		break;
//	case 1:// 1 real data set, homogenous(single) data;
//		data_final_dimension = 1;
//		assert(data_final_dimension == 1);
//		file_address_pointer = nullptr;
//		break;
//	case 2:// 2 real data set, heterogeneous(mixed) data.
//		data_final_dimension = 1;
//		assert(data_final_dimension == 1);
//		file_address_pointer = &mixed_data_file_name;
//		break;
//	case 3:// 3 multi dimension homogenous data.
//		//file_address_pointer = multi_dimension_file_name;
//		data_final_dimension = 3;
//		file_address_pointer = nullptr;
//		assert(data_final_dimension == 3);
//		break;
//	case 4:// 4 multi dimension heterogeneous data
//		assert(0);
//		break;
//	case 5://APCA
//		assert(0);
//		break;
//	default:
//		assert(0);
//		break;
//	}
//
//	string chebyshev_write_name[] = { "181206Chebyshev prun power","181206Chebyshev IO cost","181206Chebyshev buld index time","181206Chebyshev navigate time", "181206Chebyshev LB distance","181206Chebyshev euc distane", "181206Chebyshev total KNN time","181206Chebyshev result accuracy" };
//	string* chebyshev_write_pointer = chebyshev_write_name;
//
//	/*----------------------------------------------------------------------------Matlab Bar chart---------------------------------------------------------------------------------*/
//	initial_file_od = 0;
//	initial_N = 24;//12,6
//	initial_K = 2;
//
//	//int point_number0 = 20;// the number of time series. datatype: 1, 2
//	int point_number0 = 40;// the number of time series. datatype: 3
//	int max_node0 = 3;//3 5 for Rtree node, most branches in one node.  datatype: 1, 3
//	//int max_node0 = 5;//for Rtree node, most branches in one node.  datatype: 2
//	int K_total_number = INF;
//	int file_total_number = 39;// single dimension 24, 3dimension data set 13 * 3 = 39
//	int method_total_number = 7;//1 MSPLA, 2 PLA, 3 APCA, 4 PAA, 5 Chebyshev, 6 ICDE07, 7 PLALM
//
//	vector<size_t> N_coefficient_vector;//191128  N coefficient for Matlab Barchart
//	vector<int> K_coefficient_vector;//191120 K coefficient for Matlab Barchart
//	vector<double> total_prune_power_vector;//191119
//	vector<double> total_sum_deviation_vector;//191206
//	vector<double> total_accuracy_vector;//191204
//	vector<double> total_run_time_vector;//191119
//	vector<double> method_sum_prune_power_vector(method_total_number, 0);//191123 totl prune power of every method for matlab barchart
//	vector<double> method_sum_sum_deviation_vector(method_total_number, 0);//191206 sum deviation vector
//	vector<double> method_sum_accuracy_vector(method_total_number, 0);//191204 accuracy of prune power
//	vector<double> method_run_time_vector(method_total_number, 0);//191123 totl run time of every method for matlab barchart
//	vector<double> method_approximation_time_vector(method_total_number, 0);//191204 For approximation time of every method
//	vector<double> method_knn_time_vector(method_total_number, 0);//191204 For knn time of every method
//
//	vector<double> method_file_prune_power_vector(method_total_number * file_total_number, 0.0);//191125 5*23 sum prune power for every file [mehtod number, file number]
//	vector<double> method_file_sum_deviation_vector(method_total_number * file_total_number, 0.0);//191206
//	vector<double> method_file_accuracy_vector(method_total_number * file_total_number, 0.0);//191204 accuracy of prune power
//	vector<double> method_file_run_time_vector(method_total_number * file_total_number, 0.0);//191125 5*23 sum run time for every file [mehtod number, file number]
//	vector<double> method_file_approximation_time_vector(method_total_number * file_total_number, 0.0);//191204 5*23 sum approximation time for every file [mehtod number, file number]
//	vector<double> method_file_knn_time_vector(method_total_number * file_total_number, 0.0);//191204 5*23 sum knn time for every file [mehtod number, file number]
//	/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//
//	/*==========================================read & write mixed data set===============================================*/
//	vector<size_t> data_list;
//	for (int list_id = 0; list_id < file_total_number; list_id++) {
//		data_list.emplace_back(list_id);
//	}
//	typename TOOL::DATA_SOURCE data_source(data_type, data_final_dimension, data_list, file_total_number, point_number0, n0, *file_address_pointer);// 0 manipulate data set; 1 real data set, homogenous(single) data; 2 real data set, heterogeneous(mixed) data. 3 multi single data. 4 multi mixed data
//	TOOL::initial_data_source(data_source);
//	assert(data_source.point_number != INF && data_source.file_operation_number != INF && data_source.single_time_series_length != INF && data_source.multi_single_time_series_length != INF);
//	/*=======================================================================================================================*/
//	/*----------------------      Dimension    ----------------------------------------------*/
//	//for (arity_d0 = 1; arity_d0 <= data_final_dimension; arity_d0 += 2) {//dimension fo time seires. 1:one dimension, 2, two dimension. 3: three dimension. 
//		/*----------------------   File ID      ----------------------------------------------*/
//	for (file_id = initial_file_od; file_id < data_source.file_operation_number; file_id++) {//file 0-24 file_total_number = 23
//		cout << "!!!!!!!!!!!!!==========================================file id: " << file_id + 1 << "==========================================================!!!!!!!!!!!!!!!!" << endl;
//		bool change_file = true;
//		/*=======================================   Query time series  ===================================================================*/
//		query_time_series_id = TOOL::get_random_max(data_source.point_number);
//		/*....................................*/
//		query_time_series_id = data_source.point_number - 5;//19 300 195
//		/*....................................*/
//		vector<DataType> query_time_series_vector(data_source.multi_single_time_series_length, INF);
//		DataType* query_time_series = new DataType[data_source.multi_single_time_series_length];
//		TOOL::get_read_multi_file_address(data_source, file_id);
//#ifdef _DEBUG
//		//vector<DataType> test_query_time_series_vector(data_source.multi_single_time_series_length, INF);
//		//
//		//switch (data_source.data_type) {
//		//case 0: // 0 manipulate data set;    
//		//	assert(0);
//		//	break;
//		//case 1: {// 1 real data set, homogenous(single) data;
//		//	//data_source.read_file_address = data_source.file_address_vector[file_id];
//		//	assert(TOOL::getStringByID(TOOL::file_address, file_id) == data_source.read_file_address && data_source.multi_single_time_series_length == data_source.single_time_series_length);
//		//	TOOL::getFileStreamByID(data_source.read_file_address, data_source.single_time_series_length, query_time_series_id, test_query_time_series_vector);
//		//	break;
//		//}
//		//case 2: {// 2 real data set, heterogeneous(mixed) data.
//		//	assert(!data_source.read_file_address.empty() && data_source.multi_single_time_series_length == data_source.single_time_series_length);
//		//	TOOL::getFileStreamByID(data_source.read_file_address, data_source.single_time_series_length, query_time_series_id, test_query_time_series_vector);
//		//	break;
//		//}
//		//case 3: {// 3 multi dimension homogenous(single) data.
//		//	assert(data_source.read_file_address_vector.size() == data_source.time_series_dimension && data_source.file_operation_number * data_source.time_series_dimension == data_source.file_address_vector.size() && data_source.time_series_dimension == data_final_dimension);
//		//	test_query_time_series_vector.clear();
//		//	test_query_time_series_vector.shrink_to_fit();
//		//	TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, query_time_series_id, test_query_time_series_vector);
//		//	assert(test_query_time_series_vector.size() == data_source.multi_single_time_series_length && data_source.multi_single_time_series_length == data_source.time_series_dimension * data_source.single_time_series_length);
//		//	break;
//		//}
//		//case 4:// 4 multi mixed data
//		//	assert(0);
//		//	break;
//		//case 5://APCA
//		//	assert(0);
//		//	break;
//		//default:
//		//	assert(0);
//		//	break;
//		//}
//		//copy_n(test_query_time_series_vector.begin(), test_query_time_series_vector.size(), query_time_series);
//		//TOOL::normalizeStandard(test_query_time_series_vector);//z-score normalization
//		//TOOL::normalizeStandard(data_source.multi_single_time_series_length, query_time_series);//z-score normalization
//#endif
//
//
//			//already normalized in below function
//		TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, query_time_series_id, query_time_series_vector);
//		copy_n(query_time_series_vector.begin(), query_time_series_vector.size(), query_time_series);
//
//
//#ifdef _DEBUG
//		//assert(test_query_time_series_vector.size() == query_time_series_vector.size());
//		//for (int array_id = 0; array_id < data_source.multi_single_time_series_length; array_id++) {
//		//	assert(float(test_query_time_series_vector[array_id]) == float(query_time_series_vector[array_id]));
//		//	//cout << query_time_series_vector[array_id] << ",";
//		//}
//		//cout << endl;
//#endif
//	/*==================================================================================================================================*/
//		/*================================================  Base KNN  =====================================================================*/
//		std::multiset<pair<double, int>> squential_scan_result_set;
//		TOOL::SimpleBaseKNNSearchMulti(data_source, query_time_series_vector, squential_scan_result_set);
//		assert(!squential_scan_result_set.empty());
//		/*for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//			if (it->original_time_series_id == q_base_queue.top().original_time_series_id) {
//				input_argument.result_accuracy++;
//			}
//			q_base_queue.pop();
//		}
//		double accuracy = input_argument.result_accuracy / double(input_argument.K);*/
//		//cout << "KNN result accuracy: " << accuracy << endl;
//		/*==================================================================================================================================*/
//
//		/*=========================================191028 Initial if Y-Projection Method======================================================================*/
//			//TOOL::Y_PROJECTION_ARGUMENT y_projection_argument(initial_N / 3, change_file);
//			//DoublyLinkedList<APLA::AREA_COEFFICIENT> all_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//			//DoublyLinkedList<APLA::AREA_COEFFICIENT> cluster_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
//			//vector<TOOL::Y_PROJECTION_ARGUMENT> multi_y_projection_argument_vector(data_source.point_number, TOOL::Y_PROJECTION_ARGUMENT(initial_N / 3));
//		vector<TOOL::Y_PROJECTION_ARGUMENT> multi_y_projection_argument_vector(data_source.point_number, TOOL::Y_PROJECTION_ARGUMENT());
//		vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>> multi_all_linked_list(data_source.point_number, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>());//191030
//		vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>> multi_cluster_linked_list(data_source.point_number, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>());//191030
//		/*==================================================================================================================================================*/
//
//		/*----------------------K----------------------------------------------*/
//		for (K = initial_K; K <= 32; K *= 2) {//K 2 4  8 16 32 64 128;;; 16 128
//			if (file_id == initial_file_od) { K_coefficient_vector.emplace_back(K); }
//			/*------------------------------ N ----------------------------------------------*/
//			for (N = initial_N; N <= 192; N *= 2) {//N  5-80; 6 96 - 192;384 2;191119 12 24 48 96 192 258 384;;; 192 96 24, 
//
//				if (file_id == initial_file_od && K == initial_K) { N_coefficient_vector.emplace_back(N); }
//				/*-----------------   representation method   ----------------------------------------------*/
//				for (representation_option = 1; representation_option <= method_total_number; representation_option++) {// 1 MSPLA, 2 PLA, 3APCA, 4 PAA, 5 Chebyshev, 6 ICDE07, 7 PLALM
//					/*------------------------------Initial N£º PAA & Cheby 12, APCA & PLA & PLALM 6, APLA & ICDE07 4-------------------------------------------------*/
//					//representation_option = 5;
//					//change_file = true;
//					MULTI multi(data_source, n0, N, initial_N, file_id, change_file, data_source.time_series_dimension, data_source.point_number, query_time_series_id, max_node0, K, representation_option, file_address_pointer, chebyshev_write_pointer);
//					/*-----------------------------------------------------------------------------------------------------------------------------------------*/
//					//multi.project_multi_to_single(multi_y_projection_argument_vector, multi_all_linked_list, multi_cluster_linked_list);
//					//multi.buid_rtree_knn();
//
//					multi.all_approximation_build_rtree(multi_y_projection_argument_vector, multi_all_linked_list, multi_cluster_linked_list);
//					multi.all_knn(data_source, query_time_series_vector, query_time_series, multi_y_projection_argument_vector, squential_scan_result_set);
//
//					/*-----------------------------------Write Result, Matlab, bar chart----------------------------------------------------------------------*/
//					cout << "---------------------------       Result      ------------------------------\n";
//					cout << "file id : " << file_id + 1 << ", K: " << K << ", N:" << multi.input_argument.point_dimension << ", representation id: " << representation_option << ", n: " << n0 << ", dimension: " << data_source.time_series_dimension << ", point number: " << data_source.point_number << ", query time series id: " << query_time_series_id << ", max sub node: " << max_node0 << endl;
//					cout << "prune power: " << multi.input_argument.pruning_power << endl;
//					cout << "sum deviation: " << multi.output_argument.sum_deviation << endl;
//					cout << "prune accuracy: " << multi.input_argument.result_accuracy << endl;
//					cout << "representation_time: " << multi.input_argument.representation_time << endl;
//					cout << "build_rtree_time: " << multi.input_argument.build_rtree_time << endl;
//					cout << "knn_total_time: " << multi.input_argument.knn_total_time << endl;
//					cout << "accumulation time : " << multi.input_argument.representation_time + multi.input_argument.build_rtree_time + multi.input_argument.knn_total_time << endl;
//					cout << "whole running time: " << multi.input_argument.whole_run_time << endl;
//					cout << "----------------------------------------------------------------------------\n";
//
//					assert(multi.input_argument.pruning_power <= 20 && multi.input_argument.pruning_power != INF && multi.input_argument.whole_run_time != INF);
//					total_prune_power_vector.emplace_back(multi.input_argument.pruning_power);
//					total_sum_deviation_vector.emplace_back(multi.output_argument.sum_deviation);//191206
//					total_accuracy_vector.emplace_back(multi.input_argument.result_accuracy);
//					total_run_time_vector.emplace_back(multi.input_argument.whole_run_time);
//
//					method_sum_prune_power_vector[representation_option - 1] += multi.input_argument.pruning_power;
//					method_sum_sum_deviation_vector[representation_option - 1] += multi.output_argument.sum_deviation;//191206
//					method_sum_accuracy_vector[representation_option - 1] += multi.input_argument.result_accuracy;//191204 accuracy
//					method_run_time_vector[representation_option - 1] += multi.input_argument.whole_run_time;
//					method_approximation_time_vector[representation_option - 1] += multi.input_argument.representation_time;
//					method_knn_time_vector[representation_option - 1] += multi.input_argument.knn_total_time;
//
//					method_file_prune_power_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.pruning_power;
//					method_file_sum_deviation_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.output_argument.sum_deviation;//191206
//					method_file_accuracy_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.result_accuracy;//191204 accuracy
//					method_file_run_time_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.whole_run_time;
//					method_file_approximation_time_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.representation_time;
//					method_file_knn_time_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.knn_total_time;
//
//					/*-----------------------------------------------------------------------------------------------------------------------------------------*/
//					//ttool.getMultiFoldToSingleByID(multi_file_pointer,3,10,0, test_time_series);
//					//ttool.printArray(test_time_series,n1);
//					//system("pause");
//				}
//				change_file = false;
//			}
//		}
//
//		/*----------------------read file address vector---------------------------------*/
//		data_source.read_file_address_vector.clear();
//		data_source.read_file_address_vector.shrink_to_fit();
//		/*--------------------------------------------------------------------------------*/
//		/*----------------------Clear Y projection  Memory---------------------------------*/
//		for (auto&& au : multi_y_projection_argument_vector) {
//			au.whole_difference_map.clear();
//		}
//		multi_y_projection_argument_vector.clear();
//		multi_y_projection_argument_vector.shrink_to_fit();
//
//		for (auto&& au : multi_all_linked_list) {
//			au.clear();
//		}
//		multi_all_linked_list.clear();
//		multi_all_linked_list.shrink_to_fit();
//
//		for (auto&& au : multi_cluster_linked_list) {
//			au.clear();
//		}
//		multi_cluster_linked_list.clear();
//		multi_cluster_linked_list.shrink_to_fit();
//		/*--------------------------------------------------------------------------------*/
//
//		delete[] query_time_series;
//		query_time_series = nullptr;
//	}
//	//}
//
//	/*============================================Print Write result========================================================*/
//	//1 MSPLA, 2 PLA, 3APCA, 4 PAA, 5 Chebyshev, 6 ICDE07
//	cout << "Method Sum Prune Power: \n";
//	cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07\n";
//	for (auto&& au : method_sum_prune_power_vector) {
//		cout << au << ", ";
//	}
//	cout << endl;
//	cout << "Method Sum Deviation: \n";
//	cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07\n";
//	for (auto&& au : method_sum_sum_deviation_vector) {
//		cout << au << ", ";
//	}
//	cout << endl;
//	cout << "Method Sum Accuracy: \n";
//	cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07\n";
//	for (auto&& au : method_sum_accuracy_vector) {
//		cout << au << ", ";
//	}
//	cout << endl;
//	cout << "Apporximation Sum Time: \n";
//	cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07\n";
//	for (auto&& au : method_approximation_time_vector) {
//		cout << au << ", ";
//	}
//	cout << endl;
//	cout << "KNN Sum Time: \n";
//	cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07\n";
//	for (auto&& au : method_knn_time_vector) {
//		cout << au << ", ";
//	}
//	cout << endl;
//	cout << "Sum Time: \n";
//	cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07\n";
//	for (auto&& au : method_run_time_vector) {
//		cout << au << ", ";
//	}
//	cout << endl;
//	/*=======================================================================================================================*/
//	/*-------------------------------------------Matlab Write result--------------------------------------------------------*/
//	K_total_number = K_coefficient_vector.size();
//	int total_size = K_total_number * N_coefficient_vector.size() * data_source.file_operation_number * method_total_number;
//	assert(total_prune_power_vector.size() == total_run_time_vector.size() && total_prune_power_vector.size() == total_accuracy_vector.size() && total_prune_power_vector.size() == total_size);
//	assert(K_total_number != INF, !N_coefficient_vector.empty(), data_source.file_total_number != INF, method_total_number != INF, data_source.point_number != INF, query_time_series_id != INF, max_node0 != INF);
//	vector<int> barchart_coefficient_vector = { K_total_number, int(N_coefficient_vector.size()), data_source.file_operation_number, method_total_number, data_source.point_number, query_time_series_id, max_node0, n0, data_source.data_type };
//
//	const string& const str_suffix = "200101";
//	const string& const str_1 = "./191120KNNMatlab/BarChartCoefficient" + str_suffix; //191122 chart coefficient
//	const string& const str_2 = "./191120KNNMatlab/TotalPrunePower" + str_suffix;//191122 prune power
//	const string& const str_19 = "./191120KNNMatlab/TotalSumDeviation" + str_suffix;//191204 sum deviation
//	const string& const str_18 = "./191120KNNMatlab/TotalAccuracy" + str_suffix;//191205 Accuracy
//	const string& const str_3 = "./191120KNNMatlab/TotalKNNRunTime" + str_suffix;//191122 run time
//	const string& const str_4 = "./191120KNNMatlab/KCoefficient" + str_suffix;//191122 K records. eg. K = 2, 4 , 8
//	const string& const str_11 = "./191120KNNMatlab/NCoefficient" + str_suffix;//191128 N records. eg. N = 18, 36, 72, 144, 288
//	const string& const str_5 = "./191120KNNMatlab/ReadMe" + str_suffix;//191122 
//	const string& const str_6 = "./191120KNNMatlab/MethodPrunePower" + str_suffix;////191123 size is total method number. prune power of every method for matlab barchart
//	const string& const str_20 = "./191120KNNMatlab/MethodSumDeviation" + str_suffix;//191206
//	const string& const str_16 = "./191120KNNMatlab/MethodAccuracy" + str_suffix;////191204 Accuracy of prune power
//	const string& const str_7 = "./191120KNNMatlab/MethodKNNRunTime" + str_suffix;////191123 size is total method number. total run time of every method for matlab barchart
//	const string& const str_12 = "./191120KNNMatlab/MethodApproximationTime" + str_suffix;////191123 size is total method number. total run time of every method for matlab barchart
//	const string& const str_13 = "./191120KNNMatlab/MethodOnlyKNNTime" + str_suffix;////191123 size is total method number. total run time of every method for matlab barchart
//	const string& const str_8 = "./191120KNNMatlab/FileMethodPrunePower" + str_suffix;////191125 size is total file number. prune power of every file for matlab barchart
//	const string& const str_21 = "./191120KNNMatlab/FileMethodSumDeviation" + str_suffix;//191206
//	const string& const str_17 = "./191120KNNMatlab/FileMethodAccuracy" + str_suffix;////191204 File Accuracy of prune power
//	const string& const str_9 = "./191120KNNMatlab/FileMethodKNNRunTime" + str_suffix;////191125 size is total file number. total run time of every file for matlab barchart
//	const string& const str_14 = "./191120KNNMatlab/FileMethodApproximationTime" + str_suffix;////191125 size is total file number. total run time of every file for matlab barchart
//	const string& const str_15 = "./191120KNNMatlab/FileMethodOnlyKNNTime" + str_suffix;////191125 size is total file number. total run time of every file for matlab barchart
//	const string& const str_10 = "./191120KNNMatlab/StringSuffix";////191126 prefix of file,
//
//	const string& const read_me = "Data Type: " + to_string(data_source.data_type) + "Data Attribute: " + dataset_attribute[data_source.data_type] + "Total file number: " + to_string(data_source.file_operation_number) + ". Total N number: " + to_string(N_coefficient_vector.size()) + ". Total Method number: " + to_string(method_total_number)
//		+ ". Total K number: " + to_string(K_total_number) + ". point number: " + to_string(data_source.point_number)
//		+ ". query time series id: " + to_string(query_time_series_id) + ". Rtree  max sub node: " + to_string(max_node0) + ". data dimension: " + to_string(data_source.time_series_dimension) + ". mixed dataset: " + TOOL::convert_vector_to_string(data_source.data_list)
//		+ ". N coefficient: " + TOOL::convert_vector_to_string(N_coefficient_vector)
//		+ ". K coefficient: " + TOOL::convert_vector_to_string(K_coefficient_vector)
//		+ ". ***Notice:.Test Multi dimension data set. MinMax point MBR for SAPLA. Single dataset.";
//
//	TOOL::writeSingleResult(str_1, barchart_coefficient_vector);
//	TOOL::writeSingleResult(str_2, total_prune_power_vector);
//	TOOL::writeSingleResult(str_19, total_sum_deviation_vector);//191206
//	TOOL::writeSingleResult(str_18, total_accuracy_vector);//191205
//	TOOL::writeSingleResult(str_3, total_run_time_vector);
//	TOOL::writeSingleResult(str_4, K_coefficient_vector);
//	TOOL::writeSingleResult(str_11, N_coefficient_vector);
//	TOOL::writeSingleResultTimeStamp(str_5, read_me);
//	TOOL::writeSingleResult(str_6, method_sum_prune_power_vector);
//	TOOL::writeSingleResult(str_20, method_sum_sum_deviation_vector);//191206
//	TOOL::writeSingleResult(str_16, method_sum_accuracy_vector);//191204
//	TOOL::writeSingleResult(str_7, method_run_time_vector);
//	TOOL::writeSingleResult(str_12, method_approximation_time_vector);
//	TOOL::writeSingleResult(str_13, method_knn_time_vector);
//	TOOL::writeSingleResult(str_8, method_file_prune_power_vector);
//	TOOL::writeSingleResult(str_21, method_file_sum_deviation_vector);//191206
//	TOOL::writeSingleResult(str_17, method_file_accuracy_vector);//191204
//	TOOL::writeSingleResult(str_9, method_file_run_time_vector);
//	TOOL::writeSingleResult(str_14, method_file_approximation_time_vector);
//	TOOL::writeSingleResult(str_15, method_file_knn_time_vector);
//	TOOL::writeSingleResult(str_10, str_suffix);
//	/*----------------------------------------------------------------------------------------------------------------------*/
//
//	/*-----------------------------------------------Clear Memory-----------------------------------------------------------------------*/
//	N_coefficient_vector.clear();
//	N_coefficient_vector.shrink_to_fit();
//	K_coefficient_vector.clear();
//	K_coefficient_vector.shrink_to_fit();
//	total_prune_power_vector.clear();//191119
//	total_prune_power_vector.shrink_to_fit();
//	total_sum_deviation_vector.clear();//191206
//	total_sum_deviation_vector.shrink_to_fit();
//	total_accuracy_vector.clear();//191204
//	total_accuracy_vector.shrink_to_fit();
//	total_run_time_vector.clear();//191119
//	total_run_time_vector.shrink_to_fit();
//	method_sum_prune_power_vector.clear();//191123 totl prune power of every method for matlab barchart
//	method_sum_prune_power_vector.shrink_to_fit();
//	method_sum_sum_deviation_vector.clear();
//	method_sum_sum_deviation_vector.shrink_to_fit();
//	method_sum_accuracy_vector.clear();//191204 accuracy of prune power
//	method_sum_accuracy_vector.shrink_to_fit();
//	method_run_time_vector.clear();//191123 totl run time of every method for matlab barchart
//	method_run_time_vector.shrink_to_fit();
//	method_approximation_time_vector.clear();//191204 For approximation time of every method
//	method_approximation_time_vector.shrink_to_fit();
//	method_knn_time_vector.clear();//191204 For knn time of every method
//	method_knn_time_vector.shrink_to_fit();
//
//	method_file_prune_power_vector.clear();//191125 5*23 sum prune power for every file [mehtod number, file number]
//	method_file_prune_power_vector.shrink_to_fit();
//	method_file_sum_deviation_vector.clear();
//	method_file_sum_deviation_vector.shrink_to_fit();
//	method_file_accuracy_vector.clear();//191204 accuracy of prune power
//	method_file_accuracy_vector.shrink_to_fit();
//	method_file_run_time_vector.clear();//191125 5*23 sum run time for every file [mehtod number, file number]
//	method_file_run_time_vector.shrink_to_fit();
//	method_file_approximation_time_vector.clear();//191204 5*23 sum approximation time for every file [mehtod number, file number]
//	method_file_approximation_time_vector.shrink_to_fit();
//	method_file_knn_time_vector.clear();//191204 5*23 sum knn time for every file [mehtod number, file number]
//	method_file_knn_time_vector.shrink_to_fit();
//	/*--------------------------------------------------------------------------------------------------------------------------*/
//
//	system("pause");
}

//***************************************************************
// Method:evaluate_multi_KNN_speed
// Qualifier:  evaluate KNN search for multi & single time series of APLA PLA, APCA, PAA & ICDE07_APLA
// Input:
// Output: sum deviaiton, max deviation, time , I/O cost and so on.
// Note: change the order of loop to speed up algorithm
// date:200101
// author:
//***************************************************************
TEMPLATE
template<typename T>
void Evaluation::evaluate_multi_KNN_speed(const T& const evaluation_argument_struct) {

	//string file_address = "DataAddress181218.txt";
	//string read_file_name = "CinC_ECG_torso_TEST";
	//string chebyshev_write_name[] = { "181206Chebyshev prun power","181206Chebyshev IO cost","181206Chebyshev buld index time","181206Chebyshev navigate time", "181206Chebyshev LB distance","181206Chebyshev euc distane", "181206Chebyshev total KNN time","181206Chebyshev result accuracy" };
	string* chebyshev_write_pointer = nullptr;//chebyshev_write_name;
	string dataset_attribute[] = { "", "Single dataset", "Mixed dataset", "Multi dimension data set" };

	//bool change_file = false;
	int query_time_series_id = INF;
	

	int K = INF;// for KNN, the number of nearest neighbor for query time series

	/*----   201028 Probability of max deviation points are minmax, left&right endpoints   ----*/
	long double number_point_max_deviation_true = 0;
	long double number_point_max_deviation_false = 0;
	long double number_not_smaller_than_sum_deviation = 0;
	long double number_smaller_than_sum_deviation = 0;
	long double number_not_smaller_than_sum_deviation_pow = 0;
	long double number_smaller_than_sum_deviation_pow = 0;
	/*-----------------------------------------------------------------------------------------*/

	/*=========================================================         file id choose        =====================================================*/
	int initial_file_od = INF;//191124
	int file_id = INF; // the id of file(data set / file)
	vector<size_t> file_id_choose_vector;//201020 store specific file id
	vector<size_t> data_list;//210603 mixed dataset id
	/*============================================================================================================================================*/

	/*=================================================         Representation method Title vector       ==========================================*/

	size_t representation_option = INF;// 1 MSPLA, 2 PLA, 3 APCA, 4 PAA, 5 Chebyshev, 6 ICDE07, 7 PAALM, 8 Initial_200706, 9 SAX, 10 Linear Scan

	// vector<int> representation_option_vector{ int(INF), 1, 2, 3, 4, 5, 6, 7, 8, 9};
	// vector<int> representation_option_vector{ int(INF), 1, 2, 3, 4, 5, 7, 8, 9 };
	//vector<int> representation_option_vector{ int(INF), 2, 4, 7, 8, 9, 10 };// No 2 PLA 3 APCA, 5Ceby, 6 ICDE07
	//vector<int> representation_option_vector{ int(INF), 2, 3, 4, 5, 7, 8, 9, 10 };// No ICDE07
	
	//vector<int> representation_option_vector{ int(INF), 2, 4, 5, 7, 8, 9 };// 2 PLA, 4 PAA, 5 Chebyshev, 7 PAALM, 8 Initial_200706, 9 SAX
	//vector<int> representation_option_vector{ int(INF), 2, 3, 4, 5, 6, 7, 8, 9, 10 };// has SAPLA
	
	// vector<int> representation_option_vector{ int(INF), 2};// PLA
	// vector<int> representation_option_vector{ int(INF), 3};// APCA
	// vector<int> representation_option_vector{ int(INF), 4};// PAA
	// vector<int> representation_option_vector{ int(INF), 5};// CHEBY
	// vector<int> representation_option_vector{ int(INF), 6};// 6 ICDE07
	//vector<int> representation_option_vector{ int(INF), 8 };//Only SAPLA
	// vector<int> representation_option_vector{ int(INF), 9 };//Only SAX
	// vector<int> representation_option_vector{ int(INF), 10 };//Only Linear Scan
	// vector<int> representation_option_vector{ int(INF), 2, 8 };//PLA SAPLA
	// vector<int> representation_option_vector{ int(INF), 8, 10 };//SAPLA LinearScan
	// vector<int> representation_option_vector{ int(INF), 2, 8, 10 };//PLA SAPLA SAPLA LinearScan
	// vector<int> representation_option_vector{ int(INF), 2, 3, 10 };//PLA APCA SAPLA
	// vector<int> representation_option_vector{ int(INF), 4, 7 };//PAA PAALM
	// vector<int> representation_option_vector{ int(INF), 6, 8 };//APLA, SAPLA
	// vector<int> representation_option_vector{ int(INF),  3, 4, 7 };// APCA PAA PAALM
	//vector<int> representation_option_vector(evaluation_argument_struct.representation_option_vector.begin(), evaluation_argument_struct.representation_option_vector.end());
	vector<int> representation_option_vector{ int(INF), 2, 3, 4, 5, 6, 7, 8, 9, 11 };

	const vector<string> name_representation_vector_string_all_vector{ "MSPLA", "PLA", "APCA", "PAA", "Chebyshev", "APLA", "PAALM", "SAPLA", "SAX", "LTDB", "LinearScan"};
	vector<string> name_representation_vector_string_vector;
	for (int i = 1; i < representation_option_vector.size(); i++) {
		name_representation_vector_string_vector.emplace_back(name_representation_vector_string_all_vector[representation_option_vector[i] - 1]);
	}
	/*============================================================================================================================================*/

	int N = INF;//N segment number. degree_m n=m+1. APLA = 4, PAA : 12 , PLA: 6, APCA 6. Cheby 12, ICDE07: 4
	int initial_N = INF;// for y projection method in APLA, cannot be NULL;
	int final_N = INF;
	int initial_K = INF;//191204 for Yprojction
	int final_K = INF;
	int data_final_dimension = INF;// 1 or 3. Final dimension of data, is single dimension or multi dimension data
	int arity_d0 = INF;//dimension of time series, most are one dimension
	int max_node0 = INF;
	int file_total_number = INF;
	int n0 = INF;// time series length single is 512, multi dimenson is 300
	int point_number0 = INF;

	/*=========================200130 Figure Type========================================*/
	//200325 experiment of initial_N with accuracy & time.
	const int& const evaluation_type = 0;//0: for different methods comparison. 1: for differnt initial_N comparison
	/*===================================================================================*/
	/*=========================      Data Type   ========================================*/
	const size_t& const data_type = 1; //1 single dataset, 2 mixed dataset, 3 multi dimension data set, 4 each step in SAPLA
	/*===================================================================================*/

	/*-------------------------  210603 Write Mixed DataSet  -----------------------------------*/
	string data_file_name = "./191202DataSet/data_set.txt";
	/*-----------------------------------------------------------------------------------*/
	string* file_address_pointer = nullptr;

	/*------------------------------------Multi Dimension Dataset----------------------------------------*/
	/*string multi_file_name1[] = { "./191202DataSet/MultiDateSet191220/Cricket_X_TEST","./191202DataSet/MultiDateSet191220/Cricket_Y_TEST","./191202DataSet/MultiDateSet191220/Cricket_Z_TEST" };
	string multi_file_name2[] = { "./191202DataSet/MultiDateSet191220/uWaveGestureLibrary_X_TEST","./191202DataSet/MultiDateSet191220/uWaveGestureLibrary_Y_TEST","./191202DataSet/MultiDateSet191220/uWaveGestureLibrary_Z_TEST" };
	string multi_dimension_file_name[] = {
		"Cricket_X_TEST","Cricket_Y_TEST","Cricket_Z_TEST",
		"uWaveGestureLibrary_X_TEST", "uWaveGestureLibrary_Y_TEST", "uWaveGestureLibrary_Z_TEST",
	};*/
	//string multi_to_single_file_name = "./191202DataSet/MultiToSingle191220/uWaveGestureLibrary";
	//191220 project multi data set to single data set
	//TOOL::project_multi_data_to_single_data(multi_file_pointer, 3, 3582, 315, multi_to_single_file_name);
	//TOOL::project_single_data_to_multi_data_batch(TOOL::single_to_multi_file_address, 13, 3, 40, 600, TOOL::write_single_to_multi_file_address);
	/*---------------------------------------------------------------------------------------------------*/

	/*#############################################      Coefficients of Write each file each approximation methods       ####################################*/
	const string& const str_suffix = evaluation_argument_struct.str_suffix;// "211223a";0

	int tree_type = evaluation_argument_struct.tree_type;//0;// 0 Rtree, 1//1

	/* 0 has burst data 21 ; 1 no burst data 18; 2 all 24 datasets; 3 201016 UCR2018 512 41 datasets; 4 201019 UCR2018 1024 21 datasets;
	// 5 201020 UCR2018 choose specific datasets; 6 Test SAPLA distance; */
	const int option_homogenous_data_type = evaluation_argument_struct.option_homogenous_data_type;//7;// default is 5 ,6 ,7//2
	int size_file = evaluation_argument_struct.size_file;//1; //3
	int file_id_begin = evaluation_argument_struct.file_id_begin;//1; //4
	const int size_query_time_series = evaluation_argument_struct.size_query_time_series;//2; //5
	int temp_n0 = evaluation_argument_struct.n;//1024; // 6
	int temp_point_number0 = evaluation_argument_struct.point_number;//50;//default number: 50 //11
	int temp_max_node0 = evaluation_argument_struct.max_node;//5; //12
	int temp_initial_N = evaluation_argument_struct.initial_N;//12;//default 12 // 7
	int temp_final_N = evaluation_argument_struct.final_N;//20; // 8
	int temp_initial_K = evaluation_argument_struct.initial_K;//1; // 9
	int temp_final_K = evaluation_argument_struct.final_K;//6; // 10
	
	/*====================  210906 Rtree / Partition Tree  ==============================*/
	typename TOOL::template OPTION_TREE<int, int, int> option_tree_struct(tree_type, 0, 0, 0, 0);// 0: (0 R-Tree, 1 Partition Tree, 2 SAPLA MBR)
	/*===================================================================================*/

	string notice_string = ".*Notice: suffix: " + str_suffix + ". Data type: " + to_string(data_type) + " Tree type: " + to_string(option_tree_struct.type_tree) 
		+ "Test 89 files";

	//const string& const str_33 = "./191120KNNMatlab/EachFileMethodPrunePowerCombine" + str_suffix;//201221 size is total file number. prune power of every file for matlab barchart
	//const string& const str_27 = "./191120KNNMatlab/EachFileMethodPrunePower" + str_suffix;//201028 size is total file number. prune power of every file for matlab barchart
	//const string& const str_28 = "./191120KNNMatlab/EachFileMethodSumDeviation" + str_suffix;//201028
	//const string& const str_29 = "./191120KNNMatlab/EachFileMethodAccuracy" + str_suffix;//201028 File Accuracy of prune power
	//const string& const str_30 = "./191120KNNMatlab/EachFileMethodKNNRunTime" + str_suffix;//201028 size is total file number. total run time of every file for matlab barchart
	//const string& const str_31 = "./191120KNNMatlab/EachFileMethodApproximationTime" + str_suffix;//201028 size is total file number. total run time of every file for matlab barchart
	//const string& const str_32 = "./191120KNNMatlab/EachFileMethodOnlyKNNTime" + str_suffix;//201028 size is total file number. total run time of every file for matlab barchart

	
	///************** 201028 Clear Result for one File  ****************/
	const bool clear_file = true;// false true
	//if (clear_file) {
	//	TOOL::clearTXTFile(str_33);// Prune Power Combine
	//	TOOL::clearTXTFile(str_27);// Prune Power
	//	TOOL::clearTXTFile(str_28);// SumDeviation
	//	TOOL::clearTXTFile(str_29);// Accuracy
	//	TOOL::clearTXTFile(str_30);// KNNRunTime
	//	TOOL::clearTXTFile(str_31);// ApproximationTime
	//	TOOL::clearTXTFile(str_32);// OnlyKNNTime
	//}
	///******************************************************************/
	/*########################################################################################################################################################*/

	/*###################################################      DataType: homogerous heterogerous data      ###########################################################*/


	/*******************************************************************************************************************************************************************************/

	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	switch (data_type) {
	case 0: // 0 manipulate data set;    
		assert(0);
		break;
	case 1:// 1 real data set, homogenous(single) data;

		n0 = 512;
		initial_N = 12;// 201222 Default single, mixed data? 12, multi data: 24
		final_N = 192;
		point_number0 = 20;
		final_K = 16;

		max_node0 = 4;//3 for 20 point numebr .for Rtree node, most branches in one node.  datatype: 1, 3

		/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
		switch (option_homogenous_data_type) {
		case 0: // 0 has burst data 21
			file_total_number = 21; //file has 21 files for DistLB.

			/*--201020 Choose File id : Fills the range [first, last) with sequentially increasing values--*/
			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			/*---------------------------------------------------------------------------------------------*/

			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;//19 300 195

			break;
		case 1: // 1 no burst data 18
			file_total_number = 18;//No burst time series, has 18 files.

			/*--201020 Choose File id : Fills the range [first, last) with sequentially increasing values--*/
			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			/*---------------------------------------------------------------------------------------------*/

			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;//19 300 195

			break;
		case 2: // has all 24 files
			file_total_number = 24;//All 24 files.
			/*--201020 Choose File id : Fills the range [first, last) with sequentially increasing values--*/
			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			/*---------------------------------------------------------------------------------------------*/
			break;
		case 3: // 201016 UCR2018 512 41 datasets;
			file_total_number = 41;//All 41 files.
			/*--201020 Choose File id : Fills the range [first, last) with sequentially increasing values--*/
			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			/*---------------------------------------------------------------------------------------------*/

			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;//19 300 195

			break;
		case 4: // 201019 UCR2018 1024 21 datasets;
			file_total_number = 21;//All 21 files.
			n0 = 1024;
			final_N = 384;
			point_number0 = 50;
			final_K = 32;

			/*--201020 Choose File id : Fills the range [first, last) with sequentially increasing values--*/
			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);// begin from 0;
			/*---------------------------------------------------------------------------------------------*/

			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;//19 300 195

			break;
		case 5: // 201020 UCR2018 choose specific datasets; 1-21

			/*-------  201020 Choose File id : Fills the range [first, last) with sequentially increasing values  --------*/
			file_id_choose_vector.resize(size_file);// Default 20 files that data number >= 100
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), file_id_begin);// Default 1. begin from 1, not 0;
			/*------------------------------------------------------------------------------------------------------------*/

			file_total_number = file_id_choose_vector.size();
			n0 = temp_n0;// Default 1024
			initial_N = temp_initial_N;// Default 12
			final_N = temp_final_N;//Default 384
			point_number0 = temp_point_number0;//default : 100
			final_K = temp_final_K;//Default 64
			max_node0 = temp_max_node0;// Default 4

			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;//19 300 195
			query_time_series_id = 0;

			break;
		case 6: // 210305 just for each step in SAPLA. UCR2018 21 files

			/*-------  201020 Choose File id : Fills the range [first, last) with sequentially increasing values  --------*/
			file_id_choose_vector.resize(size_file);//defualt 21. 20 files that data number >= 100
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), file_id_begin);// begin from 0
			/*------------------------------------------------------------------------------------------------------------*/

			file_total_number = file_id_choose_vector.size();
			n0 = temp_n0;// Default 1024
			initial_N = temp_initial_N;// Default 12
			final_N = temp_final_N;// Default 384
			point_number0 = temp_point_number0;//Default number: 100
			final_K = temp_final_K;// Default 64
			max_node0 = temp_max_node0;// Default 4

			query_time_series_id = 0;//19 300 195

			break;
		case 7: // 211220 just for each step in SAPLA. UCR2018 89 files n <= 128, points number > 50

		    /*-------  211220 Choose File id : Fills the range [first, last) with sequentially increasing values  --------*/
			file_id_choose_vector.resize(size_file);//defualt 89. 89 files that data number >= 100
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), file_id_begin);// begin from 0
			/*------------------------------------------------------------------------------------------------------------*/

			file_total_number = file_id_choose_vector.size();
			n0 = temp_n0;// default 128
			initial_N = temp_initial_N;//default 12
			final_N = temp_final_N;//defualt 48
			point_number0 = temp_point_number0;//default number: 50
			final_K = temp_final_K;//defualt 32
			max_node0 = temp_max_node0;// defualt 4

			query_time_series_id = 0;//19 300 195

			break;
		default:
			assert(0);
		}
		/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

		data_final_dimension = 1;
		assert(data_final_dimension == 1);
		file_address_pointer = nullptr;//210603
		file_address_pointer = &data_file_name;//210603
		/*---------------------------         210603    Mixed dataset id             ---------------------*/
		for (int list_id = 1; list_id < file_total_number + 1; list_id++) { data_list.emplace_back(list_id); }// do not first file
		/*------------------------------------------------------------------------------------------------*/
		break;
	case 2:// 2 real data set, heterogeneous(mixed) data.
		n0 = 1024;
		initial_N = 12;//single, mixed data£º 12, multi data: 24
		final_N = 192;//Default 384
		point_number0 = 20;// for each file
		final_K = 128;// Default 128.
		max_node0 = 5;//3 5 for Rtree node, most branches in one node.  datatype: 1, 3
		file_total_number = 20;//210107 except  except 5 10 11 12 18. //21 //200126file has 20 files. 200116file has 21 files for DistLB
		data_final_dimension = 1;
		assert(data_final_dimension == 1);
		file_address_pointer = &data_file_name;//210603

		file_id_choose_vector.resize(1, 0);

		/*---------------------------         210603    Mixed dataset id             ---------------------*/
		for (int list_id = 1; list_id < file_total_number + 1; list_id++) { data_list.emplace_back(list_id); }// do not first file
		/*----------------------------------------------------------------------------------------------*/

		query_time_series_id = TOOL::get_random_max(point_number0 * file_total_number);
		query_time_series_id = 359;//+20         19 300 195

		break;
	case 3:// 3 multi dimension homogenous data.
		//file_address_pointer = multi_dimension_file_name;
		n0 = 295;//200
		initial_N = 12;// Not multiple 3, this is final dimension for extended time serires length:= length * d .single, mixed data£º 12, multi data: 24
		final_N = 192;
		point_number0 = 300;//40
		final_K = 128;//64,128
		max_node0 = 10;//3 5 for Rtree node, most branches in one node.  datatype: 1, 3

		/*-------  201020 Choose File id : Fills the range [first, last) with sequentially increasing values  --------*/
		file_id_choose_vector.resize(2);//20 files that data number >= 100
		std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);//
		/*------------------------------------------------------------------------------------------------------------*/

		data_final_dimension = 3;
		file_total_number = file_id_choose_vector.size() * data_final_dimension;// single dimension 24, 3 dimension data set 13 * 3 = 39

		file_address_pointer = nullptr;
		query_time_series_id = TOOL::get_random_max(point_number0);
		query_time_series_id = 89;//269 249 219
		assert(data_final_dimension == 3);

		break;
	default:
		assert(0);
		break;
	}
	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	notice_string += "Final N: " + to_string(final_N) + ". Final K: " + to_string(final_K);
	/*:::::::::::::::::::::::::::::::::::::*/
	initial_file_od = 0;// initial is 0
	//file_total_number = 1;
	const int final_file_id = file_id_choose_vector.size();// 210602 Default = file_id_choose_vector.size();//= file_total_number;//==  Default value = data_source.file_operation_number
	//const int final_file_id = 1;
	initial_K = temp_initial_K;//default 1 
	//int K_total_number = INF;
	//if (option_homogenous_data_type == 6) initial_K = 2;
	/*:::::::::::::::::::::::::::::::::::::*/

	/*##################################################################################################################################################################*/

	/*########################################################               Matlab Bar chart              #############################################################*/

	const int method_total_number = representation_option_vector.size() - 1;// 1 MSPLA, 2 PLA, 3 APCA, 4 PAA, 5 Chebyshev, 6 ICDE07, 7 PAALM, 8 Initial_200706, 9 SAX

	/*----------------------------------Figure_type == 1---------------------------------*/
	//200325
	//if (evaluation_type == 1) {
	//	method_total_number = 10; //200130 SAPLA initial_N
	//}
	/*-----------------------------------------------------------------------------------*/

	vector<size_t> N_coefficient_vector;//191128  N coefficient for Matlab Barchart
	for (int N = initial_N; N <= final_N; N += 6) {//  += 6 N  5-80; 6 96 - 192;384 2;191119 12 24 48 96 192 258 384;;; 192 96 24, 
		N_coefficient_vector.emplace_back(N);
	}

	vector<int> K_coefficient_vector;//191120 K coefficient for Matlab Barchart
	for (int K = initial_K; K <= final_K; K *= 2) {//K 2 4 8 16 32 64 128 ;;; 16 128
		/*----------------------------------------*/
		 K_coefficient_vector.emplace_back(K);
		/*----------------------------------------*/
	}

	/*******************  Each Detail results   **********************/
	vector<long double> total_prune_power_combine_vector;//201221 old p / accuracy
	vector<long double> total_prune_power_vector;//191119
	vector<long double> total_sum_deviation_vector;//191206
	vector<long double> total_max_deviation_vector;//191206
	vector<long double> total_max_deviation_av_vector;//210910
	vector<long double> total_max_width_deviation_vector;//191206
	vector<long double> total_accuracy_vector;//191204
	vector<long double> total_run_time_vector;//191119
	vector<long double> total_run_time_has_IO_vector;//210606
	vector<long double> total_approximation_time_vector;//200108
	vector<long double> total_build_tree_time_vector;//210906
	vector<long double> total_ingest_data_time_vector;//211210
	vector<long double> total_knn_time_vector;//200108
	vector<long double> total_knn_time_has_IO_vector;//210606
	vector<long double> total_IO_vector;//211202
	/***************************************************************/

	/*********************************           For Each Mehtod            ****************************************/
	vector<long double> method_sum_prune_power_combine_vector(method_total_number, 0);//201221 old p / accuracy
	vector<long double> method_sum_prune_power_vector(method_total_number, 0);//191123 totl prune power of every method for matlab barchart
	vector<long double> method_sum_sum_deviation_vector(method_total_number, 0);//191206 sum deviation vector
	vector<long double> method_sum_max_deviation_vector(method_total_number, 0);//191206 max deviation vector
	vector<long double> method_sum_max_deviation_av_vector(method_total_number, 0);//210911 max deviation av vector
	vector<long double> method_sum_max_width_deviation_vector(method_total_number, 0);//191206 sum deviation vector
	vector<long double> method_sum_accuracy_vector(method_total_number, 0);//191204 accuracy of prune power
	vector<long double> method_run_time_vector(method_total_number, 0);//191123 totl run time of every method for matlab barchart
	vector<long double> method_run_time_has_IO_vector(method_total_number, 0);//210606 totl run time of every method for matlab barchart
	vector<long double> method_approximation_time_vector(method_total_number, 0);//191204 For approximation time of every method
	vector<long double> method_build_tree_time_vector(method_total_number, 0);//191204 For approximation time of every method
	vector<long double> method_ingest_data_time_vector(method_total_number, 0);//211211 For approximation time of every method
	vector<long double> method_knn_time_vector(method_total_number, 0);//191204 For knn time of every method
	vector<long double> method_knn_time_has_IO_vector(method_total_number, 0);//210606 For knn time of every method
	vector<long double> method_IO_vector(method_total_number, 0);//211202 For I/O cost of every method
	/*************************************************************************************************************/

	/******************************      201028 Write Result by One File    ***************************************/
	vector<long double> one_file_prune_power_combine_vector(method_total_number, 0);//201221 old p / accuracy
	vector<long double> one_file_prune_power_vector(method_total_number, 0);//191123 totl prune power of every method for matlab barchart
	vector<long double> one_file_sum_deviation_vector(method_total_number, 0);//191206 sum deviation vector
	vector<long double> one_file_max_deviation_vector(method_total_number, 0);//191206 sum deviation vector
	vector<long double> one_file_max_width_deviation_vector(method_total_number, 0);//191206 sum deviation vector
	vector<long double> one_file_accuracy_vector(method_total_number, 0);//191204 accuracy of prune power
	vector<long double> one_file_run_time_vector(method_total_number, 0);//191123 totl run time of every method for matlab barchart
	vector<long double> one_file_approximation_time_vector(method_total_number, 0);//191204 For approximation time of every method
	vector<long double> one_file_knn_time_vector(method_total_number, 0);//191204 For knn time of every method
	vector<long double> one_file_knn_time_has_IO_vector(method_total_number, 0);//210606 For knn time of every method
	/*************************************************************************************************************/
	//TOOL::writeResultNoCover(const string& const write_file_name, vector<T>& array_data)

	/**********************************      For Each method, Each File      *************************************/
	//Split as method
	vector<long double> method_file_prune_power_combine_vector(size_t(method_total_number * file_total_number), 0.0);// 201221 old p / accuracy
	vector<long double> method_file_prune_power_vector(size_t(method_total_number * file_total_number), 0.0);//191125 5*23 sum prune power for every file [mehtod number, file number]
	vector<long double> method_file_sum_deviation_vector(method_total_number * file_total_number, 0.0);//191206
	vector<long double> method_file_max_deviation_vector(method_total_number * file_total_number, 0.0);//191206
	vector<long double> method_file_max_deviation_av_vector(method_total_number* file_total_number, 0.0);//210910
	vector<long double> method_file_max_width_deviation_vector(method_total_number * file_total_number, 0.0);//191206
	vector<long double> method_file_accuracy_vector(method_total_number * file_total_number, 0.0);//191204 accuracy of prune power
	vector<long double> method_file_run_time_vector(method_total_number * file_total_number, 0.0);//191125 5*23 sum run time for every file [mehtod number, file number]
	vector<long double> method_file_run_time_has_IO_vector(method_total_number * file_total_number, 0.0);//210606 5*23 sum run time for every file [mehtod number, file number]
	vector<long double> method_file_approximation_time_vector(method_total_number * file_total_number, 0.0);//191204 5*23 sum approximation time for every file [mehtod number, file number]
	vector<long double> method_file_build_tree_time_vector(method_total_number* file_total_number, 0.0);
	vector<long double> method_file_ingest_data_time_vector(method_total_number* file_total_number, 0.0);
	vector<long double> method_file_knn_time_vector(method_total_number * file_total_number, 0.0);//191204 5*23 sum knn time for every file [mehtod number, file number]
	vector<long double> method_file_knn_time_has_IO_vector(method_total_number * file_total_number, 0.0);//210606 5*23 sum knn time for every file [mehtod number, file number]
	vector<long double> method_file_IO_vector(method_total_number * file_total_number, 0.0);//211203 I/O cost
	/*************************************************************************************************************/

	/*|||||||||||||||||||||     Assign Evaluation Result 211213   ||||||||||||||||||||||||||||||||*/
	RESULT_COEFFICIENT<int> Result_Coefficient_Struct;
	METHOD_RESULT<long double> Method_Result_Struct(method_total_number);
	METHOD_FILE_RESULT<long double> Method_File_Result_Struct(method_total_number * file_total_number);
	TOTAL_RESULT<long double> Total_Result_Struct;
	//TOTAL_RESULT<long double> Total_Result(method_total_number * file_total_number * N_coefficient_vector.size() * K_coefficient_vector.size());
	/*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

	/*++++++++++++++++++++++      evaluation of Split Point   ++++++++++++++++++++++++++++++++++++++*/
	// 1 min density,2 binary search, 3 direct intersection 4 middle point, 5 best split point method.
	// default value: 1
	const int split_methods_number = 1;// Default number : 5
	/*------------------------      Local Sum, shift, time    -------------------------------*/
	vector<long double> local_total_split_id_sum_deviation(split_methods_number, 0);//200220 sum sum deviation of four spit point finding methods
	vector<long double> local_total_split_id_shift(split_methods_number, 0);//200219 sum accracy of four split point finding methods  
	vector<long double> local_total_split_id_time(split_methods_number, 0);//200219 sum time
	/*---------------------------------------------------------------------------------------*/
	/*----------------------          Global Sum, time        -------------------------------*/
	vector<long double> global_total_knn_prune_power(split_methods_number, 0);//200227
	vector<long double> global_total_approximation_sum_deviation(split_methods_number, 0);
	vector<long double> global_total_approximation_time(split_methods_number, 0);
	/*---------------------------------------------------------------------------------------*/
	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	/*+++++++++++++++++++++++++++++   200326 when evaluation_type == 1, evaluate initial_N  +++++++++++++++++++++++++++++*/
	/*--------------------------------initial coeffcients------------------------------*/
	const double N_size = 5;//12 24 48 96 192
	const double final_initial_N_number = 30;
	vector<long double> initial_number_coefficients_vector{ -2, -1 };//200330
	//initial_number_coefficients_vector.resize(final_initial_N_number);
	for (double vector_id = 0; initial_number_coefficients_vector.size() < final_initial_N_number; vector_id += 3) {
		initial_number_coefficients_vector.emplace_back(vector_id);
	}
	assert(initial_number_coefficients_vector.size() == final_initial_N_number);
	/*---------------------------------------------------------------------------------*/
	/*--------------------------------Figure_type == 1---------------------------------*/
	vector<long double> approximation_initial_N_vector;//200330 
	/*---------------------------------------------------------------------------------*/
	/*-----------------------------total initial N number------------------------------*/
	vector<long double> total_initial_N_prune_power_vector(final_initial_N_number, 0);
	vector<long double> total_initial_N_sum_deviation_vector(final_initial_N_number, 0);
	vector<long double> total_initial_N_run_time_vector(final_initial_N_number, 0);
	vector<long double> total_initial_N_approximation_time_vector(final_initial_N_number, 0);
	vector<long double> total_initial_N_knn_time_vector(final_initial_N_number, 0);
	/*---------------------------------------------------------------------------------*/
	/*-----------------------------initial N number sort by N--------------------------*/
	vector<long double> initial_N_by_N_prune_power_vector(final_initial_N_number * N_size, 0);//initial_N * N number
	vector<long double> initial_N_by_N_sum_deviation_vector(final_initial_N_number * N_size, 0);//initial_N * N number
	vector<long double> initial_N_by_N_run_time_vector(final_initial_N_number * N_size, 0);//initial_N * N number
	vector<long double> initial_N_by_N_approximation_time_vector(final_initial_N_number * N_size, 0);//initial_N * N number
	vector<long double> initial_N_by_N_knn_time_vector(final_initial_N_number * N_size, 0);//initial_N * N number
	/*---------------------------------------------------------------------------------*/
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	/*+++++++++++    210121 Evaluate Upper Bound    +++++++++++++*/
	typename TOOL::EVALUATION_BOUND evaluation_bound_whole;
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	/*##############################################################################################################################################################*/

	/*========================================       read & write mixed data set; Initial DATA_SOURCE       ===============================================*/
	typename TOOL::DATA_SOURCE data_source(data_type, data_final_dimension, data_list, file_total_number, point_number0, n0, *file_address_pointer);// 0 manipulate data set; 1 real data set, homogenous(single) data; 2 real data set, heterogeneous(mixed) data. 3 multi single data. 4 multi mixed data
	data_source.option_has_burst_data = option_homogenous_data_type;
	TOOL::initial_data_source(data_source);
	data_source.bigger_account = 0;
	assert(data_source.point_number != INF && data_source.single_time_series_length != INF && data_source.multi_single_time_series_length != INF);
	/*====================================================================================================================================================*/

	/*================================       Query Time Series Vector       ====================================*/
	
	Result_Coefficient_Struct.id_query_series_vector = { 0,59,82,73,32,16,50,91,79,91,51,64,43,92,62,31,77,57,70,34,45,42,24,93,62,39,3,28,55,97,46,7,2,85,60,69,55,21,10,40,19,28,52,82,38,72,9,81,27,46 };//211204
	Result_Coefficient_Struct.id_query_series_vector = { 0,1 };
	TOOL::get_random_vector(size_query_time_series, point_number0, Result_Coefficient_Struct.id_query_series_vector);
	/*-------  211220 Choose File id : Fills the range [first, last) with sequentially increasing values  --------*/
	/*------------------------------------------------------------------------------------------------------------*/

	cout << "Query Time Series ID: "<< endl;
	TOOL::print_vector(Result_Coefficient_Struct.id_query_series_vector);

	vector<vector<DataType>> query_series_vector_vector(size_query_time_series);//vector of time series vector
	/*==========================================================================================================*/

	/*.......................................................    Dimension    ............................................................................*/
	//for (arity_d0 = 1; arity_d0 <= data_final_dimension; arity_d0 += 2) {//dimension fo time seires. 1:one dimension, 2, two dimension. 3: three dimension. 
	/*....................................................................................................................................................*/

	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      File ID      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	for (file_id = initial_file_od; file_id < final_file_id; file_id++) {//file 0-24 file_total_number = 24, 21 or 18
		cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%   EVALUATION    file id: " << file_id_choose_vector[file_id] + 1 << "   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
		bool change_file = true;

		//data_source.point_number = data_source.time_series_number_vector[file_id];
		//query_time_series_id = data_source.point_number - 1;
		//max_node0 = data_source.point_number / 10;

		/*||||||||||||||||||||||    Normalize All Original Datasets     |||||||||||||||||||||||||||||||||*/
		//Normalze all original time series and write all in <data_set.txt>. Miss first point?
		if (data_source.data_type == 1) {//210603  1 is homogeneous datasets
			//TOOL::get_normalize_write_one_file(data_source, file_id_choose_vector[file_id]);// Old version
			TOOL::mannually_normalize_write_one_file(data_source, file_id_choose_vector[file_id]);// 211230 mannually connect
		}
		/*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

		/*-----------------------------    Original time series   ------------------------------------------------*/
		bool print_each_result = false;
		if (print_each_result == true) {
			vector<DataType> original_time_series_vector_print(data_source.multi_single_time_series_length, INF);
			TOOL::get_read_multi_file_address(data_source, file_id_choose_vector[file_id]);//210603
			//already normalized in below function
			for (int i = 0; i < data_source.point_number; i++) {
				TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, i, original_time_series_vector_print);
				cout << "  " << i << " Original time series:" << endl;
				TOOL::print_vector(original_time_series_vector_print);//190501
				original_time_series_vector_print.clear();
				original_time_series_vector_print.shrink_to_fit();
			}

			data_source.read_file_address_vector.clear();
			data_source.read_file_address_vector.shrink_to_fit();
		}
		/*--------------------------------------------------------------------------------------------------------*/

		/*=============================================    Query time series   ==================================================================*/
		//vector<DataType> query_time_series_vector(data_source.multi_single_time_series_length, INF);
		//DataType* query_time_series = new DataType[data_source.multi_single_time_series_length];
		TOOL::get_read_multi_file_address(data_source, file_id_choose_vector[file_id]);//210603
		////get data from new file. already normalized in advanced, no miss first point.
		//TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, query_time_series_id, query_time_series_vector);
		//copy_n(query_time_series_vector.begin(), query_time_series_vector.size(), query_time_series);
		//211205 query time series vector
		for (int order_query = 0; order_query < size_query_time_series ; order_query++) {
			TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, Result_Coefficient_Struct.id_query_series_vector[order_query], query_series_vector_vector[order_query]);
		}
		/*========================================================================================================================================*/

		/*===============================================  Base KNN Linear Scan  =================================================================*/
		typename TOOL::template LINEAR_SCAN_STRUCTURE<long double> linear_scan_struct(0, 0, 0);
		//multiset<pair<double, int>> squential_scan_result_set;
		//TOOL::SimpleBaseKNNSearchMulti(data_source, query_time_series_vector, squential_scan_result_set, linear_scan_struct);
		//211205 vector of squential_scan_result_set
		vector<multiset<pair<double, int>>> squential_scan_result_set_vector(size_query_time_series);
		TOOL::SimpleBaseKNNSearchMulti(data_source, query_series_vector_vector, squential_scan_result_set_vector, linear_scan_struct);
		/*========================================================================================================================================*/

		/*=========================================191028 Initial if Y-Projection Method======================================================================*/
		//TOOL::Y_PROJECTION_ARGUMENT y_projection_argument(initial_N / 3, change_file);
		//DoublyLinkedList<APLA::AREA_COEFFICIENT> all_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
		//DoublyLinkedList<APLA::AREA_COEFFICIENT> cluster_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT>();//191030
		//vector<TOOL::Y_PROJECTION_ARGUMENT> multi_y_projection_argument_vector(data_source.point_number, TOOL::Y_PROJECTION_ARGUMENT(initial_N / 3));
		vector<TOOL::Y_PROJECTION_ARGUMENT> multi_y_projection_argument_vector(data_source.point_number, TOOL::Y_PROJECTION_ARGUMENT());
		vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>> multi_all_linked_list(data_source.point_number, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>());//191030
		vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>> multi_cluster_linked_list(data_source.point_number, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>());//191030
		/*=====================================================================================================================================================*/

		/*#######################################################                  N                    ################################################################*/
		for (int od_N = 0; od_N < N_coefficient_vector.size(); od_N++) {//N  5-80; 6 96 - 192;384 2;191119 12 24 48 96 192 258 384;;; 192 96 24, 
			N = N_coefficient_vector[od_N];
			//if (file_id == initial_file_od) { N_coefficient_vector.emplace_back(N); }

			/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   Normal Operation  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
			if (evaluation_type == 0) {//For Approximation methods compariosn. 1 MSPLA, 2 PLA, 3 APCA, 4 PAA, 5 Chebyshev, 6 ICDE07, 7 PAALM

				/*|||||||||||||||||||||||||||||||||||||||||||||       Type of Approximation Methods        ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
				for (representation_option = 1; representation_option <= method_total_number; representation_option++) {// 1 MSPLA, 2 PLA, 3APCA, 4 PAA, 5 Chebyshev, 6 ICDE07

					/*!!!!!!!!!!!!!!!!!!!!!!!!!!!         Initial N£º PAA & Cheby 12, APCA & PLA 6, APLA & ICDE07 4              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
					// Split Methods Evaluation
					for (int split_method_id = 0; split_method_id < split_methods_number; split_method_id++) {
						//representation_option = 5;
						//change_file = true;
						/*=======================================            MULTI Initial       ====================================================*/

						/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
						option_tree_struct.type_representation = representation_option_vector[representation_option];
						assert(option_tree_struct.type_representation > 1);
						MULTI multi(data_source, n0, N, INF, file_id_choose_vector[file_id], change_file, data_source.time_series_dimension, data_source.point_number, query_time_series_id, max_node0, initial_K, representation_option_vector[representation_option], file_address_pointer, chebyshev_write_pointer, option_tree_struct);
						/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

						/**********************         Query Time Series Vector 211220         ************************/
						multi.input_argument.id_query_series_vector.assign(Result_Coefficient_Struct.id_query_series_vector.begin(), Result_Coefficient_Struct.id_query_series_vector.end());
						multi.input_argument.query_series_vector_vector.assign(query_series_vector_vector.begin(), query_series_vector_vector.end());
						/************************************************************************************************/

						/**********************         Linear Scan          ************************/
						if (representation_option_vector[representation_option] == 10) {// Linear scan 
							TOOL::initial_evaluation_linear_scan(multi.input_argument, multi.output_argument, linear_scan_struct);
						}
						/*****************************************************************************/

						/*******************************        Initial input_argument        *********************************/

						/*.201028 Initial Counts of True Max deviation point / False Max deviation point.*/
						multi.input_argument.number_point_max_deviation_true = 0;
						multi.input_argument.number_point_max_deviation_false = 0;
						multi.input_argument.number_not_smaller_than_sum_deviation = 0;
						multi.input_argument.number_smaller_than_sum_deviation = 0;
						multi.input_argument.number_not_smaller_than_sum_deviation_pow = 0;
						multi.input_argument.number_smaller_than_sum_deviation_pow = 0;
						/*...............................................................................*/

						/*..........     200901 option, print each time series result       .............*/
						multi.input_argument.print_each_result = print_each_result;
						/*...............................................................................*/

						/*..................        200224 split id initial = 0       ...................*/
						APLA::initial_split_coefficients(multi.input_argument, multi.output_argument);
						multi.input_argument.option_split_method = split_method_id;
						/*...............................................................................*/

						/*****************************************************************************************************/
						/*===========================================================================================================================*/

						//multi.project_multi_to_single(multi_y_projection_argument_vector, multi_all_linked_list, multi_cluster_linked_list);
						//multi.buid_rtree_knn();
						/*=======================================      Approximation Process     ================================================*/
						multi.all_approximation_build_rtree(multi_y_projection_argument_vector, multi_all_linked_list, multi_cluster_linked_list);
						/*=======================================================================================================================*/

						/*-----------------------------      Evaluate upper bound       ---------------------------*/
						APLA::count_upper_bound_whole(multi.output_argument, evaluation_bound_whole);
						evaluation_bound_whole.difference_id_max_deviation_vs_height_diff += multi.output_argument.evaluation_bound.difference_id_max_deviation_vs_height_diff;
						/*-----------------------------------------------------------------------------------------*/

						const size_t id_method = representation_option - 1;
						const size_t id_each_method_each_file = data_source.file_operation_number * id_method + file_id;
						/*===============================================================         K          ==========================================================*/
						for (int od_K = 0; od_K < K_coefficient_vector.size(); od_K++) {//K 2 4 8 16 32 64 128 ;;; 16 128
							K = K_coefficient_vector[od_K];
							/*----------------------------------------------------------------------------------------------------------------------*/
							//if (file_id == initial_file_od && N == initial_N && representation_option == 1) { K_coefficient_vector.emplace_back(K); }
							/*----------------------------------------------------------------------------------------------------------------------*/
							multi.input_argument.K = K;

							/*::::::::::::::::::::::::::::::::::::::::::::            KNN process         ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
							//multi.all_knn(data_source, query_time_series_vector, query_time_series, multi_y_projection_argument_vector, squential_scan_result_set);
							multi.all_knn(data_source, multi.input_argument, query_series_vector_vector, multi_y_projection_argument_vector, squential_scan_result_set_vector);
							/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

							/*::::::::::::::::::::::::::::::::::::::::::::   Print Result, Matlab, bar chart    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
							cout << "#||||||||||      " << name_representation_vector_string_vector[id_method] << " Each K Result     ||||||||||||\n";
							cout << "file id : " << file_id_choose_vector[file_id] + 1 << ", K: " << K << ", N:" << multi.input_argument.point_dimension << ", representation id: " << representation_option_vector[representation_option] << ", n: " << n0 << ", dimension: " << data_source.time_series_dimension << ", point number: " << data_source.point_number << ", query time series size: " << size_query_time_series << ", max sub node: " << max_node0 << endl;
							print_input_output(multi.input_argument, multi.output_argument, linear_scan_struct);
							cout << "#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n";

							if(representation_option_vector[representation_option] != 10 && representation_option_vector[representation_option] != 11){
								assert(multi.input_argument.pruning_power > 0 && multi.input_argument.pruning_power != INF && multi.input_argument.prune_power_combine > 0 && multi.input_argument.prune_power_combine != INF && multi.input_argument.whole_run_time != INF && multi.input_argument.knn_CPU_time > 0);
							}
							else if (representation_option_vector[representation_option] == 10) {
								assert(multi.input_argument.knn_total_time == linear_scan_struct.time_linear_scan_CPU && multi.input_argument.knn_total_time_has_IO == linear_scan_struct.time_linear_scan_wall_clock && multi.input_argument.IO_cost == linear_scan_struct.IO_linear_scan);
							}

							/*...................................        200219 split id evaluation        ....................................*/
							// 1 min density,2 binary search,3  direct intersection£¬ 4 middle point and 5 best split point method.
							APLA::get_split_coefficients(multi.input_argument, multi.output_argument, multi.input_argument.option_split_method, local_total_split_id_sum_deviation, local_total_split_id_shift, local_total_split_id_time, global_total_approximation_sum_deviation, global_total_approximation_time, global_total_knn_prune_power);
							/*.................................................................................................................*/

							/******************************************             Each KNN result           *****************************************************/
							/*---------     200831 Barchart: Sum deviation & representation time has no business with K     ------------*/
							if (K == initial_K) {
								total_sum_deviation_vector.emplace_back(multi.output_argument.sum_deviation);//191206
								total_approximation_time_vector.emplace_back(multi.input_argument.representation_time);//200108
								total_build_tree_time_vector.emplace_back(multi.input_argument.build_rtree_time);//200108
								total_ingest_data_time_vector.emplace_back(multi.input_argument.time_ingest_data);//200108
								total_max_deviation_vector.emplace_back(multi.output_argument.max_deviation);//191206
								total_max_deviation_av_vector.emplace_back(multi.output_argument.max_deviation_av);//191206
								total_max_width_deviation_vector.emplace_back(multi.output_argument.max_deviation_multiple_width);//191206
								Total_Result_Struct.total_internal_node_size_vector.emplace_back(multi.input_argument.argument_index_struct.count_node_internal);//211215
								Total_Result_Struct.total_leaf_node_size_vector.emplace_back(multi.input_argument.argument_index_struct.count_node_leaf);//211215
								Total_Result_Struct.total_total_node_size_vector.emplace_back(multi.input_argument.argument_index_struct.count_node_total);//211215
								Total_Result_Struct.total_index_height_vector.emplace_back(multi.input_argument.argument_index_struct.height_index);//211215
							}
							else {
								total_sum_deviation_vector.emplace_back(0.0);//191206
								total_approximation_time_vector.emplace_back(0.0);//200108
								total_build_tree_time_vector.emplace_back(0.0);//200108
								total_ingest_data_time_vector.emplace_back(0.0);//211210
								total_max_deviation_vector.emplace_back(0.0);
								total_max_deviation_av_vector.emplace_back(0.0);
								total_max_width_deviation_vector.emplace_back(0.0);
								Total_Result_Struct.total_internal_node_size_vector.emplace_back(0.0);//211215
								Total_Result_Struct.total_leaf_node_size_vector.emplace_back(0.0);//211215
								Total_Result_Struct.total_total_node_size_vector.emplace_back(0.0);//211215
								Total_Result_Struct.total_index_height_vector.emplace_back(0.0);//211215
							}
							/*-----------------------------------------------------------------------------------------------------------*/

							total_prune_power_combine_vector.emplace_back(multi.input_argument.prune_power_combine);//sum pruning power combine 201221
							total_prune_power_vector.emplace_back(multi.input_argument.pruning_power);//sum pruning power
							//total_sum_deviation_vector.emplace_back(multi.output_argument.sum_deviation);//191206
							total_accuracy_vector.emplace_back(multi.input_argument.result_accuracy);
							total_run_time_vector.emplace_back(multi.input_argument.whole_run_time);
							total_run_time_has_IO_vector.emplace_back(multi.input_argument.whole_run_time_has_IO);//210606
							//total_approximation_time_vector.emplace_back(multi.input_argument.representation_time);//200108
							Total_Result_Struct.total_knn_cpu_time_vector.emplace_back(multi.input_argument.knn_CPU_time);//211215
							total_knn_time_vector.emplace_back(multi.input_argument.knn_total_time);//200108
							total_knn_time_has_IO_vector.emplace_back(multi.input_argument.knn_total_time_has_IO);//210606
							total_IO_vector.emplace_back(multi.input_argument.IO_cost);//211203
							/***************************************************************************************************************************************/

							/***************************************         Sum for each Approxiamtion method       ***********************************************/

							//for every method
							method_sum_prune_power_combine_vector[id_method] += multi.input_argument.prune_power_combine; //each approximation method pruning power combine 201221
							method_sum_prune_power_vector[id_method] += multi.input_argument.pruning_power;
							//method_sum_sum_deviation_vector[representation_option - 1] += multi.output_argument.sum_deviation;//191206
							method_sum_accuracy_vector[id_method] += multi.input_argument.result_accuracy;//191204 accuracy
							method_run_time_vector[id_method] += multi.input_argument.whole_run_time;
							method_run_time_has_IO_vector[id_method] += multi.input_argument.whole_run_time_has_IO;//210606
							//method_approximation_time_vector[representation_option - 1] += multi.input_argument.representation_time;
							Method_Result_Struct.method_knn_cpu_time_vector[id_method] += multi.input_argument.knn_CPU_time;//211215
							method_knn_time_vector[id_method] += multi.input_argument.knn_total_time;
							method_knn_time_has_IO_vector[id_method] += multi.input_argument.knn_total_time_has_IO;//210606
							method_IO_vector[id_method] += multi.input_argument.IO_cost;//211203
							/***************************************************************************************************************************************/

							/***************************************     201028 One File for each Approximation method    ******************************************/
							if (!clear_file) {
								one_file_prune_power_combine_vector[id_method] += multi.input_argument.prune_power_combine;//each file pruning power combine 201221
								one_file_prune_power_vector[id_method] += multi.input_argument.pruning_power;
								one_file_accuracy_vector[id_method] += multi.input_argument.result_accuracy;//191204 accuracy
								one_file_run_time_vector[id_method] += multi.input_argument.whole_run_time;
								one_file_knn_time_vector[id_method] += multi.input_argument.knn_total_time;
							}
							/***************************************************************************************************************************************/

							/***********************************        Sum for each File of each Approximation method        **************************************/
							//for every method, each file, Split as method
							method_file_prune_power_combine_vector[id_each_method_each_file] += multi.input_argument.prune_power_combine;//each file pruning power combine 201221;
							method_file_prune_power_vector[id_each_method_each_file] += multi.input_argument.pruning_power;
							//method_file_sum_deviation_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.output_argument.sum_deviation;//191206
							method_file_accuracy_vector[id_each_method_each_file] += multi.input_argument.result_accuracy;//191204 accuracy
							method_file_run_time_vector[id_each_method_each_file] += multi.input_argument.whole_run_time;
							method_file_run_time_has_IO_vector[id_each_method_each_file] += multi.input_argument.whole_run_time_has_IO;//210606
							//method_file_approximation_time_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.representation_time;
							Method_File_Result_Struct.method_file_knn_cpu_time_vector[id_each_method_each_file] += multi.input_argument.knn_CPU_time;//211215
							method_file_knn_time_vector[id_each_method_each_file] += multi.input_argument.knn_total_time;
							method_file_knn_time_has_IO_vector[id_each_method_each_file] += multi.input_argument.knn_total_time_has_IO;//210606
							method_file_IO_vector[id_each_method_each_file] += multi.input_argument.IO_cost;//211203
							/***************************************************************************************************************************************/

							/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
							//ttool.getMultiFoldToSingleByID(multi_file_pointer,3,10,0, test_time_series);
							//system("pause");
						}
						/*===============================================================================================================================================================================*/

						/*:::::::::::::::::::::::::      200831 Sum Deviation &  representation time has no business with K      ::::::::::::::::::::::::::::::::::::*/
						/*--------------------------------------     Sum for each Approxiamtion method     ------------------------------*/
						method_sum_sum_deviation_vector[id_method] += multi.output_argument.sum_deviation;//191206
						method_approximation_time_vector[id_method] += multi.input_argument.representation_time;
						method_build_tree_time_vector[id_method] += multi.input_argument.build_rtree_time;
						method_ingest_data_time_vector[id_method] += multi.input_argument.time_ingest_data;//211211
						method_sum_max_deviation_vector[id_method] += multi.output_argument.max_deviation;
						method_sum_max_width_deviation_vector[id_method] += multi.output_argument.max_deviation_multiple_width;
						Method_Result_Struct.method_internal_node_size_vector[id_method] += multi.input_argument.argument_index_struct.count_node_internal;//61
						Method_Result_Struct.method_leaf_node_size_vector[id_method] += multi.input_argument.argument_index_struct.count_node_leaf;//64
						Method_Result_Struct.method_total_node_size_vector[id_method] += multi.input_argument.argument_index_struct.count_node_total;//67
						Method_Result_Struct.method_index_height_vector[id_method] += multi.input_argument.argument_index_struct.height_index;//70
						/*---------------------------------------------------------------------------------------------------------------*/

						/*---------------------------------------        201028 For One File              -------------------------------*/
						if (!clear_file) {
							one_file_sum_deviation_vector[id_method] += multi.output_argument.sum_deviation;//191206
							one_file_approximation_time_vector[id_method] += multi.input_argument.representation_time;
						}
						/*---------------------------------------------------------------------------------------------------------------*/

						/*-----------------------------         Sum Deviation for each File of each Approximation method    -----------------------*/
						method_file_sum_deviation_vector[id_each_method_each_file] += multi.output_argument.sum_deviation;//191206
						method_file_approximation_time_vector[id_each_method_each_file] += multi.input_argument.representation_time;
						method_file_build_tree_time_vector[id_each_method_each_file] += multi.input_argument.build_rtree_time;
						method_file_ingest_data_time_vector[id_each_method_each_file] += multi.input_argument.time_ingest_data;//211211
						method_file_max_deviation_vector[id_each_method_each_file] += multi.output_argument.max_deviation;//191206
						method_file_max_width_deviation_vector[id_each_method_each_file] += multi.output_argument.max_deviation_multiple_width;//191206
						Method_File_Result_Struct.method_file_internal_node_size_vector[id_each_method_each_file] += multi.input_argument.argument_index_struct.count_node_internal;
						Method_File_Result_Struct.method_file_leaf_node_size_vector[id_each_method_each_file] += multi.input_argument.argument_index_struct.count_node_leaf;
						Method_File_Result_Struct.method_file_total_node_size_vector[id_each_method_each_file] += multi.input_argument.argument_index_struct.count_node_total;
						Method_File_Result_Struct.method_file_index_height_vector[id_each_method_each_file] += multi.input_argument.argument_index_struct.height_index;
						/*-------------------------------------------------------------------------------------------------------------------------*/
						/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

						/*==========   201028 Probability of max deviation points are minmax, left&right endpoints   =============*/
						number_point_max_deviation_true += multi.input_argument.number_point_max_deviation_true;
						number_point_max_deviation_false += multi.input_argument.number_point_max_deviation_false;
						number_not_smaller_than_sum_deviation += multi.input_argument.number_not_smaller_than_sum_deviation;
						number_smaller_than_sum_deviation += multi.input_argument.number_smaller_than_sum_deviation;
						number_not_smaller_than_sum_deviation_pow += multi.input_argument.number_not_smaller_than_sum_deviation_pow;
						number_smaller_than_sum_deviation_pow += multi.input_argument.number_smaller_than_sum_deviation_pow;
						/*=========================================================================================================*/
					}//K
					/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
				}//representation
				/*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
				change_file = false;

			}/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
			else if (evaluation_type == 1) {//For different SAPLA initial N comparison. initial_N = 0,1,2,3,4,5,6,7,8,9
				assert(0);
				///*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
				//int initial_N_number = INF;
				//double initial_number_segment = INF;
				//assert(initial_number_coefficients_vector.size() == final_initial_N_number);
				//for (initial_N_number = 0; initial_N_number < final_initial_N_number; initial_N_number++) {
				//	/*----------------------------------------Initial Segment Number Caluclation-------------------------------------------------------------*/
				//	initial_number_segment = initial_number_coefficients_vector[initial_N_number] + N / 3;
				//	if (file_id == initial_file_od) {
				//		approximation_initial_N_vector.emplace_back(initial_number_segment);
				//	}
				//	/*----------------------------------------------------------------------------------------------------------------------------------------*/
				//	/*------------------------------Initial N£º PAA & Cheby 12, APCA & PLA 6, APLA & ICDE07 4-------------------------------------------------*/
				//	//representation_option = 5;
				//	//change_file = true;
				//	representation_option = 1;//only evaluate SAPLA
				//	MULTI multi(data_source, n0, N, initial_number_segment, file_id, change_file, data_source.time_series_dimension, data_source.point_number, query_time_series_id, max_node0, initial_K, representation_option, file_address_pointer, chebyshev_write_pointer);
				//	//representation_option = N + initial_number;
				//	//for (int split_method_id = 0; split_method_id < split_methods_number; split_method_id++) {
				//	/*...................200224 split id initial = 0.............................*/
				//	APLA::initial_split_coefficients(multi.input_argument, multi.output_argument);
				//	multi.input_argument.option_split_method = 0;//split_method_id
				//	/*...............................................................................*/
				//	/*-----------------------------------------------------------------------------------------------------------------------------------------*/
				//	//multi.project_multi_to_single(multi_y_projection_argument_vector, multi_all_linked_list, multi_cluster_linked_list);
				//	//multi.buid_rtree_knn();
				//	multi.all_approximation_build_rtree(multi_y_projection_argument_vector, multi_all_linked_list, multi_cluster_linked_list);

				//	/*####################################################################          K         ###############################################################*/
				//	for (K = initial_K; K <= final_K; K *= 2) {//K 2 4  8 16 32 64 128;;; 16 128
				//		if (file_id == initial_file_od && N == initial_N && initial_N_number == 0) { K_coefficient_vector.emplace_back(K); }
				//		multi.input_argument.K = K;
				//		multi.all_knn(data_source, query_time_series_vector, query_time_series, multi_y_projection_argument_vector, squential_scan_result_set);

				//		/*====================================================      Write Result, Matlab, bar chart   ==============================================*/
				//		cout << "---------------------------       Result      ------------------------------\n";
				//		cout << "file id : " << file_id + 1 << ", K: " << K << ", N:" << multi.input_argument.point_dimension << ", representation id: " << representation_option << ", n: " << n0 << ", dimension: " << data_source.time_series_dimension << ", point number: " << data_source.point_number << ", query time series id: " << query_time_series_id << ", max sub node: " << max_node0 << endl;
				//		cout << "prune power: " << multi.input_argument.pruning_power << endl;
				//		cout << "sum deviation: " << multi.output_argument.sum_deviation << endl;
				//		cout << "prune accuracy: " << multi.input_argument.result_accuracy << endl;
				//		cout << "representation_time: " << multi.input_argument.representation_time << endl;
				//		cout << "build_rtree_time: " << multi.input_argument.build_rtree_time << endl;
				//		cout << "knn_total_time: " << multi.input_argument.knn_total_time << endl;
				//		cout << "accumulation time : " << multi.input_argument.representation_time + multi.input_argument.build_rtree_time + multi.input_argument.knn_total_time << endl;
				//		cout << "whole running time: " << multi.input_argument.whole_run_time << endl;
				//		cout << "----------------------------------------------------------------------------\n";

				//		assert(multi.input_argument.pruning_power <= 20 && multi.input_argument.pruning_power != INF && multi.input_argument.whole_run_time != INF);

				//		/*........................200219 split id evaluation........................*/
				//		// 1 min density,2 binary search,3  direct intersection and 4 best split point method.
				//		APLA::get_split_coefficients(multi.input_argument, multi.output_argument, multi.input_argument.option_split_method, local_total_split_id_sum_deviation, local_total_split_id_shift, local_total_split_id_time, global_total_approximation_sum_deviation, global_total_approximation_time, global_total_knn_prune_power);
				//		/*..........................................................................*/

				//		/*-------------------------------------   200326 Evaluate initial_N  ----------------------------------------------------------------------*/
				//		//get total initial_N coefficients
				//		APLA::get_initial_N_coefficients(multi.input_argument, multi.output_argument, initial_N_number, total_initial_N_prune_power_vector, total_initial_N_sum_deviation_vector, total_initial_N_run_time_vector, total_initial_N_approximation_time_vector, total_initial_N_knn_time_vector);
				//		//get initial_N coefficients by N
				//		const int N_id = log2(N / initial_N);
				//		APLA::get_initial_N_sort_N_coefficients(multi.input_argument, multi.output_argument, int(initial_number_coefficients_vector.size()), N_id, initial_N_number, initial_N_by_N_prune_power_vector, initial_N_by_N_sum_deviation_vector, initial_N_by_N_run_time_vector, initial_N_by_N_approximation_time_vector, initial_N_by_N_knn_time_vector);
				//		/*-----------------------------------------------------------------------------------------------------------------------------------------*/

				//		/*----------------------------------------------------------------------Normal Coeffients-----------------------------------------------------------------*/
				//		//sum pruning power
				//		total_prune_power_vector.emplace_back(multi.input_argument.pruning_power);
				//		total_sum_deviation_vector.emplace_back(multi.output_argument.sum_deviation);//191206
				//		total_accuracy_vector.emplace_back(multi.input_argument.result_accuracy);
				//		total_run_time_vector.emplace_back(multi.input_argument.whole_run_time);
				//		total_approximation_time_vector.emplace_back(multi.input_argument.representation_time);//200108
				//		total_knn_time_vector.emplace_back(multi.input_argument.knn_total_time);//200108
				//		//for every method 
				//		method_sum_prune_power_vector[representation_option - 1] += multi.input_argument.pruning_power;
				//		method_sum_sum_deviation_vector[representation_option - 1] += multi.output_argument.sum_deviation;//191206
				//		method_sum_accuracy_vector[representation_option - 1] += multi.input_argument.result_accuracy;//191204 accuracy
				//		method_run_time_vector[representation_option - 1] += multi.input_argument.whole_run_time;
				//		method_approximation_time_vector[representation_option - 1] += multi.input_argument.representation_time;
				//		method_knn_time_vector[representation_option - 1] += multi.input_argument.knn_total_time;
				//		//for every method, each file
				//		method_file_prune_power_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.pruning_power;
				//		method_file_sum_deviation_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.output_argument.sum_deviation;//191206
				//		method_file_accuracy_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.result_accuracy;//191204 accuracy
				//		method_file_run_time_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.whole_run_time;
				//		method_file_approximation_time_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.representation_time;
				//		method_file_knn_time_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.knn_total_time;

				//		/*--------------------------------------------------------------------------------------------------------------------------------------------------*/

				//		/*==========================================================================================================================================================*/
				//	}
				//	/*#############################################################################################################################################################################################*/
				//}
				//change_file = false;
				////}
				////assert(approximation_initial_N_vector.size() == final_initial_N_number);
				///*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
			}
			else {
				assert(0);
			}
		}//N

		/*#############################################################################################################################################################################################################*/

		/*----------------------read file address vector-----------------------------------*/
		data_source.read_file_address_vector.clear();
		data_source.read_file_address_vector.shrink_to_fit();
		/*---------------------------------------------------------------------------------*/

		/*----------------------Clear Y projection  Memory---------------------------------*/
		for (auto&& au : multi_y_projection_argument_vector) {
			au.whole_difference_map.clear();
		}
		multi_y_projection_argument_vector.clear();
		multi_y_projection_argument_vector.shrink_to_fit();

		for (auto&& au : multi_all_linked_list) {
			au.clear();
		}
		multi_all_linked_list.clear();
		multi_all_linked_list.shrink_to_fit();

		for (auto&& au : multi_cluster_linked_list) {
			au.clear();
		}
		multi_cluster_linked_list.clear();
		multi_cluster_linked_list.shrink_to_fit();
		/*--------------------------------------------------------------------------------*/

		/*delete[] query_time_series;
		query_time_series = nullptr;*/

		/*||||||||||||||||||||||||||||||||||||||||    Write One File    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
		if (!clear_file) {
			//TOOL::writeResultNoCover(str_33, one_file_prune_power_combine_vector);//201221
			//TOOL::writeResultNoCover(str_27, one_file_prune_power_vector);
			//TOOL::writeResultNoCover(str_28, one_file_sum_deviation_vector);
			//TOOL::writeResultNoCover(str_29, one_file_accuracy_vector);
			//TOOL::writeResultNoCover(str_30, one_file_run_time_vector);
			//TOOL::writeResultNoCover(str_31, one_file_approximation_time_vector);
			//TOOL::writeResultNoCover(str_32, one_file_knn_time_vector);
			//one_file_prune_power_combine_vector.resize(method_total_number, 0.0);//201221
			//one_file_prune_power_vector.resize(method_total_number, 0.0);
			//one_file_sum_deviation_vector.resize(method_total_number, 0.0);
			//one_file_accuracy_vector.resize(method_total_number, 0.0);
			//one_file_run_time_vector.resize(method_total_number, 0.0);
			//one_file_approximation_time_vector.resize(method_total_number, 0.0);
			//one_file_knn_time_vector.resize(method_total_number, 0.0);
		}
		/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||!!!!!!!!!!!!!||||||||||||||||||||||||||||||||||||||||*/

	}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&        Print Result        &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      Final Result     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";

	/*------------------------------------------------------- 201028 Print Probability max deviation point ------------------------------------------------------------------*/
	cout << "!!!!Probability max deviation point. Total number: " << number_point_max_deviation_true + number_point_max_deviation_false << ", True number: " << number_point_max_deviation_true << ", False number: " << number_point_max_deviation_false << ", Probability: " << number_point_max_deviation_true / (number_point_max_deviation_true + number_point_max_deviation_false) << endl;
	cout << "!!!!Probability bigger than sum deviation.Total number : " << number_not_smaller_than_sum_deviation + number_smaller_than_sum_deviation << ", True number: " << number_not_smaller_than_sum_deviation << ", False number: " << number_smaller_than_sum_deviation << ", Probability: " << number_not_smaller_than_sum_deviation / (number_not_smaller_than_sum_deviation + number_smaller_than_sum_deviation) << endl;
	cout << "!!!!Probability bigger than sum deviation pow. Total number: " << number_not_smaller_than_sum_deviation_pow + number_smaller_than_sum_deviation_pow << ", True number: " << number_not_smaller_than_sum_deviation_pow << ", False number: " << number_smaller_than_sum_deviation_pow << ", Probability: " << number_not_smaller_than_sum_deviation_pow / (number_not_smaller_than_sum_deviation_pow + number_smaller_than_sum_deviation_pow) << endl;
	/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

	/*-------------------------------------------------------        Transfer time us to s        ---------------------------------------------------------------------*/
	TOOL::transfer_us_s(total_run_time_vector);
	TOOL::transfer_us_s(total_run_time_has_IO_vector);//210606
	TOOL::transfer_us_s(total_approximation_time_vector);
	TOOL::transfer_us_s(total_build_tree_time_vector);
	TOOL::transfer_us_s(total_ingest_data_time_vector);
	TOOL::transfer_us_s(Total_Result_Struct.total_knn_cpu_time_vector);//211215
	TOOL::transfer_us_s(total_knn_time_vector);
	TOOL::transfer_us_s(total_knn_time_has_IO_vector);//210606

	TOOL::transfer_us_s(method_run_time_vector);
	TOOL::transfer_us_s(method_run_time_has_IO_vector);//210606
	TOOL::transfer_us_s(method_approximation_time_vector);
	TOOL::transfer_us_s(method_build_tree_time_vector);
	TOOL::transfer_us_s(method_ingest_data_time_vector);//211211
	TOOL::transfer_us_s(Method_Result_Struct.method_knn_cpu_time_vector);//211215
	TOOL::transfer_us_s(method_knn_time_vector);
	TOOL::transfer_us_s(method_knn_time_has_IO_vector);//210606

	TOOL::transfer_us_s(one_file_run_time_vector);
	TOOL::transfer_us_s(one_file_approximation_time_vector);
	TOOL::transfer_us_s(one_file_knn_time_vector);
	TOOL::transfer_us_s(one_file_knn_time_has_IO_vector);//210606

	TOOL::transfer_us_s(method_file_run_time_vector);
	TOOL::transfer_us_s(method_file_run_time_has_IO_vector);//210606
	TOOL::transfer_us_s(method_file_approximation_time_vector);
	TOOL::transfer_us_s(method_file_build_tree_time_vector);
	TOOL::transfer_us_s(method_file_ingest_data_time_vector);//211211
	TOOL::transfer_us_s(Method_File_Result_Struct.method_file_knn_cpu_time_vector);//211215
	TOOL::transfer_us_s(method_file_knn_time_vector);
	TOOL::transfer_us_s(method_file_knn_time_has_IO_vector);//210606
	/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

	/*-------------------------------------------      Print Normal Result      -----------------------------------------*/
	////1 MSPLA, 2 PLA, 3APCA, 4 PAA, 5 Chebyshev, 6 ICDE07
	//cout << "Method Sum Prune Power Combine: \n";
	//TOOL::print_string_vector(name_representation_vector_string_vector);
	//for (auto&& au : method_sum_prune_power_combine_vector) {
	//	cout << au << ", ";
	//}
	//cout << endl;
	
	cout << "Method Max Deviation: \n";
	//cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07  PAALM\n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_sum_max_deviation_vector) {
		cout << au << ", ";
	}
	cout << endl;
	//1 MSPLA, 2 PLA, 3APCA, 4 PAA, 5 Chebyshev, 6 ICDE07
	cout << "Method Sum Prune Power: \n";
	//cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07   PAALM\n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_sum_prune_power_vector) {
		cout << au << ", ";
	}
	cout << endl;
	cout << "Method Sum Accuracy: \n";
	//cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07  PAALM\n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_sum_accuracy_vector) {
		cout << au << ", ";
	}
	cout << endl;
	cout << "Dimension Deduction Sum Time(s): \n";
	//cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07  PAALM\n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_approximation_time_vector) {
		cout << au << ", ";
	}
	cout << endl;

	cout << "Ingest Data Time(s): \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	TOOL::print_vector(method_ingest_data_time_vector);


	cout << "KNN Sum Time(s): \n";
	//cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07  PAALM\n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_knn_time_vector) {
		cout << au << ", ";
	}
	cout << endl;
	

	cout << "Tree Height: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	TOOL::print_vector(Method_Result_Struct.method_index_height_vector);

	cout << "Internal Node: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	TOOL::print_vector(Method_Result_Struct.method_internal_node_size_vector);

	cout << "Leaf Node: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	TOOL::print_vector(Method_Result_Struct.method_leaf_node_size_vector);

	cout << "Total Node: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	TOOL::print_vector(Method_Result_Struct.method_total_node_size_vector);

	cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
	cout << data_source.bigger_account << endl;
	/*-----------------------------------------------------------------------------------------------------------------------*/
	/*------------------------------------------------200219 Print split id evaluation---------------------------------------*/
	//TOOL::print_split_coefficients(local_total_split_id_sum_deviation, local_total_split_id_shift, local_total_split_id_time, global_total_approximation_sum_deviation, global_total_approximation_time, global_total_knn_prune_power);
	/*-----------------------------------------------------------------------------------------------------------------------*/
	/*-----------------------------------------------200326 Print initial_N evaluation---------------------------------------*/
	// total result
	//TOOL::print_initial_N_coefficients(approximation_initial_N_vector, total_initial_N_prune_power_vector, total_initial_N_sum_deviation_vector, total_initial_N_run_time_vector, total_initial_N_approximation_time_vector, total_initial_N_knn_time_vector);
	// sort by N
	//TOOL::print_initial_N_sort_N_coefficients(initial_number_coefficients_vector.size(), approximation_initial_N_vector, initial_N_by_N_prune_power_vector, initial_N_by_N_sum_deviation_vector, initial_N_by_N_run_time_vector, initial_N_by_N_approximation_time_vector, initial_N_by_N_knn_time_vector);
	/*-----------------------------------------------------------------------------------------------------------------------*/
	/*--------------------------------------------- 210122 Print Upper Bound evaluation -------------------------------------*/
	//APLA::print_count_of_upper_bound(evaluation_bound_whole);
	/*-----------------------------------------------------------------------------------------------------------------------*/
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&          Write Result Barchart from Matlab         &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	//K_total_number = K_coefficient_vector.size();
	int total_size = int(INF);

	switch (evaluation_type) {
	case 0:
		total_size = K_coefficient_vector.size() * N_coefficient_vector.size() * data_source.file_operation_number * method_total_number;
		break;
	case 1:
		total_size = K_coefficient_vector.size() * N_coefficient_vector.size() * data_source.file_operation_number * method_total_number * final_initial_N_number;
		break;
	default:
		assert(0);
		break;
	}

	assert(total_prune_power_vector.size() == total_size && total_prune_power_combine_vector.size() == total_run_time_vector.size() && total_prune_power_vector.size() == total_run_time_vector.size() && total_prune_power_vector.size() == total_accuracy_vector.size() && file_id_choose_vector.size() == data_source.file_operation_number);
	assert(!K_coefficient_vector.empty(), !N_coefficient_vector.empty(), data_source.file_total_number != INF, method_total_number != INF, data_source.point_number != INF, query_time_series_id != INF, max_node0 != INF);
	
	vector<int> type_index_vector(method_total_number, option_tree_struct.type_tree);// 0 is Rtree, 1 is Dist Tree
	// 1bK_total_number , 2bN_total_number , 3bfile_total_number, 4bmethod_total_number, 5bpoint number, 6bquery size, 7bmax id, 8 single time series length,
	// 9 data type, 10 evaluation type, 11 split methods number, 12 option homogenous data type, 13 index_type
	vector<int> barchart_coefficient_vector = { int(K_coefficient_vector.size()), int(N_coefficient_vector.size()), data_source.file_operation_number, method_total_number, data_source.point_number, int(Result_Coefficient_Struct.id_query_series_vector.size()), max_node0, n0, data_source.data_type, evaluation_type, split_methods_number, data_source.option_has_burst_data, int(type_index_vector.size()) };

	/*-------------------------------------------Transfer Matlab result for 3D BarChart--------------------------------------------------------*/
	vector<long double> transfer_total_prune_power_combine_vector;//201221
	vector<long double> transfer_total_prune_power_vector;//200103
	vector<long double> transfer_total_sum_deviation_vector;//200103
	vector<long double> transfer_total_max_deviation_vector;//210327
	vector<long double> transfer_total_max_deviation_av_vector;//210327
	vector<long double> transfer_total_max_width_deviation_vector;//210327
	vector<long double> transfer_total_accuracy_vector;//200103
	vector<long double> transfer_total_whole_run_time_vector;//200103
	vector<long double> transfer_total_whole_run_time_has_IO_vector;//210606
	vector<long double> transfer_total_approximation_run_time_vector;//200108
	vector<long double> transfer_total_build_tree_time_vector;//210906
	vector<long double> transfer_total_ingest_data_time_vector;//211210
	vector<long double> transfer_total_knn_cpu_time_vector;//200108
	vector<long double> transfer_total_knn_run_time_vector;//200108
	vector<long double> transfer_total_knn_run_time_has_IO_vector;//210606
	vector<long double> transfer_total_IO_vector;//211203
	vector<long double> transfer_total_internal_node_size_vector;//63
	vector<long double> transfer_total_leaf_node_size_vector;//66
	vector<long double> transfer_total_total_node_size_vector;//69
	vector<long double> transfer_total_index_height_vector;//72

	//Sort as methods < N < K < File
	//sort as file-> K -> N -> method
	for (int id_file_transfer = 0; id_file_transfer < data_source.file_operation_number; id_file_transfer++) {
		for (int K_id = 0; K_id < K_coefficient_vector.size(); K_id++) {
			for (int array_id = 0; array_id < N_coefficient_vector.size() * method_total_number; array_id++) {
				const size_t copy_id = array_id * K_coefficient_vector.size() + K_id + id_file_transfer * K_coefficient_vector.size() * N_coefficient_vector.size() * method_total_number;
				transfer_total_prune_power_combine_vector.emplace_back(total_prune_power_combine_vector[copy_id]);
				transfer_total_prune_power_vector.emplace_back(total_prune_power_vector[copy_id]);
				transfer_total_sum_deviation_vector.emplace_back(total_sum_deviation_vector[copy_id]);
				transfer_total_max_deviation_vector.emplace_back(total_max_deviation_vector[copy_id]);
				transfer_total_max_deviation_av_vector.emplace_back(total_max_deviation_av_vector[copy_id]);
				transfer_total_max_width_deviation_vector.emplace_back(total_max_width_deviation_vector[copy_id]);
				transfer_total_accuracy_vector.emplace_back(total_accuracy_vector[copy_id]);
				transfer_total_whole_run_time_vector.emplace_back(total_run_time_vector[copy_id]);
				transfer_total_whole_run_time_has_IO_vector.emplace_back(total_run_time_has_IO_vector[copy_id]);
				transfer_total_approximation_run_time_vector.emplace_back(total_approximation_time_vector[copy_id]);
				transfer_total_build_tree_time_vector.emplace_back(total_build_tree_time_vector[copy_id]);
				transfer_total_ingest_data_time_vector.emplace_back(total_ingest_data_time_vector[copy_id]);
				transfer_total_knn_cpu_time_vector.emplace_back(Total_Result_Struct.total_knn_cpu_time_vector[copy_id]);
				transfer_total_knn_run_time_vector.emplace_back(total_knn_time_vector[copy_id]);
				transfer_total_knn_run_time_has_IO_vector.emplace_back(total_knn_time_has_IO_vector[copy_id]);
				transfer_total_IO_vector.emplace_back(total_IO_vector[copy_id]);
				transfer_total_internal_node_size_vector.emplace_back(Total_Result_Struct.total_internal_node_size_vector[copy_id]);//63
				transfer_total_leaf_node_size_vector.emplace_back(Total_Result_Struct.total_leaf_node_size_vector[copy_id]);//66
				transfer_total_total_node_size_vector.emplace_back(Total_Result_Struct.total_total_node_size_vector[copy_id]);//69
				transfer_total_index_height_vector.emplace_back(Total_Result_Struct.total_index_height_vector[copy_id]);//72
			}
		}
	}
	/*--------------------------------------------------------------------------------------------------------------------------------------------*/

	/*-------------------------------------------------             String Name           --------------------------------------------------------*/
	ADDRESS_STRING address_string_struct(str_suffix);//211210
	const string& const str_1 = "./191120KNNMatlab/BarChartCoefficient" + str_suffix; //191122 chart coefficient
	const string& const str_25 = "./191120KNNMatlab/NameRepresentationMethod" + str_suffix;//200831 the name of representation method, string vector.
	const string& const str_26 = "./191120KNNMatlab/IDFileChoosed" + str_suffix;//201020 Choose File id: file_id_choose_vector
	const string& const str_4 = "./191120KNNMatlab/KCoefficient" + str_suffix;//191122 K records. eg. K = 2, 4 , 8
	const string& const str_11 = "./191120KNNMatlab/NCoefficient" + str_suffix;//191128 N records. eg. N = 18, 36, 72, 144, 288
	const string& const str_5 = "./191120KNNMatlab/ReadMe" + str_suffix;//191122
	const string& const str_10 = "./191120KNNMatlab/StringSuffix";//191126 prefix of file.

	const string& const str_2 = "./191120KNNMatlab/TotalPrunePower" + str_suffix;//191122 prune power
	const string& const str_34 = "./191120KNNMatlab/TotalPrunePowerCombine" + str_suffix;//201221 prune power combine. prune power / accuracy
	const string& const str_39 = "./191120KNNMatlab/TotalMaxDeviation" + str_suffix;//210327 max deviation
	const string& const str_40 = "./191120KNNMatlab/TotalMaxWidthDeviation" + str_suffix;//191204 sum deviation
	const string& const str_19 = "./191120KNNMatlab/TotalSumDeviation" + str_suffix;//191204 sum deviation
	const string& const str_18 = "./191120KNNMatlab/TotalAccuracy" + str_suffix;//191205 Accuracy
	const string& const str_3 = "./191120KNNMatlab/TotalKNNRunTime" + str_suffix;//191122 whole run time
	const string& const str_48 = "./191120KNNMatlab/TotalKNNRunTimeHasIO" + str_suffix;//210606 whole run time has IO
	const string& const str_22 = "./191120KNNMatlab/TotalApproximationRunTime" + str_suffix;//200108 approximation run time
	const string& const str_49 = "./191120KNNMatlab/TotalBuildTreeTime" + str_suffix;//210906 build tree time
	const string& const str_55 = "./191120KNNMatlab/TotalIngestDataTime" + str_suffix;//211210 ingest data time
	const string& const str_23 = "./191120KNNMatlab/TotalIndexRunTime" + str_suffix;//200108 index(KNN time) run time
	const string& const str_43 = "./191120KNNMatlab/TotalIndexRunTimeHasIO" + str_suffix;//210606 whole knn run time has IO
	const string& const str_52 = "./191120KNNMatlab/TotalIOCost" + str_suffix;//211203 IO cost 
	
	const string& const str_6 = "./191120KNNMatlab/MethodPrunePower" + str_suffix;//191123 size is total method number. prune power of every method for matlab barchart
	const string& const str_35 = "./191120KNNMatlab/MethodPrunePowerCombine" + str_suffix;//201221  prune power combine. prune power / accuracy
	const string& const str_41 = "./191120KNNMatlab/MethodMaxDeviation" + str_suffix;//210327
	const string& const str_42 = "./191120KNNMatlab/MethodMaxWidthDeviation" + str_suffix;//210327
	const string& const str_20 = "./191120KNNMatlab/MethodSumDeviation" + str_suffix;//191206
	const string& const str_16 = "./191120KNNMatlab/MethodAccuracy" + str_suffix;//191204 Accuracy of prune power
	const string& const str_7 = "./191120KNNMatlab/MethodKNNRunTime" + str_suffix;//191123 size is total method number. total run time of every method for matlab barchart
	const string& const str_44 = "./191120KNNMatlab/MethodKNNRunTimeHasIO" + str_suffix;//210606 size is total method number. total run time of every method for matlab barchart
	const string& const str_12 = "./191120KNNMatlab/MethodApproximationTime" + str_suffix;//191123 size is total method number. total run time of every method for matlab barchart
	const string& const str_50 = "./191120KNNMatlab/MethodBuildTreeTime" + str_suffix;//210906
	const string& const str_13 = "./191120KNNMatlab/MethodOnlyKNNTime" + str_suffix;//191123 size is total method number. total run time of every method for matlab barchart
	const string& const str_45 = "./191120KNNMatlab/MethodOnlyKNNTimeHasIO" + str_suffix;//210606 size is total method number. total run time of every method for matlab barchart
	const string& const str_53 = "./191120KNNMatlab/MethodIOCost" + str_suffix;//211203 IO cost 

	const string& const str_8 = "./191120KNNMatlab/FileMethodPrunePower" + str_suffix;//191125 size is total file number. prune power of every file for matlab barchart
	const string& const str_36 = "./191120KNNMatlab/FileMethodPrunePowerCombine" + str_suffix;//201221  prune power combine. prune power / accuracy
	const string& const str_37 = "./191120KNNMatlab/FileMethodMaxDeviation" + str_suffix;//210327 Max deviaiton 
	const string& const str_38 = "./191120KNNMatlab/FileMethodMaxWidthDeviation" + str_suffix;////210327 Max deviaiton  * width
	const string& const str_21 = "./191120KNNMatlab/FileMethodSumDeviation" + str_suffix;//191206
	const string& const str_17 = "./191120KNNMatlab/FileMethodAccuracy" + str_suffix;//191204 File Accuracy of prune power
	const string& const str_9 = "./191120KNNMatlab/FileMethodKNNRunTime" + str_suffix;//191125 size is total file number. total run time of every file for matlab barchart
	const string& const str_46 = "./191120KNNMatlab/FileMethodKNNRunTimeHasIO" + str_suffix;//210606 size is total file number. total run time of every file for matlab barchart
	const string& const str_14 = "./191120KNNMatlab/FileMethodApproximationTime" + str_suffix;//191125 size is total file number. total run time of every file for matlab barchart
	const string& const str_51 = "./191120KNNMatlab/FileMethodBuildTreeTime" + str_suffix;//210906
	const string& const str_15 = "./191120KNNMatlab/FileMethodOnlyKNNTime" + str_suffix;//191125 size is total file number. total run time of every file for matlab barchart
	const string& const str_47 = "./191120KNNMatlab/FileMethodOnlyKNNTimeHasIO" + str_suffix;//210606 size is total file number. total run time of every file for matlab barchart
	const string& const str_54 = "./191120KNNMatlab/FileMethodIOCost" + str_suffix;//211203 IO cost 

	/*----------------------------------------------------------------------------------------------------------------------------------------------*/

	/*-----------------------------------------              Split Point Coefficients           ----------------------------------------------------*/
	//200227 1 local sum deviaton, 2 global sum deviation, 3 local time, 4 global time, 5 lcoal shit, 6 global knn prune power
	local_total_split_id_sum_deviation.insert(local_total_split_id_sum_deviation.end(), global_total_approximation_sum_deviation.begin(), global_total_approximation_sum_deviation.end());
	local_total_split_id_sum_deviation.insert(local_total_split_id_sum_deviation.end(), local_total_split_id_time.begin(), local_total_split_id_time.end());
	local_total_split_id_sum_deviation.insert(local_total_split_id_sum_deviation.end(), global_total_approximation_time.begin(), global_total_approximation_time.end());
	local_total_split_id_sum_deviation.insert(local_total_split_id_sum_deviation.end(), local_total_split_id_shift.begin(), local_total_split_id_shift.end());
	local_total_split_id_sum_deviation.insert(local_total_split_id_sum_deviation.end(), global_total_knn_prune_power.begin(), global_total_knn_prune_power.end());
	const string& const str_24 = "./200227SplitPointKNNCoefficients/SplitIDCoefficients" + str_suffix;
	/*----------------------------------------------------------------------------------------------------------------------------------------------*/


	representation_option_vector.erase(representation_option_vector.begin());
	
	const string& const read_me = "Tree type: " + to_string(option_tree_struct.type_tree) + " Data Type: " + to_string(data_source.data_type) + " Data Attribute: " + dataset_attribute[data_source.data_type] + " Total file number: " + to_string(data_source.file_operation_number) + ". Total N number: " + to_string(N_coefficient_vector.size()) + " Total Method number: " + to_string(method_total_number)
		+ " Total K number: " + to_string(K_coefficient_vector.size()) + ". Point number: " + to_string(data_source.point_number)
		+ ". Rtree  max sub node: " + to_string(max_node0) + " Data dimension: " + to_string(data_source.time_series_dimension) + " Query time series size: " + to_string(Result_Coefficient_Struct.id_query_series_vector.size())
		+ ". mixed dataset: " + TOOL::convert_vector_to_string(data_source.data_list)
		+ ". file data size: " + TOOL::convert_vector_to_string(data_source.time_series_number_vector)
		+ ". N coefficient: " + TOOL::convert_vector_to_string(N_coefficient_vector)
		+ ". K coefficient: " + TOOL::convert_vector_to_string(K_coefficient_vector)
		+ ". File id: " + TOOL::convert_vector_to_string(file_id_choose_vector)
		+ ". Method id: " + TOOL::convert_vector_to_string(representation_option_vector)
		+ ". query time series vector id: " + TOOL::convert_vector_to_string(Result_Coefficient_Struct.id_query_series_vector)
		+ ". Index(Tree) type: " + TOOL::convert_vector_to_string(type_index_vector)
		+ notice_string;

	
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&        Clear Memory        &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	N_coefficient_vector.clear();
	N_coefficient_vector.shrink_to_fit();
	K_coefficient_vector.clear();
	K_coefficient_vector.shrink_to_fit();
	total_prune_power_vector.clear();//191119
	total_prune_power_vector.shrink_to_fit();
	total_prune_power_combine_vector.clear();//201221
	total_prune_power_combine_vector.shrink_to_fit();
	total_sum_deviation_vector.clear();//191206
	total_sum_deviation_vector.shrink_to_fit();
	total_accuracy_vector.clear();//191204
	total_accuracy_vector.shrink_to_fit();
	total_run_time_vector.clear();//191119
	total_run_time_vector.shrink_to_fit();
	total_run_time_has_IO_vector.clear();//210606
	total_run_time_has_IO_vector.shrink_to_fit();
	total_approximation_time_vector.clear();
	total_approximation_time_vector.shrink_to_fit();
	total_build_tree_time_vector.clear();
	total_build_tree_time_vector.shrink_to_fit();
	total_ingest_data_time_vector.clear();
	total_ingest_data_time_vector.shrink_to_fit();
	total_knn_time_vector.clear();
	total_knn_time_vector.shrink_to_fit();
	total_knn_time_has_IO_vector.clear();//210606
	total_knn_time_has_IO_vector.shrink_to_fit();

	transfer_total_prune_power_vector.clear();//200103
	transfer_total_prune_power_combine_vector.clear();//201221
	transfer_total_sum_deviation_vector.clear();//200103
	transfer_total_accuracy_vector.clear();//200103
	transfer_total_whole_run_time_vector.clear();//200103
	transfer_total_whole_run_time_has_IO_vector.clear();//210606
	transfer_total_approximation_run_time_vector.clear();//200108
	transfer_total_build_tree_time_vector.clear();//200108
	transfer_total_ingest_data_time_vector.clear();//211210
	transfer_total_knn_cpu_time_vector.clear();
	transfer_total_knn_run_time_vector.clear();//200108
	transfer_total_knn_run_time_has_IO_vector.clear();//210606
	transfer_total_prune_power_vector.shrink_to_fit();//200103
	transfer_total_prune_power_combine_vector.shrink_to_fit();//201221
	transfer_total_sum_deviation_vector.shrink_to_fit();//200103
	transfer_total_accuracy_vector.shrink_to_fit();//200103
	transfer_total_whole_run_time_vector.shrink_to_fit();//200103
	transfer_total_whole_run_time_has_IO_vector.shrink_to_fit();//210606
	transfer_total_approximation_run_time_vector.shrink_to_fit();//200108
	transfer_total_build_tree_time_vector.shrink_to_fit();
	transfer_total_ingest_data_time_vector.shrink_to_fit();//211210
	transfer_total_knn_cpu_time_vector.shrink_to_fit();
	transfer_total_knn_run_time_vector.shrink_to_fit();//200108
	transfer_total_knn_run_time_has_IO_vector.shrink_to_fit();//210606

	method_sum_prune_power_vector.clear();//191123 totl prune power of every method for matlab barchart
	method_sum_prune_power_vector.shrink_to_fit();
	method_sum_prune_power_combine_vector.clear();//201221
	method_sum_prune_power_combine_vector.shrink_to_fit();
	method_sum_sum_deviation_vector.clear();
	method_sum_sum_deviation_vector.shrink_to_fit();
	method_sum_accuracy_vector.clear();//191204 accuracy of prune power
	method_sum_accuracy_vector.shrink_to_fit();
	method_run_time_vector.clear();//191123 totl run time of every method for matlab barchart
	method_run_time_vector.shrink_to_fit();
	method_run_time_has_IO_vector.clear();//210606 totl run time of every method for matlab barchart
	method_run_time_has_IO_vector.shrink_to_fit();
	method_approximation_time_vector.clear();//191204 For approximation time of every method
	method_approximation_time_vector.shrink_to_fit();
	method_build_tree_time_vector.clear();//191204 For approximation time of every method
	method_build_tree_time_vector.shrink_to_fit();
	method_knn_time_vector.clear();//191204 For knn time of every method
	method_knn_time_vector.shrink_to_fit();
	method_knn_time_has_IO_vector.clear();//210606For knn time of every method
	method_knn_time_has_IO_vector.shrink_to_fit();

	method_file_prune_power_vector.clear();//191125 5*23 sum prune power for every file [mehtod number, file number]
	method_file_prune_power_vector.shrink_to_fit();
	method_file_prune_power_combine_vector.clear();//201221
	method_file_prune_power_combine_vector.shrink_to_fit();
	method_file_sum_deviation_vector.clear();
	method_file_sum_deviation_vector.shrink_to_fit();
	method_file_accuracy_vector.clear();//191204 accuracy of prune power
	method_file_accuracy_vector.shrink_to_fit();
	method_file_run_time_vector.clear();//191125 5*23 sum run time for every file [mehtod number, file number]
	method_file_run_time_vector.shrink_to_fit();
	method_file_run_time_has_IO_vector.clear();//210606 5*23 sum run time for every file [mehtod number, file number]
	method_file_run_time_has_IO_vector.shrink_to_fit();
	method_file_approximation_time_vector.clear();//191204 5*23 sum approximation time for every file [mehtod number, file number]
	method_file_approximation_time_vector.shrink_to_fit();
	method_file_build_tree_time_vector.clear();//191204 5*23 sum approximation time for every file [mehtod number, file number]
	method_file_build_tree_time_vector.shrink_to_fit();
	method_file_knn_time_vector.clear();//191204 5*23 sum knn time for every file [mehtod number, file number]
	method_file_knn_time_vector.shrink_to_fit();
	method_file_knn_time_has_IO_vector.clear();//210606 5*23 sum knn time for every file [mehtod number, file number]
	method_file_knn_time_has_IO_vector.shrink_to_fit();

	file_id_choose_vector.clear();//201020
	file_id_choose_vector.shrink_to_fit();

	/*========================    clear vector of split coefficents    ======================*/
	local_total_split_id_sum_deviation.clear();
	local_total_split_id_sum_deviation.shrink_to_fit();
	local_total_split_id_shift.clear();
	local_total_split_id_shift.shrink_to_fit();
	local_total_split_id_time.clear();
	local_total_split_id_time.shrink_to_fit();
	global_total_knn_prune_power.clear();
	global_total_knn_prune_power.shrink_to_fit();
	global_total_approximation_sum_deviation.clear();
	global_total_approximation_sum_deviation.shrink_to_fit();
	global_total_approximation_time.clear();
	global_total_approximation_time.shrink_to_fit();
	/*========================================================================================*/

	/*========================    201028 Vector of One File result    ======================*/
	one_file_prune_power_vector.clear();
	one_file_prune_power_vector.shrink_to_fit();
	one_file_prune_power_combine_vector.clear();
	one_file_prune_power_combine_vector.shrink_to_fit();
	one_file_sum_deviation_vector.clear();
	one_file_sum_deviation_vector.shrink_to_fit();
	one_file_accuracy_vector.clear();
	one_file_accuracy_vector.shrink_to_fit();
	one_file_run_time_vector.clear();
	one_file_run_time_vector.shrink_to_fit();
	one_file_approximation_time_vector.clear();
	one_file_approximation_time_vector.shrink_to_fit();
	one_file_knn_time_vector.clear();
	one_file_knn_time_vector.shrink_to_fit();
	/*========================================================================================*/
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	//system("pause");
}


TEMPLATE
template<typename T, typename Y, typename U>
inline void Evaluation::print_input_output(const T& const input_argument, const Y& const output_argument, const U& const linear_scan_struct) {
	
	cout << "*prune power: " << input_argument.pruning_power << endl;
	cout << "max deviation: " << output_argument.max_deviation << endl;
	cout << "prune accuracy: " << input_argument.result_accuracy << endl;
	cout << "representation_time (us): " << input_argument.representation_time << endl;
	cout << "ingest data time (us): " << input_argument.time_ingest_data << endl;
	cout << "*knn_total_time (us): " << input_argument.knn_total_time << endl;
	cout << "*Tree height: " << input_argument.argument_index_struct.height_index << endl;
	cout << "*Internal node number: " << input_argument.argument_index_struct.count_node_internal << endl;
	cout << "*Leaf node number: " << input_argument.argument_index_struct.count_node_leaf << endl;
	cout << "*Total node number: " << input_argument.argument_index_struct.count_node_total << endl;
}

#endif

