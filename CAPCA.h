#pragma once
#ifndef CAPCA_H
#define CAPCA_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "./lib/doublyLinkedList.h"

//***************************************************************
// Class:APCA & PAA Approximation
// Qualifier:
// Input:
// Output:
// date:
// author:
//***************************************************************
TEMPLATE
class CAPCA : virtual public TOOL
{
public:

	struct INPUT_ARGUMENT;
	struct APCA;
	struct COEFFICIENT_PAIR;
	struct cmpLess;
	struct APCA_MERGE;
	struct cmpMoreDistance;
	struct cmpLessDistance;
	struct cmpMoreIndex;
	struct APCA_ARRAY;

	struct TOOL::INPUT_ARGUMENT input_argument;
	struct TOOL::OUTPUT_ARGUMENT output_argument;//190402

	CAPCA(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const int& arity_d, const string& read_file_name, string*& const write_file_name);
	CAPCA(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const int& arity_d, string*& const write_file_name);
	CAPCA(const int& n, const int& N);
	CAPCA() {};
	~CAPCA();// {};

	/*========================================================================APCA======================================================================================*/
	inline unsigned int getNextPowerOf2(const unsigned int& array_length0);
	void padZero(DataType*& const temp_original_array, const DataType& old_length, const unsigned int& new_length, DataType* after_padding_array);
	void padZeroVector(DataType*& const temp_original_array, const DataType& old_length, const unsigned int& new_length, vector<DataType>& after_padding_array);//190430
	void padZeroVectorPreMemory(DataType*& const temp_original_array, const DataType& old_length, const unsigned int& new_length, vector<DataType>& after_padding_array);//190611
	bool truncateZero(const DataType& old_length, const unsigned int& new_length, APCA& apca_presentation);
	bool truncateZeroVector(const DataType& old_length, const unsigned int& new_length, vector<APCA>& apca_presentation);//190501
	void getHDWT(const unsigned int& length_power_of_2, const int& N, DataType* wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index);
	void getHDWTVector(const unsigned int& length_power_of_2, const int& N, vector<DataType>& wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index);//190430
	void reconstructApproximateAPCA(const DataType* wavelet_transform_time_series, const unsigned int& power_of_2, const int& retained_coeffs_length, priority_queue<DataType>& fq_truncate_index, APCA& apca_presentation);
	void reconstructApproximateAPCAVector(const vector<DataType>& wavelet_transform_time_series, const unsigned int& power_of_2, const int& retained_coeffs_length, priority_queue<DataType>& fq_truncate_index, vector<APCA>& apca_presentation);//190430
	void reconstructApproximateAPCAVectorPreMemory(const vector<DataType>& wavelet_transform_time_series, const unsigned int& power_of_2, const int& retained_coeffs_length, priority_queue<DataType>& fq_truncate_index, vector<APCA>& apca_presentation);//190611
	/*===================================================================================================================================================================*/

	DataType inline getAve(const DataType* originalArray, const DataType& arrayLength);
	bool getExactMeanValue(const DataType* original_time_series, APCA& apca_presentation);
	bool getExactMeanValueVector(const DataType* original_time_series, vector<APCA>& apca_presentation);//190501
	bool mergeSegmentsIterator(APCA& apca_presentation, const DataType& n, const int& N, APCA& italicC);
	bool mergeSegmentsIteratorVector(vector<APCA>& const apca_presentation, const DataType& n, const int& N, vector<APCA>& const italicC);//190501
	template<typename T>
	bool mergeSegmentsIteratorLinkedList(vector<T>& const apca_presentation, const DataType& n, const int& N, DoublyLinkedList<T>& const italicC);//191104 Use linked list ot instead vector
	bool mergeSegmentsIteratorVectorPreMemory(vector<APCA>& const apca_presentation, const DataType& n, const int& N, vector<APCA>& const italicC);//190611

	void getMergeIndex(const int& merge_frequency, APCA& apca_presentation);
	void getMergeIndexVector(const int& merge_frequency, vector<APCA>& apca_presentation);//190501
	/*==============================================APCA Approximaiton===========================================================================================*/
	void getAPCAPoint(DataType*& const original_time_series, const double& n, const int& N, APCA& italicC);
	void getAPCAPointVector(DataType*& const original_time_series, const double& n, const int& N, vector<APCA>& const italicC);//190430 vector to instead Array
	template<typename T>
	void getAPCAPointLinkedList(DataType*& const original_time_series, const double& n, const int& N, DoublyLinkedList<T>& const italicC);//191104 Linked list to instead vector
	void getAPCAPointVectorPreMemory(DataType*& const original_time_series, const double& n, const int& N, vector<APCA>& const italicC);//190611 Pre-define memory of vector
	/*===========================================================================================================================================================*/
	/*================================================PAA Approximaiton==========================================================================================*/
	void divideRemainderPAA(const DataType* original_time_series, const APCA& italicC, const DataType& n, const int& N);

	//200923 use vector instead pointer
	template<typename T, typename Y, typename U>
	void get_PAA(const vector<T>& const original_time_series_vector, const Y& const segment_number, U& const paa_container);

	void getPAAPointVector(const DataType* original_time_series, const DataType& n, const int& N, vector<APCA>& const italicC);//190515 Only use vecter. Same with divideRemainderPAA()
	template<typename T>//190515 use Linked list to instead vector. Same with divideRemainderPAA()
	void getPAAPointLinkedList(const DataType* original_time_series, const DataType& n, const int& N, DoublyLinkedList<T>& const italicC);//
	void getPAAPointVectorPreMemory(const DataType* original_time_series, const DataType& n, const int& N, vector<APCA>& const italicC);//190611 Only use vecter. Same with divideRemainderPAA(), pre-define memory
	/*===========================================================================================================================================================*/

	/*================================================200108 Lagrangian representation algorithm==========================================================================================*/
	//¡¶A New Pattern Representation Method for Time - series Data¡·
	template<typename T>//190515 use Linked list to instead vector. Same with divideRemainderPAA()
	void get_paa_lagrangian(T& const italicC); 
	/*=================================================================================================================================================================================*/

	void approximateOriginalFunction(const APCA_QUAL::APCA& italicC, DataType*& const approximated_time_series);

	void initialAPCAArray(const int& mg_d_index_point_number, const int& mg_APCA_point_dimension, APCA_ARRAY*& const APCA_array);
	void initialAPCAArray(const int& point_number, const int& point_dimension, APCA*& const APCA_array);
	void initialAPCAArray(const typename TOOL::INPUT_ARGUMENT& input_argument, APCA_ARRAY*& const APCA_array);
	//191118
	void initial_rect_vector(const int& const point_number, const int& const point_dimension, vector<RTREE::Rect>& const rtree_rectangle_vector);
	void deleteAPCAArray(const DataType& g_index_point_number, APCA_ARRAY*& const APCA_array);
	void deleteAPCAArray(const DataType& g_index_point_number, APCA*& const APCA_array);
	void initialAPCA(APCA& apca, int array_length);
	void deleteAPCA(APCA& apca);
	void delete_rect_vector(vector<RTREE::Rect>& const rtree_rectangle_vector);

	APCA& getCmax(const double* C, const APCA& italicC, APCA& Cmax);
	APCA& getCmin(const double* C, const APCA& italicC, APCA& Cmin);

	//191118 22:15 original MBR min id <= max id 
	template<typename T, typename Y>
	Y& get_minmax_original(const vector<DataType>& const time_series_vector, const int& const dimension_id, const T& const italicC, Y& const rtree_rectangle);

	//191126 22:19  APCA paper. MBR min id == max id
	template<typename T, typename Y>
	Y& get_minmax_apca(const vector<DataType>& const time_series_vector, const int& const dimension_id, const T& const italicC, Y& const rtree_rectangle);

	//191128 23:55  APCA paper. MBR min id == max id
	template<typename T, typename Y>
	Y& get_minmax_apca(const vector<DataType>& const time_series_vector, const T& const italicC, Y& const rtree_rectangle);

	APCA& getMBR(const APCA& Cmin, const APCA& Cmax, APCA& MBR);

	RTREE& buidRTreeIndex(RTREE& const APCARTree, const DataType& g_time_series_length, const DataType& g_index_point_number, APCA_ARRAY* APCALinkOriginal, const string& file_name);

	//friend class CPACA_KNN;
private:

};

//TEMPLATE
//APCA_QUAL::~CAPCA() {
//}

//TEMPLATE
//struct APCA_QUAL::APCA {
//	DataType* v = nullptr;
//	DataType* r = nullptr;
//	int segmentNum = INF; //the dimension of the index(eg. MBR, imput parameter)
//	//191108
//	APCA() :v(nullptr), r(nullptr), segmentNum(INF) {};
//
//	~APCA() {
//		/*if (v != nullptr) {
//			delete[] v;
//			v = nullptr;
//		}
//		if (r != nullptr) {
//			delete[] r;
//			r = nullptr;
//		}*/
//		segmentNum = INF;
//	}
//};

//TEMPLATE
//struct APCA_QUAL::APCA_ARRAY {
//	//DataType *originalLink = nullptr;
//	//int original_time_series_id = NULL;
//	APCA APCALink;
//};

//TEMPLATE
//void  APCA_QUAL::initialAPCAArray(const int& mg_d_index_point_number, const int& mg_APCA_point_dimension, APCA_ARRAY*& const APCA_array) {
//
//	APCA_array = new APCA_ARRAY[mg_d_index_point_number];
//
//	for (int i = 0; i < int(mg_d_index_point_number); i++) {
//		APCA_array[i].APCALink.segmentNum = mg_APCA_point_dimension;
//		APCA_array[i].APCALink.r = new DataType[mg_APCA_point_dimension];
//		APCA_array[i].APCALink.v = new DataType[mg_APCA_point_dimension];
//		//fill(APCA_array[i].APCALink.r, APCA_array[i].APCALink.r + mg_APCA_point_dimension, NULL);
//		//fill(APCA_array[i].APCALink.v, APCA_array[i].APCALink.v + mg_APCA_point_dimension, NULL);
//		fill_n(APCA_array[i].APCALink.r, mg_APCA_point_dimension, INF);
//		fill_n(APCA_array[i].APCALink.v, mg_APCA_point_dimension, INF);
//	}
//}

//#include "CAPCA.cpp"

//template class CAPCA<DataType>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///*********************************************************************************************************************************************************

TEMPLATE
struct APCA_QUAL::INPUT_ARGUMENT {
	DataType mg_file_time_series_length = NULL;
	DataType mg_d_index_point_number = NULL;
	int mg_APCA_point_dimension = NULL;
	int mg_max_nodes = NULL;
	int K = NULL;
	string read_file_name;
};

TEMPLATE
struct APCA_QUAL::APCA {
	DataType* v = nullptr;
	DataType* r = nullptr;
	int segmentNum = INF; //the dimension of the index(eg. MBR, imput parameter)
	//191108
	APCA() :v(nullptr), r(nullptr), segmentNum(INF) {};

	~APCA() {
		/*if (v != nullptr) {
			delete[] v;
			v = nullptr;
		}
		if (r != nullptr) {
			delete[] r;
			r = nullptr;
		}*/
		segmentNum = INF;
	}
};

//Haaer wavelet transformation
TEMPLATE
struct APCA_QUAL::COEFFICIENT_PAIR {
	int HDWT_id = NULL;
	DataType HDWT_original_value = NULL;
	DataType HDWT_normalized_value = NULL;
	DataType HDWT_magnitude = NULL;
};

TEMPLATE
struct APCA_QUAL::cmpLess {
	bool operator () (const COEFFICIENT_PAIR& a, const COEFFICIENT_PAIR& b) {
		return a.HDWT_magnitude < b.HDWT_magnitude;             // from big to small
	}
};

TEMPLATE
struct APCA_QUAL::APCA_MERGE {
	DataType segment_differences = NULL;
	DataType* r = nullptr;
	DataType* v = nullptr;
};

TEMPLATE
struct APCA_QUAL::cmpMoreDistance {
	bool operator () (const APCA_MERGE& a, const APCA_MERGE& b) {
		return a.segment_differences > b.segment_differences;             // small to big
	}
};

TEMPLATE
struct APCA_QUAL::cmpLessDistance {
	bool operator () (const APCA_MERGE& a, const APCA_MERGE& b) {
		return a.segment_differences < b.segment_differences;             // big to small
	}
};

TEMPLATE
struct APCA_QUAL::cmpMoreIndex {
	bool operator () (const APCA_MERGE& a, const APCA_MERGE& b) {
		return *a.r > *b.r;             // small to big
	}
};

TEMPLATE
struct APCA_QUAL::APCA_ARRAY {
	//DataType *originalLink = nullptr;
	//int original_time_series_id = NULL;
	APCA APCALink;
};

TEMPLATE
APCA_QUAL::CAPCA(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const int& arity_d, const string& read_file_name, string*& const write_file_name) {
	assert(arity_d > 0);
	assert(point_number >= K);

	input_argument.time_series_length = n;
	input_argument.degree_m = N;//point_dimension for APCA & PLA n=m+1
	input_argument.point_dimension = N;
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

	/*initialCHEBYSHEV_SHARE(input_argument, chebyshev_share);
	getCHEBYSHEV_SHARE(input_argument, chebyshev_share);*/
}

TEMPLATE
APCA_QUAL::CAPCA(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const int& arity_d, string*& const write_file_name) {
	assert(arity_d > 0);
	assert(point_number >= K);

	input_argument.time_series_length = n;
	input_argument.degree_m = N;//point_dimension for APCA & PLA n=m+1
	input_argument.point_dimension = N;
	input_argument.point_number = point_number;
	input_argument.rtree_max_nodes = rtree_max_nodes;
	input_argument.K = K;
	input_argument.arity_d = arity_d;
	//input_argument.read_file_name = read_file_name;
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

	/*initialCHEBYSHEV_SHARE(input_argument, chebyshev_share);
	getCHEBYSHEV_SHARE(input_argument, chebyshev_share);*/
}

TEMPLATE
APCA_QUAL::CAPCA(const int& n, const int& N) {
	input_argument.remainder = int(n) % int(N);
	double integerDividend = n - input_argument.remainder;
	input_argument.segment_length_second = integerDividend / N;
	input_argument.segment_length_first = input_argument.segment_length_second + 1;
	assert(input_argument.segment_length_second > 1);//l(l-1)(l+1), so l != 1
	input_argument.time_series_length = n;
	input_argument.point_dimension = N;
}

TEMPLATE
APCA_QUAL::~CAPCA() {
}

TEMPLATE
inline unsigned int APCA_QUAL::getNextPowerOf2(const unsigned int& array_length0) {
	unsigned int array_length = array_length0;

	if (array_length <= 0) return 0;

	array_length--;
	array_length |= array_length >> 1;
	array_length |= array_length >> 2;
	array_length |= array_length >> 4;
	array_length |= array_length >> 8;
	array_length |= array_length >> 16;

	return array_length + 1;
}

TEMPLATE
void APCA_QUAL::padZero(DataType*& const temp_original_array, const DataType& old_length, const unsigned int& new_length, DataType* after_padding_array) {
	assert(new_length >= old_length);
	fill_n(after_padding_array, new_length, 0.0);
	copy(temp_original_array, temp_original_array + int(old_length), after_padding_array);
}

//************************************
// Method:padZeroVector
// Qualifier:190430 vector to instead Array
// date:190430
// author:
//************************************
TEMPLATE
void APCA_QUAL::padZeroVector(DataType*& const temp_original_array, const DataType& old_length, const unsigned int& new_length, vector<DataType>& after_padding_array) {//190430
	assert(new_length >= old_length);

	for (int i = 0; i < new_length; i++) {
		if (i < old_length)
			after_padding_array.push_back(temp_original_array[i]);
		else
			after_padding_array.push_back(0.0);
	}
}

//************************************
// Method:padZeroVectorPreMemory
// Qualifier:190611 vector to instead Array
// date:190611
// author:
//************************************
TEMPLATE
void APCA_QUAL::padZeroVectorPreMemory(DataType*& const temp_original_array, const DataType& old_length, const unsigned int& new_length, vector<DataType>& after_padding_array) {//190611
#ifdef _DEBUG
	assert(new_length >= old_length);
#endif

	fill_n(after_padding_array.begin(), new_length, 0.0);
	copy(temp_original_array, temp_original_array + int(old_length), after_padding_array.begin());
}

TEMPLATE
bool APCA_QUAL::truncateZero(const DataType& old_length, const unsigned int& new_length, APCA& apca_presentation) {
	if (new_length == old_length) {
		return true;
	}
	assert(new_length > old_length);
	//assert(getNextPowerOf2(unsigned int(old_length)) == new_length);

	int f_segments_last_point = apca_presentation.segmentNum - 1;
	int f_old_last_point = old_length - 1;

	for (int i = f_segments_last_point; i >= 0; i--) {
		if (f_old_last_point <= apca_presentation.r[i] && f_old_last_point > apca_presentation.r[i - 1]) {
			//cout <<"***"<<i<<" "<< apca_presentation.r[i] << endl;
			apca_presentation.r[i] = f_old_last_point;
			apca_presentation.segmentNum = i + 1;
			//cout << apca_presentation.r[i] << endl;
			return true;
		}
	}
}

TEMPLATE
bool APCA_QUAL::truncateZeroVector(const DataType& old_length, const unsigned int& new_length, vector<APCA>& const apca_presentation) {//190501
	if (new_length == old_length) {
		return true;
	}
#ifdef _DEBUG
	assert(new_length > old_length);
#endif

	//assert(getNextPowerOf2(unsigned int(old_length)) == new_length);
	int truncated_segment_size = INF;

	//int f_segments_last_point = apca_presentation.segmentNum - 1;
	int f_segments_last_point = apca_presentation.size() - 1;
	int f_old_last_point = old_length - 1;

	for (int i = f_segments_last_point; i >= 0; i--) {
		if (f_old_last_point <= *apca_presentation[i].r && f_old_last_point > *apca_presentation[i - 1].r) {
			//cout <<"***"<<i<<" "<< apca_presentation.r[i] << endl;
			*apca_presentation[i].r = f_old_last_point;
			//apca_presentation.segmentNum = i + 1;
			truncated_segment_size = i + 1;
			//cout << apca_presentation.r[i] << endl;
		}
	}

	while (truncated_segment_size < apca_presentation.size()) {
		delete apca_presentation.back().v;
		apca_presentation.back().v = nullptr;
		delete apca_presentation.back().r;
		apca_presentation.back().r = nullptr;
		apca_presentation.pop_back();
	}

	return true;
}

//Haar wavelet transformation
TEMPLATE
void APCA_QUAL::getHDWT(const unsigned int& length_power_of_2, const int& N, DataType* wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index) {//Perform the Haar Discrete Wavelet Transform on C
	//cout << power_of_2 << endl;
	int f_interation_times = int(log2(length_power_of_2)) - 1;//3
	//cout << "f_interation_times: " << f_interation_times << endl;
	int temp_power_of_2 = length_power_of_2 >> 1;

	DataType* fd_temp = new DataType[length_power_of_2];
	priority_queue<COEFFICIENT_PAIR, vector<COEFFICIENT_PAIR>, cmpLess> fq_HDWT_coefficients;

	COEFFICIENT_PAIR temp_Coefficient;

	/*for (int i = 0; i < int(length_power_of_2); i++) {
		fd_temp[i] = wavelet_transform_time_series[i];
	}*************************************/
	copy(wavelet_transform_time_series, wavelet_transform_time_series + length_power_of_2, fd_temp);

	for (int i = f_interation_times; i >= 0; i--) {//3
		for (int j = 0; j < temp_power_of_2; j++) {
			wavelet_transform_time_series[j] = (fd_temp[j << 1] + fd_temp[1 + (j << 1)]) / 2.0;
			fd_temp[j] = (fd_temp[j << 1] + fd_temp[1 + (j << 1)]) / 2.0;
			wavelet_transform_time_series[j + temp_power_of_2] = wavelet_transform_time_series[j] - fd_temp[1 + (j << 1)];

			temp_Coefficient.HDWT_original_value = wavelet_transform_time_series[j + temp_power_of_2];
			temp_Coefficient.HDWT_normalized_value = temp_Coefficient.HDWT_original_value / pow(2.0, double(i) / 2.0);
			temp_Coefficient.HDWT_magnitude = fabs(temp_Coefficient.HDWT_normalized_value);
			temp_Coefficient.HDWT_id = j + temp_power_of_2;

			//cout << "(" << temp_Coefficient.HDWT_original_value << "," << temp_Coefficient.HDWT_magnitude << ") ";
			/*for (int i = 0; i < length_power_of_2; i++) {
			cout<<wavelet_transform_time_series[i] << " ";
			}
			cout << endl;*/
			/*for (int i = 0; i < length_power_of_2; i++) {
			cout << fd_temp[i] << " ";
			}
			cout << endl;*/

			fq_HDWT_coefficients.push(temp_Coefficient);
		}
		temp_power_of_2 >>= 1;
		//cout << endl;
	}

	temp_Coefficient.HDWT_original_value = *wavelet_transform_time_series;
	temp_Coefficient.HDWT_normalized_value = *wavelet_transform_time_series;
	temp_Coefficient.HDWT_magnitude = abs(temp_Coefficient.HDWT_normalized_value);
	temp_Coefficient.HDWT_id = 0;
	fq_HDWT_coefficients.push(temp_Coefficient);

	for (int i = 0; i < N; i++) {
		//cout << " Big magnitude(ID): " << fq_HDWT_coefficients.top().HDWT_id << endl;
		fq_truncate_index.push(fq_HDWT_coefficients.top().HDWT_id);
		fq_HDWT_coefficients.pop();
	}

	priority_queue<COEFFICIENT_PAIR, vector<COEFFICIENT_PAIR>, cmpLess>().swap(fq_HDWT_coefficients);
	delete[] fd_temp;
	fd_temp = nullptr;
}

//Haar wavelet transformation
//************************************
// Method:getHDWTVector
// Qualifier:190430 vector to instead Array
// date:190430
// author:
//************************************
TEMPLATE
void APCA_QUAL::getHDWTVector(const unsigned int& length_power_of_2, const int& N, vector<DataType>& wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index) {//190430
//cout << power_of_2 << endl;
#ifdef _DEBUG
	assert(length_power_of_2 == wavelet_transform_time_series.size());
#endif

	int f_interation_times = int(log2(length_power_of_2)) - 1;//3
	//cout << "f_interation_times: " << f_interation_times << endl;
	int temp_power_of_2 = length_power_of_2 >> 1;

	//DataType* fd_temp = new DataType[length_power_of_2];
	vector<DataType> fd_temp;
	//fd_temp.resize(length_power_of_2, 0.0);//190430
	priority_queue<COEFFICIENT_PAIR, vector<COEFFICIENT_PAIR>, cmpLess> fq_HDWT_coefficients;

	COEFFICIENT_PAIR temp_Coefficient;

	/*for (int i = 0; i < int(length_power_of_2); i++) {
		fd_temp[i] = wavelet_transform_time_series[i];
	}*************************************/
	//copy(wavelet_transform_time_series.begin(), wavelet_transform_time_series.end(), fd_temp.begin());
	for (auto&& i : wavelet_transform_time_series) {
		fd_temp.push_back(i);
	}

	for (int i = f_interation_times; i >= 0; i--) {//3
		for (int j = 0; j < temp_power_of_2; j++) {
			wavelet_transform_time_series[j] = (fd_temp[j << 1] + fd_temp[1 + (j << 1)]) / 2.0;
			fd_temp[j] = (fd_temp[j << 1] + fd_temp[1 + (j << 1)]) / 2.0;
			wavelet_transform_time_series[j + temp_power_of_2] = wavelet_transform_time_series[j] - fd_temp[1 + (j << 1)];

			temp_Coefficient.HDWT_original_value = wavelet_transform_time_series[j + temp_power_of_2];
			temp_Coefficient.HDWT_normalized_value = temp_Coefficient.HDWT_original_value / pow(2.0, double(i) / 2.0);
			temp_Coefficient.HDWT_magnitude = fabs(temp_Coefficient.HDWT_normalized_value);
			temp_Coefficient.HDWT_id = j + temp_power_of_2;

			//cout << "(" << temp_Coefficient.HDWT_original_value << "," << temp_Coefficient.HDWT_magnitude << ") ";
			/*for (int i = 0; i < length_power_of_2; i++) {
			cout<<wavelet_transform_time_series[i] << " ";
			}
			cout << endl;*/
			/*for (int i = 0; i < length_power_of_2; i++) {
			cout << fd_temp[i] << " ";
			}
			cout << endl;*/

			fq_HDWT_coefficients.push(temp_Coefficient);
		}
		temp_power_of_2 >>= 1;
		//cout << endl;
	}

	temp_Coefficient.HDWT_original_value = wavelet_transform_time_series.front();
	temp_Coefficient.HDWT_normalized_value = wavelet_transform_time_series.front();
	temp_Coefficient.HDWT_magnitude = abs(temp_Coefficient.HDWT_normalized_value);
	temp_Coefficient.HDWT_id = 0;
	fq_HDWT_coefficients.push(temp_Coefficient);

	for (int i = 0; i < N; i++) {
		//cout << " Big magnitude(ID): " << fq_HDWT_coefficients.top().HDWT_id << endl;
		fq_truncate_index.push(fq_HDWT_coefficients.top().HDWT_id);
		fq_HDWT_coefficients.pop();
	}

	priority_queue<COEFFICIENT_PAIR, vector<COEFFICIENT_PAIR>, cmpLess>().swap(fq_HDWT_coefficients);
	//delete[] fd_temp;
	//fd_temp = nullptr;
}

TEMPLATE
void APCA_QUAL::reconstructApproximateAPCA(const DataType* wavelet_transform_time_series, const unsigned int& length_power_of_2, const int& retained_coeffs_length, priority_queue<DataType>& fq_truncate_index, APCA& apca_presentation) {
	//DataType temp_apca_value = NULL;
	int loop_times = log2(retained_coeffs_length);
	int inner_loop_times = 1;
	queue<DataType> temp_original_queue;

	fill_n(apca_presentation.v, retained_coeffs_length, 0.0);
	fill_n(apca_presentation.r, retained_coeffs_length, length_power_of_2);

	for (int i = 0; i < int(fq_truncate_index.size()); i++) {
		apca_presentation.v[int(fq_truncate_index.top())] = wavelet_transform_time_series[int(fq_truncate_index.top())];
		fq_truncate_index.pop();
	}

	temp_original_queue.push(*wavelet_transform_time_series);

	for (int i = 0; i < loop_times; i++) {
		for (int j = 0; j < inner_loop_times; j++) {
			apca_presentation.v[(j << 1) + 1] = temp_original_queue.front() - apca_presentation.v[j + inner_loop_times];
			apca_presentation.v[j << 1] = 2 * temp_original_queue.front() - apca_presentation.v[(j << 1) + 1];

			apca_presentation.r[j << 1] = ((j << 1) + 1) * length_power_of_2 / (inner_loop_times << 1) - 1;
			apca_presentation.r[(j << 1) + 1] = apca_presentation.r[j << 1] + length_power_of_2 / (inner_loop_times << 1);

			temp_original_queue.pop();
			temp_original_queue.push(apca_presentation.v[j << 1]);
			temp_original_queue.push(apca_presentation.v[(j << 1) + 1]);
		}
		inner_loop_times <<= 1;
	}

	queue<DataType>().swap(temp_original_queue);
}

//************************************
// Method:reconstructApproximateAPCAVector
// Qualifier:190430 vector to instead Array
// date:190430
// author:
//************************************
TEMPLATE
void APCA_QUAL::reconstructApproximateAPCAVector(const vector<DataType>& wavelet_transform_time_series, const unsigned int& length_power_of_2, const int& retained_coeffs_length, priority_queue<DataType>& fq_truncate_index, vector<APCA>& apca_presentation) {//190430
	//DataType temp_apca_value = NULL;
	int loop_times = log2(retained_coeffs_length);
	int inner_loop_times = 1;
	queue<DataType> temp_original_queue;

	APCA temp_paca;

	/*========================Vector instead Array=======================================*/
	//fill_n(apca_presentation.v, retained_coeffs_length, 0.0);
	//fill_n(apca_presentation.r, retained_coeffs_length, length_power_of_2);
	for (int i = 0; i < retained_coeffs_length; i++) {
		temp_paca.v = new DataType;
		temp_paca.r = new DataType;
		*temp_paca.v = 0.0;
		*temp_paca.r = length_power_of_2;
		apca_presentation.push_back(temp_paca);
	}
	/*.................................................................*/

	for (int i = 0; i < int(fq_truncate_index.size()); i++) {
		*apca_presentation[int(fq_truncate_index.top())].v = wavelet_transform_time_series[int(fq_truncate_index.top())];
		fq_truncate_index.pop();
	}

	temp_original_queue.push(wavelet_transform_time_series.front());

	for (int i = 0; i < loop_times; i++) {
		for (int j = 0; j < inner_loop_times; j++) {
			*apca_presentation[(j << 1) + 1].v = temp_original_queue.front() - *apca_presentation[j + inner_loop_times].v;
			*apca_presentation[j << 1].v = 2 * temp_original_queue.front() - *apca_presentation[(j << 1) + 1].v;

			*apca_presentation[j << 1].r = ((j << 1) + 1) * length_power_of_2 / (inner_loop_times << 1) - 1;
			*apca_presentation[(j << 1) + 1].r = *apca_presentation[j << 1].r + length_power_of_2 / (inner_loop_times << 1);

			temp_original_queue.pop();
			temp_original_queue.push(*apca_presentation[j << 1].v);
			temp_original_queue.push(*apca_presentation[(j << 1) + 1].v);
		}
		inner_loop_times <<= 1;
	}

	queue<DataType>().swap(temp_original_queue);
}

//************************************
// Method:reconstructApproximateAPCAVectorPreMemory
// Qualifier:190611 vector to instead Array, pre define memory
// date:190611
// author:
//************************************
TEMPLATE
void APCA_QUAL::reconstructApproximateAPCAVectorPreMemory(const vector<DataType>& wavelet_transform_time_series, const unsigned int& length_power_of_2, const int& retained_coeffs_length, priority_queue<DataType>& fq_truncate_index, vector<APCA>& apca_presentation) {//190611
	//DataType temp_apca_value = NULL;
	int loop_times = log2(retained_coeffs_length);
	int inner_loop_times = 1;
	queue<DataType> temp_original_queue;

	/*========================Vector instead Array Predefine memory=======================================*/
	//fill_n(apca_presentation.v, retained_coeffs_length, 0.0);
	//fill_n(apca_presentation.r, retained_coeffs_length, length_power_of_2);
	APCA temp_paca;
	temp_paca.v = new DataType;
	temp_paca.r = new DataType;
	*temp_paca.v = 0.0;
	*temp_paca.r = length_power_of_2;
	/*for (int i = 0; i < retained_coeffs_length; i++) {
		temp_paca.v = new DataType;
		temp_paca.r = new DataType;
		*temp_paca.v = 0.0;
		*temp_paca.r = length_power_of_2;
		apca_presentation.push_back(temp_paca);
	}*/

	apca_presentation.resize(retained_coeffs_length, APCA());//190611
	for (int i = 0; i < retained_coeffs_length; i++) {
		//temp_paca.v = new DataType;
		//temp_paca.r = new DataType;
		apca_presentation[i].v = new DataType;
		apca_presentation[i].r = new DataType;
		//*apca_presentation[i].v = *temp_paca.v;
		//*apca_presentation[i].r = *temp_paca.r;
		*apca_presentation[i].v = 0.0;
		*apca_presentation[i].r = length_power_of_2;
	}
	/*.................................................................*/

	for (int i = 0; i < int(fq_truncate_index.size()); i++) {
		*apca_presentation[int(fq_truncate_index.top())].v = wavelet_transform_time_series[int(fq_truncate_index.top())];
		fq_truncate_index.pop();
	}

	temp_original_queue.push(wavelet_transform_time_series.front());

	for (int i = 0; i < loop_times; i++) {
		for (int j = 0; j < inner_loop_times; j++) {
			*apca_presentation[(j << 1) + 1].v = temp_original_queue.front() - *apca_presentation[j + inner_loop_times].v;
			*apca_presentation[j << 1].v = 2 * temp_original_queue.front() - *apca_presentation[(j << 1) + 1].v;

			*apca_presentation[j << 1].r = ((j << 1) + 1) * length_power_of_2 / (inner_loop_times << 1) - 1;
			*apca_presentation[(j << 1) + 1].r = *apca_presentation[j << 1].r + length_power_of_2 / (inner_loop_times << 1);

			temp_original_queue.pop();
			temp_original_queue.push(*apca_presentation[j << 1].v);
			temp_original_queue.push(*apca_presentation[(j << 1) + 1].v);
		}
		inner_loop_times <<= 1;
	}

	queue<DataType>().swap(temp_original_queue);
}

TEMPLATE
DataType inline APCA_QUAL::getAve(const DataType* originalArray, const DataType& arrayLength) {
	//cout << "accumulate: " << accumulate(originalArray, originalArray + int(arrayLength), 0) << ", " << arrayLength << endl;
	return accumulate(originalArray, originalArray + int(arrayLength), 0.0) / DataType(arrayLength);
}

TEMPLATE
bool APCA_QUAL::getExactMeanValue(const DataType* original_time_series, APCA& apca_presentation) {
	if (apca_presentation.segmentNum == 1) {
		return true;
	}

	double fd_segment_length = NULL;
	int i = 0;

	apca_presentation.v[i] = getAve(original_time_series, apca_presentation.r[0] + 1);

	for (i = 1; i < apca_presentation.segmentNum; i++) {
		apca_presentation.v[i] = getAve(original_time_series + int(apca_presentation.r[i - 1] + 1), double(apca_presentation.r[i] - apca_presentation.r[i - 1]));
		//cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;;
	}
	return true;
}

TEMPLATE
bool APCA_QUAL::getExactMeanValueVector(const DataType* original_time_series, vector<APCA>& apca_presentation) {//190501
	if (apca_presentation.size() == 1) {
		return true;
	}

	double fd_segment_length = NULL;
	int i = 0;

	*apca_presentation[i].v = getAve(original_time_series, *apca_presentation[0].r + 1);

	for (i = 1; i < apca_presentation.size(); i++) {
		*apca_presentation[i].v = getAve(original_time_series + int(*apca_presentation[i - 1].r + 1), double(*apca_presentation[i].r - *apca_presentation[i - 1].r));
		//cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;;
	}
	return true;
}

TEMPLATE
void APCA_QUAL::getMergeIndex(const int& merge_frequency, APCA& apca_presentation) {
	//void getMergeIndex(const DataType* original_time_series, const int& merge_frequency, APCA& apca_presentation) {
	DataType sefments_difference = NULL;
	int merge_start_index = NULL;
	double initial_segment_length = apca_presentation.r[0] + 1.0;

	for (int i = 0; i < merge_frequency; i++) { //merge times
		sefments_difference = DBL_MAX;

		apca_presentation.segmentNum--;
		//cout << "segment Number: " << apca_presentation.segmentNum << endl;
		for (int j = 0; j < apca_presentation.segmentNum; j++) {
			//cout << "sefments_difference: " << sefments_difference << " sgetments diff: " << fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]) << endl;
			if (sefments_difference > fabs(apca_presentation.v[j + 1] - apca_presentation.v[j])) {
				sefments_difference = fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]);
				merge_start_index = j;
			}
		}
		//cout << "merge_start_index: " << merge_start_index << endl;
		//apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] + apca_presentation.v[merge_start_index + 1]) / 2;
		if (merge_start_index != 0) {
			//double temp1 = getAve(original_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));
			//apca_presentation.v[merge_start_index] = getAve(original_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));

			apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] * (apca_presentation.r[merge_start_index] - apca_presentation.r[merge_start_index - 1]) + apca_presentation.v[merge_start_index + 1] * (apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index])) / double(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]);
		}
		else {
			apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] * double(apca_presentation.r[merge_start_index] + 1.0) + apca_presentation.v[merge_start_index + 1] * (apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index])) / double(apca_presentation.r[merge_start_index + 1] + 1);
		}

		apca_presentation.r[merge_start_index] = apca_presentation.r[merge_start_index + 1];

		merge_start_index++;
		//cout << "segment Number: " << apca_presentation.segmentNum << endl;
		//cout << "merge_start_index: " << merge_start_index << endl;
		while (merge_start_index < apca_presentation.segmentNum) {
			apca_presentation.r[merge_start_index] = apca_presentation.r[merge_start_index + 1];
			apca_presentation.v[merge_start_index] = apca_presentation.v[merge_start_index + 1];
			merge_start_index++;
		}
	}
}

TEMPLATE
void APCA_QUAL::getMergeIndexVector(const int& merge_frequency, vector<APCA>& apca_presentation) {//190501
	//void getMergeIndex(const DataType* original_time_series, const int& merge_frequency, APCA& apca_presentation) {
	DataType sefments_difference = NULL;
	int merge_start_index = NULL;
	double initial_segment_length = *apca_presentation[0].r + 1.0;

	for (int i = 0; i < merge_frequency; i++) { //merge times
		sefments_difference = DBL_MAX;

		//apca_presentation.segmentNum--;

		//cout << "segment Number: " << apca_presentation.segmentNum << endl;
		for (int j = 0; j < apca_presentation.size() - 1; j++) {
			//cout << "sefments_difference: " << sefments_difference << " sgetments diff: " << fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]) << endl;
			if (sefments_difference > fabs(*apca_presentation[j + 1].v - *apca_presentation[j].v)) {
				sefments_difference = fabs(*apca_presentation[j + 1].v - *apca_presentation[j].v);
				merge_start_index = j;
			}
		}
		//cout << "merge_start_index: " << merge_start_index << endl;
		//apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] + apca_presentation.v[merge_start_index + 1]) / 2;

		/*=============Merge average value===================*/
		if (merge_start_index != 0) {
			//double temp1 = getAve(original_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));
			//apca_presentation.v[merge_start_index] = getAve(original_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));
			*apca_presentation[merge_start_index].v = (*apca_presentation[merge_start_index].v * (*apca_presentation[merge_start_index].r - *apca_presentation[merge_start_index - 1].r) + *apca_presentation[merge_start_index + 1].v * (*apca_presentation[merge_start_index + 1].r - *apca_presentation[merge_start_index].r)) / double(*apca_presentation[merge_start_index + 1].r - *apca_presentation[merge_start_index - 1].r);
		}
		else {
			*apca_presentation[merge_start_index].v = (*apca_presentation[merge_start_index].v * double(*apca_presentation[merge_start_index].r + 1.0) + *apca_presentation[merge_start_index + 1].v * (*apca_presentation[merge_start_index + 1].r - *apca_presentation[merge_start_index].r)) / double(*apca_presentation[merge_start_index + 1].r + 1);
		}
		/*.....................................................*/

		/*=============Merge right end point===================*/
		*apca_presentation[merge_start_index].r = *apca_presentation[merge_start_index + 1].r;

		/*=============Erase trash segment===================*/
		apca_presentation.erase(apca_presentation.begin() + merge_start_index + 1);
		//merge_start_index++;
		////cout << "segment Number: " << apca_presentation.segmentNum << endl;
		////cout << "merge_start_index: " << merge_start_index << endl;
		//while (merge_start_index < segment_size) {
		//	apca_presentation.r[merge_start_index] = apca_presentation.r[merge_start_index + 1];
		//	apca_presentation.v[merge_start_index] = apca_presentation.v[merge_start_index + 1];
		//	merge_start_index++;
		//}
		/*.....................................................*/
	}
}

TEMPLATE
bool APCA_QUAL::mergeSegmentsIterator(APCA& apca_presentation, const DataType& n, const int& N, APCA& italicC) {
	//bool mergeSegmentsIterator(const DataType* original_time_series, APCA& apca_presentation, const DataType& n, const int& M, APCA& italicC) {
	assert(apca_presentation.segmentNum >= N);

	if (apca_presentation.segmentNum > N) {
		int merge_frequency = apca_presentation.segmentNum - N;

		//getMergeIndexQueue(merge_frequency, apca_presentation); //2 for 1
		getMergeIndex(merge_frequency, apca_presentation); //2 for 1

		assert(apca_presentation.segmentNum == italicC.segmentNum);
		copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);

		return true;
		//getExactMeanValue(original_time_series, apca_presentation);
	}
	else {
		copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);

		return true;
	}
}

TEMPLATE
bool APCA_QUAL::mergeSegmentsIteratorVector(vector<APCA>& const apca_presentation, const DataType& n, const int& N, vector<APCA>& const italicC) {//190501
	assert(apca_presentation.size() >= N);

	if (apca_presentation.size() > N) {
		int merge_frequency = apca_presentation.size() - N;

		//getMergeIndexQueue(merge_frequency, apca_presentation); //2 for 1
		getMergeIndexVector(merge_frequency, apca_presentation); //2 for 1

		assert(apca_presentation.size() == N);
		/*copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);*/

		for (auto&& au : apca_presentation) {
			italicC.push_back(au);
		}

		return true;
		//getExactMeanValue(original_time_series, apca_presentation);
	}
	else {
		/*copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);*/

		for (auto&& au : apca_presentation) {
			italicC.push_back(au);
		}

		return true;
	}
}

//************************************
// Method:mergeSegmentsIteratorVectorPreMemory
// Qualifier:190611 vector to instead Array, pre define memory
// date:190611
// author:
//************************************
TEMPLATE
template<typename T>
bool APCA_QUAL::mergeSegmentsIteratorLinkedList(vector<T>& const apca_presentation, const DataType& n, const int& N, DoublyLinkedList<T>& const italicC) {//191104 Use linked list ot instead vector
	assert(apca_presentation.size() >= N);

	if (apca_presentation.size() > N) {
		int merge_frequency = apca_presentation.size() - N;

		//getMergeIndexQueue(merge_frequency, apca_presentation); //2 for 1
		getMergeIndexVector(merge_frequency, apca_presentation); //2 for 1

		assert(apca_presentation.size() == N);
		/*copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);*/

		for (auto&& au : apca_presentation) {
			italicC.add(au);
		}

		return true;
		//getExactMeanValue(original_time_series, apca_presentation);
	}
	else {
		/*copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);*/

		for (auto&& au : apca_presentation) {
			italicC.add(au);
		}
		return true;
	}
}

//************************************
// Method:mergeSegmentsIteratorVectorPreMemory
// Qualifier:190611 vector to instead Array, pre define memory
// date:190611
// author:
//************************************
TEMPLATE
bool APCA_QUAL::mergeSegmentsIteratorVectorPreMemory(vector<APCA>& const apca_presentation, const DataType& n, const int& N, vector<APCA>& const italicC) {//190611
#ifdef _DEBUG
	assert(apca_presentation.size() >= N);
#endif

	if (apca_presentation.size() > N) {
		int merge_frequency = apca_presentation.size() - N;

		//getMergeIndexQueue(merge_frequency, apca_presentation); //2 for 1
		getMergeIndexVector(merge_frequency, apca_presentation); //2 for 1
#ifdef _DEBUG
		assert(apca_presentation.size() == N);
#endif
		/*copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);*/

		/*for (auto&& au : apca_presentation) {
			italicC.push_back(au);
		}*/

		copy_n(apca_presentation.begin(), N, italicC.begin());
		return true;
		//getExactMeanValue(original_time_series, apca_presentation);
	}
	else {
		/*copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);*/

		/*for (auto&& au : apca_presentation) {
			italicC.push_back(au);
		}*/
		copy_n(apca_presentation.begin(), N, italicC.begin());
		return true;
	}
}

TEMPLATE
void APCA_QUAL::getAPCAPoint(DataType*& const original_time_series, const double& n, const int& N, APCA& italicC) {
	assert(n > N);
	//TOOL::recordStartTime(TOOL::time_record[3]);

	int i = 0, j = 0, retained_coeffs_length = NULL;
	unsigned int length_power_of_2 = NULL;
	DataType* wavelet_transform_time_series;
	APCA general_apca_presentation;
	priority_queue<DataType> fq_truncate_index;

	length_power_of_2 = getNextPowerOf2(n);
	wavelet_transform_time_series = new DataType[length_power_of_2];
	//memory_account[3] = sizeof(wavelet_transform_time_series) * length_power_of_2;
	padZero(original_time_series, n, length_power_of_2, wavelet_transform_time_series);
	getHDWT(length_power_of_2, N, wavelet_transform_time_series, fq_truncate_index);

	retained_coeffs_length = getNextPowerOf2(fq_truncate_index.top() + 1);
	general_apca_presentation.segmentNum = retained_coeffs_length;
	general_apca_presentation.v = new DataType[retained_coeffs_length];
	general_apca_presentation.r = new DataType[retained_coeffs_length];

	//memory_account[3] += sizeof(general_apca_presentation.v) * retained_coeffs_length;
	//memory_account[3] += sizeof(general_apca_presentation.r) * retained_coeffs_length;

	reconstructApproximateAPCA(wavelet_transform_time_series, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	//cout << "Total Resolution: " << log2(getNextPowerOf2(n)) << ",  Resolution: " << log2(retained_coeffs_length) << endl;
	truncateZero(n, length_power_of_2, general_apca_presentation);
	getExactMeanValue(original_time_series, general_apca_presentation);
	mergeSegmentsIterator(general_apca_presentation, n, N, italicC);
	//mergeSegmentsRecursivly0(original_time_series, general_apca_presentation, n, N, italicC);

	priority_queue<DataType>().swap(fq_truncate_index);
	delete[] general_apca_presentation.r;
	general_apca_presentation.r = nullptr;
	delete[] general_apca_presentation.v;
	general_apca_presentation.v = nullptr;
	delete[] wavelet_transform_time_series;
	wavelet_transform_time_series = nullptr;

	//output_argument.run_time = TOOL::recordFinishTime(TOOL::time_record[3]);
	//cout << "APCA running Time: " << output_argument.run_time << endl;// compare percentage time
}

//************************************
// Method:getAPCAPointVector
// Qualifier:190430 vector to instead Array
// date:190430
// author:
//************************************
TEMPLATE
void APCA_QUAL::getAPCAPointVector(DataType*& const original_time_series, const double& n, const int& N, vector<APCA>& const italicC) {//190430 vector to instead Array
	assert(n > N);

	//TOOL::recordStartTime(TOOL::time_record[3]);

	int i = 0, j = 0, retained_coeffs_length = NULL;
	unsigned int length_power_of_2 = NULL;
	//DataType* wavelet_transform_time_series;
	vector<DataType> wavelet_transform_vector;//190430
	vector<APCA> general_apca_presentation;
	priority_queue<DataType> fq_truncate_index;

	length_power_of_2 = getNextPowerOf2(n);
	//wavelet_transform_time_series = new DataType[length_power_of_2];
	//memory_account[3] = sizeof(wavelet_transform_time_series)*length_power_of_2;
	//padZero(original_time_series, n, length_power_of_2, wavelet_transform_time_series);
	padZeroVector(original_time_series, n, length_power_of_2, wavelet_transform_vector);
	//getHDWT(length_power_of_2, N, wavelet_transform_time_series, fq_truncate_index);
	getHDWTVector(length_power_of_2, N, wavelet_transform_vector, fq_truncate_index);

	retained_coeffs_length = getNextPowerOf2(fq_truncate_index.top() + 1);
	//general_apca_presentation.segmentNum = retained_coeffs_length;
	//general_apca_presentation.v = new DataType[retained_coeffs_length];
	//general_apca_presentation.r = new DataType[retained_coeffs_length];

	//memory_account[3] += sizeof(general_apca_presentation.v)*retained_coeffs_length;
	//memory_account[3] += sizeof(general_apca_presentation.r)*retained_coeffs_length;

	//reconstructApproximateAPCA(wavelet_transform_time_series, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	reconstructApproximateAPCAVector(wavelet_transform_vector, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	assert(general_apca_presentation.size() == retained_coeffs_length);
	//cout << "Total Resolution: " << log2(getNextPowerOf2(n)) << ",  Resolution: " << log2(retained_coeffs_length) << endl;
	truncateZeroVector(n, length_power_of_2, general_apca_presentation);
	getExactMeanValueVector(original_time_series, general_apca_presentation);
	mergeSegmentsIteratorVector(general_apca_presentation, n, N, italicC);
	//mergeSegmentsRecursivly0(original_time_series, general_apca_presentation, n, N, italicC);

	/*cout << "****:" << endl;
	for (auto au : italicC) {
		cout << *au.v << ", ";
	}
	cout << endl;*/
	priority_queue<DataType>().swap(fq_truncate_index);
	//delete[] general_apca_presentation.r;
	//general_apca_presentation.r = nullptr;
	//delete[] general_apca_presentation.v;
	//general_apca_presentation.v = nullptr;
	//delete[] wavelet_transform_time_series;
	//wavelet_transform_time_series = nullptr;
	wavelet_transform_vector.clear();
	wavelet_transform_vector.shrink_to_fit();
	general_apca_presentation.clear();
	general_apca_presentation.shrink_to_fit();

	//output_argument.run_time = TOOL::recordFinishTime(TOOL::time_record[3]);
	//cout << "APCA running Time: " << output_argument.run_time << endl;// compare percentage time
}

//************************************
// Method:getAPCAPointLinkedList
// Qualifier: Use Linked list to instead vector
// date:191104
// author:
//************************************
TEMPLATE
template<typename T>
void APCA_QUAL::getAPCAPointLinkedList(DataType*& const original_time_series, const double& n, const int& N, DoublyLinkedList<T>& const italicC) {//191104 Linked list to instead vector
#ifdef _DEBUG
	assert(n > N);
#endif
	//TOOL::recordStartTime(TOOL::time_record[3]);

	int i = 0, j = 0, retained_coeffs_length = NULL;
	unsigned int length_power_of_2 = NULL;
	//DataType* wavelet_transform_time_series;
	vector<DataType> wavelet_transform_vector;//190430
	vector<APCA> general_apca_presentation;
	priority_queue<DataType> fq_truncate_index;

	length_power_of_2 = getNextPowerOf2(n);
	//wavelet_transform_time_series = new DataType[length_power_of_2];
	//memory_account[3] = sizeof(wavelet_transform_time_series)*length_power_of_2;
	//padZero(original_time_series, n, length_power_of_2, wavelet_transform_time_series);
	padZeroVector(original_time_series, n, length_power_of_2, wavelet_transform_vector);
	//getHDWT(length_power_of_2, N, wavelet_transform_time_series, fq_truncate_index);
	getHDWTVector(length_power_of_2, N, wavelet_transform_vector, fq_truncate_index);

	retained_coeffs_length = getNextPowerOf2(fq_truncate_index.top() + 1);
	//general_apca_presentation.segmentNum = retained_coeffs_length;
	//general_apca_presentation.v = new DataType[retained_coeffs_length];
	//general_apca_presentation.r = new DataType[retained_coeffs_length];

	//memory_account[3] += sizeof(general_apca_presentation.v)*retained_coeffs_length;
	//memory_account[3] += sizeof(general_apca_presentation.r)*retained_coeffs_length;

	//reconstructApproximateAPCA(wavelet_transform_time_series, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	reconstructApproximateAPCAVector(wavelet_transform_vector, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	assert(general_apca_presentation.size() == retained_coeffs_length);
	//cout << "Total Resolution: " << log2(getNextPowerOf2(n)) << ",  Resolution: " << log2(retained_coeffs_length) << endl;
	truncateZeroVector(n, length_power_of_2, general_apca_presentation);
	getExactMeanValueVector(original_time_series, general_apca_presentation);
	//mergeSegmentsIteratorVector(general_apca_presentation, n, N, italicC);
	mergeSegmentsIteratorLinkedList(general_apca_presentation, n, N, italicC);//191104 use linked list to instead vector
	//mergeSegmentsRecursivly0(original_time_series, general_apca_presentation, n, N, italicC);

	/*cout << "****:" << endl;
	for (auto au : italicC) {
		cout << *au.v << ", ";
	}
	cout << endl;*/
	priority_queue<DataType>().swap(fq_truncate_index);
	//delete[] general_apca_presentation.r;
	//general_apca_presentation.r = nullptr;
	//delete[] general_apca_presentation.v;
	//general_apca_presentation.v = nullptr;
	//delete[] wavelet_transform_time_series;
	//wavelet_transform_time_series = nullptr;
	wavelet_transform_vector.clear();
	wavelet_transform_vector.shrink_to_fit();
	general_apca_presentation.clear();
	general_apca_presentation.shrink_to_fit();

	//output_argument.run_time = TOOL::recordFinishTime(TOOL::time_record[3]);
	//cout << "APCA running Time: " << output_argument.run_time << endl;// compare percentage time

}

//************************************
// Method:getAPCAPointVectorPreMemory
// Qualifier:190611 vector to instead Array, Pre-define memory
// date:190611
// author:
//************************************
TEMPLATE
void APCA_QUAL::getAPCAPointVectorPreMemory(DataType*& const original_time_series, const double& n, const int& N, vector<APCA>& const italicC) {//190611 Pre-define memory of vector
#ifdef _DEBUG
	assert(n > N);
#endif
	//TOOL::recordStartTime(TOOL::time_record[3]);
	int i = 0, j = 0, retained_coeffs_length = NULL;
	unsigned int length_power_of_2 = NULL;
	//DataType* wavelet_transform_time_series;
	vector<DataType> wavelet_transform_vector;//190430
	vector<APCA> general_apca_presentation;
	priority_queue<DataType> fq_truncate_index;

	length_power_of_2 = getNextPowerOf2(n);//powe of 2

	/*======================pre-define memory===========================================*/
	wavelet_transform_vector.resize(length_power_of_2, DataType());//190611
	italicC.resize(N, APCA());//190611
	/*.................................................................*/

	//wavelet_transform_time_series = new DataType[length_power_of_2];
	//memory_account[3] = sizeof(wavelet_transform_time_series)*length_power_of_2;
	//padZero(original_time_series, n, length_power_of_2, wavelet_transform_time_series);
	//padZeroVector(original_time_series, n, length_power_of_2, wavelet_transform_vector);
	padZeroVectorPreMemory(original_time_series, n, length_power_of_2, wavelet_transform_vector);//190611
	//getHDWT(length_power_of_2, N, wavelet_transform_time_series, fq_truncate_index);
	getHDWTVector(length_power_of_2, N, wavelet_transform_vector, fq_truncate_index);

	retained_coeffs_length = getNextPowerOf2(fq_truncate_index.top() + 1);
	//general_apca_presentation.segmentNum = retained_coeffs_length;
	//general_apca_presentation.v = new DataType[retained_coeffs_length];
	//general_apca_presentation.r = new DataType[retained_coeffs_length];

	//memory_account[3] += sizeof(general_apca_presentation.v)*retained_coeffs_length;
	//memory_account[3] += sizeof(general_apca_presentation.r)*retained_coeffs_length;

	//reconstructApproximateAPCA(wavelet_transform_time_series, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	//reconstructApproximateAPCAVector(wavelet_transform_vector, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	reconstructApproximateAPCAVectorPreMemory(wavelet_transform_vector, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);//190611
#ifdef _DEBUG
	assert(general_apca_presentation.size() == retained_coeffs_length);
#endif
	//cout << "Total Resolution: " << log2(getNextPowerOf2(n)) << ",  Resolution: " << log2(retained_coeffs_length) << endl;
	truncateZeroVector(n, length_power_of_2, general_apca_presentation);
	getExactMeanValueVector(original_time_series, general_apca_presentation);
	//mergeSegmentsIteratorVector(general_apca_presentation, n, N, italicC);
	mergeSegmentsIteratorVectorPreMemory(general_apca_presentation, n, N, italicC);
	//mergeSegmentsRecursivly0(original_time_series, general_apca_presentation, n, N, italicC);

	/*cout << "****:" << endl;
	for (auto au : italicC) {
		cout << *au.v << ", ";
	}
	cout << endl;*/
	priority_queue<DataType>().swap(fq_truncate_index);
	//delete[] general_apca_presentation.r;
	//general_apca_presentation.r = nullptr;
	//delete[] general_apca_presentation.v;
	//general_apca_presentation.v = nullptr;
	//delete[] wavelet_transform_time_series;
	//wavelet_transform_time_series = nullptr;
	wavelet_transform_vector.clear();
	wavelet_transform_vector.shrink_to_fit();
	general_apca_presentation.clear();
	general_apca_presentation.shrink_to_fit();

	//output_argument.run_time = TOOL::recordFinishTime(TOOL::time_record[3]);
	//cout << "APCA running Time: " << output_argument.run_time << endl;// compare percentage time
}

TEMPLATE
void  APCA_QUAL::divideRemainderPAA(const DataType* original_time_series, const APCA& italicC, const DataType& n, const int& N) {
	//printf("APCA divideRemainderAPCA()\n");

	//double *y = new double[N];

	//cout << "APCA.segmentNum = " << italicC->segmentNum << endl;
	assert(n >= N && italicC.segmentNum == N);
	//cout << n << endl;
	//cout << N << endl;
	double Remainder = int(n) % int(N);
	//cout << "Remainder = " << Remainder << ", n = " << n << endl;
	double integerDividend = n - Remainder;
	double integerQuotient = integerDividend / N;
	double integerSegmentLength = integerQuotient + 1;

	int j = 0;
	int endOfLongSegment = 0, indexOfLongSegment = 0;

	for (int i = 1; i <= N; i++) {
		//cout << endl;
		double sum = 0;
		if (i <= Remainder) {
			for (j = int((i - 1) * integerSegmentLength); j < i * integerSegmentLength; j++) {
				sum += original_time_series[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
			}
			//y[i - 1] = sum / integerSegmentLength;
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			italicC.v[i - 1] = sum / integerSegmentLength;
			italicC.r[i - 1] = j - 1;
			//cout << endl << "v[" << i - 1 << "] = " << italicC.v[i - 1] << ", r[" << i - 1 << "] = " << italicC.r[i - 1] << endl;
			indexOfLongSegment = i;
			endOfLongSegment = int(i * integerSegmentLength);
		}
		else {
			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * integerQuotient); j < endOfLongSegment + (i - indexOfLongSegment) * integerQuotient; j++) {
				sum += original_time_series[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
				//cout << original_time_series[j] << " ";
			}
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			italicC.v[i - 1] = sum / integerQuotient;
			italicC.r[i - 1] = j - 1;
			//cout << endl << "v[" << i - 1 << "] = " << italicC.v[i - 1] << ", r[" << i - 1 << "] = " << italicC.r[i - 1] << endl;
		}
	}
}

//200923 use vector instead pointer
TEMPLATE
template<typename T, typename Y, typename U>
void APCA_QUAL::get_PAA(const vector<T>& const original_time_series_vector, const Y& const segment_number, U& const paa_container) {
	assert(original_time_series_vector.size() >= segment_number);
	//cout << n << endl;
	//cout << N << endl;
	const auto& const n = original_time_series_vector.size();
	const auto& const N = segment_number;

	double Remainder = int(original_time_series_vector.size()) % int(segment_number);
	//cout << "Remainder = " << Remainder << ", n = " << n << endl;
	double integerDividend = n - Remainder;
	double integerQuotient = integerDividend / N;
	double integerSegmentLength = integerQuotient + 1;

	int j = 0;
	int endOfLongSegment = 0, indexOfLongSegment = 0;

	for (int i = 1; i <= N; i++) {
		//cout << endl;
		double sum = 0;
		if (i <= Remainder) {
			for (j = int((i - 1) * integerSegmentLength); j < i * integerSegmentLength; j++) {
				sum += original_time_series_vector[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
			}
			//y[i - 1] = sum / integerSegmentLength;
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			paa_container.v[i - 1] = sum / integerSegmentLength;
			paa_container.r[i - 1] = j - 1;
			//cout << endl << "v[" << i - 1 << "] = " << italicC.v[i - 1] << ", r[" << i - 1 << "] = " << italicC.r[i - 1] << endl;
			indexOfLongSegment = i;
			endOfLongSegment = int(i * integerSegmentLength);
		}
		else {
			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * integerQuotient); j < endOfLongSegment + (i - indexOfLongSegment) * integerQuotient; j++) {
				sum += original_time_series_vector[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
				//cout << original_time_series[j] << " ";
			}
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			paa_container.v[i - 1] = sum / integerQuotient;
			paa_container.r[i - 1] = j - 1;
			//cout << endl << "v[" << i - 1 << "] = " << italicC.v[i - 1] << ", r[" << i - 1 << "] = " << italicC.r[i - 1] << endl;
		}
	}
}

//************************************
// Method:getPAAPointVector
// Qualifier:Use vector instead Array
// Input:
// Output:
// date:190515
// author:
//************************************
TEMPLATE
void  APCA_QUAL::getPAAPointVector(const DataType* original_time_series, const DataType& n, const int& N, vector<APCA>& const italicC) {//190515 Only use vecter. Same with divideRemainderPAA()
//printf("APCA divideRemainderAPCA()\n");
	//double *y = new double[N];
	//cout << "APCA.segmentNum = " << italicC->segmentNum << endl;
	assert(n >= N);
	//cout << n << endl;
	//cout << N << endl;
	APCA temp_paa;//190515

	double Remainder = int(n) % int(N);
	//cout << "Remainder = " << Remainder << ", n = " << n << endl;
	double integerDividend = n - Remainder;
	double integerQuotient = integerDividend / N;
	double integerSegmentLength = integerQuotient + 1;

	int j = 0;
	int endOfLongSegment = 0, indexOfLongSegment = 0;

	for (int i = 1; i <= N; i++) {
		//cout << endl;
		double sum = 0;
		temp_paa.v = new DataType;
		temp_paa.r = new DataType;
		if (i <= Remainder) {
			for (j = int((i - 1) * integerSegmentLength); j < i * integerSegmentLength; j++) {
				sum += original_time_series[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
			}
			//y[i - 1] = sum / integerSegmentLength;
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			//italicC.v[i - 1] = sum / integerSegmentLength;
			//italicC.r[i - 1] = j - 1;
			*temp_paa.v = DataType(sum / integerSegmentLength);
			*temp_paa.r = DataType(j - 1);
			//cout << "**** " << *temp_paa.v << " " << *temp_paa.r << endl;
			italicC.emplace_back(temp_paa);
			//cout << endl << "v[" << i - 1 << "] = " << italicC.v[i - 1] << ", r[" << i - 1 << "] = " << italicC.r[i - 1] << endl;
			indexOfLongSegment = i;
			endOfLongSegment = int(i * integerSegmentLength);
		}
		else {
			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * integerQuotient); j < endOfLongSegment + (i - indexOfLongSegment) * integerQuotient; j++) {
				sum += original_time_series[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
				//cout << original_time_series[j] << " ";
			}
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			//italicC.v[i - 1] = sum / integerQuotient;
			//italicC.r[i - 1] = j - 1;
			*temp_paa.v = DataType(sum / integerQuotient);
			*temp_paa.r = DataType(j - 1);
			//cout <<"**** "<< *temp_paa.v <<" "<< *temp_paa.r << endl;
			italicC.emplace_back(temp_paa);
			//cout << endl << "v[" << i - 1 << "] = " << *italicC[i - 1].v << ", r[" << i - 1 << "] = " <<*italicC[i - 1].r << endl;
		}
	}
	assert(italicC.size() == N);
	/*cout << "Test: " << endl;
	for (auto au : italicC) {
		cout << *au.v << endl;
	}*/

	/*delete temp_paa.v;
	temp_paa.v = nullptr;
	delete temp_paa.r;
	temp_paa.r = nullptr;*/
}


//************************************
// Method:getPAAPointLinkedList
// Qualifier:use Linked list to instead vector. Same with divideRemainderPAA()
// Input:
// Output:
// date:191104
// author:
//************************************
TEMPLATE
template<typename T>
void  APCA_QUAL::getPAAPointLinkedList(const DataType* original_time_series, const DataType& n, const int& N, DoublyLinkedList<T>& const italicC) {//
	//printf("APCA divideRemainderAPCA()\n");
	//double *y = new double[N];
	//cout << "APCA.segmentNum = " << italicC->segmentNum << endl;
	assert(n >= N);
	//cout << n << endl;
	//cout << N << endl;
	APCA temp_paa;//190515

	double Remainder = int(n) % int(N);
	//cout << "Remainder = " << Remainder << ", n = " << n << endl;
	double integerDividend = n - Remainder;
	double integerQuotient = integerDividend / N;
	double integerSegmentLength = integerQuotient + 1;

	int j = 0;
	int endOfLongSegment = 0, indexOfLongSegment = 0;

	for (int i = 1; i <= N; i++) {
		//cout << endl;
		double sum = 0;
		temp_paa.v = new DataType;
		temp_paa.r = new DataType;
		if (i <= Remainder) {
			for (j = int((i - 1) * integerSegmentLength); j < i * integerSegmentLength; j++) {
				sum += original_time_series[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
			}
			//y[i - 1] = sum / integerSegmentLength;
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			//italicC.v[i - 1] = sum / integerSegmentLength;
			//italicC.r[i - 1] = j - 1;
			*temp_paa.v = DataType(sum / integerSegmentLength);
			*temp_paa.r = DataType(j - 1);
			//cout << "**** " << *temp_paa.v << " " << *temp_paa.r << endl;
			italicC.emplace_back(temp_paa);
			//cout << endl << "v[" << i - 1 << "] = " << italicC.v[i - 1] << ", r[" << i - 1 << "] = " << italicC.r[i - 1] << endl;
			indexOfLongSegment = i;
			endOfLongSegment = int(i * integerSegmentLength);
		}
		else {
			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * integerQuotient); j < endOfLongSegment + (i - indexOfLongSegment) * integerQuotient; j++) {
				sum += original_time_series[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
				//cout << original_time_series[j] << " ";
			}
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			//italicC.v[i - 1] = sum / integerQuotient;
			//italicC.r[i - 1] = j - 1;
			*temp_paa.v = DataType(sum / integerQuotient);
			*temp_paa.r = DataType(j - 1);
			//cout <<"**** "<< *temp_paa.v <<" "<< *temp_paa.r << endl;
			italicC.emplace_back(temp_paa);
			//cout << endl << "v[" << i - 1 << "] = " << *italicC[i - 1].v << ", r[" << i - 1 << "] = " <<*italicC[i - 1].r << endl;
		}
	}
	assert(italicC.size() == N);
	/*cout << "Test: " << endl;
	for (auto au : italicC) {
		cout << *au.v << endl;
	}*/

	/*delete temp_paa.v;
	temp_paa.v = nullptr;
	delete temp_paa.r;
	temp_paa.r = nullptr;*/
}

//************************************
// Method:getPAAPointVectorPreMemory
// Qualifier:Use vector instead Array, Pre-define memory
// Input:
// Output:
// date:190611
// author:
//************************************
TEMPLATE//approximate original function
void APCA_QUAL::getPAAPointVectorPreMemory(const DataType* original_time_series, const DataType& n, const int& N, vector<APCA>& const italicC) {//190611 Only use vecter. Same with divideRemainderPAA(), pre-define memory
	//printf("APCA divideRemainderAPCA()\n");

	//double *y = new double[N];

	//cout << "APCA.segmentNum = " << italicC->segmentNum << endl;

#ifdef _DEBUG
	assert(n >= N);
#endif
	//cout << n << endl;
	//cout << N << endl;
	APCA temp_paa;//190515

	/*=========================Pre-define memory================================*/
	italicC.resize(N, APCA());//190611
	for (auto&& au : italicC) {
		au.v = new DataType;
		au.r = new DataType;
	}

	/*..........................................................................*/

	double Remainder = int(n) % int(N);
	//cout << "Remainder = " << Remainder << ", n = " << n << endl;
	double integerDividend = n - Remainder;
	double integerQuotient = integerDividend / N;
	double integerSegmentLength = integerQuotient + 1;

	int j = 0;
	int endOfLongSegment = 0, indexOfLongSegment = 0;

	for (int i = 1; i <= N; i++) {
		//cout << endl;
		double sum = 0;
		temp_paa.v = new DataType;
		temp_paa.r = new DataType;
		if (i <= Remainder) {
			for (j = int((i - 1) * integerSegmentLength); j < i * integerSegmentLength; j++) {
				sum += original_time_series[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
			}
			//y[i - 1] = sum / integerSegmentLength;
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			//italicC.v[i - 1] = sum / integerSegmentLength;
			//italicC.r[i - 1] = j - 1;
			*temp_paa.v = DataType(sum / integerSegmentLength);
			*temp_paa.r = DataType(j - 1);
			//cout << "**** " << *temp_paa.v << " " << *temp_paa.r << endl;
			//italicC.emplace_back(temp_paa);
			*italicC[i - 1].v = DataType(sum / integerSegmentLength);
			*italicC[i - 1].r = DataType(j - 1);
			//cout << endl << "v[" << i - 1 << "] = " << italicC.v[i - 1] << ", r[" << i - 1 << "] = " << italicC.r[i - 1] << endl;
			indexOfLongSegment = i;
			endOfLongSegment = int(i * integerSegmentLength);
		}
		else {
			for (j = int(endOfLongSegment + (i - indexOfLongSegment - 1) * integerQuotient); j < endOfLongSegment + (i - indexOfLongSegment) * integerQuotient; j++) {
				sum += original_time_series[j];
				//cout << "x[" << j << "] = " << original_time_series[j] << " ";
				//cout << original_time_series[j] << " ";
			}
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			//italicC.v[i - 1] = sum / integerQuotient;
			//italicC.r[i - 1] = j - 1;
			*temp_paa.v = DataType(sum / integerQuotient);
			*temp_paa.r = DataType(j - 1);
			//cout <<"**** "<< *temp_paa.v <<" "<< *temp_paa.r << endl;
			//italicC.emplace_back(temp_paa);
			*italicC[i - 1].v = DataType(sum / integerQuotient);
			*italicC[i - 1].r = DataType(j - 1);
			//cout << endl << "v[" << i - 1 << "] = " << *italicC[i - 1].v << ", r[" << i - 1 << "] = " <<*italicC[i - 1].r << endl;
		}
	}
#ifdef _DEBUG
	assert(italicC.size() == N);
#endif

	/*cout << "Test: " << endl;
	for (auto au : italicC) {
		cout << *au.v << endl;
	}*/

	/*delete temp_paa.v;
	temp_paa.v = nullptr;
	delete temp_paa.r;
	temp_paa.r = nullptr;*/
}

//¡¶A New Pattern Representation Method for Time - series Data¡·
//************************************
// Method:getPAAPointVectorPreMemory
// Qualifier:Use vector instead Array, Pre-define memory
// Input:
// Output:
// date:190611
// author:
//************************************
TEMPLATE
template<typename T>//190515 use Linked list to instead vector. Same with divideRemainderPAA()
void APCA_QUAL::get_paa_lagrangian(T& const italicC) {
	assert(italicC.segmentNum != INF);
	double Lagrange_multiplier = INF;
	double square_sum = 0;
	for (int segment_id = 0; segment_id < italicC.segmentNum; segment_id++) {
		assert(italicC.v[segment_id] != INF && italicC.r[segment_id] != INF);
		square_sum += italicC.v[segment_id] * italicC.v[segment_id];
	}

	Lagrange_multiplier = sqrt(square_sum / 4.0);
	Lagrange_multiplier = 1 / (Lagrange_multiplier * 2.0);

	assert(Lagrange_multiplier != 0);

	for (int segment_id = 0; segment_id < italicC.segmentNum; segment_id++) {
		italicC.v[segment_id] *= Lagrange_multiplier;
		assert(italicC.v[segment_id] != INF && italicC.r[segment_id] != INF);
	}

}

TEMPLATE//approximate original function
void APCA_QUAL::approximateOriginalFunction(const APCA_QUAL::APCA& italicC, DataType*& const approximated_time_series) {
	//double sum = 0;
	//double difference = NULL;
	//deviation_max = 0;

	for (int interval_id = 0; interval_id <= italicC.r[0]; interval_id++) {
		approximated_time_series[interval_id] = italicC.v[0];
		//cout <<"difference: "<< difference << endl;
		//deviation_max = difference > deviation_max ? difference : deviation_max;
		//sum += difference * difference;
	}

	for (int segment_id = 1; segment_id < italicC.segmentNum; segment_id++) {
		//cout<<"segment_id: "<< segment_id <<endl;
		for (int interval_id = 0; interval_id < italicC.r[segment_id] - italicC.r[segment_id - 1]; interval_id++) {
			approximated_time_series[int(interval_id + italicC.r[segment_id - 1] + 1)] = italicC.v[segment_id];
			/*cout << "difference: " << difference << endl;*/
			//deviation_max = difference > deviation_max ? difference : deviation_max;
			//sum += difference * difference;
		}
	}
}

TEMPLATE
void  APCA_QUAL::initialAPCAArray(const int& mg_d_index_point_number, const int& mg_APCA_point_dimension, APCA_ARRAY*& const APCA_array) {

	APCA_array = new APCA_ARRAY[mg_d_index_point_number];

	for (int i = 0; i < int(mg_d_index_point_number); i++) {
		APCA_array[i].APCALink.segmentNum = mg_APCA_point_dimension;
		APCA_array[i].APCALink.r = new DataType[mg_APCA_point_dimension];
		APCA_array[i].APCALink.v = new DataType[mg_APCA_point_dimension];
		//fill(APCA_array[i].APCALink.r, APCA_array[i].APCALink.r + mg_APCA_point_dimension, NULL);
		//fill(APCA_array[i].APCALink.v, APCA_array[i].APCALink.v + mg_APCA_point_dimension, NULL);
		fill_n(APCA_array[i].APCALink.r, mg_APCA_point_dimension, INF);
		fill_n(APCA_array[i].APCALink.v, mg_APCA_point_dimension, INF);
	}
}

TEMPLATE
void  APCA_QUAL::initialAPCAArray(const int& point_number, const int& point_dimension, APCA*& const APCA_array) {
	APCA_array = new APCA[point_number];

	for (int i = 0; i < int(point_number); i++) {
		APCA_array[i].segmentNum = point_dimension;
		APCA_array[i].r = new DataType[point_dimension];
		APCA_array[i].v = new DataType[point_dimension];
		//fill(APCA_array[i].APCALink.r, APCA_array[i].APCALink.r + mg_APCA_point_dimension, NULL);
		//fill(APCA_array[i].APCALink.v, APCA_array[i].APCALink.v + mg_APCA_point_dimension, NULL);
		fill_n(APCA_array[i].r, point_dimension, INF);
		fill_n(APCA_array[i].v, point_dimension, INF);
	}
}

TEMPLATE
void  APCA_QUAL::initialAPCAArray(const typename TOOL::INPUT_ARGUMENT& input_argument, APCA_ARRAY*& const APCA_array) {
	APCA_array = new APCA_ARRAY[input_argument.point_number];

	for (int i = 0; i < int(input_argument.point_number); i++) {
		APCA_array[i].APCALink.segmentNum = input_argument.point_dimension;
		APCA_array[i].APCALink.r = new DataType[input_argument.point_dimension];
		APCA_array[i].APCALink.v = new DataType[input_argument.point_dimension];
		//fill(APCA_array[i].APCALink.r, APCA_array[i].APCALink.r + mg_APCA_point_dimension, NULL);
		//fill(APCA_array[i].APCALink.v, APCA_array[i].APCALink.v + mg_APCA_point_dimension, NULL);
		fill_n(APCA_array[i].APCALink.r, input_argument.point_dimension, INF);
		fill_n(APCA_array[i].APCALink.v, input_argument.point_dimension, INF);
	}
}

//***************************************************************
// Method:initial_rect_vector
// Qualifier:
// Input:
// Output:
// date:191118
// author:
//***************************************************************
TEMPLATE
void APCA_QUAL::initial_rect_vector(const int& const point_number, const int& const point_dimension, vector<RTREE::Rect>& const rtree_rectangle_vector) {
#ifdef _DEBUG
	assert(point_number != INF && point_number > 0 && point_dimension != INF && point_dimension > 0);
#endif
	rtree_rectangle_vector.resize(point_number);
	for (auto&& au : rtree_rectangle_vector) {
		au.m_min = new DataType[point_dimension];
		au.m_max = new DataType[point_dimension];
	}
}

TEMPLATE
void  APCA_QUAL::deleteAPCAArray(const DataType& mg_d_index_point_number, APCA_ARRAY*& APCA_array) {
	for (int i = 0; i < mg_d_index_point_number; i++) {
		delete[] APCA_array[i].APCALink.r;
		APCA_array[i].APCALink.r = nullptr;

		delete[] APCA_array[i].APCALink.v;
		APCA_array[i].APCALink.v = nullptr;

		APCA_array[i].APCALink.segmentNum = NULL;
	}
	delete[] APCA_array;
	APCA_array = nullptr;
}

TEMPLATE
void  APCA_QUAL::deleteAPCAArray(const DataType& g_index_point_number, APCA*& const APCA_array) {
	if (APCA_array) {
		for (int i = 0; i < g_index_point_number; i++) {
			if (APCA_array[i].r != nullptr) {
				delete[] APCA_array[i].r;
				APCA_array[i].r = nullptr;
			}
			if (APCA_array[i].v != nullptr) {
				delete[] APCA_array[i].v;
				APCA_array[i].v = nullptr;
			}
			APCA_array[i].segmentNum = NULL;
		}
		delete[] APCA_array;
		APCA_array = nullptr;
	}
}

TEMPLATE
void  APCA_QUAL::initialAPCA(APCA& apca, int N) {
	apca.segmentNum = N;
	apca.r = new DataType[N];
	apca.v = new DataType[N];
}

TEMPLATE
void APCA_QUAL::deleteAPCA(APCA& apca) {
	apca.segmentNum = NULL;

	delete[] apca.r;
	apca.r = nullptr;

	delete[] apca.v;
	apca.v = nullptr;
}

//***************************************************************
// Method:delete_rect_vector
// Qualifier: clear rectangle vector
// Input:
// Output: 
// Note:
// date: 191118 22:15
// author:
//***************************************************************
TEMPLATE
void APCA_QUAL::delete_rect_vector(vector<RTREE::Rect>& const rtree_rectangle_vector) {
	for (auto au : rtree_rectangle_vector) {
		if (au.m_min != nullptr) {
			delete[] au.m_min;
			au.m_min = nullptr;
		}

		if (au.m_max != nullptr) {
			delete[] au.m_max;
			au.m_max = nullptr;
		}
	}
	rtree_rectangle_vector.clear();
}

TEMPLATE
typename APCA_QUAL::APCA& APCA_QUAL::getCmin(const double* C, const APCA& italicC, APCA& Cmin) {
	//cout << "getCmin() \n";

	Cmin.segmentNum = italicC.segmentNum;
	//cout << "segmentNum = " << Cmin.segmentNum << endl;

	int i = 0;
	Cmin.r[i] = italicC.r[i];
	Cmin.v[i] = *min_element(C, C + int(italicC.r[i]) + 1);

	//cout << "r[" << i << "] = " << Cmin->r[i] << ", Cmin0[" << i << "] = " << Cmin->v[i] << endl << endl;

	for (i = 1; i < italicC.segmentNum; i++) {
		Cmin.r[i] = italicC.r[i];
		Cmin.v[i] = *min_element(C + int(italicC.r[i - 1]) + 1, C + int(italicC.r[i]) + 1);//min[begin, end);

		//cout << "r[" << i << "] = " << Cmin.r[i] << ", Cmin0[" << i << "] = " << Cmin.v[i] << endl << endl;
	}
	return Cmin;
}

TEMPLATE
typename APCA_QUAL::APCA& APCA_QUAL::getCmax(const double* C, const APCA& italicC, APCA& Cmax) {
	//cout << "getCmax() \n";
	Cmax.segmentNum = italicC.segmentNum;
	//cout << "segmentNum = " << Cmax.segmentNum << endl;
	int i = 0;
	Cmax.r[i] = italicC.r[i];
	Cmax.v[i] = *max_element(C, C + int(italicC.r[i]) + 1);

	//cout << "r[" << i << "] = " << Cmax.r[i] << ", Cmax0[" << i << "] = " << Cmax.v[i] << endl << endl;

	for (i = 1; i < italicC.segmentNum; i++) {
		Cmax.r[i] = italicC.r[i];
		Cmax.v[i] = *max_element(C + int(italicC.r[i - 1]) + 1, C + int(italicC.r[i]) + 1);//max[begin, end);

		//cout << "r[" << i << "] = " << Cmax.r[i] << ", Cmax0[" << i << "] = " << Cmax.v[i] << endl << endl;
	}
	return Cmax;
}

//***************************************************************
// Method:get_minmax_original
// Qualifier: get min max point of every segment, Use original MBR, not APCA paper
// Input:
// Output: Approximated time series
// Note:
// date: 191118 22:15
// author:
//***************************************************************
TEMPLATE
template<typename T, typename Y>
Y& APCA_QUAL::get_minmax_original(const vector<DataType>& const time_series_vector, const int& const dimension_id, const T& const italicC, Y& const rtree_rectangle) {
#ifdef _DEBUG
	assert(!time_series_vector.empty() && italicC.segmentNum != INF && time_series_vector.size() == italicC.r[italicC.segmentNum - 1] + 1);
#endif
	int rectangle_length = italicC.segmentNum << 1;

	//odd is value, even is id
	//rtree_rectangle.m_min = new DataType[rectangle_length];
	//rtree_rectangle.m_max = new DataType[rectangle_length];

	const auto min_max = std::minmax_element(time_series_vector.begin(), time_series_vector.begin() + italicC.r[0] + 1);
	rtree_rectangle.m_min[0] = 0; // even min id
	rtree_rectangle.m_min[1] = *min_max.first;// odd first is min value
	rtree_rectangle.m_max[0] = italicC.r[0]; // even max id
	rtree_rectangle.m_max[1] = *min_max.second;// odd second is max value
#ifdef _DEBUG
	//cout << rtree_rectangle.m_min[0] <<" "<< rtree_rectangle.m_min[1] <<" "<< rtree_rectangle.m_max[0] <<" "<< rtree_rectangle.m_max[1] << endl;
	assert(rtree_rectangle.m_min[0] <= rtree_rectangle.m_max[0]);
#endif

	for (int segment_id = 1; segment_id < italicC.segmentNum; segment_id++) {
#ifdef _DEBUG
		assert(italicC.r[segment_id] != INF && italicC.v[segment_id] != INF);
#endif
		const auto min_max = std::minmax_element(time_series_vector.begin() + italicC.r[segment_id - 1] + 1, time_series_vector.begin() + italicC.r[segment_id] + 1);
		rtree_rectangle.m_min[segment_id * 2] = italicC.r[segment_id - 1] + 1; // even min id
		rtree_rectangle.m_min[segment_id * 2 + 1] = *min_max.first;// odd first is min value
		rtree_rectangle.m_max[segment_id * 2] = italicC.r[segment_id]; // even max id
		rtree_rectangle.m_max[segment_id * 2 + 1] = *min_max.second;// odd second is max value
#ifdef _DEBUG
		//cout << rtree_rectangle.m_min[segment_id * 2] << " " << rtree_rectangle.m_min[segment_id * 2 + 1] << " " << rtree_rectangle.m_max[segment_id * 2] << " " << rtree_rectangle.m_max[segment_id * 2 + 1] << endl;
		assert(rtree_rectangle.m_min[segment_id * 2] <= rtree_rectangle.m_max[segment_id * 2]);
#endif
		/*--------------------For high dimension dataset--------------------*/
		/*so far most data is one dimension time series*/
		if (dimension_id > 0) {
			int extend_right_endpoint = time_series_vector.size() * dimension_id;
			italicC.r[segment_id] += extend_right_endpoint;
			rtree_rectangle.m_min[segment_id * 2] += extend_right_endpoint;
			rtree_rectangle.m_max[segment_id * 2] = italicC.r[segment_id];
		}
		/*-----------------------------------------------------------------*/
	}

	return rtree_rectangle;
}

//***************************************************************
// Method:get_minmax_apca
// Qualifier: get min max point of every segment, Use original MBR, not APCA paper
// Input: 
// Output: Approximated time series APCA paper. MBR min id == max id
// Note:
// date: 191126 16:15
// author:
//***************************************************************
TEMPLATE
template<typename T, typename Y>
Y& APCA_QUAL::get_minmax_apca(const vector<DataType>& const time_series_vector, const int& const dimension_id, const T& const italicC, Y& const rtree_rectangle) {
#ifdef _DEBUG
	assert(!time_series_vector.empty() && italicC.segmentNum != INF && time_series_vector.size() == italicC.r[italicC.segmentNum - 1] + 1);
#endif
	//int rectangle_length = italicC.segmentNum << 1;
	//odd is value, even is id
	//rtree_rectangle.m_min = new DataType[rectangle_length];
	//rtree_rectangle.m_max = new DataType[rectangle_length];

	const auto min_max = std::minmax_element(time_series_vector.begin(), time_series_vector.begin() + italicC.r[0] + 1);//[begin, end)
	rtree_rectangle.m_min[0] = rtree_rectangle.m_max[0] = italicC.r[0]; // even min id, even max id
	rtree_rectangle.m_min[1] = *min_max.first;// odd first is min value
	rtree_rectangle.m_max[1] = *min_max.second;// odd second is max value

#ifdef _DEBUG
	//cout << rtree_rectangle.m_min[0] <<" "<< rtree_rectangle.m_min[1] <<" "<< rtree_rectangle.m_max[0] <<" "<< rtree_rectangle.m_max[1] << endl;
	assert(rtree_rectangle.m_min[0] == rtree_rectangle.m_max[0]);
	double test_min_value = INF;
	double test_max_value = -INF;
	for_each(time_series_vector.begin(), time_series_vector.begin() + italicC.r[0] + 1, [&](const auto& au) {
		//cout << au <<",";
		test_min_value = min(test_min_value, au);
		test_max_value = max(test_max_value, au);
		});
	//cout << rtree_rectangle.m_min[1] << " " << rtree_rectangle.m_max[1] << endl;
	assert(test_min_value == rtree_rectangle.m_min[1] && test_max_value == rtree_rectangle.m_max[1]);
#endif

	for (int segment_id = 1; segment_id < italicC.segmentNum; segment_id++) {
#ifdef _DEBUG
		assert(italicC.r[segment_id] != INF && italicC.v[segment_id] != INF);
#endif
		const auto min_max = std::minmax_element(time_series_vector.begin() + italicC.r[segment_id - 1] + 1, time_series_vector.begin() + italicC.r[segment_id] + 1);
		rtree_rectangle.m_min[segment_id * 2] = rtree_rectangle.m_max[segment_id * 2] = italicC.r[segment_id]; // even min id // even max id
		rtree_rectangle.m_min[segment_id * 2 + 1] = *min_max.first;// odd first is min value
		rtree_rectangle.m_max[segment_id * 2 + 1] = *min_max.second;// odd second is max value
#ifdef _DEBUG
		//cout << rtree_rectangle.m_min[segment_id * 2] << " " << rtree_rectangle.m_min[segment_id * 2 + 1] << " " << rtree_rectangle.m_max[segment_id * 2] << " " << rtree_rectangle.m_max[segment_id * 2 + 1] << endl;
		assert(rtree_rectangle.m_min[segment_id * 2] <= rtree_rectangle.m_max[segment_id * 2]);
		test_min_value = INF;
		test_max_value = -INF;
		for_each(time_series_vector.begin() + italicC.r[segment_id - 1] + 1, time_series_vector.begin() + italicC.r[segment_id] + 1, [&](const auto& au) {
			test_min_value = min(test_min_value, au);
			test_max_value = max(test_max_value, au);
			});
		assert(test_min_value == rtree_rectangle.m_min[segment_id * 2 + 1] && test_max_value == rtree_rectangle.m_max[segment_id * 2 + 1]);
#endif
		/*--------------------For high dimension dataset--------------------*/
		/*so far most data is one dimension time series*/
		if (dimension_id > 0) {
			int extend_right_endpoint = time_series_vector.size() * dimension_id;
			italicC.r[segment_id] += extend_right_endpoint;
			rtree_rectangle.m_min[segment_id * 2] += extend_right_endpoint;
			rtree_rectangle.m_max[segment_id * 2] = italicC.r[segment_id];
		}
		/*-----------------------------------------------------------------*/
	}

	return rtree_rectangle;
}

//***************************************************************
// Method:get_minmax_apca
// Qualifier: get min max point of every segment, Use original MBR, not APCA paper
// Input: 
// Output: Approximated time series APCA paper. MBR min id == max id
// Note:
// date: 191126 16:15
// author:
//***************************************************************
TEMPLATE
template<typename T, typename Y>
Y& APCA_QUAL::get_minmax_apca(const vector<DataType>& const time_series_vector, const T& const italicC, Y& const rtree_rectangle) {
#ifdef _DEBUG
	assert(!time_series_vector.empty() && italicC.segmentNum != INF && time_series_vector.size() == italicC.r[italicC.segmentNum - 1] + 1);
#endif
	//int rectangle_length = italicC.segmentNum << 1;
	//odd is value, even is id
	//rtree_rectangle.m_min = new DataType[rectangle_length];
	//rtree_rectangle.m_max = new DataType[rectangle_length];

	const auto min_max = std::minmax_element(time_series_vector.begin(), time_series_vector.begin() + italicC.r[0] + 1);//[begin, end)
	rtree_rectangle.m_min[0] = rtree_rectangle.m_max[0] = italicC.r[0]; // even min id, even max id
	rtree_rectangle.m_min[1] = *min_max.first;// odd first is min value
	rtree_rectangle.m_max[1] = *min_max.second;// odd second is max value

#ifdef _DEBUG
	//cout << rtree_rectangle.m_min[0] <<" "<< rtree_rectangle.m_min[1] <<" "<< rtree_rectangle.m_max[0] <<" "<< rtree_rectangle.m_max[1] << endl;
	assert(rtree_rectangle.m_min[0] == rtree_rectangle.m_max[0]);
	double test_min_value = INF;
	double test_max_value = -INF;

	for_each(time_series_vector.begin(), time_series_vector.begin() + italicC.r[0] + 1, [&](const auto& au) {
		//cout << au <<",";
		test_min_value = min(test_min_value, au);
		test_max_value = max(test_max_value, au);
		});
	//cout << rtree_rectangle.m_min[1] << " " << rtree_rectangle.m_max[1] << endl;
	assert(test_min_value == rtree_rectangle.m_min[1] && test_max_value == rtree_rectangle.m_max[1]);
#endif

	for (int segment_id = 1; segment_id < italicC.segmentNum; segment_id++) {
#ifdef _DEBUG
		assert(italicC.r[segment_id] != INF && italicC.v[segment_id] != INF);
#endif
		const auto min_max = std::minmax_element(time_series_vector.begin() + italicC.r[segment_id - 1] + 1, time_series_vector.begin() + italicC.r[segment_id] + 1);
		rtree_rectangle.m_min[segment_id * 2] = rtree_rectangle.m_max[segment_id * 2] = italicC.r[segment_id]; // even min id // even max id
		rtree_rectangle.m_min[segment_id * 2 + 1] = *min_max.first;// odd first is min value
		rtree_rectangle.m_max[segment_id * 2 + 1] = *min_max.second;// odd second is max value
//#ifdef _DEBUG
		//cout << rtree_rectangle.m_min[segment_id * 2] << " " << rtree_rectangle.m_min[segment_id * 2 + 1] << " " << rtree_rectangle.m_max[segment_id * 2] << " " << rtree_rectangle.m_max[segment_id * 2 + 1] << endl;
		/*assert(rtree_rectangle.m_min[segment_id * 2] <= rtree_rectangle.m_max[segment_id * 2]);
		test_min_value = INF;
		test_max_value = -INF;
		for_each(time_series_vector.begin() + italicC.r[segment_id - 1] + 1, time_series_vector.begin() + italicC.r[segment_id] + 1, [&](const auto& au) {
			test_min_value = min(test_min_value, au);
			test_max_value = max(test_max_value, au);
		});
		assert(test_min_value == rtree_rectangle.m_min[segment_id * 2 + 1] && test_max_value == rtree_rectangle.m_max[segment_id * 2 + 1]);*/
		//#endif
	}

	return rtree_rectangle;
}

TEMPLATE
typename APCA_QUAL::APCA& APCA_QUAL::getMBR(const APCA& Cmin, const APCA& Cmax, APCA& MBR) {//N=2M
#ifdef _DEBUG
	assert(Cmin.segmentNum == Cmax.segmentNum);	//cout << "getMBR()\n";
#endif
	MBR.segmentNum = Cmin.segmentNum * 2;
	int i = NULL;
	for (i = 0; i < MBR.segmentNum; i++) {
		//assert(Cmin.r[i / 2] == Cmax.r[i / 2]);
		if (i & 1) {							//odd odd is Cv
#ifdef _DEBUG
			assert(i % 2 != 0);
#endif
			MBR.r[i] = Cmin.v[(i - 1) / 2]; //r is min value
			MBR.v[i] = Cmax.v[(i - 1) / 2];  //v is max value
		}
		else {									//even is Cr
#ifdef _DEBUG
			assert(i % 2 == 0);
#endif
			MBR.r[i] = Cmin.r[i / 2]; // r: min id
			MBR.v[i] = Cmax.r[i / 2]; // v: max id
			//assert(Cmin.r[i / 2] == Cmax.r[i / 2]);
		}

		//cout << "L[" << i << "] = " << MBR.r[i] << ", H[" << i << "] = " << MBR.v[i] << endl;
	}
	return MBR;
}

TEMPLATE
RTREE& APCA_QUAL::buidRTreeIndex(RTREE& const APCARTree, const DataType& g_time_series_length, const DataType& g_index_point_number, APCA_ARRAY* APCALinkOriginal, const string& file_name) {
	//RTree<DataType, ElementType>& buidRTreeIndex(RTree<DataType, ElementType> &APCARTree, double* g_query_time_series, const double& g_time_series_length, const double& g_index_point_number, double(&test_d_original_time_series)[ROW][COLUMN], Link *APCALinkOriginal) {
	printf(">>>***Build RTree Index***<<<\n");
	assert(PAA_or_APCA == 0 || PAA_or_APCA == 1);
	assert(g_time_series_length > 0 && g_index_point_number > 0);
	TOOL::recordStartTime(TOOL::time_record[12]);//whole build_rtree_time

	int i = NULL, j = NULL, f_insert_count = NULL;

	DataType* original_time_series = new DataType[g_time_series_length];

	APCA CminParameter, CmaxParameter, MBRParameter; //Temp MBR

	initialAPCA(CminParameter, APCARTree.NUMDIMS / 2);
	initialAPCA(CmaxParameter, APCARTree.NUMDIMS / 2);
	initialAPCA(MBRParameter, APCARTree.NUMDIMS);

	//memory_account[0] = sizeof(*original_time_series) * int(g_time_series_length);
	//memory_account[1] = sizeof(CminParameter.r) * APCARTree.NUMDIMS / 2 + sizeof(CminParameter.v) * APCARTree.NUMDIMS / 2 + sizeof(CmaxParameter.r) * APCARTree.NUMDIMS / 2 + sizeof(CmaxParameter.v) * APCARTree.NUMDIMS / 2 + sizeof(MBRParameter.r) * APCARTree.NUMDIMS + sizeof(MBRParameter.v) * APCARTree.NUMDIMS;
	/*cout << "Query Point : ";
	for (i = 0; i < g_time_series_length; i++) cout << g_query_time_series[i] << ", ";
	cout << endl;*/

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	f_insert_count = 0;
	for (i = 0; i < g_index_point_number && (!file_stream.eof()) && file_stream.is_open() && file_stream.good(); i++) {
		//getRandomAPCAPoint(APCALinkOriginal[i].originalLink, g_time_series_length);
		//APCALinkOriginal[i].originalLink = test_d_original_time_series[i];
		//getNormalArray(APCALinkOriginal[i].originalLink, g_time_series_length);

		file_stream >> fs_row_string;
		//memory_account[2] = fs_row_string.size();
		stringstream sstr(fs_row_string);

		int string_id = -1;
		while (getline(sstr, fs_row_number, ',') && string_id < g_time_series_length) {
			if (string_id > -1)
				original_time_series[string_id] = stod(fs_row_number);
			//cout << original_time_series[string_id] << ", ";
			string_id++;
		}
		//cout << endl;

		TOOL::normalizeStandard(g_time_series_length, original_time_series);//z-score normalization

		/*cout << "Original Time series: \n";
		TOOL::printArray(original_time_series, g_time_series_length);*/

		if (PAA_or_APCA == 0) {
			//cout << "PAA!" << endl;
			divideRemainderPAA(original_time_series, APCALinkOriginal[i].APCALink, g_time_series_length, APCARTree.NUMDIMS / 2);
		}
		else if (PAA_or_APCA == 1) {
			//cout << "APCA!" << endl;
			getAPCAPoint(original_time_series, g_time_series_length, APCARTree.NUMDIMS / 2, APCALinkOriginal[i].APCALink);
		}
		else {
			assert(0);
		}

		getMBR(getCmin(original_time_series, APCALinkOriginal[i].APCALink, CminParameter), getCmax(original_time_series, APCALinkOriginal[i].APCALink, CmaxParameter), MBRParameter);

		APCARTree.Insert(MBRParameter.r, MBRParameter.v, i);

		f_insert_count++;
	}

	file_stream.close();
	fs_row_string.clear();
	fs_row_string.shrink_to_fit();
	fs_row_number.clear();
	fs_row_number.shrink_to_fit();

	/*cout << "Root Node : sub node number = " << APCARTree.m_root->m_count << " Root level = : " << APCARTree.m_root->m_level << "\n\nBegin to build a RTree:\n";
	cout << "\n RTree conclusion\n The number of RTree Data Point = : " << APCARTree.Count() << endl;*/

	deleteAPCA(CminParameter);
	deleteAPCA(CmaxParameter);
	deleteAPCA(MBRParameter);
	delete[] original_time_series;
	original_time_series = nullptr;

	//***** Delete Leaf node memory   ****//
	RTREE::Iterator it;
	for (APCARTree.GetFirst(it); !APCARTree.IsNull(it); APCARTree.GetNext(it)) {
		it.deleteLeafNodeMemory();
	}

	TOOL::recordFinishTime(TOOL::time_record[12], input_argument.build_rtree_time);
	cout << "RTree build time: " << input_argument.build_rtree_time << " us" << endl;

	if (PAA_or_APCA == 0) {
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[2], input_argument.build_rtree_time);
	}
	else if (PAA_or_APCA == 1) {
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[10], input_argument.build_rtree_time);
	}
	else {
		assert(0);
	}
	return APCARTree;
}


#endif

