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

#include "CAPCA.cpp"

#endif

