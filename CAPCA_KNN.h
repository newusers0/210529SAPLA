#pragma once
#ifndef CAPCA_KNN_H
#define CAPCA_KNN_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "CAPCA.h"

TEMPLATE
class CAPCA_KNN : public APCA_QUAL, public RTREE, virtual public TOOL {
public:

	struct TIME;
	TIME time_record[20];
	struct INPUT_ARGUMENT;
	struct REGION;
	struct APCA_NODE_PAIR;
	struct ORIGINAL_TIME_SERIES_PAIR;
	struct priorityDecreasing;
	struct priorityIncrement;
	struct priorityDistanceEUC;
	struct incrementCMP;

	struct TOOL::INPUT_ARGUMENT input_argument;
	struct TOOL::OUTPUT_ARGUMENT output_argument;//190619

	CAPCA_KNN(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const int& arity_d, const string& read_file_name, string*& const write_file_name);
	CAPCA_KNN() {};
	~CAPCA_KNN();
	void initialREGION(REGION& region_G, const DataType& N);
	void deleteREGION(REGION& region_G);
	double distanceEUC(const double* Q, const DataType& n, double*& const C, const DataType& m);
	double distanceLB(const typename APCA_QUAL::APCA& QProjection, const typename APCA_QUAL::APCA& italicC);
	//Max Sum deviaiton
	double& distanceAE(DataType*& const orginal_array, const int& arrary_length, const typename APCA_QUAL::APCA& italicC, double& distance_AE, double& deviation_max); //For APCA&PAA
	//191206
	double distanceAE(const vector<DataType>& const orginal_array, const int& const arrary_length, const typename APCA_QUAL::APCA& const italicC); //For APCA&PAA

	template<typename T, typename Y, typename U>
	long double distanceAE(const vector<T>& const orginal_array, const int& const arrary_length, const Y& const italicC, U& const result_collection);

	double& distanceAEVector(DataType*& const orginal_array, const int& arrary_length, vector<typename APCA_QUAL::APCA>& italicC, double& distance_AE, double& deviation_max);//190501 For APCA&PAA

	//200115
	//template<typename T>
	void get_APCA_reconstruction(const typename APCA_QUAL::APCA& const italicC, vector<DataType>& const reconstruct_time_series);

	void getFileStreamByID(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, DataType*& const original_time_series);
	void getTXTStreamSpace(const string& const file_name, const int& const array_length, DataType*& const original_time_series);

	REGION& getRegionG(const RTREE::Rect& MBR, REGION& G);// The boundary of a region G = {G[1],G[2],G[3],G[4]}
	//191117 for original Rtree MBR, max id != min id (int APCA paper, min id == max id)
	APCA_KNN_QUAL::REGION& get_region_G_original(const RTREE::Rect& const MBR, APCA_KNN_QUAL::REGION& G);
	double minDistQRt(const double* Q, const REGION& G, const int& i);//MINDIST(Q,R,t)

	template<typename T, typename Y>
	double minDistQRt(const vector<T>& const Q, const REGION& G, const Y& i);//MINDIST(Q,R,t) vector

	double MINDISTQR(const double* Q, const DataType& n, const REGION& G);//MINDIST(Q,R)
	//200923 add template
	template<typename T, typename Y>
	double MINDISTQR(const vector<T>& const Q, const Y& n, const REGION& G);//MINDIST(Q,R) vector
	//project query time series to one APCA approximation
	typename APCA_QUAL::APCA& QAPCAProjection(const DataType* Q, const double& n, typename CAPCA<DataType>::APCA& italicC, typename CAPCA<DataType>::APCA& QProjection);

	/*==================================================================KNN==================================================================================================*/
	bool APCAKNNSearch2(const DataType* g_query_time_series, const DataType& g_index_point_number, const DataType& g_time_series_length, const RTREE& apcaRTree, typename APCA_QUAL::APCA_ARRAY* APCALinkOriginal, const int& K, const string& file_name, list<ORIGINAL_TIME_SERIES_PAIR>& const result);
	bool APCAKNNMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const DataType* g_query_time_series, const DataType& g_index_point_number, const DataType& g_time_series_length, const int& arity_d, const RTREE& apcaRTree, typename APCA_QUAL::APCA_ARRAY* APCALinkOriginal, const int& K, string*& const multi_file_name);
	void SimpleBaseKNNSearch(const DataType* g_query_time_series, const DataType& m_file_time_series_length, const DataType& mg_d_index_point_number, const int& K, const string& file_name, priority_queue<DataType, vector<DataType>, greater<DataType> >& const q_base_queue);
	/*======================================================================================================================================================================*/

	void writeResult(INPUT_ARGUMENT& input_a, string write_file_name);
	//void recordStartTime(_LARGE_INTEGER& whole_first_time_start, double& whole_first_dqFreq);
	void recordStartTime(TIME& time);
	//double recordFinishTime(_LARGE_INTEGER& whole_first_time_start, _LARGE_INTEGER& whole_first_time_over, double& whole_first_dqFreq, double& whole_first_run_time);
	double& recordFinishTime(TIME& time, double& whole_first_run_time);

	template<typename T>
	void printArray(T*& const test_array, const int& array_length);

	template<typename T>
	void writeSingleResult(const string& write_file_name, T& result);

	template<typename T>
	void normalizeA_B(const DataType& left_endpoint, const DataType right_endpoint, T*& const original_array, const int& array_length, T*& const normalized_array);

	template<typename T>
	double& getAverage(T*& const original_array, const int& const array_length, double& average);

	template<typename T>
	double& getVariance(T*& const original_array, const int& const array_length, double& variance);

	template<typename T>
	void normalizeStandard(T*& const original_array, const int& const array_length, T*& const normalized_array);

	/*template<typename T>
	void approximateOriginalFunctionPAA(const APCA& const italicC, const int& const n, const int& const N, T*& const approximation_PAA);*/

	template<typename T>//for original and reconstruceted time series;
	double& getReconstructionError(T*& const original_time_series, T*& const approximation_time_series, const int& const time_series_length, double& const deviation_sum, double& const deviation_max);

	//For CHEBY
	void compareDiffIteration(DataType*& const orgignal_time_series, const int& const array_length, string*& const write_file_name, const int& segment_begin, const int& segment_end, const int& segment_interval);
	//for PAA
	void getDeviationIterationPAA(DataType*& const original_time_series, const int& const array_length, string*& const write_file_name, const int& segment_begin, const int& segment_end, const int& segment_interval);

	void getDeviationIterationAPCA(DataType*& const original_time_series, const int& const array_length, string*& const write_file_name, const int& segment_begin, const int& segment_end, const int& segment_interval);

private:
};

//TEMPLATE
//struct APCA_KNN_QUAL::TIME {
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

//TEMPLATE
//APCA_KNN_QUAL::~CAPCA_KNN() {
//}

#include "CAPCA_KNN.cpp"

#endif

