#pragma once
#ifndef CAPCA_KNN_H
#define CAPCA_KNN_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "CAPCA.h"
//#include "./lib/RTree.h"//210618

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

	//210618
	template<typename T, typename Y>
	Y& getRegionG(const T& MBR, Y& G);// The boundary of a region G = {G[1],G[2],G[3],G[4]}

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

//#include "CAPCA_KNN.cpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEMPLATE
struct APCA_KNN_QUAL::TIME {
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
struct APCA_KNN_QUAL::INPUT_ARGUMENT {
	int time_series_length = NULL;//n
	int point_dimension = NULL;//N
	int point_number = NULL;
	int remainder = NULL;
	int segment_length_first = NULL;//l=n/m+1
	int segment_length_second = NULL;//l=n/m

	int rtree_max_nodes = NULL;
	int K = NULL;
};

TEMPLATE
struct APCA_KNN_QUAL::REGION {
	int regionNum = NULL;
	double* G1 = nullptr;// min value 
	double* G2 = nullptr;// min id : last segment right endpoint
	double* G3 = nullptr;// max value
	double* G4 = nullptr;// max id : segemnt right endpoint

	~REGION() {
		regionNum = NULL;
		delete[] G1;
		G1 = nullptr;
		delete[] G2;
		G2 = nullptr;
		delete[] G3;
		G3 = nullptr;
		delete[] G4;
		G4 = nullptr;
	}
};

TEMPLATE
struct APCA_KNN_QUAL::APCA_NODE_PAIR {//queue data structure.
	double d_dist = NULL;//dist.
	int original_time_series_id = NULL;// APCA point ID.
	APCA_QUAL::APCA* p_APCA_point = nullptr;//data point. if NULL, this is internal node.
	RTREE::Node* p_rtree_node = nullptr;//subNode, if NULL, this is apca point.
};

TEMPLATE
struct APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR {  //for temp queue
	double d_dist = NULL;
	//double *p_original_time_series = nullptr;
	int  original_time_series_id = NULL;
};

TEMPLATE
struct APCA_KNN_QUAL::priorityDecreasing {//big to small
	bool operator ()(const APCA_NODE_PAIR& a, const APCA_NODE_PAIR& b) {
		return a.d_dist < b.d_dist;
	}
};
//for priority queue
TEMPLATE
struct APCA_KNN_QUAL::priorityIncrement {//small to big
	virtual bool operator ()(const APCA_NODE_PAIR& a, const APCA_NODE_PAIR& b) {
		return a.d_dist > b.d_dist;
	}
};

TEMPLATE
struct APCA_KNN_QUAL::priorityDistanceEUC {//small to big
	bool operator ()(const ORIGINAL_TIME_SERIES_PAIR& a, const ORIGINAL_TIME_SERIES_PAIR& b) {
		return a.d_dist > b.d_dist;
	}
};

TEMPLATE
struct APCA_KNN_QUAL::incrementCMP {
	bool operator ()(const double& a, const double& b) {
		return a > b;
	}
};

TEMPLATE
APCA_KNN_QUAL::~CAPCA_KNN() {
}

TEMPLATE
APCA_KNN_QUAL::CAPCA_KNN(const int& n, const int& N, const int& point_number, const int& rtree_max_nodes, const int& K, const int& arity_d, const string& read_file_name, string*& const write_file_name) {
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

	input_argument.result_accuracy = 0.0;// KNN result accuracy

	input_argument.IO_cost = 0.0;

	/*initialCHEBYSHEV_SHARE(input_argument, chebyshev_share);
	getCHEBYSHEV_SHARE(input_argument, chebyshev_share);*/
}

TEMPLATE
void APCA_KNN_QUAL::initialREGION(REGION& region_G, const DataType& N) {
	region_G.regionNum = N;
	region_G.G1 = new DataType[N];
	region_G.G2 = new DataType[N];
	region_G.G3 = new DataType[N];
	region_G.G4 = new DataType[N];
}

TEMPLATE
void APCA_KNN_QUAL::deleteREGION(REGION& region_G) {
	delete[] region_G.G1;
	region_G.G1 = nullptr;

	delete[] region_G.G2;
	region_G.G2 = nullptr;

	delete[] region_G.G3;
	region_G.G3 = nullptr;

	delete[] region_G.G4;
	region_G.G4 = nullptr;
}

TEMPLATE
double APCA_KNN_QUAL::distanceEUC(const double* Q, const DataType& n, double*& const C, const DataType& m) {
	//cout << "\n\n distanceEUC()" << endl;
	int i = 0, j = 0;
	double distance = 0, sum = 0;
	assert(n == m);

	//sum += italicC.r[0] * pow((QProjection.v[0] - italicC.v[0]), 2);

	for (i = 0; i < n; i++) {
		//cout << Q[i]<<" " << C[i]<<" " << Q[i] - C[i] << endl;
		sum += pow((Q[i] - C[i]), 2);
	}
	distance = sqrt(sum);
	//cout << "distanceEUC  : " << distance << endl;
	return distance;
}

TEMPLATE
double APCA_KNN_QUAL::distanceLB(const typename APCA_QUAL::APCA& QProjection, const typename APCA_QUAL::APCA& italicC) {
	//cout << QProjection.v[0]<<"   #####" << italicC.v[0] << endl;
	//cout << "\n\n distanceLB() \n";
	int i = 0, j = 0;
	double distance = 0;

	assert(QProjection.segmentNum == italicC.segmentNum);

	double sum = (italicC.r[0] + 1) * (QProjection.v[0] - italicC.v[0]) * (QProjection.v[0] - italicC.v[0]);
	//cout << QProjection.v[i] << " " << italicC.v[i] << " " << QProjection.v[i]+1 << endl;
	//cout << sum << endl;

#ifdef _DEBUG
	assert(italicC.r[0] != INF && QProjection.v[0] != INF && italicC.v[0] != INF && italicC.r[0] == QProjection.r[0]);
#endif

	for (i = 1; i < italicC.segmentNum; i++) {
		//cout<<i<<" "<< italicC.r[i] - italicC.r[i - 1]<<" " << QProjection.v[i]<<" " << italicC.v[i] <<" "<< QProjection.v[i] - italicC.v[i] << endl;
		sum += (italicC.r[i] - italicC.r[i - 1]) * (QProjection.v[i] - italicC.v[i]) * (QProjection.v[i] - italicC.v[i]);
		//cout << sum << endl;
#ifdef _DEBUG
		assert(italicC.r[i] != INF && QProjection.v[i] != INF && italicC.v[i] != INF && italicC.r[i] == QProjection.r[i]);
#endif
	}

	distance = sqrt(sum);
	//cout << "distanceLB  : " << distance << endl;
	return distance;
}

TEMPLATE
double& APCA_KNN_QUAL::distanceAE(DataType*& const orginal_array, const int& arrary_length, const typename APCA_QUAL::APCA& italicC, double& distance_AE, double& deviation_max) {
	double sum = 0;
	double difference = NULL;
	distance_AE = 0;
	deviation_max = 0;

	DataType* reconstruct_time_series = new DataType[arrary_length];

	for (int interval_id = 0; interval_id <= italicC.r[0]; interval_id++) {
		difference = fabs(italicC.v[0] - orginal_array[interval_id]);
		reconstruct_time_series[interval_id] = italicC.v[0];
		//cout <<"difference: "<< difference << endl;
		deviation_max = difference > deviation_max ? difference : deviation_max;
		sum += difference * difference;
	}

	for (int segment_id = 1; segment_id < italicC.segmentNum; segment_id++) {
		//cout<<"segment_id: "<< segment_id <<endl;
		for (int interval_id = 0; interval_id < italicC.r[segment_id] - italicC.r[segment_id - 1]; interval_id++) {
			difference = fabs(italicC.v[segment_id] - orginal_array[int(interval_id + italicC.r[segment_id - 1] + 1)]);
			reconstruct_time_series[int(interval_id + italicC.r[segment_id - 1] + 1)] = italicC.v[segment_id];
			/*cout << "difference: " << difference << endl;*/
			deviation_max = difference > deviation_max ? difference : deviation_max;
			sum += difference * difference;
		}
	}

	distance_AE = sqrt(sum);
	//cout << "APCA deviation: " << distance_AE << " " << deviation_max << endl;

	TOOL::deleteArray(reconstruct_time_series);
	return distance_AE;
}


TEMPLATE
double APCA_KNN_QUAL::distanceAE(const vector<DataType>& const orginal_array, const int& const arrary_length, const typename APCA_QUAL::APCA& const italicC) { //For APCA&PAA

	double sum = 0;
	double difference = NULL;
	double distance_AE = INF;
	double deviation_max = 0;

	vector<DataType> reconstruct_time_series(arrary_length, INF);

	for (int interval_id = 0; interval_id <= italicC.r[0]; interval_id++) {
		difference = fabs(italicC.v[0] - orginal_array[interval_id]);
		reconstruct_time_series[interval_id] = italicC.v[0];
		//cout <<"difference: "<< difference << endl;
		deviation_max = difference > deviation_max ? difference : deviation_max;
		sum += difference * difference;
	}

	for (int segment_id = 1; segment_id < italicC.segmentNum; segment_id++) {
		//cout<<"segment_id: "<< segment_id <<endl;
		for (int interval_id = 0; interval_id < italicC.r[segment_id] - italicC.r[segment_id - 1]; interval_id++) {
			difference = fabs(italicC.v[segment_id] - orginal_array[int(interval_id + italicC.r[segment_id - 1] + 1)]);
			reconstruct_time_series[int(interval_id + italicC.r[segment_id - 1] + 1)] = italicC.v[segment_id];
			/*cout << "difference: " << difference << endl;*/
			deviation_max = difference > deviation_max ? difference : deviation_max;
			sum += difference * difference;
		}
	}

	distance_AE = sqrt(sum);
	assert(distance_AE != INF);
	//cout << "APCA deviation: " << distance_AE << " " << deviation_max << endl;

	assert(float(distance_AE) == float(TOOL::getDeviation(orginal_array, reconstruct_time_series)));

	reconstruct_time_series.clear();
	reconstruct_time_series.shrink_to_fit();

	return distance_AE;
}

TEMPLATE
template<typename T, typename Y, typename U>
long double APCA_KNN_QUAL::distanceAE(const vector<T>& const orginal_array, const int& const arrary_length, const Y& const italicC, U& const result_collection) {
	long double sum = 0;
	long double difference = NULL;
	long double distance_AE = INF;
	long double deviation_max = -INF;

	result_collection.sum_deviation = 0;
	result_collection.max_deviation = 0;
	result_collection.max_deviation_multiple_width = 0;

	vector<DataType> reconstruct_time_series(arrary_length, INF);

	for (int interval_id = 0; interval_id <= italicC.r[0]; interval_id++) {
		difference = fabs(italicC.v[0] - orginal_array[interval_id]);
		reconstruct_time_series[interval_id] = italicC.v[0];
		//cout <<"difference: "<< difference << endl;
		deviation_max = difference > deviation_max ? difference : deviation_max;
		sum += difference * difference;
	}

	result_collection.max_deviation += deviation_max;
	result_collection.max_deviation_multiple_width += deviation_max * (italicC.r[0] + 1);

	for (int segment_id = 1; segment_id < italicC.segmentNum; segment_id++) {
		//cout<<"segment_id: "<< segment_id <<endl;
		deviation_max = -INF;
		for (int interval_id = 0; interval_id < italicC.r[segment_id] - italicC.r[segment_id - 1]; interval_id++) {
			difference = fabs(italicC.v[segment_id] - orginal_array[int(interval_id + italicC.r[segment_id - 1] + 1)]);
			reconstruct_time_series[int(interval_id + italicC.r[segment_id - 1] + 1)] = italicC.v[segment_id];
			/*cout << "difference: " << difference << endl;*/
			deviation_max = difference > deviation_max ? difference : deviation_max;
			sum += difference * difference;
		}
		result_collection.max_deviation += deviation_max;
		result_collection.max_deviation_multiple_width += deviation_max * (italicC.r[segment_id] - italicC.r[segment_id - 1]);
	}
	result_collection.max_deviation_av = result_collection.max_deviation / double(italicC.segmentNum);
	//result_collection.max_deviation = result_collection.max_deviation_av;

	result_collection.sum_deviation = distance_AE = sqrt(sum);
	assert(distance_AE != INF);
	//cout << "APCA deviation: " << distance_AE << " " << deviation_max << endl;

	assert(float(distance_AE) == float(TOOL::getDeviation(orginal_array, reconstruct_time_series)));

	reconstruct_time_series.clear();
	reconstruct_time_series.shrink_to_fit();

	return distance_AE;
}


TEMPLATE
double& APCA_KNN_QUAL::distanceAEVector(DataType*& const orginal_array, const int& arrary_length, vector<typename APCA_QUAL::APCA>& italicC, double& distance_AE, double& deviation_max) {//190501
	double sum = 0;
	double absolut_difference_sum = 0;
	double difference = NULL;
	distance_AE = 0;
	deviation_max = 0;

	//DataType* reconstruct_time_series = new DataType[arrary_length];
	vector<DataType> reconstruct_time_series;
	/*======================190619==================================*/
	output_argument.sum_area = 0;
	output_argument.sum_density = 0;
	double rectangle_width = 0;
	double rectangle_height = 0;
	std::vector<double> absolute_difference_vector;//190619
	std::vector<double> segment_time_series_vector;//190619
	/*..............................................................*/

	/*======================First segment=========================*/
	for (int interval_id = 0; interval_id <= *italicC[0].r; interval_id++) {
		difference = fabs(*italicC[0].v - orginal_array[interval_id]);
		absolute_difference_vector.push_back(*italicC[0].v - orginal_array[interval_id]);
		//reconstruct_time_series[interval_id] = *italicC[0].v;
		reconstruct_time_series.push_back(*italicC[0].v);
		//cout <<"difference: "<< difference << endl;
		deviation_max = difference > deviation_max ? difference : deviation_max;
		sum += difference * difference;
		absolut_difference_sum += difference;
	}
	rectangle_width = *italicC[0].r + 1;
	auto min_max = minmax_element(absolute_difference_vector.begin(), absolute_difference_vector.end());
	//cout << *min_max.first << " " << *min_max.second << endl;
	rectangle_height = fabs(*min_max.first) + fabs(*min_max.second);
	output_argument.sum_density += 1.0 / rectangle_height;
	output_argument.sum_area += rectangle_height * rectangle_width;
	absolute_difference_vector.clear();
	/*............................................................*/

	/*======================Rest segment=========================*/
	for (int segment_id = 1; segment_id < italicC.size(); segment_id++) {
		//cout<<"segment_id: "<< segment_id <<endl;
		for (int interval_id = 0; interval_id < *italicC[segment_id].r - *italicC[segment_id - 1].r; interval_id++) {
			difference = fabs(*italicC[segment_id].v - orginal_array[int(interval_id + *italicC[segment_id - 1].r + 1)]);
			absolute_difference_vector.push_back(*italicC[segment_id].v - orginal_array[int(interval_id + *italicC[segment_id - 1].r + 1)]);
			segment_time_series_vector.push_back(orginal_array[int(interval_id + *italicC[segment_id - 1].r + 1)]);
			//reconstruct_time_series[int(interval_id + *italicC[segment_id - 1].r + 1)] = *italicC[segment_id].v;
			reconstruct_time_series.push_back(*italicC[segment_id].v);
			//cout << "difference: " << difference << endl;
			deviation_max = difference > deviation_max ? difference : deviation_max;
			sum += difference * difference;
			absolut_difference_sum += difference;
		}
		rectangle_width = *italicC[segment_id].r - *italicC[segment_id - 1].r;
		auto min_max = minmax_element(absolute_difference_vector.begin(), absolute_difference_vector.end());
		auto min_max_series = minmax_element(segment_time_series_vector.begin(), segment_time_series_vector.end());
		rectangle_height = *min_max.second - *min_max.first;
		double test_height = *min_max_series.second - *min_max_series.first;
		assert(float(test_height) == float(rectangle_height));
		output_argument.sum_density += 1.0 / rectangle_height;
		output_argument.sum_area += rectangle_height * rectangle_width;
		absolute_difference_vector.clear();
		segment_time_series_vector.clear();
	}
	/*............................................................*/
	//cout << "sum: " << sum << endl;
	distance_AE = sqrt(sum);
	//cout << "PAA deviation: " << distance_AE << " " << deviation_max << endl;

	for (auto&& au : italicC) {
		delete au.v;
		au.v = nullptr;
		delete au.r;
		au.r = nullptr;
	}

	TOOL::getDeviation(orginal_array, reconstruct_time_series, arrary_length, output_argument);
	assert(distance_AE == output_argument.sum_deviation && deviation_max == output_argument.max_deviation);
	//distance_AE = absolut_difference_sum;
	//TOOL::deleteArray(reconstruct_time_series);
	reconstruct_time_series.clear();
	reconstruct_time_series.shrink_to_fit();
	return distance_AE;
}

//200115
TEMPLATE
//template<typename T>
void APCA_KNN_QUAL::get_APCA_reconstruction(const typename APCA_QUAL::APCA& const italicC, vector<DataType>& const reconstruct_time_series) {
	for (int interval_id = 0; interval_id <= italicC.r[0]; interval_id++) {
		//approximated_time_series[interval_id] = italicC.v[0];
		reconstruct_time_series.emplace_back(italicC.v[0]);
		//cout <<"difference: "<< difference << endl;
		//deviation_max = difference > deviation_max ? difference : deviation_max;
		//sum += difference * difference;
	}

	for (int segment_id = 1; segment_id < italicC.segmentNum; segment_id++) {
		//cout<<"segment_id: "<< segment_id <<endl;
		for (int interval_id = 0; interval_id < italicC.r[segment_id] - italicC.r[segment_id - 1]; interval_id++) {
			reconstruct_time_series.emplace_back(italicC.v[segment_id]);
			//approximated_time_series[int(interval_id + italicC.r[segment_id - 1] + 1)] = italicC.v[segment_id];
			/*cout << "difference: " << difference << endl;*/
			//deviation_max = difference > deviation_max ? difference : deviation_max;
			//sum += difference * difference;
		}
	}
}

TEMPLATE
void APCA_KNN_QUAL::getFileStreamByID(const string& file_name, const double& g_time_series_length, const int& const original_time_series_id, DataType*& const original_time_series) {
	string fs_row_string;
	string fs_row_number;

	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	int i = 0;
	while (!file_stream.eof() && file_stream.is_open())
	{
		file_stream.good();
		file_stream.fail();
		file_stream >> fs_row_string;

		if (i == original_time_series_id) {
			stringstream sstr(fs_row_string);
			int string_id = -1;
			while (getline(sstr, fs_row_number, ',') && string_id < g_time_series_length) {
				if (string_id > -1) {
					original_time_series[string_id] = stod(fs_row_number);
					assert(original_time_series[string_id] != NULL);
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

TEMPLATE
void APCA_KNN_QUAL::getTXTStreamSpace(const string& const file_name, const int& const array_length, DataType*& const original_time_series) {
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


TEMPLATE
typename APCA_KNN_QUAL::REGION& APCA_KNN_QUAL::getRegionG(const RTREE::Rect& MBR, REGION& G) {
#ifdef _DEBUG
	//assert(G.G2[i] <= G.G4[i] && G.G1[i] <= G.G3[i]);
	double run_time0;
	recordStartTime(time_record[11]);
#endif

	int i = 0;
	//cout << "getRegionG()\n";

	/*if (!NUMDIMS % 2) {
	cout << "NUMDIMS is not even    Wrong!!!!!";
	}*/
	G.G1[i] = MBR.m_min[1];// min value 
	G.G2[i] = 0;// min id : last segment right endpoint
	G.G3[i] = MBR.m_max[1];// max value
	G.G4[i] = MBR.m_max[0];// max id : segemnt right endpoint
#ifdef _DEBUG
	//cout << "min["<<0<< "] = " << MBR.m_min[0] << ", max[" << 0 << "] = " << MBR.m_max[0] << ", min[" << 1 << "] = " << MBR.m_min[1]  << ", max[" << 1 << "] = " << MBR.m_max[1] << endl;
	//cout << "G2[" << i << "] = " << G.G2[i] << ", G1[" << i << "] = " << G.G1[i] << ", G4[" << i << "] = " << G.G4[i] << ", G3[" << i << "] = " << G.G3[i] << endl;
#endif
	for (i = 1; i < G.regionNum; i++) {
		//cout << MBR.m_min[i * 2] << " "<< MBR.m_max[i * 2] << endl;
		//assert(MBR.m_min[i * 2] == MBR.m_max[i * 2]);
		G.G1[i] = MBR.m_min[i * 2 + 1];// min value 
		G.G2[i] = MBR.m_min[i * 2 - 2] + 1;// min id : last segment right endpoint
		G.G3[i] = MBR.m_max[i * 2 + 1];// max value
		G.G4[i] = MBR.m_max[i * 2]; //max id : segemnt right endpoint
#ifdef _DEBUG
		//cout << "min[" << i * 2 << "] = " << MBR.m_min[i * 2] << ", max[" << i * 2 << "] = " << MBR.m_max[i * 2] << ", min[" << i * 2 + 1 << "] = " << MBR.m_min[i * 2 + 1] << ", max[" << i * 2 + 1 << "] = " << MBR.m_max[i * 2 + 1] << endl;
		//cout << "G2[" << i << "] = " << G.G2[i] << ", G1[" << i << "] = " << G.G1[i] << ", G4[" << i << "] = " << G.G4[i] << ", G3[" << i << "] = " << G.G3[i] << endl;
#endif
	}

	//#ifdef _DEBUG
		//g_n_time_getMBR_count += recordFinishTime(time_record[11], run_time0);
	//#endif
		//printf("\nrun_time = %f us\n", run_time0);

	return G;
}


TEMPLATE
template<typename T, typename Y>
typename Y& APCA_KNN_QUAL::getRegionG(const T& MBR, Y& G) {
#ifdef _DEBUG
	//assert(G.G2[i] <= G.G4[i] && G.G1[i] <= G.G3[i]);
	double run_time0;
	recordStartTime(time_record[11]);
#endif

	int i = 0;
	//cout << "getRegionG()\n";

	/*if (!NUMDIMS % 2) {
	cout << "NUMDIMS is not even    Wrong!!!!!";
	}*/
	G.G1[i] = MBR.m_min[1];// min value 
	G.G2[i] = 0;// min id : last segment right endpoint
	G.G3[i] = MBR.m_max[1];// max value
	G.G4[i] = MBR.m_max[0];// max id : segemnt right endpoint
#ifdef _DEBUG
	//cout << "min["<<0<< "] = " << MBR.m_min[0] << ", max[" << 0 << "] = " << MBR.m_max[0] << ", min[" << 1 << "] = " << MBR.m_min[1]  << ", max[" << 1 << "] = " << MBR.m_max[1] << endl;
	//cout << "G2[" << i << "] = " << G.G2[i] << ", G1[" << i << "] = " << G.G1[i] << ", G4[" << i << "] = " << G.G4[i] << ", G3[" << i << "] = " << G.G3[i] << endl;
#endif
	for (i = 1; i < G.regionNum; i++) {
		//cout << MBR.m_min[i * 2] << " "<< MBR.m_max[i * 2] << endl;
		//assert(MBR.m_min[i * 2] == MBR.m_max[i * 2]);
		G.G1[i] = MBR.m_min[i * 2 + 1];// min value 
		G.G2[i] = MBR.m_min[i * 2 - 2] + 1;// min id : last segment right endpoint
		G.G3[i] = MBR.m_max[i * 2 + 1];// max value
		G.G4[i] = MBR.m_max[i * 2]; //max id : segemnt right endpoint
#ifdef _DEBUG
		//cout << "min[" << i * 2 << "] = " << MBR.m_min[i * 2] << ", max[" << i * 2 << "] = " << MBR.m_max[i * 2] << ", min[" << i * 2 + 1 << "] = " << MBR.m_min[i * 2 + 1] << ", max[" << i * 2 + 1 << "] = " << MBR.m_max[i * 2 + 1] << endl;
		//cout << "G2[" << i << "] = " << G.G2[i] << ", G1[" << i << "] = " << G.G1[i] << ", G4[" << i << "] = " << G.G4[i] << ", G3[" << i << "] = " << G.G3[i] << endl;
#endif
	}

	//#ifdef _DEBUG
		//g_n_time_getMBR_count += recordFinishTime(time_record[11], run_time0);
	//#endif
		//printf("\nrun_time = %f us\n", run_time0);

	return G;
}


TEMPLATE
typename APCA_KNN_QUAL::REGION& APCA_KNN_QUAL::get_region_G_original(const RTREE::Rect& const MBR, APCA_KNN_QUAL::REGION& const G) {
	int i = 0;
	//cout << "getRegionG()\n";

	/*if (!NUMDIMS % 2) {
	cout << "NUMDIMS is not even    Wrong!!!!!";
	}*/
	//G.G1[i] = MBR.m_min[1];// min value 
	//G.G2[i] = 0;// min id : last segment right endpoint
	//G.G3[i] = MBR.m_max[1];// max value
	//G.G4[i] = MBR.m_max[0];// max id : segemnt right endpoint
	//cout << "G1[" << i << "] = " << G.G1[i] << ", G2[" << i << "] = " << G.G2[i] << ", G3[" << i << "] = " << G.G3[i] << ", G4[" << i << "] = " << G.G4[i] << endl;
	for (i = 0; i < G.regionNum; i++) {
		G.G1[i] = MBR.m_min[i * 2 + 1];// min value 
		G.G2[i] = MBR.m_min[i * 2];// min id : last segment right endpoint
		G.G3[i] = MBR.m_max[i * 2 + 1];// max value
		G.G4[i] = MBR.m_max[i * 2]; // max id : segemnt right endpoint
		//cout << "G1[" << i << "] = " << G.G1[i] << ", G2[" << i << "] = " << G.G2[i] << ", G3[" << i << "] = " << G.G3[i] << ", G4[" << i << "] = " << G.G4[i] << endl;
#ifdef _DEBUG
		assert(G.G2[i] <= G.G4[i] && G.G1[i] <= G.G3[i]);
#endif//is internal Node
	}

	//printf("\nrun_time = %f us\n", run_time0);
	return G;
}

TEMPLATE
double APCA_KNN_QUAL::minDistQRt(const double* Q, const REGION& G, const int& i) {
#ifdef _DEBUG
	assert(Q[i] != INF);
#endif
	int j = 0, count = 0;
	double minDistanceQG = 0;

	priority_queue<double, vector<double>, incrementCMP > tempMinDist;
	//double* f_d_temp_result_array = new double[G.regionNum];

	for (j = 0; j < G.regionNum; j++) {
		//cout << "G2[" << j << "] = " << G.G2[j] << ", G4[" << j << "] = " << G.G4[j] << ", G1[" << j << "] = " << G.G1[j] << ", G3[" << j << "] = " << G.G3[j] << endl;
		if (i >= G.G2[j] && i <= G.G4[j]) {
			//cout << "Active!" << endl;
			if (Q[i] < G.G1[j]) {
				minDistanceQG = pow((G.G1[j] - Q[i]), 2);
				//minDistanceArray[count] = minDistanceQG;
				//count++;
				tempMinDist.push(minDistanceQG);
				//cout << "q[" << i << "] < G1[" << j << "], MINDIST = " << minDistanceQG << "\n\n";
			}
			else if (Q[i] > G.G3[j]) {
				minDistanceQG = pow((Q[i] - G.G3[j]), 2);
				//minDistanceArray[count] = minDistanceQG;
				//count++;
				tempMinDist.push(minDistanceQG);
				//cout << "q[" << i << "] > G3[" << j << "], MINDIST = " << minDistanceQG << "\n\n";
			}
			else {
				minDistanceQG = 0;
				//cout << "G1 >= qt <= G3, MINDIST = " << minDistanceQG << endl;
				//delete[] f_d_temp_result_array;
				return minDistanceQG;
			}
		}
		else {
			//cout << "G is not active!!!!! \n\n";
			continue;
		}
	}
	//minDistanceQG = *min_element(minDistanceArray, minDistanceArray + count);
	minDistanceQG = tempMinDist.top();
	priority_queue<double, vector<double>, incrementCMP >().swap(tempMinDist);

	return minDistanceQG;
}

TEMPLATE
template<typename T, typename Y>
double APCA_KNN_QUAL::minDistQRt(const vector<T>& const Q, const REGION& G, const Y& i) {
#ifdef _DEBUG
	assert(Q[i] != INF);
#endif
	int j = 0, count = 0;
	double minDistanceQG = 0;

	priority_queue<double, vector<double>, incrementCMP > tempMinDist;
	//double* f_d_temp_result_array = new double[G.regionNum];

	for (j = 0; j < G.regionNum; j++) {
		//cout << "G2[" << j << "] = " << G.G2[j] << ", G4[" << j << "] = " << G.G4[j] << ", G1[" << j << "] = " << G.G1[j] << ", G3[" << j << "] = " << G.G3[j] << endl;
		if (i >= G.G2[j] && i <= G.G4[j]) {
			//cout << "Active!" << endl;
			if (Q[i] < G.G1[j]) {
				minDistanceQG = pow((G.G1[j] - Q[i]), 2);
				//minDistanceArray[count] = minDistanceQG;
				//count++;
				tempMinDist.push(minDistanceQG);
				//cout << "q[" << i << "] < G1[" << j << "], MINDIST = " << minDistanceQG << "\n\n";
			}
			else if (Q[i] > G.G3[j]) {
				minDistanceQG = pow((Q[i] - G.G3[j]), 2);
				//minDistanceArray[count] = minDistanceQG;
				//count++;
				tempMinDist.push(minDistanceQG);
				//cout << "q[" << i << "] > G3[" << j << "], MINDIST = " << minDistanceQG << "\n\n";
			}
			else {
				minDistanceQG = 0;
				//cout << "G1 >= qt <= G3, MINDIST = " << minDistanceQG << endl;
				//delete[] f_d_temp_result_array;
				return minDistanceQG;
			}
		}
		else {
			//cout << "G is not active!!!!! \n\n";
			continue;
		}
	}
	//minDistanceQG = *min_element(minDistanceArray, minDistanceArray + count);
	minDistanceQG = tempMinDist.top();
	priority_queue<double, vector<double>, incrementCMP >().swap(tempMinDist);

	return minDistanceQG;
}

TEMPLATE
double APCA_KNN_QUAL::MINDISTQR(const double* Q, const DataType& n, const REGION& G) {
	//cout << "MINDISTQR()\n";
	int i = NULL;
	double minDistanceQR = 0;
	for (i = 0; i < n; i++) {
#ifdef _DEBUG
		assert(Q[i] != INF);
#endif
		minDistanceQR += minDistQRt(Q, G, i);
	}

	return sqrt(minDistanceQR);
}


TEMPLATE
template<typename T, typename Y>
double APCA_KNN_QUAL::MINDISTQR(const vector<T>& const Q, const Y& n, const REGION& G) {
	//cout << "MINDISTQR()\n";
#ifdef _DEBUG
	assert(Q.size() == n);
#endif
	int i = NULL;
	double minDistanceQR = 0;
	for (i = 0; i < n; i++) {
#ifdef _DEBUG
		assert(Q[i] != INF);
#endif
		minDistanceQR += minDistQRt(Q, G, i);
	}

	return sqrt(minDistanceQR);
}


TEMPLATE
typename APCA_QUAL::APCA& APCA_KNN_QUAL::QAPCAProjection(const DataType* Q, const double& n, typename APCA_QUAL::APCA& italicC, typename APCA_QUAL::APCA& QProjection) {
	//cout << "QAPCAProjection()" << endl;
#ifdef _DEBUG
	assert(n == italicC.r[italicC.segmentNum - 1] + 1);
#endif
	int i = 0, j = 0;
	double sum = 0;
	QProjection.segmentNum = italicC.segmentNum;

	QProjection.r[0] = italicC.r[0];
	while (j <= italicC.r[0]) {
		sum += Q[j];
		//cout << "Q[" << j << "] = " << Q[j] << " ";
		j++;
	}
	QProjection.v[i] = sum / (italicC.r[i] + 1);
	//cout << "\nQv[" << i << "] = " << QProjection.v[i] << ", Qr[" << i << "] = " << QProjection.r[i] << "\n\n";

	for (i = 1; i < italicC.segmentNum; i++) {
		QProjection.r[i] = italicC.r[i];
		sum = 0;
		for (j = 1; j <= (italicC.r[i] - italicC.r[i - 1]); j++) {
			sum += Q[j + int(italicC.r[i - 1])];
			//cout << "Q[" << j + italicC.r[i - 1] << "] = " << Q[j + int(italicC.r[i - 1])] << " ";
		}
		QProjection.v[i] = sum / (italicC.r[i] - italicC.r[i - 1]);
		//cout << "\nQv[" << i << "] = " << QProjection.v[i] << ", Qr[" << i << "] = " << QProjection.r[i] << "\n\n";
	}
	return QProjection;
}

TEMPLATE
bool APCA_KNN_QUAL::APCAKNNSearch2(const DataType* g_query_time_series, const DataType& g_index_point_number, const DataType& g_time_series_length, const RTREE& apcaRTree, typename APCA_QUAL::APCA_ARRAY* APCALinkOriginal, const int& K, const string& file_name, list<ORIGINAL_TIME_SERIES_PAIR>& const result) {
	cout << "APCA&PAA KNN Search function: " << endl;
	//g_n_account_apca_point = 0;
	input_argument.IO_cost = 0;// measure I/O cost
	input_argument.knn_total_time = 0.0;
	input_argument.navigate_index_time = 0.0;// navigate time
	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
	input_argument.distance_euc_time = 0.0;// distance euclidean time

	double whole_first_run_time = NULL;
	recordStartTime(time_record[0]);

	APCA_KNN_QUAL::recordStartTime(APCA_KNN_QUAL::time_record[13]);//for total KNN time

	//printf("<///////**  APCA KNN Begin  **////////>\n");
	int i = NULL, j = NULL;
	assert(K <= g_index_point_number);

	priority_queue<APCA_NODE_PAIR, vector<APCA_NODE_PAIR>, priorityIncrement > queue;
	APCA_NODE_PAIR f_APCA_Root, f_temp_APCA_Pair;
	list<ORIGINAL_TIME_SERIES_PAIR> temp;
	//list<ORIGINAL_TIME_SERIES_PAIR> result;
	ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair;

	f_APCA_Root.p_rtree_node = apcaRTree.m_root;
	f_APCA_Root.d_dist = 0;
	queue.push(f_APCA_Root);
	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
	//printf("<///////**    KNN Begin   **////////>\n");

	int n_data_point_count = 0;
	int n_distanceLB_count = 0;
	int n_push_leaf_node_count = 0;

	/*g_n_time_getMBR_count = 0;
	g_n_time_apca_point_count = 0;
	g_n_time_count1 = 0;
	g_n_time_count2 = 0;
	g_n_time_count_while_first_part = 0;

	g_n_time_loop_result_count = 0;
	g_n_time_leaf_node_push_count = 0;
	g_n_time_leaf_node_distanceLB_count = 0;
	g_n_time_leaf_node_assignment_count = 0;
	g_n_time_child_node_MINDIST_count = 0;
	g_n_time_child_node_push_count = 0;

	g_d_time_whole_first_run = recordFinishTime(time_record[0], whole_first_run_time);*/

	double temp_navigate_time = 0;
	TOOL::recordStartTime(TOOL::time_record[4]);//navigate index time

	int n_result_count = 0;
	while (!queue.empty()) {
		/*while (!queue.empty()|| result.size() != K) {*/
		//while (result.size() != K) {
		//cout << "KNN: the " << count << " turn : \n";

		double run_time3;
		//recordStartTime(whole_first_time_start[1], whole_first_dqFreq[1]);
		recordStartTime(time_record[1]);

		/*if (count == K) {
		cout << "////*********EMPTY  Found " << K << " result!!!!!!********///" << endl;break;}*/

		double run_time5;
		//recordStartTime(whole_first_time_start[2], whole_first_dqFreq[2]);
		recordStartTime(time_record[2]);

		APCA_NODE_PAIR m_temp_queue_top = queue.top();

		//if (queue.size() > memory_account[9])
			//memory_account[9] = queue.size();

		//cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator: ";
		//for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
		//cout << endl;

		//if (memory_account[10] < temp.size()) { memory_account[10] = temp.size(); }

		for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
			if (plist->d_dist <= m_temp_queue_top.d_dist) {
				//cout << "           <= " << endl;
				result.push_back(*plist);
				plist = temp.erase(plist);
			}
			else plist++;
			if (K == result.size()) {
				cout << "*****************Find result!!!!!!!!!!!!!  " << result.size() << endl;
				//result.sort(APCA_KNN_QUAL::compare);

				//time of whole KNN procedure
				APCA_KNN_QUAL::recordFinishTime(APCA_KNN_QUAL::time_record[13], input_argument.knn_total_time);

				cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;
				cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;
				cout << "I/O cost: " << input_argument.IO_cost << endl;
				input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
				//cout << "pruning power: "<< input_argument.pruning_power <<endl;

				if (PAA_or_APCA == 0) {
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[0], input_argument.pruning_power);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[1], input_argument.IO_cost);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[3], input_argument.navigate_index_time);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[4], input_argument.distance_lowbound_time);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[5], input_argument.distance_euc_time);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[6], input_argument.knn_total_time);
				}
				else if (PAA_or_APCA == 1) {
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[8], input_argument.pruning_power);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[9], input_argument.IO_cost);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[11], input_argument.navigate_index_time);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[12], input_argument.distance_lowbound_time);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[13], input_argument.distance_euc_time);
					TOOL::writeSingleResult(input_argument, input_argument.write_file_name[14], input_argument.knn_total_time);
				}
				else {
					assert(0);
				}

				result.sort([](const ORIGINAL_TIME_SERIES_PAIR& first, const  ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist; });//small to big

				/*cout << "Result list: ";
				for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
					cout << it->d_dist << ", " << it->original_time_series_id << "; ";
				}*/
				//*************************************Base KNN *******************************************************
/*priority_queue<DataType, vector<DataType>, greater<DataType> > q_base_queue;
SimpleBaseKNNSearch(g_query_time_series, input_argument.time_series_length, input_argument.point_number, input_argument.K, input_argument.read_file_name, q_base_queue);
*/
//input_argument.result_accuracy = 0;// KNN result accuracy
//for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
//	if (it->original_time_series_id == q_base_queue.top().original_time_series_id) {
//		input_argument.result_accuracy++;
//	}

//	q_base_queue.pop();
//}

//double accuracy = input_argument.result_accuracy / double(input_argument.K);
//cout << "KNN result accuracy: " << accuracy << endl;
//if (PAA_or_APCA == 0) {
//	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[7], accuracy);
//}
//else if (PAA_or_APCA == 1) {
//	TOOL::writeSingleResult(input_argument,input_argument.write_file_name[15], accuracy);
//}
//else {
//	assert(0);
//}
//********************************************************************************************

				TOOL::printInputArgument(input_argument);
				//cout << "  count apca point = " << g_n_account_apca_point << " times, p= " << g_n_account_apca_point / double(input_argument.point_number) << endl;
				priority_queue<APCA_NODE_PAIR, vector<APCA_NODE_PAIR>, priorityIncrement >().swap(queue);
				temp.clear();
				//result.clear();
				//list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
				//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
				return true;
			}
		}

		//g_n_time_loop_result_count += recordFinishTime(whole_first_time_start[2], whole_first_time_over[2], whole_first_dqFreq[2], run_time5);
		//g_n_time_loop_result_count += recordFinishTime(time_record[2], run_time5);

		//printf("\nrun_time = %f us\n", run_time0);
		//cout << "Before Pop: Queue.size(): " << queue.size() << " Pop: top.dist:" << queue.top().d_dist << " " << queue.top().original_time_series_id << endl;
		queue.pop();
		//cout << "After Pop: Queue.size(): " << queue.size();
		//if (!queue.empty()) cout << " Pop: top.dist:" << queue.top().d_dist << " " << queue.top().original_time_series_id << endl;
		//else cout << endl;

		//g_n_time_count_while_first_part = recordFinishTime(whole_first_time_start[1], whole_first_time_over[1], whole_first_dqFreq[1], run_time3);
		//g_n_time_count_while_first_part = recordFinishTime(time_record[1], run_time3);

		//printf("\nrun_time = %f us\n", run_time0);

		if (m_temp_queue_top.p_rtree_node == nullptr) { //is APCA data point
			TOOL::recordFinishTime(TOOL::time_record[4], temp_navigate_time);
			input_argument.navigate_index_time = temp_navigate_time;
			cout << TOOL::recordFinishTime(TOOL::time_record[4]) << endl;

			//cout << "    queue.top is data point\n";
			DataType* original_time_series = new DataType[g_time_series_length];

			double run_time0;
			//recordStartTime(whole_first_time_start[3], whole_first_dqFreq[3]);
			recordStartTime(time_record[3]);

			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;
			getFileStreamByID(file_name, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
			input_argument.IO_cost++;
			//g_n_account_apca_point++;

			double distance_euclidean_time = 0;
			TOOL::recordStartTime(TOOL::time_record[6]);//time for euclidean distance

			tempOriginalTimeSeriesPair.d_dist = distanceEUC(g_query_time_series, g_time_series_length, original_time_series, g_time_series_length);

			TOOL::recordFinishTime(TOOL::time_record[6], distance_euclidean_time);
			input_argument.distance_euc_time += distance_euclidean_time;

			delete[] original_time_series;
			original_time_series = nullptr;

			//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
			temp.push_back(tempOriginalTimeSeriesPair);

			//g_n_time_apca_point_count += recordFinishTime(whole_first_time_start[3], whole_first_time_over[3], whole_first_dqFreq[3], run_time0);
			//g_n_time_apca_point_count += recordFinishTime(time_record[3], run_time0);
			//printf("\nrun_time = %f us\n", run_time0);
		}
		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;

			typename APCA_QUAL::APCA QProjection;
			typename APCA_QUAL::initialAPCA(QProjection, apcaRTree.NUMDIMS / 2);

			double run_time1;
			//recordStartTime(whole_first_time_start[4], whole_first_dqFreq[4]);
			recordStartTime(time_record[4]);

			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {
				//g_n_account_leaf_node++;

				double assignment_run_time;
				recordStartTime(time_record[5]);

				f_temp_APCA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);
				//cout << tempAPCAPair.id_originalTimeSeries << endl;
				f_temp_APCA_Pair.p_APCA_point = &APCALinkOriginal[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)].APCALink;
				//cout << "KNN : data point ID = " << tempAPCAPair.id_originalTimeSeries << endl;
				//cout << tempAPCAPair.APCAValue << endl;

				n_distanceLB_count++;
				//g_n_time_leaf_node_assignment_count += recordFinishTime(time_record[5], assignment_run_time);

				//printf("\nrun_time = %f us\n", run_time0);

				double distance_PLA_time = 0;
				TOOL::recordStartTime(TOOL::time_record[5]);//distance index time

				//f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point), *f_temp_APCA_Pair.p_APCA_point);
				f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point, QProjection), *f_temp_APCA_Pair.p_APCA_point);

				TOOL::recordFinishTime(TOOL::time_record[5], distance_PLA_time);
				input_argument.distance_lowbound_time += distance_PLA_time;

				//printf("\nrun_time = %f us\n", run_time0);

				//cout << tempAPCAPair.key;
				f_temp_APCA_Pair.p_rtree_node = nullptr;
				//	cout << tempAPCAPair.key;

				double queue_run_time;
				recordStartTime(time_record[7]);

				//cout << "            Push apca data. Branch id: " << i << " DLB: " << f_temp_APCA_Pair.d_dist << ", time series ID: " << f_temp_APCA_Pair.original_time_series_id << endl;
				queue.push(f_temp_APCA_Pair);

				n_push_leaf_node_count++;
				//g_n_time_leaf_node_push_count += recordFinishTime(time_record[7], queue_run_time);

				//printf("\nrun_time = %f us\n", run_time0);
				//cout << "KNN : queue.size() = " << queue.size() << endl;
				//		//cout << tempAPCAPair.key;
			}
			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
			//	//cout << tempOriginalTimeSeriesPair.key;

			//g_n_time_count1 += recordFinishTime(time_record[4], run_time1);
			//printf("\nrun_time = %f us\n", run_time0);
			APCA_QUAL::deleteAPCA(QProjection);
		}
		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {															//is internal Node
			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
			double run_time2;
			recordStartTime(time_record[8]);

			APCA_NODE_PAIR tempApcaPair;
			REGION fs_region_G;
			initialREGION(fs_region_G, apcaRTree.NUMDIMS / 2);
			////cout << "regionNum = " << G.regionNum << endl;

			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
				//g_n_account_child_node++;

				double mindistQR_run_time;
				recordStartTime(time_record[9]);

				double distance_PLA_time = 0;
				TOOL::recordStartTime(TOOL::time_record[5]);//distance index time

				tempApcaPair.d_dist = MINDISTQR(g_query_time_series, g_time_series_length, getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));

				TOOL::recordFinishTime(TOOL::time_record[5], distance_PLA_time);
				input_argument.distance_lowbound_time += distance_PLA_time;

				/*for (int i = 0; i < fs_region_G.regionNum; i++) {
				cout << "            G1[" << i << "] = " << fs_region_G.G1[i] << ", G2[" << i << "] = " << fs_region_G.G2[i] << ", G3[" << i << "] = " << fs_region_G.G3[i] << ", G4[" << i << "] = " << fs_region_G.G4[i] << endl;
				}*/
				//g_n_time_child_node_MINDIST_count += recordFinishTime(time_record[9], mindistQR_run_time);

				//printf("\nrun_time = %f us\n", run_time0);

				tempApcaPair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;

				double push_run_time;
				recordStartTime(time_record[10]);
				//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
				queue.push(tempApcaPair);
				//cout << "KNN : queue.size() = " << queue.size() << endl;
				//g_n_time_child_node_push_count += recordFinishTime(time_record[10], push_run_time);
				//printf("\nrun_time = %f us\n", run_time0);
			}

			deleteREGION(fs_region_G);
			//g_n_time_count2 += recordFinishTime(time_record[8], run_time2);
			//printf("\nrun_time = %f us\n", run_time0);
		}
		else {
			assert(0);
		}
	}

	cout << "??????????????????     APCA KNN failed    ??????????????????" << endl;
	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
	ofstream outfile(input_argument.write_file_name[16] + ".txt", ios::app);
	assert(outfile.is_open());
	outfile << "??????????????????     APCA KNN failed    ??????????????????   " << PAA_or_APCA << endl;
	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
	outfile.close();
	TOOL::printInputArgument(input_argument);

	if (PAA_or_APCA == 0) {
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[0], "FAIL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[1], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[3], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[4], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[5], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[6], "FALL");
	}
	else if (PAA_or_APCA == 1) {
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[8], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[9], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[11], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[12], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[13], "FALL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[14], "FALL");
	}
	else {
		assert(0);
	}

	//assert(0);
	priority_queue<APCA_NODE_PAIR, vector<APCA_NODE_PAIR>, priorityIncrement >().swap(queue);
	list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
	//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
	return false;
}

TEMPLATE
bool APCA_KNN_QUAL::APCAKNNMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const DataType* g_query_time_series, const DataType& g_index_point_number, const DataType& g_time_series_length, const int& arity_d, const RTREE& apcaRTree, typename APCA_QUAL::APCA_ARRAY* APCALinkOriginal, const int& K, string*& const multi_file_name) {
#ifdef _DEBUG
	cout << "APCA&PAA KNN Search function: " << endl;
	//g_n_account_apca_point = 0;
	input_argument.knn_total_time = 0.0;
	input_argument.navigate_index_time = 0.0;// navigate time
	input_argument.distance_lowbound_time = 0.0; // distance chebyshev, PLA, APCA time
	input_argument.distance_euc_time = 0.0;// distance euclidean time
	TOOL::recordStartTime(TOOL::time_record[13]);//for total KNN time
	assert(K <= g_index_point_number);
#endif
	input_argument.IO_cost = 0;// measure I/O cost
	input_argument.sum_distance_euc = 0;
	TOOL::recordStartTime(TOOL::time_record[14]);
	int i = NULL, j = NULL;
	priority_queue<APCA_NODE_PAIR, vector<APCA_NODE_PAIR>, priorityIncrement > queue;
	APCA_NODE_PAIR f_APCA_Root, f_temp_APCA_Pair;
	list<ORIGINAL_TIME_SERIES_PAIR> temp;
	list<ORIGINAL_TIME_SERIES_PAIR> result;
	ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair;

	f_APCA_Root.p_rtree_node = apcaRTree.m_root;
	f_APCA_Root.d_dist = 0;
	queue.push(f_APCA_Root);
	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
	//printf("<///////**    KNN Begin   **////////>\n");

	TOOL::recordStartTime(TOOL::time_record[4]);//navigate index time

	while (!queue.empty()) {
		APCA_NODE_PAIR m_temp_queue_top = queue.top();

		//cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator: ";
		//for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
		//cout << endl;

		for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
			if (plist->d_dist <= m_temp_queue_top.d_dist) {
				//cout << "           <= " << endl;
				result.push_back(*plist);
				plist = temp.erase(plist);
			}
			else plist++;
			if (K == result.size()) {
				input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
				input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
				assert(input_argument.pruning_power != INF);
				cout << "pruning power: " << input_argument.pruning_power << endl;
#ifdef _DEBUG
				input_argument.knn_total_time = TOOL::recordFinishTime(TOOL::time_record[13]);

				cout << "Total KNN time : " << input_argument.knn_total_time << " us" << endl;

				cout << "R-tree index navigate time : " << input_argument.navigate_index_time << " us" << endl;

				cout << "R-tree Euclidean distance time : " << input_argument.distance_euc_time << " us" << endl;

				cout << "R-tree index distance time : " << input_argument.distance_lowbound_time << " us" << endl;

				assert(input_argument.IO_cost == input_argument.sum_distance_euc);
				cout << "I/O cost: " << input_argument.IO_cost << endl;

				cout << "!!!!!!!!!!!!!!!!!!!!!!!!            APCA&PAA Find result !!!!!!!!!!!!!      result list size: " << result.size() << endl;

				result.sort([](const ORIGINAL_TIME_SERIES_PAIR& first, const  ORIGINAL_TIME_SERIES_PAIR& second) {return first.d_dist < second.d_dist; });//small to big

				for (typename list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
					cout << it->d_dist << ", " << it->original_time_series_id << "; ";
				}
#endif
				priority_queue<APCA_NODE_PAIR, vector<APCA_NODE_PAIR>, priorityIncrement >().swap(queue);
				temp.clear();
				//result.clear();
				//list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
				//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
				return true;
			}
		}

		queue.pop();

		if (m_temp_queue_top.p_rtree_node == nullptr) { //is APCA data point
			input_argument.navigate_index_time += TOOL::recordFinishTime(TOOL::time_record[4]);

			//cout << "    queue.top is data point\n";
			DataType* original_time_series = new DataType[g_time_series_length * arity_d];

			double run_time0;
			//recordStartTime(whole_first_time_start[3], whole_first_dqFreq[3]);
			recordStartTime(time_record[3]);

			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;
			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;
			//getFileStreamByID(file_name, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
			//TOOL::getMultiFoldToSingleByID(multi_file_name, arity_d, g_time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);

			if (input_argument.read_multiple_file_name) {
				TOOL::getMultiFoldToSingleByID(multi_file_name, input_argument.arity_d, input_argument.time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
			}
			else {
				//file_stream = ifstream(TOOL::getStringByID(TOOL::file_address, input_argument.file_id));
				TOOL::getFileStreamByID(TOOL::getStringByID(TOOL::file_address, input_argument.file_id), input_argument.time_series_length, m_temp_queue_top.original_time_series_id, original_time_series);
			}
			TOOL::normalizeStandard(input_argument.time_series_length, original_time_series);

			input_argument.IO_cost++;
			//g_n_account_apca_point++;
			input_argument.sum_distance_euc++;
#ifdef _DEBUG
			TOOL::recordStartTime(TOOL::time_record[6]);//time for euclidean distance
#endif
			tempOriginalTimeSeriesPair.d_dist = distanceEUC(g_query_time_series, g_time_series_length * arity_d, original_time_series, g_time_series_length * arity_d);

#ifdef _DEBUG
			input_argument.distance_euc_time += TOOL::recordFinishTime(TOOL::time_record[6]);
#endif
			delete[] original_time_series;
			original_time_series = nullptr;

			//cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
			temp.push_back(tempOriginalTimeSeriesPair);
		}
		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
			//cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;

			typename APCA_QUAL::APCA QProjection;
			APCA_QUAL::initialAPCA(QProjection, apcaRTree.NUMDIMS / 2);

			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {
				//g_n_account_leaf_node++;

				f_temp_APCA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);
				//cout << tempAPCAPair.id_originalTimeSeries << endl;
				f_temp_APCA_Pair.p_APCA_point = &APCALinkOriginal[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)].APCALink;

				//cout << "KNN : data point ID = " << tempAPCAPair.id_originalTimeSeries << endl;
				//cout << tempAPCAPair.APCAValue << endl;
#ifdef _DEBUG
				TOOL::recordStartTime(TOOL::time_record[5]);
#endif
				//f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point), *f_temp_APCA_Pair.p_APCA_point);
				f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length * arity_d, *f_temp_APCA_Pair.p_APCA_point, QProjection), *f_temp_APCA_Pair.p_APCA_point);
#ifdef _DEBUG
				input_argument.distance_lowbound_time += TOOL::recordFinishTime(TOOL::time_record[5]);
#endif
				//printf("\nrun_time = %f us\n", run_time0);

				//cout << tempAPCAPair.key;
				f_temp_APCA_Pair.p_rtree_node = nullptr;
				//	cout << tempAPCAPair.key;
#ifdef _DEBUG
				double queue_run_time;
				recordStartTime(time_record[7]);
#endif
				//cout << "            Push apca data. Branch id: " << i << " DLB: " << f_temp_APCA_Pair.d_dist << ", time series ID: " << f_temp_APCA_Pair.original_time_series_id << endl;
				queue.push(f_temp_APCA_Pair);
#ifdef _DEBUG
				//g_n_time_leaf_node_push_count += recordFinishTime(time_record[7], queue_run_time);
#endif
				//printf("\nrun_time = %f us\n", run_time0);
				//cout << "KNN : queue.size() = " << queue.size() << endl;
				//		//cout << tempAPCAPair.key;
			}
			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
			//	//cout << tempOriginalTimeSeriesPair.key;

			APCA_QUAL::deleteAPCA(QProjection);
		}
		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {															//is internal Node
			//cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
#ifdef _DEBUG
			double run_time2;
			recordStartTime(time_record[8]);
#endif
			APCA_NODE_PAIR tempApcaPair;
			REGION fs_region_G;
			initialREGION(fs_region_G, apcaRTree.NUMDIMS / 2);
			////cout << "regionNum = " << G.regionNum << endl;

			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
				//g_n_account_child_node++;
#ifdef _DEBUG
				double mindistQR_run_time;
				recordStartTime(time_record[9]);
				TOOL::recordStartTime(TOOL::time_record[5]);//distance PLA & MBR time
#endif
#ifdef _DEBUG
				assert(input_argument.point_multi_single_length == g_time_series_length * arity_d);
#endif
				//APCA paper version min id == max id
				tempApcaPair.d_dist = MINDISTQR(g_query_time_series, g_time_series_length * arity_d, getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
				//191119 origianl version Update MBR structure, MBR in Rtree and get region funciton  min id < max id
				//tempApcaPair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::get_region_G_original(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
#ifdef _DEBUG
				input_argument.distance_lowbound_time += TOOL::recordFinishTime(TOOL::time_record[5]);
#endif
				/*for (int i = 0; i < fs_region_G.regionNum; i++) {
				cout << "            G1[" << i << "] = " << fs_region_G.G1[i] << ", G2[" << i << "] = " << fs_region_G.G2[i] << ", G3[" << i << "] = " << fs_region_G.G3[i] << ", G4[" << i << "] = " << fs_region_G.G4[i] << endl;
				}*/
#ifdef _DEBUG
				//g_n_time_child_node_MINDIST_count += recordFinishTime(time_record[9], mindistQR_run_time);
#endif
				//printf("\nrun_time = %f us\n", run_time0);

				tempApcaPair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
#ifdef _DEBUG
				double push_run_time;
				recordStartTime(time_record[10]);
#endif
				//cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
				queue.push(tempApcaPair);
#ifdef _DEBUG
				//cout << "KNN : queue.size() = " << queue.size() << endl;
				//g_n_time_child_node_push_count += recordFinishTime(time_record[10], push_run_time);
				//printf("\nrun_time = %f us\n", run_time0);
#endif
			}

			deleteREGION(fs_region_G);
#ifdef _DEBUG
			//g_n_time_count2 += recordFinishTime(time_record[8], run_time2);
			//printf("\nrun_time = %f us\n", run_time0);
#endif
		}
		else {
			assert(0);
		}
	}
	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
	input_argument.pruning_power = 2;
#ifdef _DEBUG
	assert(input_argument.pruning_power != INF);
	cout << "??????????????????     APCA KNN failed    ??????????????????" << endl;
	cout << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
	/*ofstream outfile(input_argument.write_file_name[16] + ".txt", ios::app);
	assert(outfile.is_open());
	outfile << "??????????????????     APCA KNN failed    ??????????????????   " << PAA_or_APCA << endl;
	outfile << "n = " << input_argument.time_series_length << ", N = " << input_argument.point_dimension << ", number = " << input_argument.point_number << ", K = " << input_argument.K << ", MAXNODES = " << input_argument.rtree_max_nodes << endl;
	outfile << "K: " << input_argument.K << ", result.size: " << result.size() << endl;
	outfile.close();*/
	TOOL::printInputArgument(input_argument);

	if (PAA_or_APCA == 0) {
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[0], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[1], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[3], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[4], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[5], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[6], "NULL");
	}
	else if (PAA_or_APCA == 1) {
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[8], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[9], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[11], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[12], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[13], "NULL");
		TOOL::writeSingleResult(input_argument, input_argument.write_file_name[14], "NULL");
	}
	else {
		assert(0);
	}
#endif
	//assert(0);
	priority_queue<APCA_NODE_PAIR, vector<APCA_NODE_PAIR>, priorityIncrement >().swap(queue);
	list <ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
	//list <ORIGINAL_TIME_SERIES_PAIR>().swap(result);
	return false;
}

TEMPLATE
void APCA_KNN_QUAL::SimpleBaseKNNSearch(const DataType* g_query_time_series, const DataType& m_file_time_series_length, const DataType& mg_d_index_point_number, const int& K, const string& file_name, priority_queue<DataType, vector<DataType>, greater<DataType> >& const q_base_queue) {
	//printf("<///////**  Base KNN Begin  **////////>\n");
	DataType* original_time_series = new DataType[int(m_file_time_series_length)];
	//priority_queue<DataType, vector<DataType>, greater<DataType> > q_base_queue;
	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	int i = 0, j = NULL;
	while (!file_stream.eof() && i < mg_d_index_point_number) {
		file_stream >> fs_row_string;
		stringstream sstr(fs_row_string);
		j = -1;
		while (getline(sstr, fs_row_number, ',') && j < m_file_time_series_length) {
			if (j > -1) {
				original_time_series[j] = stod(fs_row_number);
			}
			j++;
		}
		q_base_queue.push(distanceEUC(g_query_time_series, m_file_time_series_length, original_time_series, m_file_time_series_length));
		i++;
	}

	file_stream.close();
	delete[] original_time_series;
	original_time_series = nullptr;

	cout << "Top K: ";
	for (i = 0; i < K; i++) {
		cout << q_base_queue.top() << " ";
		q_base_queue.pop();
	}
	cout << endl;
}

TEMPLATE
void APCA_KNN_QUAL::writeResult(INPUT_ARGUMENT& input_a, string write_file_name) {
	time_t now = time(0);// system time
	char* dt = ctime(&now);// from now to string
	ofstream outfile(write_file_name + ".txt", ios::app);
	outfile << endl << dt;
	//outfile << dimension_reduction_name[PAA_or_APCA] << " n = " << input_a.mg_file_time_series_length << ", N = " << input_a.mg_APCA_point_dimension << ", number = " << input_a.mg_d_index_point_number << ", K = " << input_a.K << ", MAXNODES = " << input_a.mg_max_nodes << endl;
	//outfile << "  count apca point = " << g_n_account_apca_point << " times, p= " << g_n_account_apca_point / input_a.mg_d_index_point_number << endl;
	//outfile << "  APCA Memory = " << memory_account[0] + memory_account[1] + memory_account[2] + memory_account[3] << "  RTree Memeory = " << memory_account[4] + memory_account[5] + memory_account[6] + memory_account[7] + memory_account[8] << "  KNN Memory = " << memory_account[9] << endl;
	outfile.close();
}

TEMPLATE
void APCA_KNN_QUAL::recordStartTime(TIME& time) {
	LARGE_INTEGER whole_first_f;    //timer frequency
	QueryPerformanceFrequency(&whole_first_f);
	time.dqFrequency = (double)whole_first_f.QuadPart;
	QueryPerformanceCounter(&time.time_start);
}

TEMPLATE
double& APCA_KNN_QUAL::recordFinishTime(TIME& time, double& whole_first_run_time) {
	QueryPerformanceCounter(&time.time_over);    //Finish recording time  5
	whole_first_run_time = 1000000 * (time.time_over.QuadPart - time.time_start.QuadPart) / time.dqFrequency;
	//cout << whole_first_run_time << endl;

	return whole_first_run_time;
}

TEMPLATE
template<typename T>
void APCA_KNN_QUAL::printArray(T*& const test_array, const int& array_length) {
	assert(test_array != nullptr);

	for (int i = 0; i < array_length; i++) {
		cout << i << ":" << test_array[i] << ", ";
		//assert(test_array[i+1]- test_array[i]==1);
	}
	cout << endl;
}

TEMPLATE
template<typename T>
void APCA_KNN_QUAL::writeSingleResult(const string& write_file_name, T& result) {
	ofstream outfile(write_file_name + ".txt", ios::app);
	//ofstream outfile(write_file_name + ".txt");
	assert(outfile.is_open());
	outfile << result << endl;
	outfile.close();
}

TEMPLATE
template<typename T>
void APCA_KNN_QUAL::normalizeA_B(const DataType& left_endpoint, const DataType right_endpoint, T*& const original_array, const int& array_length, T*& const normalized_array) {
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
double& APCA_KNN_QUAL::getAverage(T*& const original_array, const int& const array_length, double& average) {
	double sum = 0;

	for (int array_id = 0; array_id < array_length; array_id++) {
		sum += original_array[array_id];
	}

	average = sum / double(array_length);

	//cout << "average: " << average << endl;

	return average;
}

TEMPLATE
template<typename T>
double& APCA_KNN_QUAL::getVariance(T*& const original_array, const int& const array_length, double& variance) {
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

TEMPLATE//z-score normalization
template<typename T>
void APCA_KNN_QUAL::normalizeStandard(T*& const original_array, const int& const array_length, T*& const normalized_array) {
	double variance = NULL;
	double average = NULL;

	getAverage(original_array, array_length, average);

	getVariance(original_array, array_length, variance);

	for (int array_id = 0; array_id < array_length; array_id++) {
		normalized_array[array_id] = (original_array[array_id] - average) / variance;
	}

	/*printArray(normalized_array, array_length);*/
}

//TEMPLATE//PAA approximate original time series
//template<typename T>
//void APCA_KNN_QUAL::approximateOriginalFunctionPAA(const APCA& const italicC, const int& const n, const int& const N, T*& const approximation_PAA) {
//	int endOfLongSegment = input_argument.remainder*input_argument.segment_length_first;//begin of second part
//	int indexOfLongSegment = input_argument.remainder;
//	//double difference = NULL;
//	double sum = 0;
//	double point_id = NULL;
//	int count[1280];
//	int i = 0;
//	int* pointer = count;
//
//	for (int segment_id = 0; segment_id < input_argument.point_dimension; segment_id++) {
//		/*for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//		point_id = input_argument.segment_length_second*segment_id + interval_id;
//		difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//		deviation_max = difference > deviation_max ? difference : deviation_max;
//		sum += difference * difference;
//		}*/
//		//cout << "i: " << segment_id <<endl;
//		if (segment_id<input_argument.remainder) {
//			for (int interval_id = 0; interval_id < input_argument.segment_length_first; interval_id++) {
//				point_id = input_argument.segment_length_first*segment_id + interval_id;// time series id
//																						//cout << "    j: " << point_id << endl;
//				count[i] = point_id;
//				i++;
//				//difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximate_array[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				//cout << "difference: "<<difference<<endl;
//			}
//			//indexOfLongSegment = segment_id+1;
//			//endOfLongSegment = int((segment_id+1) * input_argument.segment_length_first);
//		}
//		else {
//			//cout <<"remainder: " <<input_argument.remainder*input_argument.segment_length_first << endl;
//			//cout<<"endOfLongSegment: "<<endOfLongSegment<<endl;
//			assert(input_argument.remainder*input_argument.segment_length_first == endOfLongSegment);
//			//assert(indexOfLongSegment== input_argument.remainder-1);
//			for (int interval_id = 0; interval_id < input_argument.segment_length_second; interval_id++) {
//				point_id = endOfLongSegment + (segment_id - indexOfLongSegment)*input_argument.segment_length_second + interval_id;// time series id
//																																   //cout << "    j: " << point_id << endl;
//				assert(point_id >= 0);
//				//cout << "id: " << point_id << endl;
//				count[i] = point_id;
//				if (i != 0) {
//					assert(count[i] - count[i - 1] == 1);
//				}
//				i++;
//				//difference = fabs(pla.a[segment_id] * point_id + pla.b[segment_id] - original_array[int(point_id)]);
//				approximate_array[int(point_id)] = pla.a[segment_id] * interval_id + pla.b[segment_id];
//				/*if (difference>50) {
//				cout << pla.a[segment_id] * point_id + pla.b[segment_id]<<" "<< original_array[int(point_id)] << endl;
//				}*/
//				//cout << "difference: " << difference << endl;
//
//			}
//
//		}
//
//

TEMPLATE
template<typename T>
double& APCA_KNN_QUAL::getReconstructionError(T*& const original_time_series, T*& const approximation_time_series, const int& const time_series_length, double& const deviation_sum, double& const deviation_max) {
	double sum = 0;
	deviation_max = 0;
	double difference = NULL;

	for (int array_id = 0; array_id < time_series_length; array_id++) {
		difference = fabs(original_time_series[array_id] - approximation_time_series[array_id]);
		//deviation_sum += difference;
		sum += difference * difference;
		deviation_max = deviation_max < difference ? difference : deviation_max;
	}

	deviation_sum = sqrt(sum);

	return deviation_sum;
}

TEMPLATE
void APCA_KNN_QUAL::getDeviationIterationPAA(DataType*& const original_time_series, const int& const array_length, string*& const write_file_name, const int& segment_begin, const int& segment_end, const int& segment_interval) {
	double deviation_sum = 0;
	double deviation_max = 0;

	DataType* original_normalization = new DataType[array_length];
	fill_n(original_normalization, array_length, NULL);
	/*DataType* PAA_approximation = new DataType[array_length];
	fill_n(PAA_approximation, array_length, NULL);
	DataType* PAA_normalization = new DataType[array_length];
	fill_n(PAA_normalization, array_length, NULL);*/

	typename APCA_QUAL::APCA paa;

	APCA_KNN_QUAL::normalizeStandard(original_time_series, array_length, original_normalization);

	//APCA_KNN_QUAL::printArray(original_time_series, array_length);
	cout << "original_normalization " << endl;
	APCA_KNN_QUAL::printArray(original_normalization, array_length);

	for (int segment_number = segment_begin; segment_number <= segment_end; segment_number += segment_interval) {
		cout << "segment_number: " << segment_number << endl;
		APCA_QUAL::initialAPCA(paa, segment_number);
		APCA_QUAL::divideRemainderPAA(original_normalization, paa, array_length, segment_number);

		cout << "PAA: ";
		APCA_KNN_QUAL::printArray(paa.v, segment_number);

		//approximateOriginalFunction(input_argument, original_normalization, chebyshev_share, PAA_approximation);

		APCA_KNN_QUAL::distanceAE(original_normalization, array_length, paa, deviation_sum, deviation_max);

		//APCA_KNN_QUAL::normalizeStandard(PAA_approximation, array_length, PAA_normalization);
		/*APCA_KNN_QUAL::getAverage(chebyshve_normalization, array_length,average);
		APCA_KNN_QUAL::getVariance(chebyshve_normalization, array_length,variance);*/
		/*cout << "normalized Chebyshev approximation: ";
		*/
		//APCA_KNN_QUAL::getReconstructionError(original_normalization, PAA_normalization, array_length, deviation_sum, deviation_max);

		writeSingleResult(input_argument, write_file_name[2], deviation_sum);
		writeSingleResult(input_argument, write_file_name[3], deviation_max);

		//writeApproximationResult(input_argument, PAA_normalization, deviation_sum, deviation_max);

		APCA_QUAL::deleteAPCA(paa);
	}

	delete[] original_normalization;
	original_normalization = nullptr;
	/*delete[] PAA_approximation;
	PAA_approximation = nullptr;
	delete[] PAA_normalization;
	PAA_normalization = nullptr;*/
}

TEMPLATE
void APCA_KNN_QUAL::getDeviationIterationAPCA(DataType*& const original_time_series, const int& const array_length, string*& const write_file_name, const int& segment_begin, const int& segment_end, const int& segment_interval) {
	double deviation_sum = 0;
	double deviation_max = 0;

	DataType* original_normalization = new DataType[array_length];
	fill_n(original_normalization, array_length, NULL);
	/*DataType* PAA_approximation = new DataType[array_length];
	fill_n(PAA_approximation, array_length, NULL);
	DataType* PAA_normalization = new DataType[array_length];
	fill_n(PAA_normalization, array_length, NULL);*/

	typename APCA_QUAL::APCA apca;

	//nromalize time series;
	APCA_KNN_QUAL::normalizeStandard(original_time_series, array_length, original_normalization);

	//APCA_KNN_QUAL::printArray(original_time_series, array_length);
	cout << "original_normalization " << endl;
	APCA_KNN_QUAL::printArray(original_normalization, array_length);

	for (int segment_number = segment_begin; segment_number <= segment_end; segment_number += segment_interval) {
		cout << "segment_number: " << segment_number << endl;
		APCA_QUAL::initialAPCA(apca, segment_number);
		//APCA_QUAL::divideRemainderPAA(original_normalization, apca, array_length, segment_number);
		APCA_QUAL::getAPCAPoint(original_normalization, array_length, segment_number, apca);

		cout << "APCA: ";
		APCA_KNN_QUAL::printArray(apca.v, segment_number);

		//approximateOriginalFunction(input_argument, original_normalization, chebyshev_share, PAA_approximation);

		APCA_KNN_QUAL::distanceAE(original_normalization, array_length, apca, deviation_sum, deviation_max);

		//APCA_KNN_QUAL::normalizeStandard(PAA_approximation, array_length, PAA_normalization);
		/*APCA_KNN_QUAL::getAverage(chebyshve_normalization, array_length,average);
		APCA_KNN_QUAL::getVariance(chebyshve_normalization, array_length,variance);*/
		/*cout << "normalized Chebyshev approximation: ";
		*/
		//APCA_KNN_QUAL::getReconstructionError(original_normalization, PAA_normalization, array_length, deviation_sum, deviation_max);

		writeSingleResult(input_argument, write_file_name[4], deviation_sum);
		writeSingleResult(input_argument, write_file_name[5], deviation_max);

		//writeApproximationResult(input_argument, PAA_normalization, deviation_sum, deviation_max);

		APCA_QUAL::deleteAPCA(apca);
	}

	delete[] original_normalization;
	original_normalization = nullptr;
	/*delete[] PAA_approximation;
	PAA_approximation = nullptr;
	delete[] PAA_normalization;
	PAA_normalization = nullptr;*/
}

TEMPLATE
void APCA_KNN_QUAL::compareDiffIteration(DataType*& const original_time_series, const int& const array_length, string*& const write_file_name, const int& segment_begin, const int& segment_end, const int& segment_interval) {
	//PLA_QUAL pla(array_length, point_dimension, 5, 5, 3, "test.txt");
	assert(0);
	//double deviation_sum = 0;
	//double deviation_max = 0;

	//typename APCA_QUAL::APCA apca;
	//typename APCA_QUAL::APCA paa;
	////typename PLA_QUAL::PLA pla;

	//CHEBYSHEV_QUAL chebyshev(array_length, NULL, NULL, NULL, NULL, 1, "", "");
	///*double average = NULL;
	//double variance = NULL;*/

	//DataType* original_normalization = new DataType[array_length];
	//fill_n(original_normalization, array_length, NULL);
	//DataType* chebyshev_approximation = new DataType[array_length];//chebyshev
	//fill_n(chebyshev_approximation, array_length, NULL);
	//DataType* chebyshev_normalization = new DataType[array_length];//chebyshev
	//fill_n(chebyshev_normalization, array_length, NULL);

	//normalizeStandard(original_time_series, array_length, original_normalization);// normalize original time series

	//for (int segment_num = segment_begin; segment_num <= segment_end; segment_num += segment_interval) {
	//	chebyshev.input_argument.degree_m = segment_num - 1;//n=m+1;
	//	chebyshev.approximateOriginalFunction(chebyshev.input_argument, original_normalization, chebyshev.chebyshev_share, chebyshev_approximation);

	//	APCA_KNN_QUAL::normalizeStandard(chebyshev_approximation, array_length, chebyshev_normalization);
	//	APCA_KNN_QUAL::getReconstructionError(original_normalization, chebyshev_normalization, array_length, deviation_sum, deviation_max);

	//	chebyshev.writeSingleResult(input_argument, write_file_name[0], deviation_sum);
	//	chebyshev.writeSingleResult(input_argument, write_file_name[1], deviation_max);

	//	chebyshev.input_argument.write_file_name = write_file_name[2];
	//	chebyshev.writeApproximationResult(chebyshev.input_argument, chebyshev_normalization, deviation_sum, deviation_max);
	//}

	//delete[] original_normalization;
	//original_normalization = nullptr;
	//delete[] chebyshev_approximation;
	//chebyshev_approximation = nullptr;
	//delete[] chebyshev_normalization;
	//chebyshev_normalization = nullptr;
}

#endif

