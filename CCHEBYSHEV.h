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

	void rangeSearch(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue);

	void rangeSearchMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue);

	void rangeSearchMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const typename TOOL::DATA_SOURCE& const data_source, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree, const double& const radius, priority_queue <ID_DIST, vector<ID_DIST>, priorityDecrement >& queue);

	void KNNSearch(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree, priority_queue <ID_DIST, vector<ID_DIST>, priorityIncrement >& queue_result);

	std::multiset<pair<double, int>> KNNCHEBYMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree);

	//191223 add data_source
	std::multiset<pair<double, int>> KNNCHEBYMulti(typename TOOL::INPUT_ARGUMENT& const input_argument, const typename TOOL::DATA_SOURCE& const data_source, DataType*& const query_time_series, const CHEBYSHEV_SHARE& const  chebyshev_share, CHEBYSHEV*& const chebyshev_array, RTREE& const RTree);

	void writeApproximationResult(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const result_time_series, const double& deviation_sum, const double& deviation_max);

	template<typename T>
	void writeSingleResult(const string& write_file_name, T& result);

	void getSeveralRecronstructionError();

	//friend class APCA_KNN_QUAL;
	//friend class CAPCA_KNN<DataType>;


private:
};

#include "CCHEBYSHEV.cpp"
#endif

