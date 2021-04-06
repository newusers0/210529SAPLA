#pragma once
#ifndef MULTI_DIMENSION_H
#define MULTI_DIMENSION_H

#include "pch.h"
#include "CAPCA.h"
#include "CAPCA_KNN.h"
#include "CPLA.h"
#include "CCHEBYSHEV.h"
#include "CAPLA.h"
#include "APLA_ICDE07.h"
#include "./lib/saxquantizer.hpp"//200923

TEMPLATE
class MULTI_DIMENSION : public RTREE, virtual public TOOL, public APCA_QUAL, virtual public APCA_KNN_QUAL, public PLA_QUAL, public CHEBYSHEV_QUAL, public APLA, public APLA_ICDE07<DataType> {
public:
	struct MULTI_DIM;
	template<typename T>
	struct APLA_NODE_PAIR;
	struct RTREE_NODE_PAIR;
	template<typename T>
	struct priorityIncrement;
	MULTI_DIM multi_dim;
	struct TOOL::DATA_SOURCE data_source;
	struct TOOL::INPUT_ARGUMENT input_argument;
	struct TOOL::OUTPUT_ARGUMENT output_argument;
	struct TOOL::RESULT_RECORD result_record;
	RTREE RTree;
	vector<RTREE::Rect> rtree_rectangle_vector;
	struct APCA_QUAL::APCA* apca_MBR = nullptr;
	struct CHEBYSHEV_QUAL::CHEBYSHEV_SHARE chebyshev_share;
	struct APCA_QUAL::APCA_ARRAY* apca_array = nullptr;
	struct APCA_QUAL::APCA* apca_point_array = nullptr;
	struct PLA_QUAL::PLA* pla_array = nullptr;
	struct CHEBYSHEV_QUAL::CHEBYSHEV* chebyshev_array = nullptr;
	std::vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>> apla_array_vector;
	std::vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>> icde_array_vector;
	MULTI_DIMENSION(const int& n, const int& N, const int& d, string*& const multi_file_name, const int& representation_option);
	MULTI_DIMENSION(const int& const n, const int& const N, const int& const file_id, const int& const d, const int& const point_number, const int& const query_time_series_id, const int& const rtree_max_nodes, const int& const K, const int& const representation_option, string*& const multi_file_name, string*& const write_file_name);
	MULTI_DIMENSION(const typename TOOL::DATA_SOURCE& const data_source, const int& const n, const int& const N, const double& const initial_N, const int& const file_id, const bool& const change_file, const int& const d, const int& const point_number, const int& const query_time_series_id, const int& const rtree_max_nodes, const int& const K, const int& const representation_option, string*& const multi_file_name, string*& const write_file_name);
	~MULTI_DIMENSION();
	void initialMultiDimension();
	void deleteMultiDimension();
	template<typename T>
	void get_cmin_cmax(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, typename APCA_QUAL::APCA& const CminParameter, typename APCA_QUAL::APCA& const CmaxParameter);
	template<typename T, typename Y>
	void get_cmin_cmax_original(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle);
	template<typename T, typename Y>
	void get_cmin_cmax_apca(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle);
	template<typename T, typename Y, typename U>
	void get_cmin_cmax_apca(const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const linked_list, U& const rtree_rectangle);
	template<typename T, typename Y>
	void get_cmin_cmax_pla(const int& const time_series_length, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle);
	template<typename T>
	double& getPLAMBRSegmentDistance(const typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const int& segment_id, const T& pla_array_qeury, double& pla_MBR_segment_distance);
	template<typename T>
	double& getPLAMBRDistance(typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const T& pla_array_qeury, double& pla_MBR_distance);
	template<typename T, typename Y>
	void all_approximation_build_rtree(vector<T>& const multi_y_projection_argument, vector<DoublyLinkedList<Y>>& const multi_all_linked_list, vector<DoublyLinkedList<Y>>& const multi_cluster_linked_list);
	template<typename T>
	T compute_knn_accuracy(const int& const K, const multiset<pair<T, int>>& const squential_scan_result_set, const multiset<pair<T, int>>& const knn_result_set);
	template<typename T, typename Y, typename U>
	RTREE& all_knn(typename TOOL::DATA_SOURCE& const data_source, const vector<T>& const query_time_series_vector, T*& query_time_series, const vector<Y>& const multi_y_projection_argument, const std::multiset<pair<U, int>>& const squential_scan_result_set);
	template<typename T, typename Y>
	std::multiset<pair<double, int>>  all_knn_multi(typename TOOL::DATA_SOURCE& const data_source, typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<Y>& const multi_y_projection_argument, DataType*& g_query_time_series, const RTREE& const apcaRTree, T& const approximation_array);
};

#include "MULTI_DIMENSION.cpp"

#endif


