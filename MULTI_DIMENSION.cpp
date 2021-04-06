#pragma once
#ifndef _MULTI_DIMENSION_CPP_
#define _MULTI_DIMENSION_CPP_
#include "MULTI_DIMENSION.h"

TEMPLATE
struct MULTI::MULTI_DIM {
	DataType** multi_dimension;
	DataType* projected_dimension;
};
TEMPLATE
template<typename T>
struct MULTI::APLA_NODE_PAIR {
	double d_dist = INF;
	int original_time_series_id = INF;
	DoublyLinkedList<T>* approximation_pointer = nullptr;
	RTREE::Node* p_rtree_node = nullptr;
};

TEMPLATE
struct MULTI::RTREE_NODE_PAIR {
	double d_dist = INF;
	int original_time_series_id = INF;
	RTREE::Node* p_rtree_node = nullptr;
};
TEMPLATE
template<typename T>
struct MULTI::priorityIncrement {
	bool operator ()(const T& const a, const T& const b) {
		return a.d_dist > b.d_dist;
	}
};

TEMPLATE
MULTI::MULTI_DIMENSION(const int& n, const int& N, const int& d, string*& const multi_file_name, const int& representation_option) {
	input_argument.remainder = int(n) % int(N);
	double integerDividend = n - input_argument.remainder;
	input_argument.segment_length_second = integerDividend / N;
	input_argument.segment_length_first = input_argument.segment_length_second + 1;
	assert(input_argument.segment_length_second > 1);

	input_argument.time_series_length = n;
	input_argument.point_dimension = N;
	input_argument.degree_m = N - 1;
	input_argument.arity_d = d;
	input_argument.point_multi_single_dimension = d * N;
	input_argument.read_multiple_file_name = multi_file_name;
	input_argument.representation_option = representation_option;
}

TEMPLATE
MULTI::MULTI_DIMENSION(const int& const n, const int& const N, const int& const file_id, const int& const d, const int& const point_number, const int& const query_time_series_id, const int& const rtree_max_nodes, const int& const K, const int& const representation_option, string*& const multi_file_name, string*& const write_file_name) {
	input_argument.representation_option = representation_option;
	switch (input_argument.representation_option) {
	case 1: 
	case 6: {
		input_argument.point_dimension = N / 3;
		break;
	}
	case 2: 
	case 3: {
		input_argument.point_dimension = N / 2;
		break;
	}
	case 4: 
	case 7: 
	case 5: {
		input_argument.point_dimension = N;
		break;
	}
	default:
		assert(0);
		break;
	}
	input_argument.time_series_length = n;
	input_argument.degree_m = input_argument.point_dimension - 1;
	input_argument.file_id = file_id;
	input_argument.arity_d = d;
	input_argument.point_multi_single_dimension = d * input_argument.point_dimension;
	input_argument.point_multi_single_length = d * n;
	input_argument.remainder = int(n) % int(input_argument.point_dimension);
	input_argument.segment_length_second = (n - input_argument.remainder) / input_argument.point_dimension;
	input_argument.segment_length_first = input_argument.segment_length_second + 1;
	input_argument.query_time_series_id = query_time_series_id;
	input_argument.point_number = point_number;
	input_argument.rtree_max_nodes = rtree_max_nodes;
	input_argument.K = K;
	input_argument.read_multiple_file_name = multi_file_name;
	input_argument.write_file_name = write_file_name;
	initialMultiDimension();
	input_argument.pruning_power = 0.0;
	input_argument.sum_distance_euc = 0.0;
	input_argument.build_rtree_time = 0.0;
	input_argument.approximation_query_time = 0.0;
	input_argument.knn_rest_part_time = 0.0;
	input_argument.knn_total_time = 0.0;
	input_argument.navigate_index_time = 0;
	input_argument.distance_lowbound_time = 0; 
	input_argument.distance_euc_time = 0;
	input_argument.IO_cost = 0.0;
	input_argument.result_accuracy = 0;
}

TEMPLATE
MULTI::MULTI_DIMENSION(const typename TOOL::DATA_SOURCE& const data_source, const int& const n, const int& const N, const double& const initial_N, const int& const file_id, const bool& const change_file, const int& const d, const int& const point_number, const int& const query_time_series_id, const int& const rtree_max_nodes, const int& const K, const int& const representation_option, string*& const multi_file_name, string*& const write_file_name) {
	
	input_argument.change_file = change_file;
	this->data_source = data_source;
	input_argument.representation_option = representation_option;
	input_argument.time_series_length = n;
	input_argument.file_id = file_id;
	input_argument.arity_d = d;
	switch (input_argument.representation_option) {
	case 1: 
	case 6: 
	case 8: {
		if (initial_N != INF)
			input_argument.initial_N = initial_N;
		input_argument.point_dimension = N / 3;
		break;
	}
	case 2: 
	case 3: {
		input_argument.point_dimension = N / 2;
		break;
	}
	case 4:
	case 7: 
	case 5: 
	case 9: {
		input_argument.point_dimension = N;
		break;
	}
	default:
		assert(0);
		break;
	}
	input_argument.point_multi_single_dimension = input_argument.point_dimension;
	input_argument.point_multi_single_length = this->data_source.multi_single_time_series_length;
	input_argument.time_series_length = input_argument.point_multi_single_length;
	input_argument.point_dimension = input_argument.point_multi_single_dimension;
	input_argument.degree_m = input_argument.point_dimension - 1;
	input_argument.remainder = int(input_argument.point_multi_single_length) % int(input_argument.point_dimension);
	input_argument.segment_length_second = (input_argument.point_multi_single_length - input_argument.remainder) / input_argument.point_dimension;
	input_argument.segment_length_first = input_argument.segment_length_second + 1;
	input_argument.query_time_series_id = query_time_series_id;
	input_argument.point_number = data_source.point_number;
	input_argument.rtree_max_nodes = rtree_max_nodes;
	input_argument.K = K;
	input_argument.read_multiple_file_name = multi_file_name;
	input_argument.write_file_name = write_file_name;
	switch (data_source.data_type) {
	case 0:
		break;
	case 1:
		input_argument.read_file_name = data_source.read_file_address_vector.front();
		assert(input_argument.read_file_name == data_source.read_file_address_vector.front() && d == 1);
		break;
	case 2:
		input_argument.read_file_name = data_source.read_file_address;
		assert(input_argument.read_file_name == data_source.read_file_address_vector.front() && d == 1);

		break;
	case 3:
		assert(data_source.read_file_address_vector.size() == data_source.time_series_dimension && d > 1);
		break;
	case 4:
		assert(0);
		break;
	case 5:
		assert(0);
		break;
	default:
		assert(0);
		break;
	}
	initialMultiDimension();
	input_argument.representation_time = 0;
	input_argument.build_rtree_time = 0;
	input_argument.knn_total_time = 0;
	input_argument.whole_run_time = 0;
	input_argument.IO_cost = 0;
	input_argument.pruning_power = 0;
	input_argument.result_accuracy = 0;
	input_argument.sum_distance_euc = 0.0;
	input_argument.build_rtree_time = 0.0;
	input_argument.approximation_query_time = 0.0;
	input_argument.knn_rest_part_time = 0.0;
	input_argument.navigate_index_time = 0;
	input_argument.distance_lowbound_time = 0; 
	input_argument.distance_euc_time = 0;
}

TEMPLATE
MULTI::~MULTI_DIMENSION() {
	assert(input_argument.arity_d > 0);
	input_argument.time_series_length = INF;
	input_argument.degree_m = INF;
	input_argument.arity_d = INF;
	deleteMultiDimension();
	RTree.RemoveAll();
}

TEMPLATE
void MULTI::initialMultiDimension() {
	input_argument.representation_time = 0;
	input_argument.build_rtree_time = 0;
	input_argument.knn_total_time = 0;
	input_argument.whole_run_time = 0;
	input_argument.IO_cost = 0;
	input_argument.pruning_power = 0;
	input_argument.result_accuracy = 0;
	output_argument.sum_deviation = 0;

	switch (input_argument.representation_option) {
	case 1:
	case 6: 
	case 8: {
		RTree.initialRTree(input_argument.point_multi_single_dimension * 2, input_argument.rtree_max_nodes, input_argument.representation_option);
		APCA_QUAL::initial_rect_vector(input_argument.point_number, input_argument.point_multi_single_dimension * 2, rtree_rectangle_vector);
		apla_array_vector.resize(input_argument.point_number, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>());
		break;
	}
	case 2: {
		RTree.initialRTree(input_argument.point_multi_single_dimension * 2, input_argument.rtree_max_nodes, input_argument.representation_option);
		PLA_QUAL::initialPLAArray(input_argument.point_number, input_argument.point_multi_single_dimension, pla_array);
		break;
	}
	case 3:
	case 7:
	case 9:
	case 4: {
		RTree.initialRTree(input_argument.point_multi_single_dimension * 2, input_argument.rtree_max_nodes, input_argument.representation_option);
		APCA_QUAL::initialAPCAArray(input_argument.point_number, input_argument.point_multi_single_dimension, apca_point_array);
		break;
	}
	case 5: {
		RTree.initialRTree(input_argument.point_multi_single_dimension, input_argument.rtree_max_nodes, input_argument.representation_option);
		CHEBYSHEV_QUAL::initialCHEBYSHEVArray(input_argument.point_number, input_argument.time_series_length, input_argument.point_multi_single_dimension, chebyshev_array);
		CHEBYSHEV_QUAL::initialCHEBYSHEV_SHARE(input_argument, chebyshev_share);
		CHEBYSHEV_QUAL::getCHEBYSHEV_SHARE(input_argument, chebyshev_share);
		break;
	}
	default: 
		assert(0);
		break;
	}
	TOOL::printInputArgument(input_argument);
}

TEMPLATE
void MULTI::deleteMultiDimension() {

	switch (input_argument.representation_option) {
	case 1:
	case 6: {
		APCA_QUAL::delete_rect_vector(rtree_rectangle_vector);
		for (auto&& au : apla_array_vector) {
			au.clear();
		}
		apla_array_vector.clear();
		apla_array_vector.shrink_to_fit();
		APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_MBR);
		break;
	}
	case 8: {
		APCA_QUAL::delete_rect_vector(rtree_rectangle_vector);
		for (auto&& au : apla_array_vector) {
			for (auto auu : au) {
				if (auu.right_subsegment != nullptr) {
					delete auu.right_subsegment;
					auu.right_subsegment = nullptr;
				}
			}
			au.clear();
		}
		apla_array_vector.clear();
		apla_array_vector.shrink_to_fit();
		APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_MBR);
		break;
	}
	case 2: {
		PLA_QUAL::deletePLAArray(input_argument, pla_array);
		break;
	}
	case 3:
	case 7:
	case 9:
	case 4: {
		APCA_QUAL::deleteAPCAArray(input_argument.point_number, apca_point_array);
		break;
	}
	case 5: {
		CHEBYSHEV_QUAL::deleteCHEBYSHEV_SHARE(chebyshev_share);
		CHEBYSHEV_QUAL::deleteCHEBYSHEVArray(input_argument.point_number, chebyshev_array);
		break;
	}
	default:
		assert(0);
		break;
	}
}

TEMPLATE
template<typename T>
void MULTI::get_cmin_cmax(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, typename APCA_QUAL::APCA& const CminParameter, typename APCA_QUAL::APCA& const CmaxParameter) {

	for (int segment_id = 0; segment_id < linked_list.size(); segment_id++) {
		auto& const segment = linked_list[segment_id];

		CmaxParameter.r[segment_id] = segment.right_endpoint;
		CminParameter.r[segment_id] = segment.right_endpoint - segment.rectangle_width;
		CminParameter.v[segment_id] = segment.min_point.value;
		CmaxParameter.v[segment_id] = segment.max_point.value;
		
		if (dimension_id > 0) {
			int extend_right_endpoint = time_series_length * dimension_id;
			segment.right_endpoint += extend_right_endpoint;
			CminParameter.r[segment_id] += extend_right_endpoint;
			CmaxParameter.r[segment_id] = segment.right_endpoint;
		}
	}
}

TEMPLATE
template<typename T, typename Y>
void MULTI::get_cmin_cmax_original(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle) {

	rtree_rectangle.m_min[0] = 0; 
	rtree_rectangle.m_min[1] = linked_list.front().min_point.value;
	rtree_rectangle.m_max[0] = linked_list.front().right_endpoint; 
	rtree_rectangle.m_max[1] = linked_list.front().max_point.value;
	for (int segment_id = 1; segment_id < linked_list.size(); segment_id++) {
		auto& const segment = linked_list[segment_id];

		rtree_rectangle.m_min[segment_id * 2] = segment.right_endpoint - segment.rectangle_width + 1;
		rtree_rectangle.m_min[segment_id * 2 + 1] = segment.min_point.value;
		rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint; 
		rtree_rectangle.m_max[segment_id * 2 + 1] = segment.max_point.value;
		if (dimension_id > 0) {
			int extend_right_endpoint = time_series_length * dimension_id;
			segment.right_endpoint += extend_right_endpoint;
			rtree_rectangle.m_min[segment_id * 2] += extend_right_endpoint;
			rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;
		}
	}
}

TEMPLATE
template<typename T, typename Y>
void MULTI::get_cmin_cmax_apca(const int& const time_series_length, const int& const dimension_id, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle) {



	for (int segment_id = 0; segment_id < linked_list.size(); segment_id++) {
		auto& const segment = linked_list[segment_id];

		rtree_rectangle.m_min[segment_id * 2] = rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;
		rtree_rectangle.m_min[segment_id * 2 + 1] = segment.min_point.value;
		rtree_rectangle.m_max[segment_id * 2 + 1] = segment.max_point.value;

		if (dimension_id > 0) {
			int extend_right_endpoint = time_series_length * dimension_id;
			segment.right_endpoint += extend_right_endpoint;
			rtree_rectangle.m_min[segment_id * 2] += extend_right_endpoint;
			rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;
		}
	}
}


TEMPLATE
template<typename T, typename Y, typename U>
void MULTI::get_cmin_cmax_apca(const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const linked_list, U& const rtree_rectangle) {

	for (int segment_id = 0; segment_id < linked_list.size(); segment_id++) {
		const auto& const segment = linked_list[segment_id];
		const int id_coefficient = segment.right_endpoint + 1;
		const auto& const value_minmax = minmax_element(original_time_series_vector.begin() + int(id_coefficient - segment.rectangle_width), original_time_series_vector.begin() + id_coefficient);
		rtree_rectangle.m_min[segment_id * 2] = rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;
		rtree_rectangle.m_min[segment_id * 2 + 1] = *value_minmax.first;
		rtree_rectangle.m_max[segment_id * 2 + 1] = *value_minmax.second;
	}
}

TEMPLATE
template<typename T, typename Y>
void MULTI::get_cmin_cmax_pla(const int& const time_series_length, DoublyLinkedList<T>& const linked_list, Y& const rtree_rectangle) {
	for (int segment_id = 0; segment_id < linked_list.size(); segment_id++) {
		auto& const segment = linked_list[segment_id];
		auto result = minmax(segment.apla.b, segment.apla.a * (segment.rectangle_width - 1) + segment.apla.b);
		rtree_rectangle.m_min[segment_id * 2] = rtree_rectangle.m_max[segment_id * 2] = segment.right_endpoint;
		rtree_rectangle.m_min[segment_id * 2 + 1] = result.first;
		rtree_rectangle.m_max[segment_id * 2 + 1] = result.second;
	}
}

TEMPLATE
template<typename T>
double& MULTI::getPLAMBRSegmentDistance(const typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const int& segment_id, const T& pla_array_qeury, double& pla_MBR_segment_distance) {
	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);
	typename PLA_QUAL::POINT point;
	double sqrt_segment_legngth = sqrt(double(input_argument.segment_length_second));
	double perpendicular_b = (double(input_argument.segment_length_second) - 1.0) / 2.0;
	double argument_u_A = sqrt_segment_legngth * perpendicular_b;
	double perpendicular_a = sqrt((input_argument.segment_length_second * input_argument.segment_length_second - 1.0) / 12.0);
	double argument_v_B = perpendicular_a / perpendicular_b;
	point.u_A_1 = argument_u_A * (MBR.m_min[segment_id * 2] - pla_array_qeury.a[segment_id]);
	point.u_A_2 = argument_u_A * (MBR.m_max[segment_id * 2] - pla_array_qeury.a[segment_id]);
	point.v_A_1 = argument_v_B * point.u_A_1;
	point.v_A_2 = argument_v_B * point.u_A_2;
	point.u_B_1 = sqrt_segment_legngth * (pla_array_qeury.b[segment_id] - MBR.m_max[segment_id * 2 + 1]);
	point.u_B_2 = sqrt_segment_legngth * (pla_array_qeury.b[segment_id] - MBR.m_min[segment_id * 2 + 1]);

	double perpendicular_denominator = perpendicular_a * perpendicular_a + perpendicular_b * perpendicular_b;
	point.u_B_C_1 = -perpendicular_b * -perpendicular_b * point.u_B_1 / perpendicular_denominator; 
	point.v_B_C_1 = perpendicular_a * perpendicular_b * point.u_B_1 / perpendicular_denominator; 
	point.u_B_C_2 = -perpendicular_b * -perpendicular_b * point.u_B_2 / perpendicular_denominator; 
	point.v_B_C_2 = perpendicular_a * perpendicular_b * point.u_B_2 / perpendicular_denominator; 

	if (point.u_A_2 < 0) {
		if (point.u_B_2 < point.u_A_2) {
			return pla_MBR_segment_distance = PLA_QUAL::getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_2, point.v_B_2);
		}
		else if (point.u_B_2 >= point.u_A_2) {
			return pla_MBR_segment_distance = PLA_QUAL::getNearestDistance(point.u_B_1, point.v_B_1, point.u_B_2, point.v_B_2, point.u_A_2, point.v_A_2);
		}
		else assert(0);
	}
	else if (point.u_A_1 <= 0 && point.u_A_2 >= 0) {
		if (point.u_B_1 > 0) {
			if (point.u_B_C_1 >= point.u_A_2) {
				
				return pla_MBR_segment_distance = PLA_QUAL::getPointDistanceSquare(point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);
			}
			else if (point.u_B_C_1 < point.u_A_2) {
				
				return pla_MBR_segment_distance = PLA_QUAL::getPointDistanceSquare(point.u_B_1, point.v_B_1, point.u_B_C_1, point.v_B_C_1);
			}
			else assert(0);
		}
		else if (point.u_B_1 <= 0 && point.u_B_2 >= 0) {
			
			return pla_MBR_segment_distance = 0;
		}
		else if (point.u_B_2 < 0) {
			if (point.u_B_C_2 > point.u_A_1) {
				
				return pla_MBR_segment_distance = PLA_QUAL::getPointDistanceSquare(point.u_B_2, point.v_B_2, point.u_B_C_2, point.v_B_C_2);
			}
			else if (point.u_B_C_2 <= point.u_A_1) {
				return pla_MBR_segment_distance = PLA_QUAL::getPointDistanceSquare(point.u_A_1, point.v_A_1, point.u_B_2, point.v_B_2);
			}
			else assert(0);
		}
		else assert(0);
	}
	else if (point.u_A_1 > 0) {
		if (point.u_B_1 > point.u_A_1) {
			return pla_MBR_segment_distance = PLA_QUAL::getNearestDistance(point.u_A_1, point.v_A_1, point.u_A_2, point.v_A_2, point.u_B_1, point.v_B_1);
		}
		else if (point.u_B_1 <= point.u_A_1) {
			return pla_MBR_segment_distance = PLA_QUAL::getNearestDistance(point.u_B_1, point.v_B_1, point.u_B_2, point.v_B_2, point.u_A_1, point.v_A_1);
		}
		else assert(0);
	}
	else assert(0);

}

TEMPLATE
template<typename T>
double& MULTI::getPLAMBRDistance(typename TOOL::INPUT_ARGUMENT& const input_argument, const RTREE::Rect& MBR, const T& pla_array_qeury, double& pla_MBR_distance) {

	assert(input_argument.point_dimension != INF && input_argument.segment_length_second != INF && input_argument.segment_length_first != INF);

	double sum = 0;
	
	for (int i = 0; i < input_argument.point_dimension; i++) {
		sum += getPLAMBRSegmentDistance(input_argument, MBR, i, pla_array_qeury, pla_MBR_distance);
	}
	pla_MBR_distance = sqrt(sum);
	return pla_MBR_distance;
}

TEMPLATE
template<typename T, typename Y>
void MULTI::all_approximation_build_rtree(vector<T>& const multi_y_projection_argument, vector<DoublyLinkedList<Y>>& const multi_all_linked_list, vector<DoublyLinkedList<Y>>& const multi_cluster_linked_list) {
	input_argument.representation_time = 0;
	input_argument.build_rtree_time = 0;
	input_argument.whole_run_time = 0;
	output_argument.sum_deviation = 0;
	output_argument.max_deviation = 0;
	output_argument.max_deviation_multiple_width = 0;
	DataType* original_time_series = new DataType[input_argument.time_series_length];
	vector<DataType> normalized_series_vector(input_argument.time_series_length, INF);
	string row_string;
	string row_number;
	for (int point_id = 0; point_id < input_argument.point_number; point_id++) {
		fill_n(original_time_series, input_argument.time_series_length, INF);
		TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, point_id, normalized_series_vector);
		copy_n(normalized_series_vector.begin(), normalized_series_vector.size(), original_time_series);
		TOOL::recordStartTime(TOOL::time_record[14]);
		TOOL::recordStartTime(TOOL::time_record[0]);
		if (input_argument.representation_option == 1 || input_argument.representation_option == 6 || input_argument.representation_option == 8) {
			typename TOOL::RECTANGLE rectangle_insertion(input_argument.point_dimension << 1);
			switch (input_argument.representation_option) {
			case 1: {
				APLA::initial_SAPLA_200706(input_argument, normalized_series_vector, apla_array_vector[point_id], output_argument);
				break;
			}
			case 6: {
				APLA_ICDE07<DataType>::getAPLA_ICDE07(input_argument, normalized_series_vector, apla_array_vector[point_id]);
				break;
			}
			case 8: {
				APLA::initial_SAPLA_200706(input_argument, normalized_series_vector, apla_array_vector[point_id], output_argument);
				break;
			}
			default:
				assert(0);
				break;
			}
			get_cmin_cmax_apca(normalized_series_vector, apla_array_vector[point_id], rectangle_insertion);
			result_record.representation_time = TOOL::recordFinishTime(TOOL::time_record[0]);
			input_argument.representation_time += result_record.representation_time;
			TOOL::recordStartTime(TOOL::time_record[1]);
			RTree.Insert(rectangle_insertion.m_min, rectangle_insertion.m_max, point_id);
			input_argument.build_rtree_time += TOOL::recordFinishTime(TOOL::time_record[1]);
			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
			APLA::get_sum_deviation_no_ab(normalized_series_vector, apla_array_vector[point_id], result_record);
			output_argument.max_deviation += result_record.max_deviation;
			output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
			output_argument.sum_deviation += result_record.sum_deviation;
			TOOL::delete_rectangle(rectangle_insertion);
		}
		else if (input_argument.representation_option == 2) {
			typename PLA_QUAL::PLA pla_MBR; 
			PLA_QUAL::initialPLA(pla_MBR, input_argument.point_multi_single_dimension * 2);
			PLA_QUAL::getPLA(input_argument, original_time_series, pla_array[point_id]);
			PLA_QUAL::getPLAMBR(pla_array[point_id], pla_MBR);
			result_record.representation_time = TOOL::recordFinishTime(TOOL::time_record[0]);
			input_argument.representation_time += result_record.representation_time;
			TOOL::recordStartTime(TOOL::time_record[1]);
			RTree.Insert(pla_MBR.a, pla_MBR.b, point_id);
			input_argument.build_rtree_time += TOOL::recordFinishTime(TOOL::time_record[1]);
			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
			PLA_QUAL::get_pla_sum_deviation(input_argument, normalized_series_vector, pla_array[point_id], result_record);
			output_argument.max_deviation += result_record.max_deviation;
			output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
			output_argument.sum_deviation += result_record.sum_deviation;
			PLA_QUAL::deletePLA(pla_MBR);
		}
		else if (input_argument.representation_option == 3 || input_argument.representation_option == 4 || input_argument.representation_option == 7 || input_argument.representation_option == 9) {
			typename TOOL::RECTANGLE rectangle_insertion(input_argument.point_dimension << 1);
			switch (input_argument.representation_option) {
			case 3: {
				APCA_QUAL::getAPCAPoint(original_time_series, input_argument.time_series_length, input_argument.point_dimension, apca_point_array[point_id]);
				break;
			case 4:
				APCA_QUAL::divideRemainderPAA(original_time_series, apca_point_array[point_id], input_argument.time_series_length, input_argument.point_dimension);
				break;
			case 7:
				APCA_QUAL::divideRemainderPAA(original_time_series, apca_point_array[point_id], input_argument.time_series_length, input_argument.point_dimension);
				APCA_QUAL::get_paa_lagrangian(apca_point_array[point_id]);
				break;
			case 9:
				SaxQuantizer::SAX sax(input_argument.point_dimension);
				sax.get_SAX(normalized_series_vector, input_argument.point_dimension, apca_point_array[point_id]);
				break;
			}
			default:
				assert(0);
				break;
			}
			APCA_QUAL::get_minmax_apca(normalized_series_vector, apca_point_array[point_id], rectangle_insertion);
			result_record.representation_time = TOOL::recordFinishTime(TOOL::time_record[0]);
			input_argument.representation_time += result_record.representation_time;
			TOOL::recordStartTime(TOOL::time_record[1]);
			RTree.Insert(rectangle_insertion.m_min, rectangle_insertion.m_max, point_id); 
			input_argument.build_rtree_time += TOOL::recordFinishTime(TOOL::time_record[1]);
			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
			if (input_argument.representation_option == 9) {
				output_argument.sum_deviation = 0;
				output_argument.max_deviation = 0;
				output_argument.max_deviation_multiple_width = 0;
			}
			else {
				APCA_KNN_QUAL::distanceAE(normalized_series_vector, input_argument.time_series_length, apca_point_array[point_id], result_record);
				output_argument.max_deviation += result_record.max_deviation;
				output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
				output_argument.sum_deviation += result_record.sum_deviation;
			}
			TOOL::delete_rectangle(rectangle_insertion);
		}
		else if (input_argument.representation_option == 5) {
			typename CHEBYSHEV_QUAL::CHEBYSHEV chebyshev_MBR;
			CHEBYSHEV_QUAL::initialCHEBYSHEV(input_argument.time_series_length, input_argument.point_multi_single_dimension, chebyshev_MBR);
			CHEBYSHEV_QUAL::getCHEBYSHEV(input_argument, original_time_series, chebyshev_share, chebyshev_array[point_id]);
			CHEBYSHEV_QUAL::getChebyshevMBR(chebyshev_array[point_id], chebyshev_MBR);
			result_record.representation_time = TOOL::recordFinishTime(TOOL::time_record[0]);
			input_argument.representation_time += result_record.representation_time;
			TOOL::recordStartTime(TOOL::time_record[1]);
			RTree.Insert(chebyshev_MBR.f, chebyshev_MBR.coefficient, point_id);
			input_argument.build_rtree_time += TOOL::recordFinishTime(TOOL::time_record[1]);
			input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
			CHEBYSHEV_QUAL::get_sum_deviation(input_argument, normalized_series_vector, chebyshev_share, result_record);
			output_argument.max_deviation += result_record.max_deviation;
			output_argument.max_deviation_multiple_width += result_record.max_deviation_multiple_width;
			output_argument.sum_deviation += result_record.sum_deviation;
			CHEBYSHEV_QUAL::deleteCHEBYSHEV(chebyshev_MBR);
		}
		else {
		}
	}
	RTree.deleteAllLeafNodeMemory();	
	delete[] original_time_series;
	original_time_series = nullptr;
	normalized_series_vector.clear();
	normalized_series_vector.shrink_to_fit();
}


TEMPLATE
template<typename T>
T MULTI::compute_knn_accuracy(const int& const K, const multiset<pair<T, int>>& const squential_scan_result_set, const multiset<pair<T, int>>& const knn_result_set){
	assert(squential_scan_result_set.size() > K && knn_result_set.size() >= 0 && K > 0 && K != INF);
	T knn_accuracy = 0;

	set<int> sequential_time_series_id_set;
	int result_id = 0;
	for (multiset<pair<double, int>>::iterator squential_it = squential_scan_result_set.begin(); result_id < input_argument.K && squential_it != squential_scan_result_set.end(); squential_it++, result_id++) {
		sequential_time_series_id_set.emplace(squential_it->second);
	}
	assert(sequential_time_series_id_set.size() == K);
	if (knn_result_set.size() <= K) {

		for (multiset<pair<double, int>>::iterator knn_it = knn_result_set.begin(); knn_it != knn_result_set.end(); knn_it++) {
			if (sequential_time_series_id_set.find(knn_it->second) != sequential_time_series_id_set.end()) {
				knn_accuracy++;
			}
		}

		knn_accuracy /= double(K);
	}
	else {

		result_id = 0;

		for (multiset<pair<double, int>>::iterator knn_it = knn_result_set.begin(); result_id < K; knn_it++, result_id++) {
			if (sequential_time_series_id_set.find(knn_it->second) != sequential_time_series_id_set.end()) {
				knn_accuracy++;
			}
		}

		knn_accuracy /= knn_result_set.size();
	}
	
	if (knn_accuracy == 0) {
		knn_accuracy = 0.01 / double(K);
	}
	return knn_accuracy;
}

TEMPLATE
template<typename T, typename Y, typename U>
RTREE& MULTI::all_knn(typename TOOL::DATA_SOURCE& const data_source, const vector<T>& const query_time_series_vector, T*& query_time_series, const vector<Y>& const multi_y_projection_argument, const std::multiset<pair<U, int>>& const squential_scan_result_set) {
	input_argument.knn_total_time = 0;
	input_argument.result_accuracy = 0;
	std::multiset<pair<U, int>> result_set;
	switch (input_argument.representation_option) {
	case 1:
	case 6: 
	case 8: {
		result_set = all_knn_multi(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, apla_array_vector);
		break;
	}
	case 2: {
		result_set = all_knn_multi(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, pla_array);
		break;
	}
	case 3: 
	case 4:
	case 7:
	case 9: {
		result_set = all_knn_multi(data_source, input_argument, multi_y_projection_argument, query_time_series, RTree, apca_array);
		break;
	}
	case 5: {
		result_set = CHEBYSHEV_QUAL::KNNCHEBYMulti(input_argument, data_source, query_time_series, chebyshev_share, chebyshev_array, RTree);
		break;
	}
	default: assert(0);
		break;
	}
	if (input_argument.representation_option != 5) {
		assert(input_argument.pruning_power > 0 && input_argument.pruning_power <= 1);
	}
	input_argument.result_accuracy = compute_knn_accuracy(input_argument.K, squential_scan_result_set, result_set);
	assert(input_argument.result_accuracy > 0 && input_argument.result_accuracy <= 1);
	input_argument.prune_power_combine = input_argument.pruning_power / input_argument.result_accuracy;
	output_argument.run_time = input_argument.whole_run_time;
	result_set.clear();
	return RTree;
}


TEMPLATE
template<typename T, typename Y>
std::multiset<pair<double, int>> MULTI::all_knn_multi(typename TOOL::DATA_SOURCE& const data_source, typename TOOL::INPUT_ARGUMENT& const input_argument, const vector<Y>& const multi_y_projection_argument, DataType*& g_query_time_series, const RTREE& const apcaRTree, T& const approximation_array) {
	input_argument.knn_total_time = 0.0;
	input_argument.IO_cost = 0;
	input_argument.pruning_power = 0;
	int multi_to_single_series_length = input_argument.time_series_length;
	vector<DataType> query_time_series_vector;
	query_time_series_vector.resize(multi_to_single_series_length, INF);
	std::copy_n(g_query_time_series, multi_to_single_series_length, query_time_series_vector.begin());
	vector<DataType> reconstruct_query_time_series_vector;
	reconstruct_query_time_series_vector.resize(multi_to_single_series_length, INF);
	typename APCA_QUAL::APCA query_APCA;
	SaxQuantizer::SAX sax(input_argument.point_dimension);
	typename PLA_QUAL::PLA PLA_query;
	DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_linked_list = DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>();
	switch (input_argument.representation_option) {
	case 1:
	case 6:
	case 8: {
		query_linked_list.copy(apla_array_vector[input_argument.query_time_series_id]);
		APLA::getAPLAReconstructSeries(query_linked_list, reconstruct_query_time_series_vector);
		break;
	}
	case 2: {
		PLA_QUAL::initialPLA(PLA_query, input_argument.point_multi_single_dimension);
		PLA_QUAL::getPLA(input_argument.point_multi_single_length, input_argument.point_multi_single_dimension, g_query_time_series, PLA_query);
		break;
	}
	case 3:
	case 4:
		break;
	case 5:
		assert(0);
		break;
	case 7:
		break;
	case 9: {
		APCA_QUAL::initialAPCA(query_APCA, input_argument.point_dimension);
		sax.get_SAX(query_time_series_vector, input_argument.point_dimension, query_APCA);
		break;
	}
	default:
		assert(0);
		break;
	}
	int i = NULL, j = NULL;
	priority_queue<RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>> queue;
	RTREE_NODE_PAIR f_APLA_Root, f_temp_APLA_Pair;
	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> temp;
	list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR> result;
	std::multiset<pair<double, int>> result_set;
	typename APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair;
	RTREE_NODE_PAIR m_temp_queue_top;
	TOOL::recordStartTime(TOOL::time_record[14]); 
	TOOL::recordStartTime(TOOL::time_record[2]);
	f_APLA_Root.p_rtree_node = apcaRTree.m_root;
	f_APLA_Root.d_dist = 0;
	queue.push(f_APLA_Root);
	while (!queue.empty()) {
		m_temp_queue_top = queue.top();
		for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end();) {
			if (plist->d_dist <= m_temp_queue_top.d_dist) {
				result.push_back(*plist);
				plist = temp.erase(plist);
			}
			else
				plist++;
			if (input_argument.K == result.size()) {
				
				input_argument.knn_total_time += TOOL::recordFinishTime(TOOL::time_record[2]);
				input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);			
				input_argument.pruning_power = input_argument.IO_cost / double(input_argument.point_number);
				for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
					result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
				}
				assert(result_set.size() == input_argument.K);
				priority_queue< RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>>().swap(queue);
				temp.clear();
				result.clear();
				list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
				return result_set;
			}
		}
		queue.pop();
		if (m_temp_queue_top.p_rtree_node == nullptr) {
			DataType* original_time_series = new DataType[input_argument.point_multi_single_length];
			vector<DataType> original_time_series_vector;
			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;
			TOOL::read_normalized_multi_time_series(data_source, m_temp_queue_top.original_time_series_id, original_time_series_vector);
			input_argument.IO_cost++;
			tempOriginalTimeSeriesPair.d_dist = TOOL::distanceEUC(query_time_series_vector, original_time_series_vector);
			temp.push_back(tempOriginalTimeSeriesPair);
			original_time_series_vector.clear();
			original_time_series_vector.shrink_to_fit();
			delete[] original_time_series;
			original_time_series = nullptr;
		}
		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {
			typename APCA_QUAL::APCA QProjection;
			APCA_QUAL::initialAPCA(QProjection, input_argument.point_dimension);
			vector<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> query_apla_projection_vector(input_argument.point_dimension);
			f_temp_APLA_Pair.p_rtree_node = nullptr;
			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {
				f_temp_APLA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);
				switch (input_argument.representation_option) {
				case 1:
				case 6: 
				case 8: {
					f_temp_APLA_Pair.d_dist = APLA::get_distance_LB(apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], APLA::get_apla_projection(query_time_series_vector, apla_array_vector[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], query_apla_projection_vector));
					break;
				}
				case 2:
					f_temp_APLA_Pair.d_dist = PLA_QUAL::getPLADistance(input_argument.point_multi_single_length, pla_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], PLA_query, f_temp_APLA_Pair.d_dist);
					break;
				case 3:
				case 4:
				case 7: {
					f_temp_APLA_Pair.d_dist = APCA_KNN_QUAL::distanceLB(APCA_KNN_QUAL::QAPCAProjection(g_query_time_series, input_argument.time_series_length, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)], QProjection), apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
					break;
				}
				case 9: {
					f_temp_APLA_Pair.d_dist = sax.distance_LB_SAX(query_APCA, apca_point_array[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)]);
					break;
				}
				case 5:
					assert(0);
					break;
				default:
					assert(0);
					break;
				}
				queue.push(f_temp_APLA_Pair);
			}
			APCA_QUAL::deleteAPCA(QProjection);
		}
		else if (m_temp_queue_top.p_rtree_node->IsInternalNode()) {
			RTREE_NODE_PAIR temp_apla_pair;
			typename APCA_KNN_QUAL::REGION fs_region_G;
			typename APCA_KNN_QUAL::initialREGION(fs_region_G, input_argument.point_dimension);
			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {
				switch (input_argument.representation_option) {
				case 1:
				case 6:
				case 8:
					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
					break;
				case 2:
					temp_apla_pair.d_dist = PLA_QUAL::getPLAMBRDistance(input_argument, m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, PLA_query, temp_apla_pair.d_dist);
					break;
				case 3:
				case 4:
				case 7:
				case 9:
					temp_apla_pair.d_dist = APCA_KNN_QUAL::MINDISTQR(g_query_time_series, input_argument.point_multi_single_length, APCA_KNN_QUAL::getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));
					break;
				case 5:
					assert(0);
					break;
				default:
					assert(0);
					break;
				}
				temp_apla_pair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;
				queue.push(temp_apla_pair);
			}
			APCA_KNN_QUAL::deleteREGION(fs_region_G);
		}
		else {
			assert(0);
		}
	}
	input_argument.knn_total_time = TOOL::recordFinishTime(TOOL::time_record[2]);
	input_argument.whole_run_time += TOOL::recordFinishTime(TOOL::time_record[14]);
	input_argument.pruning_power = 1;
	priority_queue<RTREE_NODE_PAIR, vector<RTREE_NODE_PAIR>, MULTI::priorityIncrement<RTREE_NODE_PAIR>>().swap(queue);
	temp.clear();
	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(temp);
	for (typename list<APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it) {
		result_set.emplace(make_pair(it->d_dist, it->original_time_series_id));
	}
	assert(result_set.size() <= input_argument.K);
	result.clear();
	list <APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR>().swap(result);
	return result_set;
}
#endif
