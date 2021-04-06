#pragma once
#ifndef CAPLA_H
#define CAPLA_H

#define TIME_H

#include "pch.h"
#include "SHARE_TOOL.h"
#include "GEOMETRY_TOOL.h"
#include "lib/doublyLinkedList.h"
#include "CPLA.h"

TEMPLATE
class CAPLA : public RTREE, virtual public GEOMETRY, virtual public TOOL, virtual public PLA_QUAL {
public:
	struct TOOL::INPUT_ARGUMENT input_argument;
	struct TOOL::OUTPUT_ARGUMENT output_argument;
	struct BREAK_POINT_COEFFICIENT;
	struct BREAK_POINT_MAGNITUDE;
	struct AREA_COEFFICIENT;
	struct AREA_COEFFICIENT_SPEED;
	struct AREA_COEFFICIENT_SPEED_NO_MINMAX;
	struct AREA_COEFFICIENT_CONCISE;
	struct Y_AXIS_COEFFICIENT;
	struct SEGMENT_COEFFICIENT;
	struct APLA_COEFFICIENT;
	struct APLA_COEFFICIENT1;
	struct PLA_COEFFICIENT_CONCISE;
	struct SPLIT_COEFFICIENT;
	template<typename T>
	struct OPTIMIZATION_COEFFICIENT;
	template<typename T>
	struct UPPER_BOUND_COEFFICIENT;
	struct MagnitudeIncrease;
	struct AreaIncreasing;
	struct AreaDecreasing;
	struct ParallelogramHeightIncrease;//190104
	struct WidthIncreasing;//190104
	struct DensityIncrease;
	struct DensityIncrease_Pointer;
	struct DeviationIncreasing;//190115
	struct Width_Divide_Radius_Increase;//190315
	struct Set_Deviation_Decrease;//190624
	struct SPLIT_ID_INCREASE;//200211
	struct COMPARE_GREATER_SEGMENT_DENSITY;
	CAPLA() {};
	template<typename T>
	CAPLA(const T& const n, const T& const N);
	template<typename T>
	CAPLA(const T& const n, const T& const N, const bool& const change_file);
	template<typename T, typename Y>
	inline void initial_split_coefficients(T& const input_argument, Y& const output_argument);
	template<typename T, typename Y, typename U>
	inline void get_split_coefficients(const T& const input_argument, const U& const output_argument, const int& const split_methods_option, vector<Y>& const local_total_split_id_sum_deviation, vector<Y>& const local_total_split_id_shift, vector<Y>& const local_total_split_id_time, vector<Y>& const global_total_approximation_sum_deviation, vector<Y>& const global_total_approximation_time);
	template<typename T, typename Y, typename U>
	inline void get_split_coefficients(const T& const input_argument, const U& const output_argument, const int& const split_methods_option, vector<Y>& const local_total_split_id_sum_deviation, vector<Y>& const local_total_split_id_shift, vector<Y>& const local_total_split_id_time, vector<Y>& const global_total_approximation_sum_deviation, vector<Y>& const global_total_approximation_time, vector<Y>& const global_total_knn_prune_power);
	template<typename T, typename Y, typename U, typename U1>
	inline void get_initial_N_coefficients(const T& const input_argument, const U& const output_argument, const U1& const initial_N_number, vector<Y>& const total_initial_N_prune_power_vector, vector<Y>& const total_initial_N_sum_deviation_vector, vector<Y>& const total_initial_N_run_time_vector, vector<Y>& const total_initial_N_approximation_time_vector, vector<Y>& const total_initial_N_knn_time_vector);
	template<typename T, typename Y, typename U, typename U1>
	inline void get_initial_N_sort_N_coefficients(const T& const input_argument, const U& const output_argument, const U1& const N_size, const U1& const N_id, const U1& const initial_N_number, vector<Y>& const initial_N_by_N_prune_power_vector, vector<Y>& const initial_N_by_N_sum_deviation_vector, vector<Y>& const initial_N_by_N_run_time_vector, vector<Y>& const initial_N_by_N_approximation_time_vector, vector<Y>& const initial_N_by_N_knn_time_vector);
	void initialBREAK_POINT_COEFFICIENT(const typename TOOL::INPUT_ARGUMENT& input_argument, BREAK_POINT_COEFFICIENT& const break_point);
	void deleteBREAK_POINT_COEFFICIENT(BREAK_POINT_COEFFICIENT& const break_point);//181112
	void getAPLADifference(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, int*& const break_point_ID);//181112
	void copyRectangleCoefficient(vector<AREA_COEFFICIENT>& area_vector, const int vector_id, AREA_COEFFICIENT& temp_coefficient);//181205
	void initialRecArea(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//181211
	void initialRecArea0ForParallelogram(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//181211
	void initialRecArea0ForParallelogramNoPush(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190416

	void initialRecArea0ForParallelogramArray(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, SEGMENT_COEFFICIENT*& area_vector);//190417

	void initialRecArea0ForParallelogramArrayVector(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190429 Use vector to instead Array,

	void initialRecAreaPLAImprove(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190606 initial by PLA

	void initialAPLARightEndpoint(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190617 Use minmax point as right endpoint

	void initialAPLARightEndpoint3Sub(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190617 16:31 Use minmax point split segment into 3 subsegment

	void initialAPLARightEndpoint3SubAdaptive(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190705 11:31 Merge Use minmax point split segment into 3 subsegment, for flat stream, use few segments, for drastic change, use more segments

	void initialAPLARightEndpoint3SubAdaptiveMergeSplit(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190711 16:59 Merge Use minmax point split segment into 3 subsegment, for flat stream, use few segments, for drastic change, use more segments. Splie Merge when initial

	void initialAPLARightEndpoint3SubAdaptiveMergeSplitRecursive(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190719 Recursive merge and split segment
	//190812
	void initialAPLARightEndpoint3SubAdaptiveMergeSplitRecursiveSpeed(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190812 Recursive merge and split segment
	//190820 Use Linded list to instead Vector
	void initialAPLARightEndpoint3SubAdaptiveMergeSplitRecursiveSpeed(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//190820 Recursive merge and split segment
	//190918 Linked List
	void initialAPLARightEndpoint3SubAdaptiveMergeSplitRecursiveSpeed(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//190918 Recursive merge and split segment

	//190918 Linked List
	void initialAPLARightEndpoint3SubAdaptiveMergeSplitRecursiveSpeed1(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//190918 Recursive merge and split segment

	//191002 Linked List
	void initialAPLARightEndpoint3SubAdaptiveMergeSplitRecursiveSpeed2(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//191002 Recursive merge and split segment

	 //191209 Y projection initialization
	template<typename T>
	inline void initial_threshold_and_all(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, const int& const original_id, int& const min_max_count, const T& const segment_3, T& const temp_segment, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<T>& const all_segment_linked_list);

	//191209 Y projection initialization
	//200212 
	template<typename T, typename Y, typename U>
	inline void initial_threshold_and_all(const U& const input_argument, const vector<T>& const original_time_series_vector, const int& const original_id, int& const min_max_count, const Y& const segment_3, Y& const temp_segment, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<Y>& const all_segment_linked_list);

	//200314 Y projection initialization. delete minmax point. Concise above function. Delete redundant variable
	template<typename T, typename Y, typename U>
	inline void initial_threshold_and_all(const U& const input_argument, const T& const point_value, Y& const y_projection_argument);

	//200429 get the right endpoint of long segment_3 
	template<typename T, typename Y, typename U, typename T1>
	inline void get_segment_3_right_endpoint(const T& const segment_width_second, const Y& const remainder, const U& const pre_right_endpoint, T1& const temp_coefficient);

	//191209 get_minmax_y_projection
	template<typename T>
	void get_minmax_y_projection(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, const int& const pre_right_endpoint, T& const segment_3, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<T>& const all_segment_linked_list);

	//191209 get_minmax_y_projection
	//200212 Add template, use vector to instead pointer for time sereis 
	template<typename T, typename Y, typename U>
	void get_minmax_y_projection(U& const input_argument, const vector<T>& const original_time_series_vector, const int& const pre_right_endpoint, Y& const segment_3, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<Y>& const all_segment_linked_list);

	//200316 Add template, use vector to instead pointer for time sereis 
	template<typename T, typename Y, typename U>
	void get_minmax_y_projection(const U& const input_argument, const vector<T>& const original_time_series_vector, const int& const segment_begin_id, const int& const segment_end_id, int& const id_min_point, int& const id_max_point, Y& const y_projection_argument);

	//200424, get minmax point in long segment, use minmax point as endpoint to get three short segment. 
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void assign_three_segments_by_minmax(U& const input_argument, const vector<T>& const original_time_series_vector, const T1& const  pre_right_endpoint, int& const id_min_point, int& const id_max_point, Y& const segment_1, Y& const segment_2, Y& const segment_3, U1& const output_argument);

	//200428
	template<typename T, typename Y, typename U>
	inline void insert_segment_back_linkedlist(T& const temp_coefficient, DoublyLinkedList<Y>& const doubly_linked_list, U& const rest_segment_number);

	//191208 merge Y projection & MSPLA
	template<typename T>
	void initialMSPLA(int null, typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<T>& const all_segment_linked_list, DoublyLinkedList<T>& const cluster_linked_list, DoublyLinkedList<T>& const doubly_linked_list);

	//200212 merge Y projection & MSPLA
	//200212 Add template. Use vector to instead pointer for time series
	template<typename T, typename Y, typename U>
	void initialMSPLA(U& const input_argument, const vector<T>& const original_time_series_vector, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<Y>& const all_segment_linked_list, DoublyLinkedList<Y>& const cluster_linked_list, DoublyLinkedList<Y>& const doubly_linked_list);

	//200212 merge Y projection & MSPLA
	//200212 Add template. Use vector to instead pointer for time series
	template<typename T, typename Y, typename U, typename U1>
	void initialMSPLA(U& const input_argument, const vector<T>& const original_time_series_vector, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<Y>& const all_segment_linked_list, DoublyLinkedList<Y>& const cluster_linked_list, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//200731 loop to get split point at initial part.
	template<typename T, typename Y>
	void loop_to_get_split_point(const vector<T>& const original_time_series_vector, Y& const current_segment, Y& const accumulate_segment, Y& const next_two_points_segment, DoublyLinkedList<Y>& const doubly_linked_list);

	//210402
	template<typename T, typename Y, typename U, typename U1>
	void insert_last_segment_compute_merge_split_coefficients(const vector<T>& const original_time_series_vector, Y& const current_segment, multimap<U, DoublyListNode<Y>&, greater<U>>& const merge_segment_density_map, multimap<U, DoublyListNode<Y>&, std::greater<U>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//200803 get split point in one segment by biggest accumulation area
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void get_right_endpoint_by_accumulation_area(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//210402
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void get_right_endpoint_by_accumulation_area_no_split_merge(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//210301 No split & merge operation, no threshold
	//210301 No split & merge operation, no threshold
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void get_right_endpoint_by_accumulation_area_no_split_merge_threshold(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//200924 compute area difference of whole linked list
	template<typename T>
	long double get_whole_area_difference(const DoublyLinkedList<T>& const doubly_linked_list);

	//201005 compute area difference of whole linked list by MAP, not by linked list
	template<typename T>
	long double get_whole_area_difference_by_map(const T& const map);

	//201002 move the left / right endpoint of max area difference segment
	template<typename T, typename Y, typename T1>
	void optimize_segment_max_area_difference(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map);

	//201103 Speed up. move the left / right endpoint of max area difference segment
	template<typename T, typename Y, typename T1>
	bool optimize_segment_max_area_difference_speed(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map);

	template<typename T>
	inline void copy_optimization_coefficient(T& const segment_main, const T& const segment_copied);

	//************************************
	// Method:initial_segment_decrease_right
	// Qualifier: initial segemnt decrease right point, a&b, area difference.
	// date:201203  23:48
	// author:
	//************************************
	template<typename T, typename Y, typename U>
	inline bool initial_segment_decrease_right(const vector<T>& const original_time_series_vector, const Y& const segment_original, Y& const segment_decrease_right, U& const output_argument);

	//************************************
	// Method:initial_segment_increase_right
	// Qualifier: initial segemnt increase right point, a&b, area difference.
	// date:201205  16:58
	// author:
	//************************************
	template<typename T, typename Y, typename U>
	inline void initial_segment_increase_right(const vector<T>& const original_time_series_vector, const Y& const segment_original, Y& const segment_increase_right, U& const output_argument);

	//************************************
	// Method:initial_segment_increase_left
	// Qualifier: initial segemnt increase left point, a&b, area difference.
	// date:201204  00:39
	// author:
	//************************************
	template<typename T, typename Y, typename U>
	inline void initial_segment_increase_left(const vector<T>& const original_time_series_vector, const Y& const segment_original, Y& const segment_increase_left, U& const output_argument);

	//************************************
	// Method:initial_segment_decrease_left
	// Qualifier: initial segemnt decrease left point, a&b, area difference.
	// date:201211  20:02
	// author:
	//************************************
	template<typename T, typename Y, typename U>
	inline bool initial_segment_decrease_left(const vector<T>& const original_time_series_vector, const Y& const segment_original, Y& const segment_decrease_left, U& const output_argument);

	//************************************
	// Method:compare_segment_right_decrease
	// Qualifier: compare if segment decrease right endpoint can decrease area difference
	// date:201204  01:11
	// author:
	//************************************
	template<typename T, typename Y, typename T1>
	inline bool compare_segment_right_decrease(multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, T& const node_max_area_difference, Y& const segment_middle, const Y& const temp_segment_decrease_right, Y& const segment_right, const Y& const temp_segment_increase_left);

	//************************************
	// Method:compare_segment_left_decrease
	// Qualifier: Compare if segment decrease left endpoint can decrease area difference. Middle segment has max area difference
	// date:201213  09:19
	// author:
	//************************************
	template<typename T, typename Y, typename T1>
	inline bool compare_segment_left_decrease(multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, T& const node_max_area_difference, Y& const segment_middle, const Y& const temp_segment_decrease_left, Y& const segment_left, const Y& const temp_segment_increase_right);

	//************************************
	// Method:optimize_segment_max_area_difference_speed1
	// Qualifier: Optimization of segment adjust. Speed up. move the left / right endpoint of max area difference segment
	// date:201213  09:19
	// author:
	//************************************
	template<typename T, typename Y, typename T1>
	bool optimize_segment_max_area_difference_speed1(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map);

	//************************************
	// Method:erase_node_if_no_need_adjust
	// Qualifier:
	// date:201214  10:12
	// author:
	//************************************
	template<typename T, typename Y>
	inline bool erase_node_if_no_need_adjust(multimap<T, DoublyListNode<Y>& const, greater<T>>& const split_area_difference_map, const DoublyListNode<Y>& const node_candidate_erase);

	//************************************
	// Method:update_segment_map_by_min_area_difference
	// Qualifier:
	// date:201214  11:13
	// author:
	//************************************
	template<typename T, typename Y>
	inline void update_segment_map_by_min_area_difference(multimap<T, DoublyListNode<Y>& const, greater<T>>& const split_area_difference_map, Y& const segment_decrease, const Y& const temp_segment_decrease, Y& const segment_increase, const Y& const temp_segment_increase, DoublyListNode<Y>& const node_increase);

	//************************************
	// Method:update_segment_map_by_min_area_difference
	// Qualifier:
	// date:201223  15:41
	// author:
	//************************************
	template<typename T, typename Y>
	inline void update_segment_map_by_min_area_difference(multimap<T, DoublyListNode<Y>& const, greater<T>>& const split_area_difference_map, const Y& const temp_segment_move, DoublyListNode<Y>& const node_original);

	//************************************
	// Method:erase_if_adjacent_segments_are_zeros
	// Qualifier:
	// date:201215 16:01
	// author:
	//************************************
	template<typename T, typename Y>
	inline void erase_if_adjacent_segments_are_zeros(multimap<T, DoublyListNode<Y>& const, greater<T>>& const split_area_difference_map, DoublyListNode<Y>& const node_middle);

	//************************************
	// Method:adjust_compute_segment_right
	// Qualifier:
	// date:201215 20:31
	// author:
	//************************************
	template<typename T, typename Y, typename T1, typename U>
	bool adjust_compute_segment_right(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_middle, U& const output_argument);

	//************************************
	// Method:assert_if_finish_move
	// Qualifier:
	// date:201226 12:14
	// author:
	//************************************
	template<typename T>
	inline bool assert_if_finish_move(const T& const segment_middle, const T& const segment_right);

	//************************************
	// Method:compute_segment_decrease_right
	// Qualifier:
	// date:201226 09:34
	// author:
	//************************************
	template<typename T, typename Y, typename T1, typename T2, typename T3, typename T4, typename U>
	void compute_segment_decrease_right(const vector<T>& const original_time_series_vector, const T3 type_move, T4& const area_difference_min_other, multimap<T1, T2>& const map_new_area_difference_adjust_type, const Y& const segment_original_middle, const Y& const segment_original_right, Y& const segment_middle_decrease_right, Y& const segment_right_increase_left, U& const output_argument);

	//************************************
	// Method:compute_segment_increase_right
	// Qualifier:
	// date:201226 11:14
	// author:
	//************************************
	template<typename T, typename Y, typename T1, typename T2, typename T3, typename T4, typename U>
	void compute_segment_increase_right(const vector<T>& const original_time_series_vector, const T3 type_move, T4& const area_difference_min_other, multimap<T1, T2>& const map_new_area_difference_adjust_type, const Y& const segment_original_middle, const Y& const segment_original_right, Y& const segment_middle_increase_right, Y& const segment_right_decrease_left, U& const output_argument);

	//****************************************
	// Method:is_seg_out_map
	// Qualifier:short segment and flat segment outside split multimap
	// date:201228 06:51
	// author:
	//****************************************
	template<typename T>
	inline bool is_seg_out_map(const T& const segment);

	//****************************************
	// Method:label_left_node_erase_emplace_map
	// Qualifier: 1 Label left, middle node. 2 Erase, Emplace MAP
	// date:201228 15:23
	// author:
	//****************************************
	template<typename Y, typename T1>
	inline void label_left_node_erase_emplace_map(multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_middle);

	//****************************************
	// Method:label_right_node_erase_emplace_map
	// Qualifier: 1 Label middle, right node. 2 Erase, Emplace MAP
	// date:201228 15:36
	// author:
	//****************************************
	template<typename Y, typename T1>
	inline void label_right_node_erase_emplace_map(multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_middle);

	//****************************************
	// Method:update_left_node_label_map
	// Qualifier: 1 Update, Label left, middle node. 2 Update, Erase, Emplace MAP
	// date:201228 15:36
	// author:
	//****************************************
	template<typename Y, typename T1>
	inline void update_left_node_label_map(multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_middle, const Y& const segment_left_move, const Y& const segment_middle_move);

	//****************************************
	// Method:update_right_node_label_map
	// Qualifier: 1 Update, Label left, middle node. 2 Update, Erase, Emplace MAP
	// date:201228 15:36
	// author:
	//****************************************
	template<typename Y, typename T1>
	inline void update_right_node_label_map(multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_middle, const Y& const segment_middle_move, const Y& const segment_right_move);

	//****************************************
	// Method:adjust_compute_segment_2_sides
	// Qualifier:
	// date:201224 13:11
	// author:
	//****************************************
	template<typename T, typename Y, typename T1, typename U>
	bool adjust_compute_segment_2_sides(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_middle, U& const output_argument);

	//************************************
	// Method:optimize_segment_max_area_difference_speed2
	// Qualifier: Optimization of segment adjust, for max difference segment, adjust left / right endpoint of segment, choose the min area difference one.
	// date:201213  09:19
	// author:
	//************************************
	template<typename T, typename Y, typename T1>
	bool optimize_segment_max_area_difference_speed2(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map);

	//************************************
	// Method:optimize_segment_max_area_difference_speed3
	// Qualifier: Optimization of segment adjust, for max difference segment, adjust both left & right endpoint of segment, choose the min area difference from original one.
	// date:201224  09:19
	// author:
	//************************************
	template<typename T, typename Y, typename T1, typename U>
	bool optimize_segment_max_area_difference_speed3(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>& const, greater<T1>>& const split_area_difference_map, U& const output_argument);

	//201002 After initialization, split&merge operation, begin to optimizaiton of endpoints of segment.
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void optimization_segments_loop(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//201103 Speed up. After initialization, split&merge operation, begin to optimizaiton of endpoints of segment.
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void optimization_segments_loop_speed(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//201030 After initial part, The fastest way. Use split&merge&optimization to get tighter approximation
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void split_merge_optimization_segments_speed(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//************************************
	// Method: update_list_map_sub_left_is_merge
	// Qualifier: When candidate merge segment is splited sub left segment.
	// date:201229  14:44
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U>
	inline bool update_list_map_sub_left_is_merge(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_split_candidate, Y& const temp_segment_splited_sub_left, const Y& const temp_segment_splited_sub_right, const Y& const temp_segment_candidate_split_right, U& const output_argument);

	//************************************
	// Method: update_list_map_original_right_is_merge
	// Qualifier: When candidate merge segment is splited original right segment.
	// date:201229  21:19
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U>
	inline bool update_list_map_original_right_is_merge(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_split_candidate, const Y& const temp_segment_splited_sub_left, const Y& const temp_segment_splited_sub_right, Y& const temp_segment_candidate_split_right, U& const output_argument);

	//************************************
	// Method: update_list_map_split_is_not_merge
	// Qualifier: When candidate merge segment is splited original right segment.
	// date:201231  14:30
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U>
	inline bool update_list_map_split_is_not_merge(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, DoublyListNode<Y>& const node_split_candidate, DoublyListNode<Y>& const node_merge_candidate, DoublyListNode<Y>& const node_sub_left_new, const Y& const temp_segment_splited_sub_right, U& const output_argument);

	//************************************
	// Method: update_list_map_merge_split_same
	// Qualifier: When candidate split segment is candiate merge segment
	// date:210301  21:25
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U>
	inline bool update_list_map_merge_split_same(const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyListNode<Y>& const node_split_candidate, Y& const temp_segment_splited_sub_left, const Y& const temp_segment_splited_sub_right, U& const output_argument);

	//************************************
	// Method:split_merge_optimization_segments_speed0
	// Qualifier: After initial part, The fastest way. Use split&merge&optimization to get tighter approximation
	// date:201228  14:44
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void split_merge_optimization_segments_speed0(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//************************************
	// Method:split_merge_optimization_segments_speed1
	// Qualifier: After initial part, The fastest way. Use split&merge&optimization to get tighter approximation
	// date:201228  14:44
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void split_merge_optimization_segments_speed1(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//************************************
	// Method:split_merge_optimization_segments_speed2
	// Qualifier: After initial part, The fastest way. Use split&merge&optimization to get tighter approximation
	// date:210301  14:44
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void split_merge_optimization_segments_speed2(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//************************************
	// Method:split_merge_optimization_segments_speed3
	// Qualifier: After initial part, add speeded ms part, not only sm part.
	// date:210303  09:27
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U, typename U1>
	void split_merge_optimization_segments_speed3(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//200915 After initial part, Use split&merge&optimization to get tighter approximation
	template<typename T, typename T1, typename Y, typename U, typename U1, typename U2>
	void split_merge_optimization_segments(U& const input_argument, const vector<T>& const original_time_series_vector, const U2& const optimization_coefficients, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//************************************
	// Method:get_SAPLA_by_one_line_or_two_points
	// Qualifier: After initial part, add speeded ms part, not only sm part.
	// date:210307  16:03
	// author:
	//************************************
	template <typename T, typename T1, typename Y, typename U, typename U1>
	void get_SAPLA_by_one_line_or_two_points(U& const input_argument, const vector<T>& const original_time_series_vector, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	//200707 new initial method
	template<typename T, typename Y, typename U, typename U1>
	void initial_SAPLA_200706(U& const input_argument, const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);
	
	void initialAPLAYProjection(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<AREA_COEFFICIENT>& const all_linked_list, DoublyLinkedList<AREA_COEFFICIENT>& const cluster_linked_list, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//191028 Initial to check Y-projection & APLA
	//191115 for KNN algorithm
	template<typename T>
	void initial_apla_yprojection_knn(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<T>& const all_linked_list, DoublyLinkedList<T>& const cluster_linked_list, DoublyLinkedList<T>& const doubly_linked_list);

	void initialAPLARightEndpoint3SubNoMerge(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190627 16:38 Use minmax point split segment into 3 subsegment, when minmax points at endpoint, not merge.

	void initialRecArea0ForParallelogramArrayVector0(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector);//190605 Use vector to instead Array,

	template<typename T, typename T1, typename Y, typename U>
	void initial_y_projection_all_cluster(const U& const input_argument, const vector<T>& const original_time_series_vector, T1& const y_projection_argument, DoublyLinkedList<Y>& const all_segment_linked_list, DoublyLinkedList<Y>& const cluster_linked_list);

	inline double& getRecWidth(const vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, const int& const rec_num, AREA_COEFFICIENT& const temp_coefficient);//190115

	double& getRecWidth(vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, const int& const rec_num);//190319

	inline double& getRecWidth(SEGMENT_COEFFICIENT*& area_vector, const int& const vector_id, SEGMENT_COEFFICIENT& const temp_coefficient);//190418

	inline double& getRecWidthWidth(vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190605

	void getRightEndpointAndWidth(const vector<int>& const right_endpoint_vector, vector<AREA_COEFFICIENT>& const area_vector);//190624

	inline long double& getAandBSlopInterceptOnePoint(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190115

	template<typename T, typename Y>
	inline long double& getAandBSlopInterceptOnePoint(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);//200212

	template<typename T, typename Y>
	inline Y& getAandBSlopInterceptOnePoint(const T& const point_value, Y& const temp_coefficient);//190115

	inline AREA_COEFFICIENT& getAandBSlopInterceptTwoPoint(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190929

	template<typename T, typename Y>
	inline Y& getAandBSlopInterceptTwoPoint(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);//200212
	//210402
	template<typename T, typename Y>
	inline void getAandBSlopInterceptTwoPoint(const Y& const left_value, const Y& const right_value, T& const temp_coefficient);//191114

	AREA_COEFFICIENT& refreshSegmentCoefficient(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190318 for AREA_COEFFICIENT

	void refreshSegmentCoefficient(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector, const int& const vector_id);//190319 for vector

	AREA_COEFFICIENT& refreshSegmentCoefficient0ForParallelogram(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190404

	SEGMENT_COEFFICIENT& refreshSegmentCoefficient0ForParallelogram(DataType*& const original_time_series, SEGMENT_COEFFICIENT& const temp_coefficient);//190419

	void computeRectangleCoefficient(vector<AREA_COEFFICIENT>& area_vector, const int& const vector_id, const double& const lowest_id, const double& const lowest_value, const double& const heightest_id, const double& const heightest_value);//181205 old

	void computeRectangleCoefficient(vector<AREA_COEFFICIENT>& area_vector, const int& const vector_id, const typename GEOMETRY::POINT& const min_point, const typename GEOMETRY::POINT& const max_point);//190319 Improve input argument of function

	bool updateRectangle(vector<AREA_COEFFICIENT>& area_vector, const int& const vector_id, DataType*& const original_time_series);//181205 update adjacent rectangle after merge operation

	void mergeRectangleRange(vector<AREA_COEFFICIENT>& area_vector, int& first_id, int& last_id, DataType*& const original_time_series);//181210 Range

	//190410
	bool mergeRectangleRecursive(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector, int vector_id);

	bool isRecMonotonic(vector<AREA_COEFFICIENT>& area_vector, int& first_id, int& last_id);//181211

	bool moveRecEndpoint(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, DataType*& const original_time_series);//181211 Update Left current and reight rectangle

	//190319
	bool moveRectagnleEndpoint(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector, int vector_id);

	//190327
	bool moveRectagnleSplitpoint(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector, int vector_id);

	void mergeRectangle(vector<AREA_COEFFICIENT>& area_vector, const int vector_id);//181120

	void mergeRectangle(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, const AREA_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//181205 rcurrent rectangle merge.  No recursive.

	void mergeRectangle0(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, const AREA_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//181210 current rectangle merge right rectangle if current heightest value< right lowest value. merge. recursive

	void mergeRectangle1(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, const AREA_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//181210 current rectangle merge right rectangle and left rectangle. merge. No recursive No recursive

	void mergeRectangle2(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, const AREA_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//181210 current rectangle merge right rectangle and left rectangle. merge. recursive recursive

	void mergeRectangle0ForParallelogram(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, const AREA_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//190404 Only for ForParallelogram

	inline void mergeRectangle0ForParallelogramNoErase(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, int& const segment_number, const AREA_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//190410 Only for ForParallelogram, No erase function

	inline void mergeRectangle0ForParallelogramNoEraseArray(SEGMENT_COEFFICIENT*& area_vector, int& vector_id, const SEGMENT_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//190418 Only for ForParallelogram, No erase function,for array, speed up algorithm.

	void splitRectangle(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector, int vector_id);//190114

	int binarySearchIdEquality(DataType*& const original_time_series, int& const left_id, int& const right_id);//190709 get slope a & intercept b by 2 points

	double AdaptiveSearchIdEquality(DataType*& const original_time_series, int& const left_id, int& const right_id);//190716 get slope a & intercept b by adaptive length points

	double binarySearchIdPointDistance(DataType*& const original_time_series, int& const left_id, int& const right_id);//190715

	double AdaptiveSearchIdPointDistance(DataType*& const original_time_series, int& const left_id, int& const right_id);//190716 For every binary, compute left_a & right_a

	double AdaptiveSearchIdPointDistance2(DataType*& const original_time_series, int& const left_id, int& const right_id);//191008 For every binary, only comput changed left id or right id
	//191231 vector to 
	template<typename T>
	double AdaptiveSearchIdPointDistance2(const vector<T>& const original_time_series, int& const left_id, int& const right_id);//191231 For every binary, only comput changed left id or right id

	//200219
	template<typename T>
	double get_intersection_point_by_segment(const vector<T>& const original_time_series, const int& const left_id, const int& const right_id);

	//191207 
	double find_split_point_by_binary_gradient(DataType*& const original_time_series, const int& const left_id, const int& const right_id); //For every binary, only comput changed left id or right id by middle gradient
	//191206 By Gradient Descent
	int find_split_point_by_gradient_descent(DataType*& const original_time_series, const int& const left_id, const int& const right_id);

	double findSplitSegmentByEndpoint2length(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190709 Find split point by a&b of length 2 endpoint of segment, if no equal, pick mid point.
	//191231 vector to instead pointer
	template<typename T, typename Y>
	double findSplitSegmentByEndpoint2length(const vector<T>& const original_time_series, const Y& const temp_coefficient);//191231 Find split point by a&b of length 2 endpoint of segment, if no equal, pick mid point.

	//200219
	template<typename T, typename Y>
	int find_split_point_by_direct_intersection_point(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient);

	// 200221
	// regard the middle point of segment split id
	template<typename T>
	inline double find_split_point_by_middle_point(const T& const temp_coefficient);

	//200205 Use optimal method to find split point in one segment.
	//split_point_type: 0 best split point. 1 min triangle density
	template<typename T, typename Y>
	double findSplitSegmentBaseline(const vector<T>& const original_time_series, const Y& const temp_coefficient, const int split_point_type);

	//200210 get split coefficient: left & right a&b, triangle densitys
	// return trianle density
	template<typename T, typename Y>
	long double& get_split_coefficient(const vector<T>& const original_time_series_vector, Y& const temp_coefficient, SPLIT_COEFFICIENT& const split_coefficent);

	//210402 get split coefficient: left & right a&b, triangle densitys
	template<typename T, typename Y, typename U>
	long double& get_split_coefficient(const vector<T>& const original_time_series_vector, Y& const temp_coefficient, SPLIT_COEFFICIENT& const split_coefficent, const int& const split_id_left, const int& const split_id_right, vector<U>& const segments_density_vector, const int& const split_point_left_segment_id);

	//200301, Use linekd list to instead vector get split coefficient: left & right a&b, triangle densitys
	// return trianle density
	template<typename T, typename Y, typename U>
	long double& get_split_coefficient(const vector<T>& const original_time_series_vector, Y& const temp_coefficient, SPLIT_COEFFICIENT& const split_coefficent, const int& const split_id_left, const int& const split_id_right, DoublyLinkedList<U>& const segments_density_linked_list, const int& const split_point_left_segment_id);

	//210402 Find min density point fast use vector. prove coefficeints for sub left segment, sub irght segment, and long smegnetnr 
	template<typename T, typename Y, typename U>
	U get_segment_min_density_local_id(const vector<T>& const original_time_series_vector, Y& const original_long_segment, vector<U>& const segments_density_vector);

	// Use Linked List instead vector. Find min density point fast use vector. prove coefficeints for sub left segment, sub irght segment, and long smegnetnr 
	template<typename T, typename Y, typename U>
	U get_segment_min_density_local_id(const vector<T>& const original_time_series_vector, Y& const original_long_segment, DoublyLinkedList<U>& const segments_density_linked_list);

	//210402 Already know local mind density split id. Try to search border section. make local section as global section. Find min density point fast use vector. prove coefficeints for sub left segment, sub irght segment, and long smegnetnr 
	template<typename T, typename Y, typename U>
	U& get_segment_min_density_global_id(const vector<T>& const original_time_series_vector, Y& const long_segment, U& const min_density_segment_local, vector<U>& const segments_density_vector);

	//200302 Use linked list ot instead vector. Already know local mind density split id. Try to search border section. make local section as global section. Find min density point fast use vector. prove coefficeints for sub left segment, sub irght segment, and long smegnetnr 
	template<typename T, typename Y, typename U>
	U& get_segment_min_density_global_id(const vector<T>& const original_time_series_vector, Y& const long_segment, U& const min_density_segment_local, DoublyLinkedList<U>& const segments_density_linked_list);

	//200210 Find min density point fast use vector
	template<typename T, typename Y>
	double find_split_point_by_min_density_fast(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	//210402
	template<typename T, typename Y, typename U>
	double find_split_point_by_min_density_fast(const vector<T>& const original_time_series_vector, Y& const sub_segment_left, Y& const sub_segment_right, U& const temp_coefficient);

	//200225 contain 5 split point methods 1 min denisty , 2 BInary 3 Intersectoin point 4 middle point 5 best split id 
	template<typename T, typename Y>
	double group_find_split_point_methods(const vector<T>& const original_time_series_vector, Y& const temp_coefficient, const int& const option_split_method);

	//210402 contain 5 split point methods 1 min denisty , 2 BInary 3 Intersectoin point 4 middle point 5 best split id 
	template<typename T, typename Y, typename U>
	double group_find_split_point_methods(const vector<T>& const original_time_series_vector, const int& const option_split_method, Y& const sub_segment_left, Y& const sub_segment_right, Y& const temp_coefficient, U& const output_argument);

	//190715 get middle point of one segment
	template<typename T>
	inline int findSegmentMiddleID(T& const temp_coefficient);//190715

	inline double& getAreaDifference(AREA_COEFFICIENT& const temp_coefficient); //Min Max area - PLA area

	//190918
	inline double& getAreaDifference(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient); // Right end point area - PLA area

	//200212 change original time series from pointer to vector. Add template
	template<typename T, typename Y>
	inline long double& getAreaDifference(const vector<T>& const original_time_series_vector, Y& const temp_coefficient); // MinMax Area - PLA Area

	//201004 10:13 without original time series
	//template<typename T>
	//inline long double& get_area_difference_segment(T& const temp_coefficient); // MinMax Area - PLA Area

	//201007 02:21 Use the endpoint height difference instead area difference
	template<typename T, typename Y>
	inline long double& get_area_difference_segment(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);
	template<typename T, typename Y>
	inline long double& get_segment_max_deviation_instead(const T& const left_endpoint_value, const T& const right_endpoint_value, Y& const temp_coefficient);

	//210112 09:10 Compute upper bound from left&right min&max points
	template<typename T, typename Y>
	long double& get_segment_max_deviation_instead1(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	//190722
	void SplitSegmentByRealPLAArea(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector, const int& const split_number); //190722 find which point should be splited in one segment
	//190823
	void SplitSegmentByRealPLAArea(DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list, const int& const split_number); //190823 find which point should be splited in one segment

	//191006
	void splitSegmentBySplitedPoint(DataType*& const original_time_series, AREA_COEFFICIENT& const sub_left_segment, AREA_COEFFICIENT& const long_segment);

	//200212 Add template. Use pointer to instead vector for time series
	template<typename T, typename Y>
	void splitSegmentBySplitedPoint(const vector<T>& const original_time_series_vector, Y& const sub_left_segment, Y& const long_segment);

	//200226 Add split method option.  Add template. Use pointer to instead vector for time series
	template<typename T, typename Y>
	void splitSegmentBySplitedPoint(const vector<T>& const original_time_series_vector, const int& const option_split_method, Y& const sub_left_segment, Y& const long_segment);

	//200226 Add split method option.  Add template. Use pointer to instead vector for time series
	template<typename T, typename Y, typename U>
	void splitSegmentBySplitedPoint(const vector<T>& const original_time_series_vector, const int& const option_split_method, Y& const sub_left_segment, Y& const long_segment, U& const output_argument);

	//200821 find if multimap has key
	template<typename T, typename Y, typename U>
	bool find_if_in_multimap(const T& const key, const multimap<U, Y, greater<U>>& const multi_map);

	template<typename T, typename T1, typename Y>
	inline bool find_key_in_map(const multimap<Y, DoublyListNode<T>&, std::greater<Y>>& const sorted_multimap, const T1& const key);

	//************************************
	// Method:find_endpoint_key_in_map
	// Qualifier: if key, right endpoint in map
	// Input: map, key, right endpoint
	// Output: right endpoint, key in map
	// date:201228
	// author:
	//************************************
	template<typename T, typename T1, typename Y, typename U>
	bool find_endpoint_key_in_map(const multimap<Y, DoublyListNode<T>&, std::greater<Y>>& const sorted_multimap, const U& const unique_right_endpoint, const T1& const key);

	//190930
	template<typename T, typename T1, typename Y, typename U>
	void eraseMapByKey(multimap<Y, DoublyListNode<T>&, std::greater<Y>>& const sorted_multimap, const U& const unique_right_endpoint, const T1& const key);

	//190930
	template<typename T, typename T1, typename T2, typename Y, typename U>
	void updateMapByKey(multimap<Y, DoublyListNode<T>&, std::greater<Y>>& const sorted_multimap, const U& const unique_right_endpoint, const T1& const key, const T2& const new_value);

	//190830
	void splitSegmentBySlope(DataType*& const original_time_series, multimap<double, DoublyListNode<AREA_COEFFICIENT>&, std::greater<double>>& const merge_segment_density_map, multimap<double, DoublyListNode<AREA_COEFFICIENT>&, std::greater<double>>& const split_area_difference_map, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);

	//200212 change original time series from pointer to vector. Add templatea
	template<typename T, typename Y>
	void splitSegmentBySlope(const vector<T>& const original_time_series_vector, multimap<double, DoublyListNode<Y>&, std::greater<double>>& const merge_segment_density_map, multimap<double, DoublyListNode<Y>&, std::greater<double>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list);

	//210402 Add split metyhod option change original time series from pointer to vector. Add templatea
	template<typename T, typename Y, typename U, typename T1>
	void splitSegmentBySlope(const vector<T>& const original_time_series_vector, const int& const split_method_option, multimap<T1, DoublyListNode<Y>&, std::greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, std::greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U& const output_argument);

	//191124 for burst time series that need to split segment
	template<typename T>
	void split_burst_segment(DataType*& const original_time_series, const int& const  point_dimension, DoublyLinkedList<T>& const doubly_linked_list);

	//191124 for burst time series that need to split segment
	//200212 change original time series from pointer to vector. Add template
	template<typename T, typename Y>
	void split_burst_segment(const vector<T>& const original_time_series_vector, const int& const  point_dimension, DoublyLinkedList<Y>& const doubly_linked_list);
	inline void mergeSegmentByDensity(DataType*& const original_time_series, multimap<double, DoublyListNode<AREA_COEFFICIENT>&, std::greater<double>>& const merge_segment_density_map, multimap<double, DoublyListNode<AREA_COEFFICIENT>&, std::greater<double>>& const split_area_difference_map, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);
	//200212 change original time series from pointer to vector. Add templatea
	template<typename T, typename Y, typename U, typename T1>
	inline void mergeSegmentByDensity(const vector<T>& const original_time_series_vector, multimap<U, DoublyListNode<Y>&, std::greater<U>>& const merge_segment_density_map, multimap<U, DoublyListNode<Y>&, std::greater<U>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, T1& const output_argument);
	//190712
	inline void mergeEndSegment(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& area_vector, AREA_COEFFICIENT& const temp_coefficient);//190712 18:15 new segment is same with end vector, merge them.
	template<typename T, typename Y, typename U, typename T1>
	inline void split_merge_optimization(const vector<T>& const original_time_series_vector, const int& const split_method_option, multimap<T1, DoublyListNode<Y>&, std::greater<T1>>& const merge_segment_density_map, multimap<T1, DoublyListNode<Y>&, std::greater<T1>>& const split_area_difference_map, DoublyLinkedList<Y>& const doubly_linked_list, U& const output_argument);

	bool isSymmetry(const AREA_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//181220 Whether max point or min point in the middle

	bool isSymmetryForParallelogram(const AREA_COEFFICIENT& const temp_coefficient, DataType*& const original_time_series);//190408 Whether max point or min point in the middle

	void getMergedInfo(DataType*& const original_time_series, const vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//181129

	void getMergedInfor0ForParallelogram(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, int& const vector_id, int& const segment_number, AREA_COEFFICIENT& const temp_coefficient);//181129

	void getMergedInfor0ForParallelogramImprove(vector<AREA_COEFFICIENT>& const area_vector, int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//181129
	//190825 Use linked list to instead vector
	void getMergedInfor0ForParallelogramImprove(DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list, int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190825
	//190904
	double getMergedSubSegmentInformation(DoublyListNode<AREA_COEFFICIENT>& const node);
	//190918
	double getMergedSubSegmentInformation(DataType*& const original_time_series, DoublyListNode<AREA_COEFFICIENT>& const node);

	//210402
	template<typename T, typename Y, typename U>
	double getMergedSubSegmentInformation(const vector<T>& const original_time_series_vector, DoublyListNode<Y>& const node, U& const output_argument);

	//210402
	template<typename T, typename Y, typename U>
	double getMergedSubSegmentInformation(const vector<T>& const original_time_series_vector, Y& const left_segment, Y& const right_segment, Y& const merged_segment, U& const output_argument);

	void getMergedInforQueue(vector<AREA_COEFFICIENT>& const area_vector, int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190607

	void getMergedInfor0ForParallelogramArray(DataType*& const original_time_series, SEGMENT_COEFFICIENT*& area_vector, int& const vector_id, SEGMENT_COEFFICIENT& const temp_coefficient);//190418

	double& computeParallelogram(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190102

	double& computeParallelogramNoProjection(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190404 the funciton is the same with above, only improve the process, not project the deviaiton point to the left segment.

	double& computeParallelogramNoProjection(DataType*& const original_time_series, SEGMENT_COEFFICIENT& const temp_coefficient);//190419 the funciton is the same with above, only improve the process, not project the deviaiton point to the left segment.

	std::vector<typename GEOMETRY::POINT> getPointFromSeg(DataType*& const original_time_series, const AREA_COEFFICIENT& const temp_coefficient);//190306

	double& getPolygonArea(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190123


	double& getConvexHullArea(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190227

	
	void getSegmentMinMaxRectangleAreaDensity(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector);//190701
	//190826
	void getSegmentMinMaxRectangleAreaDensity(DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//190826

	//get a & b of segment PLA by original time series, only know width & right end point
	double& getAAndBByPLA(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190107 return a;

	// 191114 add template vector. get a & b of segment PLA by original time series, only know width & right end point
	template<typename T>
	T& getAAndBByPLA(const vector<DataType>& const original_time_series, T& const temp_coefficient);//190107 return a;

	//210402 get apla coeffcient minus advisor
	template<typename T>
	inline void get_apla_coefficients_segment(T& const temp_coefficient);

	//200504 Aleady know is flat segment. get 1 a & b£¬ 2 min&max point of segment PLA by original time series, only know width & right end point
	template<typename T, typename Y>
	inline void get_ab_minmax_flat_segment(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	//200730 compute a&b min&max of two points segment
	template<typename T, typename Y>
	inline void get_segment_two_points_ab_minmax(const T& const left_value, const T& const right_value, Y& const temp_coefficient);

	//210402 compute a&b min&max of accumulation segment from current segment.
	template<typename T, typename Y>
	inline Y& get_ab_minmax_segment_by_accumulation(const T& const accumulate_point_value, const Y& const last_segment, Y& const temp_coefficient);

	//210402 compute a&b, min&max of accumulation segment from current segment.
	template<typename T, typename Y>
	inline Y& get_ab_minmax_segment_by_accumulation(const T& const accumulate_point_value, Y& const current_segment);

	//200420 get 1 a & b£¬ 2 min&max point of segment PLA by original time series, only know width & right end point
	template<typename T, typename Y>
	Y& get_ab_minmax_segment(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	// 200221 add template vector. get a & b of segment PLA by original time series, only know width & right end point
	template<typename T, typename Y>
	Y& get_ab_segment(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	// 210402 get a&b of long segment by one point and old short segment
	template<typename T, typename Y>
	inline Y& get_ab_segment_by_accumulation(const vector<T>& const original_time_series_vector, const Y& const last_segment, Y& const temp_coefficient);

	// get a& b of long segment by right one pointand old short segment
	template<typename T, typename Y>
	inline Y& get_ab_segment_by_accumulation(const T& const accumulate_point_value, const Y& const last_segment, Y& const temp_coefficient);

	// 210402 get a&b of accumulate segment by one point and current segment
	template<typename T, typename Y>
	inline void get_ab_accumulation_segment(const T& const accumulate_point_value, Y& const current_segment);

	//201021 get a& b of long segment by left one pointand old short segment
	template<typename T, typename Y>
	inline Y& get_ab_segment_by_accumulation_left(const T& const point_value_left, const Y& const sub_right_segment, Y& const temp_coefficient);

	// 210402 get a1&b1 of short left segment by one point and old long segment
	template<typename T, typename Y>
	inline Y& get_ab_segment_by_decrement(const vector<T>& const original_time_series_vector, const Y& const last_segment, Y& const temp_coefficient);

	// 201021 get a1&b1 of short left segment by right one point and old long segment
	template<typename T, typename Y>
	inline Y& get_ab_segment_by_decrement_right(const T& const point_value_right, const Y& const long_segment, Y& const sub_left_segment);

	// 201022 get a2&b2 of short right segment by left one point and old long segment
	template<typename T, typename Y>
	inline Y& get_ab_segment_by_decrement_left(const vector<T>& const original_time_series_vector, const Y& const long_segment, Y& const sub_right_segment);

	// 201021 get a2&b2 of short right segment by left one point and old long segment
	template<typename T, typename Y>
	inline Y& get_ab_segment_by_decrement_left(const T& const point_value_left, const Y& const long_segment, Y& const sub_right_segment);

	double& getABSlopeInterceptOriginal(DataType*& const original_time_series, const int& const left_id, const int& const right_id, double& const slope_a, double& const intercept_b);//190716 get slope a & intercept b by original time series

	template<typename T>//191231 vector to instead pointer
	double& getABSlopeInterceptOriginal(const vector<T>& const original_time_series, const int& const left_id, const int& const right_id, double& const slope_a, double& const intercept_b);//191231 get slope a & intercept b by original time series

	double& getAAndBByPLAShortSeg(DataType*& const original_time_series, const vector<AREA_COEFFICIENT>& const area_vector, int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190412 return a&b from short sub segment;

	void getAAndBByPLAShortSegSpeed(const vector<AREA_COEFFICIENT>& const area_vector, int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190611 No evaluation method and sum value and a b
	
	inline void getAAndBByPLAShortSegSpeed(DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list, int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190825 No evaluation method and sum value and a b
	//210402
	template<typename T, typename Y>
	inline void getAAndBByPLAShortSegSpeed(const T& const left_segment, const T& const right_segment, Y& const merged_segment);

	template<typename T>
	inline T& getSubRightAAndBByPLA(const T& const left_segment, const T& const merged_segment, T& const right_segment);

	template<typename T>
	inline T& getSubLeftAAndBByPLA(const T& const right_segment, const T& const merged_segment, T& const left_segment);

	//190926
	// Input: 1 Original time series: long segment a&b, sum, width, right endpoint; 2 left(right) sub segment: a&b, sum, width, right endopint;
	// Output: right(left) sub segment: a&b, sum, sub left right min&max point
	template<typename T>
	void getSubAAndBByPLA(DataType*& const original_time_series, T& const left_segment, T& const right_segment, const T& const merged_segment);

	//200210
	// Input: 1 Original time series vector: long segment a&b, sum, width, right endpoint; 2 left(right) sub segment: a&b, sum, width, right endopint;
	// Output: right(left) sub segment: a&b, sum, sub left right min&max point
	template<typename T, typename Y>
	void getSubAAndBByPLA(const vector<T>& const original_time_series, Y& const left_segment, Y& const right_segment, const Y& const merged_segment);

	double& getAAndBByPLAShortSegArray(DataType*& const original_time_series, SEGMENT_COEFFICIENT*& area_vector, int& const vector_id, SEGMENT_COEFFICIENT& const temp_coefficient);//190418 return a&b from short segment;

	template<typename T, typename Y>
	inline long double get_sum_value(const T& const a, const T& const b, const Y& const segment_length);
	double& getSegmentDevByPLA(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190115 return deviation


	// 210120 Only return sum value. get segment sum deviaiton,  Know right endpoint, width, a & b of every segment
	template<typename T, typename Y>
	long double get_segment_sum_deviation(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient);

	// 210120 Only return sum value. get segment sum deviaiton,  Know right endpoint, width, a & b of every segment
	//template<typename T, typename Y>
	//long double get_segment_sum_deviation(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient);

	//210120 evaluate triangle area lower bound sum deviation of segment
	template<typename T, typename Y>
	bool evaluate_segment_sum_deviation_vs_triangle_area(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient);

	// 200220 get segment sum deviation by split point
	// 200113  get segment sum deviaiton, Know right endpoint, width, a & b of every segment
	template<typename T, typename Y, typename U>
	double get_segment_sum_deviation_by_split_id(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient, const U& const split_id);//200113

	// 201023 Return Sum & Max deviation. Know right endpoint, width, a & b of one segment
	template<typename T, typename Y, typename U>
	long double get_segment_sum_max_deviation(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient, U& const sum_deviation, U& const max_deviation);

	//210120
	template<typename T, typename Y>
	bool assert_merge_deviation(const vector<T>& const original_time_series_vector, const Y& const left_segment, const Y& const right_segment, const Y& const merged_segment);
	//210327
	template<typename T, typename Y>
	long double& get_sum_segment_absolute_difference(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	//void compareAreaPercentage(vector<AREA_COEFFICIENT>& area_vector, const int vector_id);//181122 move right_endpoint

	void compareAreaPercentage0(vector<AREA_COEFFICIENT>& area_vector, const int vector_id);//181128 Only compaer percentage

	void compareAreaPercentage(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, DataType*& const original_time_series);//181205

	//int chooseMergeSegmentID();//190315

	void compareAreaPercentageMinMaxDist(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, DataType*& const original_time_series);//181211

	void compareAreaPercentageMinMaxDistNoMove(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, DataType*& const original_time_series);//181211

	//void compareParallelogram(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, DataType *& const original_time_series);//190102

	void compareMiniMaxDist(vector<AREA_COEFFICIENT>& area_vector, const int vector_id);//181126

	void compareDeviation(vector<AREA_COEFFICIENT>& area_vector, int& vector_id, DataType*& const original_time_series);//190115

	void getAPLAArea(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//181119

	void getAPLAArea0(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//181122

	void getAPLAArea1(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//181122

	void getAPLAAreaMinMaxDist(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//181123

	void getAPLAArea3(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//181211 try initial area and compareAreaPercentageMinMaxDist()

	void getAPLAArea4LogLoop(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190326 Improve speed. Every loop, merge adjacent segment.

	void getAPLAArea5LowerBound(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190327 Improve speed. Every loop, merge adjacent segment. And because pre-merged two small segment PLA is better than merged one big PLA, this is lover bound

	void getAPLAArea5LowerBoundSpeed(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190404 Speed up above algorithm

	void getAPLAArea5LowerBoundSpeedNoErase(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190409 Speed up above algorithm, no erase() function(it cost so much time)

	void mergeOperation0(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, AREA_COEFFICIENT& const temp_coefficient);//190609

	void mergeOperationBatch2(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, vector<AREA_COEFFICIENT>& const temp_coefficient_vector, AREA_COEFFICIENT& const temp_coefficient);//190610 Improve mergeOperation0, use batch to speed up

	void mergeOperationBatch2Set(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, vector<AREA_COEFFICIENT>& const temp_coefficient_vector, AREA_COEFFICIENT& const temp_coefficient);//190610 Use set to instead priority_queue;

	void mergeOperationBatch2SetPointer(typename TOOL::INPUT_ARGUMENT& input_argument, vector<AREA_COEFFICIENT>& const area_vector, vector<AREA_COEFFICIENT>& const temp_coefficient_vector, AREA_COEFFICIENT& const temp_coefficient, vector<int>& const count_break_point_vector);//190612 Add sub segemt to speed up,avoid repeat calculate, Use set to instead priority_queue;

	void mergeOperationBatch2SetSpeed(typename TOOL::INPUT_ARGUMENT& input_argument, vector<AREA_COEFFICIENT>& const area_vector, vector<AREA_COEFFICIENT>& const temp_coefficient_vector, AREA_COEFFICIENT& const temp_coefficient);//190611 Use set to instead priority_queue,speed;

	void mergeOperationBatch2SetSpeedPointer(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, vector<AREA_COEFFICIENT>& const temp_coefficient_vector, AREA_COEFFICIENT& const temp_coefficient);//190612 Add sub segemt to speed up,avoid repeat calculate, Use set to instead priority_queue;
	//190823
	void mergeOperationBatch2SetSpeedPointer(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const temp_coefficient_vector, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//190823 Add sub segemt to speed up,avoid repeat calculate, Use set to instead priority_queue;
	//190827
	void mergeSplitOperationLink(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);

	//200212 Add template. Change time series from pointer to vector
	template<typename T, typename Y, typename U>
	void mergeSplitOperationLink(U& input_argument, const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const doubly_linked_list);

	//200212 Add template. Change time series from pointer to vector
	template<typename T, typename Y, typename U, typename U1>
	void mergeSplitOperationLink(U& input_argument, const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const doubly_linked_list, U1& const output_argument);

	 //190729 get adaptive right end point throgh y projection.
	void yProjectionAPLAInitial(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const map<double, int>& const whole_difference_map, vector<AREA_COEFFICIENT>& const area_vector);
	//191028 Check if burst time series
	bool checkIfBurstDataSet(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument);
	//191206 for threshold , max value may also be threshold
	bool check_if_burst_dataset(typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument);
	//191028 split original time series by all y projection value
	void initialYProjectionAll(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series,  const typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<AREA_COEFFICIENT>& const all_linked_list);
	//191030 split original time series by flat segment & mountain segment(unflat segment)
	void initialYProjectionCluster(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, const DoublyLinkedList<AREA_COEFFICIENT>& const all_linked_list, DoublyLinkedList<AREA_COEFFICIENT>& const cluster_linked_list);
	//191210
	template<typename T>
	void get_y_all_cluster_segment(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<T>& const all_segment_linked_list, DoublyLinkedList<T>& const cluster_linked_list);

	//191210
	template<typename T, typename Y, typename U>
	void get_y_all_cluster_segment(U& input_argument, const vector<T>& const original_time_series_vector, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<Y>& const all_segment_linked_list, DoublyLinkedList<Y>& const cluster_linked_list);

	//200213 whether is flat point
	template<typename T, typename U>
	inline bool is_flat_point(const U& const y_projection_argument, const T& const data);

	//200213 get rectangle width of current segment
	template<typename T, typename U>
	inline double& y_projection_get_segment_width(const DoublyLinkedList<U>& const all_segment_linked_list, T& const temp_coefficient);

	//200213 new Y projection method. for complicated situation, frist point is or not flat point.
	template<typename T, typename Y, typename U>
	void get_y_all_cluster_segment_200213(U& input_argument, const vector<T>& const original_time_series_vector, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<Y>& const all_segment_linked_list, DoublyLinkedList<Y>& const cluster_linked_list);

	//190730 y projection, merge initialed y projection method.
	void yProjectionAPLAMerge(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const map<double, int>& const whole_difference_map, vector<AREA_COEFFICIENT>& const area_vector);
	//191029 Linked list. y projection, merge & split. initialed y projection method.
	void YProjectionAPLALink(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);
	//191030 merge & split segment through all_linked_list & cluster_linked_list
	void YProjectionMerge(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const DoublyLinkedList<AREA_COEFFICIENT>& const all_linked_list, const DoublyLinkedList<AREA_COEFFICIENT>& const cluster_linked_list, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);

	//191209 merge & split segment by y_projection
	template<typename T>
	void y_projection_merge_split(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, const DoublyLinkedList<T>& const all_linked_list, const DoublyLinkedList<T>& const cluster_linked_list, DoublyLinkedList<T>& const doubly_linked_list);

	//200212 merge & split segment by y_projection
	template<typename T, typename Y, typename U>
	void y_projection_merge_split(U& input_argument, const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const all_linked_list, const DoublyLinkedList<Y>& const cluster_linked_list, DoublyLinkedList<Y>& const doubly_linked_list);

	//191116 when the time series only has one value
	template<typename T>
	void y_projection_merge_line(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, DoublyLinkedList<T>& const doubly_linked_list);

	//191116 when the time series only has one value
	template<typename T, typename Y, typename U>
	void y_projection_merge_line(U& input_argument, const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const doubly_linked_list);

	void mergeOperationJump1(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, AREA_COEFFICIENT& const temp_coefficient);//190610 loop +=2

	void getAPLAAreaPLAImprove(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190606, APLA based on PLA

	void getAPLAAreaPLAImproveNoPointer(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190612, APLA based on PLA, minimize struct, no pointer

	void getAPLAByMinMax(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190617  For every segment, use minmax point as endpoint. Then for every segment use same subsegemnt number to merge.
	//190911
	template<typename T>
	void getAPLAByMinMaxLinkedList(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<T>& const all_linked_list, DoublyLinkedList<T>& const cluster_linked_list, typename TOOL::OUTPUT_ARGUMENT& output_argument);// Use Linked List to instead vector, For every segment, use minmax point as endpoint. Then for every segment use same subsegemnt number to merge.

	//191108 no other evaluation method. directly get APLA approximation
	template<typename T>
	void get_APLA_point(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<T>& const all_linked_list, DoublyLinkedList<T>& const cluster_linked_list, DoublyLinkedList<T>& const doubly_linked_list);

	//200213 Use vector to instead pointer for time series. Add template
	template<typename T, typename Y, typename U, typename U1, typename U2>
	void get_APLA_point_new_y_projection(U& input_argument, const vector<T>& const original_time_series_vector, typename TOOL::Y_PROJECTION_ARGUMENT& const y_projection_argument, DoublyLinkedList<Y>& const all_linked_list, DoublyLinkedList<Y>& const cluster_linked_list, DoublyLinkedList<U2>& const doubly_linked_list, U1& output_argument);

	template<typename T>
	vector<T>& get_apla_projection(const vector<DataType>& const query_time_series_vector, const DoublyLinkedList<T>& const doubly_linked_list, vector<T>& const area_vector);
	//191114 compue distance LB for two segments.
	template<typename T>
	double get_segment_distance_LB(const T& const segment_1, const T& const segment_2);
	//191114 compue distance LB for two APLA points.
	template<typename T>
	double get_distance_LB(const DoublyLinkedList<T>& const doubly_linked_list, const vector<T>& const area_vector);
	//191114 directly get distanceLB from query time series & original APLA point, combine(get_apla_projection, get_segment_distance_LB, get_distance_LB)
	template<typename T>
	double get_distance_LB_by_series_apla(const vector<DataType>& const query_time_series_vector, const DoublyLinkedList<T>& const doubly_linked_list);

	//191129 compue PLA distance LB for two APLA points.
	template<typename T>
	double get_distance_LB_pla(const DoublyLinkedList<T>& const query_linked_list, const DoublyLinkedList<T>& const doubly_linked_list);

	//200109 compue query time series with reconstruct time series. APCA distance AE
	template<typename T>
	double get_distance_AE(const vector<DataType>& const query_time_series_vector, const DoublyLinkedList<T>& const doubly_linked_list);

	//191205 compue PLA distance LB for two APLA points fast.
	template<typename T>
	double get_distance_LB_pla_speed(const DoublyLinkedList<T>& const query_linked_list, const DoublyLinkedList<T>& const doubly_linked_list);

	//200110 compue PLA distance LB for two APLA points lower bound.
	template<typename T>
	double get_distance_lower_bound(const DoublyLinkedList<T>& const query_linked_list, const DoublyLinkedList<T>& const doubly_linked_list);

	//200113 compute adaptive segment endpoint distance
	template<typename T>
	double get_apla_endpoint_distance(const DoublyLinkedList<T>& const query_linked_list, const DoublyLinkedList<T>& const doubly_linked_list);

	//200107 compue distance LB for two APLA points by APCA(average value).
	template<typename T>
	double get_distance_apca_LB(const DoublyLinkedList<T>& const doubly_linked_list, const vector<T>& const area_vector);
	//For struct Length 190623
	void getRandomEndpoint(const int& const time_series_length, const int& const endpoint_size, set<int>& const result_endpoint_set);//190623 get random end point

	void getBestEndpoint(typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190621 compare all situation, find best result

	void combinationUtilJumpVector(DataType*& const original_time_series, const vector<int>& const arr, int index, vector<int> data, int i, set<pair<double, vector<int>>>& const endpoint_collection);//190624 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	void getAllAPLAResult(DataType*& const original_time_series, const vector<int>& const arr, const int& const r, set<pair<double, vector<int>>>& const endpoint_collection);//190624 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()

	void combinationUtilJumpVector(DataType*& const original_time_series, const vector<int>& const arr, int index, vector<int> data, int i, pair<double, vector<AREA_COEFFICIENT>>& const best_endpoint_pair);//190626 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	//200811 use vector time series
	template<typename T, typename Y, typename U>
	void combinationUtilJumpVector(const vector<T>& const original_time_series_vector, const vector<U>& const arr, U index, vector<U> data, U i, pair<long double, vector<Y>>& const best_endpoint_pair);//200811 No adjacent combination /* arr[] ---> Input Array data[] ---> Temporary array to store current combination start & end ---> Staring and Ending indexes in arr[] index ---> Current index in data[] r ---> Size of a combination to be printed */

	void getBestAPLAResult(DataType*& const original_time_series, const vector<int>& const arr, const int& const r, pair<double, vector<AREA_COEFFICIENT>>& const best_endpoint_pair);//190626 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()
	
	//200811 uset vector time series
	template<typename T, typename Y, typename U>
	void getBestAPLAResult(const vector<T>& const original_time_series_vector, const vector<U>& const arr, const U& const r, pair<long double, vector<Y>>& const best_endpoint_pair);//190626 // No adjacent combinationThe main function that prints all combinations of size r in arr[] of size n. This function mainly uses combinationUtil()

	void getAPLAArea5LowerBoundSpeedNoErase0(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190606 Failed, cannot speed up when loop from head and tail to middle.
	void getAPLAArea5LowerBoundSpeedNoEraseArray(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190416
	double getLineSegmentTriangleAreaImprove(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190604

	void getAPLADeviation(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType*& const original_time_series);//190115 try compare with deviation

	double getAPLASumDeviaitonNoAB(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector);//190624 get sum_deviation very fast, No a,b

	template<typename T, typename Y>
	long double getAPLASumDeviaitonNoAB(const vector<T>& const original_time_series_vector, vector<Y>& const area_vector);//200811 get sum_deviation very fast, No a,b

	double getAPLASumDeviaiton(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector);//190702 get sum_deviation very fast, already has a&b

	void getAPLASegmentDeviation(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector);//190626 Aleady has a&b, right_endpoint, width, sum_value
	//190826 Linked list
	void getAPLASegmentDeviation(DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//190826 Aleady has a&b, right_endpoint, width, sum_value

	//200303 Add template & use vector time series. Get segment sun deviation
	template<typename T, typename Y>
	void getAPLASegmentDeviation(const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const doubly_linked_list);//190826 Aleady has a&b, right_endpoint, width, sum_value

	//200303 get one segment  (sum difference)^2  //getSegmentSumDifferenceSquare
	template<typename T, typename Y, typename U>
	long double get_segment_difference_square(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient, U& const max_deviation_seg);

	//200303 get sume deviation of apprxomation, No a&b of each segment.  Only know right endpoint and rectangle width
	template<typename T, typename Y>
	long double get_sum_deviation_no_ab(const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const doubly_linked_list);

	//200327 get sume deviation & max deviation of apprxomation, No a&b of each segment.  Only know right endpoint and rectangle width
	template<typename T, typename Y, typename U>
	long double get_sum_deviation_no_ab(const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const doubly_linked_list, U& const result_collection);

	//201005 get sume deviation of apprxomation, No a&b of each segment.  Only know right endpoint and rectangle width
	template<typename T, typename Y>
	long double get_sum_deviation_no_ab(const vector<T>& const original_time_series_vector, const vector<Y>& const approximation_vector);
	//210304
	template<typename T>
	long double get_sum_upper_bound(const DoublyLinkedList<T>& const doubly_linked_list);

	void getAPLAReconstructSeriesNoAB(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, vector<double>& const reconstruct_time_series);//190625 No a,b, only has endpoint and witdh get reconstruct series vector

	//20303 Use template
	template<typename T, typename Y>
	void getAPLAReconstructSeries(const vector<Y>& const area_vector, vector<T>& const reconstruct_time_series);//190625 Has endpoint, width, a,b. get reconstruct series vector

	//20303 Use template
	template<typename T, typename Y>
	void getAPLAReconstructSeries(const DoublyLinkedList<Y>& const doubly_linked_list, vector<T>& const reconstruct_time_series);//190917 Has endpoint, width, a,b. get reconstruct series vector

	//191230 get reconstruct time series for one time series, not know a&b
	template<typename T>
	void get_PLA_reconstruct_series(const vector<T>& const reconstruct_time_series);

	//200108 get reconstruct time series for one time series
	template<typename T>
	void get_PLA_reconstruct_series(const vector<T>& const original_time_series_vector, vector<T>& const reconstruct_time_series_vector);

	//200206 get econstruct time series, a, b, segment_width for one time series
	template<typename T, typename Y>
	void get_PLA_reconstruct_series(const vector<T>& const original_time_series_vector, Y& const temp_coefficient, vector<T>& const reconstruct_time_series_vector);

	//200205 get part reconstruct time series for part of time series
	template<typename T>
	void get_PLA_part_reconstruct_series(const vector<T>& const original_time_series, const int& const begin_id, const int& const end_id, vector<T>& const part_reconstruct_time_series);

	//200206 get part reconstruct time series, a&b width. for part of time series
	template<typename T, typename Y>
	void get_PLA_part_reconstruct_series(const vector<T>& const original_time_series, const int& const begin_id, const int& const end_id, Y& const temp_coefficient, vector<T>& const part_reconstruct_time_series);

	//200205 get part reconstruct time series for part of time series
	template<typename T, typename Y>
	void get_PLA_reconstruct_series_by_endpoint(const vector<T>& const original_time_series, vector<Y>& const right_endpoint, vector<T>& const reconstruct_time_series);

	//200220 get reconstruction time series for one segemnt, already know a&b and right endpoint
	template<typename T, typename Y>
	void get_PLA_reconstruct_series_by_segment(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient, vector<T>& const reconstruct_time_series);
	//200110 get a & b of original time series
	template<typename T>
	void get_PLA_coefficient(const vector<T>& const original_time_series_vector, typename TOOL::APLA_COEFFICIENT& const pla_coefficient);

	//200110 get a & b of original time series
	template<typename T>
	typename APLA::APLA_COEFFICIENT get_PLA_coefficient(const vector<T>& const original_time_series_vector);

	void getAPLA(const typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector);//181213

	//200302 use template
	template<typename T, typename Y, typename U>
	void getAPLA(const U& const input_argument, T*& const original_time_series, DoublyLinkedList<Y>& const doubly_linked_list);//190917

	void getAPLA(const typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, list<AREA_COEFFICIENT>& const area_vector);//190408

	void getAPLAArray(const typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, SEGMENT_COEFFICIENT*& const area_vector);//190426

	//get slope and intercept a&b and sum value
	void getPLAByAdaptiveSegment(const typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector);//190617

	//190821 Linked list to instead vector, get slope and intercept a&b and sum value . No minmax point
	void getPLAByAdaptiveSegment(const typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//190821

	//200212  Linked list to instead vector, get slope and intercept a&b and sum value . No minmax point
	//200212 Linked list to instead vector, get slope and intercept a&b and sum value . No minmax point
	template<typename T, typename Y, typename U>
	void getPLAByAdaptiveSegment(const U& const input_argument, const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const doubly_linked_list);//200212

	//191031 Linked list to instead vector, get slope and intercept a&b and sum value & minmax point of every segment
	//template<typename T>
	void getPLAByAdaptiveSegmentAndMinMax(const typename TOOL::INPUT_ARGUMENT& const input_argument, DataType*& const original_time_series, DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list);//191031

	//200212 use vector to instead pointer for time series.Add template
	//191031 Linked list to instead vector, get slope and intercept a&b and sum value & minmax point of every segment
	template<typename T, typename Y, typename U>
	void getPLAByAdaptiveSegmentAndMinMax(const U& const input_argument, const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const doubly_linked_list);//191031

	//200212 use vector to instead pointer for time series.Add template
	//200319 Linked list to instead vector, get slope and intercept a&b, triangle density, area difference
	template<typename T, typename Y>
	void get_each_segment_ab_density_difference(const vector<T>& const original_time_series_vector, DoublyLinkedList<Y>& const doubly_linked_list);

	void printSpecialPoint(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190517
	void printSpecialPoint0(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190524
	void printSpecialPoint1(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190602

	double computeMagnitude(AREA_COEFFICIENT& const temp_coefficient, double difference);//190521
	double computeDifferenceA(AREA_COEFFICIENT& const temp_coefficient, int level);//190524
	AREA_COEFFICIENT& getMaxMagnitude0(AREA_COEFFICIENT& const temp_coefficient);//190521
	AREA_COEFFICIENT& getMaxMagnitude(AREA_COEFFICIENT& const temp_coefficient);//190521
	bool printMaxMagnitude(const AREA_COEFFICIENT& const temp_coefficient);//190521
	double printReconstructRecursive(const AREA_COEFFICIENT& const temp_coefficient, int  level);//190522

	//200907 print right endpoint, segment width, a&b, increment area, area difference, long segment triangle density
	template<typename T>
	void print_segment_coefficients(const T& const temp_coefficient);
	void testAChange(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190522
	void testAChange0(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient, int& const left_id, int& const right_id);//190522
	void testAChange1(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190522 max_deviation point in minmax point, endpoint, sub_deviation point
	void draw3Line(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);
	void testDeviationPointRecursive(AREA_COEFFICIENT& const temp_coefficient, int level, set <int>& deviation_id);//190524
	bool testDeviationPoint(AREA_COEFFICIENT& const temp_coefficient, int level, set <int>& deviation_id);//190524

	template<typename T>
	bool compareLinkListAndVector(DoublyLinkedList<T>& const doubly_linked_list, vector<T>& const area_vector);

	//210301
	template<typename T, typename Y>
	bool assert_adjacent_nodes(const vector<T>& const original_time_series_vector, const DoublyListNode<Y>& const node);

	//191101 assert segment & sub segment: non INF, right endpoint, width, minmax point,, a&b
	template<typename T, typename Y>
	bool assertLinkedListAndSubLinkedList(const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const doubly_linked_list);
	// 191101 assert segment: non INF, right endpoint, width, minmax point,, a&b
	template<typename T>
	bool assertLinkedList(const DoublyLinkedList<T>& const doubly_linked_list);

	// 201004 assert segment: non INF, right endpoint, width, minmax point, a&b, sub right endpoint
	template<typename T, typename Y>
	bool assertLinkedList(const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const doubly_linked_list);

	// 191101 assert segment right endpoint, width
	template<typename T>
	bool assertRightEndpoint_Width(const DoublyLinkedList<T>& const doubly_linked_list) const;
	
	template<typename T>
	bool assert_a_b(const DoublyLinkedList<T>& const doubly_linked_list) const;

	// 200214 assert a & b of  one  segment
	template<typename T>
	bool assert_segment_a_b(const T& const temp_coefficient) const;

	//200512 assert a&b of one segment with time series
	template<typename T, typename Y>
	inline bool assert_segment_a_b(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient);

	// 200320 assert compare a & b of two segments
	template<typename T>
	inline bool assert_two_segments_a_b(const T& const segment1, const T& const segment2) const;

	// 200713 assert a & b of two segments is equal
	template<typename T>
	inline bool assert_equal_a_b(const T& const segment1, const T& const segment2) const;

	// 200712 whether a&b of two segments are same
	template<typename T>
	inline bool if_equal_ab(const T& const segment1, const T& const segment2) const;
	
	//210402
	template<typename T>
	inline const bool if_similar_ab_two_segments(const T& const segment1, const T& const segment2) const;

	// test loop linked list
	template<typename T>
	void evaluate_linkedlist_for_loop(const DoublyLinkedList<T>& const doubly_linked_list);

	//191211 evaluate linked list right endpoint, width, right height difference
	template<typename T>
	void assert_linkedlist_rightEndpoint_Width_rightHeightDifference(const typename TOOL::INPUT_ARGUMENT& const input_argument, const DoublyLinkedList<T>& const doubly_linked_list);

	//191211 evaluate linked list right endpoint, width, right height difference
	template<typename T>
	void assert_split_coefficients(const T& const temp_coefficient);

	//200311 evaluate structure split coefficients, sub left & right, segment density. magnitude
	template<typename T>
	void assert_structure_split_coefficients(const T& const split_coeffcients);
	
	template<typename T>
	bool assert_segment_minmax(const T& const segment) const;

	// 200512 assert min max point of one segment
	template<typename T, typename Y>
	inline bool assert_segment_minmax(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient) const;

	//200313 assert min&max value
	template<typename T, typename Y>
	void assert_minmax_value(const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const doubly_linked_list);
	template<typename T, typename Y>
	inline void assert_segment_a_b_minmax(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient);
	
	template<typename T, typename Y, typename U>
	bool assert_loop_segment_density_area_difference(const DoublyLinkedList<T>& const doubly_linked_list, const multimap<U, Y, std::greater<U>>& const multi_map_1, const multimap<U, Y, std::greater<U>>& const multi_map_2);
	
	template<typename T>
	bool assert_split_map_area_difference(const T& const split_map);

	// 201227 Assert Split MAP
	template<typename T>
	bool assert_split_map_optimization(const T& const split_map);
	
	template<typename T>
	bool assert_merge_map_triangle_density(const T& const merge_map);
	
	template<typename T, typename Y>
	void evaluate_accuracy_time_split_point(const vector<T>& const original_time_series_vector, const Y& const temp_coefficient, const double& const split_id_min_density, const double& const split_time_min_density);
	
	//201028 evaluate 1 upper & lower bound: max devation * l > sum deviation > lower bound deviation, 2 The probablity of left&right, min&max point is max deviation poiont
	template<typename T, typename Y, typename U>
	bool evaluate_segment_upper_lower_bound(U& const input_argument, const vector<T>& const original_time_series_vector, const Y& const temp_coefficient);

	//201023 evaluate 1 upper & lower bound: max devation * l > sum deviation > lower bound deviation, 2 The probablity of left&right, min&max point is max deviation poiont
	template<typename T, typename Y, typename U>
	bool evaluate_upper_lower_bound(U& const input_argument, const vector<T>& const original_time_series_vector, const DoublyLinkedList<Y>& const doubly_linked_list);

	//void computeDft(const vector<double> &inreal, const vector<double> &inimag, vector<double> &outreal, vector<double> &outimag);
	//void computeDft(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType original_time_series[], DataType*& dft_time_series);

	inline typename GEOMETRY::POINT& getIntersectionPoint(const AREA_COEFFICIENT& const temp_coefficient1, const AREA_COEFFICIENT& const temp_coefficient2, typename GEOMETRY::POINT& const intersection_point);//190527
	inline typename GEOMETRY::POINT& getIntersectionPoint(const APLA_COEFFICIENT& const apla1, const APLA_COEFFICIENT& const apla2, typename GEOMETRY::POINT& const intersection_point);//190528
	inline typename GEOMETRY::POINT& getIntersectionPoint(const double& const a1, const double& const b1, const double& const a2, const double& const b2, typename GEOMETRY::POINT& const intersection_point);//190528
	inline typename GEOMETRY::POINT& segmentsIntr(typename  GEOMETRY::POINT& const point1, typename  GEOMETRY::POINT& const point2, typename  GEOMETRY::POINT& const point3, typename GEOMETRY::POINT& const point4, typename GEOMETRY::POINT& const intersection_point);//190528
	//210402
	inline typename GEOMETRY::POINT& segmentsIntrArea(typename  GEOMETRY::POINT& const point1, typename  GEOMETRY::POINT& const point2, typename  GEOMETRY::POINT& const point3, typename GEOMETRY::POINT& const point4, typename GEOMETRY::POINT& const intersection_point);//190528 //
	inline typename GEOMETRY::POINT& getLineSegmentIntersectionPoint(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient, typename GEOMETRY::POINT& const intersection_point1, typename GEOMETRY::POINT& const intersection_point2);//190528

	double getLineSegmentTriangleArea(DataType*& const original_time_series, vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190528
	//double getLineSegmentLeftArea(typename  GEOMETRY::POINT& const up_point, typename  GEOMETRY::POINT& const down_point, typename  GEOMETRY::POINT& const intersection_point);//190604
	//double getLineSegmentRightArea(typename  GEOMETRY::POINT& const intersection_point, typename  GEOMETRY::POINT& const up_point, typename  GEOMETRY::POINT& const down_point);//190604

	inline double getLineSegmentTriangleDensity(AREA_COEFFICIENT& const temp_coefficient);//190602

	double getLineSegmentTriangleAreaDensity(vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190611
	//190826 Use linkned list to instead vector
	double getLineSegmentTriangleAreaDensity(DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190826

	//190905 Use linkned list to instead vector
	inline double& getLineSegmentTriangleAreaDensity(int nothing, const AREA_COEFFICIENT& const left_segment, const AREA_COEFFICIENT& const right_segment, AREA_COEFFICIENT& const merged_segment);//190905

	//210402 get density of triangle area
	template<typename T, typename Y>
	inline double getLineSegmentTriangleAreaDensity(const T& const left_segment, const  T& const right_segment, Y& const merged_segment);

	//200728 get density of triangle area from current segment & accumulates segment
	template<typename T>
	inline long double get_triangle_density_by_accumulate_segment(T& const current_segment, const T& const accumulate_segment);

	//201007 get triangle area by line A1A2 and B1B2
	template<typename T, typename Y>
	inline long double get_triangle_area_by_points(T& const current_segment, const Y& const B1_y, const Y& const B2_y);


	//200206 get endpoint of long segment and two short segments, and comput height difference between them
	template<typename T>
	inline long double get_line_segment_height_diference(const T& const left_segment, const  T& const right_segment, const T& const merged_segment);

	//210402
	template<typename T, typename Y, typename U>
	inline long double get_line_segment_height_diference_from_accumulate(const vector<Y>& const original_time_series_vector, T& const current_segment, T& const accumulate_segment, U& const output_argument);

	//210402 get upper bound of sub left & right segment after splitting
	template<typename T, typename Y, typename U>
	inline void get_upper_bound_sub_segment_left_right(const vector<Y>& const original_time_series_vector, T& const left_segment, T& const right_segment, const T& const merged_segment, U& const output_argument);

	//210402 get endpoint of long segment and two short segments, compute height difference with original time series
	template<typename T, typename Y, typename U>
	inline long double get_line_segment_height_diference_with_original(const vector<Y>& const original_time_series_vector, const T& const left_segment, const  T& const right_segment, T& const merged_segment, U& const output_argument);
	
	//210114 After optimization. get upper bound when increase/decrease, left/right endpoint
	template<typename T, typename Y, typename U>
	inline long double get_upper_bound_move_endpoint(const vector<Y>& const original_time_series_vector, const T& const segment_original, T& const segment_move, U& const output_argument);
	
	//210121
	template<typename T, typename Y>
	inline void count_upper_bound_type(const T upper_bound_type, Y& const map_upper_bound);

	//210122
	template<typename T, typename Y>
	void merge_map_with_same_key(const T& const map_target, Y& const map_result);

	//210121
	template<typename T, typename Y>
	inline void count_upper_bound_whole(const T& const output_argument, Y& const evaluation_bound_whole);

	//210122
	template<typename T>
	void print_count_of_upper_bound(const T& const evaluation_bound);

	//210126
	template<typename T, typename Y>
	inline void get_current_max_id(const T& const map_result, Y& const temp_coefficient);

	//210115 Not INF
	template<typename T>
	inline bool assert_segment_bound(const T& const temp_coefficient);
	//210115
	template<typename T, typename Y>
	inline bool assert_segment_bound(const vector<Y>& const original_time_series_vector, const T& const temp_coefficient);

	//210115
	template<typename T, typename Y>
	inline bool assert_bound(const vector<Y>& const original_time_series_vector, const DoublyLinkedList<T>& const doubly_linked_list);

	inline double& getPLAMinMaxAreaDifference(vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190724

	// 190826
	// Input a&b, min max point, widht, right endpoint
	// Output PLA, minmax width area difference
	double& getPLAMinMaxAreaDifference(DoublyLinkedList<AREA_COEFFICIENT>& const doubly_linked_list, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190826

	// 190905
	// Input a&b, min max point, widht, right endpoint
	// Output PLA, minmax width area difference
	double& getPLAMinMaxAreaDifference(const AREA_COEFFICIENT& const left_segment, const AREA_COEFFICIENT& const right_segment, AREA_COEFFICIENT& const merged_segment);//190905


	double get3SubSegment2PointsABArea(vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, AREA_COEFFICIENT& const temp_coefficient);//190628 get a b of sub segment in 2 endpoints and middle.

	template<typename T>
	void getPLALineRectangleAreaDensity(vector<T>& const area_vector);//190627 Already know a&b right_endpoint and width

	template<typename T>
	void getPLALineRectangleAreaDensity(DoublyLinkedList<T>& const doubly_linked_list);//190826 Already know a&b right_endpoint and width

	//200419 scan segment to get min & max value height
	template<typename T, typename Y>
	long double get_segment_minmax_height(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	//200519 initialize min&max point of segmnet = INF
	template<typename T>
	inline void initialize_segment_minmax(T& const temp_coefficient);

	//200422 get min&max point of one point segment
	template<typename T, typename Y>
	inline void get_segment_one_point_minmax(const T& const point_value, Y& const temp_coefficient);

	//200422 get min&max point of two points segment
	template<typename T, typename Y>
	inline void get_segment_two_point_minmax(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	//200714 get min&max point of two points segment
	template<typename T, typename Y>
	inline void get_segment_two_point_minmax(const T& const left_value, const T& const right_value, Y& const temp_coefficient);

	//200422 get and verify min&max point of long segment
	template<typename T, typename Y, typename U>
	inline void get_segment_minmax_point(const T& const point_id, const Y& const point_value, U& const temp_coefficient);

	//scan to get minmax point of segment
	void getSegmentMinMaxPoint(DataType*& const original_time_series, AREA_COEFFICIENT& const temp_coefficient);//190701

	//200212 scan segment. time series from pointer to vector, Addd template
	template<typename T, typename Y>
	void getSegmentMinMaxPoint(const vector<T>& const original_time_series_vector, Y& const temp_coefficient);

	//200714 get accumulated segment minmaxpoint by last short segment minmax point
	template<typename T, typename Y>
	inline void get_minmax_segment_by_accumulation(const T& const end_point_value,const Y& const short_segment, Y& const temp_coefficient);

	//200731 get accumulated segment minmaxpoint by accumulate endpoint value
	template<typename T, typename Y>
	inline void get_minmax_segment_by_accumulation(const T& const end_point_value, Y& const accumulate_segment);

	//200506 Get min&max point of long segment from left & right short segments
	template<typename T>
	inline void get_long_segment_minmax_by_sub_segments(const T& const sub_left_segment, const T& const sub_right_segment, T& const long_segment);

	//191031
	//Input: oritinal time series. sub left right: 1 right enpoint, 2 width. Long segment: 1 right endpoint, 2 width, 3 min&max point
	//Output: sub left right: min&max point 
	void getSubMinMaxPoint(DataType*& const original_time_series, AREA_COEFFICIENT& const left_segment, AREA_COEFFICIENT& const right_segment, const AREA_COEFFICIENT& const merged_segment);

	//200213 use vector to instead pointer for time series
	//Input: oritinal time series vector. sub left right: 1 right enpoint, 2 width. Long segment: 1 right endpoint, 2 width, 3 min&max point
	//Output: sub left right: min&max point 
	template<typename T, typename Y>
	void getSubMinMaxPoint(const vector<T>& const original_time_series_vector, Y& const left_segment, Y& const right_segment, const Y& const merged_segment);

	void getSegmentADifference(vector<AREA_COEFFICIENT>& const area_vector);//190627 Get difference of left and right endpoint a, already know a&b

	void countBreakPoint(vector<AREA_COEFFICIENT>& const area_vector, const int& const vector_id, vector<int>& const count_break_point_vector);//190613

	template<typename T>
	void write_reconstruct_time_series(const string& const file_name, const DoublyLinkedList<T>& const doubly_linked_list);
	template<typename T>
	void write_segment_point(const string& const file_name, const DoublyLinkedList<T>& const doubly_linked_list);
};

#include "CAPLA.cpp"

#endif