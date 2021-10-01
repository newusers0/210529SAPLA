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

	void evaluate_multi_KNN_speed();
	template<typename T, typename Y>
	inline void print_input_output(const T& const input_argument, const Y& const output_argument);
};


TEMPLATE
void Evaluation::evaluate_multi_KNN_speed() {

	string* chebyshev_write_pointer = nullptr;
	string dataset_attribute[] = { "", "Single dataset", "Mixed dataset", "Multi dimension data set" };

	int query_time_series_id = INF;
	int K = INF;
	long double number_point_max_deviation_true = 0;
	long double number_point_max_deviation_false = 0;
	long double number_not_smaller_than_sum_deviation = 0;
	long double number_smaller_than_sum_deviation = 0;
	long double number_not_smaller_than_sum_deviation_pow = 0;
	long double number_smaller_than_sum_deviation_pow = 0;
	int initial_file_id = INF;
	int file_id = INF; 
	vector<size_t> file_id_choose_vector;
	vector<size_t> data_list;
	size_t representation_option = INF;
	const vector<int> representation_option_vector{ int(INF), 2, 3, 4, 5, 6, 7, 8, 9 };
	const vector<string> name_representation_vector_string_all_vector{ "MSPLA", "PLA", "APCA", "PAA", "Chebyshev", "APLA", "PAALM", "SAPLA", "SAX" };
	vector<string> name_representation_vector_string_vector;
	for (int i = 1; i < representation_option_vector.size(); i++) {
		name_representation_vector_string_vector.emplace_back(name_representation_vector_string_all_vector[representation_option_vector[i] - 1]);
	}
	int N = INF;
	int initial_N = INF;
	int final_N = INF;
	int initial_K = INF;
	int final_K = INF;
	int data_final_dimension = INF;
	int arity_d0 = INF;
	int max_node0 = INF;
	int file_total_number = INF;
	int n0 = INF;
	int point_number0 = INF;

	const int& const evaluation_type = 0;
	const size_t& const data_type = 1;
	const int option_tree = 1;
	string data_file_name = "./191202DataSet/data_set.txt";
	string* file_address_pointer = nullptr;
	const string& const str_suffix = "1";
	string notice_string = ".*Notice: suffix: ";
	const string& const str_33 = "./191120KNNMatlab/EachFileMethodPrunePowerCombine" + str_suffix;
	const string& const str_27 = "./191120KNNMatlab/EachFileMethodPrunePower" + str_suffix;
	const string& const str_28 = "./191120KNNMatlab/EachFileMethodSumDeviation" + str_suffix;
	const string& const str_29 = "./191120KNNMatlab/EachFileMethodAccuracy" + str_suffix;
	const string& const str_30 = "./191120KNNMatlab/EachFileMethodKNNRunTime" + str_suffix;
	const string& const str_31 = "./191120KNNMatlab/EachFileMethodApproximationTime" + str_suffix;
	const string& const str_32 = "./191120KNNMatlab/EachFileMethodOnlyKNNTime" + str_suffix;

	
	const int option_homogenous_data_type = 5;
	switch (data_type) {
	case 0:  
		assert(0);
		break;
	case 1:

		n0 = 512;
		initial_N = 12;
		final_N = 192;
		point_number0 = 20;
		final_K = 16;

		max_node0 = 3;

		switch (option_homogenous_data_type) {
		case 0: 
			file_total_number = 21; 

			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;

			break;
		case 1: 
			file_total_number = 18;
			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;

			break;
		case 2: 
			file_total_number = 24;
			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			break;
		case 3: 
			file_total_number = 41;
			
			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;

			break;
		case 4: 
			file_total_number = 21;
			n0 = 1024;
			final_N = 384;
			point_number0 = 50;
			final_K = 32;

			file_id_choose_vector.resize(file_total_number);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);

			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;

			break;
		case 5: 

			file_id_choose_vector.resize(20);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 1);

			file_total_number = file_id_choose_vector.size();
			n0 = 1024;
			final_N = 192;
			point_number0 = 100;
			final_K = 64;
			max_node0 = 3;

			query_time_series_id = TOOL::get_random_max(point_number0);
			query_time_series_id = point_number0 - 5;

			break;
		
		default:
			assert(0);
		}

		data_final_dimension = 1;
		assert(data_final_dimension == 1);
		file_address_pointer = nullptr;
		file_address_pointer = &data_file_name;
		for (int list_id = 1; list_id < file_total_number + 1; list_id++) { data_list.emplace_back(list_id); }
		break;
	
	default:
		assert(0);
		break;
	}
	initial_file_id = 0;
	const int final_file_id = file_id_choose_vector.size();// 
	initial_K = 2;
	if (option_homogenous_data_type == 6) initial_K = 2;

	const int method_total_number = representation_option_vector.size() - 1;
	
	vector<size_t> N_coefficient_vector;
	vector<int> K_coefficient_vector;

	vector<long double> total_prune_power_combine_vector;
	vector<long double> total_prune_power_vector;
	vector<long double> total_sum_deviation_vector;//191206
	vector<long double> total_max_deviation_vector;//191206
	vector<long double> total_max_deviation_av_vector;//210910
	vector<long double> total_max_width_deviation_vector;//191206
	vector<long double> total_accuracy_vector;//191204
	vector<long double> total_run_time_vector;//191119
	vector<long double> total_run_time_has_IO_vector;//210606
	vector<long double> total_approximation_time_vector;//200108
	vector<long double> total_build_tree_time_vector;//210906
	vector<long double> total_knn_time_vector;//200108
	vector<long double> total_knn_time_has_IO_vector;//210606

	vector<long double> method_sum_prune_power_combine_vector(method_total_number, 0);
	vector<long double> method_sum_prune_power_vector(method_total_number, 0);
	vector<long double> method_sum_sum_deviation_vector(method_total_number, 0);
	vector<long double> method_sum_max_deviation_vector(method_total_number, 0);
	vector<long double> method_sum_max_deviation_av_vector(method_total_number, 0);
	vector<long double> method_sum_max_width_deviation_vector(method_total_number, 0);
	vector<long double> method_sum_accuracy_vector(method_total_number, 0);
	vector<long double> method_run_time_vector(method_total_number, 0);
	vector<long double> method_run_time_has_IO_vector(method_total_number, 0);
	vector<long double> method_approximation_time_vector(method_total_number, 0);
	vector<long double> method_build_tree_time_vector(method_total_number, 0);
	vector<long double> method_knn_time_vector(method_total_number, 0);
	vector<long double> method_knn_time_has_IO_vector(method_total_number, 0);
	vector<long double> one_file_prune_power_combine_vector(method_total_number, 0);
	vector<long double> one_file_prune_power_vector(method_total_number, 0);
	vector<long double> one_file_sum_deviation_vector(method_total_number, 0);
	vector<long double> one_file_max_deviation_vector(method_total_number, 0);
	vector<long double> one_file_max_width_deviation_vector(method_total_number, 0);
	vector<long double> one_file_accuracy_vector(method_total_number, 0);
	vector<long double> one_file_run_time_vector(method_total_number, 0);
	vector<long double> one_file_approximation_time_vector(method_total_number, 0);
	vector<long double> one_file_knn_time_vector(method_total_number, 0);
	vector<long double> one_file_knn_time_has_IO_vector(method_total_number, 0);
	vector<long double> method_file_prune_power_combine_vector(size_t(method_total_number * file_total_number), 0.0);
	vector<long double> method_file_prune_power_vector(size_t(method_total_number * file_total_number), 0.0);
	vector<long double> method_file_sum_deviation_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_max_deviation_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_max_deviation_av_vector(method_total_number* file_total_number, 0.0);
	vector<long double> method_file_max_width_deviation_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_accuracy_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_run_time_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_run_time_has_IO_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_approximation_time_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_build_tree_time_vector(method_total_number* file_total_number, 0.0);
	vector<long double> method_file_knn_time_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_knn_time_has_IO_vector(method_total_number * file_total_number, 0.0);

	const int split_methods_number = 1;
	vector<long double> local_total_split_id_sum_deviation(split_methods_number, 0);
	vector<long double> local_total_split_id_shift(split_methods_number, 0);
	vector<long double> local_total_split_id_time(split_methods_number, 0);
	vector<long double> global_total_knn_prune_power(split_methods_number, 0);
	vector<long double> global_total_approximation_sum_deviation(split_methods_number, 0);
	vector<long double> global_total_approximation_time(split_methods_number, 0);

	const double N_size = 5;
	const double final_initial_N_number = 30;
	vector<long double> initial_number_coefficients_vector{ -2, -1 };
	for (double vector_id = 0; initial_number_coefficients_vector.size() < final_initial_N_number; vector_id += 3) {
		initial_number_coefficients_vector.emplace_back(vector_id);
	}
	assert(initial_number_coefficients_vector.size() == final_initial_N_number);
	vector<long double> approximation_initial_N_vector;
	vector<long double> total_initial_N_prune_power_vector(final_initial_N_number, 0);
	vector<long double> total_initial_N_sum_deviation_vector(final_initial_N_number, 0);
	vector<long double> total_initial_N_run_time_vector(final_initial_N_number, 0);
	vector<long double> total_initial_N_approximation_time_vector(final_initial_N_number, 0);
	vector<long double> total_initial_N_knn_time_vector(final_initial_N_number, 0);
	vector<long double> initial_N_by_N_prune_power_vector(final_initial_N_number * N_size, 0);
	vector<long double> initial_N_by_N_sum_deviation_vector(final_initial_N_number * N_size, 0);
	vector<long double> initial_N_by_N_run_time_vector(final_initial_N_number * N_size, 0);
	vector<long double> initial_N_by_N_approximation_time_vector(final_initial_N_number * N_size, 0);
	vector<long double> initial_N_by_N_knn_time_vector(final_initial_N_number * N_size, 0);
	typename TOOL::EVALUATION_BOUND evaluation_bound_whole;
	typename TOOL::DATA_SOURCE data_source(data_type, data_final_dimension, data_list, file_total_number, point_number0, n0, *file_address_pointer);
	data_source.option_has_burst_data = option_homogenous_data_type;
	TOOL::initial_data_source(data_source);
	data_source.bigger_account = 0;
	assert(data_source.point_number != INF && data_source.single_time_series_length != INF && data_source.multi_single_time_series_length != INF);

	for (file_id = initial_file_id; file_id < final_file_id; file_id++) {
		cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%   EVALUATION    file id: " << file_id_choose_vector[file_id] + 1 << "   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
		bool change_file = true;
		if (data_source.data_type == 1) {
			TOOL::get_normalize_write_one_file(data_source, file_id_choose_vector[file_id]);
		}

		bool print_each_result = false;
		if (print_each_result == true) {
			vector<DataType> original_time_series_vector_print(data_source.multi_single_time_series_length, INF);
			TOOL::get_read_multi_file_address(data_source, file_id_choose_vector[file_id]);//210603
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

		vector<DataType> query_time_series_vector(data_source.multi_single_time_series_length, INF);
		DataType* query_time_series = new DataType[data_source.multi_single_time_series_length];
		TOOL::get_read_multi_file_address(data_source, file_id_choose_vector[file_id]);//210603
		//get data from new file. already normalized in advanced, no miss first point.
		TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, query_time_series_id, query_time_series_vector);
		copy_n(query_time_series_vector.begin(), query_time_series_vector.size(), query_time_series);

		std::multiset<pair<double, int>> squential_scan_result_set;
		TOOL::SimpleBaseKNNSearchMulti(data_source, query_time_series_vector, squential_scan_result_set);

		vector<TOOL::Y_PROJECTION_ARGUMENT> multi_y_projection_argument_vector(data_source.point_number, TOOL::Y_PROJECTION_ARGUMENT());
		vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>> multi_all_linked_list(data_source.point_number, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>());//191030
		vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>> multi_cluster_linked_list(data_source.point_number, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED>());//191030
		for (N = initial_N; N <= final_N; N *= 2) {

			if (file_id == initial_file_id) { N_coefficient_vector.emplace_back(N); }

			if (evaluation_type == 0) {

				for (representation_option = 1; representation_option <= method_total_number; representation_option++) {

					for (int split_method_id = 0; split_method_id < split_methods_number; split_method_id++) {
						
						MULTI multi(data_source, n0, N, INF, file_id_choose_vector[file_id], change_file, data_source.time_series_dimension, data_source.point_number, query_time_series_id, max_node0, initial_K, representation_option_vector[representation_option], file_address_pointer, chebyshev_write_pointer, option_tree);

						multi.input_argument.number_point_max_deviation_true = 0;
						multi.input_argument.number_point_max_deviation_false = 0;
						multi.input_argument.number_not_smaller_than_sum_deviation = 0;
						multi.input_argument.number_smaller_than_sum_deviation = 0;
						multi.input_argument.number_not_smaller_than_sum_deviation_pow = 0;
						multi.input_argument.number_smaller_than_sum_deviation_pow = 0;

						multi.input_argument.print_each_result = print_each_result;

						APLA::initial_split_coefficients(multi.input_argument, multi.output_argument);
						multi.input_argument.option_split_method = split_method_id;

						multi.all_approximation_build_rtree(multi_y_projection_argument_vector, multi_all_linked_list, multi_cluster_linked_list);

						APLA::count_upper_bound_whole(multi.output_argument, evaluation_bound_whole);
						evaluation_bound_whole.difference_id_max_deviation_vs_height_diff += multi.output_argument.evaluation_bound.difference_id_max_deviation_vs_height_diff;

						const size_t id_method = representation_option - 1;
						const size_t id_each_method_each_file = data_source.file_operation_number * id_method + file_id;
						for (K = initial_K; K <= final_K; K *= 2) {

							if (file_id == initial_file_id && N == initial_N && representation_option == 1) { K_coefficient_vector.emplace_back(K); }
							multi.input_argument.K = K;

							multi.all_knn(data_source, query_time_series_vector, query_time_series, multi_y_projection_argument_vector, squential_scan_result_set);

							cout << "===================      " << name_representation_vector_string_vector[id_method] << " Each K Result        =========================\n";
							cout << "file id : " << file_id_choose_vector[file_id] + 1 << ", K: " << K << ", N:" << multi.input_argument.point_dimension << ", representation id: " << representation_option_vector[representation_option] << endl;
							print_input_output(multi.input_argument, multi.output_argument);
							cout << "======================================================================================\n";

							assert(multi.input_argument.pruning_power > 0 && multi.input_argument.pruning_power != INF && multi.input_argument.prune_power_combine > 0 && multi.input_argument.prune_power_combine != INF && multi.input_argument.whole_run_time != INF);

							APLA::get_split_coefficients(multi.input_argument, multi.output_argument, multi.input_argument.option_split_method, local_total_split_id_sum_deviation, local_total_split_id_shift, local_total_split_id_time, global_total_approximation_sum_deviation, global_total_approximation_time, global_total_knn_prune_power);

							if (K == initial_K) {
								total_sum_deviation_vector.emplace_back(multi.output_argument.sum_deviation);//191206
								total_approximation_time_vector.emplace_back(multi.input_argument.representation_time);//200108
								total_build_tree_time_vector.emplace_back(multi.input_argument.build_rtree_time);//200108
								total_max_deviation_vector.emplace_back(multi.output_argument.max_deviation);//191206
								total_max_deviation_av_vector.emplace_back(multi.output_argument.max_deviation_av);//191206
								total_max_width_deviation_vector.emplace_back(multi.output_argument.max_deviation_multiple_width);//191206
							}
							else {
								total_sum_deviation_vector.emplace_back(0.0);//191206
								total_approximation_time_vector.emplace_back(0.0);//200108
								total_build_tree_time_vector.emplace_back(0.0);//200108
								total_max_deviation_vector.emplace_back(0.0);
								total_max_deviation_av_vector.emplace_back(0.0);
								total_max_width_deviation_vector.emplace_back(0.0);
							}

							total_prune_power_combine_vector.emplace_back(multi.input_argument.prune_power_combine);
							total_prune_power_vector.emplace_back(multi.input_argument.pruning_power);
				
							total_accuracy_vector.emplace_back(multi.input_argument.result_accuracy);
							total_run_time_vector.emplace_back(multi.input_argument.whole_run_time);
							total_run_time_has_IO_vector.emplace_back(multi.input_argument.whole_run_time_has_IO);//210606
							
							total_knn_time_vector.emplace_back(multi.input_argument.knn_total_time);//200108
							total_knn_time_has_IO_vector.emplace_back(multi.input_argument.knn_total_time_has_IO);//210606
							method_sum_prune_power_combine_vector[id_method] += multi.input_argument.prune_power_combine; //each approximation method pruning power combine 201221
							method_sum_prune_power_vector[id_method] += multi.input_argument.pruning_power;
							//method_sum_sum_deviation_vector[representation_option - 1] += multi.output_argument.sum_deviation;//191206
							method_sum_accuracy_vector[id_method] += multi.input_argument.result_accuracy;//191204 accuracy
							method_run_time_vector[id_method] += multi.input_argument.whole_run_time;
							method_run_time_has_IO_vector[id_method] += multi.input_argument.whole_run_time_has_IO;//210606
							//method_approximation_time_vector[representation_option - 1] += multi.input_argument.representation_time;
							method_knn_time_vector[id_method] += multi.input_argument.knn_total_time;
							method_knn_time_has_IO_vector[id_method] += multi.input_argument.knn_total_time_has_IO;//210606

							
							method_file_prune_power_combine_vector[id_each_method_each_file] += multi.input_argument.prune_power_combine;//each file pruning power combine 201221;
							method_file_prune_power_vector[id_each_method_each_file] += multi.input_argument.pruning_power;
							//method_file_sum_deviation_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.output_argument.sum_deviation;//191206
							method_file_accuracy_vector[id_each_method_each_file] += multi.input_argument.result_accuracy;//191204 accuracy
							method_file_run_time_vector[id_each_method_each_file] += multi.input_argument.whole_run_time;
							method_file_run_time_has_IO_vector[id_each_method_each_file] += multi.input_argument.whole_run_time_has_IO;//210606
							//method_file_approximation_time_vector[data_source.file_operation_number * (representation_option - 1) + file_id] += multi.input_argument.representation_time;
							method_file_knn_time_vector[id_each_method_each_file] += multi.input_argument.knn_total_time;
							method_file_knn_time_has_IO_vector[id_each_method_each_file] += multi.input_argument.knn_total_time_has_IO;//210606
							
						}
					
						method_sum_sum_deviation_vector[id_method] += multi.output_argument.sum_deviation;//191206
						method_approximation_time_vector[id_method] += multi.input_argument.representation_time;
						method_build_tree_time_vector[id_method] += multi.input_argument.build_rtree_time;
						method_sum_max_deviation_vector[id_method] += multi.output_argument.max_deviation;
						method_sum_max_width_deviation_vector[id_method] += multi.output_argument.max_deviation_multiple_width;

						method_file_sum_deviation_vector[id_each_method_each_file] += multi.output_argument.sum_deviation;//191206
						method_file_approximation_time_vector[id_each_method_each_file] += multi.input_argument.representation_time;
						method_file_build_tree_time_vector[id_each_method_each_file] += multi.input_argument.build_rtree_time;
						method_file_max_deviation_vector[id_each_method_each_file] += multi.output_argument.max_deviation;//191206
						method_file_max_width_deviation_vector[id_each_method_each_file] += multi.output_argument.max_deviation_multiple_width;//191206
						
						number_point_max_deviation_true += multi.input_argument.number_point_max_deviation_true;
						number_point_max_deviation_false += multi.input_argument.number_point_max_deviation_false;
						number_not_smaller_than_sum_deviation += multi.input_argument.number_not_smaller_than_sum_deviation;
						number_smaller_than_sum_deviation += multi.input_argument.number_smaller_than_sum_deviation;
						number_not_smaller_than_sum_deviation_pow += multi.input_argument.number_not_smaller_than_sum_deviation_pow;
						number_smaller_than_sum_deviation_pow += multi.input_argument.number_smaller_than_sum_deviation_pow;
					}
				}
				change_file = false;

			}
			else if (evaluation_type == 1) {
				assert(0);
			}
			else {
				assert(0);
			}
		}


		data_source.read_file_address_vector.clear();
		data_source.read_file_address_vector.shrink_to_fit();
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

		delete[] query_time_series;
		query_time_series = nullptr;

	}

	cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      Final Result     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";

	TOOL::transfer_us_s(total_run_time_vector);
	TOOL::transfer_us_s(total_run_time_has_IO_vector);//210606
	TOOL::transfer_us_s(total_approximation_time_vector);
	TOOL::transfer_us_s(total_build_tree_time_vector);
	TOOL::transfer_us_s(total_knn_time_vector);
	TOOL::transfer_us_s(total_knn_time_has_IO_vector);//210606
	TOOL::transfer_us_s(method_run_time_vector);
	TOOL::transfer_us_s(method_run_time_has_IO_vector);//210606
	TOOL::transfer_us_s(method_approximation_time_vector);
	TOOL::transfer_us_s(method_build_tree_time_vector);
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
	TOOL::transfer_us_s(method_file_knn_time_vector);
	TOOL::transfer_us_s(method_file_knn_time_has_IO_vector);//210606
	
	cout << "Method Sum Prune Power: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_sum_prune_power_vector) {
		cout << au << ", ";
	}
	cout << endl;
	cout << "Method Max Deviation: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_sum_max_deviation_vector) {
		cout << au << ", ";
	}
	cout << endl;
	cout << "Method Sum Accuracy: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_sum_accuracy_vector) {
		cout << au << ", ";
	}
	cout << endl;
	cout << "Apporximation Sum Time(s): \n";
	//cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07  PAALM\n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_approximation_time_vector) {
		cout << au << ", ";
	}
	cout << endl;
	cout << "KNN Sum Time(s): \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_knn_time_vector) {
		cout << au << ", ";
	}
	cout << endl;
	system("pause");
}


TEMPLATE
template<typename T, typename Y>
inline void Evaluation::print_input_output(const T& const input_argument, const Y& const output_argument) {
	
	cout << "*prune power: " << input_argument.pruning_power << endl;
	cout << "*max deviation: " << output_argument.max_deviation << endl;
	cout << "*prune accuracy: " << input_argument.result_accuracy << endl;
	cout << "*representation_time (us): " << input_argument.representation_time << endl;
	cout << "*knn_total_time (us): " << input_argument.knn_total_time << endl;
	
}

#endif

