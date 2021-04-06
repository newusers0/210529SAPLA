#pragma once
#ifndef _EVALUATION_CPP_
#define _EVALUATION_CPP_

#include "EVALUATION.h"


TEMPLATE
void Evaluation::evaluate_multi_KNN_speed() {
	
	string read_file_name = "CinC_ECG_torso_TEST";
	string* chebyshev_write_pointer = nullptr;
	string dataset_attribute[] = { "", "Single dataset", "Mixed dataset", "" };
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
    size_t representation_option = INF;
	const vector<int> representation_option_vector{ int(INF), 2, 3, 4, 5, 6, 7, 8, 9 };// has SAPLA
	//const vector<int> representation_option_vector{ int(INF), 2, 8};//Only SAPLA
	//const vector<int> representation_option_vector{ int(INF), 5};// CHEBY
	//const vector<int> representation_option_vector{ int(INF), 2};// PLA
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
	int m_d = INF;
	int file_total_number = INF;
	int n0 = INF;
	int point_number0 = INF;
	const int& const evaluation_type = 0;
	const size_t& const data_type = 1; 
	string mixed_data_file_name = "./191202DataSet/04mixed_data_set";
	string* file_address_pointer = nullptr;
	const string& const str_suffix = "210328";
	const bool clear_file = true;
	const int option_homogenous_data_type = 0;
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
		m_d = 3;
		switch (option_homogenous_data_type) {
		case 0:
			file_id_choose_vector.resize(20);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 1);
			file_total_number = file_id_choose_vector.size();
			n0 = 1024;
			final_N = 384;
			point_number0 = 100;
			final_K = 64;
			m_d = 3;
			query_time_series_id = point_number0 - 5;
			break;
		case 1:
			file_id_choose_vector.resize(21);
			std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
			file_total_number = file_id_choose_vector.size();
			n0 = 1024;
			final_N = 384;
			point_number0 = 20;
			final_K = -1;
			m_d = 2;
			query_time_series_id = 0;
			break;
		default:
			assert(0);
		}
		data_final_dimension = 1;
		assert(data_final_dimension == 1);
		file_address_pointer = nullptr;
		break;
	case 2:
		n0 = 1024;
		initial_N = 12;
		final_N = 384;
		point_number0 = 20;
		final_K = 128;
		m_d = 5;
		file_total_number = 20;
		data_final_dimension = 1;
		assert(data_final_dimension == 1);
		file_address_pointer = &mixed_data_file_name;
		file_id_choose_vector.resize(1,0);
		query_time_series_id = 359;
		break;
	case 3:
		n0 = 295;
		initial_N = 12;
		final_N = 192;
		point_number0 = 300;
		final_K = 128;
		m_d = 10;
		file_id_choose_vector.resize(2);
		std::iota(file_id_choose_vector.begin(), file_id_choose_vector.end(), 0);
		data_final_dimension = 3;
		file_total_number = file_id_choose_vector.size() * data_final_dimension;
		file_address_pointer = nullptr;
		query_time_series_id = 89;
		break;
	default:
		assert(0);
		break;
	}
	initial_file_id = 0;
	const int final_file_id = file_id_choose_vector.size();
	initial_K = 2;
	if (option_homogenous_data_type == 6) initial_K = 0;
	const int method_total_number = representation_option_vector.size() - 1;
	vector<size_t> N_coefficient_vector;
	vector<int> K_coefficient_vector;
	vector<long double> total_prune_power_combine_vector;
	vector<long double> total_prune_power_vector;
	vector<long double> total_sum_deviation_vector;
	vector<long double> total_max_deviation_vector;
	vector<long double> total_max_width_deviation_vector;
	vector<long double> total_accuracy_vector;
	vector<long double> total_run_time_vector;
	vector<long double> total_approximation_time_vector;
	vector<long double> total_knn_time_vector;
	vector<long double> method_sum_prune_power_combine_vector(method_total_number, 0);
	vector<long double> method_sum_prune_power_vector(method_total_number, 0);
	vector<long double> method_sum_sum_deviation_vector(method_total_number, 0);
	vector<long double> method_sum_max_deviation_vector(method_total_number, 0);
	vector<long double> method_sum_max_width_deviation_vector(method_total_number, 0);
	vector<long double> method_sum_accuracy_vector(method_total_number, 0);
	vector<long double> method_run_time_vector(method_total_number, 0);
	vector<long double> method_approximation_time_vector(method_total_number, 0);
	vector<long double> method_knn_time_vector(method_total_number, 0);
	vector<long double> one_file_prune_power_combine_vector(method_total_number, 0);
	vector<long double> one_file_prune_power_vector(method_total_number, 0);
	vector<long double> one_file_sum_deviation_vector(method_total_number, 0);
	vector<long double> one_file_max_deviation_vector(method_total_number, 0);
	vector<long double> one_file_max_width_deviation_vector(method_total_number, 0);
	vector<long double> one_file_accuracy_vector(method_total_number, 0);
	vector<long double> one_file_run_time_vector(method_total_number, 0);
	vector<long double> one_file_approximation_time_vector(method_total_number, 0);
	vector<long double> one_file_knn_time_vector(method_total_number, 0);
	vector<long double> method_file_prune_power_combine_vector(size_t(method_total_number* file_total_number), 0.0);
	vector<long double> method_file_prune_power_vector(size_t(method_total_number * file_total_number), 0.0);
	vector<long double> method_file_sum_deviation_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_max_deviation_vector(method_total_number* file_total_number, 0.0);
	vector<long double> method_file_max_width_deviation_vector(method_total_number* file_total_number, 0.0);
	vector<long double> method_file_accuracy_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_run_time_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_approximation_time_vector(method_total_number * file_total_number, 0.0);
	vector<long double> method_file_knn_time_vector(method_total_number * file_total_number, 0.0);
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
	vector<size_t> data_list;
	for (int list_id = 1; list_id < file_total_number+1; list_id++) {data_list.emplace_back(list_id);}
	typename TOOL::DATA_SOURCE data_source(data_type, data_final_dimension, data_list, file_total_number, point_number0, n0, *file_address_pointer);// 0 manipulate data set; 1 real data set, homogenous(single) data; 2 real data set, heterogeneous(mixed) data. 3 multi single data. 4 multi mixed data
	data_source.option_has_burst_data = option_homogenous_data_type;
	TOOL::initial_data_source(data_source);
	data_source.bigger_account = 0;
	assert(data_source.point_number != INF && data_source.single_time_series_length != INF && data_source.multi_single_time_series_length != INF);
	for (file_id = initial_file_id; file_id < final_file_id; file_id++) {
		bool change_file = true;
		bool print_each_result = false;
		if (print_each_result == true) {
			vector<DataType> original_time_series_vector_print(data_source.multi_single_time_series_length, INF);
			TOOL::get_read_multi_file_address(data_source, file_id_choose_vector[file_id]);
			for (int i = 0; i < data_source.point_number; i++) {
				TOOL::getMultiFoldToSingleByID(data_source.read_file_address_vector, data_source.time_series_dimension, data_source.single_time_series_length, i, original_time_series_vector_print);
				TOOL::print_vector(original_time_series_vector_print);//190501
				original_time_series_vector_print.clear();
				original_time_series_vector_print.shrink_to_fit();
			}
			data_source.read_file_address_vector.clear();
			data_source.read_file_address_vector.shrink_to_fit();
		}
		vector<DataType> query_time_series_vector(data_source.multi_single_time_series_length, INF);
		DataType* query_time_series = new DataType[data_source.multi_single_time_series_length];
		TOOL::get_read_multi_file_address(data_source, file_id_choose_vector[file_id]);
		//already normalized in below function
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
						MULTI multi(data_source, n0, N, INF, file_id_choose_vector[file_id], change_file, data_source.time_series_dimension, data_source.point_number, query_time_series_id, m_d, initial_K, representation_option_vector[representation_option], file_address_pointer, chebyshev_write_pointer);
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
						evaluation_bound_whole.difference_id_max_deviation_vs_height_diff += multi.output_argument.evaluation_bound.difference_id_max_deviation_vs_height_diff;
						const size_t id_method = representation_option - 1;
						const size_t id_each_method_each_file = data_source.file_operation_number * id_method + file_id;
						for (K = initial_K; K <= final_K; K *= 2) {
							if (file_id == initial_file_id && N == initial_N && representation_option == 1) { K_coefficient_vector.emplace_back(K); }
							multi.input_argument.K = K;
							multi.all_knn(data_source, query_time_series_vector, query_time_series, multi_y_projection_argument_vector, squential_scan_result_set);
							cout << name_representation_vector_string_all_vector[representation_option_vector[representation_option]-1] << " file id : " << file_id_choose_vector[file_id] + 1 << ", K: " << K << ", N:" << multi.input_argument.point_dimension << endl;
							cout << "prune power: " << multi.input_argument.pruning_power << endl;
							cout << "sum deviation: " << multi.output_argument.sum_deviation << endl;
							cout << "max deviation: " << multi.output_argument.max_deviation << endl;
							cout << "prune accuracy: " << multi.input_argument.result_accuracy << endl;
							cout << "representation_time(us): " << multi.input_argument.representation_time << endl;
							cout << "build_rtree_time(us): " << multi.input_argument.build_rtree_time << endl;
							cout << "knn_total_time(us): " << multi.input_argument.knn_total_time << endl;
							cout << "accumulation time(us): " << multi.input_argument.representation_time + multi.input_argument.build_rtree_time + multi.input_argument.knn_total_time << endl;
							cout << "whole running time(us): " << multi.input_argument.whole_run_time << endl;
							APLA::get_split_coefficients(multi.input_argument, multi.output_argument, multi.input_argument.option_split_method, local_total_split_id_sum_deviation, local_total_split_id_shift, local_total_split_id_time, global_total_approximation_sum_deviation, global_total_approximation_time, global_total_knn_prune_power);
							if (K == initial_K) {
								total_sum_deviation_vector.emplace_back(multi.output_argument.sum_deviation);//191206
								total_approximation_time_vector.emplace_back(multi.input_argument.representation_time);//200108
								total_max_deviation_vector.emplace_back(multi.output_argument.max_deviation);//191206
								total_max_width_deviation_vector.emplace_back(multi.output_argument.max_deviation_multiple_width);//191206
							}
							else {
								total_sum_deviation_vector.emplace_back(0.0);//191206
								total_approximation_time_vector.emplace_back(0.0);//200108
								total_max_deviation_vector.emplace_back(0.0);
								total_max_width_deviation_vector.emplace_back(0.0);
							}
							total_prune_power_combine_vector.emplace_back(multi.input_argument.prune_power_combine);
							total_prune_power_vector.emplace_back(multi.input_argument.pruning_power);//sum pruning power
							total_accuracy_vector.emplace_back(multi.input_argument.result_accuracy);
							total_run_time_vector.emplace_back(multi.input_argument.whole_run_time);
							total_knn_time_vector.emplace_back(multi.input_argument.knn_total_time);//200108
							method_sum_prune_power_combine_vector[id_method] += multi.input_argument.prune_power_combine;
							method_sum_prune_power_vector[id_method] += multi.input_argument.pruning_power;
							method_sum_accuracy_vector[id_method] += multi.input_argument.result_accuracy;
							method_run_time_vector[id_method] += multi.input_argument.whole_run_time;
							method_knn_time_vector[id_method] += multi.input_argument.knn_total_time;
							if (!clear_file) {
								one_file_prune_power_combine_vector[id_method] += multi.input_argument.prune_power_combine;
								one_file_prune_power_vector[id_method] += multi.input_argument.pruning_power;
								one_file_accuracy_vector[id_method] += multi.input_argument.result_accuracy;
								one_file_run_time_vector[id_method] += multi.input_argument.whole_run_time;
								one_file_knn_time_vector[id_method] += multi.input_argument.knn_total_time;
							}
							method_file_prune_power_combine_vector[id_each_method_each_file] += multi.input_argument.prune_power_combine;
							method_file_prune_power_vector[id_each_method_each_file] += multi.input_argument.pruning_power;
							method_file_accuracy_vector[id_each_method_each_file] += multi.input_argument.result_accuracy;//191204 accuracy
							method_file_run_time_vector[id_each_method_each_file] += multi.input_argument.whole_run_time;
							method_file_knn_time_vector[id_each_method_each_file] += multi.input_argument.knn_total_time;
						}
						method_sum_sum_deviation_vector[id_method] += multi.output_argument.sum_deviation;//191206
						method_approximation_time_vector[id_method] += multi.input_argument.representation_time;
						method_sum_max_deviation_vector[id_method] += multi.output_argument.max_deviation; 
						method_sum_max_width_deviation_vector[id_method] += multi.output_argument.max_deviation_multiple_width;
						if (!clear_file) {
							one_file_sum_deviation_vector[id_method] += multi.output_argument.sum_deviation;//191206
							one_file_approximation_time_vector[id_method] += multi.input_argument.representation_time;
						}
						method_file_sum_deviation_vector[id_each_method_each_file] += multi.output_argument.sum_deviation;//191206
						method_file_approximation_time_vector[id_each_method_each_file] += multi.input_argument.representation_time;
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
		if (!clear_file) {
			one_file_prune_power_combine_vector.resize(method_total_number, 0.0);//201221
			one_file_prune_power_vector.resize(method_total_number, 0.0);
			one_file_sum_deviation_vector.resize(method_total_number, 0.0);
			one_file_accuracy_vector.resize(method_total_number, 0.0);
			one_file_run_time_vector.resize(method_total_number, 0.0);
			one_file_approximation_time_vector.resize(method_total_number, 0.0);
			one_file_knn_time_vector.resize(method_total_number, 0.0);
		}
	}
	TOOL::transfer_us_s(total_run_time_vector);
	TOOL::transfer_us_s(total_approximation_time_vector);
	TOOL::transfer_us_s(total_knn_time_vector);
	TOOL::transfer_us_s(method_run_time_vector);
	TOOL::transfer_us_s(method_approximation_time_vector);
	TOOL::transfer_us_s(method_knn_time_vector);
	TOOL::transfer_us_s(one_file_run_time_vector);
	TOOL::transfer_us_s(one_file_approximation_time_vector);
	TOOL::transfer_us_s(one_file_knn_time_vector);

	TOOL::transfer_us_s(method_file_run_time_vector);
	TOOL::transfer_us_s(method_file_approximation_time_vector);
	TOOL::transfer_us_s(method_file_knn_time_vector);
	cout << "Method Sum Prune Power: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_sum_prune_power_vector) {
		cout << au << ", ";
	}
	cout << endl;
	cout << "Method Sum Deviation: \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_sum_sum_deviation_vector) {
		cout << au << ", ";
	}
	cout << endl;
	cout << "Method Max Deviation: \n";
	//cout << "MSPLA  PLA  APCA   PAA  CHEBY   ICDE07  PAALM\n";
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
	cout << "Sum Time(s): \n";
	TOOL::print_string_vector(name_representation_vector_string_vector);
	for (auto&& au : method_run_time_vector) {
		cout << au << ", ";
	}
	cout << endl;
	system("pause");
}
#endif