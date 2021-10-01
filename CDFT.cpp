#pragma once
#ifndef _CDFT_CPP_
#define _CDFT_CPP_


//TEMPLATE
//template<typename T>
//DFT::CDFT(const T& n, const T& N) {
//	input_argument.time_series_length = n;
//	input_argument.point_dimension = N;
//}
//
//TEMPLATE
//vector<complex<double> > DFT::computeDft(const vector<complex<double> >& input) {
//	vector<complex<double> > output;
//	size_t n = input.size();
//	for (size_t k = 0; k < n; k++) {  // For each output element
//		complex<double> sum(0.0, 0.0);
//		for (size_t t = 0; t < n; t++) {  // For each input element
//			double angle = 2 * boost::math::constants::pi<double>() * t * k / n;
//			sum += input[t] * exp(-angle);
//		}
//		output.push_back(sum);
//	}
//	return output;
//}
//
//TEMPLATE
//void DFT::computeDft(const vector<double>& inreal, const vector<double>& inimag, vector<double>& outreal, vector<double>& outimag) {
//	size_t n = inreal.size();
//	for (size_t k = 0; k < n; k++) {  // For each output element
//		double sumreal = 0;
//		double sumimag = 0;
//		for (size_t t = 0; t < n; t++) {  // For each input element
//			double angle = 2 * boost::math::constants::pi<double>() * t * k / n;
//			sumreal += inreal[t] * cos(angle) + inimag[t] * sin(angle);
//			sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
//		}
//		outreal[k] = sumreal;
//		outimag[k] = sumimag;
//	}
//}
//
//TEMPLATE
//void DFT::computeDft(const typename TOOL::INPUT_ARGUMENT& input_argument, DataType original_time_series[], DataType*& dft_time_series) {
//	printf("Compute DFT\n");
//	size_t n = input_argument.time_series_length;
//	long double sum = NULL;
//	for (size_t k = 0; k < input_argument.point_dimension; k++) {  // For each output element
//		sum = 0;
//		for (size_t t = 0; t < n; t++) {  // For each input element
//			double angle = 2 * boost::math::constants::pi<double>() * t * k / n;
//			//cout << angle << endl;
//			sum += original_time_series[t] * exp(-angle);
//		}
//		dft_time_series[k] = sum;
//		cout << sum << endl;
//	}
//}

#endif
