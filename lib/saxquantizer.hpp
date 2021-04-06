

#ifndef SAXQUANTIZER_HPP
#define SAXQUANTIZER_HPP

#include "../CAPCA.h"
#include <deque>
#include <vector>
#include <cassert>
#include <boost/math/distributions/normal.hpp>

using std::deque;
using std::vector;
using namespace std;
namespace SaxQuantizer {

	inline void fill_cutpoints(const size_t& alphabet_size, vector<double> *cutpoints) {
		assert(alphabet_size > 0);
		static boost::math::normal dist(0.0, 1.0);
		//std::cout << "alphabet: " << alphabet_size << std::endl;
		cutpoints->reserve(alphabet_size);
		cutpoints->push_back(-DBL_MAX);
		//cout << "cdf: ";
		for (size_t i = 1; i < alphabet_size; ++i) {
			double cdf = ((double)i) / alphabet_size;
			//cout << quantile(dist, cdf) << " ";
			cutpoints->push_back(quantile(dist, cdf));
			//cout << cutpoints->begin();
		}
		//cout << endl;
	}

	/**
	 * Symbolic Aggregate Approximation with fractional sliding window, numerosity reduction, and scaling
	 */
	class SAX : virtual public APCA_QUAL {
	private:
		size_t m_window_size;
		size_t m_string_size;// the number of segments, N .
		size_t m_alphabet_size;

		double m_baseline_mean;
		double m_baseline_stdev;
		vector<double> m_cutpoints;

		bool m_trained;

		size_t segment_number;//200923 

		/**
		 * SAX with fractional sliding window and automatic scaling.
		 * @param <it>: start iterator
		 * @param <end>: end iterator
		 * @param <syms> (out): quantized range
		 */
		template<class Iter>
		void saxify(Iter it, const Iter end, vector<int> *syms) {
			// perform PAA using a fractional sliding window
		   // double paa[m_string_size];
			double* paa = new double[m_string_size];
			double paa_window = ((double)m_window_size) / m_string_size;

			double p = 0; // p for progress
			double w = 1; // w for weight
			size_t available = 0;

			for (size_t i = 0; i < m_string_size && it != end; ++i, ++available) {
				// normalize around baseline

				double normalized = (*it - m_baseline_mean) / m_baseline_stdev;

				paa[i] = 0;
				double j = 0;

				while (j < paa_window && it != end) {

					paa[i] += w * normalized; // sum of (partial) elements inside the window
					j += w;
					p += w;

					// window full
					if (paa_window == p) {

						if (fabs(w - 1.0) <= 0.01) {   // if last element fully consumed,
							++it;                        // then just move next.
						}
						else {                       // o.w.,
							w = 1.0 - w;                 // set remaining portion.
						}
						p = 0;                         // reset progress

					  // window not full, but next must be split
					}
					else if (paa_window - p < 1.0) {

						w = paa_window - p;            // set needed portion
						++it;                          // move to next

					  // window not full, next can be fully consumed
					}
					else {
						++it;                          // move to next
					}
				}

				paa[i] /= j; // averaging
			}

			// map to symbols. 0-based.
			//cout << "available" << available << endl;
			for (size_t i = 0; i < available; ++i) {
				int cnt = -1;
				for (const auto & cp : m_cutpoints) {
					if (paa[i] >= cp) ++cnt;
				}
				syms->push_back(cnt);
			}
			delete[] paa;
		}

	public:
		/**
		* Constructs a SAX quantizer of a given window size, string size and alphabet size.
		* @param <window_size>: sliding window size
		* @param <string_size>: output string size for each sliding window (can be greater than window_size)
		* @param <alphabet_size>: number of codewords
		*/
		SAX(size_t window_size, size_t string_size, size_t alphabet_size)
			: m_window_size(window_size), m_string_size(string_size), m_alphabet_size(alphabet_size),
			m_baseline_mean(0), m_baseline_stdev(1), m_trained(false) {

			assert(window_size > 0);
			assert(string_size > 0);
			assert(alphabet_size > 0);

			fill_cutpoints(alphabet_size, &m_cutpoints);
		}


		/**
		* Constructs a SAX quantizer of timeseries and alphabet size.
		* @param <window_size>: sliding window size
		* @param <string_size>: output string size for each sliding window (can be greater than window_size)
		* @param <alphabet_size>: number of codewords
		
		*/
		SAX(const size_t& alphabet_size): m_window_size(1), m_string_size(1), m_alphabet_size(alphabet_size),m_baseline_mean(0), m_baseline_stdev(1), m_trained(false) {

			assert(alphabet_size > 0);

			fill_cutpoints(alphabet_size, &m_cutpoints);
		}


		virtual ~SAX() {
			m_cutpoints.clear();
		}



		/**
		 * Trains the quantizer from a given sample. This sets the baseline mean and stdevs, which are used in
		 * normalizing the input.
		 *
		 * @param <samples>: list of training values
		 */
		template<typename Container>
		void train(const Container & samples) {
			double mean = 0;
			double stdev = DBL_MIN;

			assert(!samples.empty());

			if (samples.size() < 2) {
				mean = samples[0];
				stdev = DBL_MIN;

			}
			else {
				size_t n = 0;
				double M2 = 0;
				for (const auto & val : samples) {
					++n;
					double delta = val - mean;
					mean += delta / n;
					M2 += delta * (val - mean);
				}
				stdev = sqrt(M2 / (n - 1));
			}

			if (stdev == 0) stdev = DBL_MIN;

			m_baseline_mean = mean;
			m_baseline_stdev = stdev;

			m_trained = true;
		}

		template<typename Container>
		size_t quantize(const Container & seq, vector<int> *qseq, bool reduce = true) {
			if (!m_trained) train(seq);

			vector<int> buf1, buf2;
			auto *syms_buf = &buf1;
			auto *old_syms_buf = &buf2;

			size_t consumed = 0;
			for (consumed = 0; consumed < seq.size(); ++consumed) {

				if (reduce) { // run-length numerosity reduction
					syms_buf->clear();
					saxify(seq.begin() + consumed, seq.end(), syms_buf);

					// skip window if same as previous
					if (*syms_buf != *old_syms_buf) {
						qseq->insert(qseq->end(), syms_buf->begin(), syms_buf->end());
						std::swap(syms_buf, old_syms_buf);
					}

				}
				else { // no reduction
					saxify(seq.begin() + consumed, seq.end(), qseq);
				}

				// ignore excess elements, if sequence size isn't a multiple of window size
				//if (seq.size() - consumed <= m_window_size) break;
			}

			return consumed;
		}

		/**
		* Constructs a SAX quantizer of timeseries and alphabet size.
		* @param <window_size>: sliding window size
		* @param <string_size>: output string size for each sliding window (can be greater than window_size)
		* @param <alphabet_size>: number of codewords
		* @date : 2018/3/31 15:48
		* @author :  
		*/
		template<typename Container>
		size_t getSAX(const Container& seq, vector<char>& qseq) {
			if (!m_trained) train(seq);
			double normalized = NULL;
			for (auto& it : seq) {
				//cout << it << endl;
				normalized = (it - m_baseline_mean) / m_baseline_stdev;
				int cnt = -1;
				//cout << "*******************************************" << endl;
				for (const auto & cp : m_cutpoints) {
					//cout << normalized << ", " << cp << ", " << (normalized > cp ? true:false) << endl;
					if (normalized >= cp) ++cnt;
				}
				qseq.push_back(static_cast<char>(cnt + 65));
			}

			return 0;
		}

		/**
		 * Returns the order of the quantizer (here, the window size)
		 */
		inline size_t order() const {
			return m_window_size;
		}


		inline double ratio() const {
			return ((double)m_window_size) / m_string_size;
		}

		template<typename T, typename Y, typename U>
		void get_SAX(const vector<T>& const original_time_series_vector, const Y& const segment_number, const U& const SAX_container) {

			
			vector<double> paa_vector;
			vector<int> SAX_quantizer_vector;

			APCA_QUAL::get_PAA(original_time_series_vector, segment_number, SAX_container);

			for (int i = 0; i < segment_number; i++) {
				paa_vector.emplace_back(SAX_container.v[i]);
			}
		
			quantize(paa_vector, &SAX_quantizer_vector, false);
		
			copy_n(SAX_quantizer_vector.begin(), SAX_quantizer_vector.size(), SAX_container.v);
		
		}

		/**
		* @Name: distance_cell
		* @Qualifier: cell(r,c)
		* @Date: 200923
		* @author: 
		*/
		template<typename T>
		inline double get_distance_cell(const T& const vaule_1, const T& const vaule_2) {
			if (fabs(vaule_1 - vaule_2) <= 1) {
				return 0;
			}
			else {
				return m_cutpoints[max(vaule_1 + 1, vaule_2 + 1) - 1] - m_cutpoints[min(vaule_1 + 1, vaule_2 + 1)];
			}
		}

		template<typename T>
		double distance_LB_SAX(const T& const SAX_container1, const T& const SAX_container2) {

			const auto& const QProjection = SAX_container1;
			const auto& const italicC = SAX_container2;

			int i = 0, j = 0;
			double distance = 0;

			double sum = (italicC.r[0] + 1) * pow(get_distance_cell(QProjection.v[0], italicC.v[0]), 2);

			for (i = 1; i < italicC.segmentNum; i++) {
				sum += (italicC.r[i] - italicC.r[i - 1]) * pow( get_distance_cell( QProjection.v[i], italicC.v[i] ), 2 );
			}
			distance = sqrt(sum);
			return distance;
		}
	};
}

#endif
