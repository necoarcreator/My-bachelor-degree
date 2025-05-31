#include "FT.h"
#define _USE_MATH_DEFINES
const complex<double> i(0., 1.);
#include <vector>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <exception>
#include <omp.h>
#include <fftw3.h>

using namespace std;
FT::FT(vector<double> _omega, int num_threads = 1) : omega(_omega) {
	
	N = _omega.size();
	Y = vector<complex<double>>(N);
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	plan_forward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_backward = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_init_threads();
	fftw_plan_with_nthreads(num_threads);
		//omp_get_max_threads());
};

vector<complex<double>> FT::DFT(vector<complex<double>>& W, vector<double> x, long L, bool isInverse) {

	complex<double> intSum;
	double h = x[1] - x[0];

	if (!isInverse) {
		#pragma omp parallel for
		for (ptrdiff_t k = 0; k < N; k++) {

			intSum = complex<double>(0., 0.);
			for (ptrdiff_t n = 0; n < N; n++) {

				//iterate over x_n
				complex<double> w = exp(-i * omega[k] * x[n]);

				intSum += W[n] * w;
			}

			Y[k] = h / L * intSum;
		}
	}
	else {
		#pragma omp parallel for
		for (ptrdiff_t n = 0; n < N; n++) {

			intSum = complex<double>(0., 0.);
			for (ptrdiff_t k = 0; k < N; k++) {

				//iterate over w_k
				complex<double> w = exp(i * omega[k] * x[n]);

				intSum += W[k] * w;
			}

			Y[n] = intSum;
		}
	}
	W = Y;
	return W;
}
vector<complex<double>> FT::FFT(vector<complex<double>>& W, vector<double> x, long L, bool isInverse) {
	
	double h = x[1] - x[0];
	for (int n = 0; n < N; n++) {
		int shifted_idx = (n + N / 2) % N;  // Сдвигаем на N/2 из-за порядка частот в fftw
		Y[n] = W[shifted_idx];
	}
	if (!isInverse) {
		// Прямое преобразование
		for (int n = 0; n < N; n++) {
			in[n][0] = Y[n].real() ;
			in[n][1] = Y[n].imag() ;
		}
		fftw_execute(plan_forward); //O(N log N)

		
		for (int k = 0; k < N; k++) {
			// Применяем нормировку h/L
			Y[k] = complex<double>(out[k][0], out[k][1]) * (h / L);
		}
	}
	else {
		// Обратное преобразование
		for (int k = 0; k < N; k++) {
			in[k][0] = Y[k].real();
			in[k][1] = Y[k].imag();
		}

		fftw_execute(plan_backward);

		for (int n = 0; n < N; n++) {
			Y[n] = complex<double>(out[n][0], out[n][1]);
		}
	}
	//обратный сдвиг к нашему исходному порядку
	for (int k = 0; k < N; k++) {
		int shifted_idx = (k - N / 2) % N;
		W[k] = Y[shifted_idx];
	}
 
	return W;
}

FT::~FT() {
	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	fftw_free(in);
	fftw_free(out);
}
//void FT::fftRecursive(vector<complex<double>>& W, bool isInverse) {}
//	size_t N = W.size();
//	if (N <= 1) return;
//
//	// Разделяем на чётные/нечётные
//	size_t half = N / 2;
//	vector<complex<double>> even(half), odd(half);
//
//	for (size_t i = 0; i < half; ++i) {
//		even[i] = W[2 * i];
//		odd[i] = W[2 * i + 1];
//	}
//
//	// Рекурсивные вызовы
//	fftRecursive(even, isInverse);
//	fftRecursive(odd, isInverse);
//
//	// сопряжение или нет, в зав-ти от прямого/обр преобр
//	double sign = (isInverse ? +1.0 : -1.0);
//
//	for (size_t k = 0; k < half; ++k)
//	{
//		double angle = 2.0 * M_PI * k / double(N);
//
//		// twiddle = exp( i * sign * angle ) * odd
//		complex<double> twiddle = exp(i * (sign * angle)) * odd[k];
//
//		W[k] = even[k] + twiddle;
//		W[k + half] = even[k] - twiddle;
//		if (isInverse) {
//			W[k] /= 2.;
//			W[k + half] /= 2.;
//		}
//	}
//
//}
//vector<complex<double>> FT::FFT(vector<complex<double>>& W, vector<double> x, long L, bool isInverse) {
//	const size_t M = W.size();
//	double h = x[1] - x[0];
//	
//
//	complex<double> bias = exp(i * (L * M_PI / 2));
//
//	
//	if ((M & (M - 1)) != 0) {
//		throw invalid_argument("FFT size must be a power of 2.");
//	}
//
//	// вызываем рекурсивный метод
//	fftRecursive(W, isInverse);
//	
//	for (size_t k = 0; k < N / 2; k++) {
//		// swap W[k] <-> W[k+N/2]
//		// чтобы "центрировать" ноль частоты
//		size_t k2 = k + N / 2;
//		auto tmp = W[k];
//		W[k] = W[k2];
//		W[k2] = tmp;
//		if (!isInverse) {
//			W[k] *= (h / L);
//			W[k2] *= (h / L);
//		}
//	}
//
//	// Фазовый сдвиг во временной области, чтобы компенсировать -L/2
//	// Обычно: W[n] *= exp( + i * pi * n ) и т. д., но у нас симметрия:
//	for (size_t n = 0; n < N; n++) {
//
//		if (n % 2) {
//			W[n] = -W[n];
//		}
//	}
//
//	Y = W;
//	return Y;
//}


//void FT::printFileFT(string path) {}
//	ofstream file(path);
//	if (file.is_open())
//	{
//		file << "Order		" << "Re		" << "Im		" << endl;
//		for (size_t k = 0; k < N; k++)
//		{
//			file << k << "\t" << setw(12) << Y[k].real() << "\t" << setw(12) <<  Y[k].imag() << endl;
//		}
//	}
//	else
//	{
//		cerr << "Couldn't open file to write" << endl;
//		exit(1);
//	}
//	file.close();
//	return;
//}