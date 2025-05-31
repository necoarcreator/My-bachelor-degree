#ifndef _FT_
#define _FT_

#include <vector>
#include <complex>
#include <string>
#include <fftw3.h>

using namespace std;
class FT
{
	size_t N;
	vector<double> omega;
	vector<complex<double>> Y;
	fftw_plan plan_forward, plan_backward;
	fftw_complex *in, *out;

public:

	FT(vector<double> _omega, int num_threads);
	~FT();
	vector<complex<double>> FFT(vector<complex<double>>& W, vector<double> x, long L, bool isInverse);
	//void fftRecursive(vector<complex<double>>& W, bool isInverse = false);

	vector<complex<double>> DFT(vector<complex<double>>& W, vector<double> x, long L, bool isInverse = false);
	//void printFileFT(string path = "file.txt");
};

#endif