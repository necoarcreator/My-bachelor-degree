#include <iostream>
#include "FT.h"
#include <vector>
#include <complex>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <numeric>
#include <random>
#define _USE_MATH_DEFINES
const complex<double> i(0., 1.);
#include <math.h>
using namespace std;
void  printFile(vector<complex<double>> var, string path, long long int L);

double energyIntegral(vector<complex<double>> q, double h) {
	double qModuleSq = 0.;
	for (size_t k = 0; k < q.size(); ++k) {
		qModuleSq += q[k].real() * q[k].real() + q[k].imag() * q[k].imag();
	}
	return qModuleSq * h;
}
double energyFluxIntegral(vector<complex<double>> q, double h) {
	complex<double> res = 0.;
	vector<complex<double>> dq(q.size());

	dq[0] = (q[1] - q[q.size() - 1]) / (2 * h * h);
	for (size_t k = 1; k < q.size() - 1; ++k) {
		dq[k] = (q[k + 1] - q[k - 1]) / (2 * h * h);
	}
	dq[q.size() - 1] = (q[0] - q[q.size() - 2]) / (2 * h * h);


	for (size_t k = 0; k < q.size(); ++k) res += dq[k] * conj(q[k]) - conj(dq[k]) * q[k];

	return real(-i * h * res);
}
void filtrationSmooth(vector<complex<double>>& q, vector<double> omega, vector<double> x, double h, double alpha, double kCutDole, double L, size_t step) {

	size_t N = q.size();

	// вычисляем полуширину окна в индексах

	FT ft(omega);

	ft.DFT(q, x, L, false);

	int maxIt = distance(q.begin(), max_element(q.begin(), q.end(), [](auto a, auto b) { return abs(a) < abs(b); }));

	double maxK = abs(q[maxIt]);
	double kCut = maxK * kCutDole;

	double I1_before = energyIntegral(q, h);
	double I2_before = energyFluxIntegral(q, h);

	// Делаем мягкое отсечение высоких частот

	for (size_t p = 0; p < N; ++p) {
		if ((kCut < abs(q[p])) && (abs(q[p]) <= maxK)) q[p] *= exp(-alpha * pow(abs(q[p]) - maxK, 2));

		else if (abs(q[p]) > maxK) q[p] = 0;
	}
	ft.DFT(q, x, L, true);

	// Считаем, сколько "энергии" (I1, I2) потерялось
	double I1_after = energyIntegral(q, h);
	double I2_after = energyFluxIntegral(q, h);

	if (abs(I1_after - I1_before) > 1e-3) {
		cout << "On step " << step << " dI1 = " << I1_after - I1_before;

		if (abs(I2_after - I2_before) > 1e-3) cout << ", dI2 = " << I2_after - I2_before;

		cout << endl;
	}
}

void filtration(vector<complex<double>>& q, double h, double factor, double L, size_t step) {
	size_t N = q.size();

	// вычисляем полуширину окна в индексах

	int delta = static_cast<int>(floor(L / (2.0 * h)));

	if (delta >= N / 2) {
		cout << "Can't cut all window" << endl;
		return;
	}

	int maxIt = distance(q.begin(), max_element(q.begin(), q.end(), [](auto a, auto b) { return abs(a) < abs(b); }));

	double maxAmp = abs(q[maxIt]);

	// границы окна: i_left и i_right
	int i_left = maxIt - delta;
	int i_right = maxIt + delta;

	if (i_left < 0)   i_left += N;
	if (i_right >= N) i_right -= N;

	double I1_before = energyIntegral(q, h);
	double I2_before = energyFluxIntegral(q, h);

	// Фильтруем (делим на factor) все элементы ВНЕ [i_left, i_right]
	// Учитывая, что индекс может "переваливаться" через границу массива.

	if (i_right < i_left) {
		// Считаем, что окно "обернулось" через конец массива.
		// Тогда фильтруем U[i_right..i_left].
		for (int i = i_right; i < i_left; i++) {
			q[i] *= exp(-factor * pow(abs(q[i] - maxAmp), 2));
		}
	}
	else {
		// Сначала фильтруем "левую" часть [0..i_left], 
		// потом "правую" часть [i_right..N-1].
		// (Тут подразумевается, что "окно" - это (i_left..i_right),
		//  и мы, наоборот, за его пределами ослабляем решение.)
		for (int i = 0; i < i_left; i++) {
			q[i] *= exp(-factor * pow(abs(q[i] - maxAmp), 2));
		}
		for (int i = i_right; i < N; i++) {
			q[i] *= exp(-factor * pow(abs(q[i] - maxAmp), 2));
		}
	}

	// Считаем, сколько "энергии" (I1, I2) потерялось
	double I1_after = energyIntegral(q, h);
	double I2_after = energyFluxIntegral(q, h);

	if (abs(I1_after - I1_before) > 1e-3) {
		cout << "On step " << step << " dI1 = " << I1_after - I1_before;

		if (abs(I2_after - I2_before) > 1e-3) cout << ", dI2 = " << I2_after - I2_before;

		cout << endl;
	}
}

vector<complex<double>> firstStartPER(vector<double> x, double P0, double L, double C)
{
	size_t N = x.size();
	vector<complex<double>> signal(N);

	for (size_t k = 0; k < N; k++)
	{
		signal[k] = (pow(P0, 0.5) * exp(-pow((x[k] - L * floor(static_cast<double>(x[k] / L) + 0.5)) / L, 2)) * cos(C * pow((x[k] - L * floor(static_cast<double>(x[k] / L) + 0.5)) / L, 2) / 2)
			- i * pow(P0, 0.5) * exp(-pow((x[k] - L * floor(static_cast<double>(x[k] / L) + 0.5)) / L, 2)) * sin(C * pow((x[k] - L * floor(static_cast<double>(x[k] / L) + 0.5)) / L, 2) / 2));
	}
	return signal;
}
vector<complex<double>> firstStartBSOL(vector<double> x, double k, double nu, double mu, double A, double thetta0, double x0)
{	//bright solitary wave if the double root equals to zero (i.e. R4 = 0)
	size_t N = x.size();
	vector<complex<double>> signal(N);

	for (size_t n = 0; n < N; n++)
	{

		signal[n] = pow(2 * nu * A, 0.5) * exp(i * (k * x[n] - thetta0)) / pow(-mu + pow(mu * mu - 4 * nu, 0.5) * cosh(2 * pow(nu, 0.5) * (x[n] - x0)), 0.5);
	}
	return signal;
}


vector<complex<double>> firstStartBSOLSeries(vector<double> x, double k, double a, double b, double A, double alpha, double thetta0, double x0)
{	//bright solitary wave if it is derived from equations for series
	size_t N = x.size();
	vector<complex<double>> signal(N);

	for (size_t n = 0; n < N; n++)
	{

		signal[n] = A * (cos(k * x[n] - thetta0) + i * sin(k * x[n] - thetta0)) / (a * exp(alpha * (x[n] - x0)) + b * exp(-alpha * (x[n] - x0)));
	}
	return signal;
}

vector<complex<double>> firstStartBSOLSeriesWhiteNoize(vector<double> x, double k, double a, double b, double A, double alpha, double thetta0, double x0, double coeff1)
{	//bright solitary wave if it is derived from equations for series
	size_t N = x.size();
	vector<complex<double>> signal(N);
	uniform_real_distribution<double> unif(-1., 1.);
	default_random_engine re;
	
	for (size_t n = 0; n < N; n++)
	{
		double a_random = coeff1 * unif(re);
		signal[n] = A * (1 + a_random) * (cos(k * x[n] - thetta0) + i * sin(k * x[n] - thetta0)) / (a * exp(alpha * (x[n] - x0)) + b * exp(-alpha * (x[n] - x0)));
	}
	return signal;
}
vector<complex<double>> firstStartBSOLSeriesSinNoize(vector<double> x, double k, double a, double b, double A, double alpha, double thetta0, double x0, double coeff1,  double coeff2)
{	//bright solitary wave if it is derived from equations for series
	size_t N = x.size();
	vector<complex<double>> signal(N);
	uniform_real_distribution<double> unif(-1., 1.);
	default_random_engine re;
	
	for (size_t n = 0; n < N; n++)
	{
		double a_random = coeff1 * unif(re);
		signal[n] = A * (cos(k * x[n] - thetta0) + i * sin(k * x[n] - thetta0)) / (a * exp(alpha * (x[n] - x0)) + b * exp(-alpha * (x[n] - x0))) +  A * a_random * sin(coeff2 * x[n]);
	}
	return signal;
}
vector<complex<double>> firstStartBSOLSeriesAmpNoize(vector<double> x, double k, double a, double b, double A, double alpha, double thetta0, double x0, double coeff)
{	//bright solitary wave if it is derived from equations for series
	size_t N = x.size();
	vector<complex<double>> signal(N);
	
	for (size_t n = 0; n < N; n++)
	{
		signal[n] = A * (1 + coeff) * (cos(k * x[n] - thetta0) + i * sin(k * x[n] - thetta0)) / (a * exp(alpha * (x[n] - x0)) + b * exp(-alpha * (x[n] - x0)));
	}
	return signal;
}

vector<complex<double>> firstStartBSOLSerieCollision(vector<double> x, double k, double a, double b, double A, double alpha, double thetta0, double x0, double ampl, double mn, double syg)
{	//bright solitary wave if it is derived from equations for series
	size_t N = x.size();
	vector<complex<double>> signal(N);
	uniform_real_distribution<double> unif(-1., 1.);
	default_random_engine re;
	
	for (size_t n = 0; n < N; n++)
	{
		signal[n] = A * (cos(k * x[n] - thetta0) + i * sin(k * x[n] - thetta0)) / (a * exp(alpha * (x[n] - x0)) + b * exp(-alpha * (x[n] - x0))) + ampl * exp (-pow(x[n] - mn, 2.0) / pow(syg, 2.0));
	}
	return signal;
}
vector<complex<double>> firstStartDSOLRiccati(vector<double> x, double k, double b, double A, double thetta0, double x0)
{	//bright solitary wave if it is derived from Riccati special case a = 0 eq.
	size_t N = x.size();
	vector<complex<double>> signal(N);

	for (size_t n = 0; n < N; n++)
	{
		signal[n] = A * (cos(k * x[n] - thetta0) + i * sin(k * x[n] - thetta0)) * pow(b, 0.5) * tanh(pow(b, 0.5) * (x[n] - x0));
	}
	return signal;
}


vector<complex<double>> firstStartDSOL(vector<double> x, double k, double nu, double mu, double A, double thetta0, double x0)
{	//dark solitary wave if R3 = R4 neq 0 and R1 = 0
	double alpha = mu * mu - 3 * nu;
	double beta = -mu / 3 * pow(alpha, 0.5) + mu * mu / 3 - nu;

	size_t N = x.size();
	vector<complex<double>> signal(N);
	const complex<double> i(0., 1.);
	for (size_t n = 0; n < N; n++)
	{
		signal[n] = pow(beta * A, 0.5) * exp(i * (k * x[n] - thetta0)) / pow(1 / pow(alpha, 0.5) - pow(alpha, 0.5) / ((beta - alpha) * pow(cosh(pow(beta, 0.5) * (x[n] - x0)), 2) + alpha), 0.5);
	}
	return signal;
}
void linearFM(vector<complex<double>>& signal, vector<double> omega, vector<double> x, double tau, long long int L, double chi)
{
	FT fftLin(omega);
	fftLin.DFT(signal, x, L, false);
	long long int  N = signal.size();
	const complex<double> i(0., 1.);
	for (long long int k = 0; k < N; k++)
	{
		(signal)[k] *= exp(i * tau * (-chi * pow(omega[k], 2) + pow(omega[k], 3) + pow(omega[k], 4)));
	}
	//inversed for linear
	fftLin.DFT(signal, x, L, true);
	return;
}

void nonlinearFM(vector<complex<double>> U, vector<complex<double>>& V, vector<double> omega, double tau, long long int L, double eps)
{
	long long int N = U.size();
	const complex<double> i(0., 1.);
	//vector<complex<double>> signalConj(N);

	double sigMod2 = 0.;
	double sigMod4 = 0.;
	double sigMod6 = 0.;
	double sigMod8 = 0.;


	//FT fftNonLinConj;
	//vector<complex<double>> sigConjFourier = fftNonLinConj.FFT(vector<complex<double>> {signalConj}, false);
	//FT fftNonLin;
	//vector<complex<double>> sigFourier = fftNonLin.FFT(vector<complex<double>> {*signal}, false);
	//vector<complex<double>> sigConjFourier = vector<complex<double>>(0);
	//vector<complex<double>> sigFourier = vector<complex<double>>(0);

	for (size_t k = 0; k < N; k++)
	{
		sigMod2 = pow(abs((U)[k]), 2);
		sigMod4 = pow(abs((U)[k]), 4);
		//sigMod6 = pow(abs((*signal)[k]), 6);
		//sigMod8 = pow(abs((*signal)[k]), 8);
		// 
		//FT fftInv; 
		//complex<double> sigConjInv = complex<double>(0);
		//complex<double> sigConjInv = fftInv.FFT(vector<complex<double>> { i* omega[k] * sigConjFourier[k] }, true); //inversed for linear

		//FT fftInv2;
		//complex<double> sigInv = complex<double>(0);
		//complex<double> sigInv = fftInv2.FFT(vector<complex<double>> { i* omega[k] * sigFourier[k] }, true); //inversed for linear

		(V)[k] *= exp(i * tau * (-sigMod2 - eps * sigMod4));
	}

	return;
}

int main(int argc, char** argv)
{
	long long int N = stod(argv[1]), L = N * M_PI / 48;
	//2 * M_PI;
	double T = stod(argv[2]);


	//tau < h^4/100 / |a4|
	double h = double(L) / N;
	double tau;

	double omega, v, chi, eps, alpha, a, b;

	double K = -0.25;

	double A = stod(argv[3]);

	cout << "A = " << A << endl;

	double thetta0 = 0.0, x0 = 0.0;

	//double thetta0 = stod(argv[13]) * K, x0 = stod(argv[14]);

	//double ampGauss = A * stod(argv[10]), meanGauss = stod(argv[11]), sygGauss = stod(argv[12]);

	unsigned fo  = 100;
	
	double coeff = stod(argv[9]);
	cout << "coeff  = " << coeff  << endl;

	tau = pow(h, 4) / 2;

	vector <double> x(N);
	vector <double> omegaFr(N);
	vector<complex<double>> signal = vector<complex<double>>(N);

	bool isSeries = true;
	/*cout << "Insert initial condition" << endl;
	cout << "0 - Riccati kink, 1 - Series bright soliton" << endl;

	cin >> isSeries;*/

	if (isSeries) {
		//series

		alpha = stod(argv[4]);

		a = stod(argv[5]);

		b = stod(argv[6]);

		omega = (2304 * pow(alpha, 4) * 4 * a * b + 128 * alpha * alpha * pow(A,2) - 163 * alpha * alpha * 4 * a * b - 8 * pow(A,2))/(256 * alpha * alpha * 4 * a * b);

		cout << "omega = " << omega << endl;

		eps = 24 * pow(4 * a * b * alpha, 2) / pow(A, 4);

		cout << "eps = " << eps << endl;

		chi = -(4 * A * A + 332 * a * b * alpha * alpha) / (32 * a * b * alpha * alpha);

		cout << "chi = " << chi << endl;

		v = -(1 + 4 * chi) / 8; //series

		cout << "v = " << v << endl;

		if (tau > pow(h, 2) / 2 / abs(chi)) tau = pow(h, 2) / 2 / abs(chi);
		
		cout << "tau = " << tau << endl;

	}

	long long int M = T / tau;

	cout << "M = " << M << endl;

	cout << "file freq = " << fo << endl;

	vector <double> t(M);

	for (size_t k = 0; k < M; k++)
	{
		t[k] = k * tau;
	}

	for (size_t k = 0; k < N; k++)
	{
		omegaFr[k] = 2 * M_PI * (k - double(N) / 2) / L;
		//2 * M_PI * (k - N/2) / L;
		//2 * M_PI / (N * h) * (k - double(N) / 2);	 
	}

	omegaFr[N / 2] = 0.;

	for (size_t n = 0; n < N; n++)
	{
		x[n] = static_cast<double>(-L) / 2 + n * static_cast<double>(L) / static_cast<double>(N);

		//(signal)[n] = A * complex<double>(cos( K / static_cast<double>(L) * (x[n] - L*floor(static_cast<double>(x[n]/L) + 0.5)) + 0.),
		//	sin(K / static_cast<double>(L) * (x[n] - L * floor(static_cast<double>(x[n] / L) + 0.5)) + 0.));
	}


	ofstream outputFile1, outputFile2;

	if (isSeries) {
		//series

		outputFile1.open(string(argv[7]));
		outputFile2.open(string(argv[8]));
		signal = vector<complex<double>>{ firstStartBSOLSeriesAmpNoize(x, K, a, b, A, alpha, thetta0, x0, coeff) };
	}


	//запись для 3-д графика квадрата модуля сигнала


	size_t linw = 9;
	if (outputFile1.is_open() and outputFile2.is_open())
	{
		outputFile1 << "Time step\\Spacial step" << ',' << setw(linw) << "Module squared" << endl;
		outputFile1 << ',' << ',' << setw(linw);
		outputFile2 << "Time step\\Spacial step" << ',' << setw(linw) << "Module squared" << endl;
		outputFile2 << ',' << ',' << setw(linw);
		for (size_t n = 0; n < N; n++)
		{
			outputFile1 << x[n] << ',' << setw(linw);
			outputFile2 << x[n] << ',' << setw(linw);
		}
		outputFile1 << endl;
		outputFile2 << endl;
		try {
			for (size_t p = 0; p < M; p++)
			{

				if (p % fo == 0)
				{
					cout << "Step number " << p << " complete; the energy is ";

					//if (p * tau > T - 1.)
					//{
					outputFile1 << fixed << setprecision(linw) << t[p] << ',' << setw(linw);
					outputFile2 << fixed << setprecision(linw) << t[p] << ',' << setw(linw);
					for (size_t n = 0; n < N; n++)

					{
						//double stepMod = pow(abs((signal)[n]), 2);
						outputFile1 << fixed << setprecision(linw - 2) << signal[n].real() << ',' << setw(linw);
						outputFile2 << fixed << setprecision(linw - 2) << signal[n].imag() << ',' << setw(linw);

					}

					cout << energyIntegral((signal), h) << endl;
						//<< ", the energy flux is " << energyFluxIntegral((signal), h) << endl;

					//filtration((signal), h, 2.0, L / 3, p);

					outputFile1 << endl;
					outputFile2 << endl;
					//}
				}

				vector<complex<double>> V = signal;

				//linearFM(V, omegaFr, x, tau / 2, L, chi);

				//модули считаются по старым компонентам
				nonlinearFM(signal, V, omegaFr, tau, L, eps);

				linearFM(V, omegaFr, x, tau, L, chi);

				signal = V;

				/*for (ptrdiff_t k = 0; k < 5; ++k) {
					signal[5 - k - 1] = signal[N - k - 1];
					signal[N - 5 + k] =  signal[k];
				}*/

			}
		}
		catch (exception& e) {
			cerr << e.what() << endl;
			exit(1000);
		}
	}
	else
	{
		cerr << "Couldn't open file to write" << endl;
		exit(1);
	}
	outputFile1.close();
	outputFile2.close();
}

void  printFile(vector<complex<double>> var, string path, long long int L)
{
	ofstream file("log/" + path);
	size_t N = var.size();
	vector<double> x(N);
	for (size_t n = 0; n < N; n++)
	{
		x[n] = static_cast<double>(-L) / 2 + n * static_cast<double>(L) / static_cast<double>(N);
	}
	if (file.is_open())
	{
		file << "Order		" << "Re		" << "Im		" << endl;
		for (size_t k = 0; k < N; k++)
		{
			file << x[k] << "\t" << setw(12) << var[k].real() << "\t" << setw(12) << var[k].imag() << endl;
		}
	}
	else
	{
		cerr << "Couldn't open file to write" << endl;
		exit(1);
	}
	file.close();
	return;
}
