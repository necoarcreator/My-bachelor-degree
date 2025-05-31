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

#define _USE_MATH_DEFINES
const complex<double> i(0., 1.);
#include <math.h>
using namespace std;

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
		dq[k] = (q[k + 1] - q[k - 1])/ (2 * h * h);
	}
	dq[q.size() - 1] = (q[0] - q[q.size() - 2]) / (2 * h * h);


	for (size_t k = 0; k < q.size(); ++k) res += dq[k] * conj(q[k]) - conj(dq[k]) * q[k];

	return real(-i * h * res);
}

void filtration(vector<complex<double>>& q, double h, double factor, double L, size_t step) {
	size_t N = q.size();

	// ��������� ���������� ���� � ��������
	
	int delta = static_cast<int>(floor(L / (2.0 * h)));

	if (delta >= N / 2) {
		cout << "Can't cut all window" << endl;
		return;
	}

	int maxIt = distance(q.begin(), max_element(q.begin(), q.end(), [](auto a, auto b) { return abs(a) < abs(b); }));

	double maxAmp = abs(q[maxIt]);

	// ������� ����: i_left � i_right
	int i_left = maxIt - delta;
	int i_right = maxIt + delta;

	if (i_left < 0)   i_left += N;
	if (i_right >= N) i_right -= N;

	double I1_before = energyIntegral(q, h);
	double I2_before = energyFluxIntegral(q, h);

	// ��������� (����� �� factor) ��� �������� ��� [i_left, i_right]
	// ��������, ��� ������ ����� "��������������" ����� ������� �������.
	
	if (i_right < i_left) {
		// �������, ��� ���� "����������" ����� ����� �������.
		// ����� ��������� U[i_right..i_left].
		for (int i = i_right; i < i_left; i++) {
			q[i] *= exp(-factor * pow(abs(q[i] - maxAmp), 2));
		}
	}
	else {
		// ������� ��������� "�����" ����� [0..i_left], 
		// ����� "������" ����� [i_right..N-1].
		// (��� ���������������, ��� "����" - ��� (i_left..i_right),
		//  � ��, ��������, �� ��� ��������� ��������� �������.)
		for (int i = 0; i < i_left; i++) {
			q[i] *= exp(-factor * pow(abs(q[i] - maxAmp), 2));
		}
		for (int i = i_right; i < N; i++) {
			q[i] *= exp(-factor * pow(abs(q[i] - maxAmp), 2));
		}
	}

	// �������, ������� "�������" (I1, I2) ����������
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

		signal[n] = A * (cos(k * x[n] - thetta0) + i * sin (k * x[n] - thetta0)) / (a*exp(alpha * (x[n] - x0)) + b * exp(-alpha * (x[n] - x0)));
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


vector<complex<double>> firstStartVTA(vector<double> x, double k, double M0, double M1, double mu, double z_0)
{	//dark solitary wave if R3 = R4 neq 0 and R1 = 0

	size_t N = x.size();
	vector<complex<double>> signal(N);
	const complex<double> i(0., 1.);
	for (size_t n = 0; n < N; n++)
	{
		signal[n] = sqrt(complex<double>(M0 + M1 / (1.0 + exp(mu * (x[n] - z_0))) - M1 / pow(1.0 + exp(mu * (x[n] - z_0)), 2.0))) * exp(i * k * x[n]);
	}
	return signal;
}
void linearFM(vector<complex<double>>& signal, vector<double> omega,  vector<double> x, double tau, long long int L)
{
	FT fftLin(omega);
	fftLin.DFT(signal, x, L, false);
	long long int  N = signal.size();
	const complex<double> i(0., 1.);
	for (long long int k = 0; k < N; k++)
	{
		(signal)[k] *= exp(i * tau * (1. * pow(omega[k], 2) + 0 * pow(omega[k], 3) + 0 * pow(omega[k], 4)));
	}
	//inversed for linear
	fftLin.DFT(signal, x, L, true);
	return;
}

void nonlinearFM(vector<complex<double>> U, vector<complex<double>>& V, vector<double> omega, double tau, long long int L, double eps2, double eps3)
{
	long long int N = U.size();
	const complex<double> i(0., 1.);
	//vector<complex<double>> signalConj(N);

	double sigMod2 = 0.;
	double sigMod4 = 0.;
	double sigMod6 = 0.;
	double sigMod8 = 0.;

	for (size_t k = 0; k < N; k++)
	{
		sigMod2 = pow(abs((U)[k]), 2);
		sigMod4 = pow(abs((U)[k]), 4);
		sigMod6 = pow(abs((U)[k]), 6);
		(V)[k] *= exp(i * tau * (-1. *  sigMod2 - eps2 * sigMod4 - eps3 * sigMod6));
	}

	return;
}

int main(int argc, char** argv)
{
	long long int N = 256, L = N * M_PI / 70;
	double T =   2.01;
	double h = double(L) / N; 
	double tau;

	vector <double> x(N);
	vector <double> omegaFr(N);
	vector<complex<double>> signal = vector<complex<double>> (N);

	double k = 1.6, M0 = -1.48, M1 = 6.16, z_0 = 0.;

	double mu = sqrt(M1 / M0 / (M1 + 6 * M0));
	
	double eps2 = 3. / 4. / M0 * (M1 / (M1 + 6 * M0) - 2.);

	double eps3 = 4. / M0 / (M1 + 6 * M0);

	tau = pow(h, 2) / 10;

	long long int M = T / tau;

	cout << "M = " << M << endl;

	vector <double> t(M);

	for (size_t k = 0; k < M; k++)
	{
		t[k] = k * tau;
	}

	for (size_t k = 0; k < N; k++)
	{
		omegaFr[k] = 2 * M_PI * (k - double(N) / 2) / L;
	}

	omegaFr[N / 2] = 0.;

	for (size_t n = 0; n < N; n++)
	{
		x[n] = static_cast<double>(-L) / 2 + n * static_cast<double>(L) / static_cast<double>(N);
	}


	ofstream outputFile;

	outputFile.open("test1VTA.csv");

	signal = vector<complex<double>>{ firstStartVTA(x,k, M0, M1, mu, z_0) };

	//������ ��� 3-� ������� �������� ������ �������
	size_t linw = 9;
	if (outputFile.is_open())
	{
		outputFile << "Time step\\Spacial step" << ',' << setw(linw) << "Module squared" << endl;
		outputFile << ',' << ',' << setw(linw);
		for (size_t n = 0; n < N; n++)
		{
			outputFile <<  x[n] << ',' << setw(linw);
		}
		outputFile << endl;

		try {
			for (size_t p = 0; p < M; p++)
			{

				if (p % 1  == 0)
				{
					cout << "Step number " << p << " complete; the energy is ";
					outputFile << fixed << setprecision(linw) << t[p] << ',' << setw(linw);

					for (size_t n = 0; n < N; n++)

					{
						double stepMod = (pow(abs((signal)[n]), 2));

						outputFile << fixed << setprecision(linw - 2) << stepMod << ',' << setw(linw);

					}

					cout << energyIntegral((signal), h) << ", the energy flux is " << energyFluxIntegral((signal), h) << endl;
					
					//filtration((signal), h, 2.0, L / 3, p);
					
					outputFile << endl;
				}

				vector<complex<double>> V = signal;
				
				//linearFM(V, omegaFr, x, tau / 2, L, chi);

				//������ ��������� �� ������ �����������
				nonlinearFM(signal, V, omegaFr, tau, L, eps2, eps3);

				linearFM(V, omegaFr, x,  tau, L);
				 
				signal = V;

				/*for (ptrdiff_t k = 0; k < 10; ++k) {
					signal[k] = signal[N - k - 1]; 
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
	outputFile.close();
 }
