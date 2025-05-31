#include <iostream>
#include "DFT.h"
#include <vector>
#include <complex>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#define _USE_MATH_DEFINES
const complex<double> i(0., 1.);
#include <math.h>
using namespace std;
void  printFile(vector<complex<double>> var, string path, long long int L);
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


vector<complex<double>> firstStartBSOL2(vector<double> x, double k, double a, double b, double A, double alpha, double thetta0, double x0)
{	//bright solitary wave if the double root equals to zero (i.e. R4 = 0)
	size_t N = x.size();
	vector<complex<double>> signal(N);

	for (size_t n = 0; n < N; n++)
	{

		signal[n] = A * (cos(k * x[n] - thetta0) + i * sin (k * x[n] - thetta0)) / (a*exp(alpha * (x[n] - x0)) + b * exp(-alpha * (x[n] - x0)));
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
void linearFM(vector<complex<double>>* signal, double tau, long long int L, double a1, double a2, double a3, double a4)
{
	DFT dft(*signal);
	*signal = dft.Y;
	long long int  N = signal->size();
	vector<double> omega(N);
	const complex<double> i(0., 1.);
	for (long long int k = 0; k < N; k++)
	{
		omega[k] = - M_PI + 2 * M_PI * k / N;
		(*signal)[k] *= exp(i * tau * (- a1 * omega[k] - a2 * pow(omega[k], 2) + a3 * pow(omega[k], 3) + a4 * pow(omega[k], 4))); 
	}
	DFT dftInv(*signal, true); //inversed for linear
	*signal = dftInv.Y;
	return;
}
void nonlinearScheme(vector<complex<double>>* signal, double tau, long long int L, double b1, double b2, double b3)
{
	size_t N = signal->size();
	const complex<double> i(0., 1.);
	vector<double> omega(N);
	vector<complex<double>> signalConj(N);
	for (size_t k = 0; k < N; k++)
	{
		omega[k] = -M_PI + 2 * M_PI * k / N;
			//2 * M_PI * k / N;
			//-M_PI + 2 * M_PI * k / N;

		 
	}
}
void nonlinearFM(vector<complex<double>>* signal, double tau, long long int L, double b1, double b2, double b3)
{
	size_t N = signal->size();
	const complex<double> i(0., 1.);
	vector<double> omega(N);
	vector<complex<double>> signalConj(N);
	for (size_t k = 0; k < N; k++)
	{
		omega[k] = - M_PI + 2 * M_PI * k / N;

		signalConj[k] = 0;
			//conj((*signal)[k]);
	}

	double sigMod2 = 0.;
	double sigMod4 = 0.;

	for (size_t k = 0; k < N; k++)
	{
		sigMod2 += pow(abs((*signal)[k]), 2);
		sigMod4 += pow(abs((*signal)[k]), 4);
	}
	//DFT dft(vector<complex<double>> {signalConj});
	//DFT dft2(vector<complex<double>> {*signal});
	vector<complex<double>> sigConjFourier = vector<complex<double>>(0);
		//dft.Y;
	vector<complex<double>> sigFourier = vector<complex<double>>(0);
		//dft2.Y;

	for (size_t k = 0; k < N; k++)
	{
		//DFT dftInv(vector<complex<double>> { i* omega[k] * sigConjFourier[k] }, true); //inversed for linear
		complex<double> sigConjInv = complex<double>(0);
			//dftInv.Y[0];
		//DFT dftInv2(vector<complex<double>> { i* omega[k] * sigFourier[k] }, true); //inversed for linear
		complex<double> sigInv = complex<double>(0);
			//dftInv.Y[0];

		(*signal)[k] *= exp(i * tau * (-b1 * sigMod2 - b2 * sigMod4 - i * b3 * (*signal)[k] * sigConjInv - i * b3 * signalConj[k] * sigInv));
	}

	return;
}

int main(int argc, char** argv)
{
	long long int N = 400, L = 4 * M_PI;
	double T = 2.;
	
		
	//tau < h^2/2 / |a2|
	

	double sigK = 1., sigPh = M_PI/2, mod = sqrt(0.5), h = double(L) / N, tau = 9e-5;

	double R_1 = 0.5, R_2 = 0.05, R_3 = 0.0, R_4 = 0.0;
	
	double c_sq = (R_1 - R_3)* (R_2 - R_4), A = 4., nu = 0.025, mu = -0.5, thetta0 = 0., x0 = 0.;
	long long int M = T / tau;
	vector <double> t(M);
	vector <double> x(N+1);
	vector<complex<double>>* signal = new vector<complex<double>> (N+1);

	double k = 1., a = 2., b = 2., alpha = 2.;

	//k = 2.3404970685324216; //как в статье для периодич волны

	double a1 =  0.2;
		//1.,

	double a2 =  0.1;
		//1.;
	//1.

	double a4 =  0.05;

	//-0.001;
	//0.001;

	double v = 8 * pow(k, 3) * a4 + a1 + 2 * k * a2;
	
	double omega = 3 * a4 * pow(k ,4) + a2 * pow(k, 2) - 6 * a4 * pow(k, 2) + a1 * k - a2 - a4;

	double a3 = -4 * k * a4;
		//-4 * k * a4;
		//0.001;

	double b1 = -2 * 4 * a * b * (6 * a4 * pow(k, 2) + a2 + 10 * a4) / (pow(A, 2));
		//- 2 * 4 * a * b  * (6 * a4 * pow(k, 2) + a2 + 10 * a4) / (pow(A, 2));
		//5;


	double b2 = 24 * a4 * pow(4 * a * b, 2) / pow(A, 4);
	//5;
		//24 * a4 * pow(4 * a * b, 2)/ pow(A, 4);
		//5;

	double b3 = 0.00;
		//0.;
		//120 * mu * a4 / pow(A, 3);

	//L = 120 * M_PI / k;

	for (size_t n = 0; n < M; n++)
	{
		t[n] = n * tau;
			 
	}
	
	for (size_t n = 0; n < N+1; n++)
	{
		x[n] = static_cast<double>(-L) / 2 +  n * static_cast<double>(L) / static_cast<double>(N);
		 
		(*signal)[n] = mod * complex<double>(sin(1.0 * M_PI / static_cast<double>(L) * (x[n] - L*floor(static_cast<double>(x[n]/L) + 0.5)) + sigPh),
										0.);
	}
	
 
	//signal = new vector<complex<double>> { firstStartBSOL(x, k, nu, mu, A, thetta0, x0) };

	signal = new vector<complex<double>>{ firstStartBSOL2(x, k, a, b, A, alpha, thetta0, x0) };
 
	//запись для 3-д графика квадрата модуля сигнала
	ofstream file("log/bsolEnvTest7e.csv");
	size_t linw = 6;
	if (file.is_open())
	{
		file << "Time step\\Spacial step" << ',' << setw(linw) << "Module squared" << endl;
		file << ',' << ',' << setw(linw);
		for (size_t n = 0; n < N; n++)
		{
			file <<  x[n] << ',' << setw(linw);
		}
		file << endl;

		for (size_t p = 0; p < M; p++)
		{

			if (p % 1 == 0)
			{
				cout << "Step number " << p << " complete" << endl;
				//if (p * tau > T - 1.)
				//{
					file << fixed << setprecision(linw) << t[p] << ',' << setw(linw);

					for (size_t n = 0; n < N; n++)
									
					{
						file << fixed << setprecision(linw - 2) << (pow(abs((*signal)[n].real()), 2) + pow(abs((*signal)[n].imag()), 2)) << ',' << setw(linw);
						
					}

					file << endl;
				//}
			}

			linearFM(signal, tau / 2, L, a1, a2, a3, a4);

			nonlinearFM(signal, tau, L, b1, b2, b3);

			linearFM(signal, tau / 2, L, a1, a2, a3, a4);

		}
	}
	else
	{
		cerr << "Couldn't open file to write" << endl;
		exit(1);
	}
	file.close();


	/*
	signal[0] = 2;
	signal[1] = complex<double>(-2, -2);
	signal[2] = complex<double>(0, -2);
	signal[3] = complex<double>(4, 4);
	*/
	
	//DFT dft(*signal);
	//dft.printFile("dftTest1.txt");
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
