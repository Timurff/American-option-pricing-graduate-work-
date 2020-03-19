#include "pch.h"
#include "sobol.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <random>
#include <time.h> 
#include <thread>
#include <chrono>
#include <fstream>


using namespace std;

normal_distribution<double> norm_dis(0, 1);
random_device rd;
mt19937 gen(rd());

ofstream fout;

uniform_real_distribution<double> unif_dis(0, 1);

const int N = 16, M = 8192; //����� ����� // ���������� ������������ ����������
double S0 = 100, sigma = 0.08, mu = 0.05, dt = 1.0 / N, //��������� ����, �������������, ������� ���������� ������
K = 101, pf1 = (mu - sigma * sigma * 0.5) * dt, pf2 = sigma * sqrt(dt); //strike price � ������ ������������� ����������

Sobol sobol(N);
double sx[N];
double randFT = unif_dis(gen);

/*class TVR
{


public:
	double price[M][N + 1];

};*/


//�������??
//������� ������� ������� �������� �������� �����������
double basis_func(double x, int power)
{
	if (power == 0) return 1;
	if (power < 5) return  pow((x - 101), (power));
	if (power == 5) return max((x - 101), (double)0);
	else
	{
		double result = max((x - 101), (double)0);
		return result * result;
	}
}


//�������??
//������� ��������� ������� ���������
double pay_func(double x, double time)
{
	double temp = exp(-0.05 * (time / N));
	return
		temp * max(0.0, (K - x));
}


//�������
//������� ��������� ������ � 16 ��������� ��������������� ����������
vector<double> norm_dis_vec()
{
	sobol.next(sx);
	for (int i = 0; i < N; i++)
	{
		double& sxi = sx[i];
		sxi += randFT;
		if (sxi > 1) sxi -= 1;
	}
	vector<double> nDis(N);
	for (int i = 0; i < N; i += 2)
	{
		double p;
		p = sqrt(-2 * log(sx[i]));
		double teta = 2 * M_PI * sx[i + 1];
		nDis[i] = p * cos(teta);
		nDis[i + 1] = p * sin(teta);
	}
	return nDis;
}

//�������
//������� ������������ 1 ���������� ����
void prices(vector<double>& price)
{
	price.resize(N + 1);
	price[0] = S0;
	for (size_t i = 1; i < N + 1; i++)
		price[i] = price[i - 1] * exp(pf1 + pf2 * norm_dis(gen));
}

//�������
//������� ����������� 1 ���������� ���� � ������� ���
void prices_withQMC(double price)
{

	/*price.resize(N + 1);
	price[0] = S0;
	vector<double> normal_numbers(N);
	normal_numbers = norm_dis_vec();
	for (size_t i = 1; i <= N; i++)
		price[i] = price[i - 1] * exp(pf1 + pf2 * normal_numbers[i - 1]);*/
}



//�������
//������� ������������ ������� M * N ���, � ��������� ������ ��� ���������� ������ W
void total_price_modeling(double** price_matrix, vector<double>& W)
{
	vector<double> normal_numbers(N);
	for (int i = 0; i < M; i++)
	{
		price_matrix[i][0] = S0;
		normal_numbers = norm_dis_vec();
		
		for (int j = 1; j < N + 1; j++)
		{
			price_matrix[i][j] = price_matrix[i][j - 1] * exp(pf1 + pf2 * normal_numbers[j - 1]);
		}

		W[i] = pay_func(price_matrix[i][N], N);
	}
}


//������������ ���������
//������� ����������� ������ ����� ���� ��������� �������� ��������� time1 -> time
vector<double> beta(vector<double>& beta, vector<double> w, vector<vector<double>> price, int time1)
{
	//��������� ������� ���
	vector<vector<double>> mnk(7);
	time_t start, end;
	time(&start);



	time(&end);
	fout << difftime(end, start) << endl;


	//beta = gausswmec(mnk);
	return beta;
}

//������� �������� ������
void opt_price(double& x)
{
	double** price = new double*[M];
	for (int n = 0; n <= M; n++) 
		price[n] = new double[N + 1];
	//double price[M][N + 1];
	vector<double> W(M);
	vector<double> beta_j(7);

	double option_price = 0;

	total_price_modeling(price, W); //��� ��� ����� �����

	//����� ������ ����� ������� ���� ������� �� ��������� 
	for (int time = N; time >= 1; time--) //�����, ����� 0 ������ ����
	{
		beta(beta_j, W, price, time);
		for (int traectory = 0; traectory < M; traectory++)
		{
			double q = 0;
			for (int k = 0; k < 7; k++)
			{
				q += beta_j[k] * basis_func(price[traectory][time], k);
			}
			W[traectory] = max(pay_func(price[traectory][time], time), q);
			if (time == 1) //� ����� 0, �����
			{
				option_price += W[traectory];
			}
		}

	}
	x = option_price / M;
}

//������� �������� �������, �����, ���������
void stats(int numberOfRuns, double& sigma, double& var, double& mean)
{
	double sqSum = 0;
	time_t start, end;
	double allRunTime = 0;
	for (int i = 0; i < numberOfRuns; i++)
	{
		time(&start);
		double xi = 0;
		opt_price(xi);
		//cout << i + 1 << ")  " << xi << endl;
		//randFT = unif_dis(gen);
		sqSum += xi * xi;
		mean += xi;
		time(&end);
		allRunTime += difftime(start, end);
		fout << endl << "new run " << endl;
	}
	mean *= 1.0 / numberOfRuns;
	sqSum *= 1.0 / numberOfRuns;
	var = sqSum - (mean * mean);
	sigma = sqrt(var);
	cout << "Average time  8k, trajectories " << allRunTime / numberOfRuns << endl;
}

int main()
{
	setlocale(LC_ALL, "russian");
	//thread threads[2];

	fout.open("output.txt");
	double sigma = 0;
	double var = 0;
	double mean = 0;
	opt_price(var);
	//stats(1, sigma, var, mean);
	cout << "Mean:   " << mean << "   Sigma:   " << sigma << "    Variance:     " << var;
	//����� ��� 1000 �����  2.2048
	//��������� �� 40 ��
	return 0;
}


/*
for (int i = 0; i < price.size(); i++) //����� ���������� � ����
{
	for (int j = 0; j < price[i].size(); j++)
	{
		fout << price[i][j] << "  ";
	}
	fout << endl;
}
*/