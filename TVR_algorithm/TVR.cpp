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

const int N = 16, M = pow(2, 13); //Число шагов // количество моделируемых траекторий
double S0 = 100, sigma = 0.08, mu = 0.05, dt = 1.0 / N, //начальная цена, волатильность, текущая процентная ставка
K = 101, pf1 = (mu - sigma * sigma * 0.5) * dt, pf2 = sigma * sqrt(dt); //strike price и просто промежуточные вычисления

Sobol sobol(N);
double sx[N];
double randFT = unif_dis(gen);














//фиксить??
//Функция которая считает значение базисных многочленов
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


//фиксить??
//Функция считающая значени платежной
double pay_func(double x, double time)
{
	double temp = exp(-0.05 * (time / N));
	return
		temp * max(0.0, (K - x));
}


//фиксить
//Функция создающая вектор с 16 нормально распределенными величинами
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

//фиксить
//Функция генерирующая 1 траекторию цены
void prices(vector<double>& price)
{
	price.resize(N + 1);
	price[0] = S0;
	for (size_t i = 1; i < N + 1; i++)
		price[i] = price[i - 1] * exp(pf1 + pf2 * norm_dis(gen));
}

//фиксить
//Функция генерирущая 1 траекторию цены с помощью КМК
void prices_withQMC(vector<double>& price)
{
	price.resize(N + 1);
	price[0] = S0;
	vector<double> normal_numbers(N);
	normal_numbers = norm_dis_vec();
	for (size_t i = 1; i <= N; i++)
		price[i] = price[i - 1] * exp(pf1 + pf2 * normal_numbers[i - 1]);
}


//фиксить
//функция генерирующая матрицу M * N цен, и заполняет нужный для вычислений вектор W параметры start и end чтобы распараллелить этот процесс
void total_price_modeling(vector<vector<double>>& price_matrix, vector<double>& W, int start, int end)
{
	for (int i = start; i < end; i++)
	{
		prices_withQMC(price_matrix[i]); //Цены формируются с использованием чисел соболя
		//prices(price[i]); //Цены формируются просто случайными величинами
		W[i] = pay_func(price_matrix[i].back(), N);
	}
}


//Переделывать полностью
//Функция вычисляющая вектор кэфов бета Исправить название параметра time1 -> time
vector<double> beta(vector<double>& beta, vector<double> w, vector<vector<double>> price, int time1)
{
	//заполняем матрицу мнк
	vector<vector<double>> mnk(7);
	time_t start, end;
	time(&start);
	for (int i = 0; i < 7; i++)
	{
		mnk[i].resize(8);
		for (int j = 0; j < 8; j++)
		{
			double a_ij = 0;
			for (int traectory = 0; traectory < M; traectory++)
			{
				if (j != 7)
				{
					a_ij += basis_func(price[traectory][time1], i) * basis_func(price[traectory][time1], j);
				}
				else
				{
					a_ij += basis_func(price[traectory][time1], i) * w[traectory];
				}
			}
			mnk[i][j] = a_ij;
		}
	}


	time(&end);
	fout << difftime(end, start) << endl;


	//beta = gausswmec(mnk);
	return beta;
}

//Функция выдающая оценку
void opt_price(double& x)
{
	vector<vector<double>> price(M);
	vector<double> W(M);
	vector<double> beta_j(7);
	double option_price = 0;
	for (int i = 0; i < M; i++)
	{
		prices_withQMC(price[i]); //цены формируются с использованием чисел соболя
		//prices(price[i]); //цены формируются просто случайными величинами
		W[i] = pay_func(price[i].back(), N);
	}

	/*for (int i = 0; i < price.size(); i++) //вывод траекторий в файл
	{
		for (int j = 0; j < price[i].size(); j++)
		{
			fout << price[i][j] << "  ";
		}
		fout << endl;
	}*/

	for (int time = N; time >= 1; time--) //дебил, здесь 0 должен быть
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
			if (time == 1) //и здесь 0, мудак
			{
				option_price += W[traectory];
			}
		}

	}
	x = option_price / M;
}

//Функция выдающая среднее, сигму, дисперсию
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
	fout.open("output.txt");
	double sigma = 0;
	double var = 0;
	double mean = 0;
	//stats(1, sigma, var, mean);
	cout << "Mean:   " << mean << "   Sigma:   " << sigma << "    Variance:     " << var;
	//Сетка при 1000 шагов  2.2048
	//увеличить до 40 мб
	return 0;
}