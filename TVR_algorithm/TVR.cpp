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
#include <armadillo>


using namespace std;

normal_distribution<double> norm_dis(0, 1);
random_device rd;
mt19937 gen(rd());

ofstream fout;
typedef std::vector<double> stdvec;

uniform_real_distribution<double> unif_dis(0, 1);

const int N = 16, M = 8192, basis_fun_am = 7; //Число шагов // количество моделируемых траекторий
double S0 = 100, sigma = 0.08, mu = 0.05, dt = 1.0 / N, //начальная цена, волатильность, текущая процентная ставка
K = 101, pf1 = (mu - sigma * sigma * 0.5) * dt, pf2 = sigma * sqrt(dt); //strike price и просто промежуточные вычисления
double eps = 0.00001;

uniform_int_distribution<> unif_dis_int(0, M - 1);

Sobol sobol(N);
double sx[N];

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

//возвращает вектор с примененным к x базисными функциями
vector<double> basis(double x)
{
	vector<double> res(basis_fun_am);

	for (int i = 0; i < basis_fun_am; i++)
	{
		res[i] = basis_func(x, i);
	}

	return res;
}

//Функция считающая значени платежной
double pay_func(double x, double time)
{
	double temp = exp(-0.05 * (time / N));
	return
		temp * max(0.0, (K - x));
}

//Функция создающая вектор с 16 нормально распределенными величинами
vector<double> norm_dis_vec(double rand_num)
{
	sobol.next(sx);
	for (int i = 0; i < N; i++)
	{
		double& sxi = sx[i];
		//double rand_num = unif_dis(gen);
		sxi += rand_num;
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

//функция генерирующая матрицу M * N цен, и заполняет нужный для вычислений вектор W, использует квази числа
void price_modeling_QMC(arma::mat& price_matrix, arma::mat& W)
{
	vector<double> normal_numbers(N);
	int step = M / 5;
	double rand_num = unif_dis(gen);
	for (int i = 0; i < M; i++)
	{
		if (i >= step)
		{
			step += M / 5;
			rand_num = unif_dis(gen);
		}
		price_matrix(i, 0) = S0;
		normal_numbers = norm_dis_vec(rand_num);

		for (int j = 1; j < N + 1; j++)
		{
			price_matrix(i, j) = price_matrix(i, j - 1) * exp(pf1 + pf2 * normal_numbers[j - 1]);
		}

		W(i) = pay_func(price_matrix(i, N), N);
	}
}

//функция генерирующая матрицу M * N цен, и заполняет нужный для вычислений вектор W
void price_modeling(arma::mat& price_matrix, arma::mat& W)
{
	for (int i = 0; i < M; i++)
	{
		price_matrix(i, 0) = S0;

		for (int j = 1; j < N + 1; j++)
		{
			price_matrix(i, j) = price_matrix(i, j - 1) * exp(pf1 + pf2 * norm_dis(gen));
		}

		W(i) = pay_func(price_matrix(i, N), N);
	}
}

arma::mat construct_at_time(arma::mat prices, int time)
{
	arma::mat result = arma::mat(M, basis_fun_am);
	arma::rowvec temp;
	for (int i = 0; i < result.n_rows; i++)
	{
		temp = arma::rowvec(basis(prices(i, time)));
		result.row(i) = temp;
	}
	return result;
}

void least_square_regr(arma::mat& beta, arma::mat X, arma::mat W, int time)
{
	arma::mat Xt = X.t();
	beta = arma::pinv(Xt*X)*Xt*W;
}

arma::mat predict(arma::mat X, arma::mat coef) 
{
	return X * coef;
}

//Функция выдающая оценку
void opt_price(double& x)
{
	arma::mat price = arma::mat(M, N + 1);
	arma::mat W = arma::mat(M, 1);
	arma::mat beta_j = arma::mat(basis_fun_am, 1);
	double option_price = 0;

	price_modeling_QMC(price, W);
	//price_modeling(price, W);
	for (int time = N; time >= 0; time--)
	{
		arma::mat prices_at_time = construct_at_time(price, time);
		least_square_regr(beta_j, prices_at_time, W, time);
		arma::mat pred = predict(prices_at_time, beta_j);
		for (int trajectory = 0; trajectory < M; trajectory++)
		{
			W[trajectory] = max(pay_func(price(trajectory, time), time), pred(trajectory, 0));
			if (time == 0)
			{
				option_price += W[trajectory];
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
	time(&start);
	for (int i = 0; i < numberOfRuns; i++)
	{
		double xi = 0;
		opt_price(xi);
		sqSum += xi * xi;
		mean += xi;
	}
	time(&end);
	mean *= 1.0 / numberOfRuns;
	sqSum *= 1.0 / numberOfRuns;
	var = sqSum - (mean * mean);
	sigma = sqrt(var);
	allRunTime += difftime(start, end);
	cout << "Average running time  " << allRunTime / numberOfRuns << endl;
}

int main()
{
	setlocale(LC_ALL, "russian");
	//thread threads[2];

	fout.open("output.txt");
	double sigma = 0;
	double var = 0;
	double mean = 0;

	stats(20, sigma, var, mean);
	cout << "Mean:   " << mean << "   Sigma:   " << sigma << "    Variance:     " << var;
	//Сетка при 1000 шагов  2.2048
	//увеличить до 40 мб
	return 0;
}