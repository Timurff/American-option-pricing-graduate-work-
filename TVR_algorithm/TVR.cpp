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

enum { BASIS_1, BASIS_2 };
enum { QUASI, NOT_QUASI };

normal_distribution<double> norm_dis(0, 1);
random_device rd;
mt19937 gen(rd());

ofstream fout;

uniform_real_distribution<double> unif_dis(0, 1);


const int BASIS = BASIS_2, mode = QUASI;


int number_of_runs = 20;
const int N = 16, M = 8092 * 2; //Число шагов // количество моделируемых траекторий
double S0 = 100, sigma = 0.08, mu = 0.05, T = 1., dt = T / N, //начальная цена, волатильность, текущая процентная ставка
K = 101, pf1 = (mu - sigma * sigma * 0.5) * dt, pf2 = sigma * sqrt(dt); //strike price и просто промежуточные вычисления
double eps = 0.00001;
int basis_fun_am = 7;

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

double basis_func_2(double x, int power)
{
	if (power == 0) 
	{
		return 1;
	}
	else
	{
		return  pow((x - 101), (power));
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

vector<double> basis_2(double x)
{
	vector<double> res(basis_fun_am);

	for (int i = 0; i < basis_fun_am; i++)
	{
		res[i] = basis_func(x, i);
	}

	return res;
}

//Функция считающая значени платежной
double pay_func(double x)
{
	return max(0.0, (K - x));
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

		W(i) = pay_func(price_matrix(i, N));
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

		W(i) = pay_func(price_matrix(i, N));
	}
}

arma::mat construct_at_time(int rows, arma::vec prices)
{
	arma::mat result = arma::mat(rows, basis_fun_am);
	arma::rowvec temp;
	for (int i = 0; i < result.n_rows; i++)
	{
		temp = BASIS == BASIS_1 ? arma::rowvec(basis(prices(i))) : arma::rowvec(basis_2(prices(i)));
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
	basis_fun_am = BASIS == BASIS_1 ? 7 : 5;

	arma::mat price = arma::mat(M, N + 1);
	arma::mat W = arma::mat(M, 1);
	arma::mat beta_j = arma::mat(basis_fun_am, 1);
	double option_price = 0;

	if (mode == QUASI)
	{
		price_modeling_QMC(price, W);
	}
	else
	{
		price_modeling(price, W);
	}

	for (int time = N - 1; time >= 1; time--)
	{
		W = W * exp(-mu * dt);
		arma::vec price_now = price.col(time);
		arma::vec exercise(price_now.n_elem);

		for (int i = 0; i < exercise.n_elem; i++)
		{
			exercise(i) = pay_func(price_now(i));
		}

		
		arma::uvec itm_paths = arma::find(exercise);
		arma::vec price_now_slice(itm_paths.n_elem);
		arma::vec W_now(itm_paths.n_elem);

		for (int i = 0; i < itm_paths.n_elem; i++)
		{
			price_now_slice(i) = price_now(itm_paths(i));
			W_now(i) = W(itm_paths(i));
		}

		arma::mat prices_at_time = construct_at_time(itm_paths.n_elem, price_now_slice);

		least_square_regr(beta_j, prices_at_time, W_now, time);
		arma::mat continuation = predict(prices_at_time, beta_j);

		for (int trajectory = 0; trajectory < itm_paths.n_elem; trajectory++)
		{
			W(itm_paths(trajectory)) = max(exercise(itm_paths(trajectory)), continuation(trajectory, 0));
		}
	}

	W = W * exp(-mu * dt);
	for (int i = 0; i < M; i++)
	{
		option_price += W(i);
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

	stats(number_of_runs, sigma, var, mean);
	cout << "Mean:   " << mean << "   Sigma:   " << sigma << "    Variance:     " << var;
	//Сетка при 1000 шагов  2.2048
	//увеличить до 40 мб
	return 0;
}