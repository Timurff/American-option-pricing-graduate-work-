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

//считает норму разсности двух векторов
double norm_2(vector<double> v1, vector<double> v2)
{
	double res = 0;
	for (int i = 0; i < v1.size(); i++)
	{
		double temp = (v1[i] - v2[i]);
		res += temp * temp;
	}

	return sqrt(res);
}

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

double MSE(arma::mat& price, vector<double> W, vector<double> beta, int time)
{
	double res = 0;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < basis_fun_am; j++)
		{
			double temp = basis(price(i, time))[j] * beta[j] - W[i];
			res += temp * temp;
		}
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

//фиксить
//Функция создающая вектор с 16 нормально распределенными величинами
vector<double> norm_dis_vec()
{
	sobol.next(sx);
	for (int i = 0; i < N; i++)
	{
		double& sxi = sx[i];
		double rand_num = unif_dis(gen);
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
void price_modeling_QMC(arma::mat& price_matrix, vector<double>& W)
{
	vector<double> normal_numbers(N);
	for (int i = 0; i < M; i++)
	{
		price_matrix(i, 0) = S0;
		normal_numbers = norm_dis_vec();

		for (int j = 1; j < N + 1; j++)
		{
			price_matrix(i, j) = price_matrix(i, j - 1) * exp(pf1 + pf2 * normal_numbers[j - 1]);
		}

		W[i] = pay_func(price_matrix(i, N), N);
	}
}

//функция генерирующая матрицу M * N цен, и заполняет нужный для вычислений вектор W
void price_modeling(arma::mat& price_matrix, vector<double>& W)
{
	for (int i = 0; i < M; i++)
	{
		price_matrix(i, 0) = S0;

		for (int j = 1; j < N + 1; j++)
		{
			price_matrix(i, j) = price_matrix(i, j - 1) * exp(pf1 + pf2 * norm_dis(gen));
		}

		W[i] = pay_func(price_matrix(i, N), N);
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

vector<double> least_square_reg(arma::mat prices, vector<double> W, int time)
{
	arma::mat X = construct_at_time(prices, time);
	arma::vec Y = arma::vec(W);
	arma::mat Xt = X.t();
	arma::vec result = arma::pinv(Xt*X)*Xt*Y;
	return arma::conv_to<stdvec>::from(result);
}

//шаг стохастического градиентного шага
vector<double> SGD_step(vector<double> beta, double X_i, double W, double step)
{
	vector<double> x = basis(X_i);
	vector<double> res(basis_fun_am);
	double adjust = 0;
	for (int i = 0; i < basis_fun_am; i++)
	{
		adjust += beta[i] * x[i];
	}

	adjust -= W;
	adjust = adjust * 2 * step;

	for (int i = 0; i < basis_fun_am; i++)
	{
		res[i] = beta[i] - adjust * x[i];
	}
	
	return res;
}

//SGD шаг с регуляризацией
vector<double> SGD_step_L2_reg(vector<double> beta, double X_i, double W, double step, double lambda)
{
	vector<double> x = basis(X_i);
	vector<double> res(basis_fun_am);
	double adjust = 0;
	for (int i = 0; i < basis_fun_am; i++)
	{
		adjust += beta[i] * x[i];
	}

	adjust -= W;
	adjust = adjust * 2 * step;

	for (int i = 0; i < basis_fun_am; i++)
	{
		res[i] = beta[i] * (1 - step * lambda) - adjust * x[i];
	}

	return res;
}

void SGD(vector<double>& beta, vector<double> W, double** prices, int time)
{
	vector<double> next_beta(basis_fun_am);
	uniform_real_distribution<double> first_adj(-1. / (2 * basis_fun_am), 1. / (2 * basis_fun_am));
	for (int i = 0; i < next_beta.size(); i++)
	{
		beta[i] = first_adj(gen);
	}
	int iter_num = 0;
	double w_distance = 1;
	do {
		iter_num++;
		int ind = unif_dis_int(gen);
		//next_beta = SGD_step(beta, prices[ind][time], W[ind], 1. / M);
		next_beta = SGD_step_L2_reg(beta, prices[ind][time], W[ind], 1. / M, 55);
		w_distance = norm_2(beta, next_beta);
		beta = next_beta;
	} while (w_distance > eps);

	for (int i = 0; i < beta.size(); i++)
	{
		cout << beta[i] << "    ";
	}
	cout <<endl <<  " iter am   " << iter_num << endl;
}

//Функция вычисляющая вектор кэфов бета 
void linear_reg(vector<double>& beta, vector<double> W, arma::mat& prices, int time)
{
	beta = least_square_reg(prices, W, time);
	//cout << MSE(prices, W, beta, time) << endl;
}

//Функция выдающая оценку
void opt_price(double& x)
{
	/*
	double** price = new double*[M];
	for (int n = 0; n <= M; n++)
		price[n] = new double[N + 1];
	*/
	arma::mat price = arma::mat(M, N + 1);
	vector<double> W(M);
	vector<double> beta_j(basis_fun_am);

	double option_price = 0;

	price_modeling(price, W); //это without квази числами
	
	for (int time = N; time >= 0; time--)
	{
		linear_reg(beta_j, W, price, N);
		for (int traectory = 0; traectory < M; traectory++)
		{
			double q = 0;
			for (int k = 0; k < basis_fun_am; k++)
			{
				q += beta_j[k] * basis_func(price(traectory, time), k);
			}
			W[traectory] = max(pay_func(price(traectory, time), time), q);
			if (time == 0)
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
	time(&start);
	for (int i = 0; i < numberOfRuns; i++)
	{
		
		double xi = 0;
		opt_price(xi);
		sqSum += xi * xi;
		mean += xi;
		
		//allRunTime += difftime(start, end);
		//fout << endl << "new run " << endl;
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

	//opt_price(var);
	stats(1, sigma, var, mean);
	cout << "Mean:   " << mean << "   Sigma:   " << sigma << "    Variance:     " << var;
	//Сетка при 1000 шагов  2.2048
	//увеличить до 40 мб
	return 0;
}


/*
	for (int i = 0; i < price.n_rows; i++) //вывод траекторий в файл
	{
		for (int j = 0; j < price.n_cols; j++)
		{
			fout << price(i, j) << "  ";
		}
		fout << endl;
	}
*/