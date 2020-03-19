#include "stdafx.h"
#include "sobol.h"
#define _USE_MATH_DEFINES
#include <conio.h>
#include <afxwin.h>
#include <afxmt.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <ctime>
#include <time.h>
#include <algorithm>
#include <random>

using namespace std;

const double PI = 3.14159265358979323846;
const int N = 16, M = 2000, NTHREADS = 1, MSG_DONE = WM_USER + 1, quasi_points_am = 100, random_points_am = M / quasi_points_am;

Sobol sobol(2);

double a = 0.08, S0 = 100, x0 = log(S0), K = 101, T = 1, r = 0.05, sigma = 0.08; //a - ¬олатильность,S0 - начальна€ цена, K - strike price, r - текуща€ процентна€ ставка
double dt = T / N,
sqrt_dt = sqrt(dt),
mu_dt = (r - 0.5 * a * a) * dt,
alpha = exp(-r * dt),
sqrt2pi = sqrt(2 * PI);

struct Node
{
public:
	double x, q, Y;
	Node() { x = 0; q = 1; Y = 0;}
	~Node() {}
};

class MeshThread : public CWinThread
{
	int Run();
	void mesh_init();
	void mesh_fill();
	double margin_density();
	void make_quasi_points();
	void make_random_points();
	double P(double x, double y);
	double normal_density(double x, double s);
	void rnorm_with_quasi(int p, int q, double& normal_value);

public:

	int amount_of_runs;
	double S, S2;

	double quasi_points[quasi_points_am][2];
	double random_points[random_points_am][2];

	vector<vector<Node>> node;
	int thread_number;
	bool ready;
	CEvent go;

	MeshThread() : ready(false) { m_bAutoDelete = false; }
	BOOL InitInstance() { return TRUE; }
	int ExitInstance() { return 0; }
};

MeshThread thread[NTHREADS];

CCriticalSection cs;


///////////////////////////////////////////////////////////////////

void MeshThread::mesh_init()
{
	node.resize(N + 1);

	for (int i = 0; i < N + 1; i++)
	{
		node[i].resize(M);
	}
}

void MeshThread::mesh_fill()
{
	node[0][0].x = x0;
	
	make_quasi_points();

	for (int n = 0; n < N; n++)
	{
		make_random_points();
		int n1 = n + 1;
		double nmu = n1 * mu_dt,
			nsig = n1 * sqrt_dt * sigma; //«ачем здесь корень из n1? ≈сли должно быть просто n1?
		//делаю дл€ каждого  n набор точек просто рандомных
		for (int m1 = 0; m1 < M; m1++)
		{
			int p = m1 % random_points_am,
				q = m1 / random_points_am;
			double norm = 0;
			rnorm_with_quasi(p, q, norm);
			double x = x0 + nmu + nsig * norm;
			node[n1][m1].x = x;
			node[n1][m1].q = normal_density(x0 + nmu - x, nsig);
		}
	}

	for (int m = 0; m < M; m++)
	{
		node[N][m].Y = max(K - exp(node[N][m].x), 0.);
	}

}


double MeshThread::normal_density(double x, double s)
{
	return exp(-x * x / (2 * s*s)) / (sqrt2pi*s);
}

double MeshThread::P(double x, double y)
{
	return normal_density(x + mu_dt - y, a * sqrt_dt);
}

int MeshThread::Run()
{
	S = 0;
	S2 = 0;
	
	mesh_init();

	for (int i = 0; i < amount_of_runs; i++)
	{
		double price = margin_density();
		S += price;
		S2 += price * price;
	}

	return 0;
}

double MeshThread::margin_density()
{
	//сгенерировать квази точки
	//на каждом шаге генерировать обычные точки
	//рандомизировать все это дело 
	//генерировать нормальную случ величину дл€ поиска ’

	//делаю набор квази точек он будет 1 на всю работу алгоритма
	
	mesh_fill();

	for (int n = N - 1; n >= 0; n--)
	{
		int n1 = n + 1;

		if (n == 0)
		{
			double SY = 0, S1 = 0;
			for (int m1 = 0; m1 < M; m1++)
			{	
				double rho = P(node[n][0].x, node[n1][m1].x) / node[n1][m1].q;
				SY += rho * node[n1][m1].Y;
				S1 += rho;
			}
			double f = max(K - exp(node[n][0].x), 0.);
			node[n][0].Y = max(f, alpha*SY / S1);
			break;
		}

		for (int m = 0; m < M; m++)
		{
			double SY = 0, S1 = 0;
			for (int m1 = 0; m1 < M; m1++)
			{
				double rho = P(node[n][m].x, node[n1][m1].x) / node[n1][m1].q;
				SY += rho * node[n1][m1].Y;
				S1 += rho;
			}
			double f = max(K - exp(node[n][m].x), 0.);
			node[n][m].Y = max(f, alpha*SY / S1);
		}
	}

	return node[0][0].Y;
}

void stats(int numberOfRuns, double& sigma, double& var, double& mean)
{
	double sqSum = 0;
	time_t start, end;
	double allRunTime = 0;
	double xi = 0;

	time(&start);

	HANDLE handles[NTHREADS];

	for (int j = 0; j < NTHREADS; j++)
	{
		thread[j].amount_of_runs = numberOfRuns / NTHREADS;
		thread[j].CreateThread();
		handles[j] = thread[j].m_hThread;
	}
	WaitForMultipleObjects(NTHREADS, handles, true, INFINITE);

	double S = 0, S2 = 0, amount_of_runs = 0;

	for (int i = 0; i < NTHREADS; i++)
	{
		S += thread[i].S;
		S2 += thread[i].S2;
		amount_of_runs += thread[i].amount_of_runs;
	}

	S /= amount_of_runs;
	S2 /= amount_of_runs;

	sqSum = S2;
	mean = S;

	time(&end);
	allRunTime += difftime(start, end);
	
	var = sqSum - (mean * mean);
	sigma = sqrt(var);

	std::cout << "Average time  8k, trajectories " << allRunTime / numberOfRuns << std::endl;
	std::cout << "Mean " << mean << " var " << var << " sigma " << sigma << std::endl;

	_getch();
}

int _tmain(int argc, _TCHAR* argv[])
{
	setlocale(LC_ALL, "russian");
	//fout.open("output.txt");
	double sigma = 0;
	double var = 0;
	double mean = 0;
	stats(1, sigma, var, mean);
	std::cout << "Mean:   " << mean << "   Sigma:   " << sigma << "    Variance:     " << var;
	return 0;
}


//хочу функции создани€ квази точек, просто рандомныз точек и функцию делающую из них нормальную величину

void MeshThread::make_quasi_points()
{
	for (int i = 0; i < quasi_points_am; i++)
	{ 
		sobol.next(quasi_points[i]); 
	}
}

void MeshThread::make_random_points()
{
	normal_distribution<double> norm_dis(0, 1);
	uniform_real_distribution<double> unif_dis(0, 1);

	random_device rd;
	mt19937 gen(rd());

	for (int i = 0; i < random_points_am; i++)
	{
		random_points[i][0] = unif_dis(gen);
		random_points[i][1] = unif_dis(gen);
	}
}

void MeshThread::rnorm_with_quasi(int p, int q, double& normal_value)
{
	double first = random_points[p][0] + quasi_points[q][0];
	if (first > 1)
	{
		first -= 1;
	}
	double second = random_points[p][1] + quasi_points[q][1];
	if (second > 1) 
	{
		second -= 1;
	}
	double temp;
	temp = sqrt(-2 * log(first));
	double teta = 2 * M_PI * second;
	normal_value = temp * cos(teta);
}