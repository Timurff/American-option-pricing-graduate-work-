#include "stdafx.h"
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

normal_distribution<double> norm_dis(0, 1);
uniform_real_distribution<double> unif_dis(0, 1);

random_device rd;
mt19937 gen(rd());




const double PI = 3.14159265358979323846;
const int N = 16, M = 1000, NTHREADS = 2, MSG_DONE = WM_USER + 1;

double a = 0.08, S0 = 100, x0 = log(S0), K = 101, T = 1, r = 0.05, sigma = 0.08; //a - Волатильность,S0 - начальная цена, K - strike price, r - текущая процентная ставка
double dt = T / N,
sqrt_dt = sqrt(dt),
mu_dt = (r - 0.5 * a * a) * dt,
alpha = exp(-r * dt),
sqrt2pi = sqrt(2 * PI);

struct Node
{
public:
	double x, q, Y;
	double* p;
	Node() { p = new double[M]; }
	~Node() { delete[] p; }
};

//Node* node[N + 1];

class MeshThread : public CWinThread
{
	int Run();
	double margin_density();

public:
	int amount_of_runs;
	double S, S2;

	Node* node[N + 1];
	int thread_number;
	bool ready;
	CEvent go;

	MeshThread() : ready(false) { m_bAutoDelete = false; }
	BOOL InitInstance() { return TRUE; }
	int ExitInstance() { return 0; }
	void WaitForAll();
	void OnMsg(WPARAM wParam, LPARAM lParam);
};

MeshThread thread[NTHREADS];

CCriticalSection cs;

bool Set(int nom)
{
	cs.Lock();
	thread[nom].ready = true;
	for (int i = 1; i < NTHREADS; i++)
		if (thread[i].ready == false)
		{
			cs.Unlock();
			return false;
		}
	for (int i = 1; i < NTHREADS; i++)
	{
		thread[i].ready = false;
		thread[i].go.PulseEvent();
	}
	cs.Unlock();
	return true;
}

void MeshThread::WaitForAll()
{
	if (thread_number == 0)
	{
		for (;;)
		{
			MSG msg;
			GetMessage(&msg, 0, 0, 0);
			if (Set(msg.wParam)) break;
		}
	}
	else
	{
		thread[0].PostThreadMessage(MSG_DONE, thread_number, 0);
		WaitForSingleObject(go.m_hObject, INFINITE);
	}
}


///////////////////////////////////////////////////////////////////


double normal_density(double x, double s)
{
	return exp(-x * x / (2 * s*s)) / (sqrt2pi*s);
}

double P(double x, double y)
{
	return normal_density(x + mu_dt - y, a * sqrt_dt);
}

int MeshThread::Run()
{
	S = 0;
	S2 = 0;
	
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
	for (int n = 0; n <= N; n++) 
	{ 
		node[n] = new Node[M]; 
	}
	node[0][0].x = x0;

	for (int n = 0; n < N; n++)
	{
		int n1 = n + 1;
		double nmu = n1 * mu_dt,
			nsig = n1 * sqrt_dt * sigma; //Зачем здесь корень из n1? Если должно быть просто n1?
		for (int m1 = 0; m1 < M; m1++)
		{
			double x = x0 + nmu + nsig * norm_dis(rd);
			node[n1][m1].x = x;
			node[n1][m1].q = normal_density(x0 + nmu - x, nsig);
		}
	}

	for (int n = 0; n < N; n++)
	{
		int n1 = n + 1;
		for (int m = 0; m < M; m++)
		{
			for (int m1 = 0; m1 < M; m1++)
				node[n][m].p[m1] = P(node[n][m].x, node[n1][m1].x);
		}
	}

	for (int m = 0; m < M; m++)
	{
		node[N][m].Y = max(K - exp(node[N][m].x), 0.);
	}

	for (int n = N - 1; n >= 0; n--)
	{
		int n1 = n + 1;

		if (n == 0)
		{
			double SY = 0, S1 = 0;
			for (int m1 = 0; m1 < M; m1++)
			{
				double rho = node[n][0].p[m1] / node[n1][m1].q;
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
				double rho = node[n][m].p[m1] / node[n1][m1].q;
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
	stats(2, sigma, var, mean);
	std::cout << "Mean:   " << mean << "   Sigma:   " << sigma << "    Variance:     " << var;
	return 0;
}