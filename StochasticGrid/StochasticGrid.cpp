#include "stdafx.h"
#include <conio.h>
#include <afxwin.h>
#include <afxmt.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <random>

const double PI = 3.14159265358979323846;
const int N = 16, M = 1000, NTHREADS = 4, MSG_DONE = WM_USER + 1;

double a = 0.08, S0 = 100, x0 = log(S0), K = 101, T = 1, r = 0.05, sigma = 0.08;
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

Node* node[N + 1];

class MeshThread : public CWinThread
{
	int Run();
public:
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

int _tmain(int argc, _TCHAR* argv[])
{
	for (int n = 0; n <= N; n++) node[n] = new Node[M];
	node[0][0].x = x0;


	HANDLE handles[NTHREADS];
	for (int i = 0; i < NTHREADS; i++)
	{
		thread[i].thread_number = i;
		thread[i].CreateThread();
		handles[i] = thread[i].m_hThread;
	}
	WaitForMultipleObjects(NTHREADS, handles, true, INFINITE);

	printf("%f", node[0][0].Y);

	for (int n = 0; n <= N; n++)  delete[] node[n];
	_getch();
	return 0;
}

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
	std::default_random_engine rnd(2 * thread_number + 1);
	std::normal_distribution<double> normal(0, 1);
	std::uniform_real_distribution<double> unif(0, 1);
	//if(nom==0){ MSG msg; PeekMessage(&msg, NULL, WM_USER, WM_USER, PM_NOREMOVE);}

	int Nsteps = N / NTHREADS, bs = thread_number * Nsteps, es = bs + Nsteps;
	for (int n = bs; n < es; n++)
	{
		int n1 = n + 1;
		double nmu = n1 * mu_dt,
			nsig = n1 * sqrt_dt * sigma; //Зачем здесь корень из n1? Если должно быть просто n1?
		for (int m1 = 0; m1 < M; m1++)
		{
			double x = x0 + nmu + nsig * normal(rnd); //Умножение?
			node[n1][m1].x = x;
			node[n1][m1].q = normal_density(x0 + nmu - x, nsig);
		}
	}
	WaitForAll();

	for (int n = bs; n < es; n++)
	{
		int n1 = n + 1;
		for (int m = 0; m < M; m++)
		{
			for (int m1 = 0; m1 < M; m1++)
				node[n][m].p[m1] = P(node[n][m].x, node[n1][m1].x);
		}
	}
	WaitForAll();

	int Nnodes = M / NTHREADS, bn = thread_number * Nnodes, en = bn + Nnodes;
	for (int m = bn; m < en; m++)
		node[N][m].Y = max(K - exp(node[N][m].x), 0.);
	WaitForAll();

	for (int n = N - 1; n >= 0; n--)
	{
		if (n == 0)
			if (thread_number == 0) en = 1; else break;
		int n1 = n + 1;
		for (int m = bn; m < en; m++)
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
		if (n != 0) WaitForAll();
	}

	return 0;
}

