#include "pch.h"
#include "sobol.h"
#include <math.h>
#include <sstream>
#include <iostream>


using namespace std;

const unsigned int nbits=8*sizeof(int);
const double pbits = pow(2., -(double)nbits);

Sob1::Sob1()
{ 
	M.resize(nbits);
}


void Sob1::init(istream& in, unsigned int n)
{
	unsigned int a1, n1;
	while (!in.eof())
	{
		char line[128];
		in.getline(line, 128);
		stringstream sin(line);
		sin >> n1 >> s >> a1;
		if (n == n1)
		{
			a.resize(s);
			a[s - 1] = 1;
			for (int k = s - 2, b = 1; k >= 0; k--, b <<= 1) a[k] = b & a1 ? 1 : 0;
			for (unsigned int k = 0, m = 0; k < s; k++)
			{
				sin >> m;
				M[k] = m << (nbits - 1 - k);
			}
			break;
		}
	}

	for (unsigned int k = s; k < nbits; k++)
	{
		unsigned int Mk = M[k - s] >> s;
		for (unsigned int t = 0; t < s; t++)
			if (a[t]) Mk ^= M[k - t - 1];
		M[k] = Mk;
	}
	mm = 0;
	ii = 0;
}

void Sob1::next(double& x)
{
	unsigned int l = 1, j = ii - 1;
	for (unsigned int b = 1; l < nbits; l++, b <<= 1)
		if ((j & b) == 0) break;
	unsigned int mm1 = mm ^ M[l - 1];
	x = mm1 * pbits;
	mm = mm1;
	ii++;
}


Sobol::Sobol(vector<int>& nn)
{ 
	ifstream in("joe-kuo.txt");
	char line[64];
	in.getline(line, 64);
	sob1.resize(nn.size());
	for (unsigned int i = 0; i < nn.size(); i++) sob1[i].init(in, nn[i]);
}

Sobol::Sobol(int n)
{ 
	ifstream in("joe-kuo.txt");
	char line[64];
	in.getline(line, 64);
	sob1.resize(n);
	for (size_t i = 0; i < n; i++) sob1[i].init(in, i + 2);
}

void Sobol::next(double *x)
{ 
	int n=sob1.size();
	for (int i = 0; i < n; i++)  sob1[i].next(x[i]);
}


