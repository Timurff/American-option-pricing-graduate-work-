#ifndef SOBOL
#define SOBOL

#include <vector>
#include <fstream>


using namespace std;

struct Sob1
{ 
	unsigned int s, mm, ii;
	vector<unsigned int> a, M;
	Sob1();
	void init(istream&, unsigned int n);
	void next(double& x);
};


class Sobol
{   
	vector<Sob1> sob1;
	public:
		Sobol(vector<int>& nn);
		Sobol(int n);
		void next (double*);
};

#endif // SOBOL

