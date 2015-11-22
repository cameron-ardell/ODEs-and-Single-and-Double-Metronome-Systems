#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"

using namespace std;

/*struct Implicit{
	Doub lambda;
	Doub n;
	Implicit(Doub lambda_in, Doub n_in) {lambda_in = lambda; n = n_in;};
	~Implicit() {cout << "destructing Implicit...\n";};
};

struct Explicit{
	Doub coef;
	Explicit(Doub lambda_in, Doub n_in) {coef = 1.0 - (1.0/n_in) * lambda_in;};
	~Explicit() {cout << "destructing Explicit...\n";};

	Doub nextVal (Doub lastY){
		return lastY*coef;
	}
};*/



int main(){

	Doub n = 51.0;
	Doub lambda = 100.0;
	Doub h = 1.0/n;

	Doub coef = lambda * h;


	Doub curX = 0.0;
	Doub lastYExp = 1.0;
	Doub lastYImp = 1.0;
	Doub analytic = 0.0;

	ofstream outfile;
	outfile.open("data51.dat");
	outfile << setw(16) << "x val " << setw(16) << "Implicit Val ";
	outfile<< setw(16) << "Explicit Val" << setw(16) << "Analytic Val";
	outfile << setw(16) << "Relative Implicit Error " << setw(16)<< "Relative Explicit Error" << endl;
	outfile << "#=========================================" << endl;


	outfile << setw(16) << curX << setw(16) << lastYImp << setw(16) << lastYExp << endl;

	for(Doub step = 1.0; step < n; step++){

		curX = curX + h;

		lastYExp = (1.0 - coef) * lastYExp;
		lastYImp = lastYImp / (1.0 + coef);

		analytic = exp(-lambda * curX);

		Doub relExpError = abs(analytic - lastYExp)/analytic;
		Doub relImpError = abs(analytic - lastYImp)/analytic;


		outfile << setw(16) << curX << setw(16) << lastYImp;
		outfile << setw(16) << lastYExp << setw(16) << analytic;
		outfile << setw(16) << relImpError << setw(16) << relExpError << endl;


	}

	outfile.close();



}