#include <cassert>


#include <cstdlib>
using std::exit;

#include <cmath>
//#define M_PI 3.14159265358979323846
using std::sqrt;
using std::fabs;
using std::fmod;
using std::pow;

#include <limits>
using std::numeric_limits;

#include <iomanip>
using std::setprecision;

#include <string>
using std::string;

#include <complex>
using std::complex;

#include <vector>
using std::vector;

#include <algorithm>
using std::max;
using std::min;
using std::max_element;
using std::min_element;
using std::swap;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <fstream>
using std::ofstream;

#include <sstream>
using std::stringstream;

using namespace std;

int main()
{
	const int n = 200;
	const int i = 300;

	int xnot = 20;

	double sigma = 4;
	double L = i/100;
	double h = L/(i-1);
	double tau = 0.00005;
	double mu = 0;

	cout << "Input Potential Barrier Height (0.0-1.0):" << endl;
	double V;
	cin >> V;
	double a;
	cout<< "Input Potential Barrier Width (1-10):" << endl;
	cin >> a;
	double E;
	cout<< "Input Wave Packet Energy (0.1-1.0):" << endl;
	cin >> E;

	double m = 1;
	double hbar = 1;
	double psi11R[i];
	double psi12R[i];
	double psi11L[i];
	double psi12L[i];

	double k = sqrt(2*m*E/(hbar*hbar));

	double T = 1;
	if (V!= 0){
		if(E<V){
			mu = sqrt(2*m*(V-E)/(hbar*hbar));
			T = 1/((k*k+mu*mu)*(k*k+mu*mu)*sinh(2*mu*a)*sinh(2*mu*a)/(4*mu*mu*k*k)+1);
		}
		if(E>V){
			mu = sqrt(2*(E-V));
			T = 1/((k*k-mu*mu)*(k*k-mu*mu)*sin(2*mu*a)*sin(2*mu*a)/(4*mu*mu*k*k)+1);
		}
		if(E==V)
			T = 1/(1+4*a*a*V/2);
	}
	double R=1-T;

	ofstream R_psi_out;
	ofstream L_psi_out;
	ofstream psi_out;

	//initialize psi. Scaled to unity at xnot
	for (int j = 0; j<i; j++){
		psi11R[j]= exp(-pow(j-xnot,2)/(2*(sigma*sigma)));
	}

	//initialize moive file
	ofstream fgp("movieP.gp");
	fgp << "set title \"Paritcle Incidnet on Potential Barrier\"" << endl;
	fgp << "set xlabel \"Position\"" << endl;
	fgp << "set ylabel \"Potential Energy\"" << endl;
	fgp << "set yrange [0:1]" << endl;
	fgp << "set xrange [20:240]" << endl;
	fgp << "pause -1 \"Press return to start\"" << endl;

	//create Potential on plot
	psi_out.open("potential.dat");
	for (int j=0;j<i;j++){
		if (j<=i/2-30 || j >= i/2-30+a+1){
			psi_out<< j << " " << 0 << endl;
		}
		if (j>i/2-30 && j<i/2-30+a+1){
			psi_out<< j << " " << V << endl;
		}
	}
	psi_out.close();

	//fill psi
	for (int t=0; t<n; t++){
		ostringstream my_Rstringstream;
		ostringstream my_Lstringstream;

		my_Rstringstream << "R_time_step_" << t  <<".dat";
		my_Lstringstream << "L_time_step_" << t  <<".dat";

		R_psi_out.open(my_Rstringstream.str().c_str());
		L_psi_out.open(my_Lstringstream.str().c_str());

	if(t<n/2){
			for (int j=1;j<i-1;j++){
				psi12R[j+1]=psi11R[j]+tau/(2*h*h)*(psi11R[j+1]+psi11R[j-1]-psi11R[j]-psi11R[j]);

	 			R_psi_out<< scientific  << j << "	" << psi12R[j]*psi12R[j]<< endl;
			}

			for (int j=1;j<i;j++){
				psi11R[j]=psi12R[j];
			}

			fgp << "plot \"Potential.dat\" using 1:2 with filledcurves title \"Potential = "<< V <<" Width = "<< a <<"\", \
			 \"R_time_step_" << t <<".dat\" with filledcurves title \"Incident Particle = " << E << "\"" << endl;
			fgp << "pause 0.05" << endl;
		}

		if(t==n/2){
			for (int j=1;j<i;j++){
				psi11L[j]=R*psi12R[j];
				psi11R[j]=T*psi12R[j];
			}
		}

		if(t>=n/2){

			for (int j=1;j<i-1;j++){
				psi12R[j+1]=psi11R[j]+tau/(2*h*h)*(psi11R[j+1]+psi11R[j-1]-psi11R[j]-psi11R[j]);
					R_psi_out<< scientific  << j << "	" << psi12R[j]*psi12R[j]<< endl;

			}

			for (int j=i-1;j>1;j--){
				psi12L[j-1]=psi11L[j]+tau/(2*h*h)*(psi11L[j+1]+psi11L[j-1]-psi11L[j]-psi11L[j]);

					L_psi_out<< scientific  << j << "	" << psi12L[j]*psi12L[j]<< endl;

			}

			for (int j=1;j<i;j++){
				psi11R[j]=psi12R[j];
			}

			for (int j=1;j<i;j++){
				psi11L[j]=psi12L[j];
			}
			fgp << "plot \"Potential.dat\" using 1:2 with filledcurves title \"Potential = "<< V <<" Width = "<< a <<"\", \
			 \"R_time_step_" << t <<".dat\" with filledcurves title \"Transmitted Particle T = "<< setprecision(2) << T <<"\", \
			 \"L_time_step_" << t <<".dat\" with filledcurves title \"Reflected Particle R = "<< setprecision(2) << R <<"\""<< endl;
			fgp << "pause 0.05" << endl;
		}

		R_psi_out.close();
		L_psi_out.close();
			}


	return 0;
}
