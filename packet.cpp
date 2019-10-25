#include <cassert>

#include <cstddef>
using std::size_t;

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


//#ifdef __APPLE__
//#include <Accelerate/Accelerate.h>
//#include <vecLib/clapack.h>
//typedef __CLPK_doublereal doublereal;
//typedef __CLPK_integer integer;
//#else
typedef double doublereal;
typedef int integer;
extern "C" void dspevd_(char*, char*, integer*, doublereal*, doublereal*,
						doublereal*, integer*, doublereal*, integer*,
						integer*, integer*, integer*);
//#endif

size_t Lx, Ly, N;
double kx, ky, sigma0, sigma0_2, x_init, y_init, T, dt;

void build_Hamiltonian(vector<doublereal> &H)
{
	for (size_t n = 0; n < H.size(); ++n)
		H[n] = 0.0;

	for (size_t y = 0; y < Ly; ++y)
		for (size_t x = 0; x < Lx; ++x)
		{
			const size_t i = x+Lx*y;
			if (x != Lx-1)
			{
				const size_t j = (x+1)+Lx*y;
				H[i+j*(j+1)/2] = -0.5;
			}
			if (y != Ly-1)
			{
				const size_t j = x+Lx*(y+1);
				H[i+j*(j+1)/2] = -0.5;
			}
		}
}

void eigensolve(vector<doublereal> &H, vector<doublereal> &Eval, vector<doublereal> &Evec)
{
	// Solve the eigenvalue problem with LAPACK's dsepvd routine

	assert(H.size() == N*(N+1)/2);
	assert(Eval.size() == N);
	assert(Evec.size() == N*N);

	integer info;
	char jobz='V';
	char uplo='U';
	integer M = N;
	vector<doublereal> work(1+6*M+M*M);
	integer lwork = work.size();
	vector<integer> iwork(3+5*M);
	integer liwork = iwork.size();

	dspevd_(&jobz,&uplo,&M,&(H[0]),&(Eval[0]),&(Evec[0]),&M,&(work[0]),&lwork,&(iwork[0]),&liwork,&info);

	assert(info == 0);

	// check that eigenvectors are orthonormal
	for (size_t n = 0; n < N; ++n)
	{
		double sum = 0.0;
		for (size_t i = 0; i < N; ++i)
			sum += Evec[n+i*N]*Evec[n+i*N];
		assert(fabs(sum-1.0) < sqrt(numeric_limits<double>::epsilon()));
	 }

	for (size_t i = 0; i < N; ++i)
	{
		double sum = 0.0;
		for (size_t n = 0; n < N; ++n)
			sum += Evec[n+i*N]*Evec[n+i*N];
		assert(fabs(sum-1.0) < sqrt(numeric_limits<double>::epsilon()));
	}
}

void compute_dos(vector<doublereal> &Eval)
{
	const double Emin = *min_element(Eval.begin(),Eval.end());
	const double Emax = *max_element(Eval.begin(),Eval.end());
	const double Elow = Emin - 0.05*(Emax-Emin);
	const double Ehigh = Emax + 0.05*(Emax-Emin);
	const double width = Ehigh-Elow;
	vector<size_t> dos(100);
	for (size_t n = 0; n < 100; ++n)
		dos[n] = 0;
	for (size_t n = 0; n < Eval.size(); ++n)
	{
		const size_t p = size_t(100*(Eval[n]-Elow)/width);
		assert(p < dos.size());
		++dos[p];
	}

	ofstream fout("dos.dat");

	for (size_t n = 0; n < 100; ++n)
		fout << Elow+width*(n+0.5)/100 << "\t" << dos[n]*100/width/Eval.size() << endl;

	fout.close();
}

void build_packet(vector<doublereal> &Eval, vector<doublereal> &Evec,
				  vector< complex<double> > &psi, vector< complex<double> > &c,
				  double rx0, double ry0, double kx0, double ky0)
{
	double normalization = 0.0;
	for (size_t alpha = 0; alpha < N; ++alpha)
	{
		c[alpha] = 0.0;
		for (size_t y = 0; y < Ly; ++y)
			for (size_t x = 0; x < Lx; ++x)
			{
				const size_t i = x+Lx*y;
				const double phi = kx0*x + ky0*y;
				const double dx = double(x)-rx0;
				const double dy = double(y)-ry0;
				c[alpha] += Evec[i+N*alpha]*exp(-0.25*(dx*dx+dy*dy)/sigma0_2)*exp(complex<double>(0,phi));
			}
		normalization += norm(c[alpha]);
	}

	const double rescale = 1.0/sqrt(normalization);
	for (size_t alpha = 0; alpha < N; ++alpha)
		c[alpha] *= rescale;
}

void write(vector< complex<double> > &psi, ofstream &fgp, double t)
{
	static unsigned int dat_count = 0;
	stringstream ss("psi-");
	ss.fill('0');
	ss.width(5);
	ss << dat_count++ << ".dat";
	ofstream fout(ss.str().c_str());
	if (Ly == 1)
	{
		if (dat_count == 1)
		{
			fgp << "set samples 1000" << endl;
			fgp << "x0 = " << x_init << endl;
			fgp << "k0 = " << kx << endl;
			fgp << "set yrange [] writeback" << endl;
		}
		fgp << "t = " << t << endl;
		fgp << "sig0 = " << sigma0 << endl;
		fgp << "sigt2 = " << sigma0_2 + 0.25*t*t/sigma0_2 << endl;
		fgp << "plot \"" << ss.str() << "\" with lines, exp(-0.5*(x-x0-k0*t)**2/sigt2)/sqrt(2*pi*sigt2)" << endl;
		if (dat_count == 1)
			fgp << "pause -1 \"Press return to start\"" << endl;
		else
		{
			fgp << "set yrange restore" << endl;
			fgp << "pause 0.05" << endl;
		}

		for (size_t x = 0; x < Lx; ++x)
			fout << x << " " << norm(psi[x]) << endl;
	}
	else
	{
		if (dat_count == 1)
			fgp << "set zrange [] writeback" << endl;
		fgp << "splot \"" << ss.str() << "\" matrix with lines" << endl;
		if (dat_count == 1)
			fgp << "pause -1 \"Press return to start\"" << endl;
		else
		{
			fgp << "set zrange restore" << endl;
			fgp << "pause 0.05" << endl;
		}

		for (size_t y = 0; y < Ly; ++y)
		{
			for (size_t x = 0; x < Lx; ++x)
				fout << norm(psi[x+Lx*y]) << " ";
			fout << endl;
		}
	}
	fout.close();
}

double build_wavefunction(vector<doublereal> &Eval, vector<doublereal> &Evec,
						  vector< complex<double> > &psi, vector< complex<double> > &c, double t)
{
	double normalization = 0.0;
	for (size_t i = 0; i < N; ++i)
	{
		psi[i] = 0.0;
		for (size_t alpha = 0; alpha < N; ++alpha)
		{
			double theta = -Eval[alpha]*t;
			psi[i] += Evec[i+N*alpha]*c[alpha]*exp(complex<double>(0.0,theta));
		}
		normalization += norm(psi[i]);
	}
	return normalization;
}

void ave_position(const vector< complex<double> > &psi, double &ave_x, double &ave_y)
{
	ave_x = ave_y = 0.0;
	for (size_t y = 0; y < Ly; ++y)
		for (size_t x = 0; x < Lx; ++x)
		{
			const size_t i = x+Lx*y;
			const double psi2 = norm(psi[i]);
			ave_x += x*psi2;
			ave_y += y*psi2;
		}
}

void bad()
{
	cerr << "Usage:" << endl;
	cerr << "  packet -L=#,# dos" << endl;
	cerr << "  packet -L=#,# evolution -k=#,# -w=#" << endl;
	cerr << "  packet -L=#,# trajectory -w=#" << endl;
	exit(1);
}

string parse_command_line(int argc, char* argv[])
{
	if (argc < 3) bad();

	{ // parse the lattice size flag
		const string Lflag = argv[1];
		if (Lflag.substr(0,3) != string("-L=")) bad();
		string Lvalues = Lflag.substr(3);
		size_t pos = 0;
		while ((pos = Lvalues.find(',',pos)) != string::npos) Lvalues[pos] = ' ';
		stringstream ss(Lvalues);
		ss >> Lx;
		if (ss.fail()) bad();
		ss >> Ly;
		if (ss.fail()) Ly = 1;
		if (Lx == 0 or Ly == 0) bad();
		if (Ly > Lx) swap(Lx,Ly);
		N = Lx*Ly;
		x_init = 0.25*Lx;
		y_init = 0.5*Ly;
	}
	const string mode = argv[2];
	if (mode == "evolution")
	{
		if (argc != 5) bad();

		{ // parse the wavevector flag
			const string kflag = argv[3];
			if (kflag.substr(0,3) != string("-k=")) bad();
			string kvalues = kflag.substr(3);
			size_t pos = 0;
			while ((pos = kvalues.find(',',pos)) != string::npos) kvalues[pos] = ' ';
			stringstream ss(kvalues);
			ss >> kx;
			if (ss.fail()) bad();
			ss >> ky;
			if (ss.fail()) ky = 0.0;
		}

		{ // parse the width flag
			const string wflag = argv[4];
			if (wflag.substr(0,3) != string("-w=")) bad();
			stringstream ss(wflag.substr(3));
			ss >> sigma0;
		}

		T = 3.0*(Lx-x_init)/max(fabs(kx),0.1);
		dt = T/500.0;
	}
	else if (mode == "trajectory")
	{ // parse the width flag
		if (argc != 4) bad();

		const string wflag = argv[3];
		if (wflag.substr(0,3) != string("-w=")) bad();
		stringstream ss(wflag.substr(3));
		ss >> sigma0;

		if (ss.good()) bad();
	}
	else if (mode != "dos") bad();

	sigma0_2 = sigma0*sigma0;

	if (Ly == 1)
		cout << "A linear mesh of " << N << " points" << endl;
	else
		cout << "A rectangular " << Lx << "x" << Ly << " mesh of " << N << " points" << endl;
	if (mode == "evolution")
	{
		cout << "Initial wavepacket has wavevector k=";
		if (Ly == 1)
			cout << "(" << kx << "," << ky << ")";
		else
			cout << kx;
		cout << " and (half-max) width " << sigma0 << endl;
	}

	return mode;
}
int main(int argc, char* argv[])
{
	const string mode = parse_command_line(argc,argv);

	vector<doublereal> H(N*(N+1)/2);
	vector<doublereal> Eval(N);
	vector<doublereal> Evec(N*N);

	build_Hamiltonian(H);
	cout << "Diagonalizing the Hamiltonian ..." << endl;
	eigensolve(H,Eval,Evec);

	if (mode == "dos") compute_dos(Eval);
	else
	{
		vector< complex<double> > psi(N);
		vector< complex<double> > c(N);

		if (mode == "evolution")
		{
			build_packet(Eval,Evec,psi,c,x_init,y_init,kx,ky);

			ofstream fgp("movie.gp");

			for (double t = 0; t < T; t += dt)
			{
				cout << "t=" << t << ": " << build_wavefunction(Eval,Evec,psi,c,t) << endl;
				write(psi,fgp,t);
			}

			fgp.close();
		}
		else if (mode == "trajectory")
		{
			T = 0.25*(Lx-x_init)/M_PI;
			dt = T/20.0;
			ofstream fdat("traj.dat");
			for (int n = 0; n <= 12; ++n)
			{
				kx = M_PI*n/12.0;
				if (Ly == 1)
					cout << "Computing k="<< kx << endl;
				else
					cout << "Computing k=("<< kx << ",0)" << endl;
				build_packet(Eval,Evec,psi,c,x_init,y_init,kx,0.0);
				for (double t = 0; t < T; t += dt)
				{
					build_wavefunction(Eval,Evec,psi,c,t);
					double ave_x, ave_y;
					ave_position(psi,ave_x,ave_y);
					fdat << kx << "\t" << t << "\t" << ave_x-x_init << "\t" << ave_y-y_init << endl;
				}
				fdat << endl << endl;
			}
			fdat.close();
		}
	}

	return 0;
}
