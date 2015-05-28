#include <string>
#include <fstream>
#include <complex>

const int Nr = 1000;
double grid_size = 15.;
double rad = 1.;

double *r_coord;
complex<double> *field;
double *TDMA_a, *TDMA_b, *TDMA_c;
complex<double> *TDMA_u, *TDMA_v;

template <typename T>
T sqr(T x)
{
	return x * x;
}

void SetParam()
{
}

void OutParam(const std::string& file_name)
{
	std::ofstream out(file_name.c_str());
	out << "Nr\t= " << Nr << "\n";
	out << "grid_size\t= " << grid_size << " cm" << "\n";
	out << "rad\t= " << rad << " cm" << "\n";
}

void Allocate()
{
	r_coord = new double [Nr + 1];
	TDMA_a = new double [Nr];
	TDMA_b = new double [Nr];
	TDMA_c = new double [Nr];
	TDMA_u = new complex<double> [Nr];
	TDMA_v = new complex<double> [Nr];
	field = new complex<double> [Nr];
}

void InitRCoord()
{
	for (int i = 0; i <= Nr; ++i)
		r_coord[i] = (i + 0.5) * grid_size / (Nr + 0.5);
}

void InitTDMACoeff()
{
	TDMA_a[0] = (3.0 - r_coord[1] / r_coord[0]) / ((r_coord[1] + r_coord[0]) * (r_coord[0] + r_coord[0]));
	TDMA_b[0] = (3.0 + r_coord[0] / r_coord[0]) / ((r_coord[1] - r_coord[0]) * (r_coord[1] + r_coord[0]));
	TDMA_c[0] = TDMA_a[0] + TDMA_b[0];
	for (int i = 1; i < Nr; ++i)
	{
		TDMA_a[i] = (3.0 - r_coord[i + 1] / r_coord[i]) / ((r_coord[i + 1] - r_coord[i - 1]) * (r_coord[i] - r_coord[i - 1]));
		TDMA_b[i] = (3.0 - r_coord[i - 1] / r_coord[i]) / ((r_coord[i + 1] - r_coord[i - 1]) * (r_coord[i + 1] - r_coord[i]));
		TDMA_c[i] = TDMA_a[i] + TDMA_b[i];
	}
}

void SetFieldGauss()
{
	for (int i = 0; i < Nr; ++i)
		field[i] = exp(-0.5 * sqr(r_coord[i] / grid_size));
}

void InitField()
{
	SetFieldGauss();
}

void Initialize()
{
	InitRCoord();
	InitTDMACoeff();
	InitField();
}

void Linear()
{
	
}

void NonLinear()
{
}

void Advance()
{
	Linear();
	NonLinear();
}

void ChooseStep()
{
}

void Save()
{
}

void Propagate()
{
	for (int iter_num = 0; iter_num < 10; ++iter_num)
	{
		Advance();
		ChooseStep();
		Save();
	}
}

void Release()
{
	delete [] r_coord;
	delete [] TDMA_a;
	delete [] TDMA_b;
	delete [] TDMA_c;
	delete [] TDMA_u;
	delete [] TDMA_v;
}

int main()
{
	SetParam();
	OutParam();
	Allocate();
	Initialize();
	Propagate();
	Release();
}