#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>

#define N 7
#define K 11

enum ListOfSubstances {
	StainlessSteel = 0,
	Iron = 1,
	Gold = 2,
	Brass = 3,
	Copper = 4,
	Mercury = 5,
	Silver = 6,
	Diamond=7,
	Glass=8,
	Vacuum=9,
	Ones=10
};

struct MaterialsData {
	const double cv_HeatCapacity;
	const double lambda_ThermalConductivity;
	const double rho_Density;
};

namespace MaterialDataBase {
	std::vector<MaterialsData> materials = { 
		{490, 15, 7800},
		{452, 92, 7870}, 
		{128.7, 320, 19300}, 
		{378, 100, 8600},
		{383, 401, 8900 }, 
		{139, 8.3, 13540}, 
		{235, 430, 10500}, 
		{502, 1700, 3500}, 
		{500, 1, 25}, 
		{0,0,0},
		{1,1,1}
	};
}

struct IntegrateParameters {
	double delta_tStep;
	double h_xStep;
};

class HeatEquationSolver {
public:

	//HeatEquationSolver(); //class constructor
	//~HeatEquationSolver() {}; //destructor

	double InitialValueFunction(double x);
	double LeftBoundaryConditionFunction(double t);
	double RightBoundaryConditionFunction(double t);

	void SetCoordinateXGrid(double* x, const IntegrateParameters& params);
	void SetTimeGrid(double* t, const IntegrateParameters& params);

	void getSolutionExplicit(const MaterialsData& material, double(*T)[N], double* x, double* t, const IntegrateParameters& params);
	void getSolutionImplicit(const MaterialsData& material, double(*T)[N], double* x, double* t, const IntegrateParameters& params);
};


//НУ T(x,0) 
double HeatEquationSolver::InitialValueFunction(double x)
{
	return 0.9 + 2 * x * (1 - x);
}

//ГУ1 T(0, t)
double HeatEquationSolver::LeftBoundaryConditionFunction(double t)
{
	return 3 * (0.3 - 2 * t);
}

//ГУ2 T(l,t)
double HeatEquationSolver::RightBoundaryConditionFunction(double t)
{
	return 1.38;
}

//Создание сетки по координате х
void HeatEquationSolver::SetCoordinateXGrid(double* x, const IntegrateParameters& params)
{
	x[0] = 0;
	for (int i = 1; i < N; ++i) {
		x[i] = i * params.h_xStep;
	}
}

//Создание сетки по времени 
void HeatEquationSolver::SetTimeGrid(double* t, const IntegrateParameters& params)
{
	t[0] = 0;
	for (int k = 1; k < K; ++k) {
		t[k] = k * params.delta_tStep;
	}
}

void HeatEquationSolver::getSolutionExplicit(const MaterialsData& material, double (*T)[N], double* x, double* t, const IntegrateParameters& params)
{
	/*Задали начальные условия*/
	for (int i = 0; i < N; ++i)
	{
		T[0][i] = InitialValueFunction(x[i]);
	}

	/*Задали краевые условия*/
	for (int k = 0; k < K; ++k)
	{
		T[k][0] = LeftBoundaryConditionFunction(t[k]);
		T[k][N] = RightBoundaryConditionFunction(t[k]);
	}

	/*Для сокращения выражений*/
	double c_v = material.cv_HeatCapacity;
	double lambda = material.lambda_ThermalConductivity;
	double rho = material.rho_Density;

	double h = params.h_xStep;
	double delta_t = params.delta_tStep;

	/*Решение по явной схеме*/
	for (int k = 0; k < K; ++k)
	{
		for (int i = 1; i < N - 1; ++i)
		{
			T[k + 1][i] = (h * h * c_v * T[k][i] + delta_t * lambda * (T[k][i + 1] - 2 * T[k][i] + T[k][i - 1])) / (rho * c_v * h * h);
			//T[k + 1][i] = (delta_t / (h*h)) * (T[k][ i + 1] - 2 *T[k][ i] + T[k][ i - 1]) + T[k][ i];
		}
	}
}

void HeatEquationSolver::getSolutionImplicit(const MaterialsData& material, double(*T)[N], double* x, double* t, const IntegrateParameters& params)
{
	/*Задали начальные условия*/
	for (int i = 0; i < N; ++i)
	{
		T[0][i] = InitialValueFunction(x[i]);
	}

	/*Задали краевые условия*/
	for (int k = 0; k < K; ++k)
	{
		T[k][0] = LeftBoundaryConditionFunction(t[k]);
		T[k][N] = RightBoundaryConditionFunction(t[k]);
	}

	double c_v = material.cv_HeatCapacity;
	double lambda = material.lambda_ThermalConductivity;
	double rho = material.rho_Density;

	double h = params.h_xStep;
	double delta_t = params.delta_tStep;

	double A = (lambda) / (h * h);
	double C = (lambda) / (h * h);
	double B = (rho * c_v * h * h + 2 * lambda * delta_t) / (delta_t * h * h);
	double F[N], P[N], Q[N];

	for (int k = 0; k < K - 1; ++k)
	{

		for (int i = 0; i < N; i++)
		{
			F[i] = T[k][i] * ((rho * c_v) / (delta_t));
		}
		P[0] = C / B;
		Q[0] = F[0] / B;

		for (int j = 1; j < N; ++j)
		{
			P[j] = C / (B - A * P[j - 1]);
			Q[j] = (F[j] + A * Q[j - 1]) / (B - A * P[j - 1]);
		}

		for (unsigned int j = N - 2; j > 0; --j)
		{
			T[k + 1][j] = P[j] * T[k + 1][j + 1] + Q[j];
		}
	}


}


int main()
{
	using namespace MaterialDataBase;

	HeatEquationSolver mySolver;
	IntegrateParameters params;

	params.h_xStep = 0.1;
	params.delta_tStep = 0.01;

	double T[K][N] = { 0 }; double x[N] = { 0 }; double t[K] = { 0 };

	mySolver.SetCoordinateXGrid(x, params);
	mySolver.SetTimeGrid(t, params);

	/*Здесь можно задать материал (см. enum)*/
	mySolver.getSolutionImplicit(materials[Ones], T, x, t, params);

	for (int i = 0; i < N; i++)
		printf("%lf ", x[i]);
	printf("\n\n");
	for (int k = 0; k < K; k++)
		printf("%lf ", t[k]);
	printf("\n\n");

	for (int i = 0; i < K; ++i)
	{
		for (int j = 0; j < N; ++j)
			printf("%lf\t", T[i][j]);
		printf("\n\n");
	}

	return 0;
}