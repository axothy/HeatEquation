#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <array>


#define N 7
#define K 11

/*just enumeration of various materials*/
enum ListOfSubstances {
	StainlessSteel = 0,
	Iron = 1,
	Gold = 2,
	Brass = 3,
	Copper = 4,
	Mercury = 5,
	Silver = 6,
	Diamond = 7,
	Glass = 8,
	Vacuum = 9,
	Ones = 10,
	Wood = 11
};

/*struct consists of material parameters*/
struct MaterialsData {
	const double cv_HeatCapacity;
	const double lambda_ThermalConductivity;
	const double rho_Density;
};

/*prepared material parameters in container for materials from ListOfSubstances*/
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
		{0, 0, 0},
		{1, 1, 1},
		{2300, 110, 470 }
	};
}

/*integrate steps*/
struct IntegrateParameters {
	double delta_tStep;
	double h_xStep;
};

class HeatEquationSolver {
	//HeatEquationSolver(); //class constructor
	//~HeatEquationSolver() {}; //destructor

private:

	inline double InitialValueFunction(double x);
	inline double LeftBoundaryConditionFunction(double t);
	inline double RightBoundaryConditionFunction(double t);

public:

	void SetCoordinateXGrid(std::array<double, N>& x, const IntegrateParameters& params);
	void SetTimeGrid(std::array<double, K>& t, const IntegrateParameters& params);

	void getSolutionExplicit(const MaterialsData& material, double(*Temperature)[N], std::array<double, N>& x, std::array<double, K>& t,
		const IntegrateParameters& params);
	void getSolutionImplicit(const MaterialsData& material, double(*Temperature)[N], std::array<double, N>& x, std::array<double, K>& t,
		const IntegrateParameters& params);
};


/* NU T(x,0) */
inline double HeatEquationSolver::InitialValueFunction(double x)
{
	return 0.9 + 2 * x * (1 - x);
}

/* T(0, t) */
inline double HeatEquationSolver::LeftBoundaryConditionFunction(double t)
{
	return 3 * (0.3 - 2 * t);
}

/* GU2 T(l,t) */
inline double HeatEquationSolver::RightBoundaryConditionFunction(double t)
{
	return 1.38;
}

/*creating discrete grids*/
void HeatEquationSolver::SetCoordinateXGrid(std::array<double, N>& x, const IntegrateParameters& params)
{
	x[0] = 0;
	for (uint8_t i = 1; i < N; ++i) {
		x[i] = i * params.h_xStep;
	}
}

void HeatEquationSolver::SetTimeGrid(std::array<double, K>& t, const IntegrateParameters& params)
{
	t[0] = 0;
	for (uint8_t k = 1; k < K; ++k) {
		t[k] = k * params.delta_tStep;
	}
}

/*solution by explicit finite difference method*/
void HeatEquationSolver::getSolutionExplicit(const MaterialsData& material, double (*Temperature)[N], std::array<double, N>& x,
	std::array<double, K>& t, const IntegrateParameters& params)
{
	/*setting initial condition*/
	for (uint8_t i = 0; i < N; ++i)
	{
		Temperature[0][i] = InitialValueFunction(x[i]);
	}

	/*setting boundary conditions*/
	for (uint8_t k = 0; k < K; ++k)
	{
		Temperature[k][0] = LeftBoundaryConditionFunction(t[k]);
		Temperature[k][N-1] = RightBoundaryConditionFunction(t[k]);
	}
	
	/*to be short*/
	double c_v = material.cv_HeatCapacity;
	double lambda = material.lambda_ThermalConductivity;
	double rho = material.rho_Density;

	double h = params.h_xStep;
	double delta_t = params.delta_tStep;

	/*explicit formule*/
	for (uint8_t k = 0; k < K-1; ++k)
	{
		for (uint8_t i = 1; i < N - 1; ++i)
		{
			Temperature[k + 1][i] = (h * h * c_v * Temperature[k][i] + delta_t * lambda *
				(Temperature[k][i + 1] - 2 * Temperature[k][i] + Temperature[k][i - 1])) / (rho * c_v * h * h);
		}
	} 
	
}

void HeatEquationSolver::getSolutionImplicit(const MaterialsData& material, double(*Temperature)[N], std::array<double, N>& x,
	std::array<double, K>& t, const IntegrateParameters& params)
{
	/*setting initial condition*/
	for (uint8_t i = 0; i < N; ++i)
	{
		Temperature[0][i] = InitialValueFunction(x[i]);
	}

	/*setting boundary conditions*/
	for (uint8_t k = 0; k < K; ++k)
	{
		Temperature[k][0] = LeftBoundaryConditionFunction(t[k]);
		Temperature[k][N-1] = RightBoundaryConditionFunction(t[k]);
	}

	double c_v = material.cv_HeatCapacity;
	double lambda = material.lambda_ThermalConductivity;
	double rho = material.rho_Density;

	double h = params.h_xStep;
	double delta_t = params.delta_tStep;

	/*setting coefficients of finite difference equation*/
	double A = (lambda) / (h * h);
	double C = (lambda) / (h * h);
	double B = (rho * c_v * h * h + 2 * lambda * delta_t) / (delta_t * h * h);
	double F[N], P[N], Q[N];

	for (uint8_t k = 0; k < K - 1; ++k)
	{

		for (uint8_t i = 0; i < N; i++)
		{
			F[i] = Temperature[k][i] * ((rho * c_v) / (delta_t));
		}
		/*Progonka algorithm coefficients*/
		P[0] = C / B;
		Q[0] = F[0] / B;

		for (uint8_t j = 1; j < N; ++j)
		{
			P[j] = C / (B - A * P[j - 1]);
			Q[j] = (F[j] + A * Q[j - 1]) / (B - A * P[j - 1]);
		}

		/*back substitution*/
		for (uint8_t j = N - 2; j > 0; --j)
		{
			Temperature[k + 1][j] = P[j] * Temperature[k + 1][j + 1] + Q[j];
		}
	}


}


int main()
{
	using namespace MaterialDataBase;

	HeatEquationSolver mySolver;
	IntegrateParameters params;

	params.h_xStep = 0.1;
	params.delta_tStep = 0.001;

	double Temperature[K][N] = { 0 };
	std::array<double, N> x;
	std::array<double, K> t;

	mySolver.SetCoordinateXGrid(x, params);
	mySolver.SetTimeGrid(t, params);

	/*Здесь можно задать материал (см. enum)*/
	mySolver.getSolutionImplicit(materials[Ones], Temperature, x, t, params);

	for (uint8_t i = 0; i < N; i++)
		printf("%lf ", x[i]);
	printf("\n\n");
	for (uint8_t k = 0; k < K; k++)
		printf("%lf ", t[k]);
	printf("\n\n");

	for (uint8_t i = 0; i < K; ++i)
	{
		for (uint8_t j = 0; j < N; ++j)
			printf("& %lf\t", Temperature[i][j]);
		printf("\n\n");
	} 

	return 0;
}
