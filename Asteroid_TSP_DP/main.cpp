#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h> 
#include <stdio.h>
#include "Constants.h"
#include "Ephemeris.h"
#include "Lambert.h"
#include "fmincg.h"

#define NN 10 //countObjects
#define MAX 0x7fffffff
#define _CRT_SECURE_NO_WARNINGS

using namespace std;
using std::string;

FILE *fp_Bin;
double UnitM, UnitP, Unitb, UnitN;
//extern double Kep_elements[100][6], jd0_elements[100];
double Kep[NN][6], JD0[NN];
char Name[NN][19];
double KepTrajectory[NN][6]; // Момент времени всегда 0, и связан с датой пролёта
double V0Vk[NN][6];

void transfer(int& start, int& finish, double& jd, double& dt, double v0[3], double* vs, double* vk, double& dV, double dv1[3], double dv2[3]);
double TSP(int& start, double& jd0, double dVlim, int n,														// INPUT
	double& dVsum, int* PATH, double* dVpath, double* JDpath, int& ip, int path[][1 << (NN - 1)]);			// OUTPUT
void print_PATH(int* PATH, double* dVpath, double* JDpath, int ip, double jd0);

static double sign(double a, double b)
{
	return (b < 0.0) ? -fabs(a) : fabs(a);
}

int getDirectionMovement(const double* crossProductRadii)
{
	return crossProductRadii[2] >= 0.0 ? 0 : 1;
}

void costFunc(double* xVec, double* cost, double* gradVec, int* pari, double *pard)
{
//  xVec = dt
//  pard[] = jd, v1[0:3] - момент отлёта и скорость в момент отлёта
//  pari[2] = n0, nk
	static double crossProductRadii[3], radiusVectorStart[6], radiusVectorFinish[6], v11[3], v12[3], v2[3], dv11[3], dv12[3], dv1[3], a, p, theta, mu = 1., epsilon;
	static int directionMovement, iterations, countRevolutions, i, nctr = 11, inside;

	if (*xVec < 0) 
	{
		*cost = 1000;
		*gradVec = 100;
		
		return;
	}

	epsilon = *xVec*1.e-4;
	iterations = 0;
	countRevolutions = 0;

	Ephemeris::Ephemerisx6(fp_Bin, pard[0]            , pari[0], nctr, radiusVectorStart, &inside, &Kep[pari[0] - 100][0], JD0[pari[0] - 100], 0);
	Ephemeris::Ephemerisx6(fp_Bin, pard[0]+*xVec*UnitT, pari[1], nctr, radiusVectorFinish, &inside, &Kep[pari[1] - 100][0], JD0[pari[1] - 100], 0);

	getCrossProductVectors(radiusVectorStart, radiusVectorFinish, crossProductRadii);
	
	directionMovement = getDirectionMovement(crossProductRadii);
	
//	находит решение с заданным числом витков
	LambertI(radiusVectorStart, radiusVectorFinish, *xVec, mu, directionMovement, countRevolutions, 0,		// INPUT
				 v11, v2, a, p, theta, iterations);		// OUTPUT
	
	for (i = 0; i < 3; i++)
		dv11[i] = v11[i] - pard[1+i];

	Ephemeris::Ephemerisx6(fp_Bin, pard[0] + (*xVec + epsilon)*UnitT, pari[1], nctr, radiusVectorFinish, &inside, &Kep[pari[1] - 100][0], JD0[pari[1] - 100], 0);

	getCrossProductVectors(radiusVectorStart, radiusVectorFinish, crossProductRadii);

	directionMovement = getDirectionMovement(crossProductRadii);
	
	LambertI(radiusVectorStart, radiusVectorFinish, *xVec+ epsilon, mu, directionMovement, countRevolutions, 0,	// INPUT
					v12, v2, a, p, theta, iterations);	// OUTPUT

	for (i = 0; i < 3; i++)
		dv12[i] = v12[i] - pard[1 + i];

	for (i = 0; i < 3; i++)
		dv1[i] = (dv11[i] - dv12[i]) / epsilon;


	*cost = dv11[0]* dv11[0] + dv11[1] * dv11[1] + dv11[2] * dv11[2]; // первый импульс

	*gradVec = 2.0*dv11[0]*dv1[0] + 2.0*dv11[1]*dv1[1] + 2.0*dv11[2]*dv1[2];
}

void transfer(int &start, int &finish, double &jd, double &dt, double v0[3], double* vs, double* vk, double &dV, double dv1[3], double dv2[3]) {

	static double a1, a2, dt1, f1, df1, dt2, f2, df2, pard[4];
	static int pari[2];

	static double r[3], r0[6], rk[6], v1[3], v2[3], a, p, theta, mu = 1.;
	static int lw, iter, nrev, i, nctr = 11, inside;

	if (start == 13)
		a1 = 1.0;
	else
		a1 = Kep[start - 100][0];

	if (finish == 13)
		a2 = 1.0;
	else
		a2 = Kep[finish - 100][0];

	a1 = PI2 * sqrt(a1*a1*a1);
	a2 = PI2 * sqrt(a2*a2*a2);

	dt1 = (a1 + a2) * 0.4;
	dt2 = (a1 + a2) * 0.6;

	pari[0] = start;
	pari[1] = finish;

	pard[0] = jd;
	pard[1] = v0[0]; pard[2] = v0[1]; pard[3] = v0[2];
	int ret1 = fmincg(&costFunc, &dt1, pari, pard, 1, 100, f1, &df1); // nDim = 1 , maxCostFunctionCalls = 100
	int ret2 = fmincg(&costFunc, &dt2, pari, pard, 1, 100, f2, &df2); // nDim = 1 , maxCostFunctionCalls = 100

	if (f1 < f2) {
		dt = dt1;
		dV = sqrt(f1);
	}
	else {
		dt = dt2;
		dV = sqrt(f2);
	}

	iter = 0;
	nrev = 0;	// число витков

	Ephemeris::Ephemerisx6(fp_Bin, jd, start, nctr, r0, &inside, &Kep[start-100][0], JD0[start - 100], 0);
	
	Ephemeris::Ephemerisx6(fp_Bin, jd + dt * UnitT, finish, nctr, rk, &inside, &Kep[finish - 100][0], JD0[finish - 100], 0);

	getCrossProductVectors(r0, rk, r);
	if (r[2] >= 0.0) lw = 0;
	else lw = 1;
//	находит решение с заданным числом витков
	LambertI(r0, rk, dt, mu, lw, nrev, 0,			// INPUT
			 v1, v2, a, p, theta, iter);			// OUTPUT

	for (i = 0; i < 3; i++) {
		dv1[i] = v1[i] - v0[i];
		dv2[i] = rk[3 + i] - v1[i];
		vs[i] = v1[i];
		vk[i] = v2[i];
	}
}

void get_rowAM(int start, double jd0, double v0[][3], double DV[NN], double DV0[][3], double DVk[][3], double DT[NN], int parents[NN], int n, bool OK = 1)
{

	static int i, j;
	static double dv1[3], dv2[3];

	if (OK) {
		for (i = 100; i < n + 100; i++) 
		{
			if (parents[i - 100] == -1)
			{
				transfer(start, i, jd0, DT[i - 100], v0[start-100], &DV0[i - 100][0], &DVk[i - 100][0], DV[i - 100], dv1, dv2);
				for (j = 0; i < 3; i++) v0[i - 100][j] = DVk[i - 100][j];
			}
			DV[i - 100] = MAX;
		}
	}
	else 
	{
		for (i = 100; i < n + 100; i++)
		{
			if (parents[i - 100] == -1) 
			{
				transfer(i, start, jd0, DT[i - 100], v0[i - 100], &DV0[i - 100][0], &DVk[i - 100][0], DV[i - 100], dv1, dv2);
			}
			DV[i - 100] = MAX;
		}
	}
}

int removeCity(int j, int k, int n)// Remove the k city from the set of cities represented by j binary (k bit is set to 0)
{ 
	int allCities = (1 << (n - 1)) - 1;
	int onlyKCity = 1 << (k - 1);
	int mask = allCities ^ onlyKCity; // XOR (all 1 except for the k position)

	return j & mask;
}

void print_PATH(int* PATH, double* dVpath, double* JDpath, int ip, double jd0) {

	int i, j, n = 100;
	int nctr = 11, inside;
	double r0[6], a, jd;
	double year, month, day, hour, min, sec;

	FILE* fileout;
	char nameFile[20] = { "PATH.txt" };

	fileout = fopen(nameFile, "w");
	fprintf(fileout, "%25s %25s %25s %25s %25s %25s ", "x1, UnitR", "y1, UnitR", "z1, UnitR", "Vx1, UnitV", "Vy1, UnitV", "Vz1, UnitV");
	fprintf(fileout, "%25s %25s %19s %25s", "Name", "jd, day", "Cal", "dV, km/sec");
	fprintf(fileout, "\n");

	for (i = 0; i < ip; i++) {
		printf("%i", i);
		//Ephemeris::Ephemerisx6(fp_Bin, JDpath[i], PATH[i] + 100, nctr, r0, &inside, &Kep[PATH[i]][0], JD0[PATH[i]], 0);
		Ephemeris::Calendar(JDpath[i], &year, &month, &day, &hour, &min, &sec, 0, 0);

		for (j = 0; j < 20; j++)
			nameFile[j] = '\0';

		if (PATH[i] + 100 == 13) {
			a = 1.0;
			nameFile[0] = 'E';
			nameFile[1] = 'a';
			nameFile[2] = 'r';
			nameFile[3] = 't';
			nameFile[4] = 'h';
			j = 5;
		}
		else {
			a = Kep[PATH[i]][0];
			for (j = 0; Name[i][j] != 32 && Name[i][j]; j++)
				nameFile[j] = Name[i][j];
		}

		fprintf(fileout, "%25.16e %25.16e %25.16e %25.16e %25.16e %25.16e ", r0[0], r0[1], r0[2], r0[3], r0[4], r0[5]);
		fprintf(fileout, "%25s %25.16e ", nameFile, JDpath[i]);
		fprintf(fileout, "%04i.%02i.%02i %02i:%02i:%02i ", int(year), int(month), int(day), int(hour), int(min), int(sec));
		fprintf(fileout, "%25.16e", dVpath[i]*UnitV);
		fprintf(fileout, "\n");
	}
	fclose(fileout);

	for (i = 0; i < ip; i++) {

		for (j = 0; j<20; j++)
			nameFile[j] = '\0';

		if (PATH[i] + 100 == 13) {
			a = 1.0;
			nameFile[0] = 'E';
			nameFile[1] = 'a';
			nameFile[2] = 'r';
			nameFile[3] = 't';
			nameFile[4] = 'h';
			j = 5;
		}
		else {
			a = Kep[PATH[i]][0];
			for (j = 0; Name[i][j] != 32 && Name[i][j]; j++)
				nameFile[j] = Name[i][j];
		}

		a = PI2 * sqrt(a * a * a)*UnitT;
		a = a / n;

		nameFile[j] = '.';
		nameFile[j+1] = 't';
		nameFile[j+2] = 'x';
		nameFile[j+3] = 't';

		fileout = fopen(nameFile, "w");
		fprintf(fileout, "%25s %25s %25s %25s %25s %25s", "x1, UnitR", "y1, UnitR", "z1, UnitR", "Vx1, UnitV", "Vy1, UnitV", "Vz1, UnitV");
		fprintf(fileout, "\n");

		for (j = 0; j <= n; j++) {
			Ephemeris::Ephemerisx6(fp_Bin, jd0+a*j, PATH[i] + 100, nctr, r0, &inside, &Kep[PATH[i]][0], JD0[PATH[i]], 0);
			fprintf(fileout, "%25.16e %25.16e %25.16e %25.16e %25.16e %25.16e", r0[0], r0[1], r0[2], r0[3], r0[4], r0[5]);
			fprintf(fileout, "\n");
		}
		fclose(fileout);
	}

	fileout = fopen("Trajectory.txt", "w");
	fprintf(fileout, "%25s %25s %25s %25s %25s %25s %25s ", "jd, day", "x1, UnitR", "y1, UnitR", "z1, UnitR", "Vx1, UnitV", "Vy1, UnitV", "Vz1, UnitV");
	fprintf(fileout, "\n");


	for (i = 0; i < ip-1; i++) {

		a = JDpath[i+1] - JDpath[i];
		a = a / n;

		for (j = 0; j <= n; j++) {
			Ephemeris::Ephemerisx6(fp_Bin, a*j, 100, nctr, r0, &inside, &KepTrajectory[i][0], 0, 0);
			fprintf(fileout, "%25.3f %25.16e %25.16e %25.16e %25.16e %25.16e %25.16e", JDpath[i] + a * j, r0[0], r0[1], r0[2], r0[3], r0[4], r0[5]);
			fprintf(fileout, "\n");
		}
	}
	fclose(fileout);
}
							
double TSP(int& start, double& jd0, double dVlim, int n,															// INPUT
	double& dVsum, int* PATH, double* dVpath, double* JDpath, int& ip, int path[][1 << (NN - 1)])
{			// OUTPUT

	double v0[NN][3], DV[NN], DV0[NN][3], DVk[NN][3], DT[NN], jd[NN], dVmin, dt = 0;
	double temp;
	int parents[NN], i, j, k, min, nextCity;

	double** dp = new double*[n + 1];
	
	for (i = 0; i < n + 1; i++)
		dp[i] = new double[1 << (n - 1)]{ 0 };

	ip = 0;
	for (i = 0; i < n; i++)
		parents[i] = -1;

	int nctr = 11, inside;
	double r0[6];
	Ephemeris::Ephemerisx6(fp_Bin, jd0, start + 100, nctr, r0, &inside, &Kep[start - 100][0], JD0[start - 100], 0);
	for (i = 0; i < 3; i++)
		v0[0][i] = r0[3 + i];

	get_rowAM(start + 100, jd0, v0, DV, DV0, DVk, DT, parents, n);

	// Initialize the distance from all cities to city 0
	for (i = 1; i < n; i++) 
		dp[i][0] = DV[i];
		jd[i] = jd0 + DT[i];

	for (j = 1; j < (1 << (n - 1)); j++)
	{
		for (i = 1; i < n; i++)
		{
			if ((j >> (i - 1)) % 2 == 0)
			{
				// City i is not in the set of cities in j binary form
				min = MAX;
				nextCity = i;
				jd0 = jd[i];
				get_rowAM(i + 100, jd0, v0, DV, DV0, DVk, DT, parents, n);
				for (k = 1; k < n; k++) 
				{
					// City k is in the city set represented by j binary form
					if (j >> (k - 1) % 2 != 0) 
					{
						// means removing the k city from the j city set
						temp = DV[k] + dp[k][removeCity(j, k, n)];
						if (temp < min) 
						{
							min = temp;
							nextCity = k;
							dt = DT[k];
						}
					}
				}
				dp[i][j] = min;
				path[i][j] = nextCity;

				jd[nextCity] = jd0 + dt;
			}
		}
	}

	// Fill in the upper left corner element
	min = MAX;
	nextCity = 0;
	j = (1 << (n - 1)) - 1; // Represents the collection of all cities except city 0
	get_rowAM(i + 100, jd0, v0, DV, DV0, DVk, DT, parents, n);
	for (k = 1; k < n; k++) {
		temp = DV[k] + dp[k][removeCity(j, k, n)];
		if (temp < min) {
			min = temp;
			nextCity = k;
		}
	}
	dp[0][j] = min;
	path[0][j] = nextCity;

	return dp[0][j];
}

int main()
{
	int type;
	printf("Start! \n");
	printf("1. Data from MPCORB. 2. Data from file. \n");
	printf("Enter your type:");
	cin >> type;
	printf("Type = %i \n", type);

	const int de = 405;
	fp_Bin = Ephemeris::INITIAL_DE(de);
	
	int n = NN, nn[] = { 0, 10000 };// искать в первых 10000 объектах
	double aa[] = { 0.7, 3.5 }, e[] = { -0.1, 0.4 }, I[] = { -15 / raddeg, 15 / raddeg };
	//filters посмотреть в статье АиК

	if (type == 1)
	{										//name - nameObjects
		Ephemeris::Kep_elements_MPCORB_data(Name, Kep, JD0, aa, e, I, nn, n);
	}
	else if (type == 2)
	{
		n = 0;
		string nameFile;
		
		printf("Enter your name file, without extension:");
		cin >> nameFile;

		nameFile = nameFile + ".txt";
		
		ifstream filein(nameFile);
		while (filein.getline(Name[n], 19)) n++;
		filein.close();
		
		Ephemeris::Kep_elements_MPCORB_dataName(Name, Kep, JD0, n);
	}
	else
	{
		printf("Unknown type: %i \n", type);
		return 0;
	}

	int ip, PATH[50];
	double jd0 = Ephemeris::JD_epf(2023, 8, 29, 0);

	int start = 13 - 100;
	const double dVlim = 150 / UnitV;
	double dVsum = 0.0, dVpath[50] = { 0 }, JDpath[50] = { 0 };

	int path[NN][1 << (NN - 1)] = { 0 };
	
	if (n > NN) n = NN;
	TSP(start, jd0, dVlim, n,							// INPUT
		dVsum, PATH, dVpath, JDpath, ip, path);		// OUTPUT

	print_PATH(PATH, dVpath, JDpath, ip + 1, jd0);
	
	system("pause");
	return 0;
}