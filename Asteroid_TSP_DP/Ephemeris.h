#ifndef Ephemeris_H
#define Ephemeris_H

#include <vector>
#include <complex>
using namespace std;

namespace Ephemeris{

	void Calendar(double juday, double* year, double* month, double* day, double* hour, double* min, double* sec, int LST, double longitude);

	double JD_epf(double year,double month,double day,double hour);

	FILE* INITIAL_DE(int DE);
	void Ephemeris    (FILE *fp_Bin,double  jd,int targ,int cent,double *rrd,int *inside);
	void Cartesian_coordinate     (double* Kep, double JD0, double JD, double* Dec);
	void Ephemerisx6(FILE *fp_Bin, double  jd, int targ, int cent, double *rrd, int *inside, double* Kep, double jd0, int dimension);

    void NAME                (char * Name, int &n, int i);
	void Kep_elements_MPCORB (char * Name, double* Kep, double &jd0);
	void Kep_elements_MPCORB_data(char Name[][19], double Kep[][6], double* jd0, double a[2], double e[2], double i[2], int nn[2], int n);
	void Kep_elements_MPCORB_dataName(char Name[][19], double Kep[][6], double* jd0, int n);

	void Keplerian_elements(double* Dec, double  t, double mu, double* Kep, int par);

}

#endif	//Equations_H