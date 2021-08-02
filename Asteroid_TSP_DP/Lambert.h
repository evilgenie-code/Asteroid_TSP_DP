// ------------------------------------------------------------------------ //
// This source file is part of the 'KOSMOEXPORT Team's			            //
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2012-2012 FEDERAL SPACE AGENCY                             //
// ------------------------------------------------------------------------ //
#ifndef LAMBERT_H
#define LAMBERT_H

void getCrossProductVectors(const double*, const double*, double*);
void vers(const double*, double*);
void LambertI (const double *, const double *, double , const double , const int , int , const int ,//INPUT
	           double *, double *, double &, double &, double &, int &);//OUTPUT
double acosh (double );
double asinh (double);
double x2tof(const double &x,const double &s,const double &c,const int lw, int N);

void Lambert(double *r1, double *r2, double tof, double mu, int Nrev,
			 double *V1_total,double *V2_total,double *a_total,int *Nrev_total, int &iteration_total, int &count);

void theta_f    (const double *r1_in, const double *r2_in, double mu, double &theta);
void theta_f180 (const double *r1_in, const double *r2_in, double mu, double &theta);

#endif	//LAMBERT_H

