/*
Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13
%
%
% (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen
% 
% Permission is granted for anyone to copy, use, or modify these
% programs and accompanying documents for purposes of research or
% education, provided this copyright notice is retained, and note is
% made of any changes that have been made.
% 
% These programs and documents are distributed without any warranty,
% express or implied.  As the programs were written for research
% purposes only, they have not been tested to the degree that would be
% advisable in any important application.  All use of these programs is
% entirely at the user's own risk.

[ Sagar G V, 2013, sagar.writeme@gmail.com ] Changes Made:
- Ported to C

*/
#include<math.h>
#include<stdio.h>
#include<float.h>

#define COST_FUNC_DATATYPE double
#define COST_FUNC_DATATYPE_MIN (DBL_MIN*100)

#define RHO 0.01
#define SIG 0.5
#define INT 0.1
#define EXT 3.0
#define MAX 20
#define RATIO 100.0

// 1. pass in the cost function which takes in an array and gives out cost and gradient at the given input. 
// 2. xVector should contain the initial point which is will be modified to reflect the optimum point
// 3. nDim is the dimension of xVector
// 4. maxCostCalls is the maximum number of times the cost function may be called
// return value:  1 -> Num of Cost function calls exceeded max specified in the argument. 2-> line search failed

int fmincg(void (*costFunc)(COST_FUNC_DATATYPE* inputVector, COST_FUNC_DATATYPE* cost, COST_FUNC_DATATYPE* gradVector, int* pari, double *pard), COST_FUNC_DATATYPE* xVector, int* pari, double *pard, int nDim, int maxCostFuncCalls, double &f1, double *df1);

int trnlsp(void(*fun)(int *m, int *n, double *x, double *fvec), int &n, int &m, double *x, double *fvec, double *par, double eps[6], int iter1, int iter2, int &iter, int &st_cr, double &r1, double &r2);