/*
 * ARS.h
 *
 *  Created on: Jan 24, 2012
 *      Author: hoblitz
 */

#ifndef ARS_H_
#define ARS_H_

#define NMAX 10

struct ARS_WORKSPACE{
	double hwv[NMAX];
	double hpwv[NMAX];
	double scum[NMAX];
	double scum_norm[NMAX];
	double s[NMAX];
	double z[NMAX];
};

typedef struct ARS_WORKSPACE ARS_workspace;

#endif /* ARS_H_ */

#include <R.h>
#include <Rmath.h>
#include "RngStream.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Rversion.h>

#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

double sample_conditional(double* x, int* num_x, int nmax, double S_j, double c, double R, ARS_workspace *space,
		RngStream rng, double eps);

int update_hull(double *x, ARS_workspace *space, int *num_x, int nmax,
		double xnew, double hnew, int l_section,
		double R, double S_j, double c, double *huzmax);

double sample_hull(double *x, ARS_workspace *space, int* num_x, int *section, double p, double *huzmax);

void initialize_hull(double *x, ARS_workspace *space, int* num_x, double* huzmax);

double h_prime(double x, double S_j, double c, double R);

double h(double x, double S_j, double c, double R);

