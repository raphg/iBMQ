/*
 * ARS.h
 *
 *  Created on: Jan 24, 2012
 *      Author: hoblitz
 */

#ifndef ARS_H_
#define ARS_H_


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

double sample_conditional(double* x, int* num_x, int nmax, double S_j, double c, double R, double *hwv,
		double *hpwv, double *scum, double* scum_norm, double *z, double *s, RngStream rng,
		double eps);

int update_hull(double *hwv, double *hpwv, double *x, double *z, int *num_x, int nmax,
		double xnew, double hnew, int l_section,
		double R, double S_j, double c);

double sample_hull(double *hwv, double *hpwv, double *x, double *z, double *scum, double *s,
		double * scum_norm, int* num_x, int *section, double p);

void initialize_hull(double *hwv, double *hpwv, double *x, double *z, double *scum, double *s,
		double *scum_norm, int* num_x);

double h_prime(double x, double S_j, double c, double R);

double h(double x, double S_j, double c, double R);

