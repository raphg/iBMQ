/*
 * ARS.h
 *
 *  Created on: Jan 24, 2012
 *      Author: hoblitz
 */

#ifndef ARS_H_
#define ARS_H_

#define NMAX 500

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
#include "RngStream.h"
#include <math.h>
#include <float.h>

double sample_conditional(double* restrict x,
		int* restrict num_x,
		int nmax,
		double* restrict argvec,
		ARS_workspace *ws,
		RngStream rng,
		double (*h)(const double, const double *),
		double (*h_prime)(const double , const double *));

int update_hull(double* restrict x,
		ARS_workspace *ws,
		const double* restrict argvec,
		int* restrict num_x,
		const int nmax,
		const double xnew,
		const double hnew,
		const int l_section,
		double* restrict huzmax,
		double (*h)(const double, const double *),
		double (*h_prime)(const double , const double *));

double sample_hull(const double* restrict x,
		ARS_workspace *ws,
		const int* restrict num_x,
		int* restrict section,
		const double p,
		const double huzmax);

void initialize_hull(double* restrict x,
		ARS_workspace *ws,
		int num_x,
		double huzmax);

void check_sample(const double x_samp,
		const double *x, const ARS_workspace *ws,
		const int *num_x);

void print_hull(double *x,
		ARS_workspace *ws,
		int *num_x);

double log_apb(const double loga,
		const double logb);
