/*
 * iBMQ_common.h
 *
 *  Created on: Feb 28, 2012
 *      Author: hoblitz
 */

#ifndef IBMQ_COMMON_H_
#define IBMQ_COMMON_H_


#include <R.h>
#include <Rmath.h>
#include "RngStream.h"
#include "norm_gamma_generation.h"
#include "sparse.h"
#include "ARS.h"

#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <Rversion.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_check_range.h>


#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif

#define ZZERO 2.0e-308

void store_prob_include(int *n_iter, int *n_snps, int *n_genes, int** ProbSum, double* outProbs);

void update_prob_include(int* n_snps, int* n_genes, int** Gamma, int** ProbSum);

inline double log1m_from_logit(double x);

inline double log_from_logit(double x);

inline double expit(double x)
{
	double expit;
	expit = 1.0/(1.0 + exp(-x));
	return(expit);
}

inline double log_from_logit(double x)
{
	if(x > 0.0)
	{
		return(-log1p(exp(-x)));
	}
	else
	{
		return(x - log1p(exp(x)));
	}
}

inline double log1m_from_logit(double x)
{
	if(x > 0.0)
	{
		return(-x - log1p(exp(-x)));
	}
	else
	{
		return(-log1p(exp(x)));
	}
}


int update_pos_j(double* P, double* A, double* B, double* C, double** W_Logit,
		int** W_Ind, int** Gamma,
		int j, double* a_0, double* b_0, double* lambda_a, double* lambda_b,
		int* n_snps, int* n_genes, int* n_indivs, RngStream rng, int nmax,
		double *xA, double *xB, ARS_workspace *workspace, double eps);

void set_prior(double* lambda_a, double* lambda_b, double* a_0, double* b_0, double *tau_0,
		double* expr_means, double* expr_vars, double* alpha2_beta,
		gsl_matrix* X, gsl_matrix* Y, RngStream rng);

void initialize_parms(m_el **Beta, int** Gamma, double** W_Logit, int** W_Ind, int** ProbSum, double **xA, double **xB,
		double* A, double* B, double* C, double* P, double* Mu, double* Sig2,
		double* expr_means, double* expr_vars, double* alpha2_beta,
		double* lambda_a, double* a_0, double* lambda_b, double* b_0, double* tau_0,
		int* n_snps, int* n_genes, int* n_indivs, int* nmax, RngStream rng);

inline double expit(double x);
#endif /* IBMQ_COMMON_H_ */


