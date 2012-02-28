/*
 * iBMQ_common.c
 *
 *  Created on: Feb 28, 2012
 *      Author: hoblitz
 */

#include "iBMQ_common.h"
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

//________________store prob_include_____________//
void store_prob_include(int *n_iter, int *n_snps, int *n_genes, int** ProbSum, double* outProbs)
{
	int g,j;
	for(g = 0; g < *n_genes; g++)
	{
		for(j = 0; j < *n_snps; j++)
		{
			outProbs[*n_snps*g + j] = (double)(ProbSum[j][g])/(double)(*n_iter);
		}
	}
	return;
}


void update_prob_include(int* n_snps, int* n_genes, int** Gamma, int** ProbSum)
{
	int g,j;
	for(j = 0; j < *n_snps; j++)
	{
		for(g = 0; g < *n_genes; g++)
		{
			ProbSum[j][g] += Gamma[j][g];
		}
	}
	return;
}


int update_pos_j(double* P, double* A, double* B, double* C, double** W, int** Gamma,
		int j, double* a_0, double* b_0, double* lambda_a, double* lambda_b,
		int* n_snps, int* n_genes, int* n_indivs, RngStream rng, int nmax,
		double *xA, double *xB, ARS_workspace *workspace, double eps)
{
	int S_j = 0;
	int g, num_x = 2;
	double s_log_w = 0.0, s_log_1minus_w = 0.0, ratio;
	double R = 0.0;

	// count the number of genes turned on for snp j
	for(g = 0; g < *n_genes; g++)
	{
		S_j += (int)(W[j][g] > ZZERO);
	}

	if(S_j == 0)
	{
		A[j] = RngStream_GA1(1.0, rng)/(*lambda_a);
		B[j] = RngStream_GA1(1.0, rng)/(*lambda_b);
	}

	else
	{
		// update a_j and b_j using ARS; need these quantities
		for(g = 0; g < *n_genes; g++)
		{
			if(W[j][g] >= ZZERO)
			{
				s_log_w += log(W[j][g]);
			}
			if((W[j][g] >= ZZERO) & ((1.0 - W[j][g]) > ZZERO))
			{
				s_log_1minus_w += log(1.0 - W[j][g]);
			}
		}
		R = (*lambda_a - s_log_w);
		A[j] = sample_conditional(xA, &num_x, nmax, (double)S_j, B[j], R, workspace, rng, eps);

		if(A[j] == -1.0)
		{
			return(0);
		}

		num_x = 2;
		R = (*lambda_b - s_log_1minus_w);
		B[j] = sample_conditional(xB, &num_x, nmax, (double)S_j, A[j], R, workspace, rng, eps);
		if(B[j] == -1.0)
		{
			return(0);
		}
	}

	P[j] = RngStream_Beta(*a_0 + (double)(*n_genes - S_j), *b_0 + (double)(S_j), rng);

	// update W, row j (going across genes here)
	double frac;
	frac = B[j]/(A[j] + B[j]);

	for(g = 0; g < *n_genes; g++)
	{
		ratio = (P[j]*(double)(Gamma[j][g]==0))/
				(P[j]*(double)(Gamma[j][g]==0) + (1.0 - P[j])*frac);

		if(RngStream_RandU01(rng) <= ratio)
		{
			W[j][g] = 0.0;
		}
		else
		{
			W[j][g] = RngStream_Beta(A[j] + (double)(Gamma[j][g]),
					B[j] + 1.0 - (double)(Gamma[j][g]), rng);
		}
	}
	return(1);
}

void set_prior(double* lambda_a, double* lambda_b, double* a_0, double* b_0, double *tau_0,
		double* expr_means, double* expr_vars, double* alpha2_beta,
		gsl_matrix* X, gsl_matrix* Y, RngStream rng)
{
/* TODO: allow hyperparameters to be passed as a function argument
 *
 */
	*lambda_a = 10.0;
	*lambda_b = 0.1;
	*a_0 = 10.0;
	*b_0 = 0.1;
	*tau_0 = 1000.0;

	double out;
	int g,j, n_genes = Y->size2, n_indivs = Y->size1, n_snps = X->size2;
	gsl_vector_view tmp;

	//compute emprical mean and variances of gene expression phenotype
	for(g = 0; g < n_genes; g++)
	{
		tmp = gsl_matrix_column(Y, g);
		expr_means[g] = gsl_stats_mean(tmp.vector.data, tmp.vector.stride, n_indivs);
		expr_vars[g] = gsl_stats_variance_m(tmp.vector.data, tmp.vector.stride, n_indivs,
				expr_means[g]);
	}

	// compute inverse of inner product of X matrix values
	for(j = 0; j < n_snps; j++)
	{
		tmp = gsl_matrix_column(X, j);
		gsl_blas_ddot(&tmp.vector, &tmp.vector, &out);
		alpha2_beta[j] = 1.0/(out);
	}
	return;
}

void initialize_parms(m_el **Beta, int** Gamma, double** W, int** ProbSum, double **xA, double **xB,
		double* A, double* B, double* C, double* P, double* Mu, double* Sig2,
		double* expr_means, double* expr_vars, double* alpha2_beta,
		double* lambda_a, double* a_0, double* lambda_b, double* b_0, double* tau_0,
		int* n_snps, int* n_genes, int* n_indivs, int* nmax, RngStream rng)
{
	int j, g;
	for(g = 0; g < *n_genes; g++)
	{
		Mu[g] = expr_means[g] + sqrt(Sig2[g])*RngStream_N01(rng);
		Sig2[g] = expr_vars[g];
		C[g] = (double)(*n_indivs);
	}

	for(j = 0; j < *n_snps; j++)
	{
		// sample directly from prior distributions
		A[j] = RngStream_GA1(1.0, rng)/ (*lambda_a);
		B[j] = RngStream_GA1(1.0, rng)/ (*lambda_b);
		P[j] = RngStream_Beta(*a_0,*b_0, rng);

		for(g = 0; g < *n_genes; g++)
		{
			ProbSum[j][g] = 0;
			if(RngStream_RandU01(rng) <= P[j])
			{
				W[j][g] = 0.0;
			}
			else
			{
				W[j][g] = RngStream_Beta(A[j], B[j], rng);
			}
		}

		for(g = 0; g < *nmax; g++)
		{
			xA[j][g] = 0.0;
			xB[j][g] = 0.0;
		}
		xA[j][0] = .1;
		xB[j][0] = .1;
		xA[j][1] = 	2.0;
		xB[j][1] =  2.0;

		double x;
		for(g = 0; g < *n_genes; g++)
		{
			if(RngStream_RandU01(rng) <= W[j][g])
			{
				Gamma[j][g] = 1;
				x = Sig2[g]*alpha2_beta[j]*C[g];
				SV_add_el(Beta[g], j, RngStream_N01(rng)*sqrt(x));
			}
			else
			{
				Gamma[j][g]= 0;
			}
		}
	}
	return;
}


