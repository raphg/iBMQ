/*
 * main_parallel_sparse.c
 *
 *  Created on: Sep 20, 2011
 *      Author: Gregory Imholte
 */

#include <R.h>
#include <Rmath.h>
#include "RngStream.h"
#include "norm_gamma_generation.h"
#include "sparse.h"
#include "ARS.h"
#include "iBMQ_common.h"

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

void update_gene_g(m_el *beta_g, int** Gamma, double** W, gsl_matrix* X, gsl_vector* Y,
		double* A, double* B, double* C_g, double* P, double* Mu_g, double* Sig2_g,
		double* expr_means, double* expr_vars, double* alpha2_beta,
		double* lambda_a, double* a_0, double* lambda_b, double* b_0, double* tau_0,
		int* n_snps, int* n_genes, int* n_indivs, int g, gsl_vector* one, RngStream rng);

void store_mcmc_output(FILE *Afile, FILE *Bfile,
		FILE *Cfile, FILE *Pfile, FILE *Mufile, FILE *Sig2file,
		int *n_snps, int *n_genes,
		double* A, double* B, double* C, double* P,
		double* Mu, double* Sig2);

/*
void c_qtl_main_parallel_sparse(double *gene, int *n_indivs, int *n_genes, double *snp,
		int *n_snps, int *n_iter, int *burn_in, int *n_sweep, double *outProbs, int *nP, int *nmax,
		double *eps);
*/

void c_qtl_main_parallel_sparse(double *gene, int *n_indivs, int *n_genes, double *snp,
		int *n_snps, int *n_iter, int *burn_in, int *n_sweep, double *outProbs, int *nP, int *nmax,
		double *eps)
{
	/* Structures to hold various parameters.
	 *
	 * When possible variable names are analogous to those found in the paper */


	/* copying the elements of snp, gene, to matricies. R is column-major
	 * and GSL is row-major, so we send to a matrix and then
	 * transpose it to get a GSL representation of the same matrix.
	 */
	R_CStackLimit=(uintptr_t)-1;

	gsl_matrix* Y_t = gsl_matrix_calloc(*n_genes, *n_indivs);
	gsl_matrix* Y = gsl_matrix_calloc(*n_indivs, *n_genes); // worth it

	Y_t->data = gene;
	gsl_matrix_transpose_memcpy(Y, Y_t);
	gsl_matrix_free(Y_t);

	gsl_matrix* X_t = gsl_matrix_calloc(*n_snps, *n_indivs);
	gsl_matrix* X = gsl_matrix_calloc(*n_indivs, *n_snps);

	X_t->data = snp;
	gsl_matrix_transpose_memcpy(X, X_t);
	gsl_matrix_free(X_t);


	FILE *Cfile, *Pfile, *Afile, *Bfile, *Mufile, *Sig2file;
	Cfile = fopen("cfile.txt", "w");
	Afile = fopen("afile.txt", "w");
	Bfile = fopen("bfile.txt", "w");
	Pfile = fopen("pfile.txt", "w");
	Mufile = fopen("Mufile.txt", "w");
	Sig2file = fopen("Sig2file.txt", "w");

	int iter, i, j, g, th_id = 0, kill = 0;

	// statically allocated arrays for Adaptive Rejection Sampling
	ARS_workspace workspace;
	*nmax = (*nmax < NMAX) ? *nmax : NMAX;

	// initialize linked list for beta parameters
	m_el **Beta;
	Beta = malloc(*n_genes*sizeof(m_el*));
	if(Beta == NULL)
	{
		Rprintf("Memory allocation failed, exiting.\n");
		return;
	}
	for(g = 0; g < *n_genes; g++)
	{
		Beta[g] = malloc(sizeof(m_el));
		if(Beta[g] == NULL)
		{
			Rprintf("Memory allocation failed, exiting.\n");
			return;
		}
		Beta[g]->next = NULL;
	}

	double **W, **xA, **xB;
	int **Gamma, **ProbSum;

	xA = (double**) malloc(*n_snps*sizeof(double*));
	xB = (double**) malloc(*n_snps*sizeof(double*));
	W = (double**) malloc(*n_snps*sizeof(double*));
	Gamma = (int**) malloc(*n_snps*sizeof(int*));
	ProbSum = (int**) malloc(*n_snps*sizeof(int*));
	if((W == NULL) || Gamma == NULL || ProbSum == NULL)
	{
		Rprintf("Memory allocation failed, exiting.\n");
		return;
	}

	for(j = 0; j < *n_snps; j++)
	{
		W[j] = (double*) malloc(*n_genes*sizeof(double));
		xA[j] = (double*) malloc(*nmax*sizeof(double));
		xB[j] = (double*) malloc(*nmax*sizeof(double));
		Gamma[j] = (int*) malloc(*n_genes*sizeof(int));
		ProbSum[j] = (int*) malloc(*n_genes*sizeof(int));
		if((W[j] == NULL) || Gamma[j] == NULL || ProbSum[j] == NULL)
		{
			Rprintf("Memory allocation failed, exiting.\n");
			return;
		}
	}


	double *A, *B, *C, *P, *Mu, *Sig2, *expr_means, *expr_vars, *alpha2_beta;
	int *test;

	A = (double*)malloc(*n_snps*sizeof(double));
	B = (double*)malloc(*n_snps*sizeof(double));;
	P = (double*)malloc(*n_snps*sizeof(double));
	C = (double*)malloc(*n_genes*sizeof(double));
	Mu = (double*)malloc(*n_genes*sizeof(double));
	Sig2 = (double*)malloc(*n_genes*sizeof(double));
	test = (int*)malloc(*n_snps*sizeof(int));

	//empirical mean and variance of gene expression
	expr_means = (double*)malloc(*n_genes*sizeof(double));
	expr_vars = (double*)malloc(*n_genes*sizeof(double));

	// holds (X_j)^T (X_j) for each snp
	alpha2_beta = (double*)malloc(*n_snps*sizeof(double));

	if(A == NULL || B == NULL || P == NULL || C == NULL || Mu == NULL || Sig2 == NULL ||
			expr_means == NULL || expr_vars == NULL || alpha2_beta == NULL)
	{
		Rprintf("Memory allocation failed, exiting.\n");
		return;
	}

	// cheesy but effective
	double lambdaa = 0.0;
	double a0 = 0.0;
	double lambdab = 0.0;
	double b0 = 0.0;
	double tau0 = 0.0;

	double* lambda_a = &lambdaa;
	double* a_0 = &a0;
	double* lambda_b = &lambdab;
	double* b_0 = &b0;
	double* tau_0 = &tau0;

	Rprintf("%d processors running \n", *nP);
	GetRNGstate();
	unsigned long seed[6];
	for(j = 0; j < 6; j++)
	{
		seed[j] = (unsigned long)(runif(0.0,1.0)*4294944443);
	}
	PutRNGstate();

	j = RngStream_SetPackageSeed(seed);
	if(j != 0)
	{
		Rprintf("Setting Seed failed \n");
	}

	// initialize the parallel PRNGS
	RngStream rngs[*nP];
	for(i = 0; i < *nP; i++)
	{
		rngs[i] = RngStream_CreateStream("");
		//RngStream_IncreasedPrecis(rngs[i], 1);
	}

	//set hyperparameters for prior distributions
	set_prior(lambda_a, lambda_b, a_0, b_0, tau_0, expr_means, expr_vars, alpha2_beta, X, Y, rngs[0]);

	//initialize all parameters to begin MCMC
	initialize_parms(Beta, Gamma, W, ProbSum, xA, xB,
			A, B, C, P, Mu, Sig2,
			expr_means, expr_vars, alpha2_beta,
			lambda_a, a_0, lambda_b, b_0, tau_0,
			n_snps, n_genes, n_indivs, nmax, rngs[0]);

	gsl_vector_view Y_g;
	gsl_vector* one = gsl_vector_calloc(*n_indivs);
	gsl_vector_set_all(one, 1.0);

	for(iter=0; iter <= (*burn_in+(*n_sweep)*(*n_iter)); iter++)
	{
		if(iter% 1000 == 0)
		{
			Rprintf("******** iter %d *********\n",1+iter);
		}
		//if the user wants to interrupt computation (^C)
		R_CheckUserInterrupt();

		#pragma omp parallel private(th_id, Y_g, workspace) num_threads(*nP)
		{
			th_id = omp_get_thread_num();

			#pragma omp for nowait
			for(g = 0; g < *n_genes; g++)
			{
				Y_g = gsl_matrix_column(Y, g);
				update_gene_g(Beta[g], Gamma, W, X, &Y_g.vector, A, B, &(C[g]), P, &(Mu[g]), &(Sig2[g]), expr_means, expr_vars, alpha2_beta,
						lambda_a, a_0, lambda_b, b_0, tau_0, n_snps, n_genes, n_indivs, g, one, rngs[th_id]);
			}

			#pragma omp for nowait
			for(j = 0; j < *n_snps; j++)
			{
				update_pos_j(P, A, B, C, W, Gamma,
						j, a_0, b_0, lambda_a, lambda_b,
						n_snps, n_genes, n_indivs, rngs[th_id], *nmax,
						xA[j], xB[j], &workspace, *eps);
			}
		}

		if((iter > (*burn_in)) & ((iter-(*burn_in))%(*n_sweep) == 0))
		{
			update_prob_include(n_snps, n_genes, Gamma, ProbSum);
			store_mcmc_output(Afile, Bfile,
					Cfile, Pfile, Mufile, Sig2file,
					n_snps, n_genes,
					A, B, C, P, Mu, Sig2);
		}
	}
	// outputs estimate of P(gamma[j][g] = 1 | data)
	store_prob_include(n_iter, n_snps, n_genes, ProbSum, outProbs);


	gsl_matrix_free(X);
	gsl_matrix_free(Y);

	for(g = 0; g < *n_genes; g++)
	{
		SV_free(Beta[g]);
	}
	free(Beta);
	gsl_vector_free(one);

	free(Sig2);
	free(Mu);
	free(A);
	free(B);
	free(C);
	free(P);
	free(alpha2_beta);
	free(expr_means);
	free(expr_vars);
	for(j = 0; j < *n_snps; j++)
	{
		free(W[j]);
		free(Gamma[j]);
		free(ProbSum[j]);
		free(xA[j]);
		free(xB[j]);
	}

	free(W);
	free(Gamma);
	free(ProbSum);
	free(xA);
	free(xB);

	for(i = 0; i < *nP; i++)
	{
		RngStream_DeleteStream(rngs[i]);
	}

	fclose(Afile);
	fclose(Bfile);
	fclose(Cfile);
	fclose(Pfile);
	fclose(Mufile);
	fclose(Sig2file);

	return;
}


void update_gene_g(m_el* beta_g, int** Gamma, double** W, gsl_matrix* X, gsl_vector* Y_g,
		double* A, double* B, double* C_g, double* P, double* Mu_g, double* Sig2_g,
		double* expr_means, double* expr_vars, double* alpha2_beta,
		double* lambda_a, double* a_0, double* lambda_b, double* b_0, double* tau_0,
		int* n_snps, int* n_genes, int* n_indivs, int g, gsl_vector* one, RngStream rng)
{
	gsl_vector_view X_j;
	gsl_vector* Y_minus_mu_g = gsl_vector_calloc(*n_indivs);

	gsl_blas_dcopy(Y_g, Y_minus_mu_g);
	gsl_vector_add_constant(Y_minus_mu_g, -1.0*(*Mu_g));

	gsl_vector* Xbeta_sum_g = gsl_vector_calloc(*n_indivs);
	gsl_vector_set_all(Xbeta_sum_g, 0.0);

	gsl_vector* Resid_minus_j = gsl_vector_calloc(*n_indivs);
	int i, j, cur;
	double S_j, v1, v2, p1, w_jg, temp;

	// X %*% B_g   (fitted expression values for gene g), using the sparse representation of beta
	SV_gsl_dmvpy(X, beta_g, Xbeta_sum_g->data, *n_indivs);
	v1 = 1.0/sqrt(1.0 + (*C_g));

	for(j = 0; j < *n_snps; j++)
	{
		X_j = gsl_matrix_column(X, j);
		// take away the j'th component from the sum, which is non-zero only if gamma[j][g] is nonzero,
		// hence the check.
		if(Gamma[j][g] == 1)
		{
			gsl_blas_daxpy(-1.0*SV_get(beta_g, j), &X_j.vector, Xbeta_sum_g);
		}

		gsl_blas_dcopy(Y_minus_mu_g, Resid_minus_j);
		gsl_blas_daxpy(-1.0, Xbeta_sum_g, Resid_minus_j);
		gsl_blas_ddot(&X_j.vector, Resid_minus_j, &S_j);

		v2 = alpha2_beta[j]*(*C_g)/(1.0 + (*C_g));

		p1=v1*exp(0.5*v2*pow(S_j,2)/(*Sig2_g));
		w_jg = W[j][g];

		// update gamma
		cur = (int)(RngStream_RandU01(rng) <= w_jg*p1/(1 - w_jg + w_jg*p1));

		// update betas
		if(Gamma[j][g] == 1 && cur == 0)
		{
			// removing element j is the same as setting it to zero in a sparse representation
			SV_remove_el(beta_g, j);
			Gamma[j][g] = cur;
			// no need to re-update the sum for next iteration, because it's zero now and
			// we already removed it above
		}

		if(cur == 1)
		{
			// now the element beta[j][g] is nonzero, so we sample it and add it to the list
			temp = S_j*alpha2_beta[j]*(*C_g)/(1.0 + (*C_g)) + sqrt(v2*(*Sig2_g))*RngStream_N01(rng);

			// add element j to the list (if it's already there that's ok, it will just change the value then)
			SV_add_el(beta_g, j, temp);
			Gamma[j][g] = cur;

			// re-update the sum for next iteration
			gsl_blas_daxpy(SV_get(beta_g, j), &X_j.vector, Xbeta_sum_g);
		}
	}

	//update C_g
	double G_g = 0.0;
	v1=(1.0)/2.0;

	for(j = 0;j < *n_snps; j++)
	{
		if(Gamma[j][g] == 1)
		{
			v1 += (double)(Gamma[j][g])/2.0;
			G_g += gsl_pow_2(SV_get(beta_g, j))/alpha2_beta[j];
		}
	}

	v2 = (double)*n_indivs/2.0 + G_g/(2*(*Sig2_g));
	*C_g = v2/RngStream_GA1(v1, rng);

	//update Sig2_g
	G_g = G_g/(*C_g);
	G_g += gsl_pow_2((*Mu_g) - expr_means[g])/(*tau_0);

	for(i = 0; i < *n_indivs; i++)
	{
		G_g += gsl_pow_2(gsl_vector_get(Y_minus_mu_g, i) - gsl_vector_get(Xbeta_sum_g, i));
	}
	v1 += (double)*n_indivs/2.0;
	G_g = G_g/2.0;
	*Sig2_g = G_g/RngStream_GA1(v1, rng);

	//update Mu_g

	double sum_y, sum_G_g, si;
	gsl_blas_ddot(Y_g, one, &sum_y);
	gsl_blas_ddot(Xbeta_sum_g, one, &sum_G_g);
	si = sum_y - sum_G_g;

	v1=(si + expr_means[g]/(*tau_0))
			/(1.0/(*tau_0) +(double)(*n_indivs));
	v2=(*Sig2_g)/((double)(*n_indivs) + 1/(*tau_0));

	*Mu_g = v1+sqrt(v2)*RngStream_N01(rng);

	gsl_vector_free(Y_minus_mu_g);
	gsl_vector_free(Xbeta_sum_g);
	gsl_vector_free(Resid_minus_j);

	return;
}


void store_mcmc_output(FILE *Afile, FILE *Bfile,
		FILE *Cfile, FILE *Pfile, FILE *Mufile, FILE *Sig2file,
		int *n_snps, int *n_genes,
		double* A, double* B, double* C, double* P,
		double* Mu, double* Sig2)
{
	int g,j;
	for(g = 0; g < *n_genes; g++)
	{
		fprintf(Mufile, "%f\t", Mu[g]);
		fprintf(Sig2file, "%f\t", Sig2[g]);
		fprintf(Cfile, "%f\t", C[g]);
	}
	for(j = 0; j < *n_snps; j++)
	{
		fprintf(Afile, "%f\t", A[j]);
		fprintf(Bfile, "%f\t", B[j]);
		fprintf(Pfile, "%f\t", P[j]);
	}
	fprintf(Afile, "\n");
	fprintf(Bfile, "\n");
	fprintf(Cfile, "\n");
	fprintf(Pfile, "\n");
	fprintf(Mufile, "\n");
	fprintf(Sig2file, "\n");
	return;
}
