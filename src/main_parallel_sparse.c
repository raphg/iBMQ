/*
 * main_parallel_sparse_constC.c
 *
 *  Created on: Oct 13, 2011
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

void update_gene_g(ptr_m_el beta_g, int** Gamma, double** W_Logit,
		int** W_Ind, gsl_matrix* X, const gsl_vector* Y_g,
		double* C_g, double* Mu_g, double* Sig2_g,
		const double* expr_means, const double* expr_vars, const double* alpha2_beta,
		int* n_snps, int* n_genes, int* n_indivs, int g, const gsl_vector* one, RngStream rng,
		ptr_memChunk ptr_chunk_g, int variable_C);

void store_mcmc_output(FILE *Afile, FILE *Bfile, FILE *Pfile, FILE *Mufile, FILE *Sig2file,
		FILE *Cfile,
		int *n_snps, int *n_genes, double* restrict A, double* restrict B, double* restrict P,
		double* restrict Mu, double* restrict Sig2, double* restrict C, int variable_C);

void c_qtl_main_parallel_sparse(double *gene, int *n_indivs, int *n_genes, double *snp,
		int *n_snps, int *n_iter, int *burn_in, int *n_sweep, double *outProbs, int *nP, int *nmax,
		double *eps, int *write_output, int *variable_C)
{
	R_CStackLimit = (uintptr_t)-1;

	int iter, i, j, g, th_id = 0;

	// initialize a memory pool for linked list elements of the sparse matrix;
	// it's way bigger than it needs to be, just to be on the safe side,
	// and to facilitate a more simple implementation.

	memPool pool;
	ptr_memPool ptr_pool = &pool;
	initializePool(*n_genes, (*n_snps+2), ptr_pool);

	/* copying the elements of snp, gene, to matricies. R is column-major
	 * and GSL is row-major, so we send to a matrix and then
	 * transpose it to get a GSL representation of the same matrix.
	 */

	gsl_matrix* Y_t = gsl_matrix_calloc(*n_genes, *n_indivs);
	gsl_matrix* Y = gsl_matrix_calloc(*n_indivs, *n_genes);

	Y_t->data = gene;
	gsl_matrix_transpose_memcpy(Y, Y_t);
	gsl_matrix_free(Y_t);

	gsl_matrix* X_t = gsl_matrix_calloc(*n_snps, *n_indivs);
	gsl_matrix* X = gsl_matrix_calloc(*n_indivs, *n_snps);

	X_t->data = snp;
	gsl_matrix_transpose_memcpy(X, X_t);
	gsl_matrix_free(X_t);

	FILE *Pfile, *Afile, *Bfile, *Mufile, *Sig2file, *Cfile;
	// open files to store MCMC output
	if(*write_output != 0)
	{
		Afile = fopen("afile.txt", "w");
		Bfile = fopen("bfile.txt", "w");
		Pfile = fopen("pfile.txt", "w");
		Mufile = fopen("Mufile.txt", "w");
		Sig2file = fopen("Sig2file.txt", "w");

		if(*variable_C)
		{
			Cfile = fopen("cfile.txt", "w");
		}
	}
	//statically allocated workspace for Adaptive Rejection Sampling
	ARS_workspace workspace;
	*nmax = (*nmax < NMAX) ? *nmax : NMAX;

	// initialize linked list for beta parameters
	ptr_m_el *Beta;
	Beta = (ptr_m_el*) R_alloc(*n_genes, sizeof(ptr_m_el));
	for(g = 0; g < *n_genes; g++)
	{
		Beta[g] = (ptr_m_el) R_alloc(1, sizeof(m_el));
		Beta[g]->next = NULL;
	}

	// initialize double indexed arrays
	// W_logit holds logit(omega) values
	double **W_Logit, **xA, **xB;
	int **Gamma, **ProbSum, **W_Ind;

	W_Logit = (double**) R_alloc(*n_snps, sizeof(double*));
	W_Ind = (int**) R_alloc(*n_snps, sizeof(int*));
	xA = (double**) R_alloc(*n_snps, sizeof(double*));
	xB = (double**) R_alloc(*n_snps, sizeof(double*));
	Gamma = (int**) R_alloc(*n_snps, sizeof(int*));
	ProbSum = (int**) R_alloc(*n_snps, sizeof(int*));

	for(j = 0; j < *n_snps; j++)
	{
		W_Logit[j] = (double*) R_alloc(*n_genes, sizeof(double));
		W_Ind[j] = (int*) R_alloc(*n_genes, sizeof(int));
		xA[j] = (double*) R_alloc(*nmax, sizeof(double));
		xB[j] = (double*) R_alloc(*nmax, sizeof(double));
		Gamma[j] = (int*) R_alloc(*n_genes, sizeof(int));
		ProbSum[j] = (int*) R_alloc(*n_genes, sizeof(int));
	}

	// initialize single indexed arrays
	double *A, *B, *C, *P, *Mu, *Sig2, *expr_means, *expr_vars, *alpha2_beta;

	A = (double*) R_alloc (*n_snps, sizeof(double));
	B = (double*) R_alloc(*n_snps, sizeof(double));;
	P = (double*) R_alloc(*n_snps, sizeof(double));
	C = (double*) R_alloc(*n_genes, sizeof(double));
	Mu = (double*) R_alloc(*n_genes, sizeof(double));
	Sig2 = (double*) R_alloc(*n_genes, sizeof(double));

	//empirical mean and variance of gene expression
	expr_means = (double*) R_alloc(*n_genes, sizeof(double));
	expr_vars = (double*) R_alloc(*n_genes, sizeof(double));

	// holds (X_j)^T (X_j) for each snp j
	alpha2_beta = (double*) R_alloc(*n_snps, sizeof(double));

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

	// initialize random number streams
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

	RngStream rngs[*nP];
	for(i = 0; i < *nP; i++)
	{
		rngs[i] = RngStream_CreateStream("");
		//RngStream_IncreasedPrecis(rngs[i], 1);
	}

	//set hyperparameters for prior distributions, compute
	// auxiallary quantities
	set_prior(lambda_a, lambda_b, a_0, b_0, tau_0, expr_means,
			expr_vars, alpha2_beta, X, Y, rngs[0]);

	//initialize all parameters to begin MCMC
	initialize_parms(Beta, ptr_pool, Gamma,  W_Logit, W_Ind, ProbSum, xA, xB,
			A, B, C, P, Mu, Sig2,
			expr_means, expr_vars, alpha2_beta,
			lambda_a, a_0, lambda_b, b_0, tau_0,
			n_snps, n_genes, n_indivs, nmax, rngs[0]);

	// working vectors
	gsl_vector_view Y_g;
	gsl_vector* one = gsl_vector_calloc(*n_indivs);
	gsl_vector_set_all(one, 1.0);
	ptr_memChunk chunk_g;

	Rprintf("initialization successful\n");
	for(iter=0; iter <= (*burn_in+(*n_sweep)*(*n_iter)); iter++)
	{
		if(iter% 1000 == 0)
		{
			Rprintf("******** iter %d *********\n",1+iter);
		}
		//if the user wants to interrupt computation (^C)
		R_CheckUserInterrupt();
		// update genes and positions
		#pragma omp parallel private(th_id, Y_g, workspace, chunk_g) num_threads(*nP)
		{
			th_id = omp_get_thread_num();

			#pragma omp for
			for(g = 0; g < *n_genes; g++)
			{
				Y_g = gsl_matrix_column(Y, g);
				chunk_g = ptr_pool->array_head[g];
				update_gene_g(Beta[g], Gamma, W_Logit, W_Ind, X,
						&Y_g.vector, &(C[g]), &(Mu[g]), &(Sig2[g]),
						expr_means, expr_vars, alpha2_beta,
						n_snps, n_genes, n_indivs, g, one, rngs[th_id], chunk_g,
						*variable_C);
			}

			#pragma omp for
			for(j = 0; j < *n_snps; j++)
			{
				update_pos_j(P, A, B, W_Logit, W_Ind, Gamma,
						j, a_0, b_0, lambda_a, lambda_b,
						n_genes, rngs[th_id], *nmax,
						xA[j], xB[j], &workspace, *eps);
			}
		}

		if((iter > (*burn_in)) & ((iter-(*burn_in))%(*n_sweep) == 0))
		{
			update_prob_include(n_snps, n_genes, Gamma, ProbSum);
			if(*write_output != 0)
			{
				store_mcmc_output(Afile, Bfile, Pfile, Mufile, Sig2file, Cfile,
						n_snps, n_genes, A, B, P, Mu, Sig2, C, *variable_C);
			}

		}
	}
	// outputs estimate of P(gamma[j][g] = 1 | data), and free resources
	store_prob_include(n_iter, n_snps, n_genes, ProbSum, outProbs);

	gsl_matrix_free(X);
	gsl_matrix_free(Y);
	gsl_vector_free(one);

	for(i = 0; i < *nP; i++)
	{
		RngStream_DeleteStream(rngs[i]);
	}

	if(*write_output != 0)
	{
		fclose(Afile);
		fclose(Bfile);
		fclose(Pfile);
		fclose(Mufile);
		fclose(Sig2file);
		if(*variable_C)
		{
			fclose(Cfile);
		}
	}
	return;
}

void update_gene_g(ptr_m_el beta_g, int** Gamma,  double** W_Logit,
		int** W_Ind, gsl_matrix* X, const gsl_vector* Y_g,
		double* C_g, double* Mu_g, double* Sig2_g,
		const double* expr_means, const double* expr_vars, const double* alpha2_beta,
		int* n_snps, int* n_genes, int* n_indivs, int g, const gsl_vector* one, RngStream rng,
		ptr_memChunk ptr_chunk_g, int variable_C)
{
	gsl_vector* Y_minus_mu_g = gsl_vector_calloc(*n_indivs);

	gsl_blas_dcopy(Y_g, Y_minus_mu_g);
	gsl_vector_add_constant(Y_minus_mu_g, -1.0*(*Mu_g));

	gsl_vector* Xbeta_sum_g = gsl_vector_calloc(*n_indivs);
	gsl_vector_set_all(Xbeta_sum_g, 0.0);

	gsl_vector* Resid_minus_j = gsl_vector_calloc(*n_indivs);
	int i, j, cur;
	double S_j, v1, v2, p1, w_jg, temp, logodds;
	double Cg = *C_g;
	double Sg = *Sig2_g;

	// X %*% B_g   (fitted expression values for gene g), using the sparse representation of beta
	SV_gsl_dmvpy(X, beta_g, Xbeta_sum_g->data, *n_indivs);

	v1 = -.5*log1p(Cg);

	for(j = 0; j < *n_snps; j++)
	{
		if(W_Ind[j][g] == 1)
		{
			gsl_vector_view X_j = gsl_matrix_column(X, j);
			// take away the j'th component from the sum, which is non-zero only if gamma[j][g] is nonzero
			if(Gamma[j][g] == 1)
			{
				gsl_blas_daxpy(-1.0*SV_get(beta_g, j), &X_j.vector, Xbeta_sum_g);
			}

			gsl_blas_dcopy(Y_minus_mu_g, Resid_minus_j);
			gsl_blas_daxpy(-1.0, Xbeta_sum_g, Resid_minus_j);

			//calculates sum of x values times residuals, without j'th value
			gsl_blas_ddot(&X_j.vector, Resid_minus_j, &S_j);

			v2 = alpha2_beta[j]*(Cg)/(1.0 + (Cg));

			logodds = v1 + 0.5*v2*pow(S_j,2)/(Sg) + W_Logit[j][g];

			// update gammas
			cur = (int)(logit(RngStream_RandU01(rng)) <= logodds);

			// update betas
			if(Gamma[j][g] == 1 && cur == 0)
			{
				// removing element j is the same as setting it to zero in a sparse representation
				SV_remove_el(beta_g, j, ptr_chunk_g);
				Gamma[j][g] = cur;
				// no need to re-update the sum for next iteration, because it's zero now and
				// we already removed it above
			}

			if(cur == 1)
			{
				// now the element beta[j][g] is nonzero, so we simulate it and add it to the list
				temp = S_j*alpha2_beta[j]*(Cg)/(1.0 + (Cg)) + sqrt(v2*(Sg))*RngStream_N01(rng);

				// add element j to the list (if it's already there that's ok, it will just change the value then)
				SV_add_el(beta_g, j, temp, ptr_chunk_g);
				Gamma[j][g] = cur;

				// re-update the sum for next iteration
				gsl_blas_daxpy(SV_get(beta_g, j), &X_j.vector, Xbeta_sum_g);
			}
		}
	}

	//update C_g
	if(variable_C)
	{
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

		v2 = (double)*n_indivs/2.0 + G_g/(2.0*(Sg));
		Cg = v2/RngStream_GA1(v1, rng);
		*C_g = Cg;
	}

	//update Sig2_g
	double G_g = 0.0;
	v1=(1.0)/2.0;

	for(j = 0;j < *n_snps; j++)
	{
		v1+=(double)(Gamma[j][g])/2.;
		if(Gamma[j][g] == 1)
		{
			G_g += gsl_pow_2(SV_get(beta_g, j))/alpha2_beta[j];
		}
	}
	G_g = G_g/(Cg); // (C_g is fixed equal to number of subjects

	for(i = 0; i < *n_indivs; i++)
	{
		G_g += gsl_pow_2(gsl_vector_get(Y_minus_mu_g, i) - gsl_vector_get(Xbeta_sum_g, i));
	}
	v1 += (double)*n_indivs/2.0;
	G_g = G_g/2.0;

	Sg = G_g/RngStream_GA1(v1, rng);
	*Sig2_g = Sg;

	//update Mu_g

	double sum_y, sum_G_g, si;

	gsl_blas_ddot(Y_g, one, &sum_y);
	gsl_blas_ddot(Xbeta_sum_g, one, &sum_G_g);
	si = sum_y - sum_G_g;

	v1 = (si / Sg + expr_means[g] / expr_vars[g])/
			((double)*n_indivs / (Sg) + 1.0 / expr_vars[g]);
	v2 = 1.0 / ((double)*n_indivs / (Sg) + 1.0 / expr_vars[g]);

	*Mu_g = v1 + sqrt(v2) * RngStream_N01(rng);

	gsl_vector_free(Y_minus_mu_g);
	gsl_vector_free(Xbeta_sum_g);
	gsl_vector_free(Resid_minus_j);

	return;
}

void store_mcmc_output(FILE *Afile, FILE *Bfile, FILE *Pfile, FILE *Mufile, FILE *Sig2file,
		FILE *Cfile,
		int *n_snps, int *n_genes, double* restrict A, double* restrict B, double* restrict P,
		double* restrict Mu, double* restrict Sig2, double* restrict C, int variable_C)
{
	int g,j;

	for(g = 0; g < *n_genes; g++)
	{
		const double mg = Mu[g];
		const double sg = Sig2[g];
		const double cg = C[g];

		fprintf(Mufile, "%f\t", mg);
		fprintf(Sig2file, "%f\t", sg);
		if(variable_C)
		{
			fprintf(Cfile, "%f\t", cg);
		}
	}
	for(j = 0; j < *n_snps; j++)
	{
		const double aj = A[j];
		const double bj = B[j];
		const double pj = P[j];

		fprintf(Afile, "%f\t", aj);
		fprintf(Bfile, "%f\t", bj);
		fprintf(Pfile, "%f\t", pj);
	}

	fprintf(Afile, "\n");
	fprintf(Bfile, "\n");
	fprintf(Pfile, "\n");
	fprintf(Mufile, "\n");
	fprintf(Sig2file, "\n");
	if(variable_C)
	{
		fprintf(Cfile, "\n");
	}
	return;
}
