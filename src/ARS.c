/*
 * ARS.c
 *
 *  Created on: Jan 22, 2012
 *      Author: Gregory Imholte
 */

#include "ARS.h"
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

// function to test the sampler, not used in main function.
void sample_one(double* start_val, int* n_max, double* S_j, double* c, double* R, double* out,
		double* eps)
{
	double *x;
	double sample;
	int num_x = 2;

	x = (double*)malloc(*n_max * sizeof(double));

	ARS_workspace space;
	RngStream rng;
	rng = RngStream_CreateStream("");

	x[0] = start_val[0];
	x[1] = start_val[1];
	sample = sample_conditional(x, &num_x, *n_max, *S_j, *c, *R, &space, rng, *eps);
	out[0] = sample;
	out[1] = x[0];
	out[2] = x[1];

	free(x);
	RngStream_DeleteStream(rng);

	return;
}

/* Reference:
 * Gilks and Wild (1992)
 * "Adaptive Rejection Sampling for Gibbs Sampling"
 *
 * Adaptive algorithm to sample from a target density
 * of the full conditionals of a_j and b_j
 * (note, this algorithm can be adapted to any log-concave
 * density on (0, infinity), and with a little more work,
 * a density on any connected support of R)
 *
 * to make this more general, one could pass pointers to
 * target functions h, h' instead of calling them outside the function
 *
 * for arbitrary densities there may be overflow or underflow issues,
 * but for the designed usage, values are stable.
 *
 * x: starting lattice of points
 * num_x: number of starting points
 * nmax: maximum number of points in lattice
 *
 * S_j, c, R: parameters from conditional density
 *
 * ARS_workspace is a struct with elements
 * hwv: working vector for the h(x) values
 * hpwv: working vector for h'(x) values
 * scum: working vector for cumulative integral values of the hull section
 * scum_norm: working vector cumulative density of hull sections
 * z: intersection of hull tangents
 * s: working vector for calculating scum
 *
 *
 * rng: random number generator
 * eps: double machine precision value
 */


// returns -1.0 if an error is detected.
double sample_conditional(double* x, int* num_x, int nmax, double S_j, double c, double R, ARS_workspace *space,
		RngStream rng, double eps)
{
	int i, u_section, l_section;
	int accept = 0;
	int update = 0;
	double x_samp, W, l, u, hnew, hpnew, zlower, zupper, V;
	double huzmax, test1, test2, tmp0, tmp1;

	// intialize the lattice
	for(i = 0; i < *num_x; i++)
	{
		space->hwv[i] = h(x[i], S_j, c, R);
		space->hpwv[i] = h_prime(x[i], S_j, c, R);
	}

	// need last hpwv  < 0, move x_last farther out until this is so. This only works with
	// strict concavity!

	while(space->hpwv[*num_x-1] >= 0)
	{
		x[*num_x - 1] = x[*num_x - 1] + 2.0;
		space->hpwv[*num_x - 1] = h_prime(x[*num_x - 1], S_j, c, R);

	}
	space->hwv[*num_x - 1] = h(x[*num_x - 1], S_j, c, R);

	//now calculate intersection of lines, z
	// i-th element corresponds to intersection to right of x[i]
	for(i = 0; i < (*num_x - 1); i++)
	{
		if(i == 0)
		{
			// initialize the search largest y value in outer hull
			huzmax = space->hwv[0] - x[0]*space->hpwv[0];
		}
		// tangents are nearly parallel, use midpoint
		// approximation to avoid instability
		if((space->hpwv[i] - space->hpwv[i+1]) <= eps)
		{
			space->z[i] = (x[i+1] + x[i])/2.0;
			test1 = huzmax;
			test2 = (space->hwv[i] + (space->z[i] - x[i])*(space->hpwv[i]));
			huzmax = (test1 > test2) ? test1 : test2;
		}
		else
		{
			space->z[i] = (space->hwv[i+1] - space->hwv[i] - x[i+1]*(space->hpwv[i+1]) +
					x[i]*(space->hpwv[i]))/
					(space->hpwv[i] - space->hpwv[i+1]);
			test1 = huzmax;
			test2 = (space->hwv[i] + (space->z[i] - x[i])*(space->hpwv[i]));
			huzmax = (test1 > test2) ? test1 : test2;
		}

	}

	//huzmax holds the maximum value of the hull on the log scale

	//Rprintf("computed z values\n");
	// this initializes the hull, i.e. creates the cumulative distribution function of the
	// outer hull.
	initialize_hull(x, space, num_x, &huzmax);

	//Rprintf("Hull initialized\n");
	while(accept == 0)
	{
		V = RngStream_RandU01(rng);

		x_samp = sample_hull(x, space, num_x, &u_section, V, &huzmax);
		// check for sampling trouble
		if(isnan(x_samp) || isinf(x_samp) || (x_samp <= 0))
		{
			for(i = 0; i < *num_x; i++)
			{
				Rprintf("i = %d, x = %.3lf, scum = %.3lf, scum_norm = %.3lf, hpwv = %.3lf, hwv = %.3lf \n",
						i, x[i], space->scum[i], space->scum_norm[i], space->hpwv[i], space->hwv[i]);
			}
			for(i = 0; i < *num_x - 1; i++)
			{
				Rprintf("z_%d = %lf\n", i, space->z[i]);
			}
			Rprintf("R = %lf, S_j = %lf, c = %lf\n", R, S_j, c);
			Rprintf("x_samp = %lf \n", x_samp);
			return(-1.0);
		}
		W = RngStream_RandU01(rng);

		// u_section is the j corresponding to x[j]
		u = space->hwv[u_section] + (x_samp - x[u_section])*(space->hpwv[u_section]);

		// check sandwich acceptance
		l_section = 0;
		while((l_section < *num_x) & (x[l_section] < x_samp))
		{
			l_section++;
		}

		// x_samp is between x[l_section - 1] and x[l_section]
		if((l_section != 0) & (l_section != *num_x))
		{
			int j = l_section;
			// we're in the interior of the lattice

			// calculate value in lower hull.
			l =	((x[j] - x_samp)*(space->hwv[j-1]) + (x_samp - x[j-1])*(space->hwv[j]))/
					(x[j] - x[j-1]);

			//squeeze test.
			if(log(W) <= (l - u))
			{
				tmp0 = sample_hull(x, space, num_x, &u_section, .15, &huzmax);
				tmp1 = sample_hull(x, space, num_x, &u_section, .85, &huzmax);
				x[0]=tmp0;
				x[1]=tmp1;
				if(isnan(x[0]) || isnan(x[1]))
				{
					for(i = 0; i < *num_x; i++)
					{
						Rprintf("i = %d, x = %.3lf, scum = %.3lf, scum_norm = %.3lf, hpwv = %.3lf, hwv = %.3lf \n",
								i, x[i], space->scum[i], space->scum_norm[i], space->hpwv[i], space->hwv[i]);
					}
					for(i = 0; i < *num_x - 1; i++)
					{
						Rprintf("z_%d = %lf\n", i, space->z[i]);
					}
					Rprintf("R = %lf, S_j = %lf, c = %lf\n", R, S_j, c);
					Rprintf("x_samp = %lf \n", x_samp);

				}
				return(x_samp);
			}
			// failed squeeze test, need to evaluate h at x_samp
			// and compare directly to outer hull.
			else
			{
				hnew = h(x_samp, S_j, c, R);
				if(log(W) <= (hnew - u))
				{
					tmp0 = sample_hull(x, space, num_x, &u_section, .15, &huzmax);
					tmp1 = sample_hull(x, space, num_x, &u_section, .85, &huzmax);
					x[0]=tmp0;
					x[1]=tmp1;
					if(isnan(x[0]) || isnan(x[1]))
					{
						for(i = 0; i < *num_x; i++)
						{
							Rprintf("i = %d, x = %.3lf, scum = %.3lf, scum_norm = %.3lf, hpwv = %.3lf, hwv = %.3lf \n",
									i, x[i], space->scum[i], space->scum_norm[i], space->hpwv[i], space->hwv[i]);
						}
						for(i = 0; i < *num_x - 1; i++)
						{
							Rprintf("z_%d = %lf\n", i, space->z[i]);
						}
						Rprintf("R = %lf, S_j = %lf, c = %lf\n", R, S_j, c);
						Rprintf("x_samp = %lf \n", x_samp);
					}
					return(x_samp);
				}
				else
				{
					update = update_hull(x, space, num_x, nmax,
							x_samp, hnew, l_section, R, S_j, c, &huzmax);
					// if update == 1 we added a new point to the lattice
					if(update == 1)
					{
						initialize_hull(x, space, num_x, &huzmax);
					}
					//Rprintf("huzmax = %lf\n", huzmax);
				}
			}
		}

		else
		{
			// we're on either the left or right edge or the
			// hull when sampling here
			hnew = h(x_samp, S_j, c, R);
			if(log(W) <= (hnew - u))
			{
				// sample accepted
				tmp0 = sample_hull(x, space, num_x, &u_section, .15, &huzmax);
				tmp1 = sample_hull(x, space, num_x, &u_section, .85, &huzmax);
				x[0]=tmp0;
				x[1]=tmp1;
				if(isnan(x[0]) || isnan(x[1]))
				{
					for(i = 0; i < *num_x; i++)
					{
						Rprintf("i = %d, x = %.3lf, scum = %.3lf, scum_norm = %.3lf, hpwv = %.3lf, hwv = %.3lf \n",
								i, x[i], space->scum[i], space->scum_norm[i], space->hpwv[i], space->hwv[i]);
					}
					for(i = 0; i < *num_x - 1; i++)
					{
						Rprintf("z_%d = %lf\n", i, space->z[i]);
					}
					Rprintf("R = %lf, S_j = %lf, c = %lf\n", R, S_j, c);
					Rprintf("x_samp = %lf \n", x_samp);
				}
				return(x_samp);
			}
			else
			{
				update = update_hull(x, space, num_x, nmax,
						x_samp, hnew, l_section, R, S_j, c, &huzmax);
				// if update == 1 we added a new point to the lattice
				// and we must re-normalize the outer hull
				if(update == 1)
				{
					initialize_hull(x, space, num_x, &huzmax);
				}
			}
		}
	}
	//value sampled successfully
}

// insert a new abcissae in the hull, update the Z and
int update_hull(double *x, ARS_workspace *space, int *num_x, int nmax,
		double xnew, double hnew, int l_section,
		double R, double S_j, double c, double *huzmax)
{
	if(nmax == *num_x)
	{
		// not allowed to add more breakpoints;
		// do not update hull and return to sampler
		return(0);
	}

	double hpnew, test2, test1;
	int i;

	hpnew = h_prime(xnew, S_j, c, R);
	// we are past l_section, so x_samp goes between
	// x[l_section - 1] and x[l_section]
	if(l_section == 0)
	{
		//new sample fell in leftmost chamber
		*num_x = *num_x + 1;
		for(i = *num_x - 1; i > 0; i--)
		{
			x[i] = x[i-1];
			space->hwv[i] = space->hwv[i-1];
			space->hpwv[i] = space->hpwv[i-1];
		}
		x[0] = xnew;
		space->hwv[0] = hnew;
		space->hpwv[0] = hpnew;

		for(i = *num_x - 2; i > 0; i--)
		{
			space->z[i] = space->z[i-1];
		}

		space->z[0] = (space->hwv[1] - space->hwv[0] - x[1]*(space->hpwv[1]) + x[0]*(space->hpwv[0]))/
				(space->hpwv[0] - space->hpwv[1]);
	}
	else if(l_section == *num_x)
	{
		// new sample goes on the "right-hand" end of the vector.
		x[*num_x] = xnew;
		space->hwv[*num_x] = hnew;
		space->hpwv[*num_x] = hpnew;
		space->z[*num_x - 1] = (space->hwv[*num_x] - space->hwv[*num_x-1] - x[*num_x]*(space->hpwv[*num_x]) +
				x[*num_x-1]*(space->hpwv[*num_x-1]))/
				(space->hpwv[*num_x-1] - space->hpwv[*num_x]);

		*num_x = *num_x + 1;
	}
	else
	{
		// new sample goes between x[l_section - 1] and x[l_section]
		// x[l_section - 1] and below remain untouched

		for(i = *num_x - 1; i > l_section - 1; i--)
		{
			x[i+1] = x[i];
			space->hwv[i+1] = space->hwv[i];
			space->hpwv[i+1] = space->hpwv[i];
		}
		x[l_section] = xnew;
		space->hwv[l_section] = hnew;
		space->hpwv[l_section] = hpnew;

		for(i = *num_x - 2; i > l_section - 1; i--)
		{
			space->z[i+1] = space->z[i];
			//hzwv[i+1] = hzwv[i];
		}
		for(i = l_section - 1; i < l_section + 1; i++)
		{
			space->z[i] = (space->hwv[i+1] - space->hwv[i] - x[i+1]*(space->hpwv[i+1]) + x[i]*(space->hpwv[i]))/
					(space->hpwv[i] - space->hpwv[i+1]);
		}
		*num_x = *num_x + 1;
	}

	// need to update the hull-max-value
	for(i = 0; i < *num_x - 1; i++)
	{
		if(i == 0)
		{
			// initialize the search largest y value in outer hull
			*huzmax = space->hwv[0] - x[0]*(space->hpwv[0]);
		}
		test1 = *huzmax;
		test2 = (space->hwv[i] + (space->z[i] - x[i])*(space->hpwv[i]));
		*huzmax = (test1 > test2) ? test1 : test2;
	}
	/*
	Rprintf("Hull updated: \n");
	for(i = 0; i < *num_x; i++)
	{
		Rprintf("i = %d, x = %.3lf, hpwv = %.3lf, hwv = %.3lf \n",
				i, x[i], hpwv[i], hwv[i]);
	}
	for(i = 0; i < *num_x - 1; i++)
	{
		Rprintf("z_%d = %lf\n", i, z[i]);
	}
	*/
	return(1);
}

// this function doubles as the quantile function of the outer hull
double sample_hull(double *x, ARS_workspace *space, int* num_x, int *section, double p, double *huzmax)
{
	int j = 0;
	double R, y;


	// find which "chunk" we're in. returns the index of the j such that
	// our sampled x is between z[j-1] and z[j] (z[-1] = 0, z[num_x] = infinity)
	while(p > space->scum_norm[j])
	{
		j++;
	}

	if(j == 0)
	{
		R = (p*(space->scum[*num_x - 1]))*(space->hpwv[j]) + exp(space->hwv[j] - x[j]*(space->hpwv[j]) - *huzmax);
		y = (log(R) + *huzmax - space->hwv[j])/(space->hpwv[j]) + x[j];
	}
	else
	{
		R = (p*(space->scum[*num_x - 1]) - (space->scum[j-1]))*(space->hpwv[j]) +
				exp(space->hwv[j] + (space->z[j-1] - x[j])*(space->hpwv[j]) - *huzmax);
		y = (log(R) + *huzmax - space->hwv[j])/(space->hpwv[j]) + x[j];
	}

	*section = j;
	if(isnan(y))
	{
		Rprintf("hull sample failed R = %lf\n", R);
	}
	return(y);
}

// computes hull normalizing constant and updates chunk cumulative
// probabilities
void initialize_hull(double *x, ARS_workspace *space, int* num_x, double *huzmax)
{
	// compute the integrals for the upper hull sections.
	int i;
	for(i = 0; i < *num_x; i++)
	{
		// first and last chambers are special cases.
		if(i == 0)
		{
			space->s[i] = (exp(space->hwv[i] + (space->z[i] - x[i])*(space->hpwv[i]) - *huzmax) -
					exp(space->hwv[0] - x[0]*(space->hpwv[0]) - *huzmax)) / (space->hpwv[i]);
			space->scum[i] = space->s[i];
		}

		else if(i == *num_x - 1)
		{
			space->s[i] = -1.0*exp(space->hwv[i] + (space->z[i-1] - x[i])*(space->hpwv[i]) - *huzmax)/(space->hpwv[i]);
			space->scum[i] = space->scum[i-1] + space->s[i];
		}

		else
		{
			space->s[i] = (exp(space->hwv[i] + (space->z[i] - x[i])*(space->hpwv[i]) - *huzmax) -
					exp(space->hwv[i] + (space->z[i-1] - x[i])*(space->hpwv[i]) - *huzmax))/(space->hpwv[i]);
			space->scum[i] = space->scum[i-1] + space->s[i];
		}
	}

	// normalize the cumulative values
	for(i = 0; i < *num_x; i++)
	{
		space->scum_norm[i] = (space->scum[i])/(space->scum[*num_x - 1]);
	}
	return;
}


// h is the log-density of our target distribution
double h(double x, double S_j, double c, double R)
{
	double out;
	out = -1.0*x*R + S_j*lgamma(x + c) - S_j*lgamma(x) + R;

	//gamma test distribution
//	  out = .001*log(x) - x/.1;
	return(out);
}

double h_prime(double x, double S_j, double c, double R)
{
	double out;
	out = -1.0*R + S_j*gsl_sf_psi(x + c) - S_j*gsl_sf_psi(x);

	//gamma test distribution, derivative.
//	out = .001/x - 10.0;
	return(out);
}






