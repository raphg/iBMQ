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
 * x: starting lattice of points
 * num_x: number of starting points
 * nmax: maximum number of points in lattice
 *
 * ARS_workspace is a struct with elements
 * hwv: working vector for the h(x) values
 * hpwv: working vector for h'(x) values
 * scum: working vector for cumulative integral values of the hull section
 * scum_norm: working vector cumulative density of hull sections
 * z: intersection of hull tangents
 * s: working vector for calculating scum
 *
 * h: function taking a double and a vector of arguments, returns log density
 * h_prime: function taking double and vector of arguments, returns derivative of log density
 * rng: random number generator
 * eps: double machine precision value
 */


// returns -1.0 if an error is detected.
double sample_conditional(double* restrict x,
		int* restrict num_x,
		int nmax,
		double* restrict argvec,
		ARS_workspace *ws,
		RngStream rng,
		double eps, double (*h)(const double, const double *),
		double (*h_prime)(const double , const double *))
{
	int i, u_section, l_section;
	int accept = 0;
	int attempts = 0;
	int update = 0;

	double x_samp, W, l, u, hnew, V;
	double huzmax, test1, test2, tmp0, tmp1;

	double* restrict hwv = ws->hwv;
	double* restrict hpwv = ws->hpwv;
	double* restrict z = ws->z;

	// intialize the lattice
	for(i = 0; i < *num_x; i++)
	{
		const double xval = x[i];
		const double nhwv = h(xval, argvec);
		const double nhpwv = h_prime(xval, argvec);

		hwv[i] = nhwv;
		hpwv[i] = nhpwv;
	}

	// need last hpwv  < 0, move x_last farther out until this is so. This only works with
	// strict concavity!

	while(hpwv[*num_x-1] >= 0)
	{
		const double xnew = x[*num_x - 1] + 2.0;
		x[*num_x - 1] = xnew;

		const double nhpwv = h_prime(x[*num_x - 1], argvec);
		hpwv[*num_x - 1] = nhpwv;
	}
	hwv[*num_x - 1] = h(x[*num_x - 1], argvec);

	//now calculate intersection of lines, z
	// i-th element corresponds to intersection to right of x[i]
	for(i = 0; i < (*num_x - 1); i++)
	{
		if(i == 0)
		{
			// initialize the search largest y value in outer hull
			huzmax = hwv[0] - x[0]*hpwv[0];
		}

		const double x_i = x[i];
		const double x_ip1 = x[i+1];
		const double hwv_i = hwv[i];
		const double hwv_ip1 = hwv[i+1];
		const double hpwv_i = hpwv[i];
		const double hpwv_ip1 = hpwv[i+1];
		double znew;

		// tangents are nearly parallel, use midpoint
		// approximation to avoid instability

		if((hpwv_i - hpwv_ip1) <= eps)
		{

			znew = (x_i + x_ip1)/2.0;
			test1 = huzmax;
			test2 = hwv_i + (znew - x_i) * (hpwv_i);
			huzmax = (test1 > test2) ? test1 : test2;
		}
		else
		{
			znew = (hwv_ip1 - hwv_i - x_ip1*hpwv_ip1 +
					x_i*hpwv_i)/
					(hpwv_i - hpwv_ip1);
			test1 = huzmax;
			test2 = hwv_i + (znew - x_i) * hpwv_i;
			huzmax = (test1 > test2) ? test1 : test2;
		}
		z[i] = znew;
	}

	//huzmax holds the maximum value of the hull on the log scale

	// this initializes the hull, i.e. creates the cumulative distribution function of the
	// outer hull.
	initialize_hull(x, ws, *num_x, huzmax);

	while(attempts < 100000)
	{
		V = RngStream_RandU01(rng);

		x_samp = sample_hull(x, ws, num_x, &u_section, V, huzmax);
		attempts++;

		// check for sampling trouble
		check_sample(x_samp, x, ws, num_x);

		W = RngStream_RandU01(rng);

		// u_section is the j corresponding to x[j]
		u = hwv[u_section] + (x_samp - x[u_section])*(hpwv[u_section]);

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
			l =	((x[j] - x_samp)*(hwv[j-1]) + (x_samp - x[j-1])*(hwv[j]))/
					(x[j] - x[j-1]);

			//squeeze test.
			if(log(W) <= (l - u))
			{
				// acquire 15th and 85th percentile.
				tmp0 = sample_hull(x, ws, num_x, &u_section, .15, huzmax);
				tmp1 = sample_hull(x, ws, num_x, &u_section, .85, huzmax);

				x[0]=tmp0;
				x[1]=tmp1;

				check_sample(x[0], x, ws, num_x);
				check_sample(x[1], x, ws, num_x);

				return(x_samp);
			}
			// failed squeeze test, need to evaluate h at x_samp
			// and compare directly to outer hull.
			else
			{
				hnew = h(x_samp, argvec);
				if(log(W) <= (hnew - u))
				{
					// sample accepted!
					//acquire 15th and 85th percentiles
					tmp0 = sample_hull(x, ws, num_x, &u_section, .15, huzmax);
					tmp1 = sample_hull(x, ws, num_x, &u_section, .85, huzmax);

					x[0]=tmp0;
					x[1]=tmp1;

					check_sample(x[0], x, ws, num_x);
					check_sample(x[1], x, ws, num_x);

					return(x_samp);
				}
				else
				{
					// failed direct comparison, we add the new sampled value to the hull
					update = update_hull(x, ws, argvec, num_x, nmax,
							x_samp, hnew, l_section, &huzmax, h, h_prime);

					// if update == 1 we added a new point to the lattice
					if(update == 1)
					{
						initialize_hull(x, ws, *num_x, huzmax);
					}
				}
			}
		}

		else
		{
			// we're on either the left or right edge or the
			// hull when sampling here
			hnew = h(x_samp, argvec);
			if(log(W) <= (hnew - u))
			{
				// sample accepted
				// gather 15th and 85th percentiles for next sample
				tmp0 = sample_hull(x, ws, num_x, &u_section, .15, huzmax);
				tmp1 = sample_hull(x, ws, num_x, &u_section, .85, huzmax);

				x[0]=tmp0;
				x[1]=tmp1;

				check_sample(x[0], x, ws, num_x);
				check_sample(x[1], x, ws, num_x);

				return(x_samp);
			}
			else
			{
				update = update_hull(x, ws, argvec, num_x, nmax,
						x_samp, hnew, l_section, &huzmax, h, h_prime);

				// if update == 1 we added a new point to the lattice
				// and we must re-normalize the outer hull
				if(update == 1)
				{
					initialize_hull(x, ws, *num_x, huzmax);
				}
			}
		}
	}

	print_hull(x, ws, num_x);
	error("Rejection Sampler failed \n");
}

// insert a new abcissae in the hull, update the Z values, H, Hprime, max hull values
int update_hull(double* restrict x,
		ARS_workspace *ws,
		double* restrict argvec,
		int* restrict num_x,
		int nmax,
		double xnew,
		double hnew,
		int l_section,
		double* restrict huzmax,
		double (*h)(const double, const double *),
		double (*h_prime)(const double , const double *))
{
	if(nmax == *num_x)
	{
		// not allowed to add more breakpoints;
		// do not update hull and return to sampler
		return(0);
	}
	double hpnew, test2, test1;
	double* restrict hwv = ws->hwv;
	double* restrict hpwv = ws->hpwv;
	double* restrict z = ws->z;
	int i;

	hpnew = h_prime(xnew, argvec);
	// we are past l_section, so x_samp goes between
	// x[l_section - 1] and x[l_section]
	if(l_section == 0)
	{
		//new sample fell in leftmost chamber
		*num_x = *num_x + 1;

		for(i = *num_x - 1; i > 0; i--)
		{
			x[i] = x[i-1];
			hwv[i] = hwv[i-1];
			hpwv[i] = hpwv[i-1];
		}

		x[0] = xnew;
		hwv[0] = hnew;
		hpwv[0] = hpnew;

		for(i = *num_x - 2; i > 0; i--)
		{
			z[i] = z[i-1];
		}
		z[0] = (hwv[1] - hwv[0] - x[1]*(hpwv[1]) + x[0]*(hpwv[0]))/
				(hpwv[0] - hpwv[1]);
	}
	else if(l_section == *num_x)
	{
		// new sample goes on the "right-hand" end of the vector.
		x[*num_x] = xnew;
		hwv[*num_x] = hnew;
		hpwv[*num_x] = hpnew;

		//compute midpoint
		z[*num_x - 1] = (hwv[*num_x] - hwv[*num_x-1] - x[*num_x]*(hpwv[*num_x]) +
				x[*num_x-1]*(hpwv[*num_x-1]))/
				(hpwv[*num_x-1] - hpwv[*num_x]);
		*num_x = *num_x + 1;
	}
	else
	{
		// new sample goes between x[l_section - 1] and x[l_section]
		// x[l_section - 1] and below remain untouched

		for(i = *num_x - 1; i > l_section - 1; i--)
		{
			x[i+1] = x[i];
			hwv[i+1] = hwv[i];
			hpwv[i+1] = hpwv[i];
		}

		x[l_section] = xnew;
		hwv[l_section] = hnew;
		hpwv[l_section] = hpnew;

		// shift old midpoints to the right
		for(i = *num_x - 2; i > l_section - 1; i--)
		{
			z[i+1] = z[i];
			//hzwv[i+1] = hzwv[i];
		}
		// compute new midpoints
		for(i = l_section - 1; i < l_section + 1; i++)
		{
			const double x_i = x[i];
			const double x_ip1 = x[i+1];
			const double hwv_i = hwv[i];
			const double hwv_ip1 = hwv[i+1];
			const double hpwv_i = hpwv[i];
			const double hpwv_ip1 = hpwv[i+1];
			double znew;

			znew = (hwv_ip1 - hwv_i - x_ip1*hpwv_ip1 +
					x_i*hpwv_i)/
					(hpwv_i - hpwv_ip1);
			z[i] = znew;
		}
		*num_x = *num_x + 1;
	}

	// need to update the hull-max-value
	for(i = 0; i < *num_x - 1; i++)
	{
		const double x_i = x[i];
		const double hwv_i = hwv[i];
		const double hpwv_i = hpwv[i];
		const double z_i = z[i];

		if(i == 0)
		{
			// initialize the search largest y value in outer hull
			*huzmax = hwv_i - x_i * hpwv_i;
		}
		test1 = *huzmax;
		test2 = hwv_i + (z_i - x_i) * hpwv_i;
		*huzmax = (test1 > test2) ? test1 : test2;
	}

	return(1);
}

// this function doubles as the quantile function of the outer hull
double sample_hull(double* restrict x,
		ARS_workspace *ws,
		int* restrict num_x,
		int* restrict section,
		double p,
		double huzmax)
{
	int j = 0;
	double R, y;
	double a, b, pstar;
	double logp = log(p);

	double* restrict scum_norm = ws->scum_norm;
	double* restrict scum = ws->scum;
	double* restrict hwv = ws->hwv;
	double* restrict hpwv = ws->hpwv;
	double* restrict z = ws->z;

	// find which "chunk" we're in. returns the index of the j such that
	// our sampled x is between z[j-1] and z[j] (z[-1] = 0, z[num_x] = infinity)
	while(logp > scum_norm[j])
	{
		j++;
	}

	const double hwv_j = hwv[j];
	const double hpwv_j	= hpwv[j];
	const double x_j = x[j];
	const double z_j = z[j];
	const double z_jm1 = z[j - 1];
	const double scum_max = scum[*num_x - 1];
	const double scum_jm1 = scum[j-1];

	b = hwv_j - huzmax - hpwv_j*x_j;

	if(j == 0)
	{
		pstar = p*exp(scum_max);
		y = (log(pstar*hpwv_j + exp(b)) - b)/hpwv_j;
	}
	else
	{
		pstar = p*exp(scum_max) - exp(scum_jm1);
		y = (log(pstar * hpwv_j + exp(hpwv_j * z_jm1 + b)) - b)/hpwv_j;
	}

	*section = j;
	if(isnan(y) || isinf(y) || y <= 0)
	{
		Rprintf("hull sample failed j = %d, p = %.5lf\n", j, p);
		Rprintf("a = %.3lf, b = %.3lf, pstar = %.3lf\n", hpwv_j, b, pstar);
	}
	return(y);
}

// computes hull normalizing constant and updates cumulative
// probabilities by hull section
void initialize_hull(double* restrict x,
		ARS_workspace *ws,
		int num_x,
		double huzmax)
{
	// compute the integrals for the upper hull sections.
	int i;
	double tmp;
	double* restrict hpwv = ws->hpwv;
	double* restrict hwv = ws->hwv;
	double* restrict z = ws->z;
	double* restrict s = ws->s;
	double* restrict scum = ws->scum;
	double* restrict scum_norm = ws->scum_norm;

	for(i = 0; i < num_x; i++)
	{
		const double hpwv_i = hpwv[i];
		const double hwv_i = hwv[i];
		const double z_i = z[i];
		const double z_im1 = z[i - 1];
		const double x_i = x[i];
		double snew;

		// first and last chambers are special cases.
		// get the log of the integral from z[i-1] to z[i] under the hull
		if(i == 0)
		{
			if(hpwv_i > 0)
			{
				snew = (hwv_i + (z_i- x_i)*hpwv_i - huzmax) +
						log1p(-1.0*exp(-1.0*z_i*hpwv_i)) - log(hpwv_i);
			}
			else
			{
				snew = (hwv_i - x_i*hpwv_i - huzmax) +
						log1p(-1.0*exp(z_i*hpwv_i)) - log(fabs(hpwv_i));
			}
			s[i] = snew;
			scum[i] = snew;
		}
		// last chamber
		else if(i == num_x - 1)
		{
			snew = (hwv_i + (z_im1 - x_i) * hpwv_i - huzmax) -
					log(fabs(hpwv_i));
			scum[i] = log_apb(scum[i-1], snew);
			s[i] = snew;
		}
		// middle chambers
		else
		{
			if(hpwv_i > 0)
			{
				snew = (hwv_i + (z_i - x_i)*(hpwv_i) - huzmax) +
						log1p(-1.0*exp((z_im1 - z_i)*hpwv_i)) -
						 log(hpwv_i);
				s[i] = snew;
			}
			else if(hpwv_i < 0)
			{
				snew = (hwv_i + (z_im1 - x_i)*(hpwv_i) - huzmax) +
						log1p(-1.0*exp((z_i - z_im1)*hpwv_i)) -
						log(fabs(hpwv_i));
				s[i] = snew;
			}
			// tangent is zero
			else
			{
				snew = log(z_i - z_im1) + hwv_i - huzmax;
				s[i] = snew;
			}
			scum[i] = log_apb(scum[i-1], snew);
		}
	}
	// normalize the cumulative values
	for(i = 0; i < num_x; i++)
	{
		scum_norm[i] = (scum[i]) - (scum[num_x - 1]);
	}
	return;
}

void print_hull(double *x, ARS_workspace *ws, int *num_x)
{
	int i;
	for(i = 0; i < *num_x; i++)
	{
		Rprintf("i = %d, x = %.3lf, scum = %.3lf, scum_norm = %.3lf, hpwv = %.3lf, hwv = %.3lf,"
				"s = %.3lf \n",
				i, x[i], ws->scum[i], ws->scum_norm[i], ws->hpwv[i], ws->hwv[i], ws->s[i]);
	}
	for(i = 0; i < *num_x - 1; i++)
	{
		Rprintf("z_%d = %lf\n", i, ws->z[i]);
	}
	return;
}

// check that the sampler gave a valid sample x_samp
void check_sample(double x_samp, double *x, ARS_workspace *ws, int *num_x)
{
	int i;

	if(isnan(x_samp) || isinf(x_samp) || (x_samp <= 0.0))
	{
		for(i = 0; i < *num_x; i++)
		{
			Rprintf("i = %d, x = %.3lf, scum = %.3lf, scum_norm = %.3lf, hpwv = %.3lf, hwv = %.3lf,"
					"s = %.3lf \n",
					i, x[i], ws->scum[i], ws->scum_norm[i], ws->hpwv[i], ws->hwv[i], ws->s[i]);
		}
		for(i = 0; i < *num_x - 1; i++)
		{
			Rprintf("z_%d = %lf\n", i, ws->z[i]);
		}
		Rprintf("x_samp = %lf \n", x_samp);
		error("invalid x sample in function ARS\n");
	}
	return;
}

//Return log(a+b) from log(a) and log(b)
double log_apb(double loga, double logb)
{
	double tmp_diff = loga - logb;

	if(tmp_diff > 0)
	{
		return(loga + log1p(exp(-tmp_diff)));
	}
	else if(tmp_diff <= 0)
	{
		return(logb + log1p(exp(tmp_diff)));
	}
	else
	{
		return(loga + logb);
	}

}
