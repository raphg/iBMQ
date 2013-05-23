/*
 * norm_gamma_generation.h
 *
 *  Created on: Oct 6, 2011
 *      Author: hoblitz
 */

#ifndef NORM_GAMMA_GENERATION_H_
#define NORM_GAMMA_GENERATION_H_

#endif /* NORM_GAMMA_GENERATION_H_ */

#include "RngStream.h"
#include <float.h>
#include <R.h>
#include <math.h>

double RngStream_N01 (const RngStream r);
double RngStream_GA1 (const double a, RngStream r);
double RngStream_Beta (const double a, const double b, RngStream r);
double RngStream_LogitBeta(const double a, const double b, RngStream rng);
