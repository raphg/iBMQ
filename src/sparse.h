/*
 * sparse.h
 *
 *  Created on: Sep 29, 2011
 *      Author: hoblitz
 */

#ifndef SPARSE_H_
#define SPARSE_H_
#include<gsl/gsl_matrix.h>

#endif /* SPARSE_H_ */

typedef struct M_EL {
	int ind;
	double x;
	struct M_EL *next;
} m_el;

void SV_gsl_dmvpy(gsl_matrix* X, m_el *header, double* y, int nrow);
void SV_printlist(m_el *header);
void SV_remove_el(m_el *header, int j);
void SV_add_el(m_el *header, int j, double val);
void SV_dmvpy(double** X, m_el *header, double* y, int nrow);
void SV_free(m_el *header);
void SV_ddot(const double *x, m_el *header, double* out);
double SV_get(m_el *header, int j);

