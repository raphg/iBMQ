/*
 * sparse.h
 *
 *  Created on: Sep 29, 2011
 *      Author: hoblitz
 */

#ifndef SPARSE_H_
#define SPARSE_H_
#include<gsl/gsl_matrix.h>

typedef struct MemoryPool *ptr_memPool;
typedef struct MemoryChunk *ptr_memChunk;
typedef struct M_EL *ptr_m_el;

typedef struct MemoryPool {
	ptr_memChunk* array_head;
	int n_chunks;
} memPool;

typedef struct MemoryChunk {
	ptr_m_el next_avail;
	ptr_m_el last_avail;
	ptr_m_el array_head;
	int n_elements;
} memChunk;

typedef struct M_EL {
	int ind;
	int co;
	double x;
	ptr_m_el next;
} m_el;

#define N_EL_MIN 3

void deletePool(ptr_memPool ptr_pool);
void deleteChunk(ptr_memChunk ptr_chunk);
ptr_memChunk initializeChunk(int n_els);
void initializePool(int n_chunks, int n_els, ptr_memPool ptr_pool);
ptr_m_el checkoutElementFromChunk(ptr_memChunk ptr_chunk);
void returnElementToChunk(ptr_memChunk ptr_chunk, ptr_m_el ptr_to_return);

void SV_gsl_dmvpy(const gsl_matrix* X, ptr_m_el header, double* y, int nrow);
void SV_printlist(ptr_m_el header);
void SV_remove_el(ptr_m_el header, int j, ptr_memChunk ptr_chunk);
void SV_add_el(ptr_m_el header, int j, double val, ptr_memChunk ptr_chunk);
void SV_dmvpy(double** X, ptr_m_el header, double* y, int nrow);
void SV_free(ptr_m_el header);
void SV_ddot(const double *x, ptr_m_el header, double* out);
double SV_get(ptr_m_el header, int j);

#endif /* SPARSE_H_ */
