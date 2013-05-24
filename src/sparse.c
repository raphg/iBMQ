/*
 * sparse.c
 *
 *  Created on: Sep 29, 2011
 *      Author: hoblitz
 */
#include<stdio.h>
#include<stdlib.h>
#include "sparse.h"
#include<gsl/gsl_matrix.h>
#include<R.h>

// initialize a memory pool for the linked list
void initializePool(int n_chunks, int n_els, ptr_memPool ptr_pool)
{
	int i;
	ptr_pool->array_head = (ptr_memChunk*) R_alloc(n_chunks, sizeof(ptr_memChunk));
	ptr_pool->n_chunks = n_chunks;

	if (ptr_pool->array_head == NULL) {
		error("Failed to allocate memory pool\n");
	}

	for(i = 0; i < n_chunks; i++) {
		ptr_pool->array_head[i] = initializeChunk(n_els);
	}
	return;
}

// returns a pointer to a dynamically allocated block of elements
ptr_memChunk initializeChunk(int n_els)
{
	// allocate a memChunk header to the pointer

	ptr_memChunk ptr_chunk;
	ptr_chunk = (ptr_memChunk) R_alloc(1, sizeof(memChunk));
	ptr_m_el el_array;

	if (ptr_chunk == NULL) {
		error("failed to allocate chunk header\n");
	}

	// allocate memory and fill in fields of the chunk.
	el_array = (ptr_m_el) R_alloc(n_els, sizeof(m_el));

	ptr_chunk->array_head = el_array;
	ptr_chunk->next_avail = el_array;
	ptr_chunk->n_elements = n_els;

	if(el_array == NULL)
	{
		error("Failed to allocate memory chunk\n");
	}

	// link together the elements of the chunk in a list
	int i;
	for(i = 0; i < (n_els - 1); i++)
	{
		el_array[i].next = &(el_array[i+1]);
		el_array[i].co = 0;
	}
	el_array[n_els - 1].co = 0;
	el_array[n_els - 1].next = NULL;

	// put address of last element in the chunk
	ptr_chunk->last_avail = &el_array[n_els - 1];
	return(ptr_chunk);
}

//
ptr_m_el checkoutElementFromChunk(ptr_memChunk ptr_chunk)
{
	ptr_m_el new_m_el, avail_m_el;

	//get address of the first available element in the chunk
	new_m_el = ptr_chunk->next_avail;

	//update the address of the next available element in the chunk
	avail_m_el = new_m_el->next;
	new_m_el->next = NULL;
	if(new_m_el->co == 1)
	{
		error("Memory pool exhausted\n");
	}
	new_m_el->co = 1;

	ptr_chunk->next_avail = avail_m_el;
	return(new_m_el);
}

void returnElementToChunk(ptr_memChunk ptr_chunk, ptr_m_el ptr_to_return)
{
	ptr_m_el ptr_last_el = ptr_chunk->last_avail;
	// put the returned element into the list among available elements
	// i.e, tack the returned element to the end
	ptr_last_el->next = ptr_to_return;
	// delete the old address that was in the returned element
	ptr_to_return->next = NULL;
	ptr_to_return->co = 0;
	//
	ptr_chunk->last_avail = ptr_to_return;
	return;
}

void deleteChunk(ptr_memChunk ptr_chunk)
{
	//SV_printlist(ptr_chunk->array_head);
	free(ptr_chunk->array_head);
	free(ptr_chunk);
	return;
}

void deletePool(ptr_memPool ptr_pool)
{
	int i, n_chunks = ptr_pool->n_chunks;
	for(i = 0; i < n_chunks; i++)
	{
		deleteChunk(ptr_pool->array_head[i]);
	}
	free(ptr_pool->array_head);
	return;
}

// retrieve vector element j from sparse vector.
double SV_get(ptr_m_el header, int j)
{
	ptr_m_el tmp_ptr;
	double x;
	tmp_ptr = header->next;
	while((tmp_ptr != NULL) && (tmp_ptr->ind < j))
	{
		tmp_ptr = tmp_ptr->next;
	}
	if((tmp_ptr == NULL) || tmp_ptr->ind > j)
	{
		Rprintf("element %d not found \n", j);
		return 0.0;
	}
	x = tmp_ptr->x;
	return x;
}

// Sparse vector/regular vector dot product: output stored in out.
// assumes that these are reasonable things to take a dot product from,
// i.e. dimensions of y and SV *header match properly.
void SV_ddot(const double *y, ptr_m_el header, double* out)
{
	*out = 0.0;
	ptr_m_el tmp_ptr;
	tmp_ptr = header->next;
	while(tmp_ptr != NULL)
	{
		*out += y[tmp_ptr->ind]*(tmp_ptr->x);
		tmp_ptr = tmp_ptr->next;
	}
	return;
}

// diagnostic tool: prints contents of a sparse vector to R console.
void SV_printlist(ptr_m_el  header)
{
	ptr_m_el  tmp_ptr;
	tmp_ptr = header->next;
	while(tmp_ptr != NULL)
	{
		Rprintf("%d %f ",tmp_ptr->ind, tmp_ptr->x);
		Rprintf("\n");
		tmp_ptr = tmp_ptr->next;
	}
	Rprintf("\n");
	return;
}

// calculates X(vec m_el) + y = y
// Matrix vector plus y:
void SV_dmvpy(double** X, ptr_m_el header, double* y, int nrow)
{
	int indx, i;
	if(header->next == NULL)
	{
		// all elements of the sparse vector are zero, change nothing and return
		return;
	}

	ptr_m_el temp_ptr;
	temp_ptr = header->next;
	while(temp_ptr != NULL)
	{
		indx = temp_ptr->ind;
		for(i = 0; i < nrow; i++)
		{
			y[i] += X[i][indx]*(temp_ptr->x);
		}
		temp_ptr = temp_ptr->next;
	}
	return;
}

// GSL matrix times sparse vector plus Y
void SV_gsl_dmvpy(const gsl_matrix* X, ptr_m_el header, double* y, int nrow)
{
	int i;
	if(header->next == NULL)
	{
		// all elements of the sparse vector are zero, change nothing and return
		return;
	}

	ptr_m_el temp_ptr;
	temp_ptr = header->next;
	while(temp_ptr != NULL)
	{
		const int indx = temp_ptr->ind;
		const double x = temp_ptr->x;
		for(i = 0; i < nrow; i++)
		{
			const double mat_el = gsl_matrix_get(X, i, indx);
			double yi = y[i];
			yi += mat_el*x;
			y[i] = yi;
		}
		temp_ptr = temp_ptr->next;
	}
	return;
}

// this function frees the matrix elements, but not the header
void SV_free(ptr_m_el header)
{
	ptr_m_el temp_ptr, temp_ptr_cur;
	temp_ptr = header->next;
	while(temp_ptr != NULL)
	{
		temp_ptr_cur = temp_ptr;
		temp_ptr = temp_ptr->next;
		free(temp_ptr_cur);
	}
}

// this function assumes element j is in the list, and removes it
// error handling is somewhat sub-par.
void SV_remove_el(ptr_m_el header, int j, ptr_memChunk ptr_chunk)
{
	ptr_m_el temp_ptr, temp_ptr_prev;
	temp_ptr = header->next;
	temp_ptr_prev = header;
	while((temp_ptr != NULL) && (temp_ptr->ind < j)) {
		temp_ptr_prev = temp_ptr;
		temp_ptr = temp_ptr->next;
	}
	// since we assume that j is in the list, when we've stopped we are
	// sure to be AT element j
	if (temp_ptr->ind != j) {
		Rprintf("failed to locate index %d in list\n", j);
		return;
	}
	// link last item to next item
	temp_ptr_prev->next = temp_ptr->next;

	// remove j'th item.
	returnElementToChunk(ptr_chunk, temp_ptr);
	return;
}

// Add an element at index j to a sparse vector object.
void SV_add_el(ptr_m_el header, int j, double val, ptr_memChunk ptr_chunk)
{
	ptr_m_el temp_ptr, temp_ptr1;
	temp_ptr = header->next;
	temp_ptr1 = header;

	// check whether list is empty, if so insert as first element
	if (temp_ptr == NULL) {
		ptr_m_el new_element;
		new_element = checkoutElementFromChunk(ptr_chunk);
		new_element->x = val;
		new_element->ind = j;
		new_element->next = NULL;
		header->next = new_element;
	}

	// list is not empty, must traverse list until we either
	// reach the end or find the current index
	else {
		while((temp_ptr != NULL) && (temp_ptr->ind < j)) {
			temp_ptr1 = temp_ptr;
			temp_ptr = temp_ptr->next;
		}

		// we are at end of list and/or we have reached the appropriate index
		if (temp_ptr == NULL) {
			// if we are at the end of the list, tmp_ptr 1 is the address
			// of the last element, and either the index is equal,
			// or it is greater.
			if(j == temp_ptr1->ind) {
				temp_ptr1->x = val;
			}
			// if index is greater then we must insert a new element
			// at the end
			else if(j > temp_ptr1->ind) {
				ptr_m_el new_element;
				new_element = checkoutElementFromChunk(ptr_chunk);
				new_element->x = val;
				new_element->ind = j;
				new_element->next = NULL;

				//link last element to new element
				temp_ptr1->next = new_element;
			}
		}
		// did not reach end of list
		else {
			// index is equal, so element is present and we need not
			// create a new one

			if(j == temp_ptr->ind) {
				temp_ptr->x = val;
			}

			// index is greater, so we're in the middle and still must create
			// a new element
			else if(j < temp_ptr->ind) {
				ptr_m_el new_element;
				new_element = checkoutElementFromChunk(ptr_chunk);;
				new_element->x = val;
				new_element->ind = j;
				new_element->next = temp_ptr;
				temp_ptr1->next = new_element;
			}
		}
	}
	return;
}

