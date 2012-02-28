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

// retrieve vector element j from sparse vector.
double SV_get(m_el *header, int j)
{
	m_el *tmp_ptr;
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
void SV_ddot(const double *y, m_el *header, double* out)
{
	*out = 0.0;
	m_el *tmp_ptr;
	tmp_ptr = header->next;
	while(tmp_ptr != NULL)
	{
		*out += y[tmp_ptr->ind]*(tmp_ptr->x);
		tmp_ptr = tmp_ptr->next;
	}
	return;
}

// diagnostic tool: prints contents of a sparse vector to R console.
void SV_printlist(m_el *header)
{
	m_el *tmp_ptr;
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
void SV_dmvpy(double** X, m_el *header, double* y, int nrow)
{
	int indx, i;
	if(header->next == NULL)
	{
		// all elements of the sparse vector are zero, change nothing and return
		return;
	}

	m_el *temp_ptr;
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
void SV_gsl_dmvpy(gsl_matrix* X, m_el *header, double* y, int nrow)
{
	int indx, i;
	if(header->next == NULL)
	{
		// all elements of the sparse vector are zero, change nothing and return
		return;
	}

	m_el *temp_ptr;
	temp_ptr = header->next;
	while(temp_ptr != NULL)
	{
		indx = temp_ptr->ind;
		for(i = 0; i < nrow; i++)
		{
			y[i] += gsl_matrix_get(X, i, indx)*(temp_ptr->x);
		}
		temp_ptr = temp_ptr->next;
	}
	return;
}

// this function frees the matrix elements, but not the header
void SV_free(m_el *header)
{
	m_el *temp_ptr, *temp_ptr_cur;
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
void SV_remove_el(m_el *header, int j)
{
	m_el *temp_ptr, *temp_ptr_prev;
	temp_ptr = header->next;
	temp_ptr_prev = header;
	while((temp_ptr != NULL) && (temp_ptr->ind < j))
	{
		temp_ptr_prev = temp_ptr;
		temp_ptr = temp_ptr->next;
	}
	// since we assume that j is in the list, when we've stopped we are
	// sure to be AT element j
	if(temp_ptr->ind != j)
	{
		printf("failed to locate index %d in list\n", j);
		return;
	}
	// link last item to next item
	temp_ptr_prev->next = temp_ptr->next;

	// remove j'th item.
	free(temp_ptr);
	return;
}

// Add an element at index j to a sparse vector object.
void SV_add_el(m_el *header, int j, double val)
{
	int temp_ind;
	m_el *temp_ptr, *temp_ptr1;
	temp_ptr = header->next;
	temp_ptr1 = header;

	// check whether list is empty, if so insert as first element
	if(temp_ptr == NULL)
	{
		m_el *new_element;
		new_element = malloc(sizeof(m_el));
		if(new_element == NULL)
		{
			Rprintf("allocation failed \n");
			return;
		}
		new_element->x = val;
		new_element->ind = j;
		new_element->next = NULL;
		header->next = new_element;
	}

	// list is not empty, must traverse list until we either
	// reach the end or find the current index
	else
	{
		temp_ind = temp_ptr->ind;
		while((temp_ptr != NULL) && (temp_ptr->ind < j))
		{
			temp_ptr1 = temp_ptr;
			temp_ptr = temp_ptr->next;
		}


		// we are at end of list and/or we have reached the appropriate index
		if(temp_ptr == NULL)
		{
			// if we are at the end of the list, tmp_ptr 1 is the address of the last
			// element, and
			// either the index is equal,
			// or it is greater.
			if(j == temp_ptr1->ind)
			{
				temp_ptr1->x = val;
			}
			// if index is greater then we must insert a new element at the end
			else if(j > temp_ptr1->ind)
			{
				m_el *new_element;
				new_element = malloc(sizeof(m_el));
				if(new_element == NULL)
				{
					Rprintf("allocation failed \n");
					return;
				}
				new_element->x = val;
				new_element->ind = j;
				new_element->next = NULL;
				//link last element to new element
				temp_ptr1->next = new_element;
			}
		}
		// did not reach end of list
		else
		{
			// index is equal, so element is present and we need not
			// create a new one

			if(j == temp_ptr->ind)
			{
				temp_ptr->x = val;
			}

			// index is greater, so we're in the middle and still must create
			// a new element
			else if(j < temp_ptr->ind)
			{
				m_el *new_element;
				new_element = malloc(sizeof(m_el));
				if(new_element == NULL)
				{
					Rprintf("allocation failed \n");
					return;
				}
				new_element->x = val;
				new_element->ind = j;
				new_element->next = temp_ptr;
				temp_ptr1->next = new_element;
			}
		}
	}
	return;
}

