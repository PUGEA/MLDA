/*
    dmatrix.c
    an implementation of double matrix.
    $Id: dmatrix.c,v 1.3 2004/11/01 01:01:06 dmochiha Exp $

*/
#include <stdlib.h>
#include "dmatrix.h"

double **
dmatrix (int rows, int cols)
{
	double **matrix;
	int i;

	matrix = (double **)calloc(rows, sizeof(double *));
	if (matrix == NULL)
		return NULL;
	for (i = 0; i < rows; i++) {
		matrix[i] = (double *)calloc(cols, sizeof(double));
		if (matrix[i] == NULL)
			return NULL;
	}
	
	return matrix;
}

void
free_matrix (double **matrix, int rows)
{
	int i;
	for (i = 0; i < rows; i++)
		free(matrix[i]);
	free(matrix);
}

double ***
tmatrix (int tri, int rows, int cols)
{
	double ***matrix;
	int i,j;
	matrix = (double ***)calloc(tri, sizeof(double **));
	if (matrix == NULL)
		return NULL;
	for (i = 0; i < tri; i++) {
		matrix[i] = (double **)calloc(rows, sizeof(double *));
		for (j =0 ; j<rows; j++)
		{
			matrix[i][j] = (double *)calloc(cols, sizeof(double));
			if (matrix[i][j] == NULL)
				return NULL;
		}
	}
	return matrix;
}

void
free_tmatrix (double ***matrix,int tri, int rows)
{
	int i,j;
	for (i = 0; i < tri; i++)
		for (j = 0; j < rows; j++)
			free(matrix[i][j]);
	free(matrix);
}


