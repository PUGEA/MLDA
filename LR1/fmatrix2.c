/*
    fmatrix2.c
    an implementation of double matrix.
    $Id: fmatrix2.c,v 1.3 2019/7/16 06:52:22 dmochiha Exp $

*/
#include <stdlib.h>
#include "fmatrix2.h"

double **
fmatrix2 (int rows, int cols)
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
free_dmatrix (double **matrix, int rows)
{
	int i;
	for (i = 0; i < rows; i++)
		free(matrix[i]);
	free(matrix);
}

