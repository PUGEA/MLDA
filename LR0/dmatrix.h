/*
    dmatrix.h
    a header file of double matrix.
    $Id: dmatrix.h,v 1.2 2004/10/28 06:52:22 dmochiha Exp $
    
*/
#ifndef __DMATRIX_H__
#define __DMATRIX_H__

extern double **dmatrix(int rows, int cols);
extern void free_dmatrix(double **matrix, int rows);
extern double ***tmatrix(int tri, int rows, int cols);
extern void free_tmatrix(double ***matrix, int tri, int rows);
#endif
