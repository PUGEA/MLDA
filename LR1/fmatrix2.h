/*
    fmatrix2.h
    a header file of double matrix.
    $Id: fmatrix2.h,v 1.2 2019/7/16 06:52:22 dmochiha Exp $
    
*/
#ifndef __FMATRIX2_H__
#define __FMATRIX2_H__

extern double **fmatrix2(int rows, int cols);
extern void free_fmatrix2(double **matrix, int rows);

#endif