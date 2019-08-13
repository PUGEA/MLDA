#include <stdio.h>
#include <stdlib.h>
#include "mnormalize.h"

int mnormalize(double** data ,int nrows,int ncols,int arrow)
{
	int i,j;
	if (arrow==1)
	{
		double datasum[ncols];
		for (i=0; i<ncols; i++)
			datasum[i]=0;
		for (i = 0; i<nrows; i++)
			for(j = 0; j<ncols; j++)
				datasum[j] = datasum[j] + data[i][j];
		for (i = 0; i<nrows; i++)
			for(j = 0; j<ncols; j++)
			{
				if (data[i][j]!=0)
					data[i][j]= data[i][j]/datasum[j];
				else
					data[i][j] = 0;
			}
	}
	if (arrow==2)
	{
		double datasum[nrows];
		for(i =0; i<nrows; i++)
			datasum[i] = 0;
		for (i = 0; i<nrows; i++)
			for(j = 0; j<ncols; j++)
				datasum[i] = datasum[i] + data[i][j];
		for (i = 0; i<nrows; i++)
			for(j = 0; j<ncols; j++)
			{
				if (data[i][j]!=0)
					data[i][j] = data[i][j]/datasum[i];
				else
					data[i][j] = 0;
			}
	}
	return (1);
}

int clac(double** data ,int nrows,int ncols,int arrow)
{
	int i,j;
	if (arrow==1)
	{
		double datasum[ncols];
		for (i=0; i<ncols; i++)
			datasum[i]=0;
		for (i = 0; i<nrows; i++)
			for(j = 0; j<ncols; j++)
				datasum[j] = datasum[j] + data[i][j];
		for (j = 0; j<ncols; j++)
			if (datasum[j]==0)
				for(i =0; i<nrows; i++)
					data[i][j] = 1;
	}
	if (arrow==2)
	{
		double datasum[nrows];
		for(i =0; i<nrows; i++)
			datasum[i] = 0;
		for (i = 0; i<nrows; i++)
			for(j = 0; j<ncols; j++)
				datasum[i] = datasum[i] + data[i][j];
		for (i = 0; i<nrows; i++)
			if (datasum[i]==0)
				for(j = 0; j<ncols; j++)
					data[i][j] = 1;
	}
	return (1);
}