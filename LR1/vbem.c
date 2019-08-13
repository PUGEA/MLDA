/*
	vbem.c
	Latent Dirichlet Allocation, VB-EM for a document.
	$Id: vbem.c,v 1.6 2004/11/06 08:26:19 dmochiha Exp $

*/
#include <math.h>
#include "lda.h"
#include "gamma.h"
#include "feature.h"
#include "util.h"
#include "mnormalize.h"
#include "dmatrix.h"

void
vbem (const document *d, double *gamma, double **lambda, double **q,
	double *nt, double *pnt, double *ap,double *isonum,
	const double *alpha, const double **beta,
	int L, int K, int emmax)
{
	int i,j, k, l;
	double z;
	double isotemp[K],lambdas[K];
	double tempq[L];
	for(i =0; i<L;i++)
		tempq[i] = 0;
	double **b;
	if ((b = dmatrix(L, K)) == NULL) {
		fprintf(stderr, "lda:: cannot allocate b.\n");
		exit(1);
	}
	for (k = 0; k < K; k++)
	{
		nt[k] = (double) L / K;
		pnt[k] = (double) L/ K;
		lambdas[k] = 0;
		isotemp[k] = 1/isonum[k];
	}
	for (j = 0 ;j<L; j++)
		for (k =0; k<K; k++)
			lambdas[k] = lambdas[k] +lambda[j][k];
	for (j = 0 ;j<L; j++)
		for (k =0; k<K; k++)
			q[j][k] = 1/K;
	for (j = 0; j < emmax; j++)
	{
		/* vb-estep */
		for (k = 0; k < K; k++)
		{
			ap[k] = psi(gamma[k]);
			ap[k] = exp(ap[k]);
			//ap[k] = exp(psi(gamma[k]));
		}
		for (l = 0; l < L; l++)
			for (k = 0; k < K; k++)
				b[l][k] = exp((psi(lambda[l][k])-psi(lambdas[k])))*beta[l][k];
		for (l = 0; l < L; l++)
			for (k = 0; k < K; k++)
			{
				q[l][k] = b[l][k]*ap[k]*isotemp[k];
				tempq[l] = tempq[l] +q[l][k];
			}
		for (l = 0; l < L; l++)
			if (tempq[l]==0)
				for (k = 0; k < K; k++)
					q[l][k] = 1;
		mnormalize(q, L, K, 2);
		/* vb-mstep */
		for (k = 0; k < K; k++) {
			z = 0;
			for (l = 0; l < L; l++)
				z += q[l][k] * d->cnt[l];
			nt[k] = z;
		}
		/* converge? */
		if ((j > 0) && converged(nt, pnt, K, 1.0e-2))
			break;
		for (k = 0; k < K; k++)
			pnt[k] = nt[k];
	}
	for (k = 0; k < K; k++)
		gamma[k] = alpha[k] + nt[k];
	free_dmatrix(b, L);
	return;
}

