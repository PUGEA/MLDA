/*
	learn.c
	$Id: learn.c,v 1.8 2013/01/13 14:23:27 daichi Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "learn.h"
#include "vbem.h"
#include "newton.h"
#include "likelihood.h"
#include "feature.h"
#include "dmatrix.h"
#include "fmatrix2.h"
#include "util.h"

double
lda_learn (document *data, double *alpha, double *eta, double **gammas, double **beta, double *isonum, int condition, int replicate, 
		int nclass, int nlex, int dlenmax,
		int emmax, int demmax, double epsilon)
{
	double *gamma, **q, **lambda, **lambda1, **mnt, **pmnt, **betas, *nt, *pnt, *ap,**lambdatran;
	double ppl, pppl = 0;
	double z,like=0,like1=0,datasum=0;
	document *dp;
	int i, j, t, n,ii;
	int start, elapsed;
	int cons = condition*replicate;
	/*
	 *  randomize a seed
	 *
	 */
	srand(time(NULL));

	/*
	 *	count data length
	 *
	 */
	for (dp = data, n = 0; (dp->len) != -1; dp++, n++)
		;

	/*
	 *  initialize parameters
	 *
	 */

	for (i = 0; i < nclass; i++)
	{
		alpha[i] = 50/nclass;
	}
	for (i =0; i<nlex; i++)
	{
		eta[i] = 0.1*nlex;
	}

	if ((lambda = dmatrix(nlex, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate lambda.\n");
		return;
	}
	if ((lambdatran = dmatrix(nclass, nlex)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate lambdatran.\n");
		return;
	}
	if ((lambda1 = dmatrix(nlex, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate lambda1.\n");
		return;
	}
	if ((mnt = dmatrix(nlex, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate mnt.\n");
		return;
	}
	if ((pmnt = dmatrix(nlex, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate pmnt.\n");
		return;
	}
	if ((betas = dmatrix(nlex, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate betas.\n");
		return;
	}
	/*
	 *  initialize buffers
	 *
	 */
	if ((q = dmatrix(nlex, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate q.\n");
		return;
	}
	/*if ((qs = tmatrix(condition,nlex, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate qs.\n");
		return;
	}*/
	if ((gamma = (double *)calloc(nclass, sizeof(double))) == NULL)
	{
		fprintf(stderr, "lda_learn:: cannot allocate gamma.\n");
		return;
	}
	if ((ap = (double *)calloc(nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate ap.\n");
		return;
	}
	if ((nt = (double *)calloc(nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate nt.\n");
		return;
	}
	if ((pnt = (double*)calloc(nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate pnt.\n");
		return;
	}
	for (i = 0; i<nlex; i++)
		for(j = 0; j<nclass; j++)
		{
			mnt[i][j] = 0;
			pmnt[i][j] = 0;
		}

	for (t = 0; t < emmax; t++)
	{
		for (i = 0; i < condition; i++)
			for (j = 0; j< nclass; j++)
				gammas[i][j] = alpha[j] + 1/nclass;
		for (i = 0; i<nlex; i++)
			for (j = 0; j<nclass; j++)
				lambda[i][j] = eta[i] + 1/nclass;
		//printf("iteration %2d/%3d..\t", t + 1, emmax);
		fflush(stdout);
		/*
		 *  VB-E step
		 *
		 */
		/* iterate for data */
		if (t==0)
		{
			for(i =0; i<nlex; i++)
				datasum = datasum+data->cnt[i];
			if (datasum<1)
				break;
		}
		for (ii = 0; ii<demmax; ii++)
		{
			for (dp = data, i = 0; (dp->len) != -1; dp++, i++)
			{
				vbem(dp, gamma, lambda, q, nt, pnt, ap,isonum,(const double *)alpha, (const double **)beta, dp->len, nclass, demmax);
				accum_gammas(gammas, gamma, i, nclass);
//				accum_qs(qs, q, i, nclass, nlex);
				accum_betas(lambda1, q, nclass, dp);
			}
			accum_lambda(lambda, lambda1, eta, nlex, nclass);
			for(i = 0; i<nlex; i++)
				for(j =0; j<nclass; j++)
					mnt[i][j] = lambda1[i][j];
			if((t>1) && Twoconverged(mnt,pmnt,nlex,nclass,1.0e-2))
				break;
			for (i = 0; i<nlex; i++)
				for (j =0; j<nclass; j++)
				{
					pmnt[i][j] = mnt[i][j];
					lambda1[i][j] = 0;
				}
		}

		/*
		 *  VB-M step
		 *
		 */
		/* Newton-Raphson for alpha */
		for(i =0; i<nlex; i++)
			for(j =0; j<nclass; j++)
				lambdatran[j][i] = lambda[i][j];
		newton_alpha(alpha, gammas, n, nclass, 0);
		newton_alpha(eta, lambdatran, nclass, nlex, 0);

		/*  converge?
		 *
		 */
		like = lda_lik2(data,lambdatran,gammas,n,nclass,nlex);
		elapsed = myclock() - start;
		//printf("PPL = %.04f\t", ppl); fflush(stdout);
		if ((t > 1) && (fabs(like - like1)/like) < epsilon)
		{
			break;
		}
		if(t==emmax)
			break;
		like1 = like;
		/*
		 * ETA
		 *
		 */
		/*printf("ETA:%s (%d sec/step)\r",
			rtime(elapsed * ((double) emmax / (t + 1) - 1)),
			(int)((double) elapsed / (t + 1) + 0.5));*/
	}

	free_dmatrix(betas, nlex);
	free_dmatrix(q, nlex);
	free_dmatrix(lambda, nlex);
	free_dmatrix(lambdatran, nclass);
	free_dmatrix(lambda1, nlex);
	free_dmatrix(mnt, nlex);
	free_dmatrix(pmnt, nlex);
	//free_tmatrix(qs,condition,dlenmax);
	free(gamma);
	free(ap);
	free(nt);
	free(pnt);

	return like;
}

void
accum_gammas (double **gammas, double *gamma, int n, int K)
{
	/* gammas(n,:) = gamma for Newton-Raphson of alpha */
	int k;
	for (k = 0; k < K; k++)
		gammas[n][k] = gamma[k];
	return;
}

void
accum_betas (double **betas, double **q, int K, document *dp)
{
	int i, k;
	int n = dp->len;

	for (i = 0; i < n; i++)
		for (k = 0; k < K; k++)
			betas[dp->id[i]][k] += q[i][k] * dp->cnt[i];
}

void
accum_lambda(double **lambda, double **lambda1, double *eta, int L, int K)
{
	int i,j;
	for(i= 0; i< L; i++)
		for(j =0; j<K; j++)
			lambda[i][j] = eta[i]+lambda1[i][j];
}
