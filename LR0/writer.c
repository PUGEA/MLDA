/*
	writer.c
	an implementation of vector/matrix writer.
	$Id: writer.c,v 1.3 2004/11/06 09:44:42 dmochiha Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "writer.h"

void
lda_write (FILE *ap, FILE *bp, double *alpha, double **beta,
		int nclass, int nlex)
{
	//printf("writing model..\n"); fflush(stdout);
	write_vector(ap, alpha, nclass+1);
	write_matrix(bp, beta, nlex, nclass+1);
	//printf("done.\n"); fflush(stdout);
}

void
write_vector (FILE *fp, double *vector, int n)
{
	int i;
	for (i = 0; i < n; i++)
		fprintf(fp, "%.7e%s", vector[i], (i == n - 1) ? "\n" : "   ");
}

void
write_matrix (FILE *fp, double **matrix, int rows, int cols)
{
	int i, j;
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
			fprintf(fp, "%.7e%s", matrix[i][j],
				(j == cols - 1) ? "\n" : "   ");
}
void
write_log(char *exp_name,double*isoexpre,int nclass,char* gene_name)
{
	int cnt=0;
	FILE* input =fopen(exp_name,"a");
	fprintf(input,"%s\t",gene_name);
	for(cnt=0;cnt<nclass;cnt++)
	{
		fprintf(input,"%.7e\t",isoexpre[cnt]);
	}
	fprintf(input,"%s","\n");
	fclose(input);
}
void
write_OUTCOME(const char *exp_name,double *alpha, double *eta, double **gammas,double **iso_expre,double *gene_expre,int condition,int replicate,int nclass,int nlex,char* gene_name, double like)
{
	int cnt=0,i,j;

	char *alpha_path = (char*)calloc(256,sizeof(char));
	strcpy(alpha_path,exp_name);
	strcat(alpha_path,"/alpha");

	char *gammas_path = (char*)calloc(256,sizeof(char));
	strcpy(gammas_path,exp_name);
	strcat(gammas_path,"/gammas");

	char *eta_path = (char*)calloc(256,sizeof(char));
	strcpy(eta_path,exp_name);
	strcat(eta_path,"/eta");

	char *isoexp_path = (char*)calloc(256,sizeof(char));
	strcpy(isoexp_path,exp_name);
	strcat(isoexp_path,"/isoExpre");

	char *geneexp_path = (char*)calloc(256,sizeof(char));
	strcpy(geneexp_path,exp_name);
	strcat(geneexp_path,"/geneExpre");

	char *like_path = (char*)calloc(256,sizeof(char));
	strcpy(like_path,exp_name);
	strcat(like_path,"/likelihood");


	FILE* falpha = fopen(alpha_path,"a");
	FILE* feta = fopen(eta_path,"a");
	FILE* fgammas = fopen(gammas_path,"a");
	FILE* fiso = fopen(isoexp_path,"a");
	FILE* fgene = fopen(geneexp_path,"a");
	FILE* flike = fopen(like_path,"a");

	fprintf(falpha, "%s\n", gene_name);
	fprintf(feta, "%s\n", gene_name);
	fprintf(fgammas, "%s\n", gene_name);
	fprintf(fiso, "%s\n", gene_name);
	fprintf(fgene, "%s\n", gene_name);
	fprintf(flike, "%s\n", gene_name);

	for(cnt= 0; cnt<nclass; cnt++)
		fprintf(falpha,"%.7e\t",alpha[cnt]);
	fprintf(falpha, "\n");

	for(cnt=0; cnt<nlex; cnt++)
		fprintf(feta, "%.7e\t", eta[cnt]);
	fprintf(feta, "\n");

	for(i = 0; i<condition; i++)
	{
		for(j = 0; j<nclass; j++)
			fprintf(fgammas, "%.7e\t", gammas[i][j]);
		fprintf(fgammas, "\n");
	}

	fprintf(flike, "%.7e\n", like);

	for(i =0 ; i<condition*replicate; i++)
	{
		for(j =0; j<nclass; j++)
			fprintf(fiso, "%.7e\t", iso_expre[i][j]);
		fprintf(fiso, "\n");
	}

	for(i =0; i<condition*replicate; i++)
		fprintf(fgene, "%.7e\n", gene_expre[i]);

	fclose(falpha);
	fclose(feta);
	fclose(fgammas);
	fclose(fiso);
	fclose(fgene);
	fclose(flike);

	free(alpha_path);
	free(eta_path);
	free(gammas_path);
	free(isoexp_path);
	free(geneexp_path);
	free(like_path);
}
