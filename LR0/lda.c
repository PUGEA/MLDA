/*
    lda.c
    Latent Dirichlet Allocation, main driver.
    (c) 2004 Daichi Mochihashi, All Rights Reserved.
    $Id: lda.c,v 1.11 2017/7/14 04:23:06 dmochiha Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "lda.h"
#include "learn.h"
#include "writer.h"
#include "feature.h"
#include "dmatrix.h"
#include "fmatrix2.h"
#include "likelihood.h"
#include "mnormalize.h"
#include "util.h"
#include "read_beta.h"
#include "gene_expression.h"

#define SAMPLE_CNT 1
/*int main(void)
{
	printf("start LDA\n");
	gene_expre("/home/tutut/MLDA/software/MLDA/TEST_DATA/workfloder","ENSG00000001461",4,3,"/home/tutut/MLDA/software/MLDA/LR0_result");
	return(1);
}*/

int
gene_expre (char* data_path,char* gene_name,int condition, int replicate,char *output_path)
{
	printf("start \t%s\n", gene_name);
	document *data=NULL;
	document *sum_exon=NULL;//creat the space for the sum of read tunnel
	double *alpha;
	double *eta;
	//double *short_alpha;
	double **gammas;
	double *theta;
	double **beta,**betatran;
	double ** isoform;
	double *iso_expre;
	double *iso_var;
	double *isolen;
	double *isonum;
	long *readNum;
	int sumnum,i,j,k;
	double like=0;
	long SEQ_DEPTH;
	int nlex, dlenmax,lda_nclass;
	int nclass = CLASS_DEFAULT;		// default in lda.h

	int beta_size[2];//line and row of beta
	/*input data path*/
	char *gene_data_path=(char*)calloc(256,sizeof(char));
	gene_data_path=strcat( gene_data_path,data_path);
	gene_data_path=strcat( gene_data_path,"/ModelMultiGene_Data/");
	gene_data_path=strcat( gene_data_path,gene_name);

	char *gene_norm_data_path=(char*)calloc(256,sizeof(char));
	gene_norm_data_path=strcat(gene_norm_data_path,data_path);
	gene_norm_data_path=strcat(gene_norm_data_path,"/ModelMultiGene_NormData/");
	gene_norm_data_path=strcat(gene_norm_data_path,gene_name);

	char *isoNum_path=(char*)calloc(256,sizeof(char));
	isoNum_path=strcat(isoNum_path,data_path);
	isoNum_path=strcat(isoNum_path,"/isoNum/");
	isoNum_path=strcat(isoNum_path,gene_name);

	char *isoLen_path=(char*)calloc(256,sizeof(char));
	isoLen_path=strcat(isoLen_path,data_path);
	isoLen_path=strcat(isoLen_path,"/isoLen/");
	isoLen_path=strcat(isoLen_path,gene_name);

	char *map_path=(char*)calloc(256,sizeof(char));
	map_path=strcat(map_path,data_path);
	map_path=strcat(map_path,"/ModelMultiGene_Map/");
	map_path=strcat(map_path,gene_name);

	char *readNum_path=(char*)calloc(256,sizeof(char));
	readNum_path=strcat(readNum_path,data_path);
	readNum_path=strcat(readNum_path,"/seq_depth");

		/*output file path*/
	char *output_data_path=(char*)calloc(256,sizeof(char));
	output_data_path=strcat(output_data_path,output_path);

	beta_size_cnt(map_path,beta_size);
	nlex=beta_size[0];//THE COW OF BETA
	nclass=beta_size[1];//THE LINE OF BETA
	//lda_nclass=nclass+1;//LDA模型计算添加一个主题。
	lda_nclass = nclass;

	int emmax = EMMAX_DEFAULT;		// default in lda.h
	int demmax = DEMMAX_DEFAULT;	// default in lda.h
	double epsilon = EPSILON_DEFAULT;	// default in lda.h

	/* open data */
	if ((data = feature_matrix(gene_data_path, condition, replicate, &nlex, &dlenmax)) == NULL) {
		fprintf(stderr, "lda:: cannot open training data.\n");
		exit(1);
	}
	/* create the space for the sum of read*/
	 if ((sum_exon = (document*)calloc(1, sizeof(document))) == NULL)
	{
		fprintf(stderr,"gene_expression:: cannot allocate the total of read.\n");
		exit(1);
	}
	/* allocate parameters */
	if ((alpha = (double *)calloc(lda_nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda:: cannot allocate alpha.\n");
		exit(1);
	}
	if ((readNum = (long *)calloc(condition, sizeof(long))) == NULL) {
		fprintf(stderr, "lda:: cannot allocate readNum.\n");
		exit(1);
	}
	if ((eta = (double *)calloc(nlex, sizeof(double))) == NULL) {
		fprintf(stderr, "lda:: cannot allocate eta.\n");
		exit(1);
	}
	/*LDA截取前面真实读段数据*/
	/*if ((short_alpha = (double *)calloc(nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda:: cannot allocate alpha.\n");
		exit(1);
	}*/
	if ((beta = dmatrix(nlex, lda_nclass)) == NULL) {
		fprintf(stderr, "lda:: cannot allocate beta.\n");
		exit(1);
	}
	if ((betatran = dmatrix(lda_nclass,nlex)) == NULL) {
		fprintf(stderr, "lda:: cannot allocate beta.\n");
		exit(1);
	}
	if ((gammas = dmatrix(condition, lda_nclass)) == NULL) {
		fprintf(stderr, "lda:: cannot allocate gammas.\n");
		exit(1);
	}

	if((theta=(double *)calloc(nclass,sizeof(double)))==NULL)
	{
		fprintf(stderr, "lda:: cannot allocate theta.\n");
		exit(1);
	}
	if((isoform=dmatrix(nclass,nlex))==NULL)
	{
		fprintf(stderr, "lda:: cannot allocate isoform.\n");
		exit(1);
	}
	if((iso_expre=(double *)calloc(nclass,sizeof(double)))==NULL)//申请基因异构体表达值空间
	{
		fprintf(stderr, "lda:: cannot allocate iso_expre.\n");
		exit(1);
	}
	if((iso_var=(double *)calloc(nclass,sizeof(double)))==NULL)//申请基因异构体表达值空间
	{
		fprintf(stderr, "lda:: cannot allocate iso_var.\n");
		exit(1);
	}
	if((isonum=(double*)calloc(nclass,sizeof(double)))==NULL)
	{
		fprintf(stderr, "lda:: cannot allocate isonum.\n");
		exit(1);
	}
	if((isolen=(double*)calloc(nclass,sizeof(double)))==NULL)//异构体长度
	{
		fprintf(stderr, "lda:: cannot allocate isolen.\n");
		exit(1);
	}

	printf("Allocate space end!\n");
	read_isonum(isoNum_path,isonum,nclass);
	read_beta(map_path,beta);
	for(i =0; i<nlex; i++)
		for(j =0; j<nclass; j++)
			betatran[j][i] = beta[i][j];
	//normal_beta(beta,lda_nclass,nlex);
	read_readnum(readNum_path, readNum, condition,replicate);
	sumnum = sum_read_exon(data,sum_exon,nlex);//计算所有通道读段总数
	read_isolen(isoLen_path,isolen,nclass);
	double **iso_expre_sum;
	double *gene_expre_sum;
	if((iso_expre_sum=dmatrix(condition*replicate,nclass))==NULL)
	{
		fprintf(stderr, "lda:: cannot allocate iso_expre_sum.\n");
		exit(1);
	}
	if((gene_expre_sum=(double *)calloc(condition*replicate,sizeof(double)))==NULL)//申请基因异构体表达值空间
	{
		fprintf(stderr, "lda:: cannot allocate gene_expre_sum.\n");
		exit(1);
	}
	if (nclass>1 && sumnum>30)
	{
		/*lda learn */
		int t_cnt;

		like = lda_learn (data, alpha, eta, gammas,beta,isonum, condition,replicate,nclass, nlex, dlenmax,
			emmax, demmax, epsilon);
		printf("LDA end!\n");
		int cnt=0;
		double _sum=0;
		double _sumgamma [condition];
		double _sumeta=0;
		//alpha 归一化
		for(cnt=0;cnt<nclass;cnt++)
		{
			_sum=alpha[cnt]+_sum;
		}
		for(cnt=0;cnt<nclass;cnt++)
		{
			alpha[cnt] = alpha[cnt]/_sum;
		}
		for(i =0; i<condition; i++)
		{
			_sumgamma[i] = 0;
			for(j = 0; j<nclass; j++)
				_sumgamma[i]+=gammas[i][j];
		}
		for(i =0; i<condition; i++)
			for(j =0; j<nclass; j++)
				gammas[i][j] = gammas[i][j]/_sumgamma[i];
		for(i =0; i<nlex; i++)
			_sumeta+=eta[i];
		for(i =0; i<nlex; i++)
			eta[i] =eta[i]/_sumeta;
		/*calculate the gene expression*/
		free_feature_matrix(data);
		if ((data = feature_matrix1(gene_data_path, condition, replicate, &nlex, &dlenmax)) == NULL) {
			fprintf(stderr, "lda:: cannot open training data.\n");
			exit(1);
		}

		sum_read_exon(data,sum_exon,nlex);
		printf("Start calculate gene expression!\n");
		for (i = 0; i<condition*replicate; i++)
		{
			int index;
			index = i/replicate;
			for (k =0; k<nclass; k++)
			{
				theta[k] = gammas[index][k];
			}
			SEQ_DEPTH = readNum[index];
			//printf("%ld\n", SEQ_DEPTH);
			double gene_expression=0;
			read_mapping(isoform,data[i],theta,beta,nclass,nlex);
			expr_sig(isoform,nclass,nlex,isolen,iso_expre,SEQ_DEPTH);
			for(t_cnt=0;t_cnt<nclass;t_cnt++)
			{
				gene_expression=gene_expression+iso_expre[t_cnt];//异构体表达值相加求基因表达值
				iso_expre_sum[i][t_cnt] = iso_expre[t_cnt];
			}
			gene_expre_sum[i] = gene_expression;
		}
		printf("Start writing result!\n");
		write_OUTCOME(output_data_path,alpha,eta,gammas,iso_expre_sum,gene_expre_sum,condition,replicate,nclass,nlex,gene_name,like);
		free_feature_matrix(data);
	}
	if (sumnum>30 && nclass==1)
	{
		int t_cnt,index;
		free_feature_matrix(data);
		if ((data = feature_matrix1(gene_data_path, condition, replicate, &nlex, &dlenmax)) == NULL) {
			fprintf(stderr, "lda:: cannot open training data.\n");
			exit(1);
		}
		sum_read_exon(data,sum_exon,nlex);
		printf("Start calculate gene expression!\n");
		for(i =0; i<nlex; i++)
			eta[i] = (float)1.0/i;
		alpha[0] = 1;
		for(i =0 ;i<condition; i++)
			gammas[i][0]=1;
		for(i =0; i<condition*replicate; i++)
		{
			theta[0]=1;
			index = i/replicate;
			SEQ_DEPTH = readNum[index];
			double gene_expression=0;
			read_mapping(isoform,data[i],theta,beta,nclass,nlex);
			expr_sig(isoform,nclass,nlex,isolen,iso_expre,SEQ_DEPTH);
			for(t_cnt=0;t_cnt<nclass;t_cnt++)
			{
				gene_expression=gene_expression+iso_expre[t_cnt];//异构体表达值相加求基因表达值
				iso_expre_sum[i][t_cnt] = iso_expre[t_cnt];
			}
			gene_expre_sum[i] = gene_expression;
		}
		printf("Start writing result!\n");
		write_OUTCOME(output_data_path,alpha,eta,gammas,iso_expre_sum,gene_expre_sum,condition,replicate,nclass,nlex,gene_name,like);
		free_feature_matrix(data);
	}

	free(map_path);
	free(gene_data_path);
	free(gene_norm_data_path);
	free(isoLen_path);
	free(output_data_path);
	free(isoNum_path);

	free(iso_expre);
	free(iso_var);
	free(isolen);
	free(isonum);
	free(readNum);

	free_sum_exon(sum_exon);
	free_dmatrix(beta, nlex);
	free(theta);
	free_dmatrix(isoform,nclass);
	free_dmatrix(iso_expre_sum,condition*replicate);
	free(gene_expre_sum);
	free(alpha);
	free(eta);
	return(0);
}



