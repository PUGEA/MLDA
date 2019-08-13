#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "gene_expression.h"

#define SAMPLE_CNT 1
//#define SEQ_DEPTH   100//测序深度
int
create_sampling(double *alpha,double *theta,int nclass)
{
	//int n;
	//n=SAMPLE_CNT;//dirchlet sampling count
	/* set up GSL RNG */
	//!const gsl_rng_type *T;
	//!gsl_rng *r ;
	//!gsl_rng_env_setup();
	//!T=gsl_rng_default;
	//!r=gsl_rng_alloc(T);
	/* end of GSL setup */
	//!gsl_ran_dirichlet(r,nclass,alpha,theta);
	//!gsl_rng_free(r);
	return(0);
}
/*map the read to the isoform*/
int
read_mapping(double **isoform,document sum_exon,double *theta,double **Map,double nclass,double nlex)
{
	int t_nlex,t_nclass;
	double* t_scale=(double*)calloc(nclass,sizeof(double));//申请存放临时比例空间
	int i,j;
	for(t_nlex=0;t_nlex<nlex;t_nlex++)
	{
		//根据MAP排除不分配的基因项
		double t_cnt=0;
		int _zero_cnt=0;
		for(t_nclass=0;t_nclass<nclass;t_nclass++)
		{
			if(Map[t_nlex][t_nclass]!=0.0)
			{
				t_scale[t_nclass]=theta[t_nclass];
				t_cnt=t_cnt+theta[t_nclass];
				_zero_cnt++;
			}
			else
			{
				t_scale[t_nclass]=0;
			}
		}
		//计算新的分配比例
		if(t_cnt==0)
		{
			for(t_nclass=0;t_nclass<nclass;t_nclass++)
			{
				if(t_scale[t_nclass]!=0.0)
				{
					isoform[t_nclass][t_nlex]=(double)((sum_exon.cnt[t_nlex]))/(double)(_zero_cnt);
				}
				else
				{
					isoform[t_nclass][t_nlex]=(double)0.0;
				}
			}
		}
		else
		{
			for(t_nclass=0;t_nclass<nclass;t_nclass++)
			{
				isoform[t_nclass][t_nlex]=(double)((sum_exon.cnt[t_nlex]))*((double)(t_scale[t_nclass]/(double)t_cnt));
			}
		}
	}
	free(t_scale);
	return(0);
}
int
read_isolen(char *isolen_name,double* isolen,char t_lex)
{
	FILE* input;
	char line[16];
	int data_cnt=0;
	if((input=fopen(isolen_name,"r"))==NULL)
		return 0;
	while(!feof(input))//read data until the end of the document
	{
		if((fgets(line,sizeof(line),input))!=NULL)
		{
			isolen[data_cnt]=atof(line);
			//printf(line);
		}
		data_cnt++;
		if(data_cnt==t_lex)
		{
			fclose(input);
			return (1);
		}
	}
	fclose(input);
   //printf("theme count%d\n",theme_cnt);
    //printf("word count%d\n",word_cnt);
	return(1);
}
double
sum_gene_exon(double *p_exon,int len)
{
	double sum_exon=0;
	int t_cnt=0;
	for(t_cnt=0;t_cnt<len;t_cnt++)
	{
		sum_exon=sum_exon+p_exon[t_cnt];
	}
	return sum_exon;
}
void
expr_sig(double ** isoform,int nclass,int nlex, double* isolen,double *iso_expre,int SEQ_DEPTH)
{
	double total_exon=0;
	int t_cnt=0;
	for(t_cnt=0;t_cnt<nclass;t_cnt ++)
	{
		//isolen[t_cnt]=100;
		total_exon=sum_gene_exon(isoform[t_cnt],nlex);
		iso_expre[t_cnt]=(total_exon)*1.0e9/(SEQ_DEPTH*isolen[t_cnt]);
	}
}
int
read_isonum(char *isonum_name,double* isonum,char t_lex)
{
	FILE* input;
	char line[16];
	int data_cnt=0;
	if((input=fopen(isonum_name,"r"))==NULL)
		return 0;
	while(!feof(input))//read data until the end of the document
	{
		if((fgets(line,sizeof(line),input))!=NULL)
		{
			isonum[data_cnt]=atoi(line);
			//printf(line);
		}
		data_cnt++;
		if(data_cnt==t_lex)
		{
			fclose(input);
			return (1);
		}
	}
	fclose(input);
   //printf("theme count%d\n",theme_cnt);
    //printf("word count%d\n",word_cnt);
	return(1);
}
int read_readnum(char *readNum_name, long *readNum, int condition,int replicate)
{
	FILE* input;
	char line[100];
	int data_cnt=0;
	int index;
	int i=0;
	for(i=0; i<condition; i++)
		readNum[i] = 0;
	if((input=fopen(readNum_name,"r"))==NULL)
		return 0;
	while(!feof(input))//read data until the end of the document
	{
		if((fgets(line,sizeof(line),input))!=NULL)
		{
			index = data_cnt/replicate;
			readNum[index]=readNum[index] + atol(line);
		}
		data_cnt++;
		if(data_cnt==condition*replicate)
		{
			fclose(input);
			return (1);
		}
	}
	fclose(input);
   //printf("theme count%d\n",theme_cnt);
    //printf("word count%d\n",word_cnt);
	return(1);
}
