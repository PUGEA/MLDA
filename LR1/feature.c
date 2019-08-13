/*
    feature.c
    an implementation of feature matrix.
    $Id: feature.c,v 1.10 2019/7/16 02:12:03 dmochiha Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "feature.h"
#define  BUFSIZE  65536
static int n_fields (char *line);
static int file_lines (FILE *fp);
static int isspaces (char *s);

document *
feature_matrix (char *file_name, int condition, int replicate, int *maxid, int *maxlen)
{
	document *d;
	int  n,j, m,i,index,flag=0;
	FILE *fp;
	char line[BUFSIZE];
	char lineNo[10];
	*maxid  = -1;
	*maxlen =  0;

	if ((fp = fopen(file_name, "r")) == NULL)//open the document data
		return NULL;
	//m = file_lines(fp);
	m = condition;
	if ((d = (document *)calloc(m + 1, sizeof(document))) == NULL)//initialise the data space
		return NULL;
	d[m].len = -1;

	n = 0;
	flag = 0;
	while (fgets(lineNo, sizeof(lineNo), fp))//read a line from the document
	{
		index = atoi(lineNo);
		n = index/replicate;
		fgets(line, sizeof(line), fp);
		int i, len;
		char *cp, *sp, *lp = line;

		if (isspaces(line))
			continue;

		len = n_fields(line);
		if (flag==0)
		{
			for(j =0; j<m;j++)
			{
				d[j].len = len;
				d[j].id  = (int *)calloc(len, sizeof(int));
				d[j].cnt = (double *)calloc(len, sizeof(double));
			}
			flag = 1;
		}
		if (len > *maxlen)
			*maxlen = len;
		if (!(len > 0)) {
			fprintf(stderr, "feature_matrix: suspicious line:\n%s",
				line);
			exit(1);
		}
		if ((d[n].id == NULL) || (d[n].cnt == NULL))
		{
			printf( "id or cnt is null\n");
			return NULL;
		}

		i = 0;
		while (*lp)
		{
			int id;
			double cnt;
			if ((cp = strchr(lp, ':')) == NULL)
				break;
			if ((sp = strpbrk(cp + 1, " \t\n")) == NULL)
				break;
			*cp = '\0';
			*sp = '\0';
			id = atoi(lp) - 1;	/* zero origin */
			cnt = atof(cp + 1);
			if (id >= *maxid)
				*maxid = id + 1;
			d[n].id[i]  = id;
			d[n].cnt[i] = d[n].cnt[i]+cnt;
			lp  = sp + 1;
			i++;
		}
		//n++;
	}
	fclose(fp);
	return(d);
}
document *
feature_matrix1 (char *file_name, int condition, int replicate, int *maxid, int *maxlen)
{
	document *d;
	int  n,j,flag=0, m;
	FILE *fp;
	char line[BUFSIZE];
	char lineNo[10];
	*maxid  = -1;
	*maxlen =  0;

	if ((fp = fopen(file_name, "r")) == NULL)//open the document data
		return NULL;
	//m = file_lines(fp);
	m = condition*replicate;
	if ((d = (document *)calloc(m + 1, sizeof(document))) == NULL)//initialise the data space
		return NULL;
	d[m].len = -1;

	n = 0;
	flag =0;
	while (fgets(lineNo, sizeof(lineNo), fp))//read a line from the document
	{
		n = atoi(lineNo);
		fgets(line, sizeof(line), fp);
		int i, len;
		char *cp, *sp, *lp = line;

		if (isspaces(line))
			continue;

		len = n_fields(line);
		if (flag==0)
		{
			for(j =0; j<m;j++)
			{
				d[j].len = len;
				d[j].id  = (int *)calloc(len, sizeof(int));
				d[j].cnt = (double *)calloc(len, sizeof(double));
			}
			flag = 1;
		}
		if (len > *maxlen)
			*maxlen = len;
		if (!(len > 0)) {
			fprintf(stderr, "feature_matrix: suspicious line:\n%s",
				line);
			exit(1);
		}
		if ((d[n].id == NULL) || (d[n].cnt == NULL))
		{
			printf( "id or cnt is null\n");
			return NULL;
		}

		i = 0;
		while (*lp)
		{
			int id;
			double cnt;
			if ((cp = strchr(lp, ':')) == NULL)
				break;
			if ((sp = strpbrk(cp + 1, " \t\n")) == NULL)
				break;
			*cp = '\0';
			*sp = '\0';
			id = atoi(lp) - 1;	/* zero origin */
			cnt = atof(cp + 1);
			if (id >= *maxid)
				*maxid = id + 1;
			d[n].id[i]  = id;
			d[n].cnt[i] = cnt;
			lp  = sp + 1;
			i++;
		}
		//n++;
	}
	fclose(fp);
	return(d);
}
double
sum_read_exon(document* data,document *sum_exon,int nlex)//count the total read of every exon
{
	int i=0,j=0;
	int doc_cnt=0;
	int word_cnt=0;
	double sumnum =0 ;
	(*sum_exon).id  = (int *)calloc(nlex, sizeof(int));
	(*sum_exon).cnt = (double *)calloc(nlex, sizeof(double));
	(*sum_exon).len=nlex;
	for(word_cnt=0; word_cnt<nlex; word_cnt++)
	{
		doc_cnt=0;
		while(data[doc_cnt].len!=-1)
		{
			sum_exon->id[word_cnt]=data[doc_cnt].id[word_cnt];
			sum_exon->cnt[word_cnt]=sum_exon->cnt[word_cnt] + data[doc_cnt].cnt[word_cnt];
			sumnum = sumnum+data[doc_cnt].cnt[word_cnt];
			doc_cnt++;
		}
	}
	return(sumnum);
}

void
free_sum_exon(document* sum_exon)
{
	free(sum_exon->id);
	free(sum_exon->cnt);
	free(sum_exon);
}
void
free_feature_matrix (document *matrix)
{
	document *dp;
	for (dp = matrix; (dp->len) != -1; dp++)
	{
		free(dp->id);
		free(dp->cnt);
	}
	free(matrix);
}

static int
file_lines (FILE *fp)
{
	int n = 0;
	char buf[BUFSIZE];

	while (fgets(buf, sizeof(buf), fp))
	{
		if (!isspaces(buf))
			n++;
	}
	rewind(fp);

	return n;
}

static int
n_fields (char *line)//
{
	int n = 0;
	char *cp, *sp, *lp = line;

	while (*lp)
	{
		if ((cp = strchr(lp, ':')) == NULL)
			break;
		if ((sp = strpbrk(cp + 1, " \t\n")) == NULL)
			break;
		lp = sp + 1;
		n++;
	}

	return n;
}

static int
isspaces (char *s)	/* return 1 if s consists of only white spaces */
{
	char *c = s;
	while (*c)
	{
		if (!isspace(*c))
			return 0;
		c++;
	}
	return 1;
}


/* $Id: feature.c,v 1.10 2019/7/16 02:12:03 dmochiha Exp $ */
