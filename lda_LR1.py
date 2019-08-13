#!/usr/bin/env python
#coding=utf-8
import os
import os.path
import shutil
import glob
import ctypes
import multiprocessing as mp
########################################
#function:copy file from one dir to another
#input: res file path, target file path
#output:none
########################################
def copy_file(res_file,tar_file):
	shutil.copy(res_file,tar_file)
########################################
#function:remove all file in the dir and the floder itself
#input: dir to remove
#output:none
########################################
def removeFileInFirstDir(targetDir): 
	for file_path in os.listdir(targetDir): 
		targetFile = os.path.join(targetDir,file_path) 
		if os.path.isfile(targetFile): 
			os.remove(targetFile)

########################################
#function: calcultate the abundances of gene and isoform 
#input:EM parement
#output:count of gene
########################################
def pp_lda(gene_list,isoLen_path,Map_path,Data_path,output_path,data_path,condition,replicate):
	for item in gene_list:
		doc_str=item
		isoLen_exists =os.path.exists(os.path.join(isoLen_path,doc_str))
		map_exists =os.path.exists(os.path.join(Map_path,doc_str))
		normdata_exists =os.path.exists(os.path.join(Data_path,doc_str))
		if(isoLen_exists and map_exists and normdata_exists):
			so.gene_expre(data_path,doc_str,condition,replicate,output_path)

def gene_expre(geneFile,data_path,condition,replicate,output_path):
	isoLen_path = os.path.join(data_path,'isoLen')
	Data_path = os.path.join(data_path,'ModelMultiGene_Data')
	Map_path = os.path.join(data_path,'ModelMultiGene_Map')
	NormData_path = os.path.join(data_path,'ModelMultiGene_NormData')
	output_path=os.path.join(output_path,'LR1_result')
	if(os.path.exists(output_path)):
		shutil.rmtree(output_path)
	os.mkdir(output_path)
	res_list = []
	#res_list=os.listdir(Data_path)
	with open(geneFile,'r') as fr:
		for line in fr:
			line = line.strip()
			line = line.split('\t')
			genename = line[0]
			if genename not in res_list:
				res_list.append(genename)
	block_size=1000
	start=0
	end=len(res_list)
	endi=0;
	starti=0;
	pool=mp.Pool(processes=1)
	#pool=mp.Pool(processes=mp.cpu_count()*2)
	while(endi<end):
		endt=starti+block_size;
		endi=min(endt,end);
		pool.apply_async(pp_lda,(res_list[starti:endi],isoLen_path,Map_path,Data_path,output_path,data_path,condition,replicate))
		starti=starti+block_size;
	pool.close()
	pool.join()
########################################
#function:lda main 
#input: data input ,Sequence depth ,expression output path
#output:none
########################################
so=ctypes.CDLL('./LR1/libgene_expre.so')
def lda_expre(geneFile,data_path,condition,replicate,output_path):
	gene_expre(geneFile,data_path,condition,replicate,output_path)
#lda_expre('/home/tutut/MLDA/software/MLDA/geneFileP','/home/tutut/MLDA/software/MLDA/Result/workfloder',3,2,'/home/tutut/MLDA/software/MLDA/Result')



