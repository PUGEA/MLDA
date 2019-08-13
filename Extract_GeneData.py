#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  

import PE_Extract_BowtieData
import SE_Extract_BowtieData
import static_exonJunc
import static_readOnExon
import geneExonLen
import getIsoLength
import getIsonum

import plusGeneName
import PE_Extract_CalcAbsLoc_my
import SE_Extract_CalcAbsLoc_my


import os
import os.path

import pp





def startProcess(readFileList,annotation_type,AnnotationDir,ReadLen,readType):
	rep=readFileList
	lanFile=[]
	CatFile = ['ExtractMultiGeneData','ModelMultiGene_AbsLoc','ModelMultiGene_Data','ModelMultiGene_NormData']
	for i in range(len(readFileList)):
		lanFile.append('Lane'+str(i+1))

	targetDir = AnnotationDir
	workfloder = os.path.join(targetDir,'workfloder')
	annotationPath = os.path.join(targetDir,'Annotation')

	####Annocation Files Path 
	GeneExonSplit = os.path.join(annotationPath,annotation_type+'.Exon.Split')
	GeneFile =  os.path.join(annotationPath,annotation_type+'.Gene.Info')
	#####workfloder path
	geneMapout = os.path.join(workfloder,'ModelMultiGene_Map')
	geneExonout = os.path.join(workfloder,'GeneExonFile')
	GeneData_Path = os.path.join(workfloder,'ExtractMultiGeneData')
	GeneAbsLoc_path = os.path.join(workfloder,'ModelMultiGene_AbsLoc')

	exonLen = os.path.join(workfloder,'exonLen')	 
	exonjunfile = os.path.join(workfloder,'GeneExonFile')	 
	spandExonsLen = os.path.join(workfloder,'exonLen')
	readNumFile = os.path.join(workfloder,'seq_depth')

	static_exonJunc.staticGeneNewExon(GeneExonSplit,geneMapout,geneExonout)
	getIsonum.getIsonum(GeneExonSplit,targetDir)
	getIsoLength.getIsoLen(GeneFile,GeneExonSplit,targetDir)#get isoform length
	geneExonLen.staticGeneNewExon(GeneExonSplit,exonjunfile,spandExonsLen,ReadLen)	

	plusGeneName.plusGeneName(GeneFile,readFileList,targetDir,readType,ReadLen)#covert bowtie

	DictIsofGene = {}
	gene_list=[]
	isoNo_list=[]
	length_list=[]
	isoName_list=[]
	f = open(GeneFile,'r')
	for line in f:
		line = line.rstrip()
		line = line.split('\t')
		gene = line[0]
		IsoNo = int(line[1])
		IsoNameList=  line[3:]
		for i in xrange(IsoNo):
			DictIsofGene[IsoNameList[i]]=gene
		gene_list.append(line[0])
		isoNo_list.append(int(line[1]))
		length_list.append(int(line[2]))
		isoName_list.append(line[3:])
	f.close()

	targetDir=AnnotationDir
	LocInputPath = os.path.join(targetDir,'workfloder','ModelMultiGene_AbsLoc')
	exonInputFile = os.path.join(targetDir,'workfloder','GeneExonFile')
	exonLenFile  = os.path.join(targetDir,'workfloder','exonLen')
	DataOutputPath = os.path.join(targetDir,'workfloder','ModelMultiGene_Data')
	NormDataOutput = os.path.join(targetDir,'workfloder','ModelMultiGene_NormData')
	if readType=='Paired':
		for j in range(len(rep)):
			InputPath = os.path.join(workfloder,'readInput')
			InputFile = os.path.join(InputPath,'Lane'+str(j)+'.plusGene')
			OutputFile= os.path.join(GeneData_Path,lanFile[j])
			print InputFile
			print OutputFile
			PE_Extract_BowtieData.ExtractGeneBowtieData(InputFile,GeneFile,OutputFile,readNumFile)
		for l in range(len(rep)):
			CalLocFileInPath = os.path.join(GeneData_Path,lanFile[l])
			CalLocFileOutPath = os.path.join(GeneAbsLoc_path,lanFile[l])
			print CalLocFileInPath
			print CalLocFileOutPath
			PE_Extract_CalcAbsLoc_my.CalculateAbsoluteLocation(GeneFile,GeneExonSplit,CalLocFileInPath,CalLocFileOutPath,ReadLen)

	if readType=='Single':
		for j in range(len(rep)):
			InputPath = os.path.join(workfloder,'readInput')
			InputFile = os.path.join(InputPath,'Lane'+str(j)+'.plusGene')
			OutputFile= os.path.join(GeneData_Path,lanFile[j])
			print InputFile
			print OutputFile
			SE_Extract_BowtieData.ExtractGeneBowtieData(InputFile,GeneFile,OutputFile,readNumFile)
		for l in range(len(rep)):
			CalLocFileInPath = os.path.join(GeneData_Path,lanFile[l])
			CalLocFileOutPath = os.path.join(GeneAbsLoc_path,lanFile[l])
			print CalLocFileInPath
			print CalLocFileOutPath
			SE_Extract_CalcAbsLoc_my.CalculateAbsoluteLocation(GeneFile,GeneExonSplit,CalLocFileInPath,CalLocFileOutPath,ReadLen)
	f_in = open(GeneFile,'r')
	for line in f_in:
		line = line.rstrip();
		line = line.split('\t');
		gene = line[0];
		isoNo = int(line[1]);
		length = int(line[2]);
		isoName = line[3:];
		for j in range(len(rep)):
			LocInputPath = os.path.join(targetDir,'workfloder')
			LocInputFile = os.path.join(LocInputPath,CatFile[1])
			LocInputFile = os.path.join(LocInputFile,lanFile[j])
			exonInputFile=os.path.join(targetDir,'workfloder')
			exonInputFile=os.path.join(exonInputFile,'GeneExonFile')
			exonLenFile = os.path.join(targetDir,'workfloder')
			exonLenFile = os.path.join(exonLenFile,'exonLen')
			DataOutputPath = os.path.join(targetDir,'workfloder')
			DataOutputFile = os.path.join(DataOutputPath,CatFile[2])
			NormDataOutput=os.path.join(DataOutputPath,CatFile[3])
			static_readOnExon.ModelMultiGeneDataScale(j,gene,length,isoNo,isoName,LocInputFile,exonInputFile,exonLenFile,DataOutputFile,NormDataOutput);
	f_in.close()
