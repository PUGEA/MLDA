from __future__ import division
import os
import numpy as np
import math
import shutil

def getLrt(result_path,annotype):
    annoFile = os.path.join(result_path,'Annotation',annotype+'.Gene.Info')
    LR0_file = os.path.join(result_path,'LR0_result','likelihood')
    LR1_file = os.path.join(result_path,'LR1_result','likelihood')
    LRT_path = os.path.join(result_path,'LRT_result')
    if( os.path.exists(LRT_path)):
        shutil.rmtree(LRT_path)
    os.makedirs(LRT_path)
    DicIso={}
    fr = open(annoFile,'r')
    for line in fr:
        line = line.strip()
        line = line.split('\t')
        genename = line[0]
        isonum = int(line[1])
        DicIso[genename] = isonum
    fr.close()
    DicDiff = {}
    DicSame = {}
    fr = open(LR0_file,'r')
    lines = fr.readlines()
    genenum = int(len(lines)/2)
    for i in range(genenum):
        index = 2*i
        genename = lines[index].strip()
        likeli = float(lines[index+1].strip())
        if likeli<0:
            DicDiff[genename] = likeli
    fr.close()
    fr = open(LR1_file,'r')
    lines = fr.readlines()
    genenum = int(len(lines)/2)
    for i in range(genenum):
        index = 2*i
        genename = lines[index].strip()
        likeli = float(lines[index+1].strip())
        if likeli<0:
            DicSame[genename] = likeli
    fr.close()
    LRT_file = os.path.join(LRT_path,'LR_value')
    fw = open(LRT_file,'a')
    fr = open(annoFile,'r')
    for line in fr:
        line = line.strip()
        line = line.split('\t')
        genename = line[0]
        if genename in DicDiff and genename in DicSame:
            L1 = DicDiff[genename]
            L0 = DicSame[genename]
            num = abs(L1-L0)
            LR = abs(2*num)
            isonum = DicIso[genename]
            fw.write(genename+'\t'+str(isonum)+'\t'+str(LR)+'\n')
    fw.close()
    fr.close()
#getLrt('/home/tutut/MLDA/software/MLDA/Result','Ensembl')