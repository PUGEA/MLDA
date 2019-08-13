#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os 
import os.path
import shutil
import ctypes
import getLrt
import testLRT

if not os.path.exists("./LR0/libgene_expre.so"):
    print "Build software filed.Please rebuild again use the install_LR0.sh script"

if os.path.exists("./TEST_DATA/LR0_result"):
    shutil.rmtree("./TEST_DATA/LR0_result")
os.mkdir("./TEST_DATA/LR0_result")
if not os.path.exists("./TEST_DATA/LR0_result"):
        print "Create workfloder failed .Please get the access control"
else:
        so=ctypes.CDLL("./LR0/libgene_expre.so")
        so.gene_expre("./TEST_DATA/workfloder","ENSG00000000419",4,1,"./TEST_DATA/LR0_result")

if os.path.exists("./TEST_DATA/LR1_result"):
    shutil.rmtree("./TEST_DATA/LR1_result")
os.mkdir("./TEST_DATA/LR1_result")
if not os.path.exists("./TEST_DATA/LR1_result"):
        print "Create workfloder failed .Please get the access control"
else:
        so=ctypes.CDLL("./LR1/libgene_expre.so")
        so.gene_expre("./TEST_DATA/workfloder","ENSG00000000419",4,1,"./TEST_DATA/LR1_result")
getLrt.getLrt('./TEST_DATA','Ensembl')
testLRT.testLRT('./TEST_DATA/LRT_result','./TEST_DATA/Annotation/gene',4,0.05)
if os.path.exists('./TEST_DATA/LRT_result/LR_value'):
    fr = open('./TEST_DATA/LRT_result/LR_value','r')
    line = fr.readline()
    line = line.strip()
    line = line.split('\t')
    genename = line[0]
    if genename=='ENSG00000000419':
        print "Test passed! The software has been installed correctly"
    else:
        print "Test filed .Please try again."
else:
    print "Test filed .Please try again."
