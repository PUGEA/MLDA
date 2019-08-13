#! /usr/bin/env python

import os
import os.path
def plusGeneName(refFlat_Check_GeneUnionExonLenFile,readFileList,targetDir,readType,ReadLen): 
	f_isfInfor=open(refFlat_Check_GeneUnionExonLenFile,'r')
	isoList={}
	for line in f_isfInfor: 
		line = line.rstrip();
		line = line.split('\t');
		gene=line[0]
		i=0
		while i<int(line[1]):
			isoList[line[3+i]]=gene
			i=i+1
	f_isfInfor.close();

	targetDir= os.path.join(targetDir,'workfloder')
	readFile = os.path.join(targetDir,'readInput')
	if readType == 'Paired':
		for j in range(len(readFileList)):
			f_reads=open(readFileList[j],'r')  
			f_output = open(readFile+'/'+'Lane'+str(j)+'.plusGene','w')
			lineNum=0;
			while True:
				line=f_reads.readline()
				line=line.strip();
				line = line.split('\t')
				if len(line) == 1:
					f_reads.close()
					f_output.close()
					break
				else:
					lineNum=lineNum+1;
					try:
						isfName1=line[1]
					except IndexError as ie:
						print(j)
						print(line)
						print(lineNum)
					read1=line[0] 
					loc1=line[2]
					if lineNum%2 !=0:
						line=f_reads.readline()
						line=line.rstrip()
						line = line.split('\t')
						lineNum=lineNum+1
						try:
							isfName2=line[1]
						except IndexError as ie:
							print(j)
							print(line)
							print(lineNum)
						read2=line[0] 
						loc2=line[2]
						if isfName1 in isoList:
							gene=isoList[isfName1]
							f_output.write(read1+'\t'+isfName1+'\t'+loc1+'\t'+loc2+'\t'+str(int(loc2)+ReadLen-int(loc1))+'\t'+gene+'\n')
	else:
		for i in range(len(readFileList)):
			f_reads=open(readFileList[i],'r')  
			f_output = open(readFile+'/'+'Lane'+str(i)+'.plusGene','a')
			i=0;
			for line in f_reads:
				line=line.rstrip();
				line = line.split('\t');
				if len(line)==3:
					i=i+1;
					isfName=line[1]
					read=line[0]
					loc=line[2]
					if isfName in isoList:
						gene=isoList[isfName]
						f_output.write(read+'\t'+isfName+'\t'+loc+'\t'+gene+'\t'+'\n')
				else:
					pass
	f_output.close()
	f_reads.close()

