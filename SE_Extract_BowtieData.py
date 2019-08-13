#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  


import os
import os.path

#在读段后面添加读段所在基因的比例，忽略映射到多基因的读段，概率为1.0
def ExtractGeneBowtieData(InputFile,TargetGeneFile,OutputFile,readNumFile):
    print 'Processing and Reassigning Read for'+ OutputFile
    f = open(InputFile,'r')
    ReadNo = {}
    ReadGene = {}

    for line in f:                   
      line = line.rstrip()
      line = line.split('\t')
      if int (len(line))==4:
        read = line[0]
        gene = line[3]
        if read not in ReadNo:
           ReadNo[read]=1
           genelist=[gene]
           ReadGene[read] = genelist
        else:
            genelist = ReadGene[read]
            if gene in genelist:
                pass
            else:
              ReadNo[read]=ReadNo[read]+1
              genelist=ReadGene[read]
              genelist.append(gene)
              ReadGene[read] = genelist
    f.close()
    seq_depth=open(readNumFile,'a')
    print 'ReadNo.' + str(len(ReadNo))
    seq_depth.write(str(len(ReadNo))+'\n')

    ReadLen = len(ReadNo)          
    ReadList = ReadNo.items()
    #print 'delete ReadNo...'
    del ReadNo 

    GeneNo = {}
    for i in range(ReadLen):
        read = ReadList[i][0]
        readNo = ReadList[i][1]
        
        if readNo ==1:
           gene = ReadGene[read][0]
           if gene not in GeneNo:
              GeneNo[gene] = 1
           else:
              GeneNo[gene] = GeneNo[gene]+1
        else:
           #pass
            for tt in range(readNo):
                gene = ReadGene[read][tt]
                if gene not in GeneNo:
                    GeneNo[gene] = 1.0/readNo
                else:
                    GeneNo[gene] = GeneNo[gene]+1.0/readNo
    print 'GeneNo.'+str(len(GeneNo))

    print 'Calculate Rate of Read...'
    ReadRate = {}

    for i in range(ReadLen):
        read = ReadList[i][0]
        readNo = ReadList[i][1]
        gene = ReadGene[read]


        if len(gene) ==1:          
           geneRate =[1.0]
           ReadRate[read] = gene,geneRate  
           geneNo = GeneNo[gene[0]]
#           fo.write(read+'\t'+str(gene)+'\t'+str(geneNo)+'\t'+str(geneRate)+'\n')    
        else:
           geneNo = []
           for j in range(len(gene)):
               genetemp = gene[j]
               if genetemp in GeneNo:
                  tempNo = GeneNo[genetemp]
               else:
                  tempNo = 0
               geneNo.append(tempNo)
       
           sumNo = sum(geneNo)
           geneRate = []
           if sumNo == 0:
              geneRate = [1/(len(geneNo)) for i in range(len(geneNo))]
           else:
               for j in range(len(gene)):
                   tempRate = geneNo[j]/sumNo 
                   geneRate.append(tempRate) 

#           print read,gene,geneNo,geneRate
#           fo.write(read+'\t'+str(gene)+'\t'+str(geneNo)+'\t'+str(geneRate)+'\n')    
           ReadRate[read] = gene,geneRate
           
           
#    fo.close()

    #print len(ReadRate)
    #print 'Delete ReadList, ReadGene and GeneNo...'
    del ReadList                  
    del ReadGene                 
    del GeneNo 
    
    f_ref = open(TargetGeneFile,'r')       
    dic_ref = {}
    for line in f_ref:
        line = line.rstrip();
        line = line.split('\t');
        gene = line[0];
        isono = line[1]
        if gene not in dic_ref:
           dic_ref[gene] = isono;
        else:
           pass;
    f_ref.close()
    dicGeneout={}
    print 'Beginning extracting bowtie data....' 
    f = open(InputFile,'r')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')	
        if len(line)==4:
        	read, iso, loc, gene = line;
        #print line
        
        	if gene in dic_ref:
           		genelist,geneRate = ReadRate[read]
           		index=genelist.index(gene)
           		rate = geneRate[index]
#           print genelist,geneRate,index,rate
           		isono = int(dic_ref[gene])          
          # if isono ==1 :
          #   pass;
          # else:
          		tem=[read,iso,loc,gene,rate]
          		if gene not in dicGeneout:
          			dicGeneout[gene]=[]
          			dicGeneout[gene].append(tem)
          		else:
          			dicGeneout[gene].append(tem)

#           		f_output = open(OutputFile+'/'+gene,'a')
#           		f_output.write(read+'\t'+loc+'\t'+iso+'\t'+gene+'\t'+str(rate)+'\n')
#           		f_output.close();
    f.close();
    for gene_name in dicGeneout.keys():#a.keys得到字典a中键的list
        f=open(OutputFile+'/'+gene_name,'w')
        for line in dicGeneout[gene_name]:
                for data in line:
                        f.write(str(data)+'\t')
                f.write('\n')
        f.close()
	del dicGeneout[gene_name]
    del dicGeneout

