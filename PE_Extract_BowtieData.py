#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  


import os
import os.path

#提取每个通道下每个基因对应的读段信息，一个基因一个文件，里面含有对应到多异构体的读段
def ExtractGeneBowtieData(InputFile,TargetGeneFile,OutputFile,readNumFile):
    print 'Processing and Reassigning Read for'+ OutputFile
    f = open(InputFile,'r')
    ReadNo = {}
    ReadGene = {}

    for line in f:                   # 统计每个读段映射的基因个数和名称
      line = line.rstrip()
      line = line.split('\t')
      if int (len(line))==6:
        read = line[0]
        gene = line[5]
    
        if read not in ReadNo:
           ReadNo[read]=1
           genelist=[gene]
           ReadGene[read] = genelist
        else:
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
    
    ReadLen = len(ReadNo)           # 统计每个基因上唯一映射的读段个数
    ReadList = ReadNo.items()
    #print 'delete ReadNo...'
    del ReadNo                    # 删除ReadNo字典  
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



#    fo = open('./ReadRate/'+OutputFile+'.GeneNo','w')
#    GeneNoList = GeneNo.items()
#    for i in range(len(GeneNoList)):
#        fo.write(GeneNoList[i][0]+'\t'+str(GeneNoList[i][1])+'\n')
#    fo.close()


    print 'Calculate Rate of Read...'
    ReadRate = {}                    
# 计算每个读段的分配概率，根据其对应基因的上的唯一映射读段数目，如果对应基因上的唯一映射读段都是0，这此读段平均分配。
#    fo = open('./ReadRate/'+OutputFile+'.ReadRate','w')
    #print ReadLen
    for i in range(ReadLen):
        
        read = ReadList[i][0]
        readNo = ReadList[i][1]
        gene = ReadGene[read]


        if len(gene) ==1:           # 对于唯一映射的读段，即只对应一个基因，其概率为1
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
           ReadRate[read] = gene,geneRate
    del ReadList                   # 删除ReadList list数据
    del ReadGene                  # 删除ReadGene字典
    del GeneNo 

    f_ref = open(TargetGeneFile,'r')        # 根据注释文件中需求基因来提起此基因的读段数据
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
    print 'Extracting bowtie data....' 
    f = open(InputFile,'r')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        
        read,iso,loc1,loc2,readlen,gene = line;
        #print line
        if gene in dic_ref:
           genelist,geneRate = ReadRate[read]
           index=genelist.index(gene)
           rate = geneRate[index]
           isono = int(dic_ref[gene])          # 输出注释文件中的基因读段数据，包括读段名，相对位置，isoform名，基因名，分配概率
           tem=[read,loc1,loc2,iso,rate,readlen]
           if gene not in dicGeneout:
                dicGeneout[gene]=[]
                dicGeneout[gene].append(tem)
           else:
                dicGeneout[gene].append(tem)
    f.close()
    for gene_name in dicGeneout.keys():
        f=open(OutputFile+'/'+gene_name,'w')
        for line in dicGeneout[gene_name]:
                for data in line:
                        f.write(str(data)+'\t')
                f.write('\n')
        f.close()
        del dicGeneout[gene_name]
    del dicGeneout
    del dic_ref



