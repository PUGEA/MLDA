#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  



import os
import os.path
Sample = ['A']

#将每个通道中的读段和外显子进行比对，得到每个读段对应的外显子信息，一个外显子时加一个，跨外显子时显示两个外显子的链接1-2
def CalculateAbsoluteLocation(SubTargetGeneFile,TargetGeneExonSplitFile,CalLocFileInPath,CalLocFileOutPath,ReadLen):
    readLen=ReadLen
    f_ref = open(TargetGeneExonSplitFile,'r');            
    infoAll={}
    while True:
    	line=f_ref.readline();
    	if not line:
    		break;
    	else:
    		line=line.split('\t')
    		iso_no=line[1]
    		if iso_no.isdigit():
    			gene_info = [[]for i in range(3)];
                        iso_info = [[[]for j in range(3)]for i in range(int(iso_no))]
                        iso_nm_info = []
	    		gene=line[0]

	    		split_no = line[2];

		        split_index = line[3];
		        split_index = split_index.rstrip();
		        split_index = split_index.split(',')

		        split_len = line[4];
		        split_len = split_len.rstrip();
		        split_len = split_len.split(',');

		        gene_info[0].append(int(split_no));

		        for i in range(int(split_no)):
				gene_info[1].append(int(split_index[i]));
				gene_info[2].append(int(split_len[i]));
			gene_info[2].insert(0,0)
			flag=1;
                	while flag != int(iso_no)+1:                        #gene信息  对应的isoforms信息
				        line = f_ref.readline();
				        line = line.rstrip();
				        line = line.split('\t'); 

				        if line[0] != gene:
				            continue;
				        else:                    
				                              
				                iso_nm = line[1]; 
				                split_no = line[2];
		
				                split_index = line[3];
				                split_index = split_index.rstrip();
				                split_index = split_index.split(',')  
		
				                split_len = line[4];
				                split_len = split_len.rstrip();
				                split_len = split_len.split(',');
		
				                iso_nm_info.append(iso_nm);                                
				                iso_info[flag-1][0].append(int(split_no));        
				                    
				                for i in range(int(split_no)):
				                    iso_info[flag-1][1].append(int(split_index[i]));
				                    iso_info[flag-1][2].append(int(split_len[i]));
				                            
				        flag = flag + 1;    
	        	infoAll[gene]=[gene_info,iso_info,iso_nm_info]        #iso_nm_info是gene所有isoforms名字
    f_ref.close();			

    f_info = open(SubTargetGeneFile,'r')  
    for info in f_info:
        info = info.rstrip();
        info = info.split('\t');
        gene = info[0];
        iso_no = int(info[1]);

        #if int(iso_no) ==1:              
         #  pass;
        #else:
        gene_info =infoAll[gene][0]     #   [[split_no],[split_index],[split_len]]     split_len=0,23,34,57....
        iso_info = infoAll[gene][1]     #[ [[split_no],[split_index],[split_len]],[[split_no],[split_index],[split_len]].... ]
        iso_nm_info = infoAll[gene][2]  #  [iso1,iso2....]
               

         

        for sam in range(len(Sample)):             
               if sam == 1: 
                   rep1 = 1
               else:
                   rep1 = 1
               for rep in range(rep1):
                   targetFile = os.path.join(CalLocFileInPath, gene);
            
                   if os.path.isfile(targetFile): 
                     f_gene = open(CalLocFileInPath+'/'+gene,'r');
                     f_out = open(CalLocFileOutPath+'/'+ gene,'w');
    
                     for line in f_gene:
                         line = line.rstrip();
                         line = line.split('\t');
    
                         readid = line[0];
                         loc1 = int(line[1]);
                         iso_na = line[3];
                         prob = line[4] #rate
                         
                         loc2=int(line[2])


                         if iso_na in iso_nm_info:
                            ind = iso_nm_info.index(iso_na);
                            exon_no = iso_info[ind][0][0];    
                            exon_ind = iso_info[ind][1];
                            exon_len=[]
                            exon_len_temp = iso_info[ind][2];
                            for i in range(len(exon_len_temp)):
                                exon_len.append((iso_info[ind][2])[i])
                            #exon_len = iso_info[ind][2];
                            exon_len.insert(0,0)
                           # print exon_len    [0,23,34,57.]
                            #print 'loc'+str(loc)
                            abs_loc1 = 0;
                            relat_loc1 = 0;
                            relat_ind1 = 0;
                            abs_loc2=0
                            relat_loc2=0
                            relat_ind2=0
                            flag1=0;#denote find startExon
                            flag2=0#denote find endExon
                            for i in range(exon_no):
                               if loc1 >= exon_len[i] and loc1<= exon_len[i+1]:#start location   第i个exon上
                                   flag1=1
                                   relat_loc1 = loc1 - exon_len[i];
                                   relat_ind1 = exon_ind[i];  #开始的exon
                                   starExon=exon_ind.index(relat_ind1)
                                 #  print 'relat_ind:'+str(relat_ind)
                                   abs_loc1 = gene_info[2][relat_ind1-1]+relat_loc1;
                                   if  i< exon_no-1:
                                     relat_conjuction=str(exon_ind[i])+'-'+str(exon_ind[i+1])
                                   
                                   break;
                            for i in range(exon_no):
                            	if loc2+readLen> exon_len[i] and loc2+readLen<= exon_len[i+1]:#end location
                                   flag2=1
                                   relat_loc2 = loc2+readLen - exon_len[i];
                                   relat_ind2 = exon_ind[i]; #结束的exon
                                   endExon=exon_ind.index(relat_ind2)
                                 #  print 'relat_ind:'+str(relat_ind)
                                   abs_loc2 = gene_info[2][relat_ind2-1]+relat_loc2;  #另一个末端结束的地址
                                   #if loc+readLen > exon_len[i+1] and i< exon_no-1:
                                    #   relat_ind=str(exon_ind[i])+'-'+str(exon_ind[i+1])
                                   break;
                            if relat_ind1==relat_ind2:
                            	f_out.write(readid+'\t'+iso_na+'\t'+str(abs_loc1)+'\t'+str(relat_ind1)+'\t'+prob+'\n')
                            else:
                            	f_out.write(readid+'\t'+iso_na+'\t'+str(abs_loc1)+'\t'+str(relat_conjuction)+'\t'+prob+'\n')
                            
                         else:
                           pass
    
                     f_gene.close();
                     f_out.close();

                   else:
                       pass;
    f_info.close();






 
