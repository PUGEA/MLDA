import sys
import os
import numpy

def extract(filename):
    num = len(filename)
    print('Start extract sam file:')
    for i in range(num):
        print(filename[i])
        name1 = filename[i]+'.sam'
        name2 = filename[i]
        with open(name1,'rb')as fr:
            with open(name2,'w')as fw:
                for line in fr:
                    line=line.rstrip();
                    line=line.split("\t")
                    read=line[0];
                    iso=line[2];
                    loc=line[3];
                    fw.write(str(read)+'\t'+str(iso)+'\t'+str(loc)+'\n')
