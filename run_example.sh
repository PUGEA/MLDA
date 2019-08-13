#..........SE................
python MLDA_calculation.py -s -c 4 -r 1 -o /home/tutut/MLDA/MLDA/MLDA/Result \
-g /home/tutut/MLDA/software/MLDA/geneFileS -b 0.05 \
-i /home/tutut/MLDA/software/bowtie/SE/sam/example1,\
/home/tutut/MLDA/software/bowtie/SE/sam/example4,\
/home/tutut/MLDA/software/bowtie/SE/sam/example7,\
/home/tutut/MLDA/software/bowtie/SE/sam/example10 \
-t 4 -a /home/tutut/MLDA/software/bowtie/Homo_sapiens.GRCh37.71.gtf
#........PE................
#python MLDA_calculation.py -c 3 -r 2 -o /home/tutut/MLDA/software/MLDA/Result \
#-g /home/tutut/MLDA/software/MLDA/geneFileP -b 0.05 \
#-i /home/tutut/MLDA/software/bowtie/PE/sam/example1,\
#/home/tutut/MLDA/software/bowtie/PE/sam/example2,\
#/home/tutut/MLDA/software/bowtie/PE/sam/example3,\
#/home/tutut/MLDA/software/bowtie/PE/sam/example4,\
#/home/tutut/MLDA/software/bowtie/PE/sam/example5,\
#/home/tutut/MLDA/software/bowtie/PE/sam/example6 \
#-t 4 -a /home/tutut/MLDA/software/bowtie/Homo_sapiens.GRCh37.71.gtf