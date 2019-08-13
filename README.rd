#MLDA

#What is MLDA?

MLDA is the software for detecting DTU genes and obtaining gene or isoform expressions under multiple conditions from RNA-seq data given a reference transcriptome,  The program obtains results using the alignment from Bowtie 2.

The software is an open-source software. You can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation.

#How MLDA works?

##Features

>* MLDA detects DTU genes under mulptiple conditions by likelihood ratio test based on maximum likelihood of two models(LR0 and LR1).
>* MLDA provides an approach to accurately estimate gene and isoform expression from RNA-Seq data by modeling the isoform- and exon-specific read sequencing biases.
>* MLDA provides an approach to accurately estimate isoform abundances under different conditions.

#Contact
Users can ask technical questions by sending emails to Xuejun Liu (xuejun.liu@nuaa.edu.cn).

#Installation

##Installation requirements:

>* Operating system :
	Linux( Ubuntu13.10 or Fedora 17…)
	Mac OS X (Mac OS X 10.8 5 or higher)
>* Python 2.7 
>* GUN Scientific Library GCC 4.8.1 or higher(http://www.gnu.org/software/gsl/)
>* Python package: parallel python (pp)(http://www.parallelpython.com/) 
>* Bowtie2

###NOTE: 
 ***PP*** is a python module which provides mechanism for parallel execution of python code on SMP(systems with multiple processors or cores) and clusters (computers connected via network). It is light, easy to install and integrate with other python software. PP is an open-source and cross-platform module written in pure python.

The software depends on the free and open-source software, the ***GNU Scientific Library (GSL)*** (http://www.gnu.org/software/gsl/), so the GSL needs to be installed on the user’s system. NLDMseq requires the GSL installed in /usr/local, which is the default location. The GSL can be compiled by the user. Users who are compiling NLDMseq from the source code should install GSL in the standard location (/usr/local). If you are not sure if GSL is already installed, at the Terminal prompt $ type:`$gsl-config --cflags --libs-without-cblas`
If GSL is installed, the command above should return the following information:
```shell
-I/sw/include
-L/sw/lib -lgsl –lm
```
GSL can be found in the gsl subdirectory on your nearest GNU mirror ( http://ftpmirror.gnu.org/gsl/) or the Main GNU ftp site (ftp://ftp.gnu.org/gnu/gsl/). The users can download gsl-1.6.tar.gz or a higher version for NLDMseq. Please follow the instructions in the included file “INSTALL” to guide  the installation of this library.

We recommend using the Linux operating system. 
##Installation

>* First, go to the location of NLDMseq

```
    $cd ~
```

>* Second, install the software by the install.sh shell script.

```
$./install_LR0.sh
$./install_LR1.sh
```

NOTE:The user may have to give the script execute permissions by `chmod u+x install_LR0.sh and install_LR1.sh`. This script will use the gcc complier to build C codes in the software.
>* Test the software

```shell
$python test.py
```

To test this software, the user can run the python script *test.py* .This will check whether the software has been built and installed correctly. It will detect DTU gene and calculate the expression of genes and isoforms using the data from the TEST_DATA folder.

##Running MLDA
 run example to test:
$ bash run_example.sh

Use the example(Single-End Reads/Paired-End Reads) from aligning sequenced reads with Bowtie 2 to test this software, the user can run the script *run_example.sh* .This will check whether the software can run on data of reads under multiple conditions. It will detect DTU gene and calculate the expression of genes and isoforms using the data from the bowtie folder.

>* Step 1. Aligning sequenced reads with Bowtie 2
The following commands are used to align sequenced reads to a reference transcriptome.
```shell
$ bowtie2-build -f ref_transcript. Fasta ref_transcript.index
$ bowtie2 –t –f -k 20 -p 4 --no-hd --no-unal -x ref_transcript.index raw_data.fasta -S example.sam
```

If the paired-end reads are processed, the Bowtie command should be listed as

```shell
$ bowtie2 –t –f -k 20 -p 4 --no-hd --no-unal --no-mixed --no-discordant -x ref_transcript.index -1 raw_data.fasta -2 raw_data.fasta -S example.sam
```

The above transcriptome reference sequence can be downloaded from UCSC or Ensembl website.

>* Step 2. Detect DTU gene and calculate expression values with MLDA

The software uses the alignment SAM file from Bowtie 2 as input data.

```
$ python MLDA_calculation.py -s single -c condition -r replicate -o output_path -g GeneName_path -b cutoff -i align_file(s) -t AnnotationFile Type -a AnnotationFile.txt
```

Options:
>* -h/--Help #get help info
>* -s/--SingleEnd: Input data are Single-end reads alignment. (Default: off)
>* -c/--Condition: <int> The value of condition
>* -r/--Replicate: <int > The value of replicate
>* -o/--OutputPath: All output files are in this direction
>* -g/--GeneFile: The path of genename file including all genenames for detecting DTU, the genename is the beginning of each line.(format:"ENSG00000000003\t...\n")
>* -b/--Bound: <float>  The cutoff of pvalue(qvalue) for DTU gene(e.g,0.05)'
>* -i/--Input: File Path+Prefix name of alignment File(s) under all conditions,use comma separating the input files if more than one files are provided.(e.g,/home/tutut/example1,/home/tutut/example2)'
>* -t/--AnnotationType   <int>    supported four annotation types: refGene: 1, ensGene: 2, knownGene: 3 and Ensembl: 4'
>* -a/--AnnotationFile: The path of annotation file, the file includes the gene and isoform information. eg: refGene, knownGene, and ensGene, which all can be download UCSC website. If you use the Ensembl dataset, you may use the .gtf file.'

Output files:
The MLDA produces qvalue of each gene for detecting DTU and gene/isoform expression output files, which are put under the path given in the option -o/--OutputPaht


Description of output files:
1 LR0_result/LR1_result:
(1) alpha: Dirichlet distribution hyperparameter of isoform abundances
(2) eta: Dirichlet distribution hyperparameter of exon abundances
(3) gammas: isoform abundances under all conditions
(4) geneExpre: FPKM expression value of gene under all samples
(5) isoExpre: FPKM expression value of isoform under all samples
(6) likelihood: Maximum log likelihood of each gene
2 LRT_result
(1) LR_value:Chi-square value and isoform num of each gene
(2) DTU_value:p_value and q_value of each gene for detecting DTU


#Example
Here, we use a simple example to show the usage of MLDA. The alignment files from Bowtie 2 and the annotation file are also supplied.

>* Annotation file: Homo_sapiens.GRCh37.71.gtf
>* Alignment from Bowtie2: example1.sam...example8.sam (single-end,4 condition,2 replicate)
>* Alignment from Bowtie2: example1.sam...example6.sam (paired-end,3 condition, 2 replicate)

Since the Bowtie2 output has been supplied, so you can skip step1 and just run the following command.

For single-end data:
```
$ python MLDA_calculation.py -s -c 4 -r 1 -o /home/tutut/MLDA/software/MLDA/Result \
-g /home/tutut/MLDA/software/MLDA/geneFileS -b 0.05 \
-i /home/tutut/MLDA/software/bowtie/SE/sam/example1,\
/home/tutut/MLDA/software/bowtie/SE/sam/example4,\
/home/tutut/MLDA/software/bowtie/SE/sam/example7,\
/home/tutut/MLDA/software/bowtie/SE/sam/example10 \
-t 4 -a /home/tutut/MLDA/software/bowtie/Homo_sapiens.GRCh37.71.gtf
```
After running the above command, you will obtain three output folders of this single-end data, ***LR0_result*** and ***LR1_result*** ***LRT_result*** under -o path.

For paired-end data :
```
$ $ python MLDA_calculation.py -c 3 -r 2 -o /home/tutut/MLDA/software/MLDA/Result \
-g /home/tutut/MLDA/software/MLDA/geneFileP -b 0.05 \
-i /home/tutut/MLDA/software/bowtie/PE/sam/example1,\
/home/tutut/MLDA/software/bowtie/PE/sam/example2,\
/home/tutut/MLDA/software/bowtie/PE/sam/example3,\
/home/tutut/MLDA/software/bowtie/PE/sam/example4,\
/home/tutut/MLDA/software/bowtie/PE/sam/example5,\
/home/tutut/MLDA/software/bowtie/PE/sam/example6 \
-t 4 -a /home/tutut/MLDA/software/bowtie/Homo_sapiens.GRCh37.71.gtf
```
After running the above command, you will obtain three output folders of this paired-end data, ***LR0_result*** and ***LR1_result*** ***LRT_result*** under -o path.

#Authors

Xuejun Liu (xuejun.liu@nuaa.edu.cn).College of Computer Science and Technology, Nanjing University of Aeronautics and Astronautics, 29Jiangjun Rd., Jiangning District, 211106Nanjing China.

The MLDA algorithm is developed by Xuejun Liu and Jing Li. The MLDA software is mainly implemented by Jing Li.

#Related software

Bowtie2: for alignment of RNA-Seq reads to transcriptome (optionally with the precomputed set of splice junctions).

