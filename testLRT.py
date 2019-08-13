import os
from scipy import stats
import statsmodels.stats.multitest as smt

def testLRT(LR_Path,Gene_File,condition,cutoff):
	LR_File = os.path.join(LR_Path,'LR_value')
	DTU_File = os.path.join(LR_Path,'DTU_value')
	if os.path.exists(DTU_File):
		os.remove(DTU_File)
	genelist = []
	with open(Gene_File,'r') as fgene:
		for line in fgene:
			line = line.strip()
			line = line.split('\t')
			genename = line[0]
			if genename not in genelist:
				genelist.append(genename)
	geneFDR = []
	pvalueFDR = []
	qvalueFDR = []
	reject = []

	with open(LR_File,'r') as flr:
		for line in flr:
			line = line.strip()
			line = line.split('\t')
			genename = line[0]
			if genename in genelist:
				isonum = int(line[1])
				LR = float(line[2])
				freedom = (condition-1)*isonum
				p_value = 1.0 - stats.chi2.cdf(LR, freedom)
				geneFDR.append(genename)
				pvalueFDR.append(p_value)
	print(len(pvalueFDR))
	reject, qvalueFDR,c,d = smt.multipletests(pvalueFDR,0.05,'fdr_bh')
	with open(DTU_File,'a') as fdtu:
		fdtu.write('GeneName\tp_value\tq_value\tDTU\n')
		for i in range(len(geneFDR)):
			genename = geneFDR[i]
			pvalue = pvalueFDR[i]
			qvalue = qvalueFDR[i]
			flag = 'FALSE'
			if qvalue<cutoff:
				flag = 'TRUE'
			else:
				flag = 'FALSE'
			fdtu.write(genename+'\t'+str(pvalue)+'\t'+str(qvalue)+'\t'+flag+'\n')
#testLRT('/home/tutut/MLDA/software/MLDA/Result/LRT_result','/home/tutut/MLDA/software/MLDA/Result/Annotation/Ensembl.Gene.Info',4,0.05)
