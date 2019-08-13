import os

def getIsonum(refFlat_Check_MultiGeneExonSplitFile,targetDir):
	print'Begin get isoform num...'
	targetDir  = os.path.join(targetDir,'workfloder')
	isoNumFile = os.path.join(targetDir,'isoNum')
	if not os.path.exists(isoNumFile):
		 os.mkdir(isoNumFile)
	with open(refFlat_Check_MultiGeneExonSplitFile) as fr:
		while 1:
			line = fr.readline()
			line = line.strip()
			if(line == ''):
				break
			line = line.split('\t')
			genename = line[0]
			iso = int(line[1])
			if os.path.exists(isoNumFile+'/'+genename):
				os.remove(isoNumFile+'/'+genename)
			fw = open(isoNumFile+'/'+genename,'a')
			for i in range(iso):
				line = fr.readline()
				line = line.strip()
				line = line.split('\t')
				fw.write(line[2]+'\n')
			fw.close()