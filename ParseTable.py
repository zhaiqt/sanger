# this will parse the output table, containing GERMline, CDR, DNA information of individual reads
import os
import IsolateClone

def ParseTable(filePath):
	AbDict={}
	keyList=['RID','LID','DNAlen','FV-PRO','FV-DNA',"GERMLINE-V",'GERMLINE-D','GERMLINE-J','FR1-PRO','CDR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO','PRO',"DNA",'FR1-PRO','CDR1-DNA','FR2-PRO','CDR2-DNA','FR3-DNA','CDR3-DNA','FR4-DNA','FR4-DNA']

	#print "how many key? " + str(len(keyList))
	count_seq=0
	for filename in os.listdir(filePath):
		if not filename.endswith("all.xls"):
			continue
		currentFile=os.path.join(filePath,filename)
		with open (currentFile) as tableObject:
			for line in tableObject:	
				if line.startswith("#"):
					continue
				count_seq +=1
				line=line.rstrip("\n")
				item=line.split("\t")

				#print "how many item?" + str(len(item))
				AbDict[item[0]]={}
				#print item
				for i in range(0,len(keyList)):
					#print AbDict[item[0]]
					#print keyList[i]
					#print item[i]
					try:
						AbDict[item[0]].update({keyList[i]:item[i+1]})
					except:
						AbDict[item[0]].update({keyList[i]:""})
				AbDict[item[0]]["DNA"] = AbDict[item[0]]["DNA"]
	return (AbDict,count_seq)							

##################
'''
file="/dlab/NGS/usem-seqanalysis/160314_zhaiqi1_miseq_HBx52-60DNA.20160214_AN2N4/test/RESULTS"
AbDict=ParseTable(file)
group = IsolateClone.identifyClone(AbDict)
#print group
#print "Unique comibnation:" + str(len(group))
#print "nubmer of individiual" + str(len(AbDict))
#print AbDict
'''
