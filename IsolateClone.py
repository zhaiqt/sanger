import re

def identifyClone (inputDict, keywords=['CDR3-DNA','RID','DNAlen']):
	#keyList=["DNAlen","GERMLINE-V","CDR3-DNA","RID","DNA"]
 	#keyList=["Name","DNAlen","GERMLINE-V","GERMLINE-D","GERMLINE-J","PRODUCT","CHAIN","LID","RID","DNA","PRO",'FR1-PRO','CDR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO','FR1-DNA','CDR1-DNA','FR2-DNA','CDR2-DNA',"FR3-DNA",'CDR3-DNA','FR4-DNA']
	#keywords.sort()
 	keywords = tuple(keywords)
 	groupDict = {}
	# example of final groupDict:
	#{('', 'DFL16.1', 'JH1', 'J558.40'): {'M00680:164:000000000-AN2N4:1:2119:22686:25114': 'GGGCCCATGAGGTCCGGCTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATAAACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAAATATTTATCCTTCTGATAGTTATACTAACTACAATCAAAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCCGACATCTGAGGATTCTGCGGTCTATTACTTTATTACTACGGTAGTAGCTACTGCTGGTACTTCGATGTCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCACATTCAAG'}, ('AREGGNYHYFDY', 'DSP2.5', 'JH2', '3:3.9'): {'M00680:164:000000000-AN2N4:1:2119:20823:25147': 'TAAAGTGGGAGGTGCAGCTTCCGGAGTCTGGGGGAGACTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTTCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGGCATGTCTTGGGTTCGCCAGACTCCAGACAAGAGGCTGGAGTGGGTCGCAACCATTAGTAGTGGTGGTAGTTACACCTACTATCCAGACAGTGTGAAGGGGCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGAGAGGGGGGTAACTACCACTACTTTGACTACTGGGGCCAAGGCACCACTCTCACCGTCTCCTCAACATTCGT'}}

	for abName, abInfo in inputDict.iteritems():
		combinedKey = []
 		for key in keywords:
			combinedKey.append(abInfo[key])
		#if "" in combinedKey or "N/A" in combinedKey:
		#	continue 
		combinedKey = tuple(combinedKey)
		if groupDict.has_key(combinedKey):
			groupDict[combinedKey][abName]= abInfo
		else:
			groupDict[combinedKey]={abName:abInfo}
	#writeCount(groupDict)
	return groupDict

def writeCount(groupDict,outputfilename,keywordlist):
	Outfile=open(outputfilename,'w')
	#print groupDict
	line = '' 
	for keyword in keywordlist:
		line +=keyword +'\t'	
	line += 'Counts \n'
	Outfile.write(line)

	for keys,value in groupDict.iteritems():
		line = ''
		for i,key in enumerate(keys):
			#print key
			line += key + '\t'
		line += str(len(value)) + '\n'
		Outfile.write(line)
	return



