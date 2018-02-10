#import pdb
import os
import NameGermline
def write2FastQ(Q2,outpath):
	Outfile_trimed2FastQ = open(outpath, "a+")
	for fastq in Q2:
		for element in fastq:
			Outfile_trimed2FastQ.write("%s\n" % element)
	Outfile_trimed2FastQ.close()
	return 

def writeFastA(FastA,outpath):
	Outfile_FastA = open(outpath, "a+")
	Outfile_tmpFastA = open(outpath+"_tmp", "w")
	entry= ">"+FastA[0]+"\n"+FastA[1]+"\n"
	Outfile_FastA.write (entry) 
	Outfile_tmpFastA.write (entry)
	Outfile_FastA.close()
	Outfile_tmpFastA.close()
	return

def writeDictFastA(inDict,outpath):
        Outfile_FastA = open(outpath+"/merged.fasta", "a+")
        Outfile_tmpFastA = open(outpath+"/_tmp_merged.fasta", "w")
	for fasta in iter(inDict.keys()):
                entry= ">"+fasta[0]+"\n"+fasta[1]+"\n"
               	Outfile_FastA.write (entry)
               	Outfile_tmpFastA.write (entry)
        Outfile_FastA.close()
        Outfile_tmpFastA.close()
	return

def writeDictFastQ(inDict,outpath,prefix):
        Outfile_trimed2FastQ = open(os.path.join(outpath,prefix+"_trimed_fastq_pair.txt"), "a+")
        for pairQ in inDict:
		for fastq in pairQ:
                	for element in fastq:
                        	if element == '+':
					continue
				Outfile_trimed2FastQ.write("%s\n" % element)
        Outfile_trimed2FastQ.close()
        return


def writeDict_ProDNA(inDict,outpath, prefix):
	Outfile_tmpPro = open(os.path.join(outpath,prefix+"_tmp_protein.fasta"), "w")
	Outfile_tmpDNA = open(os.path.join(outpath,prefix+ "_tmp_DNA.fasta"), "w")
	Outfile_FastA = open(os.path.join(outpath, prefix+"_merged.fasta"), "a+")
	for seqName in iter(inDict.keys()):	
		dnaFasta=">"+seqName+'\n'+inDict[seqName]["DNA"]+"\n"
		proFasta=">"+seqName+'\n'+inDict[seqName]["PRO"]+"\n"
		if inDict[seqName]["DNA"] :
			Outfile_tmpDNA.write(dnaFasta)
		if inDict[seqName]["PRO"] :
			Outfile_tmpPro.write(proFasta)
		Outfile_FastA.write(dnaFasta)
	return os.path.join(outpath,prefix+ "_tmp_DNA.fasta")

def writeDict_keys(inDict,outpath,prefix):
	Outfile_all = open(os.path.join(outpath,prefix+"_all.xls"), "a+")
	Outfile_error= open(os.path.join(outpath,prefix+"_all_error.xls"), "a+")
        keyList=['RID','LID','DNAlen','FV-PRO','FV-DNA',"GERMLINE-V",'GERMLINE-D','GERMLINE-J','FR1-PRO','CDR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO','PRO',"DNA",'FR1-PRO','CDR1-DNA','FR2-PRO','CDR2-DNA','FR3-DNA','CDR3-DNA','FR4-DNA','FR4-DNA']

	Outfile_all.write("# Name\t")

	for keyword in keyList:
		Outfile_all.write(keyword+'\t') 
	Outfile_all.write("\n")

	for ID,info in inDict.iteritems():
		line=ID+"\t"
		#print info.values()
		if '' in info.values():
			print "\t space was found in" + ID 
			continue
		for keyword in keyList:
			#pdb.set_trace()
			#print keyword + ":" + str(info[keyword])
			try:
				#line +=str(info[keyword])
				#line +="\t"
				#print "found it !"
				tmp= info[keyword]
				#if tmp=="" :
				#	Outfile_error.write(line+info["DNA"]+"\n")
				#	line=""
				#	break
				#else:
				#	line += (str(tmp)+"\t")
				line += (str(tmp)+"\t")
				
			except:
				print "IN EXCEPT, key is  "+ keyword
				print "info is :" + str(info[keyword])
				Outfile_error.write(line+str(info['DNAlen'])+info["DNA"]+"\n")
				line=""
				break
	
		if line != "" :
			Outfile_all.write(line+"\n")

	return
			

def writeDict_keys2(inDict,prefix,outpath):

        Outfile_all = open(os.path.join(outpath,prefix+"_individual.xls"), "a+")
        Outfile_error = open(os.path.join(outpath,prefix+"error_individual.xls"),"a+")

        keyList=['RID','LID','DNAlen','FV-PRO','FV-DNA',"GERMLINE-V",'GERMLINE-D','GERMLINE-J','FR1-PRO','CDR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO','PRO',"DNA",'FR1-PRO','CDR1-DNA','FR2-PRO','CDR2-DNA','FR3-DNA','CDR3-DNA','FR4-DNA','FR4-DNA']
        # keyList=["COUNT","GERMLINE-V","DNA","PRO",'FR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO']

        ##Write the keyword title

        Outfile_all.write("# Name\t")

        for keyword in keyList:
                Outfile_all.write(keyword+'\t')
        Outfile_all.write("\n")

        for ID,info in inDict.iteritems():
                line='ID_'+ID+"\t"
                for keyword in keyList:
                        try:
                                tmp= info[keyword]
                                line += (str(tmp)+"\t")
                        except:
                                line=""
                                break
                if line != "" :
                        Outfile_all.write(line+"\n")
        return


def writeDict_all(inDict,outpath,prefix):
        Outfile_all = open(os.path.join(outpath,prefix+"_Final_corrected.xls"), "w")
	Outfile_error = open(os.path.join(outpath,prefix+"_error_aggregate.xls"),"w")

	keyList=['CDR3-PRO','RIDcount','RIDcount%','Normcount','Normcount%','FV-PRO','FV-DNA',"GERMLINE-V",'GERMLINE-D','GERMLINE-J','FR1-PRO','CDR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO','PRO',"DNA",'FR1-DNA','CDR1-DNA','FR2-DNA','CDR2-DNA','FR3-DNA','CDR3-DNA','FR4-DNA','FR4-DNA'] 
        ##Write the keyword title

        Outfile_all.write("# Name\t")

        for keyword in keyList:
                Outfile_all.write(keyword+'\t')
        Outfile_all.write("\n")

        for ID,info in inDict.iteritems():
                line='ID_'+ID+"\t"
                for keyword in keyList:
                        #pdb.set_trace()
			#print "$$$$ keyword"
			#print keyword
			#print "info: "
			#print info
                        try:
                                #line +=str(info[keyword])
                                #line +="\t"
                                #print "found it !"
                                tmp= info[keyword]
				if keyword=='GERMLINE-V':
					tmp =NameGermline.IMGT_germline(tmp)
                                line += (str(tmp)+"\t")

                        except:
                                #print "IN EXCEPT, key is  "+ keyword
                                #print "info is :" + info[keyword]
                                #Outfile_error.write(line+str(info['DNAlen'])+'\t'+info["DNA"]+"\n")
                                line=""
                                break
		#print line
                if line != "" :
                        Outfile_all.write(line+"\n")

        return

	
#################
'''
a={'dadfdgfda':{'LID': 'ACTATAAA', 'DNAlen': 360, 'FR1tail_pos': 75, 'DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAATATTAATCCTAGCAATGGTGGTACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTATATAATTGTGCAAGAAGGGGCAGTAACTACGGGACTTGGTTTGCATACTGGGGCCAAGGGACTCTGGTCACTGTCTCTGCA', 'CDR3head_pos': 289, 'CDR2tail_pos': 174, 'FR1-PRO': 'QVQLQQPGTELVKPGASVKLSCKAS', 'CDR2-PRO': 'INPSNGGT', 'CDR3-DNA': 'GCAAGAAGGGGCAGTAACTACGGGACTTGGTTTGCATAC', 'CDR1tail_pos': 99, 'FR3-DNA': 'ACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTATATAATTGT', 'FR4-PRO': 'WGQGTLVTVS', 'FR1-DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCT', 'FR2-PRO': 'MHWVKQRPGQGLEWIGN', 'FR3tail_pos': 288, 'CDR1-DNA': 'GGCTACACCTTCACCAGCTACTGG', 'CDR1-PRO': 'GYTFTSYW', 'GERMLINE-J': 'JH3', 'FR1head_pos': 1, 'CDR3tail_pos': 294, 'FR4-DNA': 'TGGGGCCAAGGGACTCTGGTCACTGTCTCT', 'FR2tail_pos': 150, 'CDR2-DNA': 'ATTAATCCTAGCAATGGTGGTACT', 'PRO': 'QVQLQQPGTELVKPGASVKLSCKASGYTFTSYWMHWVKQRPGQGLEWIGNINPSNGGTNYNEKFKGKATLTVDKSSSTAYMQLSSLTSEDSAVYNCARRGSNYGTWFAYWGQGTLVTVSA', 'CDR1head_pos': 76, 'GERMLINE-D': 'DSP2.x', 'FR3head_pos': 175, 'FR2head_pos': 100, 'CDR2head_pos': 151, 'CDR3-PRO': 'ARRGSNYGTWFAY', 'FV-DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAATATTAATCCTAGCAATGGTGGTACTACTAACTACAATGAGAAGTTCAAGGGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTATATAATTGTGCAAGAAGGGGCAGTAACTACGGGACTTGGTTTGCATACTGGGGCCAAGGGACTCTGGTCACTGTCTCT', 'FV-PRO': 'QVQLQQPGTELVKPGASVKLSCKASGYTFTSYWMHWVKQRPGQGLEWIGNINPSNGGTTNYNEKFKGKATLTVDKSSSTAYMQLSSLTSEDSAVYNCARRGSNYGTWFAYWGQGTLVTVS', 'RID': 'ATTATACT', 'FR3-PRO': 'TNYNEKFKGKATLTVDKSSSTAYMQLSSLTSEDSAVYNC', 'GERMLINE-V': 'J558.53.146', 'FR2-DNA': 'ATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAAT'},'M00680:178:000000000-AM^YVA:1:1101:18669:2957':{'LID': 'ATATGTAT', 'DNAlen': 295, 'FR1tail_pos': 75, 'DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGTAATATTAATCCTAGCAATGGTGGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGAAATATGAGGACAATGCGGTGTAGTAATGAGAAAAGAGGGATCTTAATACTACGGTAGAAGATACGTACTACTTTGAAAACTGGGGCCAAGGCACCACACTCACAGTCTCCTCA', 'CDR3head_pos': 289, 'CDR2tail_pos': 174, 'FR1-PRO': 'QVQLQQPGTELVKPGASVKLSCKAS', 'CDR2-PRO': 'INPSNGGT', 'CDR3-DNA': 'GAAA', 'CDR1tail_pos': 99, 'FR3-DNA': 'ACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGAAATATGAGGACAATGCGGTGTAGTAATGA', 'FR4-PRO': '', 'FR1-DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCT', 'FR2-PRO': 'MHWVKQRPGQGLEWIGN', 'FR3tail_pos': 288, 'CDR1-DNA': 'GGCTACACCTTCACCAGCTACTGG', 'CDR1-PRO': 'GYTFTSYW', 'GERMLINE-J': 'JH2', 'FR1head_pos': 1, 'CDR3tail_pos': 292, 'FR4-DNA': '', 'FR2tail_pos': 150, 'CDR2-DNA': 'ATTAATCCTAGCAATGGTGGTACT', 'PRO': '*GDCECGALAPVFKVVRIFYRSIKIPLFSLLHRIVLIFQAAELHVGCAGGFVYSQCGLALELLIVVSTTIARINITNPLKALSRPLLHPVHPVAGEGVARSLAGQLH*SPRLHQFSPRLLQLDL', 'CDR1head_pos': 76, 'GERMLINE-D': 'DFL16.1', 'FR3head_pos': 175, 'FR2head_pos': 100, 'CDR2head_pos': 151, 'CDR3-PRO': 'E', 'FV-DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGTAATATTAATCCTAGCAATGGTGGTACTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGAAATATGAGGACAATGCGGTGTAGTAATGAGAAA', 'FV-PRO': 'QVQLQQPGTELVKPGASVKLSCKASGYTFTSYWMHWVKQRPGQGLEWIGNINPSNGGTLALELLIVVSTTIARINITNPLKALSRPLLHPVHPVAGEE', 'RID': 'TGGAAAAG', 'FR3-PRO': 'LALELLIVVSTTIARINITNPLKALSRPLLHPVHPVAGE', 'GERMLINE-V': 'J558.53.146', 'FR2-DNA': 'ATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGTAAT'}}

writeDict_keys(a,'test2/','writeonly')
'''
