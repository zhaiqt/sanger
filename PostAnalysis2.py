# example to run this program
#qsub -b y python PostAnalysis.py  -d /dlab/NGS/usem-seqanalysis/160314_zhaiqi1_miseq_HBx52-60DNA.20160214_AN2N4/HC/RESULTS -c H  -s mouse

import pdb
import argparse
import os
import IsolateClone
import ParseTable
import ClusterClone
import translator
import AnnotateProtein
import WriteFast
import ReadIgBlastn
import math
parser= argparse.ArgumentParser(prog='cat all.xls files',description="python PostAnalysis.py -d path -s species -c chain",epilog='')
parser.add_argument ('-d','--directory',help='input file directory',default='/home/zhaiqi1/NGS/mycode/Ab_NGS_4/test/results',action='store')
#parser.add_argument ('-d','--directory',help='input file directory',default='/Users/zhaiqi1/Documents/Novartis/my_code/NGS/mycode/Ab_NGS_4/test/results',action='store')
parser.add_argument('-s', '--species', help='mouse, rabbit or human', default="mouse")
parser.add_argument('-c', '--chain', help="folder", default="H")

args=parser.parse_args()

#print args.directory
#print args.chain

############### read the table from the Fastq2fastA################
raw_AbDict,count_seq=ParseTable.ParseTable(args.directory)
print "There are total %s sequences in the table." % str(count_seq)
print "Total number of sequences meets the keywords requirement\t:%s\n" %  str(len(raw_AbDict))

Outfile_summary=open(os.path.join(args.directory,"Summary.txt"),'w')

Outfile_summary.write("Total number of sequences meets the keywords requirement\t:%s\n" %  str(len(raw_AbDict)))
#print raw_AbDict
############################## cluster the clone based on the keywords_3, and then correct the pcr error ########
keywords_3=['CDR3-PRO','RID','DNAlen']
groupDict = IsolateClone.identifyClone(raw_AbDict,keywords_3)
Outfile_keywords3=os.path.join(args.directory,"uniqueclone.txt")
IsolateClone.writeCount(groupDict,Outfile_keywords3,keywords_3)  #this output has not been corrected

Outfile_summary.write("There are  DNA sequences by same CDR3-DNA, RID, DNAlen : %s \n " % str(len(groupDict)))
#print groupDict
print ("There are  DNA sequences by same CDR3-PRO,  RID, DNAlen : %s \n" % str(len(groupDict)))

'''
i=0
for key,info in groupDict.iteritems():
        print key
        print info
        i +=1
        if i>2:
                break

'''
##################### cluster dna with the same CDR3 and GERMLine-V, count the same protein with barcode
Ab_dict={}
Outfilename_protein=os.path.join(args.directory,"_tmp_protein.txt")
Outfile_tmpprotein=open(Outfilename_protein,'w')
Num_seq_inCluster=0
DNAconsensusDict ={}
PROconsensusDict={}
for keyword,tmp_AbDict in groupDict.iteritems():
	if len(tmp_AbDict)<3:
		continue
	cluster_handle = ClusterClone.ClusterClone(tmp_AbDict)
	cluster_handle.mergeList()
	Num_seq_inCluster  += cluster_handle.memberCount_inGroup()
	# combine all the RID and LID which generate the same DNA or protein
	for dna,info in cluster_handle._consensusDict.iteritems():
		protein= translator.choose_translation(dna)
		if dna in DNAconsensusDict.keys():
			DNAconsensusDict[dna]['RID'].extend(info['RID'])
			DNAconsensusDict[dna]['LID'].extend(info['LID'])
		else:
			DNAconsensusDict[dna]=info
		DNAconsensusDict[dna]['PRO']=protein

		if protein in PROconsensusDict.keys():
			PROconsensusDict[protein]['RID'].extend(info['RID'])
			PROconsensusDict[protein]['LID'].extend(info['LID'])
		else:
			PROconsensusDict[protein]=info
		PROconsensusDict[protein]['DNA']=dna

	#example of consensusDict
	#{'GGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAATATTAATCCTAGCAATGGTGGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGATGGGGCAGTAACTACGGGACTTGGTTTGCTTACTGGGGCCAAGGGACTCTGGTCACTGTCTCTGCAATACCATATACCCATACGATGTTCCA': {'LID': ['AATGGAGA', 'ACCGCCGT', 'TGTGGCGA'], 'RID': 'ATACCATA'}}
	#{'GGGGCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGTTGTCCTGTAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAATATTAATCCTAGCAATGGTGGTACTAATTACATTGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGATCCGGATACTACGGTAGTAGCTACAAGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAGTTACTAATACCCATACGATGTTCCA': {'LID': ['TGCTGCGA', 'CAAATCGT', 'CAGATCGT'], 'RID': 'GTTACTAA'}}
'''
i=0
for key,info in PROconsensusDict.iteritems():
	print key
	print info
	i +=1
	if i>2:
		break
'''
######### convert to stabndard library format########

Ab_dict_dna={}
ID = 0
RIDtotal=0
Normtotal=0
for key,info in DNAconsensusDict.iteritems(): # generate Ab_dict based on same dna seq
	RIDcount = len(set(info['RID']))
	LIDcount = len(set(info['LID']))
	Normcount = math.pow(RIDcount,2)/LIDcount 
	Ab_dict_dna[str(ID)]={'RIDcount':RIDcount, 'Normcount':Normcount,'DNA':key, 'PRO':info['PRO']}
	ID +=1
	RIDtotal += RIDcount
	Normtotal += Normcount
for key,info in Ab_dict_dna.iteritems(): 
	Ab_dict_dna[key].update({'RIDcount%':info['RIDcount']/float(RIDtotal), 'Normcount%': float(info['Normcount'])/float(Normtotal)})
print ("Total unique dna  sequences after error correction and count by unique barcode number\t:%s\n" % str(len(Ab_dict_dna)))
Outfile_summary.write("Total unique DNA  sequences after error correction and count by unique barcode number\t:%s\n" % str(len(Ab_dict_dna)))

print("The sum of the RID counts of dna is %d , the sum of theNormalized RID counts is %d. " % (RIDtotal, Normtotal) )
Outfile_summary.write( "The sum of the RID counts of dna is %d , the sum of theNormalized RID counts is %d." % (RIDtotal, Normtotal) )


Ab_dict_pro={}
ID = 0
RIDtotal=0
Normtotal=0
for key,info in PROconsensusDict.iteritems(): # generate Ab_dict based on same dna seq
	RIDcount = len(set(info['RID']))
	LIDcount = len(set(info['LID']))
	Normcount = math.pow(RIDcount,2)/LIDcount
	#pdb.set_trace()
	Ab_dict_pro[str(ID)]={'RIDcount':RIDcount, 'Normcount':Normcount,'DNA':info['DNA'], 'PRO':key}
	ID +=1
	RIDtotal += RIDcount
	Normtotal += Normcount
for key,info in Ab_dict_pro.iteritems():
	Ab_dict_pro[key].update({'RIDcount%':info['RIDcount']/float(RIDtotal), 'Normcount%': float(info['Normcount'])/float(Normtotal)} )

print ("Total unique DNA  sequences after error correction and count by unique barcode number\t:%s\n" % str(len(Ab_dict_dna)))
Outfile_summary.write("Total unique protein  sequences after error correction and count by unique barcode number\t:%s\n" % str(len(Ab_dict_pro)))

print("The sum of the RID counts of dna is %d , the sum of theNormalized RID counts is %d. " % (RIDtotal, Normtotal) )
Outfile_summary.write( "The sum of the RID counts of dna is %d , the sum of theNormalized RID counts is %d." % (RIDtotal, Normtotal) )
Outfile_summary.write('%d sequences have are in cluster' % Num_seq_inCluster)
print ('%d sequences have are in cluster' % Num_seq_inCluster)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~generate output files from Ab_dict_dna~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### submit 1000 dna.fasta files to Igblasn  ----------##
for prefix in ['dna','pro']:
	dictName='Ab_dict_'+prefix
	print '#####'
	print eval(dictName)
	print "There are total %d unique %s sequences." % (len(eval(dictName)), prefix)
	consensusFastA_filename=WriteFast.writeDict_ProDNA(eval(dictName),args.directory,prefix+"_Consensus_")
	os.system("python /home/zhaiqi1/NGS/mycode/Ab_NGS_4/WrapIgBlastn.py -s %s -i %s" %( args.species,consensusFastA_filename ))
	igblastnFilename=consensusFastA_filename.rstrip('.fasta')+".igblastn"

	# extract results from Igblastn, the results are returned as dictionary {name: }
	foo=ReadIgBlastn.ReadIgBlastn(igblastnFilename)
	foo.readIgBlastn()
	#print igblastn_results
	for key in foo._dict:
		eval(dictName)[key].update(foo._dict[key])

	##### Anotate Ab Protein 1000 sequences using PWM      #######
	foo = AnnotateProtein.AnnotateProtein(eval(dictName),args.species,args.chain)
	foo.AnnotateDict()
	print "~~~~~~~~~~~~~~~~error corrected Ab_dict"
	########################
	print "write all the information of raw_AbDict_%s into %s_Final_corrected.xls" %(prefix,prefix)
	Outfile_summary.write("write all the information of raw_AbDict_%s into %s_Final_corrected.xls" %(prefix,prefix))
	WriteFast.writeDict_all(eval(dictName),args.directory,prefix)
