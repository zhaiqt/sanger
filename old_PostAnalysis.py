# example to run this program
#qsub -b y python PostAnalysis.py  -d /dlab/NGS/usem-seqanalysis/160314_zhaiqi1_miseq_HBx52-60DNA.20160214_AN2N4/HC/RESULTS -c H  -s mouse


import argparse
import os
import IsolateClone
import ParseTable 
import ClusterClone
import translator
import AnnotateProtein
import WriteFast
import ReadIgBlastn
parser= argparse.ArgumentParser(prog='cat all.xls files',description="python PostAnalysis.py -d path -s species -c chain",epilog='')
parser.add_argument ('-d','--directory',help='input file directory',default='/home/zhaiqi1/NGS/mycode/Ab_NGS_4/test/results',action='store')
parser.add_argument('-s', '--species', help='mouse, rabbit or human', default="mouse")
parser.add_argument('-c', '--chain', help="folder", default="H")

args=parser.parse_args()

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

Outfile_summary.write("There are  DNA sequences by same CDR3-DNA, GERMLINE-V, RID, DNAlen : %s \n " % str(len(groupDict)))
print ("There are  DNA sequences by same CDR3-DNA, GERMLINE-V, RID, DNAlen : %s \n" % str(len(groupDict)))

#print groupDict
        # example of final groupDict:
        #{('', 'DFL16.1', 'JH1', 'J558.40'): {'M00680:164:000000000-AN2N4:1:2119:22686:25114': 'GGGCCCATGAGGTCCGGCTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATAAACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAAATATTTATCCTTCTGATAGTTATACTAACTACAATCAAAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCCGACATCTGAGGATTCTGCGGTCTATTACTTTATTACTACGGTAGTAGCTACTGCTGGTACTTCGATGTCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCACATTCAAG'}, ('AREGGNYHYFDY', 'DSP2.5', 'JH2', '3:3.9'): {'M00680:164:000000000-AN2N4:1:2119:20823:25147': 'TAAAGTGGGAGGTGCAGCTTCCGGAGTCTGGGGGAGACTTAGTGAAGCCTGGAGGGTCCCTGAAACTCTTCTGTGCAGCCTCTGGATTCACTTTCAGTAGCTATGGCATGTCTTGGGTTCGCCAGACTCCAGACAAGAGGCTGGAGTGGGTCGCAACCATTAGTAGTGGTGGTAGTTACACCTACTATCCAGACAGTGTGAAGGGGCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCAGTCTGAAGTCTGAGGACACAGCCATGTATTACTGTGCAAGAGAGGGGGGTAACTACCACTACTTTGACTACTGGGGCCAAGGCACCACTCTCACCGTCTCCTCAACATTCGT'}}

'''
Outfile_keywords4=os.path.join(args.directory,"Aggregateclone.txt")
Outfile_Aggregate= open(Outfile_keywords4, "w")
for seq,id in mergedDict.iteritems():
	Outfile_Aggregate.write(seq + '\t'+id+'\n')
'''


##################### cluster dna with the same CDR3 and GERMLine-V, count the same protein with barcode
Ab_dict={}
Outfilename_protein=os.path.join(args.directory,"_tmp_protein.txt")
Outfile_tmpprotein=open(Outfilename_protein,'w')
ID=0
Num_seq_inCluster=0
for keyword,tmp_AbDict in groupDict.iteritems():
	cluster_handle = ClusterClone.ClusterClone(tmp_AbDict.values())
	cluster_handle.mergeList()
	Num_seq_inCluster  +=cluster_handle.memberCount_inGroup()
	consensusList = cluster_handle._consensusList
	for consensus in consensusList:
		found_flag = False
		protein_consensus = translator.choose_translation(consensus)
		for abID, info in Ab_dict.iteritems():
			if protein_consensus  ==info['PRO']:
				Ab_dict[abID]['COUNT'] +=1
				found_flag = True
				break
		if found_flag ==False:	
			ID += 1
			Ab_dict[str(ID)]={'DNA':consensus,'COUNT':1, "PRO": protein_consensus,"GERMLINE-V":keyword[1]}
print ("Total unique protein  sequences after error correction and count by unique barcode number\t:%s\n" % str(len(Ab_dict)))
Outfile_summary.write("Total unique protein  sequences after error correction and count by unique barcode number\t:%s\n" % str(len(Ab_dict)))
Outfile_summary.write('%d sequences have are in cluster' % Num_seq_inCluster)
print ('%d sequences have are in cluster' % Num_seq_inCluster)



##### submit 1000 dna.fasta files to Igblasn  ----------##
consensusFastA_filename=WriteFast.writeDict_ProDNA(Ab_dict,args.directory,"Consensus_")
os.system("python /home/zhaiqi1/NGS/mycode/Ab_NGS_3/WrapIgBlastn.py -s %s -i %s" %( args.species,consensusFastA_filename ))
igblastnFilename=consensusFastA_filename.rstrip('.fasta')+".igblastn"

# extract results from Igblastn, the results are returned as dictionary {name: }
foo=ReadIgBlastn.ReadIgBlastn(igblastnFilename)
foo.readIgBlastn()
#print igblastn_results
for key in foo._dict:
	Ab_dict[key].update(foo._dict[key])
#print Ab_dict
##### Anotate Ab Protein 1000 sequences using PWM      #######
foo = AnnotateProtein.AnnotateProtein(Ab_dict,args.species,args.chain)
foo.AnnotateDict()
print "~~~~~~~~~~~~~~~~error corrected Ab_dict"
#print Ab_dict
########################
#keyList=['COUNT',"GERMLINE-V","DNA","PRO",'FR1-PRO','CDR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO','FR1-DNA','CDR1-DNA','FR2-DNA','CDR2-DNA',"FR3-DNA",'CDR3-DNA','FR4-DNA']
print "write all the information of raw_AbDict into Final_corrected.xls"
WriteFast.writeDict_all(Ab_dict,args.directory)

