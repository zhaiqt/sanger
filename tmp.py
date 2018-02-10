import argparse
import os
import Parse2FastQ 
import ExtractBarcodeQ
import TrimEnds
import Merge2FastQ
import translator
import WriteFast
import ReadIgBlastn
import AnnotateProtein
############ set up directory for input and output#####

parser = argparse.ArgumentParser(prog='process and assemble Ab FastQ, translate and annotate CDR and germline',
                                 description="python workflow.py -s human -c h", epilog='')
parser.add_argument('-s', '--species', help='mouse, rabbit or human', default="mouse")
parser.add_argument('-c', '--chain', help="folder", default="H")
parser.add_argument('-k', '--keylist',nargs='+', help="DNA,FV,CDR1,CDR2,CDR3, GERMLINE", default="DNAlen GERMLINE-V CDR3-DNA RID")
parser.add_argument ('-d',"--directory",type=str, help="input files directory", default='/home/zhaiqi1/NGS/mycode/Ab_NGS_3/test',action='store')
parser.add_argument ('-o',"--outputpath",type=str, help="input files directory", default='/home/zhaiqi1/NGS/mycode/Ab_NGS_3/test/results/')

args = parser.parse_args()

print ("input file: %s" % args.species)
print ("string size: %s" % args.chain)
print ("input files directory")
print  args.directory
if isinstance(args.keylist, basestring): keylist=[args.keylist]
else:
	keylist=args.keylist
#try:
#	keylist.remove("DNA")
#except:
#	pass
#keylist.append("DNA")

print ("count abundance of the same keywords: %s" %keylist)
input_path= args.directory

## Simultaneously open R1 and R2 fastq file####
FileName1=""
FileName2=""
for filename in os.listdir(input_path):

	print "filename:::" +filename
	if  filename.endswith("_R1_001.fastq"):
		R1fastqFile=os.path.join(input_path,filename)
		FileName1=filename.rstrip("_R1_001.fastq") 
	elif filename.endswith("_R2_001.fastq"):
		R2fastqFile=os.path.join(input_path,filename)
		FileName2=filename.rstrip("_R2_001.fastq")
if FileName1 == FileName2: 
	prefix=FileName1.split("_")[0]
 
print "file1 is : %s" % R1fastqFile 
print "file2 is : %s" % R2fastqFile

#import shutil
#try:
#	shutil.rmtree(args.outputpath) #remove the RESUTLS folder 
#except:
	#os.makedirs(args.outputpath) #recreate the RESULTS folder
#	print
try:
	os.makedirs(args.outputpath) 
except:
	pass
'''
Outfilename1 = "trimed_fastq_pair.txt"
args.outputpath_trimed2FastQ = os.path.join(args.outputpath, Outfilename1)
open(args.outputpath_trimed2FastQ,'w').close()  # clear the content of this file 
Outfilename2 = "mergedFastA.txt"
args.outputpath_FastA = os.path.join(args.outputpath, Outfilename2)
open(args.outputpath_FastA,'w').close()
'''

############ read one fastq from each file, Trim the ends, and Assemble before Reading the next fastq#####

print "#########################################################################################"
print "######### READing FASTQ###################################################################"
c=Parse2FastQ.Parse2FastQ(R1fastqFile,R2fastqFile)
#print "c.next"
#print type(c) 
#tmp = c.next()
#print tmp

end_of_file_flag=False
readCounts=0
IDCounts=0
emptyCounts =0
while not end_of_file_flag:
	pairDict=[]
	# read 1000 fastq from read 1 and read2##
	for i in range(1000):
		try:
			pairDict.append(c.next())
			readCounts +=1
		except:
			end_of_file_flag = True
			break

	#Trim from both end and Merge###	
	trimPairDictQ=[]	
	AbDict={}
	
	for pair in pairDict:
		# obtain barcode from the FastQ pair	
		foo=ExtractBarcodeQ.ExtractBarcode(pair,('',''))
		LRID=foo.extractBarcode()

		# remove low quality region 
		trimedQ2=[]
		for element in pair:	
			#Trim from 3 end###
			fastq=TrimEnds.TrimEnds(element)
			trimedQ2.append(fastq.trim3End())
       		trimPairDictQ.append(trimedQ2)
		#WriteFast.write2FastQ(Q2,args.outputpath_trimed2FastQ )  # write and accumulate the fastq file

		
		#Merger 2 fastq files and return fasta
		objectFastA=Merge2FastQ.Merge2FastQ(trimedQ2)
		fasta=objectFastA.merge2FastQ()
		fasta[1]=fasta[1][8:-8]

		if not fasta[1]:
			emptyCounts +=1
		elif len(fasta[1]) < 300:
			continue
                #translate after barcode

		protein=translator.choose_translation(fasta[1])
		
		# counts if both end barcodes suvived
		if len(LRID[0]) >4 and len(LRID[1])>4 :
			IDCounts +=1
			
		# Generate AbDict library: Key 
		AbDict[fasta[0]]={}
		AbDict[fasta[0]].update({"DNA":fasta[1]})	
		AbDict[fasta[0]].update({"DNAlen":len(fasta[1])})
		AbDict[fasta[0]].update({"LID":LRID[0]})
		AbDict[fasta[0]].update({"RID":LRID[1]})		
		AbDict[fasta[0]].update({"PRO":protein})


	#WriteFast.write2FastQ(Q2,args.outputpath_trimed2FastQ )  # write and accumulate the fastq file
	WriteFast.writeDictFastQ(trimPairDictQ,args.outputpath,prefix)			
	print "The trimed and processed FastQ files are in \n %s/trimed_fastq_pair.txt" % args.outputpath
	
	#This command write 3 files, and return the 1000 dna fasta file name##
	#3 files including: merged.fasta, _tmp_DNA.fasta, _tmp_protein.fasta	
        print "~~~~~~~~~~~~~~~~~~~Writing Merged Fasta files~~~~~~~~~~~~~"
        print "all %d fasta are merged and stored in \n %s/merged.fasta. merged.fasta contains all, including those with no DNA sequence." % (readCounts, args.outputpath)
        print "The current batch of merged non-empty PROTEIN fasta are in \n %s/_tmp_protein.fasta" % args.outputpath
        print "The current batch of merged non-empty DNA fasta are in \n %s/_tmp_DNA.fasta" % args.outputpath
	tmpDNA_filename=WriteFast.writeDict_ProDNA(AbDict,args.outputpath,prefix )
	
	# submit 1000 dna.fasta files to Igblasn  ----------##
	os.system("python /home/zhaiqi1/NGS/mycode/Ab_NGS_3/WrapIgBlastn.py -s %s -i %s" %( args.species, tmpDNA_filename))
	igblastnFilename=tmpDNA_filename.rstrip('.fasta')+".igblastn" 
	
	
	# extract results from Igblastn, the results are returned as dictionary {name: }	
	#'M00680:164:000000000-AN2N4:1:2119                                                     :11126:25130': {'GERMLINE-J': 'VK', 'FR1head_pos': 1, 'CDR3head_pos': 277, 'FR2tail_pos': 159, 'CDR2tail_pos': 168, 'CDR1head_pos': 7                                                     9, 'GERMLINE-D': 'JK1', 'CDR1tail_pos': 108, 'FR2head_pos': 109, 'CDR2head_pos': 160, 'FR1tail_pos': 78, 'GERMLINE-V': '21-4', 'CDR3t                                                     ail_pos': 296},
	foo=ReadIgBlastn.ReadIgBlastn(igblastnFilename)	
	foo.readIgBlastn()
	#print igblastn_results	
	for key in foo._dict: 	
		#print key
		#print igblastn_results[key]
		AbDict[key].update(foo._dict[key])
	#print AbDict

        ##### Anotate Ab Protein 1000 sequences using PWM      ####### 
	foo = AnnotateProtein.AnnotateProtein(AbDict,args.species,args.chain)
	foo.AnnotateDict()
	print AbDict		

	######## write all the information of AbDict into all.txt######
	### keyList=["GERMLINE-V","GERMLINE-D","GERMLINE-J","PRODUCT","CHAIN","LID","RID","DNA","PRO",'FR1-PRO','CDR1-PRO','FR2-PRO','CDR2-PRO',"FR3-PRO",'CDR3-PRO','FR4-PRO','FR1-DNA','CDR1-DNA','FR2-DNA','CDR2-DNA',"FR3-DNA",'CDR3-DNA','FR4-DNA']
	print "write all the information of AbDict into all.txt"
	#WriteFast.writeDict_keys(AbDict,args.outputpath,prefix,keylist)
	WriteFast.writeDict_keys(AbDict,args.outputpath,prefix)

'''
	
for key in protein_fasta_dict:
    #print protein_fasta_dict[key] + "\n"
    annotated_Seq = annotate_Ab(protein_fasta_dict[key], args.species.lower(), args.chain.upper())
    if annotated_Seq:
        antibody_dict[key]=annotated_Seq
        antibody_dict[key].update({"PROTEIN":protein_fasta_dict[key]})
        antibody_dict[key].update({"DNA":dna_fasta_dict[key]}) 	
'''	

print "There are %d FastQ Seqs in each file." % readCounts
print "There are %d FastQ seqs have 5' ID barcode" % IDCounts 
print "Empty Sequences after assembly are %d." % emptyCounts 
