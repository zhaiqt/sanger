import zipfile
import os
import TrimEnds
import Merge2FastQ
import ReadSanger
import MergeMultiFastQ
import logging
import csv
import argparse

#python IDAP384screen.py -i /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/IDAP384screen/IDAP_5_star/rawdata/ensingle.csv.zip -c /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/IDAP384screen/IDAP_5_star/rawdata/barcode_Table.csv -o /Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/IDAP384screen/IDAP_5_star/result
parser = argparse.ArgumentParser( prog='SangerPairSeqContig',description="Read paired-ends sanger sequneces, convert to fastq, and merge the Forward and the reverse reads, -> 1) fasta contains the assembled sequences, 2) fasta with non-assembled sequences", epilog='python Sanger2Fastq -i inputfile.zip  -o outputpath')
parser.add_argument ('-i','--inputzip',help='Input zipfile', default='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/raw/TillerLC.zip')
parser.add_argument('-o', '--outputpath',help='outputpath for output',default='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results')
parser.print_help()
args=parser.parse_args()
print args

log_file = args.outputpath+'/runlog.txt'
log_level = logging.DEBUG
logging.basicConfig(filename=log_file, level=log_level, format='%(asctime)s %(message)s')


####unzip the inputfile to the output filepath#####
abi_path = ReadSanger.extract_zip(args.inputzip,args.outputpath)

#generate directories for the All, F , R fastq files
localinfilename = os.path.basename (abi_path)
local_Allname = localinfilename+"_All"
local_Fname = localinfilename+"_F"
local_Rname = localinfilename+"_R"
outfilepath_All = os.path.join(args.outputpath,local_Allname)
outfilepath_F = os.path.join(args.outputpath,local_Fname )
outfilepath_R = os.path.join(args.outputpath,local_Rname )

for dirname in [outfilepath_All,outfilepath_F,outfilepath_R]:
    try:
        os.makedirs(dirname)
    except OSError:
        if os.path.exists(dirname):
            # We are nearly safe
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise


#####  convert Abi/ab1 to fastq files#############


ReadSanger.read_sanger(abi_path,outfilepath_All,outfilepath_F,outfilepath_R )

F_fastq_path = ReadSanger.concentrate2single (outfilepath_F,args.outputpath)
R_fastq_path = ReadSanger.concentrate2single (outfilepath_R,args.outputpath)


######## generate
(All_fastq_path,All_fasta_dict)=MergeMultiFastQ.read_FastQs_dict(outfilepath_All,args.outputpath)

#print len(All_fasta_dict)

All_Fasta_dict = MergeMultiFastQ.Contig_2Fastq_Fasta(All_fasta_dict, args.outputpath)
#print All_Fasta_dict

####### write out the fasta
contig_fasta_name = os.path.join(args.outputpath,'all.fasta')
MergeMultiFastQ.writeFasta (All_Fasta_dict, contig_fasta_name )
