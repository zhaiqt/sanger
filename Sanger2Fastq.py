import zipfile
import os
import TrimEnds
import Merge2FastQ
import logging
import csv

from Bio import SeqIO

outfilePath='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/'

log_file = outfilePath+'runlog.txt'
log_level = logging.DEBUG
logging.basicConfig(filename=log_file, level=log_level, format='%(asctime)s %(message)s')

def extract_zip(infile,outfilePath):
    if infile.endswith(".zip"):
        with zipfile.ZipFile(infile, 'r') as z:
            print "zipfile zipfile"
            z.extractall(outfilePath)

###############
def read_sanger(infilepath,outfilepath):
    print infilepath
    outfilepath0 = outfilepath+'all/'
    outfilepath1=outfilepath+'F/'
    outfilepath2=outfilepath+'R/'

    i = 0
    for filename in os.listdir(infilepath):
        if filename.endswith('.ab1'):
            filename_list= filename.split(";")

            fastq_name =filename_list[0]+".fastq"
            fastq_name=filename.replace('.ab1','.fastq')
            SeqIO.convert(infilepath+filename,'abi',outfilepath0+fastq_name,'fastq')
            #fastq_name_inall = outfilepath0 + fastq_name_inall
            #print fastq_name_inall
            #filename=infilepath+filename
            #SeqIO.convert(infilepath+filename,'abi',fastq_name_inall,'fastq')

            if filename_list[1].endswith("QB5505.ab1"):
                fastq_name =filename.split(";")[0]+ '.fastq'
                SeqIO.convert(infilepath+filename,'abi',outfilepath1+fastq_name,'fastq')
            else:
                fastq_name =filename.split(";")[0]+ '.fastq'
                SeqIO.convert(infilepath+filename,'abi',outfilepath2+fastq_name,'fastq')

            i = i+1
    print "There are total %d sequences" %i

    return
#############

###### read all csv files in the file
def merge_FASTQ_files(infilepath,outfile):
    combinedfile = open(outfile, 'wb')
    all_fastq = {}
    i=0
    for filename in os.listdir(infilepath):
        if filename.endswith('.fastq'):
            with open(infilepath+'/'+filename) as f:
                single_fastq = []
                for row in f:
                    single_fastq.append(row.rstrip('\n'))
                    # if row.startswith('@'):
                    #     combinedfile.write(row.replace(';QB5505', '').replace(';QB5506', ''))
                    combinedfile.write(row)
                # fix name of fastq
                name = single_fastq[0].split(';')[0]
                name.lstrip('@')
                #print name
                single_fastq[0]=name
                # cover fastq score to the numbers
                single_fastq[-1] = covert_Qscore(single_fastq[-1])

                all_fastq[name] = all_fastq.get(name, [])

                all_fastq[name].append(single_fastq)
            f.close()
            i =i+1
    print "combined %d fastq."   %i
    combinedfile.close()

    fastadict= convert_fastqdict_fastadict(all_fastq)

    return fastadict

##########
def convert_fastqdict_fastadict(fastqdict):
    outfilename1= outfilePath + "count.csv"
    outfile1 = open(outfilename1,'wb')
    outfile1.write ("ID,seq1,seq2,seq1_length, seq2_length, assembed_seq, assembed_length\n")
    #outfilename2= outfilePath+"short_fastq.txt"
    #outfile1 = open(outfilename1,'wb')

    all_fasta_dict={}
    for key,value in fastqdict.iteritems():

        objectFastA=Merge2FastQ.Merge2FastQ(value) #apply Merge2FastQ class
        objectFastA.anneal2fastq()
        single_fasta= objectFastA.output_fasta()
        print "^^^^^^^^"
        print single_fasta
        all_fasta_dict[key]= single_fasta

        output ="%s,%s,%d,%s,%d,%s,%d\n" %(key,value[0][1],len(value[0][1]),value[1][1],len(value[1][1]),single_fasta,len(single_fasta))
        outfile1.write(output)

        if len(single_fasta) < 300:
            logging.info(value)

    return all_fasta_dict




###########
def writeFasta(infastadict,outfile):
    #print infastadict
    combinefile = open(outfile, 'wb')
    for key, value in infastadict.iteritems():
        try:
            output= ">" +key+"\n"+value+"\n"
            combinefile.write(output)
        except:
            logging.info(key )
            logging.info(value)



########
def covert_Qscore(String):
	Q_score=[]
	for i in String:
		Q_score.append(int((ord(i))-33))

        #Q_score ([ord(x)-33 for x in String])
	return Q_score

########
def concentrate2single(infilepath,outputfilename):
    combinedfile = open(outputfilename, 'wb')
    print "$$$$$$$$$$$$$"
    print outputfilename
    for filename in os.listdir(infilepath):
        if filename.endswith('.fastq'):
            with open(infilepath+'/'+filename) as f:
                for row in f:
                    combinedfile.write(row)
                f.close()
    combinedfile.close()


###########

from Bio import SeqIO
import itertools
import sys
import os

def merge_fastq(fastq_path1, fastq_path2, outpath):
    outfile = open(outpath,"w")
    fastq_iter1 = SeqIO.parse(open(fastq_path1),"fastq")
    fastq_iter2 = SeqIO.parse(open(fastq_path2),"fastq")
    for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
        rec1.id = rec1.id.split(';')[0]
        rec2.id = rec2.id.split(';')[0]
        rec1.name = rec1.name.split(';')[0]
        rec2.name = rec2.name.split(';')[0]
        rec1.description = rec1.description.split(';')[0]
        rec2.description = rec2.description.split(';')[0]
        SeqIO.write([rec1,rec2], outfile, "fastq")
    outfile.close()


##########################
'''
#Setup variables (could parse command line args instead)
file_f = "/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/F/W1136_HC_Ti.fastq"
file_r = "/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/R/W1136_HC_Ti.fastq"
file_out = "/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/W1136_HC.fastq"
format = "fastq" #or "fastq-illumina", or "fasta", or ...

def interleave(iter1, iter2) :
    for (forward, reverse) in itertools.izip(iter1,iter2):
        #assert forward.id  == reverse.id
        #forward.id += "/1"
        #reverse.id += "/2"
        yield forward
        yield reverse

records_f = SeqIO.parse(open(file_f,"rU"), format)
records_r = SeqIO.parse(open(file_r,"rU"), format)

handle = open(file_out, "w")
count = SeqIO.write(interleave(records_f, records_r), handle, format)
handle.close()
print "%i records written to %s" % (count, file_out)
'''
##########



#F_abi_file='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/raw/F.zip'
#R_abi_file='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/raw/R.zip'
all_abi_file='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/raw/Ab1_all/'
#all_abi_file='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/raw/all/'
outfilePath='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/'


# log_file = outfilePath+'runlog.txt'
# log_level = logging.DEBUG
# logging.basicConfig(filename=log_file, level=log_level, format='%(asctime)s %(message)s')


read_sanger(all_abi_file, outfilePath)
concentrate2single('/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/F/',outfilePath+'F.fastq')
concentrate2single('/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/R/',outfilePath+'R.fastq')
#extract_zip(F_abi_file,outfilePath+'/F/')
#extract_zip(R_abi_file,outfilePath+'/R/')

#merge_FASTQ_files(outfilePath+'QB5505',outfilePath+'F.fastq' )
#merge_FASTQ_files(outfilePath+'QB5506',outfilePath+'R.fastq' )

fastadict=merge_FASTQ_files(outfilePath+'all',outfilePath+'all.fastq' )

writeFasta(fastadict,outfilePath+'all.fasta')
print outfilePath+'all.fasta'

# merge_fastq(outfilePath+'F.fastq',outfilePath+'R.fastq',outfilePath+'merge.fastq')
# read_sanger('/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/raw/F.zip','/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/')
