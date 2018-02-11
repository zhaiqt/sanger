
from Bio import SeqIO
import zipfile
import os
import re

##################
def extract_zip(infile,outfilePath):
    if infile.endswith(".zip"):
        with zipfile.ZipFile(infile, 'r') as z:
            print "zipfile zipfile"
            z.extractall(outfilePath)
        localfilename = os.path.basename (z.filename)
        abi_path = os.path.join(outfilePath,localfilename.rsplit('.zip')[0])
        print abi_path
    return abi_path



###############
def read_sanger(infilepath,outfilepath_All,outfilepath_F,outfilepath_R):
    count_infile = 0
    for filename in os.listdir(infilepath):
        #print "infilepath :" +infilepath
        #print filename
        if filename.endswith('.ab1'):
            abi_filename=os.path.join(infilepath,filename)
            fastq_filename= os.path.join(outfilepath_All,filename.replace('.ab1','.fastq'))
            #print abi_filename
            ##print fastq_filename
            SeqIO.convert(abi_filename,'abi',fastq_filename,'fastq')

            if filename.endswith("QB5505.ab1"):
                F_fastq_filename =os.path.join(outfilepath_F,filename.replace('.ab1','.fastq'))
                SeqIO.convert(abi_filename,'abi',F_fastq_filename,'fastq')
            else:
                R_fastq_filename =os.path.join(outfilepath_R,filename.replace('.ab1','.fastq'))
                SeqIO.convert(abi_filename,'abi',R_fastq_filename,'fastq')

            count_infile += 1
    print "There are total %d sequences" %count_infile

    return

##########



############################################
def concentrate2single(infilepath,outputfilename):
    combinedfile_name = os.path.basename(infilepath)
    print "Xxxxx"
    print combinedfile_name
    combinedfile_name =os.path.join(outputfilename,combinedfile_name+'.fastq')
    combinedfile = open(combinedfile_name , 'wb')
    print "$$$$$$$$$$$$$"
    print outputfilename
    for filename in os.listdir(infilepath):
        if filename.endswith('.fastq'):
            seq_name = filename.split(';')[0]
            with open(os.path.join(infilepath,filename)) as f:
                for row in f:
                    if row.startswith("@"):
                        row = row.split(";")
                        row= row[0]+"\n"
                    combinedfile.write(row)
                f.close()
    combinedfile.close()
    return combinedfile_name
