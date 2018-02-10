
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
#Setup variables (could parse command line args instead)
file_f = "/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/F/W1136_HC_Ti.fastq"
file_r = "/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/R/W1136_HC_Ti.fastq"
file_out = "/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/W1136_HC.fastq"
format = "fastq" #or "fastq-illumina", or "fasta", or ...

def interleave(iter1, iter2) :
    for (forward, reverse) in itertools.izip(iter1,iter2):
        print forward
        print reverse
        assert forward.id  == reverse.id
        forward.id += "/1"
        reverse.id += "/2"
        yield forward
        yield reverse

records_f = SeqIO.parse(open(file_f,"rU"), format)
records_r = SeqIO.parse(open(file_r,"rU"), format)

handle = open(file_out, "w")
count = SeqIO.write(interleave(records_f, records_r), handle, format)
handle.close()
print "%i records written to %s" % (count, file_out)
