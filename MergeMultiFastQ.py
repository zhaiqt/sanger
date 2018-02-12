import os
import Merge2FastQ


###### read all csv files in the file
def read_FastQs_dict(infilepath,outputfilename):

    combinedfile_name = os.path.basename(infilepath)
    combinedfile_name =os.path.join(outputfilename,combinedfile_name+'.fastq')
    combinedfile = open(combinedfile_name , 'wb')

    all_fastq = {}
    count_file =0
    for filename in os.listdir(infilepath):
        if filename.endswith('.fastq'):
            seq_name = filename.split(';')[0]
            with open(os.path.join(infilepath,filename)) as f:
                single_fastq = []
                name =""
                count_file  +=1

                for row in f:
                    row = row.rstrip('\n')
                    if row.startswith("@"):
                        row = row.split(";")[0]
                        name = row.lstrip('@')
                    single_fastq.append(row)
                    combinedfile.write(row+'\n')

                single_fastq[-1] = covert_Qscore(single_fastq[-1])
                all_fastq[name] = all_fastq.get(name, [])
                all_fastq[name].append(single_fastq)
                #print single_fastq
            f.close()
    combinedfile.close()
    #print all_fastq
    print "combined %d fastq."   %count_file
    print len(all_fastq)
    return (combinedfile_name,all_fastq)


    #fastadict= convert_fastqdict_fastadict(all_fastq)


##########
def Contig_2Fastq_Fasta(fastqdict,outfilePath):
    print "running Contig_2Fastq_Fasta"
    outfilename1= os.path.join (outfilePath , "count.csv")
    outfile1 = open(outfilename1,'wb')
    outfile1.write ("ID,seq1,seq1_length,seq2,seq2_length, assembed_seq, assembed_length, mismatch in overlap \n")

    outfilename2 = os.path.join (outfilePath ,"assembled.fasta")
    outfile2 = open(outfilename2,'wb')

    outfilename3 = os.path.join (outfilePath ,"non_assembled.fasta")
    outfile3 = open(outfilename3,'wb')
    print outfilename1


    all_fasta_dict={}
    for key,value in fastqdict.iteritems():

        objectFastA=Merge2FastQ.Merge2FastQ(value) #apply Merge2FastQ class
        objectFastA.anneal2fastq()
        single_fasta= objectFastA.output_fasta()
        mismach_count=objectFastA.output_mismatch_count()

        all_fasta_dict[key]= single_fasta

        output ="%s,%s,%d,%s,%d,%s,%d,%d\n" %(key,value[0][1],len(value[0][1]),value[1][1],len(value[1][1]),single_fasta,len(single_fasta),mismach_count)
        outfile1.write(output)

        if len(single_fasta) > 300:
            output= ">" +key+"\n"+single_fasta+"\n"
            outfile2.write(output)
        else:
            output = ">" + key +"_1\n" + value[0][1]+"\n"
            outfile3.write(output)
            output = ">" + key +"_2\n" + value[1][1]+'\n'
            outfile3.write(output)
    outfile1.close()
    outfile2.close()
    outfile3.close()

    return all_fasta_dict


#contig_fasta_name

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




'''
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

    #fastadict= convert_fastqdict_fastadict(all_fastq)

    return fastadict

'''
###########
'''
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
'''

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
