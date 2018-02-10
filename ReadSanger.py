
from Bio import SeqIO
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
