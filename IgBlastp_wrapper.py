


# Purpse: To execute IgBlastp without having to type the whole command each time
#
# Functionality: Replace certain variable parts of the commands with user input and execute IgBlast
# Note: This script has to be withing the IgBlast folder so the folder structure has to follow the IgBlast structure
#
# Usage: python Igblast_wrapper.py -s species -i Input_filename

#from argparse import ArgumentParser
from argparse import ArgumentParser
import os

os.chdir('/home/zhaiqi1/NGS')
print "new diectory:"+ os.getcwd()
# Parsing arguments
parser=ArgumentParser()
parser.add_argument("-s","--species",default='human',help='Enter the species you are interested. Example: human, mouse, rabbit. Default: human')
parser.add_argument("-i","--input",default="protein.txt", help='Input FASTA file')
args=parser.parse_args()

species=args.species
inputname=args.input

print "inputname :"+inputname 
# Executing IgBlast command
#os.system("./igblastp -germline_db_V database/%s_gl_V  -domain_system imgt -num_alignments_V 1 -outfmt 3 -num_threads 2 -query %s -organism %s > %s-igblast.txt" %(species,species,inputname,inputname))
print "$$$$$$$$$$$$$$$$$$$4 running Igblastp"
os.system("/home/zhaiqi1/NGS/IgBlast/1.4.0/ncbi-igblast-1.4.0/bin/igblastp -germline_db_V /home/zhaiqi1/NGS/IgBlast/database/%s_gl_V  -domain_system imgt -num_alignments_V 1 -query %s -organism %s > %s.igblastp" %(species,inputname,species,inputname))
print "IgBlastp has been performed. the results are in %s.igblastp" % inputname

#os.system("igblastp -germline_db_V database/%s_gl_V  -domain_system imgt -num_alignments_V 1 -query %s -organism %s > %s-igblast.txt" %(species,inputname,species,inputname))


#            igblastp -germline_db_V database/human_gl_V -domain_system kabat -num_alignments_V 1 -query protein.txt -outfmt 3 -organism human > output_kabat.txt
