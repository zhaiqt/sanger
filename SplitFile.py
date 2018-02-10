import os
import argparse

def splitFile(input_path,lineChunk):
	os.chdir(input_path)
	chucksize=lineChunk
	fout=None
	for filename in os.listdir(input_path):
		if not filename.endswith(".fastq"):
			continue
		with open(filename,"rt") as file:
			for (i, line) in enumerate(file):
				if i % chucksize ==0:
					if fout: fout.close()
					output_path=input_path+"/part"+ str(i/chucksize)
					try:
						os.mkdir(output_path)
					except:
						pass
					OutFileName=os.path.join(output_path,"part%d_%s" %(i/chucksize,filename))
					fout= open(OutFileName,'wb')
				fout.write (line)
			fout.close()
	return

#############
parser = argparse.ArgumentParser(prog='split large files into small files', description="python SplitFile.py -d inputpath -s size", epilog='')
#parser.add_argument("-d","--directory",type=str, help="large fastq files directory", default="/dlab/NGS/usem-seqanalysis/160314_zhaiqi1_miseq_HBx52-60DNA.20160214_AN2N4/HC", action="store")
parser.add_argument("-s",'--size',help="the line number of the desired small file, should be 4N", type=int, default=400)
parser.add_argument("-d","--directory",type=str, help="large fastq files directory", default="/home/zhaiqi1/NGS/mycode/Ab_NGS_2_protein/test", action="store")
args=parser.parse_args()


splitFile(args.directory, args.size)
