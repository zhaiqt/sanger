
import os


class ReadIgBlastn(object):
    """ return dictionary with germline, CDR boundaries"""
    def __init__(self,inputName):
        self._file = open(inputName,'r')
        self._dict ={}


    def readIgBlastn(self):
            seqID=""
	    germFlag=False
            while True:
                line = self._file.readline()
                if not line:
                    break
		line.strip()
		if line.startswith("# Query: "):
			line=line.split()
			seqID=line[2]
			self._dict[seqID]={'GERMLINE-V':"",'GERMLINE-D':"",'GERMLINE-J':""}
                elif line.startswith("# V-(D)-J rearrangement summary"):
			germFlag=True 
		elif any([line.startswith(x) for x in ("FR1","CDR1","FR2",'CDR2','FR3','CDR3') ]):
			line = line.split()
			y = line[0].split('-')
			self._dict[seqID].update({y[0]+"head_pos":int(line[-7]),y[0]+"tail_pos":int(line[-6])} )
		elif germFlag:
                        if line.startswith("#"):
                                self._dict[seqID]={'GERMLINE-V':"",'GERMLINE-D':"",'GERMLINE-J':""}
                        else:
                                germline=line.split()
                                self._dict[seqID].update({"GERMLINE-V":germline[0],"GERMLINE-D":germline[1],'GERMLINE-J':germline[2]})
                        germFlag=False



'''
test=ReadIgBlastn("./test/results/HC_tmp_DNA.igblastn")
test.readIgBlastn()
print test._dict
for key in test._dict:
    print key
    for key2 in test._dict[key]:
	print test._dict[key]

'''
