
import os


class IgBlastp(object):
    """ return dictionary with germline, CDR boundaries"""
    def __init__(self,inputName):
        self._file = open(inputName,'r')
        self._dict ={}


    def extract_IgBlastp(self):
            seqID=""
            while True:
                line = self._file.readline()
                if line=="":
                    break
                line=line.strip()
                if line.startswith("Query="):
                    line=line.split(" ")
                    #self._IgBlastp_boundary["name"] = line[1]
                    seqID=line[1]
                    self._dict[seqID]={}
                elif line.startswith("lcl|"):
                    line=line.split()
                    self._dict[seqID].update({"germline": line[0][4:]})
                    #self._dict[seqID]["germline"]=line[0][4:]

'''
        #extract the CDR position from IgBlast
                elif line.startswith("FR1"):
                    line = line.split()
                    self._dict[seqID]["FR1"]=line[1:7]
                elif line.startswith('CDR1'):
                    line = line.split()
                    self._dict[seqID]["CDR1"] = line[1:7]
                elif line.startswith("FR2"):
                    line = line.split()
                    self._dict[seqID]["FR2"] = line[1:7]
                elif line.startswith("CDR2"):
                    line = line.split()
                    self._dict[seqID]["CDR2"] = line[1:7]
                elif line.startswith("FR3"):
                    line = line.split()
                    self._dict[seqID]["FR3"] = line[1:7]
    '''


"""test=IgBlastp("./results/assembled_protein-igblastp")
test.extract_IgBlastp()
print test._dict
for key in test._dict:
    print key
    for key2 in test._dict[key]:
        print key2
        print test._dict[key][key2]"""



