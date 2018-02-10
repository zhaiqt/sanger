import json
import pdb
import os
import translator
#parser = argparse.ArgumentParser( prog='calculate PMW, and find the most match string',description="", epilog='python get_PMW_table.py -i inputfile')
#parser.add_argument ('-i','--input',help='Input File Name', default="./data/mFR1.txt")

class AnnotateProtein(object):

	def __init__(self,inDict,species,chaintype):

        	self._dict=inDict

        	if species.lower() in ["mouse","human","rabbit"]:
			self._species=species.lower()
		else:
			print "Error: the species can't be recognized by AnnotateProtein"
	
		if chaintype.upper() in ["H",'HC',"VH"]:
                        self._chain='H'
		elif chaintype.upper() in ["K","KAPPA","KL","LK","VK","KV","L","VL","LV","LC"]:
                        self._chain='K'
                elif chaintype.upper() in ["LAMBDA"]:
                        self._chain='L'

		databasePath='/home/zhaiqi1/NGS/PWM'

		PMW_FR1head_name=databasePath+'/'+self._species+'/'+self._chain+'/FR1_head_PMW_'+self._species+'.json'
		PMW_FR1tail_name=databasePath+'/'+self._species+'/'+self._chain+'/FR1_tail_PMW_'+self._species+'.json'
		PMW_FR2head_name=databasePath+'/'+self._species+'/'+self._chain+'/FR2_head_PMW_'+self._species+'.json'
		PMW_FR2tail_name=databasePath+'/'+self._species+'/'+self._chain+'/FR2_tail_PMW_'+self._species+'.json'
		PMW_FR3head_name=databasePath+'/'+self._species+'/'+self._chain+'/FR3_head_PMW_'+self._species+'.json'
		PMW_FR3tail_name=databasePath+'/'+self._species+'/'+self._chain+'/FR3_tail_PMW_'+self._species+'.json'
		PMW_FR4head_name=databasePath+'/'+self._species+'/'+self._chain+'/FR4_head_PMW_'+self._species+'.json'
		#PMW_FR4tail_name=databasePath+'/'+self._species+'/'+self._chain+'/FR4_tail_PMW_'+self._species+'.json'
	
		with open(PMW_FR1head_name,'r') as fp:
			self._pwmFR1head=json.load(fp)
		with open(PMW_FR1tail_name,'r') as fp:
			self._pwmFR1tail=json.load(fp)
		with open(PMW_FR2head_name,'r') as fp:
			self._pwmFR2head=json.load(fp)
		with open(PMW_FR2tail_name,'r') as fp:
			self._pwmFR2tail=json.load(fp)
		with open(PMW_FR3head_name,'r') as fp:
			self._pwmFR3head=json.load(fp)
		with open(PMW_FR3tail_name,'r') as fp:
			self._pwmFR3tail=json.load(fp)
		with open(PMW_FR4head_name,'r') as fp:
			self._pwmFR4head=json.load(fp)
		#with open(PMW_FR4tail_name,'r') as fp:
			#self._pwmFR4tail=json.load(fp)

	def AnnotateDict(self):	
		#print self._dict
		for seqID,seqInfo in self._dict.iteritems():
			#print "ANnotate"
			#print seqInfo
			self.AnnotateSingleAb(seqID,seqInfo["PRO"],seqInfo["DNA"]) 	
			#print "seqID : "+seqID
			#print "PRO :" + seqInfo['PRO']
			#print "DNA :" + seqInfo['DNA']
		return
	
	def AnnotateSingleAb(self,inSeqID, inProteinSeq,inDNAseq):
                if not inProteinSeq and  not inDNAseq:
        	        return ""

               #initiate the FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4 to emtry, just in case any can't be found
                self._dict[inSeqID].update({'FR1-PRO':"",'CDR1-PRO':'','FR2-PRO':'','CDR2-PRO':'',"FR3-PRO":'','CDR3-PRO':'','FR4-PRO':'','FV-PRO':''})
                self._dict[inSeqID].update({'FR1-DNA':"",'CDR1-DNA':'','FR2-DNA':'','CDR2-DNA':'',"FR3-DNA":'','CDR3-DNA':'','FR4-DNA':'','FV-DNA':''})

		count_missing_fragment =0
		for x in ("FR1","CDR1","FR2",'CDR2','FR3','CDR3'):
			try: 
				#pdb.set_trace()
				dna_frag= inDNAseq[self._dict[inSeqID][x+"head_pos"]-1:self._dict[inSeqID][x+'tail_pos']]
				#print "DNA_frag : "  + dna_frag
				#print "PRO_Frag : " + translator.translate_dna_single(dna_frag)
				self._dict[inSeqID].update({x+'-DNA': dna_frag})
				self._dict[inSeqID].update({x+'-PRO':translator.translate_dna_single(dna_frag)})	
			except: 
				print "IgBlastn couldn't find: "+ inSeqID + ':' + x
				count_missing_fragment += 1
				if count_missing_fragment > 2:
					print "Warining! missing the some fragmentation from IgBlastn in this sequence:"
					print inDNAseq
					return 

							
		FR1head=self.find_high_PMW_score(self._pwmFR1head,inProteinSeq)
		FR1head_pos=inProteinSeq.find(FR1head)

		#print "start to find FR2"
		try :
			FR1tail_pos=self._dict[inSeqID]["FR1tail_pos"]/3 -1
		except: 
			FR1tail=self.find_high_PMW_score(self._pwmFR1tail,inProteinSeq[FR1head_pos:])
			FR1tail_pos=inProteinSeq.find(FR1tail)+len(FR1tail)-1

		new_FR1protein=inProteinSeq[FR1head_pos:FR1tail_pos+1]
		new_FR1dna= translator.convertPROtoDNA(inDNAseq,FR1head_pos,FR1tail_pos+1)
		if len(new_FR1dna)>=len(self._dict[inSeqID]['FR1-DNA']):
			self._dict[inSeqID].update({'FR1-PRO':new_FR1protein})	
			self._dict[inSeqID].update({'FR1-DNA':new_FR1dna})
		'''
		try:
			FR2head=self.find_high_PMW_score(self._pwmFR2head,inProteinSeq[FR1tail_pos:])
			FR2head_pos=inProteinSeq.find(FR2head)
			self._dict[inSeqID].update({'CDR1-PRO':inProteinSeq[FR1tail_pos+1:FR2head_pos]})
			self._dict[inSeqID].update({'CDR1-DNA':translator.convertPROtoDNA(inDNAseq,FR1tail_pos+1,FR2head_pos)})
		except:
			return
		'''

		try :
                        FR2head_pos=(self._dict[inSeqID]["FR2head_pos"]-1)/3
		except:
			FR2head=self.find_high_PMW_score(self._pwmFR2head,inProteinSeq[FR1tail_pos:])
			FR2head_pos=inProteinSeq.find(FR2head)
		if self._dict[inSeqID]["FR2-DNA"] !='':
			pass	
                elif FR1tail_pos and FR2head_pos :
			self._dict[inSeqID].update({'CDR1-PRO':inProteinSeq[FR1tail_pos+1:FR2head_pos]})
                	self._dict[inSeqID].update({'CDR1-DNA':translator.convertPROtoDNA(inDNAseq,FR1tail_pos+1,FR2head_pos)})

		try : 
			FR2tail_pos= (self._dict[inSeqID]["CDR2head_pos"]-1)/3 -1
		except:
			FR2tail=self.find_high_PMW_score(self._pwmFR2tail,inProteinSeq[FR2head_pos:])
			FR2tail_pos=inProteinSeq.find(FR2tail)+len(FR2tail)-1
		if self._dict[inSeqID]["FR2-DNA"] != '' :
                        pass
		elif FR2head_pos and FR2tail_pos:
			self._dict[inSeqID].update({'FR2-PRO':inProteinSeq[FR2head_pos:FR2tail_pos+1]})
			self._dict[inSeqID].update({'FR2-DNA':translator.convertPROtoDNA(inDNAseq,FR2head_pos,FR2tail_pos+1)})
		#print "FR2 tail:" + FR2tail
		#print FR2tail_pos

		try :
			FR3head_pos= (self._dict[inSeqID]["FR3head_pos"]-1)/3 -1
		except:
			FR3head=self.find_high_PMW_score(self._pwmFR3head,inProteinSeq[FR2tail_pos:])
			FR3head_pos=inProteinSeq.find(FR3head)
		if self._dict[inSeqID]["CDR2-DNA"] !='' :
			pass
		elif FR2tail_pos and FR3head :
			self._dict[inSeqID].update({'CDR2-PRO':inProteinSeq[FR2tail_pos+1:FR3head_pos]})
			self._dict[inSeqID].update({'CDR2-DNA':translator.convertPROtoDNA(inDNAseq,FR2tail_pos+1,FR3head_pos)})
		#print "FR3head :" + FR3head
		#print FR3head_pos

		try:
			FR3tail_pos= (self._dict[inSeqID]["CDR3head_pos"]-1)/3 -1
		except: 
			FR3tail=self.find_high_PMW_score(self._pwmFR3tail,inProteinSeq[FR3head_pos:])
			FR3tail_pos=inProteinSeq.find(FR3tail)+len(FR3tail)-1
		if self._dict[inSeqID]["FR3-DNA"] !='' !='':
			pass
		elif FR3head_pos and FR3tail_pos :
			#new_FR3protein= inProteinSeq[FR3head_pos:FR3tail_pos]
			#new_FR3dna = translator.convertPROtoDNA(inDNAseq,,FR3head_pos, FR3tail_pos])
			self._dict[inSeqID].update({'FR3-PRO':inProteinSeq[FR3head_pos:FR3tail_pos+1]})
			self._dict[inSeqID].update({'FR3-DNA':translator.convertPROtoDNA(inDNAseq,FR3head_pos,FR3tail_pos+1)})

		try:
			FR4head_pos = self._dict[inSeqID]["CD3tail_pos"]/3
		except:
			FR4head=self.find_high_PMW_score(self._pwmFR4head,inProteinSeq[FR3tail_pos:])
			FR4head_pos=inProteinSeq.find(FR4head)
		new_CDR3protein= inProteinSeq[FR3tail_pos+1:FR4head_pos]
		new_FR3dna= translator.convertPROtoDNA(inDNAseq,FR3tail_pos+1,FR4head_pos)
		if FR3tail_pos and FR4head_pos and len(new_FR3dna)> len(self._dict[inSeqID]["CDR3-DNA"]):	
			self._dict[inSeqID].update({'CDR3-PRO':inProteinSeq[FR3tail_pos+1:FR4head_pos]})
			self._dict[inSeqID].update({'CDR3-DNA':translator.convertPROtoDNA(inDNAseq,FR3tail_pos+1,FR4head_pos)})
		#pdb.set_trace()
		FR4tail=self.find_high_PMW_score(self._pwmFR4head,inProteinSeq[FR3tail_pos:])
		FR4tail_pos=inProteinSeq.find(FR4tail)+len(FR4tail)-1
		if FR4head and FR4tail:
			self._dict[inSeqID].update({'FR4-PRO':inProteinSeq[FR4head_pos:FR4tail_pos+1]})
			self._dict[inSeqID].update({'FR4-DNA':translator.convertPROtoDNA(inDNAseq,FR4head_pos,FR4tail_pos+1)})
		elif FR4head:
			self._dict[inSeqID].update({'FR4-PRO':FR4head})
			self._dict[inSeqID].update({'FR4-DNA':translator.convertPROtoDNA(inDNAseq,FR4head_pos,len(inProteinSeq))})
		else:
			self._dict[inSeqID].update({'FR4-PRO':''})
			self._dict[inSeqID].update({'FR4-DNA':''})

		FVpro =self._dict[inSeqID]["FR1-PRO"]+self._dict[inSeqID]["CDR1-PRO"]+self._dict[inSeqID]["FR2-PRO"]+self._dict[inSeqID]["CDR2-PRO"]+self._dict[inSeqID]["FR3-PRO"]+self._dict[inSeqID]["CDR3-PRO"]+self._dict[inSeqID]["FR4-PRO"]
		self._dict[inSeqID].update({'FV-PRO': FVpro })
		FVdna=self._dict[inSeqID]["FR1-DNA"]+self._dict[inSeqID]["CDR1-DNA"]+self._dict[inSeqID]["FR2-DNA"]+self._dict[inSeqID]["CDR2-DNA"]+self._dict[inSeqID]["FR3-DNA"]+self._dict[inSeqID]["CDR3-DNA"]+self._dict[inSeqID]["FR4-DNA"]
                self._dict[inSeqID].update({'FV-DNA':FVdna})
		#self._dict[inSeqID].update({'DNAlen':len(FVdna)})


        def find_high_PMW_score (self,PMW_dict,fragment):
		#print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
		#print fragment
		if len(fragment)<len(PMW_dict):
			return "" 
		input_seq=fragment
		motif_size=len(PMW_dict)
                find_flag=False
                while (not find_flag) :
                        stop= min(len(input_seq),motif_size)
                        frag_window=input_seq[0:stop]
			backup_frag= None
			backup_score=  -1 
                        score=0
                        i=0
			while (i< motif_size) :
				if input_seq[i]in ["*","x"]:
					i+=1
					continue	
                                #print "PMW_dict[str(i)]" + str(PMW_dict[str(i)])
				#print "[input_seq[i]]" + input_seq[i]
				score= score+ float(PMW_dict[str(i)][input_seq[i]])
				#tmp_table=PMW_dict[str(i)]['E']
				#print "the current accumulate  score is at position"+str(i) +"is :" + str(score)	
				i+=1
                        #print "score of %s is %d" %(frag_window,score)
			if score >5:
				return frag_window
				break
				#return frag_window
			elif score <=5 and score > -1 and score >backup_score:
				find_flag =True
				backup_score =score
				backup_frag = frag_window
			elif len(input_seq)<=motif_size:
				break
			elif len(input_seq)>motif_size and score <=5:
                        	input_seq=input_seq[1:]

               	#print find_flag 
		#print frag_window
                #print fragment.find(frag_window)
		if find_flag == True:
			#print "Final Score of " + frag_window +' : ' + str(score)
			return  backup_frag 
		else: 
			#print "can't find"
			return "" 


###----------Test----------
#a={ 'M00680:164:000000000-AN2N4:1:2119:18838:25113': {'LID': 'TTCTGTAC', 'DNAlen': 357, 'FR1tail_pos': 78, 'DNA': 'GACATTCTGATGACCCAGTCTCCAGCTTCTTTGGCTGTGTCTCTAGGGCAGAGGGCCACCATCTCCTGCAAGGCCAGCCAAAGTGTTGATTATGATGGTGATAGTTATATGAACTGGTACCAACAGAAACCAGGACAGCCACCAAAACTCCTCATCTATGCTGCATCCAATCTAGAATCTGGGATCCCAGCCAGGTTTAGTGGCAGTGGGTCTGGGACAGACTTCACCCTCAACATCCATCCTGTGGAGGAGGAGGATGCTGCAACCTATTACTGTCAGCAAAGTAATGAGGATCCGTGGACGTTCGGTGGAGGCACCAAGCTGGAAATCAAACGGGCTGATGCTGCACCAACAGCATCCATCTTCCC', 'CDR3head_pos': 277, 'CDR2tail_pos': 168, 'CDR1tail_pos': 108, 'GERMLINE-J': 'JH3', 'FR1head_pos': 1, 'CDR3tail_pos': 296, 'FR2tail_pos': 159, 'PRO': 'DILMTQSPASLAVSLGQRATISCKASQSVDYDGDSYMNWYQQKPGQPPKLLIYAASNLESGIPARFSGSGSGTDFTLNIHPVEEEDAATYYCQQSNEDPWTFGGGTKLEIKRADAAPTASIF', 'CDR1head_pos': 79, 'GERMLINE-D': 'DSP2.9', 'FR2head_pos': 109, 'CDR2head_pos': 160, 'RID': 'GAGGCCCC', 'GERMLINE-V': 'J558.89pg.195','FR3head_pos':169,'FR3tail_pos':276}}

'''
a= {'test': { "DNA": 'CAGGTCCAACTGCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGATCCGGTATTACTACGGTAGTAGCTACCCCCCTGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAGCCAAAACAACACCCCCATCAGTCTATCCACTGGCCCCTGGGTGTGGAGATACAACTGGTTCCTCTGTGACTCTGGGATGCCTGGTCAAGGGCTACTTCCCTGAGTCAGTGACTGTGACTTGGAACTCTGGATCCCTGTCCAGCAGTGTGCACACCTTCCCAGCTCTCCTGCAGTCTGGACTCTACACTATGAGCAGCTCAGTGACTGTCCCCTCCAGCACCTGGCCAAGTCAGACCGTCACCTGCAGCGTTGCTCACCCAGCCAGCAGCACCACTCTCACAGTCTCCTCATTACAGGTTACCCATACGATGTTCCAGATTACGCT', 'LID': 'TTCTGTAC','CDR3head_pos': 53, 'CDR3tail_pos': 58, 'FR3head_pos':13, 'FR3tail_pos':52, 'PRO':'RSNCSSLTSEDSAVYYCARSGITTVVATPLYFDVWGTGTTVTVSSAKTTPPSVYPLAPGCGDTTGSSVTLGCLVKGYFPESVTVTWNSGSLSSSVHTFPALLQSGLYTMSSSVTVPSSTWPSQTVTCSVAHPASSTTLTVSSLQVTHTMFQIT'}}

b=AnnotateProtein(a,'mouse',"H")
b.AnnotateDict()
for ID,infor in a.iteritems():
	print "~~~ID is :" + ID 
	for key,value in infor.iteritems():
		if key.endswith("PRO") or key.endswith("DNA"):
			print key +":" + value
'''
