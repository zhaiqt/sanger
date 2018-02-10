import re
import pdb
#sys.setrecursionlimit(1500)

class ClusterClone(object):

	def __init__(self, inAbdict):
		self._inlist = []
		#[seq1,seq2,seq3]
		for name,info in inAbdict.iteritems():
			self._inlist.append({
				'DNA': info['DNA'],
				'TRUNC_DNA': info['DNA'][22:],
				'RID': info['RID'],
				'LID': info['LID']
			})
		self._grouplist = []
		self._finallist = [] #[[seq2,seq1],[seq3,seq4]]
		self._consensusDict={}
		# {"consensus": {ab1:{'DNA':'','Pro':'','RID:'',LID:''}, ab2{}}}

	def mergeList(self):
		if len(self._inlist)<3:
			return
		self.divideGroup(self._inlist,'TRUNC_DNA')

		for ele in self._grouplist:
			self.checkConsensus(ele)

		for list in self._finallist:
			dna_list=[]
			LID_list=[]
			RID_list=[]
			for info in list:
				LID_list.append(info['LID'])
				RID_list.append(info['RID'])
			each_consensus =  self.calculate_consensus(list,)
			#print "each consensus :" + each_consensus
			self._consensusDict[each_consensus]={"RID":RID_list,'LID':LID_list}

	def checkConsensus(self, info_list):
		if not info_list or len(info_list) < 3: return

		consensus, second_position, second_freq_base = self.find_consensus_outlier(info_list)
		#consensus = self.calculate_consensus(info_list)
		if not second_freq_base:
			self._finallist.append(info_list)
			return
		else:
			new_info_list=self.outlier_consensus(consensus,info_list)
			self.checkConsensus(new_info_list)

	def outlier_consensus(self,consensus_seq,info_list):
		if not info_list: return
		if len(info_list) <= 5:
			self._finallist.append(info_list)
			return
		new_info_list=[]
		group = []
		for info in info_list:
			if self.diff_filter2(consensus_seq, info['TRUNC_DNA'], 3):
				group.append(info)
			else:
				new_info_list.append(info)
		if len(new_info_list)>=3 and len(info_list)-len(new_info_list) >=3:
			self._finallist.append(group)
			return new_info_list
		else:
			self._finallist.append(info_list)
			return

	def diff_filter2(self,consensus,s2,max_count):
		count=0
		i=0
		while i < len(consensus) and count < max_count:
			if consensus[i] != s2[i]:
				count +=1
			i +=1
		return count< max_count


	def freq_dict_of_lists(self,info_list,key):  # find_consensus uses this function
		n=max([len(info[key]) for info in info_list])
		#frequency_matrix={base:defaultdict(lambda:0) for base in 'ATGCN'  }
		frequency_matrix={base: {index: 0 for index in range(n)} for base in 'ATCGN' }
		for info in info_list:
			for i, base in enumerate(info[key]):
				if base in 'ATGCN':
					frequency_matrix[base][i] +=1
				elif base !='' or base != None:
					frequency_matrix['N'][i] +=1
		return frequency_matrix

	def calculate_consensus(self,info_list):
		#pdb.set_trace()
		frequency_matrix=self.freq_dict_of_lists(info_list,'DNA')
		consensus=''
		dna_length=len(frequency_matrix['A'])
		for i in range(dna_length):
			#if counter > 1: return (None, None, None)
			max_freq = -1
			max_freq_base= None
			for base in "ATGCN":
				if frequency_matrix[base][i]>max_freq:
					max_freq=frequency_matrix[base][i]
					max_freq_base=base
			consensus +=max_freq_base
		return consensus

	def find_consensus_outlier(self,info_list):
		frequency_matrix=self.freq_dict_of_lists(info_list,'TRUNC_DNA')
		consensus=''
		dna_length=len(frequency_matrix['A'])
		counter = 0
		second_freq = -1
		second_freq_base=None
		second_position= -1
		for i in range(dna_length):
			#if counter > 1: return (None, None, None)
			max_freq = -1
			max_freq_base= None
			base_outlier={}
			for base in "ATGCN":
				if frequency_matrix[base][i]>=3:
					base_outlier = {frequency_matrix[base][i]:base}
				if frequency_matrix[base][i]>max_freq:
					max_freq=frequency_matrix[base][i]
					max_freq_base=base
				#elif frequency_matrix[base][i] == max_freq:
				#	max_freq_base = "-"
			consensus += max_freq_base
			if max_freq_base == "-":	counter += 1

			try:
				clean_baseoutlier={freq: base_outlier[freq] for freq in base_outlier if base_outlier[freq] !=max_freq_base}
			except:
				pass
			if base_outlier:
				for k,v in clean_baseoutlier.iteritems():
					if second_freq< k:
						second_freq=k
						second_freq_base=v
						second_position=i
		return (consensus,second_position, second_freq_base)
	        # return consensus

	def divideGroup(self,info_list,keyword):
		if not info_list: return
		if len(info_list) <= 5:
			self._grouplist.append(info_list)
			return
		consensus,outlier_position, outlier_base=self.find_consensus_outlier(info_list)
		if consensus and (outlier_base is None or outlier_position ==-1):
			self._grouplist.append(info_list)
			return
		elif  consensus is None:
			return
		dna_length=len(consensus)
		outlier_list=[]
		left_list=list(info_list)
		for info in info_list:
				if info[keyword]=='':
					left_list.remove(info)
					continue
				#pdb.set_trace()
				try:
					if info[keyword][outlier_position] == outlier_base:
						outlier_list.append(info)
						left_list.remove(info)
				except:
					print "The problem dna sequence is :" + info[keyword]
					print "The length of the problem dna is :" + str(len(info[keyword]))
					print "the outlier_position is :" + str( outlier_position )
					print "the outlier_base is : " + outlier_base
		# return (left_list,outlier_list)
		if len(outlier_list)<3 or len(left_list)<3 :
			self._grouplist.append(info_list)
			return
		self.divideGroup(outlier_list,keyword)
		self.divideGroup(left_list,keyword)

	def memberCount_inGroup(self):
		count =0
		for list in self._finallist:
			count += len(list)
		return count

##############################################################################
'''
a={'test1':{'LID': 'ATCTAACT', 'DNAlen': 375, 'FR1tail_pos': 75, 'DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGTTTGGAAATATTAATCCAAGCAATGGTGGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGATCAGCAGCCTGACATATGAGGACTCTGCGGTCTATTATTGTGCAAAGAGGGAACTAATTACTACGGTAGTAGCTACGTACTACTTTGAAAACTGGGGCCAAGGAACCACACTCACAGTCTCCTCA', 'CDR3head_pos': 289, 'CDR2tail_pos': 174, 'FR1-PRO': 'QVQLQQPGTELVKPGASVKLSCKAS', 'CDR2-PRO': 'INPSNGGT', 'CDR3-DNA': 'GCAAAGAGGGAACTAATTACTACGGTAGTAGCTACGTACTACTTTGAAAAC', 'CDR1tail_pos': 99, 'FR3-DNA': 'ACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGATCAGCAGCCTGACATATGAGGACTCTGCGGTCTATTATTGT', 'FR4-PRO': 'WGQGTTLTVSS', 'FR1-DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCT', 'FR2-PRO': 'MHWVKQRPGQGLEWFGN', 'FR3tail_pos': 288, 'CDR1-DNA': 'GGCTACACCTTCACCAGCTACTGG', 'CDR1-PRO': 'GYTFTSYW', 'GERMLINE-J': 'JH2', 'FR1head_pos': 1, 'CDR3tail_pos': 292, 'FR4-DNA': 'TGGGGCCAAGGAACCACACTCACAGTCTCCTCA', 'FR2tail_pos': 150, 'CDR2-DNA': 'ATTAATCCAAGCAATGGTGGTACT', 'PRO': 'QVQLQQPGTELVKPGASVKLSCKASGYTFTSYWMHWVKQRPGQGLEWFGNINPSNGGTNYNEKFKSKATLTVDKSSSTAYMQISSLTYEDSAVYYCAKRELITTVVATYYFENWGQGTTLTVSS', 'CDR1head_pos': 76, 'GERMLINE-D': 'DFL16.1', 'FR3head_pos': 175, 'FR2head_pos': 100, 'CDR2head_pos': 151, 'CDR3-PRO': 'AKRELITTVVATYYFEN', 'FV-DNA': 'CAGGTCCAACTGCAGCAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGTTTGGAAATATTAATCCAAGCAATGGTGGTACTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGATCAGCAGCCTGACATATGAGGACTCTGCGGTCTATTATTGTGCAAAGAGGGAACTAATTACTACGGTAGTAGCTACGTACTACTTTGAAAACTGGGGCCAAGGAACCACACTCACAGTCTCCTCA', 'FV-PRO': 'QVQLQQPGTELVKPGASVKLSCKASGYTFTSYWMHWVKQRPGQGLEWFGNINPSNGGTTNYNEKFKSKATLTVDKSSSTAYMQISSLTYEDSAVYYCAKRELITTVVATYYFENWGQGTTLTVSS', 'RID': 'ACAGAAGA', 'FR3-PRO': 'TNYNEKFKSKATLTVDKSSSTAYMQISSLTYEDSAVYYC', 'GERMLINE-V': 'J558.53.146', 'FR2-DNA': 'ATGCACTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGTTTGGAAAT'}, 'test2':{'LID': 'TAACCTTA', 'DNAlen': 363, 'FR1tail_pos': 75, 'DNA': 'CAGGTCCAACTGCAGCAGCCTGGGCCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCTGACTCCATTTTCTCCACCTACTGGGTGCACTGGGTGAAGCAGAGGCCTGGTCATGGCCTTGTGTGTATTGGTTTGATTCATCCTAATTGTGGTAGTACTAACTACCATTAGAATTTCAAGAGCAAGGGCACAATTACTTTAGACAAATCGTCCAGCACAGACTACATGCAACTCAGCAGCATGACATCTGACGACTCAGAGGTAAAAAACTGTGGGTATAATTACAACGGAAGTAGAGGAGATATGGAATACTGGGGTCAAGGAACCACACTCACAGTCACCTCA', 'CDR3head_pos': 289, 'CDR2tail_pos': 174, 'FR1-PRO': 'QVQLQQPGPELVKPGASVKLSCKAS', 'CDR2-PRO': 'IHPNCGST', 'CDR3-DNA': 'GGGTATAATTACAACGGAAGTAGAGGAGATATGGAATAC', 'CDR1tail_pos': 99, 'FR3-DNA': 'ACTAACTACCATTAGAATTTCAAGAGCAAGGGCACAATTACTTTAGACAAATCGTCCAGCACAGACTACATGCAACTCAGCAGCATGACATCTGACGACTCAGAGGTAAAAAACTGT', 'FR4-PRO': 'WGQGTTLTVTS', 'FR1-DNA': 'CAGGTCCAACTGCAGCAGCCTGGGCCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCT', 'FR2-PRO': 'VHWVKQRPGHGLVCIGL', 'FR3tail_pos': 288, 'CDR1-DNA': 'GACTCCATTTTCTCCACCTACTGG', 'CDR1-PRO': 'DSIFSTYW', 'GERMLINE-J': 'JH4', 'FR1head_pos': 1, 'CDR3tail_pos': 289, 'FR4-DNA': 'TGGGGTCAAGGAACCACACTCACAGTCACCTCA', 'FR2tail_pos': 150, 'CDR2-DNA': 'ATTCATCCTAATTGTGGTAGTACT', 'PRO': 'QVQLQQPGPELVKPGASVKLSCKASDSIFSTYWVHWVKQRPGHGLVCIGLIHPNCGSTNYH*NFKSKGTITLDKSSSTDYMQLSSMTSDDSEVKNCGYNYNGSRGDMEYWGQGTTLTVTS', 'CDR1head_pos': 76, 'GERMLINE-D': 'DFL16.1', 'FR3head_pos': 175, 'FR2head_pos': 100, 'CDR2head_pos': 151, 'CDR3-PRO': 'GYNYNGSRGDMEY', 'FV-DNA': 'CAGGTCCAACTGCAGCAGCCTGGGCCTGAGCTGGTAAAGCCTGGGGCTTCAGTGAAGTTGTCCTGCAAGGCTTCTGACTCCATTTTCTCCACCTACTGGGTGCACTGGGTGAAGCAGAGGCCTGGTCATGGCCTTGTGTGTATTGGTTTGATTCATCCTAATTGTGGTAGTACTACTAACTACCATTAGAATTTCAAGAGCAAGGGCACAATTACTTTAGACAAATCGTCCAGCACAGACTACATGCAACTCAGCAGCATGACATCTGACGACTCAGAGGTAAAAAACTGTGGGTATAATTACAACGGAAGTAGAGGAGATATGGAATACTGGGGTCAAGGAACCACACTCACAGTCACCTCA', 'FV-PRO': 'QVQLQQPGPELVKPGASVKLSCKASDSIFSTYWVHWVKQRPGHGLVCIGLIHPNCGSTTNYH*NFKSKGTITLDKSSSTDYMQLSSMTSDDSEVKNCGYNYNGSRGDMEYWGQGTTLTVTS', 'RID': 'TGGAAAAA', 'FR3-PRO': 'TNYH*NFKSKGTITLDKSSSTDYMQLSSMTSDDSEVKNC', 'GERMLINE-V': 'J558.67.166', 'FR2-DNA': 'GTGCACTGGGTGAAGCAGAGGCCTGGTCATGGCCTTGTGTGTATTGGTTTG'}}
print len(a)
b=ClusterClone(a)
b.mergeList()
for l in b._finallist:
	print l
	print '\n'
print "Consensus are :"
print b._consensusDict
'''
