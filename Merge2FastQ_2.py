import collections
import translator
import logging

outfilePath='/Users/zhaiqi1/Documents/Novartis/my_code/ToolBox/sanger/results/'

log_file = outfilePath+'runlog.txt'
log_level = logging.DEBUG
logging.basicConfig(filename=log_file, level=log_level, format='%(asctime)s %(message)s')

class Merge2FastQ():

    def __init__(self,pairQ2,overlap_baseNumber=10):
	self._merged=[]
	self._fastQ1=pairQ2[0]
	self._fastQ2=translator.reverse_complement_FastQ(pairQ2[1]) # from now on, fastQ are sense strand
	self._overlapNumber =overlap_baseNumber
    self._assemble_seq=""
    #def __iter__(self):
        #return self

    #def next(self):
    def merge2FastQ(self):
        # assemble forward and backward sequences without template file
        # input are forward sequence, backwardSeq, minimum overlap base, default is different sense strands

        assemble_seq = ''
        assemble_seq = self.anneal_1direction(self._fastQ1, self._fastQ2)

        #print "length of sequence F: %d" % len(self._fastQ1[1])

        #print "length of sequence R: %d" % len(self._fastQ2[1])
        #print "length of assembled: %d" % len(assemble_seq)
        if len(assemble_seq)<max(len(self._fastQ1[1]),len(self._fastQ2[1])):
            logging.info(self._fastQ1)
            logging.info(self._fastQ2)
            logging.info ("assembled seq: " +assemble_seq)

        #print '========'
        #print assemble_seq
        '''
        if self.compare_end_score():
    		assemble_seq= self.anneal_1direction(self._fastQ1, self._fastQ2)
        else:
            antiSense_fastQ1= translator.reverse_complement_FastQ(self._fastQ1)
            antiSense_fastQ2= translator.reverse_complement_FastQ(self._fastQ2)
            antiSense_assemble_seq= self.anneal_1direction(antiSense_fastQ2, antiSense_fastQ1)
            assemble_seq=translator.reverse_complement(antiSense_assemble_seq)
        '''
        return	[self._fastQ1[0],assemble_seq]

    '''
    def compare_end_score(self):
	sumq1=sum(self._fastQ1[3][-self._overlapNumber:])
	sumq2=sum(self._fastQ2[3][:self._overlapNumber])
	if sumq2 >= sumq1:
		return True
	else:
		return False
    '''
    def anneal_1direction(self, inputfastq1, inputfastq2):
        print "$$$$$$"
        print inputfastq1
        print inputfastq2
        firstSeq=inputfastq1[1]
        secondSeq=inputfastq2[1]
        fastQ1= inputfastq1
        fastQ2 = inputfastq2
        print len(fastQ1[1])
        print len(fastQ2[1])
        if not firstSeq and not secondSeq:
                print " firstSeq and not secondSeq empty"
                return ' '
        elif not firstSeq or len(firstSeq) <= self._overlapNumber or firstSeq in secondSeq:
                print "lif not firstSeq or len(firstSeq) <= self._overlapNumber"
                return secondSeq
        elif not secondSeq or len(secondSeq) <= self._overlapNumber or secondSeq in firstSeq:
                print "firtseq> 2"
                return firstSeq


        couple_list=[]
        print "-------- caculate assemb1--------"
        (overlap_seq1,assembleSeq1) = self.find_overlap(fastQ1,fastQ2)
        couple_list.append((overlap_seq1,assembleSeq1))
        print "-------- caculate assemb2--------"
        (overlap_seq2,assembleSeq2) =self.find_overlap(fastQ2,fastQ1)
        couple_list.append((overlap_seq2,assembleSeq2))
        print "-------- caculate assemb3--------"
        antiSense_fastQ2= translator.reverse_complement_FastQ(inputfastq2)
        (overlap_seq3,assembleSeq3) = self.find_overlap(fastQ1,antiSense_fastQ2)
        couple_list.append((overlap_seq3,assembleSeq3))

        longest_overlap=('','')
        print couple_list
        for tmp_couple in couple_list:
            if len(tmp_couple[0]) > len(longest_overlap):
                longest_overlap = tmp_couple
        assembled_fasta = longest_overlap[1]

        print "assembled_fasta: "+ assembled_fasta
        self._assemble_seq = assembled_fasta
        return assembled_fasta



    def shrinkID(self):
        longName1=self._fastQ1[0].split()
        self._fastQ1[0]=longName1[0].lstrip('@')

        longName2=self._fastQ2[0].split()
        self._fastQ2[0]=longName2[0].lstrip('@')
        return

    def find_overlap(self,inputfastq1,inputfastq2):
        firstSeq=inputfastq1[1]
        secondSeq=inputfastq2[1]
        fastQ1= inputfastq1
        fastQ2 = inputfastq2

        overlap_position=firstSeq.rfind(secondSeq[-self._overlapNumber:])
        #print secondSeq[-self._overlapNumber:]
        #print overlap_position   # Seq A

        overlap_seq = ''
        if (overlap_position > -1):
            overlap_length=overlap_position+self._overlapNumber
            overlapAQ=[fastQ1[1][:overlap_length],fastQ1[3][:overlap_length]]
            overlapBQ=[fastQ2[1][-overlap_length:],fastQ2[3][-overlap_length:]]
            overlap_seq=self.overlap_consensus(overlapAQ, overlapBQ)
            print "AQ: " + overlapAQ[0]
            print "BQ: " + overlapBQ[0]
            print "overlap" +overlap_seq
            #overlap_seq:
            assembleSeq=fastQ2[1][:-overlap_length]+ overlap_seq + fastQ1[1][overlap_length:]
            print "5': " +  fastQ2[1][:-overlap_position]
            print "3': " +fastQ1[1][overlap_length:]
            print len(assembleSeq)
            print "length of assembled sequence is %d."  %len(assembleSeq)
            print "assembled seq is " + assembleSeq
            return (overlap_seq,assembleSeq)
        else:
            return ('','')



    def overlap_consensus(self,overlapQA,overlapQB):	  # QA is (DNAsequence1, Qscore1), QB[DNA sequence2, Qscore]
        mismatchCounts=0
        matchCounts=0
        overlapSeq=""
    	AQ=overlapQA
    	BQ=overlapQB
        #print "AQ[0]: "+ AQ[0]
        #print "BQ[0]: " + BQ[0]
        for i in range(min(len(AQ[0]),len(BQ[0]))):
            if AQ[0][i]==BQ[0][i]:
                overlapSeq +=AQ[0][i]
                matchCounts +=1
            else:
                mismatchCounts += 1

                if AQ[1][i]>BQ[1][i]:
                    overlapSeq +=AQ[0][i]
                else:
                    overlapSeq +=BQ[0][i]
        if matchCounts >= mismatchCounts *3:
            return overlapSeq
        else:
            #print "Warning!!! can't find overlap."
            #return overlapQB[0]

            return " "



########################
'''
a=["M00680:164:000000000-AN2N4:1:2119:14442:25148","CAGGTCCACTTGCAGCAGCCTGGGGCTGAGCTGGTGAGGCCTGGGGTTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACATTCACCAGCTACTGGATGCACTGGATTAAGCAGAGGCCTGAGCAAGGCCTTGAGAGGATTGGAGAGATTAATCCTAGCAATGGTGGTACTAACTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAATCCTCCAGCACAGCCTACCTGTAGCTCAGCAGCCTGACATCTGAGGACACTGCCG","+",[37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 36, 37, 38, 38, 37, 38, 38, 37, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 34, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 34, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 37, 38, 37, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 38, 38, 37, 37, 38, 38, 38, 38, 38, 38, 38, 35, 38, 36, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 32, 37, 38, 35, 27, 37, 37, 38, 38, 35, 38, 38, 38, 38, 38, 38, 38, 38, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 34, 37, 34, 37, 37, 38, 38, 37, 37, 38, 37, 37, 38, 38, 38, 38, 34, 38, 38, 37, 38, 38, 38, 37, 28, 37, 37, 38, 37, 37, 37, 38, 38, 37, 38, 38, 34, 38, 34, 38, 37, 37, 34]]
#b=["M00680:164:000000000-AN2N4:1:2119:14442:25148","AGACAAATCCTCCAGCACAGCCTACCTGTAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTGCTAGAAAACGGGTATCCCATGCTATGGACTACTGGGGTCAAGGAACCTCAGTCACCGTCTCCTCATCAGCAAT","+",[38, 38, 37, 37, 37, 37, 30, 38, 38, 38, 38, 37, 37, 29, 27, 34, 38, 38, 38, 38, 38, 36, 32, 30, 34, 38, 38, 38, 37, 37, 38, 37, 35, 37, 23, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 32, 37, 36, 32, 37, 36, 27, 27, 36, 11, 37, 37, 31, 37, 33, 34, 38, 37, 38, 37, 36, 34, 38, 38, 38, 38, 38, 38, 37, 36, 37, 24, 38, 38, 34, 38, 38, 38, 38, 38, 38, 38, 35, 34, 38, 37, 38, 37, 36, 34, 34, 38, 38, 38, 38, 38, 38, 38, 37, 37, 38, 38, 38, 37, 37, 37, 37, 38, 37, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 37, 34, 36, 36, 38, 34, 38, 38, 38, 38, 35, 38, 38, 38, 37, 38, 34, 34, 34, 34, 34]]

b=['M00680:164:000000000-AN2N4:1:2119:14442:25148', 'ATTGCTGATGAGGAGACGGTGACTGAGGTTCCTTGACCCCAGTAGTCCATAGCATGGGATACCCGTTTTCTAGCACAGTAATAGACGGCAGTGTCCTCAGATGTCAGGCTGCTGAGCTACAGGTAGGCTGTGCTGGAGGATTTGTCT', '+', [34, 34, 34, 34, 34, 38, 37, 38, 38, 38, 35, 38, 38, 38, 38, 34, 38, 36, 36, 34, 37, 38, 38, 38, 38, 38, 38, 38, 38, 37, 38, 38, 38, 37, 38, 37, 37, 37, 37, 38, 38, 38, 37, 37, 38, 38, 38, 38, 38, 38, 38, 34, 34, 36, 37, 38, 37, 38, 34, 35, 38, 38, 38, 38, 38, 38, 38, 34, 38, 38, 24, 37, 36, 37, 38, 38, 38, 38, 38, 38, 34, 36, 37, 38, 37, 38, 34, 33, 37, 31, 37, 37, 11, 36, 27, 27, 36, 37, 32, 36, 37, 32, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 23, 37, 35, 37, 38, 37, 37, 38, 38, 38, 34, 30, 32, 36, 38, 38, 38, 38, 38, 34, 27, 29, 37, 37, 38, 38, 38, 38, 30, 37, 37, 37, 37, 38, 38]]
'''

'''
a=['@W1180_HC_Ti', 'TTGTCCTGTAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATAAATTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAAATATTTATCCTCCTGAAAGTTATACTAACTACAATCAAAAGTTCAAGGACAAGGTCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCCGACATCTGAGGACTCTGCGGTCTATTACTGTACAAGAGAGTATGGTAACCTTGACTACTGGGGCCAAGGCACCACTCTCACAGTCTCCTCA', '+', [62, 54, 47, 44, 62, 62, 62, 62, 62, 62, 59, 49, 43, 59, 62, 62, 62, 62, 62, 59, 62, 59, 59, 59, 51, 62, 62, 62, 62, 59, 41, 62, 49, 59, 62, 62, 59, 62, 59, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 49, 62, 59, 59, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 59, 59, 62, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 54, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 42, 54, 62, 62, 54, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 54, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 56, 49, 38, 54, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 49, 49, 62, 62, 62, 62, 62, 49, 49, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62]]


b=['@W1180_HC_Ti', 'TCTGGGGCTGAGCTGGTGAGGCCTGGGGCTTCAGTGAAGTTGTCCTGTAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATAAATTGGGTGAAGCAGAGGCCTGGACAAGGCCTTGAGTGGATCGGAAATATTTATCCTCCTGAAAGTTATACTAACTACAATCAAAAGTTCAAGGACAAGGTCACATTGACTGTAGACAAATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCCGACATCTGAGG', '+', [30, 62, 62, 62, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 49, 49, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 59, 62, 62, 51, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 54, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 54, 62, 62, 62, 62, 62, 62, 59, 62, 62, 59, 62, 62, 62, 62, 62, 59, 59, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 59, 59, 62, 62, 62]]

c= ['@W1208_HC_Ti', 'GCCTGGGGGTTCTCTGCGACTCTCCTGTGCAACTTCTGGGTTC', '+', [39, 62, 56, 44, 62, 62, 62, 62, 36, 43, 62, 59, 46, 62, 49, 59, 49, 34, 59, 41, 51, 62, 62, 62, 62, 59, 54, 62, 62, 59, 34, 59, 59, 62, 62, 62, 59, 62, 62, 59, 62, 62, 59]]

d=['@W1208_HC_Ti', 'CAATAATAAGTGGCACTGTCCTCAGCTCTCAGGGTGTTCATTTGAAGATAGAGGATGCTTTGGGAATTATCTCTGGAGATGGTGAACCGACCCTTCACAGAAGGACTGTATTCTGTTGTGTAACCATCAGCTATGTTTCGAATAAAACCCACCCACTCAAGTGTCTTTCCTGGAGGCTGGCGGACCCAGGTCATGTAGTAATCAGTGAAGGTGAACCCAGAAGTTGCACAGGAGAGTCGCAGAGAACCCCCAGGCTGTACCAAGTCTCCTCCAGACTCC', '+', [44, 39, 59, 59, 62, 59, 62, 62, 59, 59, 59, 62, 59, 56, 62, 59, 62, 54, 49, 42, 62, 62, 59, 59, 59, 62, 62, 62, 62, 59, 59, 59, 62, 50, 34, 56, 34, 34, 39, 30, 43, 62, 54, 43, 59, 49, 46, 62, 59, 59, 59, 62, 62, 62, 59, 36, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 54, 54, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 59, 62, 62, 59, 59, 62, 62, 59, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 51, 62, 54, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 59, 56, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 49, 62, 62, 43, 62, 51, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 49, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 56, 62, 51, 56, 38, 31]]
'''

'''
c=['@W1193_HC_Ti', 'GCTGCTGGATACACCTTCACTAACTACTGGATAGGTTGGATATACCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGATATTTACCCTGGAGGTGATTATACTAACTACAATGAGAAATTCAAGGTCAAGGCCACACTGACTGCAGACACATCGTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCCATCTATTACTGTGCAAGAGGTGGCCCTTCCTACGCTTACTGGGGCCAGGGG', '+', [41, 27, 39, 54, 62, 62, 62, 49, 44, 51, 56, 59, 51, 62, 62, 62, 49, 59, 43, 62, 59, 62, 59, 62, 59, 54, 62, 62, 62, 59, 59, 62, 62, 59, 59, 62, 62, 62, 62, 54, 49, 62, 62, 62, 62, 59, 59, 49, 54, 38, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 59, 62, 62, 62, 59, 54, 42, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 59, 62, 62, 59, 54, 62, 62, 62, 62, 62, 62, 62, 62, 59, 51, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 54, 62, 62, 62, 59, 62, 62, 54, 54, 62, 62, 62, 56, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 54, 62, 62, 62, 62, 59, 43, 59, 62, 62, 62, 62, 54, 54, 50, 59, 62, 62, 42, 54, 62, 62, 62, 62, 50, 62, 62, 62, 62, 56, 56, 62, 62, 62, 62, 59, 62, 62, 62, 62, 49, 62, 62, 51, 62, 62, 62, 62, 59, 62, 59, 62, 62, 62, 62, 62, 56, 62, 59, 59, 62, 62, 62, 62, 62, 62, 59, 62, 54, 42, 54, 62, 62, 62, 59, 62, 62, 62, 62, 46, 59, 50, 45, 62, 48, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 56, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 59, 59, 59, 59, 62, 62, 62]]

d=['@W1193_HC_Ti', 'CTGGAGCTGAGCTGGTAAGGCCTGGGACTTCAGTGAAAATGTCCTGCAAGGCTGCTGGATACACCTTCACTAACTACTGGATAGGTTGGATATACCAGAGGCCTGGACATGGCCTTGAGTGGATTGGAGATATTTACCCTGGAGGTGATTATACTAACTACAATGAGAAATTCAAGGTCAAGGCCACACTGACTGCAGACACATCGTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCCATCTATTACTGTGCAAGAGG', '+', [56, 46, 59, 56, 56, 59, 43, 59, 62, 62, 62, 51, 62, 62, 54, 62, 62, 62, 62, 62, 62, 51, 59, 62, 59, 62, 62, 51, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 54, 62, 51, 59, 62, 62, 62, 62, 62, 62, 51, 62, 62, 59, 51, 59, 62, 59, 62, 59, 62, 62, 62, 54, 62, 62, 62, 59, 49, 59, 59, 62, 51, 62, 62, 62, 59, 54, 54, 62, 62, 59, 62, 59, 62, 59, 62, 62, 62, 62, 59, 62, 59, 62, 54, 56, 62, 59, 59, 62, 62, 62, 59, 59, 59, 56, 62, 59, 56, 62, 62, 62, 62, 62, 59, 62, 54, 62, 62, 54, 49, 49, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 59, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 59, 56, 59, 62, 62, 56, 62, 59, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 56, 54, 62, 62, 56, 59, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 59, 59, 59, 62, 54, 62, 62, 62, 62, 62, 54, 54, 62, 51, 59, 59, 62, 62, 59, 62, 59, 62, 62, 51, 62, 62, 62, 62, 62, 62, 62, 62, 49, 49, 59, 46, 59, 59, 59, 49, 59, 59, 62, 62, 62, 62, 56, 62, 51, 59, 59, 54, 62, 59, 59, 49, 59, 59, 59, 62, 59, 56, 36, 43, 56, 44, 54, 62, 56, 62, 62, 62, 59]]




me=Merge2FastQ([c,d])
# print len(a[-1])
# print len(b[-1])
#c.find_overlap(a,b)
tmp= me.merge2FastQ()

print tmp
print len(tmp[1])
'''
