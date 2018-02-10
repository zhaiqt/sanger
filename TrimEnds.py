import collections
import translator

class TrimEnds:

    def __init__(self,element):
        self._entry=list(element)
        # print self._entry[0]
        # print len(self._entry)

    # def __iter__(self):
    #     return self

    # def next(self):
    def output_trimed_fastq(self):

        return self._entry

        
    def trimEnds(self, sliding_window_size=10, minimum_quality=10):
        # trim 5' and 3' end
		# assert len(self._entry[1])==len(self._entry[3]),\
		# "**ERROR: The length of the sequence file does not match the length of the Phred file.\n\
		# Please check FastQ file near #%s and try again**" % (self._entry)
        remove_flag=True
        stop = min(len(self._entry[3]), sliding_window_size)+1
        while remove_flag:
            # Sliding window to detect poor quality region
            # remove 5' reads
            window_phred=self._entry[3][0:stop]
            bad_read_count=len(filter(lambda x: x- minimum_quality<0, window_phred))
            if (bad_read_count>=1 or "N" in self._entry[1][0:stop] or 'n' in self._entry[1][0:stop]):
                self._entry[1]=self._entry[1][1:]
                self._entry[3]=self._entry[3][1:]
            else:
				remove_flag=False

        remove_flag=True
        while remove_flag:
            # stop = max(0,len(self._entry[3])-sliding_window_size)
            window_phred=self._entry[3][(-stop+1):]
            bad_read_count=len(filter(lambda x: x- minimum_quality<0, window_phred))
            if bad_read_count>=1 or "N" in self._entry[1][-stop:] or 'n' in self._entry[1][-stop:]:
                self._entry[1]=self._entry[1][:-1]
                self._entry[3]=self._entry[3][:-1]
            else:
                remove_flag=False
		#print ">" + self._entry[0] +"\n"+self._entry[1]
		return
		#return (self._entry[1],self._entry[3])

    def trim3End(self, sliding_window_size=6, minimum_quality=10):

                #only trim 3' end
                assert len(self._entry[1])==len(self._entry[3]),\
                "**ERROR: The length of the sequence file does not match the length of the Phred file.\n\
                Please check FastQ file near #%s and try again**" % (self._entry)
                #print "IT is trying to remove low score"
                remove_flag=True

                while remove_flag:
                        stop = max(0,len(self._entry[3])-sliding_window_size)
                        window_phred=self._entry[3][stop:]
                        bad_read_count=len(filter(lambda x: x- minimum_quality<0, window_phred))
                        if bad_read_count>=1:
                                self._entry[1]=self._entry[1][:-1]
                                self._entry[3]=self._entry[3][:-1]
                        else:
                                remove_flag=False
                #print ">" + self._entry[0] +"\n"+self._entry[1]
                return self._entry
                #return (self._entry[1],self._entry[3])


################
'''
a = ['@W1135_HC_Ti', 'GNNNNGNNNNNTCTGGNGGAGGCTTAGTGAGCCTGGAGGGTCCCTGAAACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTGACTATGCCCTGTCTTGGGTTCGCCAGACTCCGGAGAAGAGGCTGGAGTGGGTCGCAGCCATTAGTGATGGTGGTAGTTCCACCTACTACCCAGACACTGTGAAGGGTCGATTCACCATCTCCAGAGACAATGCCAAGAACACCCTGTACCTGCAAATGAGCCGTCTGAGGTCTGAGGACACGGCCATGTATTACTGTGCAAGAAATCTGGATTATAGATACGACGGCCCCGCCTGGTTTACTTACTGGGGCCAAGGGACTCTGGTCACTGTCTCTGCAGCGTCGACTTCGCAA', '+', [8, 4, 4, 6, 4, 9, 5, 4, 5, 5, 4, 12, 13, 10, 46, 30, 5, 12, 15, 15, 7, 25, 17, 52, 23, 20, 21, 58, 54, 38, 9, 43, 62, 36, 62, 62, 47, 29, 62, 37, 26, 62, 62, 47, 62, 62, 59, 62, 56, 56, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 43, 62, 62, 62, 62, 62, 62, 62, 62, 59, 43, 62, 59, 13, 54, 49, 62, 62, 44, 59, 56, 62, 62, 59, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 59, 59, 62, 59, 59, 62, 62, 62, 62, 62, 59, 59, 59, 62, 59, 62, 62, 59, 59, 62, 62, 49, 42, 42, 49, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 54, 62, 62, 54, 59, 62, 62, 59, 59, 62, 62, 49, 62, 62, 62, 62, 62, 62, 62, 54, 62, 62, 62, 62, 59, 62, 62, 59, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 59, 54, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 49, 49, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 47, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 49, 62, 62, 10, 12, 35, 62, 62, 62, 56, 62, 62, 31]]

t = TrimEnds(a)
print t.trimEnds()
'''
