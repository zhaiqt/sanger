class Parse2FastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath1,filePath2,headerSymbols=['@','+']):
        '''
	Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
	if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        elif filePath.endswith('.fastq'):
            self._file = open(filePath, 'rU')
        '''
	self._currentLineNumber1 = 0
	self._currentLineNumber2 = 0
        self._hdSyms = headerSymbols
	if filePath1.endswith('.gz'):
		self._file1=gzip.open(filePath1)
	elif filePath1.endswith('.fastq'):
		self._file1 = open(filePath1, 'rU')
	else:
	#	print filePath1
	#	print filePath1.endswith('.fastq')
		print "Error: The 1 FASTQ can't be found."

	if filePath2.endswith('.gz'):
                self._file2=gzip.open(filePath1)
        elif filePath2.endswith('.fastq'):
                self._file2 = open(filePath2, 'rU')
	else:
		print "Error: The 2 FASTQ can't be found."
	
 
    def __iter__(self):
        return self
     
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList1= []
        for i in range(4):
            line = self._file1.readline()
            self._currentLineNumber1 += 1 ## increment file position
            if line:
                elemList1.append(line.strip('\n'))
            else: 
                elemList1.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList1].count(True)
        nones = elemList1.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber1)
        # -- Make sure we are in the correct "register" --
        assert elemList1[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber1) 
        assert elemList1[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber1) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList1[1]) == len(elemList1[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber1) 
        
	########  read the file of the 2nd 
       # ++++ Get Next Four Lines ++++
        elemList2= []
        for i in range(4):
            line = self._file2.readline()
            self._currentLineNumber2 += 1 ## increment file position
            if line:
                elemList2.append(line.strip('\n'))
            else:
                elemList2.append(None)

        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList2].count(True)
        nones = elemList2.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber2)
        # -- Make sure we are in the correct "register" --
        assert elemList2[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber2)
        assert elemList2[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber2)
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList2[1]) == len(elemList2[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber2)

	scoreQ1=self.covert_Qscore(elemList1[3])	
 	elemList1[3]=scoreQ1
	
        scoreQ2=self.covert_Qscore(elemList2[3])
        elemList2[3]=scoreQ2
	#print (elemList1,elemList2)

        # ++++ Return fatsQ data as tuple ++++
	return (elemList1,elemList2)



    def covert_Qscore(self,String):
	Q_score=[]
	for i in String:
		Q_score.append(int((ord(i))-33))
		
        #Q_score ([ord(x)-33 for x in String])
	return Q_score

####################
'''
a='/dlab/NGS/usem-seqanalysis/160314_zhaiqi1_miseq_HBx52-60DNA.20160214_AN2N4/test/R1.fastq'

b='/dlab/NGS/usem-seqanalysis/160314_zhaiqi1_miseq_HBx52-60DNA.20160214_AN2N4/test/R2.fastq'

c=Parse2FastQ(a,b)
print c.next()
print "~~~~~~~ 2nd next~~~~~~~"
print c.next()
print "~~~~~~~ 3nd next~~~~~~~"
print c.next()

'''
