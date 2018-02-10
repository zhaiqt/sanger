import translator
class ExtractBarcode:

    def __init__(self,pairfastQ,barcodePair):
	self._fastq1=pairfastQ[0]
	self._fastq2=pairfastQ[1]
	self._LID=barcodePair[0]
	self._RID=barcodePair[1]

    def extractBarcode(self, barcodeLen=8):
        if len(self._fastq1[1])>barcodeLen:
		self._LID=self._fastq1[1][:barcodeLen]
        else:
		self._LID=''

	if len(self._fastq2[1])>8:
		RID=self._fastq2[1][35:35+barcodeLen]
		self._RID=translator.reverse_complement(RID)
	else: 
		self._RID=''
        return (self._LID, self._RID)		

    
    def updateBarcode(self):
	if not self._fastq1[1].startswith(self._LID):
		self._LID=""

	if not self._fastq2[1].endswith(self._RID):
		self._RID=""
	return (self._LID, self._RID)
 
