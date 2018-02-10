gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}


# a function to translate a single codon
def translate_codon(codon):
    return gencode.get(codon.upper(), 'x')

# a function to split a sequence into codons
def split_into_codons(dna, frame):
    if dna is None or not dna:
	return ''
    codons = []
    for i in range(frame - 1, len(dna)-2, 3):
        codon = dna[i:i+3]
        codons.append(codon)
    return codons

# a function to translate a dna sequence in a single frame
def translate_dna_single(dna, frame=1):
    if dna is  None  or not dna:
	return  ''
    codons = split_into_codons(dna, frame)
    amino_acids_seperate = ''
    for codon in codons:
        amino_acids_seperate = amino_acids_seperate +","+ translate_codon(codon)
    amino_acids_whole=amino_acids_seperate.replace(",","")
    #return amino_acids_whole+','+amino_acids_seperate
    return amino_acids_whole

# a function to translate a dna sequence in 3 forward frames
def translate_dna_3frame(dna):
    if dna is None  or not dna:
	return ''
    all_translations = []
    for frame in range(1,4):
        all_translations.append(translate_dna_single(dna, frame))
    return all_translations

def reverse_complement(dna):
	if dna is  None or not dna:
		return ''
	baseComplement={'A':"T", 'C':'G', 'T':'A','G':'C','N':'N', ' ':' ' }
	return ''.join(baseComplement[base] for base in reversed(dna))

def reverse_complement_FastQ(element):
	element[1]=reverse_complement(element[1])
	element[3].reverse()
	return element

def choose_translation(dna):
    if dna is None or not dna:
	return ''
    translate_protein_list=translate_dna_3frame(dna)
    translate_protein_list.extend(translate_dna_3frame((reverse_complement(dna))))
    protein_winner="*******"
    for protein in translate_protein_list:
        if (protein_winner.count("*")>= protein.count("*")):
            protein_winner=protein
    return protein_winner


def assembly(forwardSeq, backwardSeq,overlap_baseNumber=3,different_sense=True):
        # assemble forward and backward sequences without template file
	# input are forward sequence, backwardSeq, minimum overlap base, default is different sense strands
        """

        :rtype : object
        """
        if different_sense==True:
		secondSeq = reverse_complement(backwardSeq)
	else:
		secondSeq = backwardSeq
        overlap_position=forwardSeq.rfind(secondSeq[:overlap_baseNumber])
        while (overlap_position > -1):
                overlap_length=len(forwardSeq[overlap_position:])
                if secondSeq.startswith(forwardSeq[overlap_position:]):
                        return forwardSeq+secondSeq[overlap_length:]
                else:
                        overlap_position=forwardSeq.rfind(secondSeq[:overlap_baseNumber],0,overlap_position)
        if  (len(forwardSeq)>=len(secondSeq)):
                return forwardSeq
        else:
                return secondSeq


# find corresponding dna fragment, given the protein start and end position.#
# The DNA and Protein should be identical
def convertPROtoDNA(inDNAseq,pStat_pos,pEnd_pos):
                dnaStat_pos = pStat_pos * 3
                dnaEnd_pos  = pEnd_pos *3
                outputDNA = inDNAseq[dnaStat_pos:dnaEnd_pos]
                return outputDNA

###############
'''

a=["M00680:164:000000000-AN2N4:1:2119:14442:25148","AGACAAATCCTCCAGCACAGCCTACCTGTAGCTCAGCAGCCTGACATCTGAGGACACTGCCGTCTATTACTGTGCTAGAAAACGGGTATCCCATGCTATGGACTACTGGGGTCAAGGAACCTCAGTCACCGTCTCCTCATCAGCAAT","+",[38, 38, 37, 37, 37, 37, 30, 38, 38, 38, 38, 37, 37, 29, 27, 34, 38, 38, 38, 38, 38, 36, 32, 30, 34, 38, 38, 38, 37, 37, 38, 37, 35, 37, 23, 38, 38, 38, 38, 38, 38, 38, 38, 38, 37, 32, 37, 36, 32, 37, 36, 27, 27, 36, 11, 37, 37, 31, 37, 33, 34, 38, 37, 38, 37, 36, 34, 38, 38, 38, 38, 38, 38, 37, 36, 37, 24, 38, 38, 34, 38, 38, 38, 38, 38, 38, 38, 35, 34, 38, 37, 38, 37, 36, 34, 34, 38, 38, 38, 38, 38, 38, 38, 37, 37, 38, 38, 38, 37, 37, 37, 37, 38, 37, 38, 38, 38, 37, 38, 38, 38, 38, 38, 38, 38, 38, 37, 34, 36, 36, 38, 34, 38, 38, 38, 38, 35, 38, 38, 38, 37, 38, 34, 34, 34, 34, 34]]

print a
reverse_complement_FastQ(a)
print a
'''
