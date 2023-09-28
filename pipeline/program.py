"""given:
amino acid string
sequence
"""

inputSeq = """
augaccggagtatcagattcatagtaagaaacgcacgaacaccgcggattccgcgctagacgt
gtcggctactcgcttaatctgataatggcaactggttcatgtcgagaaccgctcctatgc
tatgcattaaaatcacttagtttcgccatatactgcagatctcgtcctcatctgagagaa
acggacgttgtaatgccctagccgttcggtagatggtaatttaatgttcaaaaacctaag
gaggtgagacggaggttatccaccctcaagcgtgttgggatggtgcccggcaccaagaaa
ttgtctccgatcaatggtaccccgcaacttcgggtacccacaccctactccgcaggacgg
agacctgatcgaaacatagtactcgggtcgtgatgcggcttgtcgacgaggcgggcgctc
ggttagatccgcaatgtgggctggaggcaattaaatacgaagttcactattaaggcatgg
tgtgtcaggacgctccccgactaagcgtggagtattagaacggatagactaacctctcgg
gaacgaggatcagctcatcttgaccgacatcaggattggggcctattaacggtcggttca
ctgacggtttccgcagcgtaggacagttgtccaattgaaccacgaccgtgttgccacgaa
aagcagatccggacagtgagccgcagcagagttgaccattagcggatgggactccaagta
ttctccgtaacatacagctactcgatgggggtgtgaccctaccgcgctgcggcatcttgg
accgggaaaattatgtatcatgtacctgttcctctggcggggcctggtgagtagtccaaa
ttatattgaggcagtacgggtatgtacggccccctaatttgaaagccactcacgacatgg
atatggaaaggagggtcccccacgtataatcgatggatatactgttaactcatgataagt
gtcactgctctttcgaatcgaggtaagcgtacagtctaa
"""

protein = input()

stopCodons = {'UAA', 'UAG', 'UGA'}

revCodonTable = {
    'F': {'UUU', 'UUC'},  # Phe, phenylalanine
    'L': {'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'},  # Leu, leucine
    'S': {'UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'},  # Ser, serine
    'Y': {'UAU', 'UAC'},  # Tyr, tyrosine
    'C': {'UGU', 'UGC'},  # Cys, cysteine
    'W': {'UGG'},  # Trp, tryptophan
    'P': {'CCU', 'CCC', 'CCA', 'CCG'},  # Pro, proline
    'H': {'CAU', 'CAC'},  # His, histidine
    'Q': {'CAA', 'CAG'},  # Gln, glutamine
    'R': {'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'},  # Arg, arginine
    'I': {'AUU', 'AUC', 'AUA'},  # Ile, isoleucine
    'M': {'AUG'},  # Met, methionine
    'T': {'ACU', 'ACC', 'ACA', 'ACG'},  # Thr, threonine
    'N': {'AAU', 'AAC'},  # Asn, asparagine
    'K': {'AAA', 'AAG'},  # Lys, lysine
    'V': {'GUU', 'GUC', 'GUA', 'GUG'},  # Val, valine
    'A': {'GCU', 'GCC', 'GCA', 'GCG'},  # Ala, alanine
    'D': {'GAU', 'GAC'},  # Asp, aspartic acid
    'E': {'GAA', 'GAG'},  # Glu, glutamic acid
    'G': {'GGU', 'GGC', 'GGA', 'GGG'},  # Gly, glycine
}

sequences = []

# populate sequences with the starting codon 
for codon in revCodonTable[protein[0]]:
    sequences.append(codon)

for aa in protein:
    temp = []
    for codon in revCodonTable[aa]:
        sequencesCopy = [seq + codon for seq in sequences]
        temp.append(sequencesCopy)
    sequences = []
    for seqList in temp:
        sequences.extend(seqList)



print(len(sequences))




