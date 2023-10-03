class SeqGraph:
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

    def __init__(self, aaSeq):
        self.aaSeq = aaSeq
        if self.aaSeq[0] != 'M':
            self.aaSeq = 'M' + self.aaSeq

        self.codons = [SeqGraph.revCodonTable[aa] for aa in self.aaSeq]
        self.nodes = [[GraphNode(codon) for codon in codonSet] for codonSet in self.codons]
        self.firstNode = self.nodes[0][0]
        self.setNodeParents()
        self.setNodeChildren()

    def __str__(self):
        """Prints each layer as a row of GraphNode objects"""
        output = ''
        for layer in self.nodes:
            output += str(layer) + '\n'
        return output

    def setNodeParents(self):
        """Iterates through every codonnode and sets their parents (the nodes that point to them)"""
        for i in range(1, len(self.nodes)):
            prevLayerNodes = {node: False for node in self.nodes[i-1]}
            for node in self.nodes[i]:
                node.parentsVisited = prevLayerNodes

    def setNodeChildren(self):
        """Iterates through every codon node and sets their children nodes (the nodes they point to)"""
        for i in range(len(self.nodes)-1):  # iterate through layers
            nextLayerNodes = [node for node in self.nodes[i+1]]  # get all nodes in latter layer
            for node in self.nodes[i]:  # go through every node in current layer
                node.nextCodons = nextLayerNodes

    def validateNodes(self, inSeq: str, inNode=None, pathLen=0):
        """Passes an input sequence and validates each GraphNode's nextCodons list. An impossible
        to code for codon in a nextCodons list will get removed from said list. Any node no longer
        pointing to any nodes will be deleted and removed from all its former children's parentsVisited
        dictionary."""
        if not inNode:
            inNode = self.nodes[0][0]
        for nextNode in inNode.nextCodons:  # iterate over each child node
            nextSeq = nextNode.navigateSeq(inSeq)  # generate its return sequence

            print('codon:', inNode.codon, 'nextSeq:', nextSeq)
            if not nextSeq and pathLen != len(self.nodes):  # nucleotide seq has ran out, not end of graph
                inNode.nextCodons.remove(nextNode)  # remove nextNode as inNode's child
                nextNode.parentsVisited.pop(inNode)  # remove inNode as nextNode's parent
                print()
                continue

            nextNode.parentsVisited[inNode] = True

            self.validateNodes(inSeq=nextSeq,
                               inNode=nextNode,
                               pathLen=pathLen+1)
            
            



class GraphNode:
    def __init__(self, codon):
        self.codon = codon
        self.nextCodons = []
        self.parentsVisited = dict()

    def __repr__(self):
        return f'GraphNode({self.codon})'
    
    def navigateSeq(self, inSeq: str) -> str:
        """Returns the sequence after all characters for """
        for base in self.codon:
            try:
                inSeq = inSeq[inSeq.index(base) + 1:]
            except ValueError:
                return ''
        return inSeq


if __name__ == '__main__':
    aaSeq = 'MCGK'
    seq = SeqGraph(aaSeq=aaSeq)
    print()
    print(seq)
    for layer in seq.nodes:
        for node in layer:
            print(node, node.parentsVisited)

    seq.validateNodes('AUGUGUGUGUGUGUCGGAGGAGGA')


