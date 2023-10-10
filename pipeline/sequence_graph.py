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
        self.setNodeChildren()

    def __str__(self):
        """Prints each layer as a row of GraphNode objects"""
        output = ''
        for layer in self.nodes:
            output += str(layer) + '\n'
        return output

    def setNodeChildren(self):
        """Iterates through every codon node and sets their children nodes (the nodes they point to)"""
        for i in range(len(self.nodes)-1):  # iterate through layers
            nextLayerNodes = [node for node in self.nodes[i+1]]  # get all nodes in latter layer
            for node in self.nodes[i]:  # go through every node in current layer
                node.nextNodes = nextLayerNodes

    def validateNodes(self, inSeq: str, inNode=None, pathLen=1):
        """Passes an input sequence and validates each GraphNode's nextCodons list. An impossible
        to code for codon in a nextCodons list will get removed from said list. Any node no longer
        pointing to any nodes will be deleted and removed from all its former children's parentsVisited
        dictionary.""" 
        if not inNode:
            inNode = self.nodes[0][0]
            inNode.visited = True
            inSeq = inNode.navigateSeq(inSeq)

        for nextNode in inNode.nextNodes:  # iterate over each child node
            if not nextNode:
                continue
            nextSeq = nextNode.navigateSeq(inSeq)  # generate its return sequence
            # print(f'Attempting {inNode} -> {nextNode} << {nextSeq} >> \t\t {pathLen}')

            if nextSeq == "END":
                nextNode.visited = True
                # print(f'END: {inNode} -> {nextNode}')
                continue

            if len(nextSeq) < 3 and pathLen < len(self.nodes):  # nucleotide seq has ran out, not end of graph
                # print(f'{inNode} -X-> {nextNode}, {pathLen}, {nextNode.nextNodes}')
                inNode.nextNodes[inNode.nextNodes.index(nextNode)] = None  # remove nextNode as inNode's child
                
            elif pathLen < len(self.nodes):
                nextNode.visited = True
                self.validateNodes(inSeq=nextSeq,
                                inNode=nextNode,
                                pathLen=pathLen+1)            

    def clearUnvisitedNodes(self):
        for i in range(len(self.nodes)):
            layer = self.nodes[i]
            validNodes = []
            for node in layer:
                if not node.visited:
                    continue
                validNodes.append(node)
                node.nextNodes = [nextNode for nextNode in node.nextNodes if nextNode]
            self.nodes[i] = validNodes

    def getNumCombinations(self):
        combos = 1
        for layer in self.codons:
            combos *= len(layer)
        return combos


class GraphNode:
    def __init__(self, codon):
        self.codon = codon
        self.nextNodes = []
        self.seq = []
        self.visited = False

    def __repr__(self):
        return f'GraphNode({self.codon})'
    
    def navigateSeq(self, inSeq: str) -> str:
        """Returns the sequence after all characters have been found """
        for base in self.codon:
            try:
                inSeq = inSeq[inSeq.index(base) + 1:]
            except ValueError:
                return ''
        if not inSeq:
            return 'END'
        return inSeq 


if __name__ == '__main__':
    with open('pipeline/seq.txt', 'r') as f:
        rnaSeq = f.read()
        print(rnaSeq)

    aaSeq = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
    
    aaSeq = "MVHL"

    seq = SeqGraph(aaSeq=aaSeq)
    print(True)

    seq.validateNodes(rnaSeq)
    seq.clearUnvisitedNodes()

    print(seq.getNumCombinations())

    print('POST VALIDATION:\n', seq)
