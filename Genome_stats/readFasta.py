import re
from collections import Counter

class gene:
    header_parser = re.compile(r'^>.*; (\d*)-(\d*); (\S); (\S*)')
    # group 1 = gene start
    # group 2 = gene end
    # group 3 = strand (+/-)
    # group 4 = gene name

    def __init__(self, header, sequence):

        self.gene_name = gene.header_parser.search(header).group(4)
        self.gene_start = int(gene.header_parser.search(header).group(1))
        self.gene_end = int(gene.header_parser.search(header).group(2))
        self.gene_strand = gene.header_parser.search(header).group(3)

        self.seq = sequence

    def getStart(self):
        return self.gene_start
    def getEnd(self):
        return self.gene_end
    def getLength(self):
        return self.gene_end - self.gene_start

    def startCodon(self):
        return self.seq[0:3]
    def stopCodon(self):
        return self.seq[len(self.seq)-3: len(self.seq)]

    def getGC(self):
        g = self.seq.count('G')
        c = self.seq.count('C')
        return (g+c)/len(self.seq) *100


class readFasta:

    def __init__(self, inputFile, name):
        self.species_name = name
        self.genes = []
        try:
            with open(inputFile, 'r') as f:
                sequence = ''
                created_header = False
                for line in f:
                    line = line.rstrip()

                    if line.startswith('>'):
                        if created_header:
                            #print("creating gene: ", header, "seq: ", sequence)
                            self.genes.append(gene(header, sequence))
                            sequence = ''
                        header = line
                        created_header = True
                    else:
                        sequence += line
                # last gene
                if created_header:
                    #print("creating last gene:", header, "seq:", sequence)
                    self.genes.append(gene(header, sequence))

        except FileNotFoundError:
            print("File could not be opened")

        #print(readFasta.genes)

