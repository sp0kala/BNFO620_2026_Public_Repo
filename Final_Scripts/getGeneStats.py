import re
import pandas as pd
from collections import Counter

class gene:
    header_parser = re.compile(r'^>.*; (\d*)-(\d*); (\S); (\S*)')
    # group 1 = gene start
    # group 2 = gene end
    # group 3 = strand (+/-)
    # group 4 = gene name
    protein_coding = ["atp6", "atp8", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4", "nad4l",
                           "nad5", "nad6"]

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
        return self.gene_end - self.gene_start +1
    def getStrand(self):
        return self.gene_strand
    def startCodon(self):
        if self.gene_name in gene.protein_coding:
            return self.seq[0:3]
        else:
            return ""
    def stopCodon(self):
        if self.gene_name in gene.protein_coding:
            return self.seq[len(self.seq) - 3: len(self.seq)]
        else:
            return ""
    def getGC(self):
        g = self.seq.count('G')
        c = self.seq.count('C')
        return (g+c)/len(self.seq) *100
    def getAT(self):
        a = self.seq.count('A')
        t = self.seq.count('T')
        return (a+t)/len(self.seq) *100
    def getATSkew(self):
        a = self.seq.count('A')
        t = self.seq.count('T')
        return (a-t)/(a+t)
    def getGCSkew(self):
        g = self.seq.count('G')
        c = self.seq.count('C')
        return (g-c)/(g+c)

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

#### read in fasta files + genes
species = {
'arcuata': readFasta("T_arcuata_annotated.fasta", 'arcuata'),
'acuta': readFasta("T_acuta_annotated.fasta", 'acuta'),
'alcocki': readFasta("T_alcocki_annotated.fasta", 'alcocki'),
'capricornis': readFasta("T_capricornis_annotated.fasta", 'capricornis'),
'paradussumieri': readFasta("T_paradussumieri_annotated.fasta", 'paradussumieri'),
'polita': readFasta("T_polita_annotated.fasta", 'polita'),
'rosea': readFasta("T_rosea_annotated.fasta", 'rosea'),
'borealis': readFasta("G_borealis_annotated.fasta", 'borealis')
}

def gather_stats(organism):
    organism_stats = {}
    for g in organism.genes:
        organism_stats[g.gene_name] = {'location': (str(g.getStart()) + "-" + str(g.getEnd())),
                                       #'end': g.getEnd(),
                                       'strand': g.getStrand(),
                                       'length': g.getLength(),
                                           "start codon": g.startCodon(),
                                           "stop codon": g.stopCodon(),
                                           "GC%": g.getGC(),
                                       'AT%': g.getAT(),
                                       'AT Skew': g.getATSkew(),
                                       'GC skew': g.getGCSkew()}
    #print(organism.species_name, len(organism_stats))
    return pd.DataFrame.from_dict(organism_stats).transpose()

#species_files = ["arcuata", "acuta", "alcocki","capricornis", "paradussumieri", "polita", "rosea", "borealis"]
# = ["atp6", "atp8", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5","nad6"]

###########
# get stat files for each organism
###########
info = {}
for o,s in species.items():
    #print(o,s)
    data = gather_stats(s)
    data.to_csv(o+"_stats.csv")
