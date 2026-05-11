from readFasta import *
import pandas as pd


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

def gather_stats(organism, protein_coding):
    organism_stats = {}
    for g in organism.genes:
        #if g.gene_name in protein_coding:
        organism_stats[g.gene_name] = {'start': g.getStart(),
                                       'end': g.getEnd(),
                                       'length': g.getLength(),
                                           "start codon": g.startCodon(),
                                           "stop codon": g.stopCodon(),
                                           "GC%": g.getGC()}
    print(organism.species_name, len(organism_stats))
    return pd.DataFrame.from_dict(organism_stats).transpose()

species_files = ["arcuata", "acuta", "alcocki","capricornis", "paradussumieri", "polita", "rosea", "borealis"]

protein_coding = ["atp6", "atp8", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5",
                  "nad6"]



info = {}
for o,s in species.items():
    print(o,s)
    data = gather_stats(s, protein_coding)
    data.to_csv(o+"_stats.csv")

print(info)
#print(df)
#pd.DataFrame.from_dict(info).to_csv("all_stats.csv")






