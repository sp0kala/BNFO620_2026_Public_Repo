from collections import Counter
import pandas as pd

l = ['T_acuta_sequence.fasta','T_alcocki_sequence.fasta','T_arcuata_sequence.fasta','T_capricornis_sequence.fasta', 'T_paradussumieri_sequence.fasta',
'T_polita_sequence.fasta','T_rosea_sequence.fasta']

overall_table = {'Species':[], 'Length': [], "A":[], "C":[], "T":[], "G":[]}
for file in l:
    with open("FASTA_Sequence/" + file) as f:
        seq = ""
        nucleotide_perc = ['A', 'C', 'G','T']
        for line in f:
            if line.startswith(">"):
                header = file.split(".")[0]
                overall_table['Species'].append(header)
            else:
                seq += line.strip()
        # compute a, t, c, g percentage
        for nuc in nucleotide_perc:
            try:
                overall_table[nuc].append(seq.count(nuc))
            except:
                print("weird key", nuc)
        # add length
        overall_table["Length"].append(len(seq))

# switch to panda and export
pd.DataFrame.from_dict(overall_table).to_csv("Whole_genome_stats1.csv")
print(overall_table)