from collections import Counter
import pandas as pd

#----------------------------------------
# Program to generate nucleotide composition statistics
# for whole mitogenome per species
#---------------------------------------

# list of raw fasta files
l = ['T_acuta_sequence.fasta','T_alcocki_sequence.fasta','T_arcuata_sequence.fasta','T_capricornis_sequence.fasta', 'T_paradussumieri_sequence.fasta',
'T_polita_sequence.fasta','T_rosea_sequence.fasta']

overall_table = {'Species':[], 'Length': [], "A":[], "C":[], "T":[], "G":[], "AT":[], "GC":[], "AT Skew": [], "GC Skew":[]}

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
                overall_table[nuc].append(seq.count(nuc)/len(seq) * 100)
            except:
                print("weird key", nuc)
        # add AT + GC percentages
        overall_table['AT'].append(( seq.count('A') + seq.count('T'))/len(seq) * 100)
        overall_table['GC'].append(( seq.count('G') + seq.count('C'))/len(seq) * 100)
        # add skews
        overall_table['AT Skew'].append((seq.count('A') - seq.count('T'))/(seq.count('A') + seq.count('T')))
        overall_table['GC Skew'].append((seq.count("G") - seq.count("C"))/(seq.count("G") + seq.count("C")))
        # add length
        overall_table["Length"].append(len(seq))
        print(seq)
# switch to panda and export
pd.DataFrame.from_dict(overall_table).to_csv("Whole_genome_stats3.csv")
print(overall_table)