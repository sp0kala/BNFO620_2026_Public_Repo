# Program to create a fasta file of one gene from multiple species

# specify which gene to get
gene = "nad6"
# working directory
wd = "/Users/stellacastro/Downloads/BNFO620_2/BNFO620_2026_3/Aligned_Sequences/"
# path to output file
output_file = wd + "reannotated_genes/nad6_sequences.fasta"
# text file with file names of each species annotated fasta file
file_list = wd + "file_list.txt"

# clear output file if re-running the same gene
open(output_file, "w").close()
# loop through each file listed in the file_list
with open(file_list, "r") as fl:
    for file in fl:
        file_name = file.strip()
        if file_name == "":
            continue

        # get rid of _annotated.fasta part of file name (leaves just the species name)
        species_name = file_name.replace("_annotated.fasta", "")

        next = False
        # loop through individual file + write to output
        with open(wd + file_name, "r") as f, open(output_file, "a") as out:
            for line in f:
                line = line.rstrip()

                # isolate header line
                if line.startswith(">"):
                    # if specific gene is in headerline
                    if gene.lower() in line.lower():
                        out.write(">" + species_name + "_" + gene + "\n")
                        next = True # signal we found a header with the specified gene + can move on to get the sequence
                    else:
                        next = False # header not of specified gene = don't get the sequence below
                else:
                    if next: # just got the header with specified gene, now collect sequence
                        out.write(line + "\n")

