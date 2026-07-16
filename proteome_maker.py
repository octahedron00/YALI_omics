import re

annot_file = "C:/Github/YALI_omics/Genome/Po1f_annot.gff3"
genome_file = "C:/Github/YALI_omics/Genome/Po1f.gff3"


exons = {}

with open("C:/Github/YALI_omics/Genome/Po1f_annot.gff3", "r") as f:
    
    for i, line in enumerate(f):
        if line.startswith('#'): continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9: continue
        if c[2] == 'exon':
            lt = re.search(r'locus_tag=([^;]+)', c[8])
            note = re.search(r'Note=([^;]+)', c[8])
            description = re.search(r'description=([^;]+)', c[8])
            print(f"{lt.group(1)}\t{note.group(1)}\t{description.group(1)}")





