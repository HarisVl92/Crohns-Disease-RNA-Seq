#!/usr/bin/env python3
# Χτίζει gene_id -> gene_name/biotype table από το ΙΔΙΟ GTF που μετρήσαμε (offline, ακριβές).
import re
GTF="/home/haris/Desktop/Bioinformatics_Project/data/reference/gtf/GRCh38_r110/Homo_sapiens.GRCh38.110.chr.gtf"
OUT="/home/haris/Desktop/Bioinformatics_Project/data/reference/gene_annotation.tsv"
gid_re=re.compile(r'gene_id "([^"]+)"'); gnm_re=re.compile(r'gene_name "([^"]+)"')
gbt_re=re.compile(r'gene_biotype "([^"]+)"')
n=0
with open(GTF) as f, open(OUT,"w") as o:
    o.write("gene_id\tgene_name\tbiotype\n")
    for line in f:
        if line.startswith("#"): continue
        c=line.split("\t")
        if len(c)<9 or c[2]!="gene": continue
        g=gid_re.search(c[8])
        if not g: continue
        gid=g.group(1)
        nm=gnm_re.search(c[8]); nm=nm.group(1) if nm else gid
        bt=gbt_re.search(c[8]); bt=bt.group(1) if bt else "NA"
        o.write(f"{gid}\t{nm}\t{bt}\n"); n+=1
print("genes:", n, "->", OUT)
