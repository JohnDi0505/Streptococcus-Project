import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Create a Python dictionary for reference genome
genome_ncbi = SeqIO.to_dict(SeqIO.parse("NC_004368.gb", "genbank"))
gene_dict = {}
for gene in genome_ncbi['NC_004368.1'].features:
    if gene.type != 'gene' and gene.type != 'source' and gene.type != 'misc_binding' and gene.type != 'STS' and gene.type != 'repeat_region':
        try:
            locus = gene.qualifiers['locus_tag'][0]
            gene_dict[locus] = str()
            gene_dict[locus] = str(genome_ncbi['NC_004368.1'][int(gene.location.start):int(gene.location.end)].seq)
        except:
            pass

file = open("upload_nt.txt", "w")
for key in gene_dict.keys():
    file.write("update reforf set ntseq='%s' where locus='%s';\n" % (gene_dict[key], key))
file.close()
