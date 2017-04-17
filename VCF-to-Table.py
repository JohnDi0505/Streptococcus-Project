import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


# 1. Create a dictionary for the VCF file data
genomes_ncbi = SeqIO.to_dict(SeqIO.parse("NC_004368.gb", "gb")) # Load load reference genome downloaded from NCBI database
f = open("bubbles-in-sample.decomp.vcf")
genomes = {}
for line in f:
    if not re.match(r'^#', line, re.M|re.I):
        SNP = line.rstrip().split("\t")
        SNP[0] = SNP[0] + '.1'
        if genomes.has_key(SNP[0]):
            genomes[SNP[0]][int(SNP[1])] = SNP
        else:
            genomes[SNP[0]] = {}
            genomes[SNP[0]][int(SNP[1])] = SNP
f.close()


# 2. Parse VCF data to a varlist table
output_file = open('to-varlist.table', 'w')
header = ['var_id', 'acc', 'pos', 'var_type', 'refNT', 'altNT', 'refCodon', 'altCodon', 'refAA', 'altAA', 'locus', 'protein_accession', 'proj_id', 'gene_desc', 'call_batch']
output_file.write('\t'.join(header) + '\n')
for acc in genomes.keys():
    for gene in genomes_ncbi[acc].features:
        if re.search(r'RNA', gene.type, re.M|re.I) or re.search(r'CDS', gene.type, re.M|re.I):
           for snp in genomes[acc].keys():
                index = snp - 1
                var_type = re.search(r'SVTYPE=(.*)', genomes[acc][snp][7], re.M|re.I)
                List_output = [genomes[acc][snp][2], #0, var_id
                               acc[:-2], #1, acc
                               genomes[acc][snp][1], #2, pos
                               var_type.group(1), #3, var_type
                               genomes[acc][snp][3], #4, refNT
                               genomes[acc][snp][4], #5, altNT
                               '\N', #6, refCodon
                               '\N', #7, altCodon
                               '\N', #8, refAA
                               '\N', #9, altAA
                               '\N', #10, locus
                               '\N', #11, protein_accession
                               '2',  #12, proj_id
                               '\N', #13, gene_desc
                               '2']  #14, call_batch
                if index in gene: # if coding genes
                    List_output[10] = gene.qualifiers['locus_tag'][0]
                    List_output[13] = gene.qualifiers['product'][0]
                    if re.search(r'RNA', gene.type, re.M|re.I): # if RNA
                        List_output[13] = gene.type
                    else: # not RNA (if ORF)  
                        if re.search(r'SNP', genomes[acc][snp][7], re.M|re.I): # if ORF SNP (simple & complex)
                            if (index - gene.location.start) % 3 == 0:
                                refSeq = genomes_ncbi[acc][index] + genomes_ncbi[acc][index + 1] + genomes_ncbi[acc][index + 2]
                                altSeq = genomes[acc][snp][4] + genomes_ncbi[acc][index + 1] + genomes_ncbi[acc][index + 2]
                            elif (index - gene.location.start) % 3 == 1:
                                refSeq = genomes_ncbi[acc][index - 1] + genomes_ncbi[acc][index] + genomes_ncbi[acc][index + 1]
                                altSeq = genomes_ncbi[acc][index - 1] + genomes[acc][snp][4] + genomes_ncbi[acc][index + 1]
                            else:
                                refSeq = genomes_ncbi[acc][index - 2] + genomes_ncbi[acc][index - 1] + genomes_ncbi[acc][index]
                                altSeq = genomes_ncbi[acc][index - 2] + genomes_ncbi[acc][index - 1] + genomes[acc][snp][4]
                            refCodon = Seq(refSeq, IUPAC.unambiguous_dna)
                            altCodon = Seq(altSeq, IUPAC.unambiguous_dna)
                            if gene.location.strand == 1:
                                pass
                            else:
                                refCodon = refCodon.reverse_complement()
                                altCodon = altCodon.reverse_complement()
                            refAA = refCodon.translate()
                            altAA = altCodon.translate()
                            List_output[6] = str(refCodon)
                            List_output[7] = str(altCodon)
                            List_output[8] = str(refAA)
                            List_output[9] = str(altAA)
                            try:
                                List_output[11] = gene.qualifiers['protein_id'][0][:-2]
                            except:
                                pass
                        else: # not SNP (INDEL)
                            if re.search(r'INDEL_FROM_COMPLEX', genomes[acc][snp][7], re.M|re.I): # if COMPLEX INDEL
                                pass
                            else:
                                pass
                    output_file.write('\t'.join(List_output) + '\n')
                    del genomes[acc][snp]
    # Output IGS SNPs
    for snp_igs in genomes[acc]:
        var_type = re.search(r'SVTYPE=(.*)', genomes[acc][snp_igs][7], re.M|re.I)
        List_output = [genomes[acc][snp_igs][2], #0, var_id
                       acc[:-2], #1, acc
                       genomes[acc][snp_igs][1], #2, pos
                       var_type.group(1), #3, var_type
                       genomes[acc][snp_igs][3], #4, refNT
                       genomes[acc][snp_igs][4], #5, altNT
                       '\N', #6, refCodon
                       '\N', #7, altCodon
                       '\N', #8, refAA
                       '\N', #9, altAA
                       '\N', #10, locus
                       '\N', #11, protein_accession
                       '2', #12, proj_id
                       'IGS', #13, gene_desc
                       '2'] #14, call_batch
        output_file.write('\t'.join(List_output) + '\n')
output_file.close()


# 3. Parse VCF data to reforf table
output_file2 = open('to-reforf.table', 'w')
header = ['acc', 'locus', 'start', 'stop', 'strand', 'project', 'gene_desc', 'sequence']
output_file2.write('\t'.join(header) + '\n')

for acc in genomes.keys():
    for gene in genomes_ncbi[acc].features:
        if gene.type != 'gene' and gene.type != 'source' and gene.type != 'misc_binding' and gene.type != 'STS' and gene.type != 'repeat_region':
            List_output2 = [acc[:-2], #0, acc
                               '', #1, locus
                               '', #2, start
                               '', #3, stop
                               '', #4, strand
                               '2', #5, project
                               '', #6, gene_desc
                               ''] #7, sequence
            List_output2[1] = gene.qualifiers['locus_tag'][0]
            start_ex = re.search(r'(\d+)', str(gene.location.start + 1), re.M|re.I)
            end_ex = re.search(r'(\d+)', str(gene.location.end), re.M|re.I)
            List_output2[2] = start_ex.group(1)
            List_output2[3] = end_ex.group(1)
            if gene.location.strand == 1:
                List_output2[4] = 't'
            else:
                List_output2[4] = 'f'
            List_output2[6] = gene.qualifiers['product'][0]
            try:
                List_output2[7] = gene.qualifiers['translation'][0]
            except:
                pass
            output_file2.write('\t'.join(List_output2) + '\n')
output_file2.close()


# 4. Parse VCF data to var table
f = open("bubbles-in-sample.decomp.vcf")
genomes = {}
genome_id_list = {}
output_dic = {}
for line in f:
    if re.match(r'^#CHROM', line, re.M|re.I):
        gid_header = line.rstrip().split("\t")
        genome_id = 20
        field_id = 9
        for item in gid_header[9:]:
            gid = re.search(r'(^[a-zA-Z0-9]*)', item, re.M|re.I)
            genome_id_list[gid.group()] = [genome_id, field_id]
            genome_id += 1
            field_id += 1
    if not re.match(r'^#', line, re.M|re.I):
        SNP = line.rstrip().split("\t")
        SNP[0] = SNP[0] + '.1'
        if genomes.has_key(SNP[0]):
            genomes[SNP[0]][int(SNP[1])] = SNP
        else:
            genomes[SNP[0]] = {}
            genomes[SNP[0]][int(SNP[1])] = SNP
f.close()
for element in genome_id_list.keys():
    output_dic[element] = {}
output_var = open('to-var.table2', 'w')
header = ['genome_id', 'var_id', 'status', 'coverage', 'conf', 'call_batch']
output_var.write('\t'.join(header) + '\n')
for acc in genomes.keys():
    for snp in genomes[acc].keys():
        for gid in output_dic.keys():
            ln_gid_list = genome_id_list[gid]
            variant = re.search(r'(.*):(\d+),(\d+):(\d+.\d+)', genomes[acc][snp][ln_gid_list[1]], re.M | re.I)
            status = variant.group(1) # Append status
            if status != './.':
                output_dic[gid][genomes[acc][snp][2][4:]] = []
                output_dic[gid][genomes[acc][snp][2][4:]].append(str(ln_gid_list[0])) # Append genome_id
                output_dic[gid][genomes[acc][snp][2][4:]].append(genomes[acc][snp][2]) # Append var_id
                variant = re.search(r'(.*):(\d+),(\d+):(\d+.\d+)', genomes[acc][snp][ln_gid_list[1]], re.M | re.I)
                status = variant.group(1) # Append status
                if status == '0/0':
                    output_dic[gid][genomes[acc][snp][2][4:]].append('t') # '0/0' as 't'
                elif status == '1/1':
                    output_dic[gid][genomes[acc][snp][2][4:]].append('f') # '1/1' as 'f'
                if variant.group(2) == '0': # Append coverage
                    output_dic[gid][genomes[acc][snp][2][4:]].append(variant.group(3))
                else:
                    output_dic[gid][genomes[acc][snp][2][4:]].append(variant.group(2))
                output_dic[gid][genomes[acc][snp][2][4:]].append(variant.group(4)) # Append confidence
                output_dic[gid][genomes[acc][snp][2][4:]].append("2") # Append call_batch
                output_var.write('\t'.join(output_dic[gid][genomes[acc][snp][2][4:]]) + '\n')
output_var.close()
