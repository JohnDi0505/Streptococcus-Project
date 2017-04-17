import re
import json
import psycopg2
import pandas as pd
import numpy as np
from pandas import DataFrame
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

con = None
con = psycopg2.connect(host, dbname, user, password) 
cur = con.cursor()
cur.execute("select b.genome_id, a.var_id, a.locus, a.var_type, a.refnt, a.altnt, a.refcodon, a.altcodon, a.refaa, a.altaa, a.pos, b.status, c.start, c.stop, c.strand from ((varlist2 as a inner join var2 as b on (a.var_id=b.var_id)) inner join reforf as c on (a.locus=c.locus)) where a.gene_desc != 'IGS' and a.var_type='SNP'")
x = cur.fetchall()
con.close()
df = DataFrame(x, columns=['genome_id', 'var_id', 'locus', 'var_type', 'refnt', 'altnt', 'refcodon', 'altcodon', 'refaa', 'altaa', 'pos', 'status', 'start', 'stop', 'strand'])
df = df.sort(columns='pos', ascending=True)
df = df.reset_index(drop=True)

genome_id_list = ["var_id", "pos", "locus", "strand", "snp_site", "gid_ref", "gid_22", "gid_23", "gid_20", "gid_21", "gid_24", "gid_25", "gid_29", "gid_31", "gid_26", "gid_27", "gid_28", "gid_30", "gid_32"]

SNP_table = DataFrame()
var_id_list = list(df.var_id.unique())
for var_id in var_id_list:
    temp = DataFrame(columns=genome_id_list, index=np.arange(1))
    samples = df[df.var_id == var_id]
    temp.var_id = samples.var_id.unique()[0]
    temp.pos = samples.pos.unique()[0]
    temp.locus = samples.locus.unique()[0]
    if samples.strand.unique()[0] == True:
        temp.strand = 'forward'
    else:
        temp.strand = 'reverse'
    temp.snp_site = (samples.pos.unique()[0] - samples.start.unique()[0]) % 3
    temp_ref = []
    temp_ref.append(samples.refcodon.unique()[0])
    temp_ref.append(samples.refaa.unique()[0])
    temp["gid_ref"][0] = temp_ref
    for i in samples.index:
        temp_sample = []
        if samples.status[i] == True:
            temp_sample.append('.')
            temp_sample.append('.')
        elif samples.status[i] == False:
            temp_sample.append(samples.altnt.unique()[0])
            if samples.refaa[i] == samples.altaa[i]:
                temp_sample.append('.')
            else:
                temp_sample.append(samples.altaa[i])
        temp["gid_" + str(samples.genome_id[i])][0] = temp_sample
    SNP_table = pd.concat([SNP_table, temp], axis=0)
SNP_table = SNP_table.reset_index(drop=True)
for col in ["gid_22", "gid_23", "gid_20", "gid_21", "gid_24", "gid_25", "gid_29", "gid_31", "gid_26", "gid_27", "gid_28", "gid_30", "gid_32"]:
    for row in SNP_table.loc[SNP_table[col].isnull(), col].index:
        SNP_table.at[row, col] = ["?", "?"]
SNP_to_json = SNP_table[genome_id_list].set_index('var_id')
SNP_to_json = SNP_to_json.T

SNP_to_json.to_json("/home/john/public_html/Streptococcus/js/SNP_table.json")
