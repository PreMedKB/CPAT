import pandas as pd
import numpy as np
import re, json

gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",
             "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15", "IFNL3", 
             "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]
######################
##### definition #####
######################
output = {}
for gene in gene_list:
  allele_definition_table = "./assets/definition/%s_allele_definition_table.txt" % gene
  define_df = pd.read_csv(allele_definition_table, sep='\t', index_col=0)
  single_gene = {}
  single_gene['reference'] = 'GRCh38'
  single_gene['chr'] = re.findall('chromosome (\w+),', define_df.index[1])[0]
  ng = define_df.iloc[1,]
  rs = define_df.iloc[3,].to_list()
  pos_rs = []
  for i in range(0, len(ng)):
    ng_item = ng[i]
    matchobj = re.search(r'\w\.(\d+)(\w*)', ng_item)
    if matchobj:
      pos = matchobj.group(1)
    else:
      matchobj = re.search(r'\w\.(\d+)(\_\d+)?(del|ins)(\w*)', ng_item)
      if matchobj:
        pos = matchobj.group(1)
      else:
        pos = ''
    pos_rs.append('%s:%s' % (pos, rs[i]))
  
  # Definition of Haplotype
  haplotype = {}
  reference = define_df.iloc[5,].to_list()
  for index, row in define_df.iloc[5:,].iterrows():
    hap_name = index
    defined_base = []
    for i in range(0, len(row)):
      if row[i] is np.nan:
        defined_base.append(reference[i])
      else:
        defined_base.append(row[i])
    haplotype[index] = dict(zip(pos_rs, defined_base))
  
  single_gene['haplotype'] = haplotype
  output[gene] = single_gene


with open("./assets/allele_definition_merge.json", 'w', encoding='utf-8') as f:
  json.dump(output, f, ensure_ascii=False, indent=4)

json.loads(open("./assets/allele_definition_merge.json").read())
###########################
##### check indel/dup #####
###########################
for gene in gene_list:
  output[gene]['haplotype']