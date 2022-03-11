import predict_diplotype
import pandas as pd
import numpy as np
import os, re

gene_list = ("ABCG2", "CACNA1S", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "MT-RNR1", "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1")

########### Test the accuracy
pharmgkb_100genome = {'AFR': 'African_American_Afro_Caribbean', 'AMR': 'Latino', 'EAS': 'East_Asian', 'EUR': 'European'}
metadata = pd.read_csv('./test/population.txt', sep='\t')
test_vcfs = os.listdir('./test/88samples/')
result = {}
for vcf in test_vcfs:
  if vcf.endswith('.vcf'):
    vcf_file = './test/88samples/%s' % vcf
    # Get sample and race info
    sample = vcf.split('.')[0]
    race = metadata.loc[metadata['Get-RM 137 Samples']==sample, 'Superpopulation'].to_list()[0]
    race = pharmgkb_100genome[race]
    print(sample, race)
    result[sample] = predict_diplotype.main(vcf_file, race, gene_list)

test_df = pd.DataFrame(result)

### Validation data
getrm = pd.read_csv('./test/getrm_raw.txt', sep='\t', index_col=0)
# Filter samples and genes
overlap_gene = list(set(gene_list).intersection(set(getrm.columns.to_list())))
samples = list(result.keys())
getrm_new = getrm[getrm.index.isin(samples)][overlap_gene].T.replace(np.nan, '-')


# # samples used in PharmCAT
# pharmcat_samples = []
# for i in metadata['PharmCAT 59 Samples'].to_list():
#   if i != '-' and pd.isna(i) is False:
#     pharmcat_samples.append(i)


# Accuracy by gene
for g in overlap_gene:
  test = test_df.loc[g, samples].to_list()
  validate = getrm_new.loc[g, samples].to_list()
  # remove the nan
  index = [x for x, y in list(enumerate(validate)) if y != '-']
  right = 0
  for i in index:
    if g in ['ABCG2', 'VKORC1']:
      x = test[i].split('/')
      genotype = validate[i].replace('G', 'C').replace('A', 'T')
      yy = [genotype[0], genotype[1]]
    else:
      if g == 'DPYD':
        test[i] = test[i].replace('Reference', '*1')
      x = re.findall(r'\*\w+', test[i])
      y = re.findall(r'\*\w+', validate[i])
      yy = [j.replace('*0', '*') for j in y]
      # if g == 'DPYD':
      #   print(x, yy, test[i], validate[i])
    # if g == 'VKORC1':
    #   print(x, yy, test[i], validate[i])
    if set(x) <= set(yy):
      right = right+1
    else:
      # Has the overlap of haplotype is passed
      if 'no consensus' in validate[i]:
        y = re.findall(r'\*\w+', validate[i])
        yy = [j.replace('*0', '*') for j in y]
        if len(set(x).intersection(set(yy))) > 0:
          right = right+1
  
  print(g, round(right/len(index)*100, 2))

