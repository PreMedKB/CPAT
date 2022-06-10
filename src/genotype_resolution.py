from src import predict_diplotype
import re
import pandas as pd
import numpy as np
from pybedtools import BedTool

def resolution(race, germline_vcf):
  ## Filter loci based on PharmGKB's bed file: delete all loci in the user's vcf that are not in the cpat.bed file
  ## Extract the user's bed file: skip the lines starting with '##'
  vcf = []
  vcf_bed = []
  with open(germline_vcf, "r", encoding = "utf-8") as file:
    for line in file:
      if 'CHROM' in line:
        colnames = line.strip().split('\t')
      if line[0] != '#':
        info = line.strip().split('\t')
        vcf.append(info)
        vcf_bed.append([info[0], info[1], info[1]])
  
  vcf_bed_df = pd.DataFrame(vcf_bed, columns=['chrom', 'start', 'end'])
  cpat_bed = pd.read_csv('./assets/bed/cpat.hla.bed', sep="\t", names=['chrom', 'start', 'end', 'rsid'])
  # Replace chr
  vcf_bed_df['chrom'] = vcf_bed_df['chrom'].map(lambda x: re.sub('chr|Chr|CHR', '', x)).astype('str')
  cpat_bed['chrom'] = cpat_bed['chrom'].map(lambda x: re.sub('chr|Chr|CHR', '', x)).astype('str')
  filter_bed = BedTool.from_dataframe(vcf_bed_df).intersect(BedTool.from_dataframe(cpat_bed.iloc[:,0:3])).to_dataframe().drop_duplicates()
  
  ## Save the user_vcf to a vcf file for the diplotype function to call
  vcf_df = pd.DataFrame(vcf, columns=colnames)
  vcf_df.loc[:,'#CHROM'] = vcf_df['#CHROM'].map(lambda x: re.sub('chr|Chr|CHR', '', x)).astype('str')
  vcf_df.loc[:,colnames[1]] = vcf_df[colnames[1]].astype('int64')
  filter_bed = filter_bed.iloc[:, 0:2].rename(columns={'chrom': colnames[0], 'start': colnames[1]})
  filtered_vcf = pd.merge(vcf_df, filter_bed, how='inner', on=colnames[:2])
  
  ## Class 1: Diplotype
  gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",\
              "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15", "IFNL3", \
              "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]
  dic_diplotype = predict_diplotype.predict(filtered_vcf, race, gene_list)
  ## Class 2: HLA genes
  hla_subtypes = []
  ## Class 3: Genotypes of detected positions
  dic_rs2gt = {}
  format_index = colnames.index("FORMAT")
  # Reload cpat_bed
  cpat_bed_rsid = cpat_bed.dropna()
  for index, row in filtered_vcf.iterrows():
    info = row.to_list()
    format = info[format_index].split(":")
    gt_index = format.index("GT")
    gt = info[-1].split(":")[gt_index]
    if re.findall('0', gt) == ['0', '0']:
      genotype = 0
    elif re.findall('0', gt) == ['0']:
      genotype = 1
    else:
      genotype = 2
    if row[0].startswith('HLA') and genotype != 0:
      hla_subtypes.append(row[0])
    # If the variant was within the clinical relevant list, add it into dis_rs2gt
    tmp = cpat_bed_rsid[(cpat_bed_rsid.chrom == info[0]) & (cpat_bed_rsid.start == info[1])].rsid.to_list()
    if tmp != []:
      rsids = tmp
    elif info[2] in cpat_bed_rsid.rsid.to_list(): # The chromosome positions of a rsID may not complete.
      rsids = [info[2]]
    else:
      rsids = None
    
    if rsids:
      ref = info[3]
      alts = info[4].split(',')
      alleles = [ref] + alts
      var = []
      for g in re.split('\||/', gt):
        var.append(alleles[int(g)])
      for rsid in rsids:
        dic_rs2gt[rsid] = tuple(var)
  
  return(dic_diplotype, dic_rs2gt, set(hla_subtypes))