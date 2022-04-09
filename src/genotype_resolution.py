from src import predict_diplotype
import re, os
import pandas as pd
import numpy as np
from pybedtools import BedTool

def resolution(race, germline_vcf, outdir):
  ## 基于PharmGKB的bed文件过滤位点：将用户vcf中不在cpat.bed文件中的位点都删掉
  ## 提取用户的bed文件： 跳过以##开头的行
  vcf = []
  vcf_bed = []
  with open(germline_vcf, "r", encoding = "utf-8") as file:
    detected_pos = []
    for line in file:
      if 'CHROM' in line:
        colnames = line.strip().split('\t')
      if line[0] != '#':
        info = line.strip().split('\t')
        vcf.append(info)
        vcf_bed.append([info[0], info[1], info[1]])

  vcf_bed_df = pd.DataFrame(vcf_bed, columns=['chrom', 'start', 'end'])
  cpat_bed = pd.read_csv('./assets/bed/cpat.hla.bed', usecols=[0,1,2], sep="\t", names=['chrom', 'start', 'end'])
  filter_bed = BedTool.from_dataframe(vcf_bed_df).intersect(BedTool.from_dataframe(cpat_bed)).to_dataframe().drop_duplicates()

  ## 将user_vcf保存到vcf文件，以便diplotype函数调用
  vcf_df = pd.DataFrame(vcf, columns=colnames)
  vcf_df[colnames[1]] = vcf_df[colnames[1]].astype('int64')
  filter_bed = filter_bed.iloc[:, 0:2].rename(columns={'chrom': colnames[0], 'start': colnames[1]})
  filtered_vcf = pd.merge(vcf_df, filter_bed, how='inner', on=colnames[:2])
  fp = '%s/filtered.vcf' % outdir.strip('/')
  filtered_vcf.to_csv(fp, sep='\t', index=0)
  
  ## Class 1: Diplotype
  gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C19", "CYP2C9", "CYP2D6",
               "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15",
               "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]
  dic_diplotype = predict_diplotype.predict(fp, race, gene_list)
  os.system('rm %s' % fp)
  
  ## Class 2: HLA genes
  # hla_genes = ["HLA-B", "HLA-A", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1", "HLA-DQA1", "HLA-DRB3"]
  # hla_subtypes = vcf_df[vcf_df[colnames[0]].str.contains('^HLA', regex=True)][colnames[0]].drop_duplicates().to_list()
  # print(len(hla_subtypes))
  hla_subtypes = []
  
  ## Class 3: Genotypes of detected positions
  dic_rs2gt = {}
  format_index = colnames.index("FORMAT")
  # Reload cpat_bed
  cpat_bed = pd.read_csv('./assets/bed/cpat.hla.bed', sep="\t", names=['chrom', 'start', 'end', 'rsid'])
  cpat_bed_rsid = cpat_bed.dropna()
  for index, row in filtered_vcf.iterrows():
    info = row.to_list()
    format = info[format_index].split(":")
    gt_index = format.index("GT")
    gt = info[-1].split(":")[gt_index]
    genotype = predict_diplotype.phase_genotype(gt)
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


