from matplotlib.pyplot import annotate
from src import genotype_resolution, clinical_annotation, pgx_report
import getopt, sys
import pandas as pd
import numpy as np
import os, re

gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",# "CYP2D6",
             "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15",
             "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]
########### Test the accuracy
race_dic = {'AFR': 'African_American_Afro_Caribbean', 'AMR': 'Latino', 'EAS': 'East_Asian', 'EUR': 'European'}
metadata = pd.read_csv('./test/population.txt', sep='\t')
test_vcfs = os.listdir('./test/88samples/')
outdir = './test/88samples_res/'

annotation_df = pd.DataFrame()
merged_score_df = pd.DataFrame()
for vcf in test_vcfs:
  if vcf.endswith('.vcf'):
    germline_vcf = './test/88samples/%s' % vcf
    # Get sample and race info
    basename = vcf.split('.')[0]
    race = metadata.loc[metadata['Get-RM 137 Samples']==basename, 'Superpopulation'].to_list()[0]
    race = race_dic[race]
    print(basename, race)
    ### Run CPAT
    dic_diplotype, dic_rs2gt, hla_subtypes = genotype_resolution.resolution(race, germline_vcf, outdir)
    # pgx_summary, clinical_anno_table = clinical_annotation.annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
    pgx_summary, clinical_anno_table, dosing_guideline_table = annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
    clinical_anno_table['Sample'] = basename
    #pgx_report.report(race, pgx_summary, clinical_anno_table, dosing_guideline_table, outdir, basename)
    annotation_df = pd.concat([annotation_df, clinical_anno_table], axis=0)
    merged_score_df = pd.concat([merged_score_df, pgx_summary], axis=0)

annotation_df.EvidenceLevel.drop_duplicates()
merged_score_df.groupby('EvidenceLevel').count()
merged_score_df[merged_score_df.EvidenceLevel == 'D']