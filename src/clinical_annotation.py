import sqlite3
from matplotlib.pyplot import axis
import pandas as pd
import numpy as np
from pybedtools import BedTool


def annotation(dic_diplotype, dic_rs2gt, hla_subtypes):
  
  ## Connected database
  conn = sqlite3.connect("./assets/dbs/CPAT_v20220409.db")
  cursor = conn.cursor()
  
  ## Find therapies, extract the table
  ann = cursor.execute("SELECT ID, Gene, VariantOrHaplotype, GenotypeOrAllele FROM ClinAnn;")
  ann = cursor.fetchall()
  ann_df = pd.DataFrame(ann, columns=['ID', 'Gene', 'Variant', 'Alleles'])
  
  # Split the df by ' + ', and then all the genotype will be single-locus paired alleles or multi-locus genotype joint by '/'
  ann_df = ann_df.drop(['Alleles'], axis=1).join(ann_df['Alleles'].str.split(' \+ ', expand=True).stack().reset_index(level=1, drop=True).rename('Alleles')).reset_index(drop = True)
  
  # Part 1: Diplotype and haplotypes
  anno_ids_multi = []
  for gene in dic_diplotype.keys():
    #cpat_dip = set(dic_diplotype[gene])
    for cpat_dip in dic_diplotype[gene]:
      res_df = ann_df[ann_df.Gene == gene]
      for index, row in res_df.iterrows():
        if len(row.Alleles) == 2 and '*' not in row.Alleles:
          genotype = {row.Alleles[0], row.Alleles[1]}
        else:
          genotype = set(row.Alleles.split("/"))
        # Compare the genotype in PharmGKB with CPAT predicted.
        if genotype == set(cpat_dip):
          anno_ids_multi.append(row.ID)
  
  # Part 2: Single-locus based on rsIDs
  anno_ids_single = []
  for rsid in dic_rs2gt.keys():
    cpat_gt = set(dic_rs2gt[rsid])
    sub_df = ann_df[ann_df.Variant == rsid]
    for index, row in sub_df.iterrows():
      # HLA subtyping
      if set(row.Variant) <= set(hla_subtypes):
        anno_ids_single.append(row.ID)
      # SNPs/INDELs
      if len(row.Alleles) == 2:
        genotype = {row.Alleles[0], row.Alleles[1]}
      else:
        genotype = set(row.Alleles.split("/"))
      if genotype == cpat_gt:
        anno_ids_single.append(row.ID)
  
  #########################################
  # Fetch data from CPAT database
  res1 = cursor.execute("SELECT Gene, VariantOrHaplotype, Drug, Phenotypes, EvidenceLevel, Score, PhenotypeCategoryID, GenotypeOrAllele, Annotation, Function, URL, SpecialtyPopulation, GeneID, DrugID FROM ClinAnn WHERE EvidenceLevel != 3 AND EvidenceLevel != 4 AND ID IN (%s);" % ','.join([str(i) for i in anno_ids_multi]))
  res1 = cursor.fetchall()
  res1_df = pd.DataFrame(res1, columns=["Gene", "Variant", "Drug", "Phenotypes", "EvidenceLevel", "EvidenceScore", "PhenotypeCategoryID", "Alleles", "Annotation", "Function", "URL", "Pediatric", "GeneID", "DrugID"])
  res1_df['Class'] = 'Diplotype'

  res2 = cursor.execute("SELECT Gene, VariantOrHaplotype, Drug, Phenotypes, EvidenceLevel, Score, PhenotypeCategoryID, GenotypeOrAllele, Annotation, Function, URL, SpecialtyPopulation, GeneID, DrugID FROM ClinAnn WHERE EvidenceLevel != 3 AND EvidenceLevel != 4 AND ID IN (%s);" % ','.join([str(i) for i in anno_ids_single]))
  res2 = cursor.fetchall()
  res2_df = pd.DataFrame(res2, columns=["Gene", "Variant", "Drug", "Phenotypes", "EvidenceLevel", "EvidenceScore", "PhenotypeCategoryID", "Alleles", "Annotation", "Function", "URL", "Pediatric", "GeneID", "DrugID"])
  res2_df['Class'] = 'Single'
  res_df = pd.concat([res1_df, res2_df])
  res1_df.shape; res2_df.shape
  
  # Summary the data based on drug-phenotype category
  category = cursor.execute("SELECT * FROM PhenotypeCategoryDic;")
  category = cursor.fetchall()
  category_df = pd.DataFrame(category, columns=["PhenotypeCategoryID", "PhenotypeCategory"])
  res_df = pd.merge(res_df, category_df).drop(columns=['PhenotypeCategoryID'])
  
  response_score = {'Unknown': np.nan, 'Uncertain': np.nan, 'Unrelated': np.nan, 'Decreased': 0.5, 'Normal': 1, 'Moderate': 1, 'Increased': 2}
  score_df = pd.DataFrame({'Function': response_score.keys(), 'ResponseScore': response_score.values()})
  res_df = pd.merge(res_df, score_df)
  res_df.EvidenceScore = res_df.EvidenceScore.astype('float')
  
  # Output table 1: Original clinical annotation of PharmGKB
  clinical_anno_table = res_df[['Gene', 'Variant', 'Drug', 'Phenotypes', 'EvidenceLevel', 'Alleles', 'PhenotypeCategory', 'Annotation', 'Function', 'URL', 'Pediatric', 'Class']]

  # Summary Clinical Dosing Guideline
  dosing = cursor.execute("SELECT Source, Annotation, RelatedGeneID, RelatedDrugID, URL FROM ClinDosingGuideline;")
  dosing = cursor.fetchall()
  dosing = pd.DataFrame(dosing, columns=["DosingSource", "DosingAnnotation", "GeneID", "DrugID", "DosingURL"])
  dosing_df = pd.merge(res_df[['Gene', 'Variant', 'Alleles', 'Drug', 'GeneID', 'DrugID']].drop_duplicates(), dosing, on=["GeneID", "DrugID"])
  # Merge with drug table to get PAID
  drug = cursor.execute("SELECT ID, PAID FROM Drug;")
  drug = cursor.fetchall()
  drug = pd.DataFrame(drug, columns=["DrugID", "DrugPAID"])
  dosing_df = pd.merge(dosing_df, drug, on=["DrugID"]).drop(columns=["GeneID", "DrugID"])
  # Output table 2: Dosing Guidelines
  dosing_guideline_table = dosing_df#[['Gene', 'Variant', 'Drug', 'DrugPAID', 'DosingSource', 'DosingAnnotation', 'DosingURL']]
  
  ## Level A
  a1 = cursor.execute("SELECT DISTINCT GeneID, DrugID from ClinAnn WHERE EvidenceLevel = '1A';")
  a1 = cursor.fetchall(); a1 = pd.DataFrame(a1, columns = ['GeneID', 'DrugID'])
  a2 = cursor.execute("SELECT DISTINCT RelatedGeneID,RelatedDrugID FROM ClinDosingGuideline;")
  a2 = cursor.fetchall(); a2 = pd.DataFrame(a2, columns = ['GeneID', 'DrugID'])
  level_a = pd.concat([a1, a2], axis=0).drop_duplicates()
  a1.shape; a2.shape; level_a.shape # 119, 220, 237
  level_a['A'] = 1
  
  # Categorize by phenotypes and drugs
  pgx_summary = pd.DataFrame()
  categories = ['Toxicity', 'Dosage', 'Efficacy', 'Metabolism/PK', 'Other']
  for cat in categories:
    cat_df = res_df[res_df.PhenotypeCategory == cat]
    cat_df = pd.merge(cat_df, level_a, on=['GeneID', 'DrugID'], how='left')
    # cat_df['Drug+Gene'] = cat_df.Drug.str.cat(cat_df.Gene, sep='+') #cat_df.groupby("Drug")['Gene'].agg(lambda x: x.str.cat(sep=';'))
    # Calculating cat_pgx
    cat_pgx = cat_df.groupby("Drug")['ResponseScore', 'EvidenceScore', 'A'].mean()
    cat_pgx['PhenotypeCategory'] = cat
    cat_pgx['EvidenceLevel'] = ''; cat_pgx['Response'] = ''
    for index, row in cat_pgx.iterrows():
      if row.A > 0 and row.EvidenceScore >= 80:
        cat_pgx.loc[index, 'EvidenceLevel'] = 'A'
      else:
        cat_pgx.loc[index, 'EvidenceLevel'] = 'B'
      x = row.ResponseScore
      if x > 1:
        cat_pgx.loc[index, 'Response'] = 'Increased'
      elif x < 1:
        cat_pgx.loc[index, 'Response'] = 'Decreased'
      else:
        cat_pgx.loc[index, 'Response'] = 'Moderate'
    pgx_summary = pd.concat([pgx_summary, cat_pgx])
  
  cursor.close()
  conn.close()
  
  return(pgx_summary, clinical_anno_table, dosing_guideline_table)

