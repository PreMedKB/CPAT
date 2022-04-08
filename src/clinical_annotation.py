import sqlite3
import pandas as pd
import numpy as np
from pybedtools import BedTool


def annotation(dic_diplotype, dic_rs2gt, hla_subtypes):
  
  ## Connected database
  conn = sqlite3.connect("./assets/dbs/CPAT_v20220317.db")
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
  res1 = cursor.execute("SELECT Gene, VariantOrHaplotype, Drug, Phenotypes, EvidenceLevel, Score, PhenotypeCategoryID, GenotypeOrAllele, Annotation, Function, URL, SpecialtyPopulation FROM ClinAnn WHERE EvidenceLevel != 3 AND EvidenceLevel != 4 AND ID IN (%s);" % ','.join([str(i) for i in anno_ids_multi]))
  res1 = cursor.fetchall()
  res1_df = pd.DataFrame(res1, columns=["Gene", "Variant", "Drug", "Phenotypes", "EvidenceLevel", "EvidenceScore", "PhenotypeCategoryID", "Alleles", "Annotation", "Function", "URL", "Pediatric"])
  res1_df['Class'] = 'Diplotype'

  res2 = cursor.execute("SELECT Gene, VariantOrHaplotype, Drug, Phenotypes, EvidenceLevel, Score, PhenotypeCategoryID, GenotypeOrAllele, Annotation, Function, URL, SpecialtyPopulation FROM ClinAnn WHERE EvidenceLevel != 3 AND EvidenceLevel != 4 AND ID IN (%s);" % ','.join([str(i) for i in anno_ids_single]))
  res2 = cursor.fetchall()
  res2_df = pd.DataFrame(res2, columns=["Gene", "Variant", "Drug", "Phenotypes", "EvidenceLevel", "EvidenceScore", "PhenotypeCategoryID", "Alleles", "Annotation", "Function", "URL", "Pediatric"])
  res2_df['Class'] = 'Single'
  res_df = pd.concat([res1_df, res2_df])
  
  # Summary the data based on drug-phenotype category
  category = cursor.execute("SELECT * FROM PhenotypeCategoryDic;")
  category = cursor.fetchall()
  category_df = pd.DataFrame(category, columns=["PhenotypeCategoryID", "PhenotypeCategory"])
  res_df = pd.merge(res_df, category_df).drop(columns=['PhenotypeCategoryID'])
  
  function_score = {'Unknown': np.nan, 'Uncertain': np.nan, 'Unrelated': np.nan, 
                    'Decreased': 0.5, 'Normal': 1, 'Moderate': 1, 'Increased': 2}
  score_df = pd.DataFrame({'Function': function_score.keys(), 'FunctionScore': function_score.values()})
  res_df = pd.merge(res_df, score_df)
  res_df.EvidenceScore = res_df.EvidenceScore.astype('float')
  res_df['CPATScore'] = res_df.EvidenceScore * res_df.FunctionScore
  res_df.to_csv('test.txt', sep='\t')
  
  # Output table 1
  clinical_anno_table = res_df[['Gene', 'Variant', 'Drug', 'Phenotypes', 'EvidenceLevel', 'Alleles', 'PhenotypeCategory', 'Annotation', 'Function', 'URL', 'Pediatric', 'Class']]
  
  # Categorize by phenotypes and drugs
  drug_score_summary = {}
  categories = ['Toxicity', 'Dosage', 'Efficacy', 'Metabolism/PK', 'Other']
  for cat in categories:
    cat_df = res_df[res_df.PhenotypeCategory == cat]
    # Calculating drug scores
    drug_score = cat_df.groupby("Drug").CPATScore.sum().sort_values(ascending = False)
    # Calculating ratio
    ratio = cat_df.groupby("Drug")['CPATScore', 'EvidenceScore'].mean()
    ratio['Ratio'] = ratio.CPATScore/ratio.EvidenceScore
    drug_score = ratio.CPATScore
    correspond_level =[]
    for v in drug_score.values:
      if v >= 80:
        correspond_level.append('1A')
      elif v >= 25 and v < 80:
        correspond_level.append('1B')
      elif v >= 8 and v < 25:
        correspond_level.append('2')
      elif v >= 0 and v < 8:
        correspond_level.append('3')
      elif v >= 0:
        correspond_level.append('4')
    drug_score_summary[cat] = drug_score.to_dict()
    drug_score_summary['%s.Level' % cat] = dict(zip(drug_score.keys(), correspond_level))
  
  # Output table 2
  drug_score_df = pd.DataFrame(drug_score_summary)
  
  # Output table 3
  levels = ['1A', '1B', '2', '3', '4'] # Row
  level_categories = ['Toxicity.Level', 'Dosage.Level', 'Efficacy.Level', 'Metabolism/PK.Level', 'Other.Level'] # Column
  
  table_res = []
  for i in levels:
    line = []
    for j in level_categories:
      tmp = drug_score_df[drug_score_df[j] == i].index.to_list()
      line.append('; '.join(tmp))
    table_res.append(line)
  
  pgx_summary = pd.DataFrame(table_res, columns=categories, index=levels)
  
  cursor.close()
  conn.close()
  
  return(clinical_anno_table, drug_score_df, pgx_summary)