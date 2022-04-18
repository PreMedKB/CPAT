"""
Evalute the diplotype population frequency by using Hardy-Weinberg equilibrium.
"""

import re
import pandas as pd

gene = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C19", "CYP2C9", "CYP2D6", 
        "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "IFNL3", 
        "MT-RNR1", "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]

race = ["African_American_Afro_Caribbean", "American", "Central_South_Asian", "East_Asian", "European", "Latino", "Near_Eastern", "Oceanian", "Sub_Saharan_African"]

for g in gene:
  print (g + "\n")
  # Definition
  gene_definition_file = "./assets/definition/%s_allele_definition_table.txt" % (g)
  definition = pd.read_csv(gene_definition_file, sep='\t', skiprows=5)
  haplotypes = definition.iloc[:, 0].to_list()
  # Frequency
  gene_frequency_file = "./assets/population_frequency/%s_frequency_table.txt" % (g)
  diplotype_frequency_file = "./assets/diplotype_frequency/%s_diplotype_frequency_table.txt" % (g)
  dic_race2dtfq = {}
  for i in range(9):
    locals() ['x' + str(i)] = {}
  
  with open(gene_frequency_file,'r',encoding = 'utf-8') as file:
    for line in file:
      if re.search('allele', line, re.IGNORECASE):
        line = line.replace(" ","_").replace("/","_").replace("-","_")
        #print(line.strip(), file = dt_fq)
        continue
      line = line.replace("\n","")
      info = line.split("\t")
      hap = info.pop(0)
      if hap not in haplotypes:
        continue
      
      for i in range(9):
        if info[i]:
          fre = float(info[i])
        else:
          fre = 0.000001
        # if fre > 0.000001:
        #   locals() ['x' + str(i)][hap] = fre
        # else:
        #   locals() ['x' + str(i)][hap] = 0.000001
        locals() ['x' + str(i)][hap] = fre
  
  for i in range(9):
    locals() ['d' + str(i)] = {}
  
  # Hardy–Weinberg equilibrium
  import itertools
  combination = list(itertools.combinations_with_replacement(sorted(locals() ['x' + str(i)]), 2))
  for i in range(9):
    for j in combination:
      diplotype = j[0] + '/' + j[1]
      if j[0]!= j[1]:
        # heterozygotes
        locals() ['d' + str(i)] [diplotype] = 2 * locals() ['x' + str(i)][j[0]] *  locals() ['x' + str(i)][j[1]]
      else:
        # homozygotes
        locals() ['d' + str(i)] [diplotype] = locals() ['x' + str(i)][j[0]] *  locals() ['x' + str(i)][j[1]]
  
  dic_race2dtfq = {
    "African_American_Afro_Caribbean" : d0,
    "American" : d1,
    "Central_South_Asian" : d2,
    "East_Asian" : d3,
    "European" : d4,
    "Latino" : d5,
    "Near_Eastern" : d6,
    "Oceanian" : d7,
    "Sub_Saharan_African" : d8,
  }
  
  race2dtfq = pd.DataFrame(dic_race2dtfq)
  race2dtfq['Global'] = race2dtfq.T.mean()
  race2dtfq.to_csv(diplotype_frequency_file, sep='\t')

