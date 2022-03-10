"""
Evalute the diplotype population frequency by using Hardy-Weinberg equilibrium.
"""

import re

gene = ["ABCG2", "CACNA1S", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "MT-RNR1", "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]

race = ["African_American_Afro_Caribbean", "American", "Central_South_Asian", "East_Asian", "European", "Latino", "Near_Eastern", "Oceanian", "Sub_Saharan_African"]

for g in gene:
  print (g + "\n")
  #gene_definition_file = "%s_allele_definition_table.txt" % (g)
  gene_frequency_file = "./population_frequency/%s_frequency_table.txt" % (g)
  diplotype_frequency_file = "./diplotype_frequency/%s_diplotype_frequency_table.txt" % (g)
  dic_race2dtfq = {}
  for i in range(9):
    locals() ['x' + str(i)] = {}
  dt_fq = open(diplotype_frequency_file, "w", encoding='utf-8')
  with open(gene_frequency_file,'r',encoding = 'utf-8') as file:
    for line in file:
      if re.search('allele', line, re.IGNORECASE):
        line = line.replace(" ","_").replace("/","_").replace("-","_")
        print(line.strip(), file = dt_fq)
        continue
      line = line.replace("\n","")
      info = line.split("\t")
      hap = info.pop(0)
      
      for i in range(9):
        if info[i]:
          fre = float(info[i])
        else:
          fre = 0.000001
        if fre > 0.000001:
          locals() ['x' + str(i)][hap] = fre
        else:
          locals() ['x' + str(i)][hap] = 0.000001

  for i in range(9):
    locals() ['d' + str(i)] = {}

  for i in range(9):
    for j in locals() ['x' + str(i)]:
      for k in locals() ['x' + str(i)]:
        diplotype = j + '/' + k
        locals() ['d' + str(i)] [diplotype] = locals() ['x' + str(i)][j] *  locals() ['x' + str(i)][k]

  #print (d1)
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

  for i in race:
    locals()[str(i)] = dic_race2dtfq[i]
  
  for diplotype in American.keys():
    f = []
    for r in race:
      f.append(str((locals()[str(r)][diplotype])))
    print(diplotype + "\t" + "\t".join(f), file = dt_fq)
  
  dt_fq.close()


