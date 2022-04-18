import pandas as pd
import numpy as np
import re, itertools


##################
##### Assets #####
##################
dic_race2col = {
  "African_American_Afro_Caribbean": 0,
  "American": 1,
  "Central_South_Asian": 2,
  "East_Asian": 3,
  "European": 4,
  "Latino": 5,
  "Near_Eastern": 6,
  "Oceanian": 7,
  "Sub_Saharan_African": 8
}

# Degenerate base means more than one base possibility at a particular position, 
# this is usually the case when a DNA sequence is derived from amino acid sequence with codon based sequence.
dic_degenerate_bases = {
  'R': ['A', 'G'],
  'Y': ['C', 'T'],
  'M': ['A', 'C'],
  'K': ['G', 'T'],
  'S': ['G', 'C'],
  'W': ['A', 'T'],
  'H': ['A', 'T', 'C'],
  'B': ['G', 'T', 'C'],
  'V': ['G', 'A', 'C'],
  'D': ['G', 'A', 'T'],
  'N': ['A', 'T', 'C', 'G']
}

dic_gene_ref_haplotype = {
  'ABCG2': 'G',
  'CACNA1S': 'Reference',
  'CYP2B6': '*1',
  'CYP2C19': '*38',
  'CYP2C9': '*1',
  'CYP2C8': '*1',
  'CYP2D6': '*1',
  'CYP3A4': '*1',
  'CYP3A5': '*1',
  'CYP4F2': '*1',
  'CFTR': 'ivacaftor non-responsive CFTR sequence',
  'DPYD': 'Reference',
  'MT-RNR1': 'Reference',
  'G6PD': 'B (wildtype)',
  'NUDT15': '*1',
  'RYR1': 'Reference',
  'SLCO1B1': '*1',
  'TPMT': '*1',
  'UGT1A1': '*1',
  'VKORC1': 'C',
  'IFNL3': 'C'
}


######################
##### Funcitions #####
######################
def phase_genotype(gt):
  # Sometimes, there are more than 0 and 1 codes
  ref_count = len(re.findall('0', gt))
  if ref_count == 2:
    genotype = 0
  elif ref_count == 1:
    genotype = 1
  elif ref_count == 0:
    genotype = 2
  return (genotype)


def index2alt_gt(pos, var, gt, non_gt):
  pos_index = pos.index(var)
  non_gt[pos_index] = gt
  return (non_gt)


def phase_alle(allele_definition_table, gene):
  pos = []
  list_alt = []
  list_ref = []
  index_nonsnp = []
  dic_alle2genotype = {}
  with open(allele_definition_table, 'r', encoding='utf-8') as file:
    flag = 0
    for line in file:
      line = line.replace("\n", "")
      info = line.split('\t')
      flag_judge = info[0]
      if 'GRCh38' in line:
        for i in info:
          matchobj = re.search(r'\w\.(\d+)(\w*)', i)
          # The traditional variant will match the following pattern.
          #matchobj = re.search(r'\w\.(\d+)(\w)\>(\w)', i)
          if matchobj:
            pos.append(matchobj.group(1))
            list_alt.append(matchobj.group(2))
          else:
            matchobj = re.search(r'\w\.(\d+)(\_\d+)?(del|ins)(\w*)', i)
            #matchobj = re.search(r'g\.(\d+)del(\w*)', i)
            if matchobj:
              pos.append(matchobj.group(1))
              list_alt.append(matchobj.group(4))
            elif "GRCh38" not in i:
              index_nonsnp.append(info.index(i)-1)
              pos.append(info.index(i)-1)
              list_alt.append(info.index(i)-1)
          print(pos, list_alt)
      
      else:
        if "rsID" in line:
          for i in index_nonsnp:
            pos[i] = info[i+1]
            list_alt[i] = info[i+1]
      
      # Start to process the haplotypes
      if flag == 1:
        alle = info.pop(0)
        
        if alle == dic_gene_ref_haplotype[gene]:
          list_ref = info
          gt = [0 for genotype in range(len(info))]
        else:
          gt = []
          for j in range(0, len(info)):
            if info[j] == '' or info[j] == list_ref[j]:
              genotype = 0
            elif info[j] in list(dic_degenerate_bases.keys()):
              genotype = []
              for wobble in dic_degenerate_bases[info[j]]:
                if wobble == list_ref[j]:
                  genotype.append(0)
                else:
                  genotype.append(1)
            else:
              genotype = 1
            gt.append(genotype)
        dic_alle2genotype[alle] = gt
      
      if re.search('%s Allele' % gene, flag_judge, re.IGNORECASE):
        flag = 1
      
  return (pos, list_alt, dic_alle2genotype)


def read_vcf(vcf_file, pos):
  dic_rs2gt = {}
  gene_gt = [0 if i else i for i in pos]
  with open(vcf_file, "r", encoding="utf-8") as file:
    for line in file:
      if line[1] != '#':
        line = line.replace("\n", "")
        info = line.split('\t')
        if "CHROM" in line:
          format_index = info.index("FORMAT")
          continue
        
        alt_pos, rsid = info[1], info[2]
        format = info[format_index].split(":")
        gt_index = format.index("GT")
        gt = info[-1].split(":")[gt_index]
        
        genotype = phase_genotype(gt)
        # Replace the genotype of gene_gt
        if alt_pos in pos:
          gene_gt = index2alt_gt(pos, alt_pos, genotype, gene_gt)
        elif rsid in pos:
          gene_gt = index2alt_gt(pos, rsid, genotype, gene_gt)
        
        if info[2] != ".":
          if genotype == 2:
            var = info[4] + info[4]
          elif genotype == 1:
            var = info[3] + info[4]
          else:
            var = info[3] + info[3]
          dic_rs2gt[info[2]] = var
    
    return(gene_gt, dic_rs2gt)


def import_fre(fre_file):
  dic_diplotype2fre = {}
  with open(fre_file, 'r', encoding='utf-8') as file:
    for line in file:
      if "alle" not in line:
        line = line.replace('\n', '')
        info = line.split('\t')
        dip = info.pop(0)
        dic_diplotype2fre[dip] = info
  return (dic_diplotype2fre)


def less_than_0(x):
  if x < 0:
    return x


def mix_rank(candidate, dic_diplotype2fre, race):
  # Rank the diplotype frequency
  dip_fre = []
  for key in candidate.keys():
    dip_fre.append(float(dic_diplotype2fre[key][dic_race2col[race]]))
  dip_fre_sorted = sorted(set(dip_fre))
  # Rank the difference of genotype
  cut_gt_sorted = sorted(set(candidate.values()), reverse=True)
  # Population rank (bigger is better) - Genotype rank (smaller is better)
  dic_diplotype2score = {}
  for key in candidate.keys():
    fre = float(dic_diplotype2fre[key][dic_race2col[race]])
    # Return index(rank)
    dic_diplotype2score[key] = dip_fre_sorted.index(fre) + cut_gt_sorted.index(candidate[key])
  
  return(dic_diplotype2score)


def phase_diplotype(dic_alle2genotype, gene_genotype, dic_diplotype2fre, race):
  candidate = {}
  combination = list(itertools.combinations_with_replacement(sorted(dic_alle2genotype), 2))
  loci_count = len(list(dic_alle2genotype.values())[0])
  for comb in combination:
    i = comb[0]; j = comb[1]
    cut_gt = []
    for index in range(0, loci_count):
      hap1 = dic_alle2genotype[i][index]
      hap2 = dic_alle2genotype[j][index]
      if type(hap1) is int:
        hap1 = [int(h) for h in str(hap1)]
      if type(hap2) is int:
        hap2 = [int(h) for h in str(hap2)]
      # cut_gt_tmp
      cut_gt_tmp = []
      for item in itertools.product(hap1, hap2):
        cut_gt_tmp.append(float(gene_genotype[index]) - sum(item))
      # cut_gt
      if 0 in cut_gt_tmp:
        cut_gt.append(0)
      elif 1 in cut_gt_tmp:
        cut_gt.append(1)
      else:
        cut_gt.append(max(cut_gt_tmp))
    print(cut_gt)
    if len(list(filter(less_than_0, cut_gt))):
      pass
    else:
      candidate[(i + '/' + j)] = sum(cut_gt)
  
  min_score = candidate[min(candidate, key=candidate.get)]
  # if min score is 0 means the exact map
  if min_score == 0:
    sub_candidate = {}
    for i in candidate:
      if candidate[i] == min_score:
        sub_candidate[i] = min_score
    # diplotype Score
    dic_diplotype2score = mix_rank(sub_candidate, dic_diplotype2fre, race)
  else:
    dic_diplotype2score = mix_rank(candidate, dic_diplotype2fre, race)
  
  # Sum result has equal rank
  diplotype = []
  max_rank = max(dic_diplotype2score.values())
  for key in dic_diplotype2score.keys():
    if dic_diplotype2score[key] == max_rank:
      # Store in the form of tuple
      diplotype.append(tuple(key.split('/')))
  return (diplotype)


def predict(vcf_file, race, gene_list):
  dic_gene2diplotype = {}
  for gene in gene_list:
    allele_definition_table = "./assets/definition/%s_allele_definition_table.txt" % (gene)
    fre_table = "./assets/diplotype_frequency/%s_diplotype_frequency_table.txt" % (gene)
    
    pos, list_alt, dic_alle2genotype = phase_alle(allele_definition_table, gene)
    dic_diplotype2fre = import_fre(fre_table)
    gene_genotype, dic_rs2gt = read_vcf(vcf_file, pos)
    diplotype = phase_diplotype(dic_alle2genotype, gene_genotype, dic_diplotype2fre, race)
    # Save the result in the format of list(tuple1, tuple2, ...)
    dic_gene2diplotype[gene] = diplotype
  
  return(dic_gene2diplotype)