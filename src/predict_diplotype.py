import pandas as pd
import numpy as np
import re, itertools, json, os


def parse_input_allele(vcf_file, info):
  hap_define = info['haplotype_definition']
  ref_hap = info['reference_haplotype']
  
  ## Read input VCF
  with open(vcf_file, "r", encoding="utf-8") as file:
    lines = file.readlines()
  
  lines_info = []
  for line in lines:
    if line.startswith('##') is False:
      lines_info.append(line.strip().split('\t'))
  
  vcf_df = pd.DataFrame(lines_info)
  vcf_df.columns = vcf_df.iloc[0]
  vcf_df = vcf_df.iloc[1:,]
  vcf_df['#CHROM'] = vcf_df['#CHROM'].map(lambda x: re.sub('chr|Chr|CHR', '', x)).astype('str')
  vcf_df = vcf_df[vcf_df['#CHROM'] == info['chrom']]
  vcf_df['POS'] = vcf_df['POS'].astype('int')
  cols = vcf_df.columns.to_list()
  
  vcf_alleles = {}
  hap_pos = list(hap_define[ref_hap].keys())
  for source_pos in hap_pos:
    # ref_hap_base, Only two loci of CYP2D6 gene will have more than one ref_hap_base
    ref_hap_base = hap_define[ref_hap][source_pos]
    
    # Get all possible bases at this position
    defined = []
    for key in hap_define.keys():
      defined.extend(hap_define[key][source_pos])
    defined = list(set(defined))
    
    # Transfer the positions into the format of list
    pos_rs = source_pos.split(':')
    pos = []
    tmp = pos_rs[0].split('-')
    if len(tmp) > 1:
      for p in range(int(tmp[0]), int(tmp[1])+1):
        pos.append(p)
    else:
      pos = [int(pos_rs[0])]
    
    ###### Start to parse the input vcf
    # Filter by pos and rs_id
    mat = vcf_df[(vcf_df['POS'].isin(pos)) | (vcf_df['ID'] == pos_rs[1])]
    if mat.empty:
      vcf_alleles[source_pos] = (ref_hap_base, ref_hap_base)
    else:
      is_wild_type = 1
      for index, row in mat.iterrows():
        cols = mat.columns.to_list()
        format = row[cols.index("FORMAT")].split(":")
        gt = row[-1].split(":")[format.index("GT")]
        if re.findall('0', gt) == ['0', '0']:
          continue
        else:
          ## Only process the first line which genotype is not wild type
          ## Therefore, the end of 'else' is break
          tuple_res = ()
          is_wild_type = 0
          ref = row[cols.index("REF")]
          alts = row[cols.index("ALT")].split(",")
          opts = [ref]; opts.extend(alts)
          gts = re.split('/|\|', gt)
          if len(gts) == 1: # chrX
            #print('GT info is not standard!'); print(gts, row.POS)
            gts.append(gts[0])
          for gt_index in gts:
            base = None
            if gt_index == '0':
              base = ref_hap_base
            else:
              alt = opts[int(gt_index)]
              # SNP #
              if len(alt) == len(ref):
                base = alt
              ### ! If alt was not match the definition, try to combine with the following position base and test again
              # Del #
              elif len(alt) < len(ref):
                base = 'del%s' % ref[len(alt):]
                if base not in defined:
                  pos_of_mat = mat.index.to_list()
                  pos_order = pos_of_mat.index(index)
                  if pos_order < len(pos_of_mat)-1 and mat.loc[pos_of_mat[pos_order+1], 'POS'] == row.POS+1:
                    new_row = mat.loc[pos_of_mat[pos_order+1],]
                    base = 'del%s' % ref[len(alt)+1:] + new_row.REF
                  elif int(row.POS) == 42128173 and base == 'delCTT':
                    base = 'delTCT'
                ## A smooth judge part
                if base not in defined:
                  # small indel
                  if ref_hap_base[0].startswith('ref') and base[-1] == ref_hap_base[0][-2]:
                    base = 'del' + base[4:] + ref_hap_base[0][-1]
                  # long del
                  else:
                    modd = ref_hap_base[0].replace('ref', '')
                    if re.search(modd, base):
                      matched_span = re.search(modd, base).span()
                      if modd.startswith(base[matched_span[1]:]):
                        base = 'del' + modd * int((matched_span[1] - matched_span[0])/len(modd) + 1)
                      else:
                        base = 'del' + modd * int((matched_span[1] - matched_span[0])/len(modd))
                if base not in defined:
                  print('Warning in Del!'); print(row.T); print(pos); print(base)
              # Ins #
              elif len(alt) > len(ref):
                # if alt[0:len(ref)] == ref:
                base = 'ins%s' % alt[len(ref):]
                if base not in defined:
                  pos_of_mat = mat.index.to_list()
                  pos_order = pos_of_mat.index(index)
                  if pos_order < len(pos_of_mat)-1 and mat.loc[pos_of_mat[pos_order+1], 'POS'] == row.POS+1:
                    new_row = mat.loc[pos_of_mat[pos_order+1],]
                    base = 'ins%s' % alt[len(ref)+1:] + new_row.REF
                ## A smooth judge part
                if base not in defined:
                  # small indel
                  if ref_hap_base[0].startswith('ref') and base[-1] == ref_hap_base[0][-2]:
                    base = 'ins' + base[4:] + ref_hap_base[0][-1]
                  # long dup
                  else:
                    flag = 0
                    if ref_hap_base[0].startswith('ref') is False:
                      for r in defined:
                        if r.startswith('delins'):
                          modd = 'C'; flag = 1; break
                    else:
                      modd = ref_hap_base[0].replace('ref', '')
                    print(modd, base)
                    if re.search(modd, base):
                      matched_span = re.search(modd, base).span()
                      if modd.startswith(base[matched_span[1]:]):
                        base = 'ins' + modd * int((matched_span[1] - matched_span[0])/len(modd) + 1)
                      else:
                        base = 'ins' + modd * int((matched_span[1] - matched_span[0])/len(modd))
                    if flag == 1:
                      base = 'del' + base + 'C'
                if base not in defined:
                  print(base)
                  print('Warning in Ins or Dup!'); print(row.T); print(pos)
                # else:
                #   print('Warning in Ins or Dup!'); print(row.T); print(pos)
            ## Add the result into tuple_res
            tuple_res = tuple_res + (base,)
          break
      
      if is_wild_type == 1:
        vcf_alleles[source_pos] = (ref_hap_base, ref_hap_base)
      else:
        vcf_alleles[source_pos] = tuple_res
        if None in tuple_res:
          print('Warning! None was reported.')
  
  # Check the format of output
  fine_vcf_alleles = {}
  for source_pos in vcf_alleles.keys():
    tuple_res = ()
    for base in vcf_alleles[source_pos]:
      if type(base) != list:
        tuple_res = tuple_res + ([base],)
      else:
        tuple_res = tuple_res + (base,)
    fine_vcf_alleles[source_pos] = tuple_res
  
  return(fine_vcf_alleles)


def predict_diplotype(vcf_alleles, info, race):
  hap_define = info['haplotype_definition']
  hap_mutated_loci = info['haplotype_mutated_loci']
  diplotype_candidates = list(itertools.combinations_with_replacement(hap_define.keys(), 2))
  ### 1st Ranking by haplotype definition
  candidate = {}
  for dip in diplotype_candidates:
    hap1 = hap_define[dip[0]]
    hap2 = hap_define[dip[1]]
    ## Calculate difference
    difference_step1 = dict(zip(vcf_alleles.keys(), [np.nan] * len(hap1)))
    for pos_rs in vcf_alleles.keys():
      defined_alleles = list(itertools.product(hap1[pos_rs], hap2[pos_rs]))
      tuple_res = vcf_alleles[pos_rs]
      alleles = list(itertools.product(tuple_res[0], tuple_res[1]))
      # Must sort the results
      defined_alleles = [sorted(ele) for ele in defined_alleles]
      alleles = [sorted(ele) for ele in alleles]
      # Extact matched
      score1 = []; score2 = []
      for allele in alleles:
        if allele in defined_alleles:
          #difference_step1[pos_rs] = 0
          score1 = [0]; score2 = [0]
          break
        else:
          diff = []
          ### Part 1: Input vcf can have more variants than definition but less is better.
          for define in defined_alleles:
            in_vcf_notin_defined = len(set(allele).difference(set(define)))
            if in_vcf_notin_defined == 1:
              score1.append(-1)
            else:
              score1.append(-2)
          ### Part 2: Haplotype muated loci. We need to make a case-by-case judgment.
          # When only one haplotype is involved in a locus, we do not need to ensure that both strands meet the requirements.
          allele_cp = allele
          for hap in dip:
            score2_hap = 0
            # It is possible to delete both items of allele for degenerate bases.
            if len(allele_cp) == 0:
              allele_cp = allele
            if pos_rs in hap_mutated_loci[hap]:
              star_define = hap_define[hap][pos_rs]
              star_matched = set(allele_cp).intersection(set(star_define))
              for m in star_matched:
                allele_cp.remove(m)
              if len(allele_cp) == 2:
                score2_hap = -99
            score2.append(score2_hap)
      # Add the results of score1 and score2
      difference_step1[pos_rs] = max(score1) + min(score2)
      # print(pos_rs, score1, score2)
    if min(difference_step1.values()) > -99:
    # if '*80' in dip:
    #   pass
    # else:
      candidate['%s/%s' % (dip[0], dip[1])] = difference_step1
  
  # print(candidate.keys())
  # Whether is 0
  extract_match_res = []
  for dip in candidate.keys():
    # print(dip, candidate[dip].values())
    if set(candidate[dip].values()) == {0}:
      extract_match_res.append(dip)
  
  # 1. select diplotypes which 1st rank is max
  rank_step1 = {}
  for dip in candidate.keys():
    rank_step1[dip] = sum(candidate[dip].values())
  
  uniq_diff = sorted(set(rank_step1.values()))
  rank_step1_res = [k for k,v in rank_step1.items() if v == max(uniq_diff)]#; rank_step1_res
  # 2. only ranked the above diplotypes with population frequency
  dip_freq = info['diplotype_frequency']
  rank_step2 = {}
  for dip in rank_step1_res:
    rank_step2[dip] = dip_freq[dip][race]
  uniq_freq = sorted(set(rank_step2.values()))
  final_rank_res = [k for k,v in rank_step2.items() if v == max(uniq_freq)]#; final_rank_res
  return("; ".join(extract_match_res), "; ".join(rank_step1_res), "; ".join(final_rank_res))


def main(vcf_file, race, gene_list):
  cpat_dip_base = json.loads(open("./assets/cpat_dip_base.json").read())
  dic_gene2diplotype = {}
  for gene in gene_list:
    # print(gene)
    info = cpat_dip_base[gene]
    vcf_alleles = parse_input_allele(vcf_file, info)
    extract_match_res, rank_step1_res, final_rank_res = predict_diplotype(vcf_alleles, info, race)
    dic_gene2diplotype[gene] = {'extract_match_res': extract_match_res, 'rank_step1_res': rank_step1_res, 'final_rank_res': final_rank_res}
    #dic_gene2diplotype[gene] = final_rank_res
  
  return(dic_gene2diplotype)






########################################################################
#######################################################################
cpat_dip_base = json.loads(open("./assets/cpat_dip_base.json").read())
gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",\
             "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15", "IFNL3", \
             "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]
########### Test the accuracy
pharmgkb_100genome = {'AFR': 'African American/Afro-Caribbean', 'AMR': 'Latino', 'EAS': 'East Asian', 'EUR': 'European'}
metadata = pd.read_csv('../cpat_data/test/population.txt', sep='\t')
test_vcfs = os.listdir('../cpat_data/test/88samples/')
extract_match_res = {}
rank_step1_res = {}
final_rank_res = {}
n = 0
# test_vcfs = ['NA18945.cpat.vcf']
gene_list = ['CYP2D6']
for vcf in test_vcfs:
  if vcf.endswith('.vcf'):
    n = n+1
    print("%i/88" % n)
    vcf_file = '../cpat_data/test/88samples/%s' % vcf
    # Get sample and race info
    sample = vcf.split('.')[0]
    race = metadata.loc[metadata['Get-RM 137 Samples']==sample, 'Superpopulation'].to_list()[0]
    race = pharmgkb_100genome[race]
    print(sample, race)
    res = main(vcf_file, race, gene_list)
    print(res)
    # Extract match res
    extract_match_res[sample] = {}
    for gene in gene_list:
      extract_match_res[sample][gene] = res[gene]['extract_match_res']
    # Step1 rank res
    rank_step1_res[sample] = {}
    for gene in gene_list:
      rank_step1_res[sample][gene] = res[gene]['rank_step1_res']
    # Final rank res
    final_rank_res[sample] = {}
    for gene in gene_list:
      final_rank_res[sample][gene] = res[gene]['final_rank_res']


extract_df = pd.DataFrame(extract_match_res).T
step1_df = pd.DataFrame(rank_step1_res).T
final_df = pd.DataFrame(final_rank_res).T

## Merge the above dfs into one df
extract_df.columns = ['%s.raw' % i for i in extract_df.columns]
step1_df.columns = ['%s.step1' % i for i in step1_df.columns]
final_df.columns = ['%s.step2' % i for i in final_df.columns]
df = pd.concat([extract_df, step1_df, final_df], axis=1)
df[sorted(df.columns)].to_csv('../manuscript/diplotype_concordance/detail_merge_v20220422_cyp2d6.txt', sep='\t')


