#!/usr/bin/python
# -*- coding: UTF-8 -*-

from src import genotype_resolution, clinical_annotation, pgx_report
import getopt, sys

def main(argv):
  version = '1.0.0'
  help = '''
  Usage: python cpat.py [OPTIONS]
  CPAT takes the variant calling format (VCF) file and population information as input
  and outputs an HTML report of drug responses with prescription recommendations.
  Options:
    -s, --sample_id                 Sample ID.
    -i, --germline_vcf              Unannotated vcf file, preferably germline variant.
    -p, --population                The three-letter abbreviation for biogeographic groups:
                                    AAC (African American/Afro-Caribbean), AME (American),
                                    EAS (East Asian), EUR (European), LAT (Latino), NEA (Near Eastern),
                                    OCE (Oceanian), SAS (Central/South Asian), SSA (Sub-Saharan African).
    -o, --outdir TEXT               Create report in the specified output path.
    -v, --version                   Show the version and exit.
    -h, --help                      Show this message and exit.
  '''
  
  try:
    opts, args = getopt.getopt(argv, "hvs:i:p:o:", ["help", "version", "sample_id=", "germline_vcf=", "population=", "outdir="])
  except getopt.GetoptError:
    print(help)
  
  for opt, arg in opts:
    if opt in ("-h", "--help"):
      print(help)
      sys.exit()
    elif opt in ("-v", "--version"):
      print(version)
      sys.exit()
    elif opt in ("-s", "--sample_id"):
      sample_id = arg
    elif opt in ("-i", "--germline_vcf"):
      germline_vcf = arg
    elif opt in ("-p", "--population"):
      population = arg
    elif opt in ("-o", "--output"):
      outdir = arg
  
  try:
    population = population.upper()
    pop_dic = {'AAC': 'African American/Afro-Caribbean', 'AME': 'American', 'SAS': 'Central/South Asian', 'EAS': 'East Asian', 'EUR': 'European', 'LAT': 'Latino', 'NEA': 'Near Eastern', 'OCE': 'Oceanian', 'SSA': 'Sub-Saharan African'}
    if population not in pop_dic.keys():
      print('The input population is not included in CPAT. Please check if the abbreviation is used correctly.')
    
    print('  - Parsing PGx related genotypes ...')
    dic_diplotype, dic_rs2gt, hla_subtypes = genotype_resolution.resolution(pop_dic[population], germline_vcf)
    
    print('  - Annotating clinical information ...')
    pgx_summary, clinical_anno_table, dosing_guideline_table = clinical_annotation.annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
    
    print('  - Generating CPAT report ...')
    pgx_report.report(pop_dic[population], pgx_summary, dic_diplotype, clinical_anno_table, dosing_guideline_table, outdir, sample_id)
  except:
    print('Error appeared.')


if __name__ == "__main__":
  main(sys.argv[1:])