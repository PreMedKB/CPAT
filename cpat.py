#!/usr/bin/python
# -*- coding: UTF-8 -*-

from src import genotype_resolution, clinical_annotation, pgx_report
import getopt, sys

def main(argv):
  version = '1.0.0'
  help = '''
  Usage: python cpat.py [OPTIONS]
  CPAT takes the germline variant calling format (VCF) file and population information as input and 
  outputs an HTML report of drug responses with prescription recommendations.
  Options:
    -s, --sample_id                 Sample ID.
    -r, --population                Biogeographic groups: African American/Afro-Caribbean, American, Central/South Asian,
                                    East Asian, European, Latino, Near Eastern, Oceanian, Sub-Saharan African.
    -i, --germline_vcf              Unannotated vcf file, preferably germline variant.
    -o, --outdir TEXT               Create report in the specified output path.
    -v, --version                   Show the version and exit.
    -h, --help                      Show this message and exit.
  '''
  
  try:
    opts, args = getopt.getopt(argv, "hvs:i:p:o", ["help", "version", "sample_id=", "germline_vcf=", "population=", "outdir="])
  except getopt.GetoptError:
    print(help)
  
  for opt, arg in opts:
    if opt in ("-h", "--help"):
      print(help)
      sys.exit()
    elif opt in ("-v", "--version"):
      print(version)
      sys.exit()
    elif opt in ("-r", "--population"):
      population = arg
    elif opt in ("-i", "--germline_vcf"):
      germline_vcf = arg
    elif opt in ("-o", "--output"):
      outdir = arg
    elif opt in ("-s", "--sample_id"):
      sample_id = arg
  
  try:
    print('  - Parsing PGx related genotypes ...')
    dic_diplotype, dic_rs2gt, hla_subtypes = genotype_resolution.resolution(population, germline_vcf)
    print('  - Annotating clinical information ...')
    pgx_summary, clinical_anno_table, dosing_guideline_table = clinical_annotation.annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
    print('  - Generating CPAT report ...')
    pgx_report.report(population, pgx_summary, dic_diplotype, clinical_anno_table, dosing_guideline_table, outdir, sample_id)
  except:
    print('Error appeared.')


if __name__ == "__main__":
  main(sys.argv[1:])