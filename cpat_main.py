#!/usr/bin/python
# -*- coding: UTF-8 -*-

from src import genotype_resolution, clinical_annotation#, pgx_report
import getopt, sys
# import pandas as pd
# import numpy as np
# from pybedtools import BedTool


def main(argv):
  race = ''
  germline_vcf = ''
  outdir = ''
  
  version = '0.1.0'
  help = '''
  Usage: python pgi.py [OPTIONS]
  Personal multi-omics data, including somatic variants, germline variants, gene expression 
  status, copy number variants, msi and tmb are combined for drug recommendations.
  Options:
    -r, --race                      Race of the input vcf. There are 9 selections, i.e., African American/Afro-Caribbean, American, 
                                    Central/South Asian, East Asian, European, Latino, Near Eastern, Oceanian, Sub-Saharan African.
    -i, --germline_vcf              Unannotated vcf file, preferably germline variant.
    -o, --outdir TEXT               Create report in the specified output path.
    --version                       Show the version and exit.
    -h, --help                      Show this message and exit.
  '''

  try:
    opts, args = getopt.getopt(
      argv, "hvr:i:o", ["help", "version", "race=", "germline_vcf=", "outdir="])
  except getopt.GetoptError:
    print(help)
  for opt, arg in opts:
    if opt == '-h':
      print(help)
      sys.exit()
    elif opt in ("-v", "--version"):
      print(version)
      sys.exit()
    elif opt in ("-r", "--race"):
      race = arg
    elif opt in ("-i", "--germline_vcf"):
      germline_vcf = arg
    elif opt in ("-o", "--output"):
      outdir = arg
  
  germline_vcf = './test/88samples/HG00276.pgx.vcf'
  race = 'European'
  outdir = './test/'
  dic_diplotype, dic_rs2gt, hla_subtypes = genotype_resolution.resolution(race, germline_vcf, outdir)
  pgx_summary, clinical_anno_table = clinical_annotation.annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
  pgx_report.report(race, pgx_summary, clinical_anno_table, outdir)