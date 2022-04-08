#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pandas as pd
import numpy as np
import os, re, sqlite3, time

"""
There three columns of the input dataframe:
1. Gene Symbol
2. Variant (rsID) or Haplotype
3. Genotype or Diplotype
"""

def report(race, pgx_summary, clinical_anno_table, outdir):
  with open(outdir, 'w+') as f:
    ## Style
    style="""
    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <title>CPAT Report</title>
    <link rel="stylesheet" href="./css/darkdown.css">
    </head>
    <body>
    """
    print(style, file=f)

    ## Part 0: Basic information
    time_info = time.asctime(time.localtime(time.time()))
    print('<h1>CPAT Report</h1>\n<blockquote>', file=f)
    print('<p>Report Time: %s</p>' % time_info, file=f)
    print('<p>Biogeographic Group: %s</p>' % race, file=f)
    
    ## Part 1: Diplotype
    # time_info = time.asctime(time.localtime(time.time()))
    # print('<h1>CPAT Report</h1>\n<blockquote>', file=f)
    # print('<p>Report Time: %s</p>' % time_info, file=f)
    # print('<p>Biogeographic Group: %s</p>' % race, file=f)
    
    ## Part 2: Pharmacogenomics Annotation
    part2_header = """<h2>Pharmacogenomics Summary</h2>\n<p>CPAT calculated the effect scores of the drugs of interest in the five phenotype categories based on the genotypes obtained from the analysis. Evidence levels 1A, 1B, 2, 3, and 4 are used to indicate the degree of influence of a drug on a particular drug response.</p>\n</blockquote>
    """
    print(part2_header, file=f)
    
    table_html = """
    <table border="1" cellspacing="0">
      <tr><th>Evidence</th><th>Drug</th></tr>
      <tr><td bgcolor="#55B979" width="80px"><font color="white"><b>Level 1A</b></font></td><td>%s</td></tr>
      <tr><td bgcolor="#55B979" width="80px"><font color="white"><b>Level 1B</b></font></td><td>%s</td></tr>
      <tr><td bgcolor="#3F72D8" width="80px"><font color="white"><b>Level 2</b></font></td><td>%s</td></tr>
      <tr><td bgcolor="#F5C344" width="80px"><font color="white"><b>Level 3</b></font></td><td>%s</td></tr>
      <tr><td bgcolor="#B64641" width="80px"><font color="white"><b>Level 4</b></font></td><td>%s</td></tr>
    </table>"""
    
    for index, col in pgx_summary.iteritems():
      print('<h3>%s</h3>' % index, file=f)
      print(table_html%(col[0], col[1], col[2], col[3], col[4]), file=f)
    
    
    ## Part 3: Genotype
    print('<h2>Pharmacogenomics Details</h2>', file=f)
    print('<h3>Diplotypes predicted by CPAT</h3>', file=f)
    detail1 = clinical_anno_table[clinical_anno_table.Class == 'Diplotype'].drop(columns=['Class']).reset_index(drop = True)
    header = '<table border="1" cellspacing="0">\n<tr><th>ID</th><th>Gene</th><th>Variant</th><th>Drug</th><th>Phenotypes</th><th>Evidence</th><th>Alleles</th><th>Category</th><th>Function</th></tr>'
    for index, row in detail1.iterrows():
      header = header + '\n<tr><td><a href=%s>%s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (row[9], index, row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[8])
    header = header + '\n</table>'
    print(header, file=f)

    print('<h3>Genotypes called by VCF</h3>', file=f)
    detail2 = clinical_anno_table[clinical_anno_table.Class == 'Single'].drop(columns=['Class']).reset_index(drop = True)
    header = '<table border="1" cellspacing="0">\n<tr><th>ID</th><th>Gene</th><th>Variant</th><th>Drug</th><th>Phenotypes</th><th>Evidence</th><th>Alleles</th><th>Category</th><th>Function</th></tr>'
    for index, row in detail2.iterrows():
      header = header + '\n<tr><td><a href=%s>%s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (row[9], index, row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[8])
    header = header + '\n</table>'
    print(header, file=f)
    
    # Part 4: Disclaimers
    disclaimer = """
    <h2>Disclaimers</h2>
    <p>The report incorporates analyses of peer-reviewed studies and other publicly available information identified by CPAT by State Key Laboratory of Genetic Engineering from the School of Life Sciences and Human Phenome Institute, Fudan University, Shanghai, China. These analyses and information may include associations between a molecular alteration (or lack of alteration) and one or more drugs with potential clinical benefit (or potential lack of clinical benefit), including drug candidates that are being studied in clinical research.</p>
    <p>Note: A finding of biomarker alteration does not necessarily indicate pharmacologic effectiveness (or lack thereof) of any drug or treatment regimen; a finding of no biomarker alteration does not necessarily indicate lack of pharmacologic effectiveness (or effectiveness) of any drug or treatment regimen.</p>
    <p>No Guarantee of Clinical Benefit: This Report makes no promises or guarantees that a particular drug will be effective in the treatment of disease in any patient. This report also makes no promises or guarantees that a drug with a potential lack of clinical benefit will provide no clinical benefit.</p>
    <p>Treatment Decisions are Responsibility of Physician: Drugs referenced in this report may not be suitable for a particular patient. The selection of any, all, or none of the drugs associated with potential clinical benefit (or potential lack of clinical benefit) resides entirely within the discretion of the treating physician. Indeed, the information in this report must be considered in conjunction with all other relevant information regarding a particular patient, before the patient's treating physician recommends a course of treatment. Decisions on patient care and treatment must be based on the independent medical judgment of the treating physician, taking into consideration all applicable information concerning the patient's condition, such as patient and family history, physical examinations, information from other diagnostic tests, and patient preferences, following the standard of care in a given community. A treating physician's decisions should not be based on a single test, such as this test or the information contained in this report.</p>
    <p>When using results obtained from CPAT, you agree to cite CPAT.</p>
    \n</body>\n</html>"""
    print(disclaimer, file=f)

