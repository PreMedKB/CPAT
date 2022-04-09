#!/usr/bin/python
# -*- coding: UTF-8 -*-

import time

"""
There three columns of the input dataframe:
1. Gene Symbol
2. Variant (rsID) or Haplotype
3. Genotype or Diplotype
"""

def report(race, pgx_summary, clinical_anno_table, outdir, basename):
  fp = "%s/%s.cpat.html" % (outdir, basename)
  with open(fp, 'w+') as f:
    ## Style
    style="""
    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <title>CPAT Report</title>
    <script src="https://kit.fontawesome.com/e540049a97.js" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://raw.githack.com/lyaqing/CPAT/main/html/css/custom.css">
    <ul>
      <p></p>
      <li><a href="#summary"><b>PGx Summary</b></a></li>
      <li><a href="#toxicity">&nbsp;&nbsp;Toxicity</a></li>
      <li><a href="#dosage">&nbsp;&nbsp;Dosage</a></li>
      <li><a href="#efficacy">&nbsp;&nbsp;Efficacy</a></li>
      <li><a href="#metabolism/pk">&nbsp;&nbsp;Metabolism/PK</a></li>
      <li><a href="#other">&nbsp;&nbsp;Other</a></li>
      <li><a href="#guideline"><b>Dosing Guideline</b></a></li>
      <li><a href="#detail"><b>Genotype Detail</b></a></li>
      <li><a href="#haplotype/pk">&nbsp;&nbsp;Haplotype/PK</a></li>
      <li><a href="#snp/indel">&nbsp;&nbsp;SNP/Indel</a></li>
      <li><a href="#about"><b>About</b></a></li>
    </ul>
    <div style="margin-left:15rem;margin-right:0rem;padding:1px 16px;height:1000px;">
    """
    print(style, file=f)

    ## Part 0: Basic information
    basic_info = """
  <h1>CPAT Report</h1>
  <blockquote>
    <p><b>Report Time:</b> %s<br><b>Biogeographic Group:</b> %s</p>
  </blockquote>
    """
    print(basic_info%(time.asctime(time.localtime(time.time())), race), file=f)

    ## Part 1: Sort disclaimer
    disclaimer_short = """
  <div class="alert alert-info">
    <em>Disclaimer:</em> CPAT reports are subject to iteration as the release version changes. In the current release you should only use this software to assess whether the CPAT executable will compile and run properly on your system. CPAT can only generate recommendations based on information from the imported software, so all information in the reports section is interpreted directly from the uploaded vcf file. Users recognise that they are using CPAT at their own risk.
  </div>
    """
    print(disclaimer_short, file=f)

    ## Part 2: Pharmacogenomics Annotation
    part2_header = """
  <h2 id="summary">Pharmacogenomics Summary</h2>
  <p>CPAT calculated the effect scores of the drugs of interest in the five phenotype categories based on the genotypes obtained from the analysis. Evidence levels 1A, 1B, 2, 3, and 4 are used to indicate the degree of influence of a drug on a particular drug response.</p>
    """
    print(part2_header, file=f)
    
    table_html = """
  <table id="customers" border="1" cellspacing="0">
    <tr>
    <th bgcolor="#ffffff"></th>
    <th><i class="fa-solid fa-square-caret-down"></i> Decreased</th>
    <th><i class="fa-solid fa-square-minus"></i> Moderate</th>
    <th><i class="fa-solid fa-square-caret-up"></i> Increased</th>
    </tr>
    <tr>
      <td bgcolor="#10304E" width="80px"><font color="white"><b>Level A</b></font></td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
    </tr>
    <tr>
      <td bgcolor="#4F7C58" width="80px"><font color="white"><b>Level B</b></font></td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
    </tr>
    <tr>
      <td bgcolor="#DAB156" width="80px"><font color="white"><b>Level C</b></font></td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
    </tr>
    <tr>
      <td bgcolor="#761621" width="80px"><font color="white"><b>Level D</b></font></td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
    </tr>
  </table>
  """
    
    categories = ['Toxicity', 'Dosage', 'Efficacy', 'Metabolism/PK', 'Other']
    fa = dict(zip(categories, ['fa-skull-crossbones', 'fa-pills', 'fa-vial-circle-check', 'fa-disease', 'fa-feather-pointed']))
    for cat in categories:
      i = 0
      html_input = []
      cat_pgx = pgx_summary[pgx_summary.PhenotypeCategory == cat]
      response = ['Decreased', 'Moderate', 'Increased']
      levels = ['A', 'B', 'C', 'D']
      for level in levels:
        for res in response:
          tmp = cat_pgx[(cat_pgx.Response == res) & (cat_pgx.EvidenceLevel == level)]
          if tmp.empty:
            html_input.append('')
          else:
            html_input.append('; '.join(tmp.index.to_list()))
          i = i+1
      # print
      print('<h3 id="%s"><i class="fa-solid %s"></i> %s</h3>' % (cat.lower(), fa[cat], cat), file=f)
      print(table_html%(html_input[0], html_input[1], html_input[2],
                        html_input[3], html_input[4], html_input[5],
                        html_input[6], html_input[7], html_input[8],
                        html_input[9], html_input[10], html_input[11]), file=f)
    
    ## Part 3: Dosing Guideline
    print('<h2 id="guideline">Dosing Guidelines</h2>', file=f)
    

    ## Part 4: Genotype Details
    print('<h2 id="detail">Genotype Details</h2>', file=f)
    print('<h3 id="haplotype">Haplotype/Diplotype predicted by CPAT</h3>', file=f)
    print('<p>This subsection provides CPAT predicted haplotypes based on the VCF calls, the haplotype definition table, and the haplotype population frequency table. PharmGKB and CPIC together summarized pharmacogenomic genes with explicit star(*) or named allele information are involved, specifically ABCG2, CACNA1S, CFTR, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP3A4, CYP3A5, CYP4F2, DPYD, G6PD, MT-RNR1, NUDT15, RYR1, SLCO1B1, TPMT, UGT1A1, VKORC1.</p>', file=f)
    detail1 = clinical_anno_table[clinical_anno_table.Class == 'Diplotype'].drop(columns=['Class']).reset_index(drop = True)
    header = '<table id="customers" border="1" cellspacing="0">\n<tr><th>ID</th><th>Gene</th><th>Variant</th><th>Alleles</th><th>Drug</th><th>Evidence</th><th>Category</th><th>Function</th></tr>'
    for index, row in detail1.iterrows():
      header = header + '\n<tr><td><a href=%s>%s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (row.URL, index, row.Gene, row.Variant, row.Alleles, row.Drug, row.EvidenceLevel, row.PhenotypeCategory, row.Function)
    header = header + '\n</table>'
    print(header, file=f)
    
    print('<h3 id="snp/indel">SNPs/Indels called by VCF</h3>', file=f)
    print('<p>This subsection depends on the VCF calls for the query positions provided. CPAT does not assume any reference calls for missing positions in the submitted VCF; all missing query positions are not considered in the allele determination process.</p>', file=f)
    detail2 = clinical_anno_table[clinical_anno_table.Class == 'Single'].drop(columns=['Class']).reset_index(drop = True)
    header = '<table id="customers" border="1" cellspacing="0">\n<tr><th>ID</th><th>Gene</th><th>Variant</th><th>Alleles</th><th>Drug</th><th>Evidence</th><th>Category</th><th>Function</th></tr>'
    for index, row in detail2.iterrows():
      header = header + '\n<tr><td><a href=%s>%s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (row.URL, index, row.Gene, row.Variant, row.Alleles, row.Drug, row.EvidenceLevel, row.PhenotypeCategory, row.Function)
    header = header + '\n</table>'
    print(header, file=f)
    
    # Part 4: Disclaimers
    disclaimer = """
  <h2 id="about">Disclaimers</h2>
  <p>The report incorporates analyses of peer-reviewed studies and other publicly available information identified by CPAT by State Key Laboratory of Genetic Engineering from the School of Life Sciences and Human Phenome Institute, Fudan University, Shanghai, China. These analyses and information may include associations between a molecular alteration (or lack of alteration) and one or more drugs with potential clinical benefit (or potential lack of clinical benefit), including drug candidates that are being studied in clinical research.<br>
  Note: A finding of biomarker alteration does not necessarily indicate pharmacologic effectiveness (or lack thereof) of any drug or treatment regimen; a finding of no biomarker alteration does not necessarily indicate lack of pharmacologic effectiveness (or effectiveness) of any drug or treatment regimen.<br>
  No Guarantee of Clinical Benefit: This Report makes no promises or guarantees that a particular drug will be effective in the treatment of disease in any patient. This report also makes no promises or guarantees that a drug with a potential lack of clinical benefit will provide no clinical benefit.<br>
  Treatment Decisions are Responsibility of Physician: Drugs referenced in this report may not be suitable for a particular patient. The selection of any, all, or none of the drugs associated with potential clinical benefit (or potential lack of clinical benefit) resides entirely within the discretion of the treating physician. Indeed, the information in this report must be considered in conjunction with all other relevant information regarding a particular patient, before the patient's treating physician recommends a course of treatment. Decisions on patient care and treatment must be based on the independent medical judgment of the treating physician, taking into consideration all applicable information concerning the patient's condition, such as patient and family history, physical examinations, information from other diagnostic tests, and patient preferences, following the standard of care in a given community. A treating physician's decisions should not be based on a single test, such as this test or the information contained in this report.<br>
  When using results obtained from CPAT, you agree to cite CPAT.
  </p>
  </div>
  </body>
  </html>"""
    print(disclaimer, file=f)

