#!/usr/bin/python
# -*- coding: UTF-8 -*-

import time

"""
There three columns of the input dataframe:
1. Gene Symbol
2. Variant (rsID) or Haplotype
3. Genotype or Diplotype
"""

def report(race, pgx_summary, clinical_anno_table, outdir):
  #fp = "%stest_v20220407.html" % outdir
  fp = "./html/air_v20220408.html"
  with open(fp, 'w+') as f:
    ## Style
    style="""
    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <title>CPAT Report</title>
    <script src="https://kit.fontawesome.com/e540049a97.js" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="./css/custom.css">
    <style>
      #customers {
      font-family: Arial, Helvetica, sans-serif;
      border-collapse: collapse;
      width: 100%;
      }

      #customers td, #customers th {
        border: 1px solid #ddd;
        padding: 8px;
      }

      #customers tr:nth-child(even){background-color: #f2f2f2;}

      #customers tr:hover {background-color: #ddd;}

      #customers th {
        padding-top: 12px;
        padding-bottom: 12px;
        text-align: center;
        color: #44308d;
      }

      /* Full-height Fixed Vertical Navbar */
      body {
        margin: 0;
      }

      ul {
        list-style-type: none;
        margin: 0;
        padding: 0;
        width: 10%;
        background-color: #f1f1f1;
        position: fixed;
        height: 100%;
        overflow: auto;
      }

      li a {
        display: block;
        color: #000;
        padding: 8px 16px;
        text-decoration: none;
      }

      li a.active {
        background-color: #04AA6D;
        color: white;
      }

      li a:hover:not(.active) {
        background-color: #555;
        color: white;
      }
    </style>
    </head>
    <body>
    <ul>
      <li><a href="#summary">Summary</a></li>
        <ul>
          <li><a href="#toxicity">Toxicity</a></li>
          <li><a href="#dosage">Dosage</a></li>
          <li><a href="#efficacy">Efficacy</a></li>
          <li><a href="#metabolism/pk">Metabolism/PK</a></li>
          <li><a href="#other">Other</a></li>
        </ul>
      <li><a href="#detail">Detail</a></li>
      <li><a href="#about">About</a></li>
    </ul>

    <div style="margin-left:20%;padding:1px 16px;height:1000px;">
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
      <td bgcolor="#55B979" width="80px"><font color="white"><b>Level 1</b></font></td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
    </tr>
    <tr>
      <td bgcolor="#3F72D8" width="80px"><font color="white"><b>Level 2</b></font></td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
    </tr>
    <tr>
      <td bgcolor="#F5C344" width="80px"><font color="white"><b>Level 3</b></font></td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
      <td width="250px">%s</td>
    </tr>
  </table>
  """
    
    categories = ['Toxicity', 'Dosage', 'Efficacy', 'Metabolism/PK', 'Other']
    for cat in categories:
      i = 0
      html_input = []
      cat_pgx = pgx_summary[pgx_summary.PhenotypeCategory == cat]
      response = ['Decreased', 'Moderate', 'Increased']
      levels = ['1', '2', '3']
      for level in levels:
        for res in response:
          tmp = cat_pgx[(cat_pgx.Response == res) & (cat_pgx.EvidenceLevel == level)]
          if tmp.empty:
            html_input.append('')
          else:
            html_input.append('; '.join(tmp.index.to_list()))
          i = i+1
      # print
      print('<h3 id="%s">%s</h3>' % (cat.lower(), cat), file=f)
      print(table_html%(html_input[0], html_input[1], html_input[2],
                        html_input[3], html_input[4], html_input[5],
                        html_input[6], html_input[7], html_input[8]), file=f)
    
    ## Part 3: Genotype
    print('<h2 id="detail">Pharmacogenomics Details</h2>', file=f)
    print('<h3>Diplotypes predicted by CPAT</h3>', file=f)
    detail1 = clinical_anno_table[clinical_anno_table.Class == 'Diplotype'].drop(columns=['Class']).reset_index(drop = True)
    header = '<table id="customers" border="1" cellspacing="0">\n<tr><th>ID</th><th>Gene</th><th>Variant</th><th>Drug</th><th>Phenotypes</th><th>Evidence</th><th>Alleles</th><th>Category</th><th>Function</th></tr>'
    for index, row in detail1.iterrows():
      header = header + '\n<tr><td><a href=%s>%s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (row[9], index, row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[8])
    header = header + '\n</table>'
    print(header, file=f)

    print('<h3>Genotypes called by VCF</h3>', file=f)
    detail2 = clinical_anno_table[clinical_anno_table.Class == 'Single'].drop(columns=['Class']).reset_index(drop = True)
    header = '<table id="customers" border="1" cellspacing="0">\n<tr><th>ID</th><th>Gene</th><th>Variant</th><th>Drug</th><th>Phenotypes</th><th>Evidence</th><th>Alleles</th><th>Category</th><th>Function</th></tr>'
    for index, row in detail2.iterrows():
      header = header + '\n<tr><td><a href=%s>%s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (row[9], index, row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[8])
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

