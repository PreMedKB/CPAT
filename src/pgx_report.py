#!/usr/bin/python
# -*- coding: UTF-8 -*-

import time

def report(race, pgx_summary, dic_diplotype, clinical_anno_table, dosing_guideline_table, outdir, sample_id):
  fp = "%s/%s.cpat.html" % (outdir, sample_id)
  with open(fp, 'w+') as f:
    ## Style
    style="""
    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <title>CPAT Report</title>
    <script src="https://kit.fontawesome.com/e540049a97.js" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://raw.githack.com/premedkb/cpat/main/assets/css/custom.css">
    <ul>
      <p></p>
      <li><a href="#home"><b>Home</b></a></li>
      <li><a href="#summary"><b>PGx Summary</b></a></li>
      <li><a href="#toxicity">&nbsp;&nbsp;Toxicity</a></li>
      <li><a href="#dosage">&nbsp;&nbsp;Dosage</a></li>
      <li><a href="#efficacy">&nbsp;&nbsp;Efficacy</a></li>
      <li><a href="#metabolism/pk">&nbsp;&nbsp;Metabolism/PK</a></li>
      <li><a href="#guideline"><b>Dosing Guideline</b></a></li>
      <li><a href="#diplotype"><b>Diplotype Detail</b></a></li>
      <li><a href="#annotation"><b>PGx Annotation</b></a></li>
      <li><a href="#about"><b>About</b></a></li>
    </ul>
    <div style="margin-left:15rem;margin-right:0rem;padding:1px 16px;height:1000px;">
    """
    print(style, file=f)

    ## Part 0: Basic information
    basic_info = """
    <h1 id="home">CPAT Report</h1>
    <blockquote>
      <p><b>Sample ID:</b> %s<br><b>Biogeographic Group:</b> %s<br><b>Report Time:</b> %s</p>
    </blockquote>
    """
    print(basic_info%(sample_id, race, time.asctime(time.localtime(time.time()))), file=f)
    
    ## Part 1: Sort disclaimer
    disclaimer_short = """
    <div class="alert alert-info-yellow">
      <em><b>Disclaimer:</b></em> The CPAT report iterates as the release version changes. In the current release, you should only use it to evaluate whether CPAT will compile and run properly on your system. All information in the CPAT report is interpreted directly from the uploaded VCF file. Users recognise that they are using CPAT at their own risk.
    </div>
    """
    print(disclaimer_short, file=f)

    ## Part 2: Pharmacogenomics Annotation
    part2_header = """
    <h2 id="summary">Pharmacogenomics Summary</h2>
    <p>CPAT combines genotypes with PharmGKB's clinical annotations to summarise a patient's response to a specific drug across four dimensions, including <font color="#000"><b>toxicity</b></font>, <font color="#000"><b>dosage</b></font>, <font color="#000"><b>efficacy</b></font>, <font color="#000"><b>metabolism/PK</b></font>. Response levels are indicated using <font color="#8E529A"><b>decreased</b></font>, <font color="#653F92"><b>moderate</b></font> and <font color="#44308d"><b>increased</b></font>.<br>CPAT labels results with variant-specific prescribing guidance in professional clinical guidelines or FDA-approved drug label annotations as <font color="#54A052"><b>Level A</b></font>, and otherwise as <font color="#3978B1"><b>Level B</b></font>.</p>
    """
    print(part2_header, file=f)
    
    table_html = """
    <table id="customers" border="1" cellspacing="0">
      <tr>
      <th bgcolor="#ffffff"></th>
      <th><font color="#8E529A"><i class="fa-solid fa-square-caret-down"></i> Decreased</font></th>
      <th><font color="#653F92"><i class="fa-solid fa-square-minus"></i> Moderate</font></th>
      <th><font color="#44308d"><i class="fa-solid fa-square-caret-up"></i> Increased</font></th>
      </tr>
      <tr>
        <td bgcolor="#54A052" width="80px"><font color="white"><b>Level A</b></font></td>
        <td width="250px">%s</td>
        <td width="250px">%s</td>
        <td width="250px">%s</td>
      </tr>
      <tr>
        <td bgcolor="#3978B1" width="80px"><font color="white"><b>Level B</b></font></td>
        <td width="250px">%s</td>
        <td width="250px">%s</td>
        <td width="250px">%s</td>
      </tr>
    </table>
    """
    
    categories = ['Toxicity', 'Dosage', 'Efficacy', 'Metabolism/PK']
    fa = dict(zip(categories, ['fa-skull-crossbones', 'fa-pills', 'fa-vial-circle-check', 'fa-disease', 'fa-feather-pointed']))
    for cat in categories:
      print('<h3 id="%s"><i class="fa-solid %s"></i> %s</h3>' % (cat.lower(), fa[cat], cat), file=f)
      i = 0
      html_input = []
      cat_pgx = pgx_summary[pgx_summary.PhenotypeCategory == cat]
      if cat_pgx.empty:
        print('<div class="alert alert-info-red">No clinical annotations matched for this category.</div>', file=f)
      else:
        response = ['Decreased', 'Moderate', 'Increased']
        levels = ['A', 'B']
        for level in levels:
          for res in response:
            tmp = cat_pgx[(cat_pgx.Response == res) & (cat_pgx.EvidenceLevel == level)]
            if tmp.empty:
              html_input.append('')
            else:
              html_input.append('; '.join(tmp.index.to_list()))
            i = i+1
        print(table_html%(html_input[0], html_input[1], html_input[2],
                          html_input[3], html_input[4], html_input[5]), file=f)
      
    ## Part 3: Dosing Guideline
    print('<h2 id="guideline">Dosing Guideline</h2>', file=f)
    print('<p>CPAT integrates brief annotations of genotype-based dosing recommendations after PharmGKB processing. Original PGx-based drug dosing guidelines include the <a href="http://cpicpgx.org/">Clinical Pharmacogenetics Implementation Consortium</a> (CPIC), the <a href="https://www.knmp.nl/dossiers/farmacogenetica/">Royal Dutch Association for the Advancement of Pharmacy - Pharmacogenetics Working Group</a> (DPWG), the <a href="https://cpnds.ubc.ca/">Canadian Pharmacogenomics Network for Drug Safety</a> (CPNDS), the French National Network for Pharmacogenetics (RNPGx), The Australian and New Zealand consensus guidelines (AusNZ), the Spanish Pharmacogenetics and Pharmacogenomics Society (SEFF), the Cystic Fibrosis Foundation (CFF), the American College of Rheumatology.</p>', file=f)
    # Drug - Detected variant or alleles - Dosing guidelines
    for drug in list(dosing_guideline_table.Drug.drop_duplicates()):
      print('<h3>%s</h3>' % drug, file=f)
      drug_sub = dosing_guideline_table[dosing_guideline_table.Drug == drug]
      for gene in list(drug_sub.Gene.drop_duplicates()):
        drug_by_gene = drug_sub[drug_sub.Gene == gene]
        # drug_by_gene.insert(drug_by_gene.shape[1], 'Report', drug_by_gene.Variant.str.cat(drug_by_gene.Alleles, sep=": "))
        print('<p><font color="#444">Related Gene: %s<br>Detected Alleles: %s</font></p>' % (gene, "; ".join(list(drug_by_gene.Report.drop_duplicates()))), file=f)
        drug_guide = drug_by_gene[['DosingURL', 'DosingSource', 'DosingAnnotation']].copy()
        drug_guide = drug_guide.drop_duplicates()
        for index, row in drug_guide.iterrows():
          if row.DosingAnnotation.startswith('There are currently no'):
            #print('<div class="alert alert-info-yellow"><a href=%s target="_blank" style="color: #8a6d3b"><i class="fa-solid fa-circle-info"></a></i><b> %s: </b>%s</div>' % (row.DosingURL, row.DosingSource, row.DosingAnnotation), file=f)
            print('<div class="alert alert-info-red"><a href=%s target="_blank" style="color: #7C3A37"><i class="fa-solid fa-circle-info"></a></i><b> %s: </b>%s</div>' % (row.DosingURL, row.DosingSource, row.DosingAnnotation), file=f)
          else:
            print('<div class="alert alert-info-blue"><a href=%s target="_blank"><i class="fa-solid fa-circle-info"></a></i><b> %s: </b>%s</div>' % (row.DosingURL, row.DosingSource, row.DosingAnnotation), file=f)
    
    ## Part 4: Diplotype Detail
    print('<h2 id="diplotype">Diplotype Detail</h2>', file=f)
    print('<p>CPAT predicts the diplotypes for PGx genes with star (*) or named allele information, which include ABCG2, CACNA1S, CFTR, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP3A4, CYP3A5, CYP4F2, DPYD, G6PD, MT-NR1, NUDT15, RYR1, SLCO1B1, TPMT, UGT1A1, and VKORC1. CPAT assumes that no variation occurs for the missing positions in the submitted VCFs, and the final inferred diplotypes are a composite ranking result.</p>', file=f)
    for gene in list(dic_diplotype.keys()):
      print('<h3>%s: %s</h3>' % (gene, dic_diplotype[gene]['step2_res']), file=f)
      gene_detail = dic_diplotype[gene]['detail']
      header = '<table id="customers" border="1" cellspacing="0">\n<tr><th>Position</th><th>rsID</th><th>Effect on Protein</th><th>Definition of Alleles</th><th>Detected Alleles</th></tr>'
      for pos_res in gene_detail:
        position = pos_res[0] + ':' + pos_res[1] + ':' + pos_res[2]
        header = header + '\n<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (position, pos_res[3], pos_res[4], pos_res[5], pos_res[6])
      header = header + '\n</table>'
      print(header, file=f)
    
    ## Part 5: PGx Annotation
    print('<h2 id="annotation">PGx Annotation</h2>', file=f)
    print('<p>CPAT annotates PGx-related single-locus genotypes and predicted diplotypes using the knowledge from PharmGKB, which are the basis for calculating the PGx summary results.</p>', file=f)
    detail = clinical_anno_table.drop(columns=['Class']).reset_index(drop = True)
    header = '<table id="customers" border="1" cellspacing="0">\n<tr><th></th><th>Gene</th><th>Variant</th><th>Alleles</th><th>Drug</th><th>Evidence</th><th>Category</th><th>Function</th></tr>'
    for index, row in detail.iterrows():
      header = header + '\n<tr><td>%s <a href=%s target="_blank" style="color: #ffb30e"><i class="fa-solid fa-up-right-from-square"></i></a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (index+1, row.URL, row.Gene, row.Variant, row.Alleles, row.Drug, row.EvidenceLevel, row.PhenotypeCategory, row.Function)
    header = header + '\n</table>'
    print(header, file=f)
    
    # Part 4: Disclaimers
    disclaimer = """
    <h2 id="about">Disclaimers</h2>
    <p>The report incorporates analyses of peer-reviewed studies and other publicly available information identified by CPAT by State Key Laboratory of Genetic Engineering from the School of Life Sciences and Human Phenome Institute, Fudan University, Shanghai, China. These analyses and information may include associations between a molecular alteration (or lack of alteration) and one or more drugs with potential clinical benefit (or potential lack of clinical benefit), including drug candidates that are being studied in clinical research.<br>
    Note: A finding of biomarker alteration does not necessarily indicate pharmacologic effectiveness (or lack thereof) of any drug or treatment regimen; a finding of no biomarker alteration does not necessarily indicate lack of pharmacologic effectiveness (or effectiveness) of any drug or treatment regimen.<br>
    No Guarantee of Clinical Benefit: This Report makes no promises or guarantees that a particular drug will be effective in the treatment of disease in any patient. This report also makes no promises or guarantees that a drug with a potential lack of clinical benefit will provide no clinical benefit.<br>
    Treatment Decisions are Responsibility of Physician: Drugs referenced in this report may not be suitable for a particular patient. The selection of any, all, or none of the drugs associated with potential clinical benefit (or potential lack of clinical benefit) resides entirely within the discretion of the treating physician. Indeed, the information in this report must be considered in conjunction with all other relevant information regarding a particular patient, before the patient's treating physician recommends a course of treatment. Decisions on patient care and treatment must be based on the independent medical judgment of the treating physician, taking into consideration all applicable information concerning the patient's condition, such as patient and family history, physical examinations, information from other diagnostic tests, and patient preferences, following the standard of care in a given community. A treating physician's decisions should not be based on a single test, such as this test or the information contained in this report.<br>
    When using results obtained from CPAT, you agree to cite CPAT.<br><br>
    </p>
    </div>
    </body>
    </html>
    """
    print(disclaimer, file=f)

