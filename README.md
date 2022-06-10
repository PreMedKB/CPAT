# CPAT

Annotation of pharmacogenomics relevant genotypes is a key component of clinical genomic testing. However, the inherent limitations of the currently widely used next generation sequencing make accurate calling of some genotypes still challenging, especially diplotypes. Therefore, <span style="color:#20376D"><b>we developed an automated annotation tool, CPAT, which takes the germline variant calling format (VCF) file and population information as input and outputs an HTML report of drug responses with prescription recommendations.</b></span> CPAT implements a ranking model based on the allele definition and population frequency for the inference of diplotype, and its performance is validated in comparison with GeT-RM consensus and four other similar tools. CPAT further constructed an annotation model based on PharmGKB knowledge to summarize the <b>level of drug response</b> (<b>decreased</b>, <b>moderate</b>, <b>increased</b>) and the <b>level of evidence</b> (<b>level A</b>, <b>level B</b>) for the resolved genotypes. In summary, CPAT is able to resolve, annotate and report germline variants of an individual, providing an end-to-end solution for clinical pharmacogenomics decision support.

<p align="center">
<img src="./assets/cpat_architecture.png" width="40%" />
</p>

## Status
CPAT is still under _active development_.
## Prerequisite
- Bash
- Python3 >= 3.6
## Input data
### VCF file
As the diplotype definitions only match to the human genome GRCh38 and given its increasing generality, CPAT requires that the VCF file is based on GRCh38.

CPAT directly uses the NGS-derived VCF file as input and assumes that it has undergone quality control. Therefore, if the VCF file is of poor quality, inaccurate genotype resolution results and inappropriate clinical recommendations may be reported.

### Population
There are nine biogeographic groups, i.e., African American/Afro-Caribbean, American, Central/South Asian, East Asian, European, Latino, Near Eastern, Oceanian, Sub-Saharan African.


## CPAT Models
### CPAT ranking model for diplotype inference
<p align="center">
<img src="./assets/cpat_genotype_resolution.png" width="100%" />
</p>

### CPAT annotation model for predicting drug response at individual level
<p align="center">
<img src="./assets/cpat_clinical_annotation.png" width="60%" />
</p>