# CPAT

Annotation of pharmacogenomics relevant genotypes is a key component of clinical genomic testing. However, the inherent limitations of the currently widely used next generation sequencing make accurate calling of some genotypes still challenging, especially diplotypes. Therefore, <span style="color:#20376D"><b>we developed an automated annotation tool, CPAT, which takes the germline variant calling format (VCF) file and population information as input and outputs an HTML report of drug responses with prescription recommendations.</b></span> CPAT implements a ranking model based on the allele definition and population frequency for the inference of diplotype, and its performance is validated in comparison with GeT-RM consensus and four other similar tools. CPAT further constructed an annotation model based on PharmGKB knowledge to summarize the <b>level of drug response</b> (<b>decreased</b>, <b>moderate</b>, <b>increased</b>) and the <b>level of evidence</b> (<b>level A</b>, <b>level B</b>) for the resolved genotypes. In summary, CPAT is able to resolve, annotate and report germline variants of an individual, providing an end-to-end solution for clinical pharmacogenomics decision support.
<p align="center">
<img src="./assets/cpat_architecture.png" width="40%" />
</p>

## Status
CPAT is still under _active development_. In the current release, you should only use it to evaluate whether CPAT will compile and run properly on your system. All information in the CPAT report is interpreted directly from the uploaded VCF file. Users recognize that they are using CPAT at their own risk.
## Prerequisite
- Bash
- Python3 >= 3.6
## Usage
```Bash
git clone https://github.com/premedkb/cpat.git
cd cpat
python cpat.py -s sample_id -i germline_vcf -p population -o outdir
```
## Examples
The test VCF files of 1000 Genomes Project are stored in _./data/vcf_ directory, and the corresponding CPAT reports are stored in _./data/report_ directory.
## Input data
### VCF file
As the diplotype definitions only match to the human genome _GRCh38_ and given its increasing generality, CPAT requires that the VCF file is based on _GRCh38_.

CPAT directly uses the NGS-derived VCF file as input and assumes that it has undergone quality control. Therefore, if the VCF file is of poor quality, inaccurate genotype resolution results and inappropriate clinical recommendations may be reported.
### Population
There are nine biogeographic groups provided by CPAT: **AAC** (African American/Afro-Caribbean), **AME** (American), **EAS** (East Asian), **EUR** (European), **LAT** (Latino), **NEA** (Near Eastern), **OCE** (Oceanian), **SAS** (Central/South Asian), **SSA** (Sub-Saharan African). More information is available at https://www.pharmgkb.org/page/biogeographicalGroups.

Please use the *three-letter abbreviation* as input. This is to prevent errors caused by special symbols such as spaces.
## CPAT Models
### CPAT ranking model for diplotype inference
The aim of genotype resolution is to extract the alleles of small variants and the diplotypes related to PGx from the user-submitted VCF file. CPAT processes the “GT” information to obtain all relevant single-locus genotypes. Afterwards, the genotypes of small variants will be passed to clinical annotation directly, while the genotypes related diplotype definitions will be passed to CPAT ranking model. The output diplotypes with the highest ranking will then be annotated.
<p align="center">
<img src="./assets/cpat_genotype_resolution.png" width="80%" />
</p>

### CPAT annotation model for predicting drug response at individual level
The aim of this component is to resolve the "drug-genotype-response-evidence" relationship. CPAT annotation model translates the literal PGx knowledge about genotypes into quantitative scores. The association between multiple genotypes and a single drug is then further translated into an individual-level association with this drug. Ultimately, individual responses to specific drugs are reported in terms of the strength of the response and the reliability of the evidence.
<p align="center">
<img src="./assets/cpat_clinical_annotation.png" width="60%" />
</p>