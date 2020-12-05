# PUMA resequencing data
This is a repository for scripts I used in analyzing Purple Martin (*Progne subis*) resequencing data.
These scripts were executed using the Texas A&M High Preformance Research Computing resources (https://hprc.tamu.edu/).
<br/>
# Filtering SNPs
The _filter_ folder contains code for steps in filtering snps. The order is based on this particular snp set and on the first run through. Quality filtering is in `02-filterSNPS_round1.sh`. Determining sex-linked scaffolds was done through a couple methods in this project, first by mapping reference genome to chicken genome using satsuma (not listed here), and using fst between males and females to find snp positions driving large divide in pca between male/female. For population-related analyses, linkage-disequilibrium (LD) pruning was done prior to fst and pcas. For gwas and other analyses the snps were not pruned for LD, and positions were ID'ed again for sex-linked scaffolds in order to filter for analyses needing autosomes only. 
<br/>
# Analyses

**Population genomics**: 
* _admixture_ for looking at ancestry 
* _fst_ to examine population differentiation
* _pca_ for population structure<br/>
**Demographic history**:
* _fastphase_ to phase data to prep for msmc
* _msmc_ to run demographic history analysis<br/>
**Genome-wide associations**:
* _gwas_ for running lmms and bslmms with GEMMA, ordering chromosomes and plotting, running ETI, and regressing out covariates in phenotype for specific models to include covariates in bslmm
* _heritability_ for plotting PVE, PGE, n_gamma outputs from GEMMA
* _lfmm_ for running lfmms with LEA <br/>
**Other**:
* _ld_ for characterizing linkage disequilirium
