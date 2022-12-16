# PUMA resequencing data
This is a repository for scripts I used in analyzing Purple Martin (*Progne subis*) resequencing data.
These scripts were executed using the Texas A&M High Preformance Research Computing resources (https://hprc.tamu.edu/).
<br/>
# Filtering SNPs
The _filter_ folder contains code for steps in filtering snps. The order is based on this particular data set. Quality filtering is in `02-filterSNPS_round1.sh`. Determining sex-linked scaffolds was done through a couple methods in this project, first by mapping reference genome to chicken genome using satsuma (not listed here), and using fst between males and females to find snp positions driving large divide in pca between male/female. For population-related analyses, linkage-disequilibrium (LD) pruning was done prior to fst and pcas. For gwas and other analyses the snps were not pruned for LD, and positions were ID'ed again for sex-linked scaffolds in order to filter for analyses requiring autosomes only. 
<br/>
# Analyses

**Population genomics**: 
* _admixture_ for looking at ancestry 
* _fst_ to examine population differentiation
* _pca_ for population structure  

**Demographic history**:
* _fastphase_ to phase data to prep for msmc to test runs with groups of multiple individuals
* _msmc_ to run demographic history analysis, incuding bootstrapping 

**Genome-wide associations**:
* _gwas_ for running lmms and bslmms with GEMMA, ordering chromosomes and plotting, running ETI, and regressing out covariates in phenotype for specific models to include covariates in bslmm
* _fst_ for differentiation across genome
* _heritability_ for plotting PVE, PGE, n_gamma outputs from GEMMA. The actual code for caluclating heritability comes with the bslmms models in _gwas_.The `01-plot_pves.R` is really the only code relavent here for plotting, but the folder contains other code when I was trying out GCTA
* _lfmm_ for running lfmms with LEA  
* Polygenic score code and polygenic score accuracy code by Matt Thorstensen: https://github.com/BioMatt/PUMA_PGS

**Other**:
* _ld_ for characterizing linkage disequilibrium
