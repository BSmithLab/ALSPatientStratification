# ALSPatientStratification
R and Python scripts used in the ALS Patient Stratification analysis

Author: Jarrett Eshima
PI: Dr. Barbara Smith
Institute: Arizona State University, School of Biological and Health Systems Engineering

Short description: This file contains the instructions to replicate the ALS Patient Stratification analysis with GEO Accession GSE153960.

Note:
1. Raw RNA-seq files can be downloaded from NCBI or the EMBL-EBI mirror.
2. The NIH's Globus software was employed to leverage high-speed download capabilities.
3. SQuIRE (1) was run on Arizona State University's Agave High Performance Computing Cluster

Script Order:
1. SQuIRE for Transposable Element quantification using RNA-seq .BAM or .FASTQ files
2. ALSPatientStratification_SQuIRE_PostProcessing.R
3. ALSPatientStratification_SQuIRE_Meta.R
4. ALSPatientStratification_DifferentialExpression_ALSPatients.R
5. ALSPatientStratification_ClusterPrep.R
6. ALSPatientStratification_Estimate_NMF_Rank.R
7. ALSPatientStratification_SAKE.R (This file also contains download instructions to obtain the SAKE (2) version used in this study)
8. ALSPatientStratification_FeatureScore.R 



(1) Yang, W. R., Ardeljan, D., Pacyna, C. N., Payer, L. M., & Burns, K. H. (2019). SQuIRE reveals locus-specific regulation of interspersed repeat expression. Nucleic acids research, 47(5), e27-e27.
(2) Ho, Y. J., Anaparthy, N., Molik, D., Mathew, G., Aicher, T., Patel, A., ... & Hammell, M. G. (2018). Single-cell RNA-seq analysis identifies markers of resistance to targeted BRAF inhibitors in melanoma cell populations. Genome research, 28(9), 1353-1363.
