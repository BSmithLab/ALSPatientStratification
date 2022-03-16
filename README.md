## ALS Patient Stratification Analysis
R and Python scripts used in the ALS Patient Stratification analysis by Eshima et al.

Author: Jarrett Eshima

PI: Dr. Barbara Smith

Institute: School of Biological and Health Systems Engineering, Arizona State University

Short description: This file contains the instructions to replicate the ALS Patient Stratification analysis with GEO Accession GSE153960.

Note:
1. Raw RNA-seq files can be downloaded from NCBI or the EMBL-EBI mirror.
2. The NIH's Globus software was employed to leverage high-speed download capabilities.
3. SQuIRE (1) was run on Arizona State University's Agave High Performance Computing Cluster

## ALS Patient Stratification Analysis - Code Order

1) SQuIRE[1] for Transposable Element quantification using RNA-seq .BAM or .FASTQ files
2) ALSPatientStratification_SQuIRE_PostProcessing.R
3) ALSPatientStratification_SQuIRE_Meta.R
4) ALSPatientStratification_DifferentialExpression_ALSPatients.R
5) ALSPatientStratification_ClusterPrep.R
6) ALSPatientStratification_Estimate_NMF_Rank.R
7) ALSPatientStratification_SAKE.R (This file also contains download instructions to obtain the SAKE [2] version used in this study)
8) ALSPatientStratification_PostSAKE.R
9) ALSPatientStratification_FeatureScore.R (Feature selection for machine learning)
10) ALSPatientStratification_Pheno4loom.R (Generates phenotype dataframe for .loom file development)
11) ALSPatientStratification_ReadLoom.R (Generates .loom files for compatibility with plaisier-lab [3] supervised classification scripts)
12) Plaisier Lab Python scripts* for development of random forest, support vector machine, and k-nearest neighbor classifiers using Scikit-learn [4]
13) Smith Lab Python script* for development of multilayer perceptron classifier [4]
14) ALSPatientStratification_MLP_CV_to_F1.R* 
15) ALSPatientStratification_Supervised_Classification_PostProcessing.R
16) ALSPatientStratification_TEQuant_DifferentialExpression_Controls.R (Cleans up control subject count information for GSEA and univariate analysis)
17) ALSPatientStratification_GSEA_Prep.R
18) Run GSEA [5] (This analysis can be paired with Cytoscape [6] for construction of gene network)
19) ALSPatientStratification_GSEA_CustomHeatmap.R 
20) ALSPatientStratification_WGCNA.R [7]
19) ALSPatientStratification_MetaData_Extractor.R
20) ALSPatientStratification_survival.R
21) ALSPatientStratification_UnivariateAnalysis.R


A few notes on supervised classification:
* Required Plaisier Lab scripts [3]: classifiersV3.py ; cvClassification_FullAnalysis.py (ACTINN Neural Network not used in this analysis)

* Required Smith Lab scripts: SVM.py, RF.py, KNN.py, and MLP_CV.py (very basic python script for multilayer perceptron development with CV - does not generate a classification "report" file)

* Supplemental Smith Lab scripts: MLP_CV_to_F1.R (Generates a "report" file with Precision, Recall, and F1 metrics, based on MLP classification results)


## Supplemental Code
1) ALSPatientStratification_CustomTEGeneSet_Prep.R (aids with the development of a custom gene set for TEs, using RepBase [8])
2) ALSPatientStratification_PatientDemographics.R
3) ALSPatientStratification_TEcluster_WGCNA.R (aids with the decision to collapse locus-specific TEs to subfamily for gene enrichment)
4) Subtype_DiscriminatoryFeature_Assignment.R (compares clustering-assigned subtype genes between the NovaSeq and HiSeq cohorts)

References:
[1] Yang, W. R., Ardeljan, D., Pacyna, C. N., Payer, L. M., & Burns, K. H. (2019). SQuIRE reveals locus-specific regulation of interspersed repeat expression. Nucleic acids research, 47(5), e27-e27.
[2] Ho, Y. J., Anaparthy, N., Molik, D., Mathew, G., Aicher, T., Patel, A., ... & Hammell, M. G. (2018). Single-cell RNA-seq analysis identifies markers of resistance to targeted BRAF inhibitors in melanoma cell populations. Genome research, 28(9), 1353-1363.
[3] https://github.com/plaisier-lab/U5_hNSC_Neural_G0
[4] Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.
[5] Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences, 102(43), 15545-15550.
[6] Shannon, P., Markiel, A., Ozier, O., Baliga, N. S., Wang, J. T., Ramage, D., ... & Ideker, T. (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome research, 13(11), 2498-2504.
[7] Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), 1-13.
[8] Bao, W., Kojima, K. K., & Kohany, O. (2015). Repbase Update, a database of repetitive elements in eukaryotic genomes. Mobile Dna, 6(1), 1-6.
