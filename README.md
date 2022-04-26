## ALS Patient Stratification Analysis
R and Python scripts used in the ALS Patient Stratification analysis by Eshima, O'Connor, Marschall, NYGC ALS Consortium, Bowser, Plaisier, and Smith

Author: Jarrett Eshima

PI: Dr. Barbara Smith

Institute: School of Biological and Health Systems Engineering, Arizona State University

Short description: This file contains the instructions to replicate the ALS Patient Stratification analysis with GEO Accession GSE153960.

Note:
1. Raw RNA-seq files can be downloaded from NCBI or the EMBL-EBI mirror.
2. The NIH's Globus software was employed to leverage high-speed download capabilities.
3. SQuIRE [1] was run on Arizona State University's Agave High Performance Computing Cluster

## ALS Patient Stratification Analysis - Code Order

1) SQuIRE[1] for Transposable Element quantification using RNA-seq .BAM or .FASTQ files
2) ALSPatientStratification_SQuIRE_PostProcessing.R
3) ALSPatientStratification_SQuIRE_Meta.R
4) ALSPatientStratification_DifferentialExpression_ALSPatients.R
5) ALSPatientStratification_ClusterPrep.R
6) ALSPatientStratification_Estimate_NMF_Rank.R
7) ALSPatientStratification_SAKE.R (This file also contains download instructions to obtain the SAKE [2] version used in this study)
8) ALSPatientStratification_PostSAKE.R
9) ALSPatientStratification_EnrichmentFeatureScore.R
10) ALSPatientStratification_TEQuant_DifferentialExpression_Controls.R (Cleans up control subject count information for GSEA and univariate analysis)
11) ALSPatientStratification_GSEA_Prep.R
12) Run GSEA [3] (This analysis can be paired with Cytoscape [4] for construction of gene network)
13) ALSPatientStratification_GSEA_CustomHeatmap.R 
14) ALSPatientStratification_WGCNA.R [5]
15) ALSPatientStratification_MetaData_Extractor.R
16) ALSPatientStratification_survival.R
17) ALSPatientStratification_UnivariateAnalysis.R
18) ALSPatientStratification_BootstrapClassification.R

## Supplemental Code

Supervised Classification - Code Order
1) ALSPatientStratification_SupervisedClassificationFeatureScore.R (Feature selection for machine learning)
2) ALSPatientStratification_Pheno4loom.R (Generates phenotype dataframe for .loom file development)
3) ALSPatientStratification_ReadLoom.R (Generates .loom files for compatibility with plaisier-lab [6] supervised classification scripts)
4) Python scripts* for development of k-nearest neighbor, multilayer perceptron, random forest, and support vector machine classifiers using Scikit-learn framework [7]
5) ALSPatientStratification_CVresults_to_F1report.R* 
6) ALSPatientStratification_Supervised_Classification_PostProcessing.R

Supporting
1) ALSPatientStratification_PatientDemographics.R
2) ALSPatientStratification_TEcluster_WGCNA.R (aids with the decision to avoid collapsing locus-specific TEs to subfamilies for custom gene enrichment)
3) ALSPatientStratification_Subtype_DiscriminatoryFeature_Assignment.R (compares clustering-assigned subtype genes between the NovaSeq and HiSeq cohorts)


*A few notes on supervised classification:

* Plaisier Lab scripts [6]: classifiersV3.py ; cvClassification_FullAnalysis.py (ACTINN Neural Network not used in this analysis)

* Smith Lab scripts: KNN_CV_TrainTest.py, KNN_Holdout.py, MLP_CV_TrainTest.py, MLP_Holdout.py, RF_CV_TrainTest.py, RF_Holdout.py, SVM_CV_TrainTest.py, and SVM_Holdout.py (basic python scripts for supervised classifier development with cross validation - does not generate a classification "report" file)

* Supplemental Smith Lab scripts: ALSPatientStratification_CVresults_to_F1report.R (Generates a "report" file with Precision, Recall, and F1 metrics, based on supervised classification CV results)


References:
[1] Yang, W. R., Ardeljan, D., Pacyna, C. N., Payer, L. M., & Burns, K. H. (2019). SQuIRE reveals locus-specific regulation of interspersed repeat expression. Nucleic acids research, 47(5), e27-e27.
[2] Ho, Y. J., Anaparthy, N., Molik, D., Mathew, G., Aicher, T., Patel, A., ... & Hammell, M. G. (2018). Single-cell RNA-seq analysis identifies markers of resistance to targeted BRAF inhibitors in melanoma cell populations. Genome research, 28(9), 1353-1363.
[3] Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences, 102(43), 15545-15550.
[4] Shannon, P., Markiel, A., Ozier, O., Baliga, N. S., Wang, J. T., Ramage, D., ... & Ideker, T. (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome research, 13(11), 2498-2504.
[5] Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), 1-13.
[6] https://github.com/plaisier-lab/U5_hNSC_Neural_G0
[7] Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.
