# Menegollo_Bentham
Code and data in Menegollo et al., Cancer Research, 2024


# 1. Identifying biclusters from METABRIC and TCGA using the MCbiclust package

## 1.1. 
METABRIC data (https://doi.org/10.7303/syn1688369) was analysed using METABRIC_DATA.RData, METABRIC_gene_lists.RData as inputs.

`METABRIC_mito_legion.R` finds the highly correlated seeds.

`METABRIC_legion_analysis.R` was used to load the results from the above script and use the silhouette method to identify similar biclusters and then calculate the average clusters which can then be used to sort the samples. (output: `METABRIC_sort_data.RData`)

`METABRIC_mito_sort.R` loads a set of average correlation vectors generated from `METABRIC_legion_analysis.R` and sorts the samples by this bicluster.

`METABRIC_sort_analysis.R` is used to analyse the resulting biclusters. It plots the CVplot comparing the different outputs and also runs a GSEA on the correlation vectors. (output: `METABRIC_PC1_GSEA.RData`)

`compile_all_METABRIC_data_and_derived_parameters_all_biclusters.R` was used to generate the final METABRIC dataset, used for scripts generating the figure panels. 

As additional input, `Pommier_dev_genesets.xlsx`, `GS metabolic genes list.xlsx`, [`Mitocarta 2.0`](https://www.broadinstitute.org/files/shared/metabolism/mitocarta/human.mitocarta2.0.html), `METABRIC_estimate_score.gct`, `CL_CoreCL_METABRIC_Fougner.xlsx` were used. (output: `METABRIC_starting_data_final_corr_groups.Rdata`).

#### Large files are shared on Gdrive: [METABRIC_DATA.RData](https://drive.google.com/file/d/171f94kgKDdWKvPSNAYQX4ORyOZ4RVzGX/view?usp=sharing), [METABRIC_starting_data_final_corr_groups.Rdata](https://drive.google.com/file/d/17A6JRdsL07RIZtiNpYS5HgokbBpdG6h6/view?usp=sharing)

#### *Notes:*
*`METABRIC_mito_legion.R` and `METABRIC_mito_sort.R` need to be run on a high preformance cluster.*

*`METABRIC_mito_legion.R` accepts one argument which is for set.seed() to set the random number generator.*

*The initial bicluster nomenclature is maintained in these scripts:*

*ICT1 biclusters resulted from an ICT1 correlated geneset, Mito1-3 biclusters resulted from the entire Mitocarta1.0 geneset, Rand1-2 resulted from using random genesets as starting points. The resulting biclusters were shown to have the following relationships:*

*`ICT1 ~ Mito2`*
*`Rand1 = Mito1`*
*`Rand2 = Mito3`*

*During further analysis (as in METABRIC_starting_data_final_corr_groups.Rdata) a new nomenclature was introduced for clarity: the following three biclusters were further analysed:* 

*`MB1 = ICT1`* 
*`MB2 = Mito2`* 
*`MB3 = Mito1`*

### 1.2 
Biclusters in TCGA data (NCBI dbGaP at phs000178.v11.p8) were identified using the `TCGA_RNASeq_Analysis_MB1.R`, `TCGA_RNASeq_Analysis_MB2.R`, `TCGA_RNASeq_Analysis_MB3.R` scripts.

These take `Legion_BRCA_RNAseq1.RData` as input and generate `TCGA_MB1_RNASeq_data_nonorm.RData`, `TCGA_MB2_RNASeq_data_nonorm.RData`, `TCGA_MB3_RNASeq_data_nonorm.RData`. These are datafiles by biclusters, and are unified in `TCGA.all.biclusters.RNAseq.Rdata` by the `Fig_7ABC_S7AB_Histological_and_Genetic_Analysis.R script`.

#### Large files are shared on Gdrive:[Legion_BRCA_RNAseq1.RData](https://drive.google.com/file/d/17RD79BPvJQlTYTPXq3lUMHzLu3yG_SBW/view?usp=sharing)

# 2. Further analysis of the METABRIC and TCGA biclusters, scripts generating figure panels
### 2.1 METABRIC further analysis - figures
		`Fig_1DE_3B_Three_way_forks.R`
		`Fig_2ABCD_S2B_pathway_CV_column_plots.R`
		`Fig_2EFGH_pathway_activities_in_forks_vs_PAM50.R`
		`Fig_1GF_CV_Complex_Heatmap_Geneset_Enrichment.R`
		`Fig_1B_Two_way_forks.R`
### 2.2 TCGA further analysis - figures
		`Fig_3A_S3_TCGA_ParadigmAnalysis.R`
  		`Fig_7ABC_S7AB_Histological_and_Genetic_Analysis.R'

