
## R scripts and inputs to reproduce the main results and figures in Perofsky _et al._ 2021 "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes." _Molecular Ecology_
---
#### Microbiome and statistical analyses are performed with the statistical computing software [R](https://www.r-project.org/).

#### The analysis is split into eleven steps:

- 1_PCOA_figs.R

	- Principal coordinates analysis of sifaka gut microbiome Bray-Curtis dissimilarities and weighted Unifrac distances.

- 2_CST_figs.R

	- Partitions sifaka gut microbiome taxonomic profiles into community state types (CSTs).

- 3_CST_alpha_diversity.R

	- Compares Chao1 richness and Simpson's diversity across microbial CSTs.

- 4_CST_samr_differential_abundance.R

	- Estimates differential abundance of bacterial phyla and genera across microbial CSTs.

- 5_Microbiome_Random_Forest.R

	- Random forest classifier models that use microbial ASVs to predict the year and social group associated with individual microbiome samples.

- 6_PERMANOVA.R

	- Permutational multivariate analysis of variance (PERMANOVA) tests to assess variation in gut microbiome beta diversity according to study year, social group affiliation, host age, host sex, and individual identity.

- 7_Fig4_within_and_between_individual_distances.R

	- Assesses differences in pairwise microbial dissimilarity among six categories of intra- and inter-individual sample comparisons.

- 8_between_individual_GLMMs.R

	- Uses Beta generalized linear mixed models (GLMMs) to estimate the effects of social relationships, group membership, genetic relatedness, and shared diet on pairwise gut microbial similarity.

- 9_within_host_GLMMs.R

	- Uses Beta GLMMs to estimate the effects of demographic predictors (age, sex, and dispersal history) on gut microbiome dynamics in individual animals.

- 10_immigrant_resident_GLMMs.R

	- Uses Beta GLMMs to measure the effects of group tenure (i.e., length of time spent in the social group) and social relationships on microbial dissimilarity between immigrants and long-term group residents.

- 11_social_partner_stability_GLMMs.R

	- Uses a Beta GLMM to test the effects of social partner stability on within-host gut microbiome dynamics.

### Data dictionary
- perofsky_2021_mol_ecol_data_dictionary.xlsx
	- Tabs define variables ("Column"), allowable values for variables ("Value"), whether variables contain sample or host information ("Category"), and variable definitions ("Explanation").

### Folders

#### The _figures_ folder contains figures created with the analysis scripts.

#### The _Rdata_ folder contains all data inputs necessary to run the analysis scripts.

##### Pre-processed phyloseq objects with ASV count tables, ASV taxonomic classifications, ASV phylogenies, and sample metadata:
- sifaka_allotus_phyloseq_openref_2016.RData (raw ASV counts)
- sifaka_trimmed_phyloseq_normalized_dada2.RData (normalized ASV counts)
- sifaka_trimmed_phyloseq_normalized_dada2_CSTs.RData (normalized ASV counts with CST assignments)

##### Input dataframes for analysis of immigrant-resident microbial dissimilarity:
- imm_group_tenure_df_lim.rds
- imm_res_diss_df_lim.rds

##### Input data frame for analysis of withinand between-individual microbial dissimilarity:
- within_and_bw_host_df_lim.rds

##### Input data frame for analysis of between-individual microbial dissimilarity:
- pairwise_social_behavior_and_microbiome_df_lim.rds

##### Input data frame for analysis of within-host microbial dynamics:
- pairwise_predictors_same_ind_lim.RData

##### Input data frame for analysis of social partner stability and within-host gut microbiome dynamics:
- PSI_vs_microbiome_df_lim.rds

##### Outputs of random forest models:
- sifaka_RF_fit_group_classification.Rdata
- sifaka_RF_fit_year_classification.Rdata