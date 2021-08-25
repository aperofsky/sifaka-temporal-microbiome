Microbiome and statistical analyses are performed with the statistical computing software R (https://www.r-project.org/). The dataset is comprised of custom R scripts and inputs to reproduce the main statistical analyses and figures in Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes." 

The analysis starts with pre-processed phyloseq objects and sifaka metadata and is split into eleven steps: 

1_PCOA_figs.R

Principal coordinates analysis of sifaka gut microbiome Bray-Curtis and weighted Unifrac dissimilarities. Sources Fig2_home_range_maps.R to make Figure 2 in the publication. 

2_CST_figs.R

Partitions sifaka gut microbiome taxonomic profiles into community state types (CSTs). 

3_CST_alpha_diversity.R

Compares Chao1 richness and Simpson's diversity across microbial CSTs. 

4_CST_samr_differential_abundance.R

Uses the samr R package to estimate differential abundance of bacterial phyla and genara across microbial CSTs.

5_Microbiome_Random_Forest.R

Random forest classifier models that use microbial ASVs to predict the year and social group associated with individual microbiome samples. 

6_PERMANOVA.R

Permutational multivariate analysis of variance (PERMANOVA) tests to assess variation in gut microbiome beta diversity according to study year, social group affiliation, host age, host sex, and individual identity. 

7_Fig4_within_and_between_individual_distances.R 

Assesses differences in pairwise microbial dissimilarity among six categories of intra- and inter-individual sample comparisons.

8_between_individual_GLMMs.R

Uses Beta generalized linear mixed models (GLMMs) to estimate the effects of social relationships, group membership, genetic relatedness, and shared diet on pairwise gut microbial similarity.

9_within_host_GLMMs.R 

Uses Beta GLMMs to estimate the effects of demographic predictors (age, sex, and dispersal history) on gut microbiome dynamics in individual animals. 

10_immigrant_resident_GLMMs.R

Uses Beta GLMMs to measure the effects of group tenure (i.e., length of time in the group) and social relationships on microbial dissimilarity between immigrants and long-term group residents.

11_social_partner_stability_GLMMs.R

Uses a Beta GLMM to test the effects of social partner stability on within-host gut microbiome dynamics. 

Other scripts

Fig2_home_range_maps.R Plots pre-computed home ranges of six sifaka social groups. 

Folders

The figures folder contains the output of the analysis scripts. 

The Rdata folder contains all data inputs necessary to run the analysis scripts. 

Pre-processed phyloseq objects with ASV count tables, ASV taxonomic classifications, ASV phylogenies, and sample metadata:
sifaka_allotus_phyloseq_openref_2016.RData (raw ASV counts)
sifaka_trimmed_phyloseq_normalized_dada2.RData (normalized ASV counts)
sifaka_trimmed_phyloseq_normalized_dada2_CSTs.RData (normalized ASV counts with CST assignments)

Pre-computed home ranges for sifaka social groups I to VI in each study year:
home_range_gps_df_2012.Rdata
home_range_gps_df_2015.Rdata
home_range_gps_df_2016.Rdata

Input dataframes for analysis of immigrant-resident microbial dissimilarity
imm_group_tenure_df_lim.rds
imm_res_diss_df_lim.rds

Input dataframe for analysis of within- and between-individual microbial dissimilarity (Figure 4 in manuscript): within_and_bw_host_df_lim.rds

Input for analysis of between-individual microbial dissimilarity: pairwise_social_behavior_and_microbiome_df_lim.rds

Input dataframes for analysis of within-host microbial dynamics: pairwise_predictors_same_ind_lim.RData

Input for analysis of social partner stability and within-host gut microbiome dynamics: PSI_vs_microbiome_df_lim.rds

Outputs of random forest models: 
sifaka_RF_fit_group_classification.Rdata
sifaka_RF_fit_year_classification.Rdata

