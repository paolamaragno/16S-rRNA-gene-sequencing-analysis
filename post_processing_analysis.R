setwd("/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/delong_covid/analisi_micca_nuovi_vecchi_insieme/")

library('phyloseq')
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(tikzDevice)
library(psych) # corr.test
library(network)
library(GGally) # ggnet
library(corrplot)
library(vegan)
library(ape)
library(dplyr)
library(cowplot)
library(DESeq2)
library(ccrepe)
library(patchwork)
library(calibrate)
library(ggrepel)
library("randomForest")
library("plyr") # for the "arrange" function
library("rfUtilities") # to test model significance
library("caret")
library(pROC)
library(rstatix)
library(seqinr)
library(pheatmap)

#install.packages('Rcpp', dependencies = TRUE)
#require(devtools)
#install_version("ggplot2", version = "0.9.1", repos = "http://cran.us.r-project.org")

divprof_to_phyloseq <- function(otufilename, taxonomyfilename, samplefilename, treefilename) {
  argumentlist = list()
  if (!is.null(otufilename)) {
    otu_frame = read.csv(otufilename, sep="\t", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
    # otu table with samples as columns and ASV as rows
    otu = otu_table(otu_frame, taxa_are_rows=TRUE)
    # the table is converted in otu table object that is added to the list
    argumentlist <- c(argumentlist, list(otu))
  }
  if (!is.null(taxonomyfilename)) {
    taxa_frame = read.table(taxonomyfilename, sep="\t", header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)
    # for each ASV the taxonomy up to the genus
    taxa_split = strsplit(taxa_frame[,2], ";")
    names(taxa_split) = taxa_frame[,1]
    # list reporting for each ASV all its taxonomy separated as vector
    taxa = build_tax_table(sapply(taxa_split, parse_taxonomy_default))
    # the information is converted in a taxonomyTable object saved inside the list
    colnames(taxa) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:dim(taxa)[[2]]]
    argumentlist <- c(argumentlist, list(taxa))
  }
  if (!is.null(samplefilename)) {
    sample_frame = read.csv(samplefilename, sep="\t", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
    # table with the samples as rows and all the info including the depth (number of reads obtained for that sample) as columns
    sample = sample_data(sample_frame)
    argumentlist <- c(argumentlist, list(sample))
  }
  if (!is.null(treefilename)) {
    tree = read_tree(treefilename)
    if (is.null(tree)) {
      warning("treefilename failed import. It will not be included.")
    }
    else {
      argumentlist <- c(argumentlist, list(tree))
    }
  }
  
  do.call("phyloseq", argumentlist)
  # do.call to apply the same function (phyloseq) to each element of the list (argumentlist)
}
###################
ps_16S <- divprof_to_phyloseq("./denovo_unoise_sv/otutable.txt", "./denovo_unoise_sv/taxa.txt", "./denovo_unoise_sv/sampledata_new.txt","./denovo_unoise_sv/tree_rooted.tree")
# phyloseq object

all_samples_id <- unique(sample_data(ps_16S)$sample_id)
sample_id_run2 <- unique(sample_data(ps_16S)[sample_data(ps_16S)$run == 'run2']$sample_id)
sample_id_run1 <- unique(sample_data(ps_16S)[sample_data(ps_16S)$run == 'run1']$sample_id)
sample_id_only_both_run <- intersect(sample_id_run1,sample_id_run2) 

# remove all the resequenced samples from the second run (inlcuding S54780)
remove <- which((sample_data(ps_16S)$sample_id %in% sample_id_only_both_run) & (sample_data(ps_16S)$run =='run2'))
sample_data(ps_16S) <- sample_data(ps_16S)[-remove,]

all_samples_id <- unique(sample_data(ps_16S)$sample_id)
# 48
sample_id_run2 <- unique(sample_data(ps_16S)[sample_data(ps_16S)$run == 'run2']$sample_id)
# 10
sample_id_run1 <- unique(sample_data(ps_16S)[sample_data(ps_16S)$run == 'run1']$sample_id)
# 38

# save the otu table to do the last commands on the cluster
write.table(otu_table(ps_16S), file = './otutable_filtered.txt', quote = FALSE, sep = '\t')

dim(otu_table(ps_16S))
# 4077   96
dim(sample_data(ps_16S))
# 96 5

###### tax4fun placebo vs vsl
fasta_ASV <- read.fasta(file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/delong_covid/analisi_micca_nuovi_vecchi_insieme/denovo_unoise_sv/otus.fasta',as.string = TRUE,
                        forceDNAtolower = FALSE)

unname(unlist(fasta_ASV))
df <- data.frame(ASV_name = names(fasta_ASV), ref = unname(unlist(fasta_ASV)))
rownames(df) <- df[,1]
df[,1] <- NULL

otu_table_copy <- otu_table(ps_16S)
rownames(otu_table_copy) <- df[rownames(otu_table_copy),1]
write.table(otu_table_copy, file = './otutable_filtered_renamed.txt', quote = FALSE, sep = '\t')

taxa <- read.table('./denovo_unoise_sv/taxa.txt', sep = '\t')
taxa_split <- tidyr::separate(data = taxa, col =  "V2",into =  c("Kingdom", "Phylum", 'Class', 'Order', 'Family', 'Genus'), sep = ";", extra = "merge")

rownames(taxa_split) <- taxa_split[,1]
taxa_split[,1] <- NULL
rownames(taxa_split) <- df[rownames(taxa_split),1]
write.table(taxa_split, file = './taxa_split_renamed.txt', quote = FALSE, sep = '\t')

tax4fun_pathways <- read.csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/delong_covid/analisi_micca_nuovi_vecchi_insieme/tax4fun/pathwayprof_tax4fun2.csv')

rownames(tax4fun_pathways) <- tax4fun_pathways[,1]
tax4fun_pathways[,1] <- NULL

tax4fun_pathways_summary <- data.frame(taxa = rep(rownames(tax4fun_pathways), ncol(tax4fun_pathways)),
                                        sample = rep(colnames(tax4fun_pathways), each=nrow(tax4fun_pathways)),
                                        frequency = as.vector(as.matrix(tax4fun_pathways)))

tax4fun_pathways_summary$treatment <- sample_data(ps_16S)[tax4fun_pathways_summary$sample, 'treatment']$treatment
tax4fun_pathways_summary$time <- sample_data(ps_16S)[tax4fun_pathways_summary$sample, 'time']$time
tax4fun_pathways_summary$outcome <- sample_data(ps_16S)[tax4fun_pathways_summary$sample, 'outcome']$outcome

tax4fun_pathways_summary_t0 <- tax4fun_pathways_summary[tax4fun_pathways_summary$time =='T0',]
tax4fun_pathways_summary_t4 <- tax4fun_pathways_summary[tax4fun_pathways_summary$time =='T4',]

tax4fun_pathways_summary_t0$taxa <- as.factor(tax4fun_pathways_summary_t0$taxa)
tax4fun_pathways_summary_t0$sample <- as.factor(tax4fun_pathways_summary_t0$sample)
tax4fun_pathways_summary_t0$treatment <- as.factor(tax4fun_pathways_summary_t0$treatment)

stats_less_t0 <- tax4fun_pathways_summary_t0  %>% group_by(taxa) %>%
  wilcox_test(frequency ~ treatment, alternative = 'less', ref.group = 'VSL3') %>%
  adjust_pvalue(method = 'fdr') 

stats_less_sig_t0 <- stats_less_t0[stats_less_t0$p.adj <= 0.2,]

stats_greater_t0 <- tax4fun_pathways_summary_t0  %>% group_by(taxa) %>%
  wilcox_test(frequency ~ treatment, alternative = 'greater', ref.group = 'VSL3') %>%
  adjust_pvalue(method = 'fdr') 

stats_greater_sig_t0 <- stats_greater_t0[stats_greater_t0$p.adj <= 0.2,]

tax4fun_pathways_summary_t0_vsl <- tax4fun_pathways_summary_t0[tax4fun_pathways_summary_t0$treatment =='VSL3',]

stats_less_resp_t0 <- tax4fun_pathways_summary_t0_vsl  %>% group_by(taxa) %>%
  wilcox_test(frequency ~ outcome, alternative = 'less') %>%
  adjust_pvalue(method = 'fdr') 

stats_less_resp_sig_t0 <- stats_less_resp_t0[stats_less_resp_t0$p.adj <= 0.2,]

stats_greater_resp_t0 <- tax4fun_pathways_summary_t0_vsl  %>% group_by(taxa) %>%
  wilcox_test(frequency ~ outcome, alternative = 'greater') %>%
  adjust_pvalue(method = 'fdr') 

stats_greater_resp_sig_t0 <- stats_greater_resp_t0[stats_greater_resp_t0$p.adj <= 0.2,]

tax4fun_pathways_summary_t4$taxa <- as.factor(tax4fun_pathways_summary_t4$taxa)
tax4fun_pathways_summary_t4$sample <- as.factor(tax4fun_pathways_summary_t4$sample)
tax4fun_pathways_summary_t4$treatment <- as.factor(tax4fun_pathways_summary_t4$treatment)

stats_less_t4 <- tax4fun_pathways_summary_t4  %>% group_by(taxa) %>%
  wilcox_test(frequency ~ treatment, alternative = 'less', ref.group = 'VSL3') %>%
  adjust_pvalue(method = 'fdr') 

stats_less_sig_t4 <- stats_less_t4[stats_less_t4$p.adj <= 0.1,]
stats_less_sig_t4$taxa <- gsub('ko','map',stats_less_sig_t4$taxa)

effect_size_less_t4 <- qnorm(stats_less_sig_t4$p, lower.tail = TRUE)/sqrt(stats_less_sig_t4$n1+stats_less_sig_t4$n2)
names(effect_size_less_t4) <- stats_less_sig_t4$taxa

stats_greater_t4 <- tax4fun_pathways_summary_t4  %>% group_by(taxa) %>%
  wilcox_test(frequency ~ treatment, alternative = 'greater', ref.group = 'VSL3') %>%
  adjust_pvalue(method = 'fdr') 

stats_greater_sig_t4 <- stats_greater_t4[stats_greater_t4$p.adj <= 0.1,]
stats_greater_sig_t4$taxa <- gsub('ko','map',stats_greater_sig_t4$taxa)

effect_size_greater_t4 <- qnorm(1-stats_greater_sig_t4$p, lower.tail = TRUE)/sqrt(stats_greater_sig_t4$n1+stats_greater_sig_t4$n2)
names(effect_size_greater_t4) <- stats_greater_sig_t4$taxa

enriched_map_t4 <- data.frame(c(stats_greater_sig_t4$taxa, stats_less_sig_t4$taxa), 
                              c(rep("#ef857c",length(stats_greater_sig_t4$taxa)),
                                rep("#28658f",length(stats_less_sig_t4$taxa))),
                              'W18')

write.table(enriched_map_t4, file ='./tax4fun/pathways_fdr_0.1_t4.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

write.csv(stats_greater_sig_t4, file ='./tax4fun/pathways_fdr_0.1_t4.csv', quote = FALSE, row.names = FALSE)

# tax4fun responders vs non
tax4fun_pathways_summary_t4_vsl <- tax4fun_pathways_summary_t4[tax4fun_pathways_summary_t4$treatment =='VSL3',]

stats_less_resp_t4 <- tax4fun_pathways_summary_t4_vsl  %>% group_by(taxa) %>%
  t_test(frequency ~ outcome, alternative = 'less') %>%
  adjust_pvalue(method = 'fdr') 

stats_less_sig_resp_t4 <- stats_less_resp_t4[stats_less_resp_t4$p.adj <= 0.1,]

stats_greater_resp_t4 <- tax4fun_pathways_summary_t4_vsl  %>% group_by(taxa) %>%
  wilcox_test(frequency ~ outcome, alternative = 'greater') %>%
  adjust_pvalue(method = 'fdr') 

stats_greater_sig_resp_t4 <- stats_greater_resp_t4[stats_greater_resp_t4$p.adj <= 0.1,]

## Prima usi tax4fun per predire il contenuto metabolico sulla base delle sequenze 16s. 
# Poi identifichi quali pathway sono arricchiti in vsl/placebo al T4 e li rappresenti con iPath. Da ipath ho ricostruito quali KO (enzimi) 
# appartengono a questo pathway che è stato identificato come arricchito e ho fatto il wilcoxon per capire quali di questi KO sono arricchiti al 
# t4 in vsl/placebo e ho calcolato l’effect size

KO_functions <- read.csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/delong_covid/analisi_micca_nuovi_vecchi_insieme/tax4fun/functionalprof_tax4fun2.csv')

KO_functions_enriched <- KO_functions[KO_functions$X.NAME %in% c('K00830','K01620','K01733','K00872','K02203', 'K02204','K00003', 'K12524', 'K12525','K00133',
                                                                 'K00928', 'K12524', 'K12525', 'K12526','K00836','K13745','K06718','K06720','K10674','K00600',
                                                                 'K01079', 'K02203','K00831','K00058','K01834', 'K01837', 'K15633', 'K15634','K15635',
                                                                 'K00865', 'K11529','K00018','K00830'),]
rownames(KO_functions_enriched) <- KO_functions_enriched[,1]
KO_functions_enriched[,1] <- NULL

KO_functions_enriched_summary <- data.frame(KO = rep(rownames(KO_functions_enriched), ncol(KO_functions_enriched)),
                                       sample = rep(colnames(KO_functions_enriched), each=nrow(KO_functions_enriched)),
                                       frequency = as.vector(as.matrix(KO_functions_enriched)))

KO_functions_enriched_summary$treatment <- sample_data(ps_16S)[KO_functions_enriched_summary$sample, 'treatment']$treatment
KO_functions_enriched_summary$time <- sample_data(ps_16S)[KO_functions_enriched_summary$sample, 'time']$time
KO_functions_enriched_summary$outcome <- sample_data(ps_16S)[KO_functions_enriched_summary$sample, 'outcome']$outcome

KO_functions_enriched_summary_t4 <- KO_functions_enriched_summary[KO_functions_enriched_summary$time =='T4',]

stats_less_t4 <- KO_functions_enriched_summary_t4  %>% group_by(KO) %>%
  wilcox_test(frequency ~ treatment, alternative = 'less', ref.group = 'VSL3') %>%
  adjust_pvalue(method = 'fdr') 

stats_less_sig_t4 <- stats_less_t4[stats_less_t4$p.adj <= 0.1,]

stats_greater_t4 <- KO_functions_enriched_summary_t4  %>% group_by(KO) %>%
  wilcox_test(frequency ~ treatment, alternative = 'greater', ref.group = 'VSL3') %>%
  adjust_pvalue(method = 'fdr') 

stats_greater_sig_t4 <- stats_greater_t4[stats_greater_t4$p <= 0.1,]

stats_greater_sig_t4$KO <- c('K00058 - D-3-phosphoglycerate dehydrogenase/\n2-oxoglutarate reductase',
                             'K00133 - aspartate-semialdehyde dehydrogenase',
                             'K00831 - phosphoserine\naminotransferase',
                             'K00872 - homoserine kinase',
                             'K00928 - aspartate kinase',
                             'K01079 - phosphoserine\nphosphatase',
                             'K01733 - threonine synthase',
                             'K01834 - 2,3-bisphosphoglycerate-dependent\nphosphoglycerate mutase',
                             'K15634 - 2,3-bisphosphoglycerate-dependent\nphosphoglycerate mutase')

stats_greater_sig_t4 <- stats_greater_sig_t4[-9,]

effect_size_greater_t4 <- qnorm(1-stats_greater_sig_t4$p, lower.tail = TRUE)/sqrt(stats_greater_sig_t4$n1+stats_greater_sig_t4$n2)
names(effect_size_greater_t4) <- stats_greater_sig_t4$KO

effect_size_db <- data.frame(KO = stats_greater_sig_t4$KO, 
                             eff_size = effect_size_greater_t4[stats_greater_sig_t4$KO],
                             p = stats_greater_sig_t4$p)

effect_size_db$treatment <- ifelse(effect_size_db$eff_size > 0, 'VSL3', 'placebo')

effect_size_db$treatment <- factor(effect_size_db$treatment)

effect_size_db$KO <- factor(effect_size_db$KO, levels=effect_size_db[order(effect_size_db$eff_size),'KO'])
effect_size_db$p <- round(effect_size_db$p, 3)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/riccione_2024/effect_size_KO_t4.pdf', height = 3.5, width = 4.5)
ggplot(effect_size_db, aes(y=eff_size,x=KO, fill=treatment)) +
  geom_bar(stat = "identity", width = 0.7) +  coord_flip() +
  geom_text(aes(label = p), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(effect_size_db$eff_size >= 0, 1.2, -0.2)) +
  scale_fill_manual(values = c("#ef857c")) +
  theme_classic() +
  theme(axis.title.y = element_blank(), legend.title = element_blank(), legend.position = 'top') +
  ylab('Effect size') 
dev.off()

enriched_pathways_metabolic <- readxl::read_excel('./tax4fun/pathways_fdr_0.1_t4_only_metabolic.xlsx')

enriched_pathways_metabolic$taxa <- gsub('map','ko', enriched_pathways_metabolic$taxa)
names(effect_size_greater_t4) <- gsub('map','ko', names(effect_size_greater_t4))
stats_greater_sig_t4$taxa <- gsub('map','ko', stats_greater_sig_t4$taxa)

effect_size_db <- data.frame(taxa = paste(enriched_pathways_metabolic$taxa, '-', enriched_pathways_metabolic$pathway), 
                             eff_size = effect_size_greater_t4[enriched_pathways_metabolic$taxa],
                             p_adj = stats_greater_sig_t4[stats_greater_sig_t4$taxa %in% enriched_pathways_metabolic$taxa, 'p.adj'])

effect_size_db$treatment <- ifelse(effect_size_db$eff_size > 0, 'VSL3', 'placebo')

effect_size_db$treatment <- factor(effect_size_db$treatment)

effect_size_db$taxa <- factor(effect_size_db$taxa, levels=effect_size_db[order(effect_size_db$eff_size),'taxa'])
effect_size_db$p.adj <- round(effect_size_db$p.adj, 3)

pdf('./tax4fun/effect_size_metabolic_t4.pdf', height = 1.5, width = 6)
ggplot(effect_size_db, aes(y=eff_size,x=taxa, fill=treatment)) +
  geom_bar(stat = "identity", width = 0.7) +  coord_flip() +
  geom_text(aes(label = p.adj), position=position_dodge(width=0.9),size=3, 
            hjust = ifelse(effect_size_db$eff_size >= 0, 1.2, -0.2)) +
  scale_fill_manual(values = c("#ef857c")) +
  theme_classic() +
  theme(axis.title.y = element_blank(), legend.title = element_blank(), legend.position = 'top') +
  ylab('Effect size') 
dev.off()

#########

t <- otu_table(ps_16S)
class(t) <- "matrix"
tab <- t(t)
pdf("./R_analysis_all_run1/rarecurve.pdf")
rarecurve(tab, step=50, cex=0.5)
dev.off()
# rarefaction is a form of normalization of the samples so that they all have the same depth. each line is a sample

ps_16S_placebo <- subset_samples(ps_16S, treatment == 'placebo')
ps_16S_vsl <- subset_samples(ps_16S, treatment == 'VSL3')
ps_16S_T0 <- subset_samples(ps_16S, time == 'T0')
ps_16S_T4 <- subset_samples(ps_16S, time == 'T4')
ps_16S_T4_vsl <- subset_samples(ps_16S, time == 'T4' & treatment == 'VSL3')
ps_16S_T0_vsl <- subset_samples(ps_16S, time == 'T0' & treatment == 'VSL3')

ps.rarefied = rarefy_even_depth(ps_16S, rngseed=1, sample.size=min(sample_sums(ps_16S)), replace=F)
ps.rarefied_placebo = subset_samples(ps.rarefied, treatment == "placebo")
ps.rarefied_vsl = subset_samples(ps.rarefied, treatment == "VSL3")
ps.rarefied_T0 = subset_samples(ps.rarefied, time == "T0")
ps.rarefied_T4 = subset_samples(ps.rarefied, time == "T4")
ps.rarefied_T4_vsl = subset_samples(ps.rarefied, treatment == "VSL3" & time == 'T4')
ps.rarefied_T0_vsl = subset_samples(ps.rarefied, treatment == "VSL3" & time == 'T0')


# stacked barplot
taxtable6 = read.table('./denovo_unoise_sv/taxtables/taxtable6.txt', header = TRUE, sep='\t', row.names = 1)
sums <- colSums(taxtable6)

for (c in 1:ncol(taxtable6)) {
  for (r in 1:nrow(taxtable6)) {
    taxtable6[r,c] <- taxtable6[r,c]*100/sums[1]
  }
}

sums <- colSums(taxtable6)
taxtable6$mean_per_genus <- apply(data.frame(taxtable6),1,mean)
taxtable6 <- taxtable6[order(taxtable6$mean_per_genus, decreasing = TRUE),]
taxtable6$mean_per_genus <- round(taxtable6$mean_per_genus, 3)
taxtable6_filtered <- taxtable6[taxtable6$mean_per_genus >= 0.8,]

mean_genus_relative_frequency_VSL_T4 <- unlist(lapply(seq_along(1:nrow(taxtable6_filtered)), function(x) {
  mean(as.numeric(taxtable6_filtered[x,rownames(sample_data(ps_16S)[(sample_data(ps_16S)$treatment == 'VSL3' & sample_data(ps_16S)$time == 'T4'),])]))
}))
mean_genus_relative_frequency_VSL_T0 <- unlist(lapply(seq_along(1:nrow(taxtable6_filtered)), function(x) {
  mean(as.numeric(taxtable6_filtered[x,rownames(sample_data(ps_16S)[(sample_data(ps_16S)$treatment == 'VSL3' & sample_data(ps_16S)$time == 'T0'),])]))
}))
mean_genus_relative_frequency_placebo_T4 <- unlist(lapply(seq_along(1:nrow(taxtable6_filtered)), function(x) {
  mean(as.numeric(taxtable6_filtered[x,rownames(sample_data(ps_16S)[(sample_data(ps_16S)$treatment == 'placebo' & sample_data(ps_16S)$time == 'T4'),])]))
}))
mean_genus_relative_frequency_placebo_T0 <- unlist(lapply(seq_along(1:nrow(taxtable6_filtered)), function(x) {
  mean(as.numeric(taxtable6_filtered[x,rownames(sample_data(ps_16S)[(sample_data(ps_16S)$treatment == 'placebo' & sample_data(ps_16S)$time == 'T0'),])]))
}))

final_df <- data.frame(genus = rownames(taxtable6_filtered), VSL3_T4 = mean_genus_relative_frequency_VSL_T4,
                       VSL3_T0 = mean_genus_relative_frequency_VSL_T0,
                       placebo_T4 = mean_genus_relative_frequency_placebo_T4,
                       placebo_T0 = mean_genus_relative_frequency_placebo_T0)

others <- unlist(lapply(seq_along(1:4), function(x) {
  100 - sum(final_df[,x+1])
}))

final_df <- rbind(final_df, c('others',others))

df <- data.frame(genus = rep(final_df$genus, 4), frequency = c(final_df$VSL3_T4, final_df$VSL3_T0, final_df$placebo_T4, final_df$placebo_T0), 
                 info = c(rep('VSL3 - T4', length(final_df$genus)),rep('VSL3 - T0', length(final_df$genus)),rep('placebo - T4', length(final_df$genus)),rep('placebo - T0', length(final_df$genus))))

write.table(df, file = './denovo_unoise_sv/taxtables/taxtable6_stacked_barplot.txt', row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')

df <- data.frame(genus = rep(c('Bacteroides','Blautia', 'Faecalibacterium', 'Lachnospiracea_incertae_sedis', 'Ruminococcaceae uncl.',
                               'Lachnospiraceae uncl.', 'Alistipes', 'Prevotella', 'Roseburia', 'Ruminococcus', 'Gemmiger', 'Clostridiales uncl.',
                               'Streptococcus'	, 'Bifidobacterium', 'Anaerostipes',	'Dorea', 'Ruminococcus2',	'Fusicatenibacter',	'Escherichia/Shigella',
                               'Akkermansia', 'Romboutsia',	'Coprococcus', 'Parabacteroides',	'Barnesiella',	'Dialister',	'Erysipelotrichaceae uncl.',	
                               'Clostridium XlVa', 'Clostridium sensu stricto', 'others'), 4), frequency = c(final_df$VSL3_T4, final_df$VSL3_T0, final_df$placebo_T4, final_df$placebo_T0), 
                 info = c(rep('VSL3 - T4', length(final_df$genus)),rep('VSL3 - T0', length(final_df$genus)),rep('placebo - T4', length(final_df$genus)),rep('placebo - T0', length(final_df$genus))))

df$genus <- factor(df$genus, levels= c('Bacteroides','Blautia', 'Faecalibacterium', 'Lachnospiracea_incertae_sedis', 'Ruminococcaceae uncl.',
                                       'Lachnospiraceae uncl.', 'Alistipes', 'Prevotella', 'Roseburia', 'Ruminococcus', 'Gemmiger', 'Clostridiales uncl.',
                                       'Streptococcus'	, 'Bifidobacterium', 'Anaerostipes',	'Dorea', 'Ruminococcus2',	'Fusicatenibacter',	'Escherichia/Shigella',
                                       'Akkermansia', 'Romboutsia',	'Coprococcus', 'Parabacteroides',	'Barnesiella',	'Dialister',	'Erysipelotrichaceae uncl.',	
                                       'Clostridium XlVa', 'Clostridium sensu stricto', 'others'))
df$info <- factor(df$info, levels = c('placebo - T0', 'placebo - T4','VSL3 - T0', 'VSL3 - T4'))
df$frequency <- as.numeric(df$frequency)

n <- 31
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
colori <- col_vector
palette(colori)

pdf('./R_analysis_all_run1/stacked_barplot.pdf')
ggplot(df, aes(fill=genus, y=frequency, x=info)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(values = palette(colori)) + labs(y = "% Abundances") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=8), axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio=3) 
dev.off()

###alphas

alphas<-estimate_richness(otu_table(ps.rarefied), split = TRUE, measures = NULL)

df <- data.frame(treatment = sample_data(ps.rarefied)[rownames(alphas),]$treatment, 
                 outcome = sample_data(ps.rarefied)[rownames(alphas),]$outcome, 
                 alpha=alphas$Shannon,time=sample_data(ps.rarefied)[rownames(alphas),]$time)
df$info <- unlist(lapply(seq_along(1:nrow(df)), function(x) paste0(df[x,]$treatment, ' - ', df[x,]$time)))

df$treatment <- factor(df$treatment, levels = c('placebo','VSL3'))
df$time <- factor(df$time, levels = c('T0', 'T4'))
df$info <- factor(df$info, levels=c('placebo - T0', 'placebo - T4','VSL3 - T0', 'VSL3 - T4'))

pdf('./R_analysis_all_run1/alpha_diversity.pdf', width = 5, height =7)
ggplot(df,aes(x = info, y=alpha, colour=treatment)) +
  geom_boxplot(size=0.5, outlier.shape = NA) +
  theme_classic() +
  theme(aspect.ratio = 2.5, axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  geom_point(aes(colour=treatment,shape = time),position = position_jitter(width = 0.2, height = 0),alpha=1,stroke=0.6,size=4) +
  #facet_wrap(~type, nrow=1, ncol=2, scales = 'free') +
  ylab('Shannon index') +
  scale_colour_manual(values = c("#28658f", "#ef857c")) 
  #stat_pvalue_manual(stat.test,label = "p.adj", tip.length=0.01 )
dev.off()

df_vsl_t4 <- df[df$time=='T4' & df$treatment == 'VSL3',]
stat.test <-  df_vsl_t4 %>%
  wilcox_test(alpha ~ outcome) %>%
  adjust_pvalue(method = 'fdr') 
# not significant

df_vsl <- df[df$treatment == 'VSL3',]
stat.test <-  df_vsl %>%
  wilcox_test(alpha ~ outcome) %>%
  adjust_pvalue(method = 'fdr') 


###beta

bray_16S <- ordinate(ps.rarefied, method = "PCoA", distance = "bray")  ##bray-curtis

pdf('./R_analysis_all_run1/beta_bray_curtis.pdf')
plot_ordination(ps.rarefied, bray_16S, color="treatment", shape = 'time', title = 'Bray-Curtis, PCoA permANOVA') +
  geom_point(size =4) + theme(aspect.ratio=1) +
  xlab('PC1 [8.9%]') +
  ylab('PC2 [6.1%]') +
  theme_classic() +
  scale_colour_manual(values = c("#28658f", "#ef857c"))
dev.off()

dist_bray_T0 <-phyloseq::distance(ps.rarefied_T0, method="bray")

set.seed(1)
# distanza al T0 dipende dal trattamento?
permANOVA<-adonis2(dist_bray_T0 ~ sample_data(ps.rarefied_T0)$treatment, permutations = 999)  
# 0.576

dist_bray_T4 <-phyloseq::distance(ps.rarefied_T4, method="bray")

set.seed(1)
# distanza al T4 dipende dal trattamento?
permANOVA<-adonis2(dist_bray_T4 ~ sample_data(ps.rarefied_T4)$treatment, permutations = 999)  
# 0.049 *

dist_bray_vsl <-phyloseq::distance(ps.rarefied_vsl, method="bray")
# distanza dei trattati dipende dal time point?
set.seed(1)
permANOVA<-adonis2(dist_bray_vsl ~ sample_data(ps.rarefied_vsl)$time, permutations = 999)  
# 0.993

dist_bray_placebo <-phyloseq::distance(ps.rarefied_placebo, method="bray")
# distanza dei placebo dipende dal time point?
set.seed(1)
permANOVA<-adonis2(dist_bray_placebo ~ sample_data(ps.rarefied_placebo)$time, permutations = 999)  
# 1

dist_bray_vsl_t4 <-phyloseq::distance(ps.rarefied_T4_vsl, method="bray")
# distanza dei vsl al t4 dipende dall'outcome
set.seed(1)
permANOVA<-adonis2(dist_bray_vsl_t4 ~ sample_data(ps.rarefied_T4_vsl)$outcome, permutations = 999)  
# 0.713

dist_bray_vsl <-phyloseq::distance(ps.rarefied_vsl, method="bray")
# distanza dei vsl dipende dall'outcome
set.seed(1)
permANOVA<-adonis2(dist_bray_vsl ~ sample_data(ps.rarefied_vsl)$outcome, permutations = 999)  
# 0.028

dist_bray_vsl_t0 <-phyloseq::distance(ps.rarefied_T0_vsl, method="bray")
# distanza dei vsl dipende dall'outcome
set.seed(1)
permANOVA<-adonis2(dist_bray_vsl_t0 ~ sample_data(ps.rarefied_T0_vsl)$outcome, permutations = 999)  
# 0.703

bray_16S_vsl <- ordinate(ps.rarefied_vsl, method = "PCoA", distance = "bray")  ##bray-curtis

pdf('./R_analysis_all_run1/beta_bray_curtis.pdf')
plot_ordination(ps.rarefied_vsl, bray_16S_vsl, color="outcome", shape = 'time', title = 'Bray-Curtis, PCoA permANOVA') +
  geom_point(size =4) + theme(aspect.ratio=1) +
  theme_classic() +
  scale_colour_manual(values = c("#28658f", "#ef857c"))
dev.off()

###### otu table for metabolomcs with microbiome analyst. the name of each asv is replaced by its representative fasta sequence

otu_t4_vsl <- data.frame(otu_table(ps_16S_T4_vsl))

ref_seq <- read.fasta('./denovo_unoise_sv/otus.fasta',as.string = TRUE,
                      forceDNAtolower = FALSE)

rownames(otu_t4_vsl) <- unname(unlist(ref_seq[rownames(otu_t4_vsl)]))

write.table(otu_t4_vsl, './metabolomica/otu_t4_vsl.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)

####
# DESEQ vsl3 vs placebo al T4 

sample_data(ps_16S_T4)$treatment <- as.factor(sample_data(ps_16S_T4)$treatment) 
# convertire i valori in sample data in factor

# Convert the phyloseq object to a DESeqDataSet and run DESeq2:
ds_T4 = phyloseq_to_deseq2(ps_16S_T4, ~ treatment)
#far girare con otu table non rarefatta
ds_T4 = DESeq(ds_T4)

alpha = 0.05
res = results(ds_T4, contrast=c("treatment", "VSL3", "placebo"), alpha=alpha)
# order the results table based on the adjusted p-value from the smaller to the bigger
res = res[order(res$padj, na.last=NA), ]
res = cbind(as(res, "data.frame"), as(tax_table(ps_16S_T4)[rownames(res), ], "matrix"))
head(res)
write.table(res, "./R_analysis_all_run1/deseq2_res_T4.txt", sep="\t")

#res <- read.table("deseq2_res.txt", header=TRUE)

pdf('./R_analysis_all_run1/volcano_differential_ASV_T4_VSL_vs_placebo.pdf', height = 7, width = 7)
res2 <- res
res2$diff <- "NO"
res2$diff[res2$log2FoldChange > 1 & res2$padj<0.05] <- "VSL3"
res2$diff[res2$log2FoldChange < (-1) & res2$padj<0.05] <- "placebo"
res2$diff <- factor(res2$diff)

enriched_denovo <- rownames(res2[res2$diff != 'NO',])
ref_seq <- read.fasta('./denovo_unoise_sv/otus.fasta',as.string = TRUE,
                      forceDNAtolower = FALSE)

write.fasta(sequences = as.list(unname(unlist(ref_seq[enriched_denovo]))), as.string = TRUE, names = names(unlist(ref_seq[enriched_denovo])),file.out = './denovo_unoise_sv/denovo_enriched_t4_vsl_placebo.fasta')

rownames(res2) <- gsub('DENOVO', 'ASV_',rownames(res2))
res2$asv_names <- ifelse(is.na(res2$Genus),rownames(res2),paste0(res2$Genus, '_', rownames(res2)))

ggplot2::ggplot(data=res2, aes(x=log2FoldChange, y=-log10(pvalue), fill=diff)) +
  geom_point(size=6, shape=21) + 
  scale_fill_manual(values=c(placebo="#28658f", not="grey",VSL3="#ef857c")) +
  geom_point(data = subset(res2, diff != 'NO'), col = "black", shape=21, size=6) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(padj<0.01 & abs(log2FoldChange)>1,asv_names,'')), size =4)  +
  geom_vline(xintercept=c(-1, 1), col="red", linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  xlim(c(-30,30)) +
  theme(legend.title = element_blank())

dev.off()

### correlation between the ASV significantly enriched in placebo vs vsl at t4 and the immune system 

sist_imm <- data.frame(readxl::read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/delong_covid/sist_immunitario/values.xlsx'))

#intersect(sist_imm$sample,sample_data(ps_16S)$sample_id)

sist_imm$treatment <- unlist(lapply(as.list(sist_imm$sample), function(x) unique(sample_data(ps_16S)[sample_data(ps_16S)$sample_id == x, 'treatment'])))

sist_imm_t4 <- sist_imm[sist_imm$time == 'T4',]

otu_rar <- otu_table(ps.rarefied_T4)
rownames(otu_rar) <- gsub('DENOVO', 'ASV_',rownames(otu_rar))
otu_rar <- otu_rar[rownames(res2),]
rownames(otu_rar) <- res2$asv_names

counts <- data.frame(t(otu_rar[res2[res2$diff != 'NO','asv_names'],]))
rownames(counts) <- sample_data(ps.rarefied_T4)[rownames(counts),'sample_id']$sample_id
counts <- counts[unique(sist_imm_t4$sample),]

sist_imm_t4$value <- as.numeric(sist_imm_t4$value)
sist_imm_t4_new <- data.frame(sample = unique(sist_imm_t4$sample))
for (c in unique(sist_imm_t4$molecule)) {
  sist_imm_t4_new <- cbind(sist_imm_t4_new, sist_imm_t4[sist_imm_t4$molecule == c, 'value'])
}
colnames(sist_imm_t4_new) <- c('sample',unique(sist_imm_t4$molecule))

sist_imm_t4_new <- sist_imm_t4_new %>% drop_na()
counts <- counts[unique(sist_imm_t4_new$sample),]
rownames(sist_imm_t4_new) <- sist_imm_t4_new[,1]
sist_imm_t4_new[,1] <- NULL

ct1 = corr.test(sist_imm_t4_new, counts,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p <- ifelse(clinic_p <= 0.05, '*','')

distance = dist(clinic_r, method = "euclidean")
rowclus <- hclust(dist(clinic_r, method = "euclidean"))
colclus <- hclust(dist(t(clinic_r), method = "euclidean"))

my_palette <- colorRampPalette(brewer.pal(10, "RdYlBu"))(n = 256)

ann <- res2[res2$diff != 'NO', c('diff','asv_names')]
colnames(clinic_r)[colclus$order][colnames(clinic_r)[colclus$order] == 'Clostridium.sensu.stricto_ASV_46'] <- 'Clostridium sensu stricto_ASV_46'
colnames(clinic_r)[colclus$order][colnames(clinic_r)[colclus$order] == 'Clostridium.XlVa_ASV_182'] <- 'Clostridium XlVa_ASV_182'
ann <- ann[match(colnames(clinic_r)[colclus$order], ann$asv_names),]
rownames(ann) <- ann$asv_names
ann$asv_names <- NULL
colnames(ann) <- 'Group'

ann_colors = list(
  Group = c(VSL3="#ef857c", placebo="#28658f"))

pdf('../sist_immunitario/heatmap_correlation_T4_citochine_ASV_enriched.pdf', height = 8,width = 7)
pheatmap(t(clinic_r), cluster_rows = TRUE, cluster_cols = TRUE, color = my_palette, cellwidth = 15, cellheight = 15, angle_col = 315, 
         annotation_row  = ann, annotation_colors = ann_colors, border_color = 'white', display_numbers = t(clinic_p), fontsize_number = 15, number_color = 'black')
dev.off()

### correlation between the ASV significantly enriched in placebo vs vsl at t4 and metabolites

association <- data.frame(name = c('acetic acid',
                                   'propionic acid',
                                   'isobutyric acid',
                                   'butyric acid',
                                   'isovaleric acid',
                                   'valeric acid',
                                   'hexanoic acid',
                                   'dopamine',
                                   'serotonine'), compounds = c('C00033',
                                                                'C00163',
                                                                'C02632',
                                                                'C00246',
                                                                'C08262',
                                                                'C00803',
                                                                'C01585',
                                                                'C03758',
                                                                'C00780'))

metabolomic <- read.table('./metabolomica/values_t4 copia.txt', sep = '\t', header = TRUE)
rownames(metabolomic) <- metabolomic$Compound
metabolomic$Compound <- NULL
rownames(metabolomic) <- unlist(lapply(as.list(rownames(metabolomic)), function(x) association[association$compounds ==x, 1]))

metabolomic <- t(metabolomic)
metabolomic <- data.frame(metabolomic)
metabolomic <- metabolomic[rownames(metabolomic) %in% rownames(sample_data(ps_16S)[sample_data(ps_16S)$time =='T4',]),]

otu_rar <- otu_table(ps.rarefied_T4)
rownames(otu_rar) <- gsub('DENOVO', 'ASV_',rownames(otu_rar))
otu_rar <- otu_rar[rownames(res2),]
rownames(otu_rar) <- res2$asv_names

counts <- data.frame(t(otu_rar[res2[res2$diff != 'NO','asv_names'],]))
counts <- counts[rownames(metabolomic) ,]

ct1 = corr.test(metabolomic, counts,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p <- ifelse(clinic_p <= 0.05, '*','')

distance = dist(clinic_r, method = "euclidean")
rowclus <- hclust(dist(clinic_r, method = "euclidean"))
colclus <- hclust(dist(t(clinic_r), method = "euclidean"))

my_palette <- colorRampPalette(brewer.pal(10, "RdYlBu"))(n = 256)

ann <- res2[res2$diff != 'NO', c('diff','asv_names')]
colnames(clinic_r)[colclus$order][colnames(clinic_r)[colclus$order] == 'Clostridium.sensu.stricto_ASV_46'] <- 'Clostridium sensu stricto_ASV_46'
colnames(clinic_r)[colclus$order][colnames(clinic_r)[colclus$order] == 'Clostridium.XlVa_ASV_182'] <- 'Clostridium XlVa_ASV_182'
ann <- ann[match(colnames(clinic_r)[colclus$order], ann$asv_names),]
rownames(ann) <- ann$asv_names
ann$asv_names <- NULL
colnames(ann) <- 'Treatment'

ann_colors = list(
  Treatment = c(VSL3="#ef857c", placebo="#28658f"))

pdf('./metabolomica/heatmap_correlation_T4_treatment_metabolites_ASV_enriched.pdf', height = 5,width = 10)
pheatmap(clinic_r, cluster_rows = TRUE, cluster_cols = TRUE, color = my_palette, cellwidth = 15, cellheight = 15, angle_col = 315, 
         annotation_col  = ann, annotation_colors = ann_colors, border_color = 'white', display_numbers = clinic_p, fontsize_number = 15, number_color = 'black')

dev.off()

                                       
# comparison between the -log2 of the ratio of the sum of SCFA in placebo at t4 and the sum of the SCFA in placebo at t0 and
# the -log2 of the ratio of the sum of SCFA in vsl at t4 and the sum of the SCFA in vsl at t0

metabolomic_t4 <- data.frame(readxl::read_excel('./metabolomica/all_values.xlsx',sheet = 11))
metabolomic_t0 <- data.frame(readxl::read_excel('./metabolomica/all_values.xlsx',sheet = 10))
rownames(metabolomic_t0) <- metabolomic_t0$X.NAME
metabolomic_t0$X.NAME <- NULL
rownames(metabolomic_t4) <- metabolomic_t4$X.NAME
metabolomic_t4$X.NAME <- NULL

all_values_t0 <- data.frame(sample=rep(colnames(metabolomic_t0), each=nrow(metabolomic_t0)),
                            molecule = rep(rownames(metabolomic_t0), ncol(metabolomic_t0)),
                         time = 'T0',
                         values = as.vector(as.matrix(metabolomic_t0)))

all_values_t0$patient <- unlist(lapply(as.list(all_values_t0$sample), function(x) {
  sample_data(ps_16S)[x, 'sample_id']
}))

all_values_t4 <- data.frame(sample=rep(colnames(metabolomic_t4), each=nrow(metabolomic_t4)),
                            molecule = rep(rownames(metabolomic_t4), ncol(metabolomic_t4)),
                            time = 'T4',
                            values = as.vector(as.matrix(metabolomic_t4)))
all_values_t4$patient <- unlist(lapply(as.list(all_values_t4$sample), function(x) {
  sample_data(ps_16S)[x, 'sample_id']
}))

all <- rbind(all_values_t0,all_values_t4)

all <- all[-which(all$molecule == 'hexanoic acid'),]
all <- all[-which(all$molecule == 'valeric acid'),]

all[which(is.na(all$values)), 'values'] <- 0

for (p in unique(all$patient)) {
  for (m in unique(all[all$patient==p, 'molecule'])) {
    if (all[all$patient == p & all$molecule == m & all$time =='T0','values'] == 0) {
      all[all$patient == p & all$molecule == m & all$time =='T0','values'] <- mean(all[all$molecule == m & all$time == 'T0', 'values'])
    }
    if (all[all$patient == p & all$molecule == m & all$time =='T4','values'] == 0) {
      all[all$patient == p & all$molecule == m & all$time =='T4','values'] <- mean(all[all$molecule == m & all$time == 'T4', 'values'])
    }
  }
}

all_scfa <- all[!all$molecule %in% c('dopamine', 'serotonine'),]

sum_SCFA <- matrix(ncol=2)

for (p in unique(all_scfa$patient)) {
  sum_T0 <- sum(all_scfa[all_scfa$patient == p & all_scfa$time =='T0', 'values'])
  sum_T4 <- sum(all_scfa[all_scfa$patient == p & all_scfa$time =='T4', 'values'])
  transf <- -log2(sum_T4/sum_T0)
  sum_SCFA <- rbind(sum_SCFA, c(p, transf))
  }

sum_SCFA <- data.frame(sum_SCFA[-1,])
colnames(sum_SCFA) <- c('patient', 'log2')

sum_SCFA$treatment <- unlist(lapply(as.list(sum_SCFA$patient), function(x) {
  unique(sample_data(ps_16S)[sample_data(ps_16S)$sample_id == x, 'treatment'])
}))

wilcox.test(as.numeric(sum_SCFA[sum_SCFA$treatment=='VSL3','log2']), as.numeric(sum_SCFA[sum_SCFA$treatment=='placebo','log2']))
# p = 0.79

# comparison between the -log2 of the ratio of the sum of SCFA in vsl responders at t4 and the sum of the SCFA in vsl responders at t0
# and the -log2 of the ratio of the sum of SCFA in vsl non responders at t4 and the sum of the SCFA in vsl non responders at t0

sum_SCFA$outcome <- unlist(lapply(as.list(sum_SCFA$patient), function(x) {
  unique(sample_data(ps_16S)[sample_data(ps_16S)$sample_id == x, 'outcome'])
}))

sum_SCFA_treated <- sum_SCFA[sum_SCFA$treatment =='VSL3',]

wilcox.test(as.numeric(sum_SCFA_treated[sum_SCFA_treated$outcome=='Responder','log2']), as.numeric(sum_SCFA_treated[sum_SCFA_treated$outcome=='Non responder','log2']))
# p = 0.44

## for each metabolite, comparison between the -log2 of the ratio the quantity of that metabolite in vsl at t4 and in vsl at t0 and 
# the -log2 of the ratio the quantity of that metabolite in placebo at t4 and in placebo at t0 
mat <- matrix(ncol=3)

for (p in unique(all$patient)) {
  for (m in unique(all[all$patient == p, 'molecule' ])) {
    transf = -log2(all[all$patient == p & all$molecule == m & all$time =='T4', 'values' ]/all[all$patient == p & all$molecule == m & all$time =='T0', 'values'])
    mat <- rbind(mat, c(p, m, transf))
  }
}

mat <- data.frame(mat[-1,])
colnames(mat) <- c('patient', 'molecule', 'variable')

mat$treatment <- unlist(lapply(as.list(mat$patient), function(x) {
  unique(sample_data(ps_16S)[sample_data(ps_16S)$sample_id == x, 'treatment'])
}))

mat$variable <- as.numeric(mat$variable)

stats <- mat %>% group_by(molecule) %>%
  wilcox_test(variable ~ treatment)

## for each metabolite, comparison between the -log2 of the ratio the quantity of that metabolite in vsl responders at t4 and in vsl responders at t0 and 
# the -log2 of the ratio the quantity of that metabolite in vsl non responders at t4 and in vsl non responders at t0 
mat_treated <- mat[mat$treatment =='VSL3',]
mat_treated$outcome <- unlist(lapply(as.list(mat_treated$patient), function(x) {
  unique(sample_data(ps_16S)[sample_data(ps_16S)$sample_id == x, 'outcome'])
}))

stats_treated <- mat_treated %>% group_by(molecule) %>%
  wilcox_test(variable ~ outcome)

####

# DESEQ responders e non responders al T4 

sample_data(ps_16S_T4_vsl)$outcome <- as.factor(sample_data(ps_16S_T4_vsl)$outcome) 
# convertire i valori in sample data in factor

# Convert the phyloseq object to a DESeqDataSet and run DESeq2:
ds_T4_vsl = phyloseq_to_deseq2(ps_16S_T4_vsl, ~ outcome)
#far girare con otu table non rarefatta
ds_T4_vsl = DESeq(ds_T4_vsl)

alpha = 0.05
res = results(ds_T4_vsl, contrast=c("outcome", "Responder", "Non responder"), alpha=alpha)
# order the results table based on the adjusted p-value from the smaller to the bigger
res = res[order(res$padj, na.last=NA), ]
res = cbind(as(res, "data.frame"), as(tax_table(ps_16S_T4_vsl)[rownames(res), ], "matrix"))
head(res)
write.table(res, "./R_analysis_all_run1/deseq2_res_T4_vsl.txt", sep="\t")

#res <- read.table("deseq2_res.txt", header=TRUE)

pdf('./R_analysis_all_run1/volcano_differential_ASV_T4_VSL_responders_vs_non.pdf', height = 7, width = 7)
res2 <- res
res2$diff <- "NO"
res2$diff[res2$log2FoldChange > 1 & res2$padj<0.05] <- "Responder"
res2$diff[res2$log2FoldChange < (-1) & res2$padj<0.05] <- "Non responder"

enriched_denovo <- rownames(res2[res2$diff != 'NO',])
ref_seq <- read.fasta('./denovo_unoise_sv/otus.fasta',as.string = TRUE,
                      forceDNAtolower = FALSE)

write.fasta(sequences = as.list(unname(unlist(ref_seq[enriched_denovo]))), as.string = TRUE, names = names(unlist(ref_seq[enriched_denovo])),file.out = './denovo_unoise_sv/denovo_enriched_t4_vsl_resp_non_resp.fasta')


res2$diff <- factor(res2$diff)
rownames(res2) <- gsub('DENOVO', 'ASV_',rownames(res2))
res2$asv_names <- ifelse(is.na(res2$Genus),rownames(res2),paste0(res2$Genus, '_', rownames(res2)))

ggplot(data=res2, aes(x=log2FoldChange, y=-log10(pvalue), fill=diff)) +
  geom_point(size=6, shape=21, stroke =0.5) +
  scale_fill_manual(values=c("grey60","#ffb703","#fb8500"), labels = c('not','Non responder','Responder')) +
  geom_point(data = subset(res2, diff != 'NO'), col = "black", shape=21, size=6) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(padj<0.01 & abs(log2FoldChange)>1,asv_names,'')), size =4, max.overlaps = 100)  +
  geom_vline(xintercept=c(-1, 1), col="red", linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  xlim(c(-30,30)) +
  theme(legend.title = element_blank())

dev.off()

### correlation between the ASV significantly enriched in vsl responder vs vsl non responders at t4 and the immune system 

sist_imm <- data.frame(readxl::read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/delong_covid/sist_immunitario/values.xlsx'))
#intersect(sist_imm$sample,sample_data(ps_16S)$sample_id)

sist_imm$treatment <- unlist(lapply(as.list(sist_imm$sample), function(x) unique(sample_data(ps_16S)[sample_data(ps_16S)$sample_id == x, 'treatment'])))

sist_imm_t4_vsl <- sist_imm[sist_imm$time == 'T4' & sist_imm$treatment == 'VSL3',]

otu_rar <- otu_table(ps.rarefied_T4_vsl)
rownames(otu_rar) <- gsub('DENOVO', 'ASV_',rownames(otu_rar))
otu_rar <- otu_rar[rownames(res2),]
rownames(otu_rar) <- res2$asv_names

counts <- data.frame(t(otu_rar[res2[res2$diff != 'NO','asv_names'],]))
rownames(counts) <- sample_data(ps.rarefied_T4)[rownames(counts),'sample_id']$sample_id
counts <- counts[unique(sist_imm_t4_vsl$sample),]

sist_imm_t4_vsl$value <- as.numeric(sist_imm_t4_vsl$value)
sist_imm_t4_vsl_new <- data.frame(sample = unique(sist_imm_t4_vsl$sample))
for (c in unique(sist_imm_t4_vsl$molecule)) {
  sist_imm_t4_vsl_new <- cbind(sist_imm_t4_vsl_new, sist_imm_t4_vsl[sist_imm_t4_vsl$molecule == c, 'value'])
}
colnames(sist_imm_t4_vsl_new) <- c('sample',unique(sist_imm_t4_vsl$molecule))

sist_imm_t4_vsl_new <- sist_imm_t4_vsl_new %>% drop_na()
counts <- counts[unique(sist_imm_t4_vsl_new$sample),]
rownames(sist_imm_t4_vsl_new) <- sist_imm_t4_vsl_new[,1]
sist_imm_t4_vsl_new[,1] <- NULL

ct1 = corr.test(sist_imm_t4_vsl_new, counts,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p <- ifelse(clinic_p <= 0.05, '*','')

distance = dist(clinic_r, method = "euclidean")
rowclus <- hclust(dist(clinic_r, method = "euclidean"))
colclus <- hclust(dist(t(clinic_r), method = "euclidean"))

my_palette <- colorRampPalette(brewer.pal(10, "RdYlBu"))(n = 256)

ann <- res2[res2$diff != 'NO', c('diff','asv_names')]
colnames(clinic_r)[colclus$order][colnames(clinic_r)[colclus$order] == 'Clostridium.sensu.stricto_ASV_2374'] <- 'Clostridium sensu stricto_ASV_2374'
ann <- ann[match(colnames(clinic_r)[colclus$order], ann$asv_names),]
rownames(ann) <- ann$asv_names
ann$asv_names <- NULL
colnames(ann) <- 'Outcome'

ann_colors = list(
  Outcome = c("Responder"="#fb8500", "Non responder"="#ffb703"))

pdf('../sist_immunitario/heatmap_correlation_T4_outcome_citochine_ASV_enriched.pdf', height = 5,width = 7)
pheatmap(t(clinic_r), cluster_rows = TRUE, cluster_cols = TRUE, color = my_palette, cellwidth = 15, cellheight = 15, angle_col = 315, 
         annotation_row  = ann, annotation_colors = ann_colors, border_color = 'white', display_numbers = t(clinic_p), fontsize_number = 15, number_color = 'black')

dev.off()

### correlation between the ASV significantly enriched in vsl resp and vsl non resp at t4 and metabolites

metabolomic <- read.table('./metabolomica/values_t4 copia.txt', sep = '\t', header = TRUE)
rownames(metabolomic) <- metabolomic$Compound
metabolomic$Compound <- NULL
rownames(metabolomic) <- unlist(lapply(as.list(rownames(metabolomic)), function(x) association[association$compounds ==x, 1]))

metabolomic <- t(metabolomic)
metabolomic <- data.frame(metabolomic)
metabolomic <- metabolomic[rownames(metabolomic) %in% rownames(sample_data(ps_16S)[sample_data(ps_16S)$treatment =='VSL3' & sample_data(ps_16S)$time =='T4',]),]

otu_rar <- otu_table(ps.rarefied_T4_vsl)
rownames(otu_rar) <- gsub('DENOVO', 'ASV_',rownames(otu_rar))
otu_rar <- otu_rar[rownames(res2),]
rownames(otu_rar) <- res2$asv_names

counts <- data.frame(t(otu_rar[res2[res2$diff != 'NO','asv_names'],]))
counts <- counts[rownames(metabolomic) ,]
counts <- counts[,-which(colSums(counts) == 0)]

ct1 = corr.test(metabolomic, counts,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p <- ifelse(clinic_p <= 0.05, '*','')

distance = dist(clinic_r, method = "euclidean")
rowclus <- hclust(dist(clinic_r, method = "euclidean"))
colclus <- hclust(dist(t(clinic_r), method = "euclidean"))

my_palette <- colorRampPalette(brewer.pal(10, "RdYlBu"))(n = 256)

ann <- res2[res2$diff != 'NO', c('diff','asv_names')]
colnames(clinic_r)[colclus$order][colnames(clinic_r)[colclus$order] == 'Clostridium.sensu.stricto_ASV_2374'] <- 'Clostridium sensu stricto_ASV_2374'
ann <- ann[match(colnames(clinic_r)[colclus$order], ann$asv_names),]
rownames(ann) <- ann$asv_names
ann$asv_names <- NULL
colnames(ann) <- 'Outcome'

ann_colors = list(
  Outcome = c("Responder"="#fb8500", "Non responder"="#ffb703"))

pdf('./metabolomica/heatmap_correlation_T4_outcome_metabolites_ASV_enriched.pdf', height = 5,width = 7)
pheatmap(clinic_r, cluster_rows = TRUE, cluster_cols = TRUE, color = my_palette, cellwidth = 15, cellheight = 15, angle_col = 315, 
         annotation_col  = ann, annotation_colors = ann_colors, border_color = 'white', display_numbers = clinic_p, fontsize_number = 15, number_color = 'black')

dev.off()

######
# predict responders and not at t4

metadata_T4_VSL <- sample_data(ps.rarefied_T4_vsl)
otu_table_T4_VSL <- data.frame(otu_table(ps.rarefied_T4_vsl))
#apply(otu_table_T4_VSL,2,sum)
dim(otu_table_T4_VSL)
# 4072   23
dim(metadata_T4_VSL)
# 23  5

otu_nonzero_counts <- apply(otu_table_T4_VSL, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  # ceiling returns the smallest integer value, which is greater than or equal the number you give in input (in thid case
  # the input is 8.8 and the function returns 9)
  for ( i in 1:nrow(table) ) {
    # for each genus identify the number of columns (samples) in which the genus has a count major than 0
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    # if the number of columns (samples) in which the genus has a count major than 0 is higher than the set cut off
    # keep that genus
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep ,])
}

# Based on the above histogram I removed OTUs that have non-zero values in <= 20% of samples:
otu_table_rare_removed <- remove_rare(table=otu_table_T4_VSL, cutoff_pro=0.2)
dim(otu_table_rare_removed)
# 1110   23

otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)

# Running Model

otu_table_scaled_T4_vsl_outcome <- data.frame(t(otu_table_scaled))  
# add a new column "variable" to the scaled and centered otu table reporting the statud (HS/CBP) of the samples (rows),
# this info is extracted from metadata object (imported before from a txt file)
otu_table_scaled_T4_vsl_outcome$outcome <- factor(metadata_T4_VSL[rownames(otu_table_scaled_T4_vsl_outcome), "outcome"]$outcome)

set.seed(88)

RF_state_classify <- randomForest(x=otu_table_scaled_T4_vsl_outcome[,1:(ncol(otu_table_scaled_T4_vsl_outcome)-1)], y=otu_table_scaled_T4_vsl_outcome[ , ncol(otu_table_scaled_T4_vsl_outcome)] , ntree=501, importance=TRUE, proximities=TRUE )

# Assessing Model Fit

RF_state_classify
#ntree = 501, importance = TRUE, proximities = TRUE) 
# Type of random forest: classification
# Number of trees: 501
# No. of variables tried at each split: 33
# 
# OOB estimate of  error rate: 17.39%
# Confusion matrix:
#   Non_responder Responder class.error
# Non_responder             0         4           1
# Responder                 0        19           0

# Permutation Test
RF_state_classify_sig <- rf.significance( x=RF_state_classify ,  xdata=otu_table_scaled_T4_vsl_outcome[,1:(ncol(otu_table_scaled_T4_vsl_outcome)-1)] , nperm=1000 , ntree=501 )  
# Number of permutations:  1000 
# p-value:  0.024 
# Model signifiant at p = 0.024 
# Model OOB error:  0.173913 
# Random OOB error:  0.173913 
# min random global error: 0.08695652 
# max random global error:  0.3478261 
# min random within class error: NA 
# max random within class error:  NA

#Accuracy Estimated by Cross-validation
fit_control <- trainControl(method = "LOOCV",)    
RF_state_classify_loocv <- train(otu_table_scaled_T4_vsl_outcome[,1:(ncol(otu_table_scaled_T4_vsl_outcome)-1)], y=otu_table_scaled_T4_vsl_outcome[, ncol(otu_table_scaled_T4_vsl_outcome)], 
                                 method="rf", ntree=501 , tuneGrid=data.frame( mtry=25), trControl=fit_control)
RF_state_classify_loocv$

#mtry Accuracy Kappa
# 1   25 0.826087     0

control <- trainControl(method="cv", classProbs = TRUE)
set.seed(55)
fit <- train(otu_table_scaled_T4_vsl_outcome[,1:(ncol(otu_table_scaled_T4_vsl_outcome)-1)], y=otu_table_scaled_T4_vsl_outcome[, ncol(otu_table_scaled_T4_vsl_outcome)], 
             method="rf",metric="ROC", trControl=control, ntree=501 , tuneGrid=data.frame(mtry=25))
print(fit)
#  Accuracy   Kappa    
# 0.85  0

#Identifying Important Features
RF_outcome_classify_imp <- as.data.frame( RF_state_classify$importance )
# add a column "feature" reporting the genus of each row (that is each feature of the model, so each genus)
RF_outcome_classify_imp$features <- rownames( RF_outcome_classify_imp )
# reorder the rows of the object on the base of decreasing values of the variable "MeanDecreaseGini"
RF_outcome_classify_imp_sorted <- arrange( RF_outcome_classify_imp  , desc(MeanDecreaseGini)  )

otu_to_taxa <- read.table('./denovo_unoise_sv/taxa_split.txt',fill = , sep = '\t', row.names = 1)
otu_to_taxa[otu_to_taxa == ''] <- NA

v <- unlist(lapply(seq_along(1:nrow(otu_to_taxa)), function(x) {otu_to_taxa[x,6]}))
names(v) <- rownames(otu_to_taxa)    

# order the genus for increasing value of gini coefficient. convert the features variable in factor
zzz <- RF_outcome_classify_imp[order(RF_outcome_classify_imp$MeanDecreaseGini, decreasing = TRUE),]
zzz$features <- ifelse(is.na(v[zzz$features]), gsub('DENOVO','ASV_',zzz$features), paste0(v[zzz$features], '_',gsub('DENOVO','ASV_',zzz$features)))
zzz2 <- zzz[10:1,]

pdf('./R_analysis_all_run1/genus_gini_coef_T4_predict_outcome.pdf', height = 2, width = 3.5)
# plot the genus in order of increasing gini coefficient
zzz2$name2 <- 1:nrow(zzz2)
ggplot(zzz2, aes(x=factor(name2), y=MeanDecreaseGini)) +
  geom_segment( aes(x=factor(name2), xend=factor(name2), y=0, yend=MeanDecreaseGini), color="black") +
  geom_point( color="orange", size=3, alpha=0.7) +
  theme_bw() +
  scale_x_discrete(labels=zzz2$features) +
  coord_flip() +
  ggtitle('OOB=17.39%;p=0.024;\nAccuracy=0.85;Kappa=0')+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.title.x = element_text(size=7),
    title = element_text(size=7))
dev.off()


####
metadata_T4_VSL_responders <- rownames(metadata_T4_VSL[metadata_T4_VSL$outcome=='Responder',])
metadata_T4_VSL_non_responders <- rownames(metadata_T4_VSL[metadata_T4_VSL$outcome=='Non_responder',])

otu_to_taxa <- read.table('./denovo_unoise_sv/taxa_split.txt',fill = , sep = '\t', row.names = 1)
otu_to_taxa[otu_to_taxa == ''] <- NA

v <- unlist(lapply(seq_along(1:nrow(otu_to_taxa)), function(x) {otu_to_taxa[x,6]}))
names(v) <- rownames(otu_to_taxa)    

df <- data.frame(samples = rep(c(metadata_T4_VSL_responders,metadata_T4_VSL_non_responders), each=10),
                 genus = rep(RF_outcome_classify_imp_sorted$features[1:10],length(c(metadata_T4_VSL_responders,metadata_T4_VSL_non_responders))),
                 outcome=c(rep('Responder', length(metadata_T4_VSL_responders)*10), rep('Non responder', length(metadata_T4_VSL_non_responders)*10)))

df$abundance <- NA

for (r in 1:nrow(df)) {
  df$abundance[r] <- otu_table_T4_VSL[df[r,2],df[r,1]]
}

df$outcome <- as.factor(df$outcome)
df$genus <- ifelse(is.na(v[df$genus]), gsub('DENOVO','ASV_',df$genus), paste0(v[df$genus], '_',gsub('DENOVO','ASV_',df$genus)))
df$genus <- as.factor(df$genus)

stat.test <-  df %>%
  group_by(genus) %>%
  rstatix::wilcox_test(abundance~outcome) %>% 
  adjust_pvalue(method = 'fdr') %>%
  add_significance("p.adj") 

p_top10 <- round(stat.test$p,3)
names(p_top10) <- stat.test$genus

stat.test_less <-  df %>%
  group_by(genus) %>%
  rstatix::wilcox_test(abundance~outcome, alternative = 'less') %>% 
  adjust_pvalue(method = 'fdr') %>%
  add_significance("p.adj") 

effect_size_top10 <- qnorm(stat.test_less$p, lower.tail = FALSE)/sqrt(stat.test_less$n1+stat.test_less$n2)
names(effect_size_top10) <- stat.test_less$genus

df <- data.frame(effsize = effect_size_top10, p = p_top10, asv = names(effect_size_top10))
df$asv <- factor(df$asv, levels=df[order(df$effsize, decreasing = TRUE),'asv'])

df$outcome <- "Non responder"
df$outcome[df$effsize > 0] <- "Responder"
df$outcome <- factor(df$outcome, levels=c("Responder","Non responder"))

pdf('./R_analysis_all_run1/effect_size_responders_nr_T4.pdf', width = 5, height = 3)
ggplot(df, aes(y=effsize,x=asv, fill=outcome)) +
  geom_bar(stat = "identity", width = 0.7) +  coord_flip() +
  geom_text(aes(label = p), position=position_dodge(width=0.9),size=3,
            hjust = ifelse(df$effsize >= 0, 1.2, -0.2)) +
  scale_fill_manual(values = c("#fb8500","#ffb703")) +
  theme_classic() +
  theme(axis.title.y = element_blank(), legend.title = element_blank(), legend.position = 'top') +
  ylab('Effect size') +
  scale_x_discrete(limits = rev(levels(df$asv)))
dev.off()

###AUROC
# AUROC to understand if the created model classifies the samples in a good way
# Load required packages

# Get the true class labels
true_labels <- metadata_T4_VSL$outcome
predicted_probs <- RF_state_classify$votes[, 2]
# Calculate AUROC
roc_obj <- roc(true_labels, predicted_probs)
auc_value <- auc(roc_obj)
# auc function computes the numeric value of area under the ROC curve (
# Plot ROC curve
pdf('./R_analysis_all_run1/ROC_T4_vsl_predict_outcome.pdf')
plot(roc_obj, main = "ROC Curve", print.auc = TRUE,)
# Print AUROC value
cat("AUROC:", auc_value, "\n")
dev.off()

######
# predict responders and not at t0

metadata_T0_VSL <- sample_data(ps.rarefied_T0_vsl)
otu_table_T0_VSL <- data.frame(otu_table(ps.rarefied_T0_vsl))
apply(otu_table_T0_VSL, 2, sum)
dim(otu_table_T0_VSL)
# 4072   23
dim(metadata_T0_VSL)
# 23  5

otu_nonzero_counts <- apply(otu_table_T0_VSL, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  # ceiling returns the smallest integer value, which is greater than or equal the number you give in input (in thid case
  # the input is 8.8 and the function returns 9)
  for ( i in 1:nrow(table) ) {
    # for each genus identify the number of columns (samples) in which the genus has a count major than 0
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    # if the number of columns (samples) in which the genus has a count major than 0 is higher than the set cut off
    # keep that genus
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F])
}

# Based on the above histogram I removed OTUs that have non-zero values in <= 20% of samples:
otu_table_rare_removed <- remove_rare(table=otu_table_T0_VSL, cutoff_pro=0.2)
dim(otu_table_rare_removed)
#1189   23

otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
#apply(otu_table_rare_removed_norm,2,sum)

# One approach is to standardize the data by subtracting each sample's mean (center) and then dividing by the 
# sample's standard deviation (scale). In other words, each value is converted into a Z-score.
otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)

# Running Model

otu_table_scaled_T0_vsl_outcome <- data.frame(t(otu_table_scaled))  
# add a new column "variable" to the scaled and centered otu table reporting the statud (HS/CBP) of the samples (rows),
# this info is extracted from metadata object (imported before from a txt file)
otu_table_scaled_T0_vsl_outcome$outcome <- factor(metadata_T0_VSL[rownames(otu_table_scaled_T0_vsl_outcome),]$outcome)

set.seed(88)

RF_state_classify <- randomForest(x=otu_table_scaled_T0_vsl_outcome[,1:(ncol(otu_table_scaled_T0_vsl_outcome)-1)], y=otu_table_scaled_T0_vsl_outcome[ , ncol(otu_table_scaled_T0_vsl_outcome)] , ntree=501, importance=TRUE, proximities=TRUE )

# Assessing Model Fit

RF_state_classify
#ntree = 501, importance = TRUE, proximities = TRUE) 
#  Type of random forest: classification
# Number of trees: 501
# No. of variables tried at each split: 34
# 
# OOB estimate of  error rate: 17.39%
# Confusion matrix:
#   Non_responder Responder class.error
# Non_responder             0         4           1
# Responder                 0        19           0

# Permutation Test
RF_state_classify_sig <- rf.significance( x=RF_state_classify ,  xdata=otu_table_scaled_T0_vsl_outcome[,1:(ncol(otu_table_scaled_T0_vsl_outcome)-1)] , nperm=1000 , ntree=501 )  
#Number of permutations:  1000 
# p-value:  0.01 
# Model signifiant at p = 0.01 
# Model OOB error:  0.173913 
# Random OOB error:  0.173913 
# min random global error: 0.08695652 
# max random global error:  0.3478261 
# min random within class error: NA 
# max random within class error:  NA 

#Accuracy Estimated by Cross-validation
fit_control <- trainControl(method = "LOOCV")    
RF_state_classify_loocv <- train(otu_table_scaled_T0_vsl_outcome[,1:(ncol(otu_table_scaled_T0_vsl_outcome)-1)], y=otu_table_scaled_T0_vsl_outcome[, ncol(otu_table_scaled_T0_vsl_outcome)], 
                                 method="rf", ntree=501 , tuneGrid=data.frame( mtry=25), trControl=fit_control, metric="Accuracy")
RF_state_classify_loocv$results
#mtry Accuracy Kappa
# 1   25 0.826087     0

control <- trainControl(method="cv", classProbs = TRUE)
set.seed(55)
fit <- train(otu_table_scaled_T0_vsl_outcome[,1:(ncol(otu_table_scaled_T0_vsl_outcome)-1)], y=otu_table_scaled_T0_vsl_outcome[, ncol(otu_table_scaled_T0_vsl_outcome)], 
             method="rf", metric="ROC", trControl=control, ntree=501 , tuneGrid=data.frame(mtry=25))
print(fit)
#  Accuracy   Kappa    
# 0.85  0

#Identifying Important Features
RF_outcome_T0_classify_imp <- as.data.frame( RF_state_classify$importance )
# add a column "feature" reporting the genus of each row (that is each feature of the model, so each genus)
RF_outcome_T0_classify_imp$features <- rownames( RF_outcome_T0_classify_imp )
# reorder the rows of the object on the base of decreasing values of the variable "MeanDecreaseGini"
RF_outcome_T0_classify_imp_sorted <- arrange( RF_outcome_T0_classify_imp  , desc(MeanDecreaseGini)  )

otu_to_taxa <- read.table('./denovo_unoise_sv/taxa_split.txt',fill = , sep = '\t', row.names = 1)
otu_to_taxa[otu_to_taxa == ''] <- NA

v <- unlist(lapply(seq_along(1:nrow(otu_to_taxa)), function(x) {otu_to_taxa[x,6]}))
names(v) <- rownames(otu_to_taxa)    

# order the genus for increasing value of gini coefficient. convert the features variable in factor
zzz <- RF_outcome_T0_classify_imp[order(RF_outcome_T0_classify_imp$MeanDecreaseGini, decreasing = TRUE),]
zzz$features <- ifelse(is.na(v[zzz$features]), gsub('DENOVO','ASV_',zzz$features), paste0(v[zzz$features], '_',gsub('DENOVO','ASV_',zzz$features)))
zzz2 <- zzz[10:1,]

pdf('./R_analysis_all_run1/genus_gini_coef_T0_predict_outcome.pdf', height = 2, width = 3.5)
# plot the genus in order of increasing gini coefficient
zzz2$name2 <- 1:nrow(zzz2)
ggplot(zzz2, aes(x=factor(name2), y=MeanDecreaseGini)) +
  geom_segment( aes(x=factor(name2), xend=factor(name2), y=0, yend=MeanDecreaseGini), color="black") +
  geom_point( color="orange", size=3, alpha=0.7) +
  theme_bw() +
  scale_x_discrete(labels=zzz2$features) +
  coord_flip() +
  ggtitle('OOB=17.39%;p=0.01;\nAccuracy=0.85;Kappa=0')+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.title.x = element_text(size=7),
    title = element_text(size=7))
dev.off()

##
metadata_T0_VSL_responders <- rownames(metadata_T0_VSL[metadata_T0_VSL$outcome=='Responder',])
metadata_T0_VSL_non_responders <- rownames(metadata_T0_VSL[metadata_T0_VSL$outcome=='Non_responder',])

otu_to_taxa <- read.table('./denovo_unoise_sv/taxa_split.txt',fill = , sep = '\t', row.names = 1)
otu_to_taxa[otu_to_taxa == ''] <- NA

v <- unlist(lapply(seq_along(1:nrow(otu_to_taxa)), function(x) {otu_to_taxa[x,6]}))
names(v) <- rownames(otu_to_taxa)    

df <- data.frame(samples = rep(c(metadata_T0_VSL_responders,metadata_T0_VSL_non_responders), each=10),
                 genus = rep(RF_outcome_T0_classify_imp_sorted$features[1:10],length(c(metadata_T0_VSL_responders,metadata_T0_VSL_non_responders))),
                 outcome=c(rep('Responder', length(metadata_T0_VSL_responders)*10), rep('Non responder', length(metadata_T0_VSL_non_responders)*10)))

df$abundance <- NA

for (r in 1:nrow(df)) {
  df$abundance[r] <- otu_table_T0_VSL[df[r,2],df[r,1]]
}

df$outcome <- as.factor(df$outcome)
df$genus <- ifelse(is.na(v[df$genus]), gsub('DENOVO','ASV_',df$genus), paste0(v[df$genus], '_',gsub('DENOVO','ASV_',df$genus)))
df$genus <- as.factor(df$genus)

stat.test <-  df %>%
  group_by(genus) %>%
  rstatix::wilcox_test(abundance~outcome) %>% 
  adjust_pvalue(method = 'fdr') %>%
  add_significance("p.adj") 

p_top10 <- round(stat.test$p,3)
names(p_top10) <- stat.test$genus

stat.test_less <-  df %>%
  group_by(genus) %>%
  rstatix::wilcox_test(abundance~outcome, alternative = 'less') %>% 
  adjust_pvalue(method = 'fdr') %>%
  add_significance("p.adj") 

effect_size_top10 <- qnorm(stat.test_less$p, lower.tail = FALSE)/sqrt(stat.test_less$n1+stat.test_less$n2)
names(effect_size_top10) <- stat.test_less$genus

df <- data.frame(effsize = effect_size_top10, p = p_top10, asv = names(effect_size_top10))
df$asv <- factor(df$asv, levels=df[order(df$effsize, decreasing = TRUE),'asv'])

df$outcome <- "Non responder"
df$outcome[df$effsize > 0] <- "Responder"
df$outcome <- factor(df$outcome, levels=c("Responder","Non responder"))

pdf('./R_analysis_all_run1/effect_size_responders_nr_T0.pdf', width = 5, height = 3)
ggplot(df, aes(y=effsize,x=asv,fill=outcome)) +
  geom_bar(stat = "identity", width = 0.7) +  coord_flip() +
  geom_text(aes(label = p), position=position_dodge(width=0.9),size=3,
          hjust = ifelse(df$effsize >= 0, 1.2, -0.2)) +
  theme_classic() +
  scale_fill_manual(values = c("#fb8500","#ffb703")) +
  theme(axis.title.y = element_blank(), legend.title = element_blank(), legend.position = 'top') +
  ylab('Effect size') +
  scale_x_discrete(limits = rev(levels(df$asv)))
dev.off()

###AUROC
# AUROC to understand if the created model classifies the samples in a good way
# Load required packages

# Get the true class labels
true_labels <- metadata_T0_VSL$outcome
# Get the predicted class probabilities for the positive class (the second column
# corresponds to the probability that the sample is HS)
predicted_probs <- RF_state_classify$votes[, 2]
# Calculate AUROC
roc_obj <- roc(true_labels, predicted_probs)
auc_value <- auc(roc_obj)
# auc function computes the numeric value of area under the ROC curve (
# Plot ROC curve
pdf('./R_analysis_all_run1/ROC_T0_vsl_predict_outcome.pdf')
plot(roc_obj, main = "ROC Curve", print.auc = TRUE)
# Print AUROC value
cat("AUROC:", auc_value, "\n")
dev.off()


#### metabolomic

sample <- data.frame(readxl::read_excel('./metabolomica/metadata_t0.xlsx'))
sample$treatment <- sample_data(ps_16S)[sample$X.NAME, 4]$treatment

write.table(sample, file='./metabolomica/metadata_t0.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

t <- data.frame(readxl::read_excel('./metabolomica/taxa.xlsx'))
rownames(t) <- t[,1]
t[,1] <- NULL

fasta_ASV <- read.fasta(file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/delong_covid/analisi_micca_nuovi_vecchi_insieme/denovo_unoise_sv/otus.fasta',as.string = TRUE,
                        forceDNAtolower = FALSE)

unname(unlist(fasta_ASV))
df <- data.frame(ASV_name = names(fasta_ASV), ref = unname(unlist(fasta_ASV)))
rownames(df) <- df[,1]
df[,1] <- NULL

ot <- otu_table(ps_16S)[,sample$X.NAME]
rownames(ot) <- df[rownames(ot),1]
ot2 <- ot[,colSums(ot) != 0]

write.table(ot2,file='./metabolomica/otu_t0.txt', quote = FALSE, sep = '\t')

taxa <- read.table('./denovo_unoise_sv/taxa.txt', sep = '\t')
rownames(taxa) <- taxa[,1]
taxa[,1] <- NULL
rownames(taxa) <- df[rownames(taxa),1]
taxa <- taxa[rownames(ot2),]
colnames(taxa) <- c("Kingdom", "Phylum", 'Class', 'Order', 'Family', 'Genus')
taxa[taxa==''] <- NA
write.table(taxa, file = './metabolomica/taxa_renamed.txt', quote = FALSE, sep = '\t')


# t4 
sample <- data.frame(readxl::read_excel('./metabolomica/metadata_t4.xlsx'))
sample$treatment <- sample_data(ps_16S)[sample$X.NAME, 4]$treatment
sample$outcome <- sample_data(ps_16S)[sample$X.NAME, 5]$outcome

write.table(sample, file='./metabolomica/metadata_t4.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

t <- data.frame(readxl::read_excel('./metabolomica/taxa.xlsx'))
rownames(t) <- t[,1]
t[,1] <- NULL

fasta_ASV <- read.fasta(file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/delong_covid/analisi_micca_nuovi_vecchi_insieme/denovo_unoise_sv/otus.fasta',as.string = TRUE,
                        forceDNAtolower = FALSE)

df <- data.frame(ASV_name = names(fasta_ASV), ref = unname(unlist(fasta_ASV)))
rownames(df) <- df[,1]
df[,1] <- NULL

ot <- otu_table(ps_16S)[,sample$X.NAME]
rownames(ot) <- df[rownames(ot),1]
ot2 <- ot[,colSums(ot) != 0]

write.table(ot2,file='./metabolomica/otu_t4.txt', quote = FALSE, sep = '\t')

taxa <- read.table('./metabolomica/taxa.txt', sep = '\t')
rownames(taxa) <- taxa[,1]
taxa[,1] <- NULL
rownames(taxa) <- df[rownames(taxa),1]
taxa <- taxa[rownames(ot2),]
colnames(taxa) <- c("Kingdom", "Phylum", 'Class', 'Order', 'Family', 'Genus')
write.table(taxa, file = './metabolomica/taxa_t4.txt', quote = FALSE, sep = '\t')

###
metabolomic <- read.table('./metabolomica/values_t4 copia.txt', sep = '\t', header = TRUE)
rownames(metabolomic) <- metabolomic$Compound
metabolomic$Compound <- NULL

df <- data.frame(compound = rep(rownames(metabolomic), ncol(metabolomic)),
                 sample = rep(colnames(metabolomic), each = nrow(metabolomic)),
                 values = as.vector(as.matrix(metabolomic)))

df$treatment <- unlist(lapply(as.list(df$sample), function(x) unique(sample_data(ps_16S)[x, 'treatment'])))
df$outcome <- unlist(lapply(as.list(df$sample), function(x) unique(sample_data(ps_16S)[x, 'outcome'])))

association <- data.frame(name = c('acetic acid',
                      'propionic acid',
                      'isobutyric acid',
                      'butyric acid',
                      'isovaleric acid',
                      'valeric acid',
                      'hexanoic acid',
                      'dopamine',
                      'serotonine'), compounds = c('C00033',
                                                    'C00163',
                                                    'C02632',
                                                    'C00246',
                                                    'C08262',
                                                    'C00803',
                                                    'C01585',
                                                    'C03758',
                                                    'C00780'))

df$compound <- unlist(lapply(as.list(df$compound), function(x) association[association$compounds == x, 1]))

stats_vsl_placebo_w <- df %>%
  group_by(compound) %>%
  wilcox_test(values~treatment) %>%
  adjust_pvalue(method = 'fdr')

stats_vsl_placebo_w$y.position <- unlist(lapply(seq_along(1:length(stats_vsl_placebo_w$compound)), function(x) {
  max(df[df$compound == stats_vsl_placebo_w$compound[x],'values'][! is.na(df[df$compound == stats_vsl_placebo_w$compound[x],'values'])]) +max(df[df$compound == stats_vsl_placebo_w$compound[x],'values'][! is.na(df[df$compound == stats_vsl_placebo_w$compound[x],'values'])])/14
}))

df_vsl <- df[df$treatment == 'VSL3',]
stats_resp_non_w <- df_vsl %>%
  group_by(compound) %>%
  wilcox_test(values~outcome) %>%
  adjust_pvalue(method = 'fdr')

stats_resp_non_w$y.position <- unlist(lapply(seq_along(1:length(stats_resp_non_w$compound)), function(x) {
  max(df_vsl[df_vsl$compound == stats_resp_non_w$compound[x],'values'][! is.na(df_vsl[df_vsl$compound == stats_resp_non_w$compound[x],'values'])]) +max(df_vsl[df_vsl$compound == stats_resp_non_w$compound[x],'values'][! is.na(df_vsl[df_vsl$compound == stats_resp_non_w$compound[x],'values'])])/14
}))

stats_vsl_placebo_t <- df %>%
  group_by(compound) %>%
  t_test(values~treatment) %>%
  adjust_pvalue(method = 'fdr')

stats_vsl_placebo_t$y.position <- unlist(lapply(seq_along(1:length(stats_vsl_placebo_t$compound)), function(x) {
  max(df[df$compound == stats_vsl_placebo_t$compound[x],'values'][! is.na(df[df$compound == stats_vsl_placebo_t$compound[x],'values'])]) +max(df[df$compound == stats_vsl_placebo_t$compound[x],'values'][! is.na(df[df$compound == stats_vsl_placebo_t$compound[x],'values'])])/14
}))

stats_resp_non_t <- df_vsl %>%
  group_by(compound) %>%
  t_test(values~outcome) %>%
  adjust_pvalue(method = 'fdr')

stats_resp_non_t$y.position <- unlist(lapply(seq_along(1:length(stats_resp_non_t$compound)), function(x) {
  max(df_vsl[df_vsl$compound == stats_resp_non_t$compound[x],'values'][! is.na(df_vsl[df_vsl$compound == stats_resp_non_t$compound[x],'values'])]) +max(df_vsl[df_vsl$compound == stats_resp_non_t$compound[x],'values'][! is.na(df_vsl[df_vsl$compound == stats_resp_non_t$compound[x],'values'])])/14
}))

df$treatment <- factor(df$treatment)
df$outcome <- factor(df$outcome)

pdf('./metabolomica/boxplot_metabolomics_t4_vsl_placebo_t_test.pdf', width = 10, height = 10)
ggplot(df, aes(x= treatment, y=values,color=treatment)) +
  geom_boxplot(aes(colour=treatment),width=0.4, lwd=1) +
  scale_colour_manual(values = c(VSL3="#ef857c",placebo="#28658f")) +
  facet_wrap(~compound,nrow = 3, ncol = 3,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  geom_point(aes(colour=treatment, shape=treatment), position=position_jitter(seed=1,width = .3),alpha = 1, size=2) +
  ylab("") +
  xlab("") +
  theme(aspect.ratio=1.5, strip.text=element_text(size = 15, face = 'bold'), axis.text.x = element_text(size=15, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=15),
        legend.text = element_text(size = 15),legend.title =  element_text(size = 15), strip.text.x = element_text(size=10)) +
  stat_pvalue_manual(stats_vsl_placebo_t, label='p', size=3)
dev.off()

pdf('./metabolomica/boxplot_metabolomics_t4_resp_non_t_test.pdf', width = 10, height = 10)
ggplot(df_vsl, aes(x= outcome, y=values,color=outcome)) +
  geom_boxplot(aes(colour=outcome),width=0.4, lwd=1) +
  scale_colour_manual(values = c("#ffb703","#fb8500")) +
  facet_wrap(~compound,nrow = 3, ncol = 3,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  geom_point(aes(colour=outcome, shape=outcome), position=position_jitter(seed=1,width = .3),alpha = 1, size=2) +
  ylab("") +
  xlab("") +
  theme(aspect.ratio=1.5, strip.text=element_text(size = 15, face = 'bold'), axis.text.x = element_text(size=15, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=15),
        legend.text = element_text(size = 15),legend.title =  element_text(size = 15), strip.text.x = element_text(size=10)) +
  stat_pvalue_manual(stats_resp_non_t, label='p', size=3)
dev.off()

pdf('./metabolomica/boxplot_metabolomics_t4_vsl_placebo_wilcoxon_test.pdf', width = 10, height = 10)
ggplot(df, aes(x= treatment, y=values,color=treatment)) +
  geom_boxplot(aes(colour=treatment),width=0.4, lwd=1) +
  scale_colour_manual(values = c(VSL3="#ef857c",placebo="#28658f")) +
  facet_wrap(~compound,nrow = 3, ncol = 3,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  geom_point(aes(colour=treatment, shape=treatment), position=position_jitter(seed=1,width = .3),alpha = 1, size=2) +
  ylab("") +
  xlab("") +
  theme(aspect.ratio=1.5, strip.text=element_text(size = 15, face = 'bold'), axis.text.x = element_text(size=15, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=15),
        legend.text = element_text(size = 15),legend.title =  element_text(size = 15), strip.text.x = element_text(size=10)) +
  stat_pvalue_manual(stats_vsl_placebo_w, label='p', size=3)
dev.off()

pdf('./metabolomica/boxplot_metabolomics_t4_resp_non_wilcoxon_test.pdf', width = 10, height = 10)
ggplot(df_vsl, aes(x= outcome, y=values,color=outcome)) +
  geom_boxplot(aes(colour=outcome),width=0.4, lwd=1) +
  scale_colour_manual(values = c("#ffb703","#fb8500")) +
  facet_wrap(~compound,nrow = 3, ncol = 3,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  geom_point(aes(colour=outcome, shape=outcome), position=position_jitter(seed=1,width = .3),alpha = 1, size=2) +
  ylab("") +
  xlab("") +
  theme(aspect.ratio=1.5, strip.text=element_text(size = 15, face = 'bold'), axis.text.x = element_text(size=15, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=15),
        legend.text = element_text(size = 15),legend.title =  element_text(size = 15), strip.text.x = element_text(size=10)) +
  stat_pvalue_manual(stats_resp_non_w, label='p', size=3)
dev.off()


#### microbiome analysis results
external_results <- read.csv('./metabolomica/limma_output_T4_vsl_placebo.csv', sep = ';')
rownames(external_results) <- external_results[,1]
external_results[,1] <- NULL

rownames(external_results) <- unlist(lapply(as.list(rownames(external_results)), function(x) association[association$compounds == x, 1]))

external_results2 <- as_tibble(external_results)

external_results2$compound <- rownames(external_results)
external_results2$group2 <- 'VSL3'
external_results2$group1 <- 'placebo'

external_results2$y.position <- unlist(lapply(seq_along(1:length(external_results2$compound)), function(x) {
  max(df[df$compound == external_results2$compound[x],'values'][! is.na(df[df$compound == external_results2$compound[x],'values'])]) +max(df[df$compound == external_results2$compound[x],'values'][! is.na(df[df$compound == external_results2$compound[x],'values'])])/14
}))

pdf('./metabolomica/boxplot_metabolomics_t4_vsl_placebo_external_results.pdf', width = 10, height = 10)
ggplot(df, aes(x= treatment, y=values,color=treatment)) +
  geom_boxplot(aes(colour=treatment),width=0.4, lwd=1) +
  scale_colour_manual(values = c(VSL3="#ef857c",placebo="#28658f")) +
  facet_wrap(~compound,nrow = 3, ncol = 3,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  geom_point(aes(colour=treatment, shape=treatment), position=position_jitter(seed=1,width = .3),alpha = 1, size=2) +
  ylab("") +
  xlab("") +
  theme(aspect.ratio=1.5, strip.text=element_text(size = 15, face = 'bold'), axis.text.x = element_text(size=15, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=15),
        legend.text = element_text(size = 15),legend.title =  element_text(size = 15), strip.text.x = element_text(size=10)) +
  stat_pvalue_manual(external_results2, label='P_value', size=3)
dev.off()


external_results_outcome <- read.csv('./metabolomica/limma_output_T4_resp_non.csv', sep = ';')
rownames(external_results_outcome) <- external_results_outcome[,1]
external_results_outcome[,1] <- NULL

rownames(external_results_outcome) <- unlist(lapply(as.list(rownames(external_results_outcome)), function(x) association[association$compounds == x, 1]))

external_results_outcome2 <- as_tibble(external_results_outcome)

external_results_outcome2$compound <- rownames(external_results_outcome)
external_results_outcome2$group2 <- 'Responder'
external_results_outcome2$group1 <- 'Non responder'

external_results_outcome2$y.position <- unlist(lapply(seq_along(1:length(external_results_outcome2$compound)), function(x) {
  max(df_vsl[df_vsl$compound == external_results_outcome2$compound[x],'values'][! is.na(df_vsl[df_vsl$compound == external_results_outcome2$compound[x],'values'])]) +max(df_vsl[df_vsl$compound == external_results_outcome2$compound[x],'values'][! is.na(df_vsl[df_vsl$compound == external_results_outcome2$compound[x],'values'])])/14
}))

pdf('./metabolomica/boxplot_metabolomics_t4_responder_non_external_results.pdf', width = 10, height = 10)
ggplot(df_vsl, aes(x= outcome, y=values,color=outcome)) +
  geom_boxplot(aes(colour=outcome),width=0.4, lwd=1) +
  scale_colour_manual(values = c("#ffb703","#fb8500")) +
  facet_wrap(~compound,nrow = 3, ncol = 3,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  geom_point(aes(colour=outcome, shape=outcome), position=position_jitter(seed=1,width = .3),alpha = 1, size=2) +
  ylab("") +
  xlab("") +
  theme(aspect.ratio=1.5, strip.text=element_text(size = 15, face = 'bold'), axis.text.x = element_text(size=15, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=15),
        legend.text = element_text(size = 15),legend.title =  element_text(size = 15), strip.text.x = element_text(size=10)) +
  stat_pvalue_manual(external_results_outcome2, label='P_value', size=3)
dev.off()



### rf

# vsl vs placebo t4

metabolites <-  metabolomic[c('C00163', 'C03758'),]
# 30 samples
cytochines <- sist_imm_t4[sist_imm_t4$molecule %in% c('IFNa', 'IL10', 'IL7','GCS-F'),c(1,2,4)]
cytochines$sample <- unlist(lapply(as.list(cytochines$sample), function(x) {rownames(sample_data(ps_16S)[sample_data(ps_16S)$time == 'T4' & sample_data(ps_16S)$sample_id == x,])}))
cytochines2 <- matrix(nrow=length(unique(cytochines$molecule)), ncol=length(unique(cytochines$sample)))
rownames(cytochines2) <- unique(cytochines$molecule)
colnames(cytochines2) <- unique(cytochines$sample)

for (r in rownames(cytochines2)) {
  for (c in colnames(cytochines2)) {
    cytochines2[r,c] <- cytochines[cytochines$sample == c & cytochines$molecule == r,'value']
  }
}

ot <- otu_table(ps_16S)
ot <- ot[,rownames(sample_data(ps_16S)[sample_data(ps_16S)$time =='T4' & sample_data(ps_16S)$treatment == 'VSL3',])]
rownames(ot) <- gsub('DENOVO', 'ASV_', rownames(ot))
asv <- ot[rownames(res2[res2$diff != 'NO',])]
rownames(asv) <- res2[rownames(asv),'asv_names']

s <- intersect(intersect(colnames(asv), colnames(cytochines2)),colnames(metabolites))

s <- intersect(colnames(asv),colnames(metabolites))

df_total <- rbind(metabolites[,colnames(metabolites) %in% s],asv[,colnames(asv) %in% s])

metadata <- sample_data(ps_16S_T4_vsl)[colnames(df_total), 'outcome']

df_total <- data.frame(t(df_total))
# 6 responders + 1 non responder

######
# predict responders and not at t0

df_total$outcome <- metadata$outcome
df_total$outcome <- factor(df_total$outcome)

set.seed(88)

RF_state_classify <- randomForest(x=df_total[,1:(ncol(df_total)-1)], y=df_total[ , ncol(df_total)] , ntree=501, importance=TRUE, proximities=TRUE )

# Assessing Model Fit

RF_state_classify
#ntree = 501, importance = TRUE, proximities = TRUE) 
#  Type of random forest: classification
# Number of trees: 501
# No. of variables tried at each split: 34
# 
# OOB estimate of  error rate: 17.39%
# Confusion matrix:
#   Non_responder Responder class.error
# Non_responder             0         4           1
# Responder                 0        19           0

# Permutation Test
RF_state_classify_sig <- rf.significance( x=RF_state_classify ,  xdata=df_total[,1:(ncol(df_total)-1)] , nperm=1000 , ntree=501 )  
#Number of permutations:  1000 
# p-value:  0.01 
# Model signifiant at p = 0.01 
# Model OOB error:  0.173913 
# Random OOB error:  0.173913 
# min random global error: 0.08695652 
# max random global error:  0.3478261 
# min random within class error: NA 
# max random within class error:  NA 

#Accuracy Estimated by Cross-validation
fit_control <- trainControl(method = "LOOCV")    
RF_state_classify_loocv <- train(df_total[,1:(ncol(df_total)-1)], y=df_total[, ncol(df_total)], 
                                 method="rf", ntree=501 , tuneGrid=data.frame( mtry=25), trControl=fit_control, metric="Accuracy")
RF_state_classify_loocv$results
#mtry Accuracy Kappa
# 1   25 0.826087     0

control <- trainControl(method="cv", classProbs = TRUE)
set.seed(55)
fit <- train(otu_table_scaled_T0_vsl_outcome[,1:(ncol(otu_table_scaled_T0_vsl_outcome)-1)], y=otu_table_scaled_T0_vsl_outcome[, ncol(otu_table_scaled_T0_vsl_outcome)], 
             method="rf", metric="ROC", trControl=control, ntree=501 , tuneGrid=data.frame(mtry=25))
print(fit)
#  Accuracy   Kappa    
# 0.85  0

#Identifying Important Features
RF_outcome_T0_classify_imp <- as.data.frame( RF_state_classify$importance )
# add a column "feature" reporting the genus of each row (that is each feature of the model, so each genus)
RF_outcome_T0_classify_imp$features <- rownames( RF_outcome_T0_classify_imp )
# reorder the rows of the object on the base of decreasing values of the variable "MeanDecreaseGini"
RF_outcome_T0_classify_imp_sorted <- arrange( RF_outcome_T0_classify_imp  , desc(MeanDecreaseGini)  )

otu_to_taxa <- read.table('./denovo_unoise_sv/taxa_split.txt',fill = , sep = '\t', row.names = 1)
otu_to_taxa[otu_to_taxa == ''] <- NA

v <- unlist(lapply(seq_along(1:nrow(otu_to_taxa)), function(x) {otu_to_taxa[x,6]}))
names(v) <- rownames(otu_to_taxa)    

# order the genus for increasing value of gini coefficient. convert the features variable in factor
zzz <- RF_outcome_T0_classify_imp[order(RF_outcome_T0_classify_imp$MeanDecreaseGini, decreasing = TRUE),]
zzz$features <- ifelse(is.na(v[zzz$features]), gsub('DENOVO','ASV_',zzz$features), paste0(v[zzz$features], '_',gsub('DENOVO','ASV_',zzz$features)))
zzz2 <- zzz[10:1,]

pdf('./R_analysis_all_run1/genus_gini_coef_T0_predict_outcome.pdf', height = 2, width = 3.5)
# plot the genus in order of increasing gini coefficient
zzz2$name2 <- 1:nrow(zzz2)
ggplot(zzz2, aes(x=factor(name2), y=MeanDecreaseGini)) +
  geom_segment( aes(x=factor(name2), xend=factor(name2), y=0, yend=MeanDecreaseGini), color="black") +
  geom_point( color="orange", size=3, alpha=0.7) +
  theme_bw() +
  scale_x_discrete(labels=zzz2$features) +
  coord_flip() +
  ggtitle('OOB=17.39%;p=0.01;\nAccuracy=0.85;Kappa=0')+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=7),
    axis.text.x = element_text(size=7),
    axis.title.x = element_text(size=7),
    title = element_text(size=7))
dev.off()

##
metadata_T0_VSL_responders <- rownames(metadata_T0_VSL[metadata_T0_VSL$outcome=='Responder',])
metadata_T0_VSL_non_responders <- rownames(metadata_T0_VSL[metadata_T0_VSL$outcome=='Non_responder',])

otu_to_taxa <- read.table('./denovo_unoise_sv/taxa_split.txt',fill = , sep = '\t', row.names = 1)
otu_to_taxa[otu_to_taxa == ''] <- NA

v <- unlist(lapply(seq_along(1:nrow(otu_to_taxa)), function(x) {otu_to_taxa[x,6]}))
names(v) <- rownames(otu_to_taxa)    

df <- data.frame(samples = rep(c(metadata_T0_VSL_responders,metadata_T0_VSL_non_responders), each=10),
                 genus = rep(RF_outcome_T0_classify_imp_sorted$features[1:10],length(c(metadata_T0_VSL_responders,metadata_T0_VSL_non_responders))),
                 outcome=c(rep('Responder', length(metadata_T0_VSL_responders)*10), rep('Non responder', length(metadata_T0_VSL_non_responders)*10)))

df$abundance <- NA

for (r in 1:nrow(df)) {
  df$abundance[r] <- otu_table_T0_VSL[df[r,2],df[r,1]]
}

df$outcome <- as.factor(df$outcome)
df$genus <- ifelse(is.na(v[df$genus]), gsub('DENOVO','ASV_',df$genus), paste0(v[df$genus], '_',gsub('DENOVO','ASV_',df$genus)))
df$genus <- as.factor(df$genus)

stat.test <-  df %>%
  group_by(genus) %>%
  rstatix::wilcox_test(abundance~outcome) %>% 
  adjust_pvalue(method = 'fdr') %>%
  add_significance("p.adj") 

p_top10 <- round(stat.test$p,3)
names(p_top10) <- stat.test$genus

stat.test_less <-  df %>%
  group_by(genus) %>%
  rstatix::wilcox_test(abundance~outcome, alternative = 'less') %>% 
  adjust_pvalue(method = 'fdr') %>%
  add_significance("p.adj") 

effect_size_top10 <- qnorm(stat.test_less$p, lower.tail = FALSE)/sqrt(stat.test_less$n1+stat.test_less$n2)
names(effect_size_top10) <- stat.test_less$genus

df <- data.frame(effsize = effect_size_top10, p = p_top10, asv = names(effect_size_top10))
df$asv <- factor(df$asv, levels=df[order(df$effsize, decreasing = TRUE),'asv'])

df$outcome <- "Non responder"
df$outcome[df$effsize > 0] <- "Responder"
df$outcome <- factor(df$outcome, levels=c("Responder","Non responder"))

pdf('./R_analysis_all_run1/effect_size_responders_nr_T0.pdf', width = 5, height = 3)
ggplot(df, aes(y=effsize,x=asv,fill=outcome)) +
  geom_bar(stat = "identity", width = 0.7) +  coord_flip() +
  geom_text(aes(label = p), position=position_dodge(width=0.9),size=3,
            hjust = ifelse(df$effsize >= 0, 1.2, -0.2)) +
  theme_classic() +
  scale_fill_manual(values = c("#fb8500","#ffb703")) +
  theme(axis.title.y = element_blank(), legend.title = element_blank(), legend.position = 'top') +
  ylab('Effect size') +
  scale_x_discrete(limits = rev(levels(df$asv)))
dev.off()

###AUROC
# AUROC to understand if the created model classifies the samples in a good way
# Load required packages

# Get the true class labels
true_labels <- metadata_T0_VSL$outcome
# Get the predicted class probabilities for the positive class (the second column
# corresponds to the probability that the sample is HS)
predicted_probs <- RF_state_classify$votes[, 2]
# Calculate AUROC
roc_obj <- roc(true_labels, predicted_probs)
auc_value <- auc(roc_obj)
# auc function computes the numeric value of area under the ROC curve (
# Plot ROC curve
pdf('./R_analysis_all_run1/ROC_T0_vsl_predict_outcome.pdf')
plot(roc_obj, main = "ROC Curve", print.auc = TRUE)
# Print AUROC value
cat("AUROC:", auc_value, "\n")
dev.off()

######### correlation matrix
# vsl t4 resp vs non resp
#metabolite <- t(data.frame(ko00260 = tax4fun_pathways_summary_t4[tax4fun_pathways_summary_t4$taxa == 'ko00260',c(3)]))
#colnames(metabolite) <- tax4fun_pathways_summary_t4[tax4fun_pathways_summary_t4$taxa == 'ko00260',c(2)]

metabolomic <- read.table('./metabolomica/values_t4 copia.txt', sep = '\t', header = TRUE)
rownames(metabolomic) <- metabolomic$Compound
metabolomic$Compound <- NULL
rownames(metabolomic) <- unlist(lapply(as.list(rownames(metabolomic)), function(x) association[association$compounds ==x, 1]))

metabolomics <-  metabolomic[c('isovaleric acid'),]
cytochines <- sist_imm_t4_vsl[sist_imm_t4_vsl$molecule %in% c('IL10', 'IL7','GCS-F'),c(1,2,4)]
cytochines$sample <- unlist(lapply(as.list(cytochines$sample), function(x) {rownames(sample_data(ps_16S)[sample_data(ps_16S)$time == 'T4' & sample_data(ps_16S)$sample_id == x,])}))
cytochines2 <- matrix(nrow=length(unique(cytochines$molecule)), ncol=length(unique(cytochines$sample)))
rownames(cytochines2) <- unique(cytochines$molecule)
colnames(cytochines2) <- unique(cytochines$sample)

for (r in rownames(cytochines2)) {
  for (c in colnames(cytochines2)) {
    cytochines2[r,c] <- cytochines[cytochines$sample == c & cytochines$molecule == r,'value']
  }
}

ot <- otu_table(ps_16S)
ot <- ot[,rownames(sample_data(ps_16S)[sample_data(ps_16S)$time =='T4',])]
rownames(ot) <- gsub('DENOVO', 'ASV_', rownames(ot))
asv <- ot[rownames(res2[res2$diff != 'NO',])]
rownames(asv) <- res2[rownames(asv),'asv_names']

quest <- data.frame(readxl::read_excel('./../questionari/valuesT4.xlsx'))
#quest_sign <- quest[,c('SF36.Role.limit.emotional','CFS', 'SF36.Phys.func','SF36.Health.change')]
rownames(quest) <- quest$sample
quest$sample <- NULL
quest$sample_id <- NULL
quest$...2 <- NULL
quest_t4_vsl3 <- quest[rownames(quest) %in% rownames(sample_data(ps_16S)[sample_data(ps_16S)$treatment =='VSL3' & sample_data(ps_16S)$time =='T4',]),]

s <- intersect(intersect(intersect(rownames(quest_t4_vsl3), colnames(asv)), colnames(cytochines2)), colnames(metabolomics))

cyto_metabolites_asv <- rbind(asv[,colnames(asv) %in% s], 
                  metabolomics[,colnames(metabolomics) %in% s],cytochines2[,colnames(cytochines2) %in% s])

cyto_metabolites_asv <- t(cyto_metabolites_asv)

cyto_metabolites_asv <- cyto_metabolites_asv[,-which(colSums(cyto_metabolites_asv) == 0)]

#rownames(df_total)[rownames(df_total) == '35'] <- 'ko00260'

quest_t4_vsl3 <- quest_t4_vsl3[rownames(quest_t4_vsl3) %in% s,]

ct1 = corr.test(cyto_metabolites_asv,quest_t4_vsl3,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_r <- clinic_r[,-12]
clinic_p<-ct1$p
clinic_p <- ifelse(clinic_p <= 0.05, '*','')
clinic_p <- clinic_p[,-12]

distance = dist(clinic_r, method = "euclidean")
rowclus <- hclust(dist(clinic_r, method = "euclidean"))
colclus <- hclust(dist(t(clinic_r), method = "euclidean"))

my_palette <- colorRampPalette(brewer.pal(10, "RdYlBu"))(n = 256)

pdf('../immagini_articolo/heatmap_correlation_T4_resp_non_quest_vs_all.pdf', height = 10, width = 10)
pheatmap(clinic_r, cluster_rows = TRUE, cluster_cols = TRUE, color = rev(my_palette), cellwidth = 15, cellheight = 15, angle_col = 315, 
        border_color = 'white', display_numbers = clinic_p, fontsize_number = 17, number_color = 'white', fontsize_row = 15, fontsize_col = 15, )
# annotation_col  = ann, annotation_colors = ann_colors, 
dev.off()

######### correlation matrix
# vsl vs placebo t4 
#metabolite <- t(data.frame(ko00260 = tax4fun_pathways_summary_t4[tax4fun_pathways_summary_t4$taxa == 'ko00260',c(3)]))
#colnames(metabolite) <- tax4fun_pathways_summary_t4[tax4fun_pathways_summary_t4$taxa == 'ko00260',c(2)]

metabolomic <- read.table('./metabolomica/values_t4 copia.txt', sep = '\t', header = TRUE)
rownames(metabolomic) <- metabolomic$Compound
metabolomic$Compound <- NULL
rownames(metabolomic) <- unlist(lapply(as.list(rownames(metabolomic)), function(x) association[association$compounds ==x, 1]))

metabolomics <-  metabolomic[c('serotonine', 'dopamine','acetic acid', 'butyric acid'),]
cytochines <- sist_imm_t4[sist_imm_t4$molecule %in% c('IFNa', 'GCS-F', 'IL-1a','IL10','IL7','IL6'),c(1,2,4)]
cytochines$sample <- unlist(lapply(as.list(cytochines$sample), function(x) {rownames(sample_data(ps_16S)[sample_data(ps_16S)$time == 'T4' & sample_data(ps_16S)$sample_id == x,])}))
cytochines2 <- matrix(nrow=length(unique(cytochines$molecule)), ncol=length(unique(cytochines$sample)))
rownames(cytochines2) <- unique(cytochines$molecule)
colnames(cytochines2) <- unique(cytochines$sample)

for (r in rownames(cytochines2)) {
  for (c in colnames(cytochines2)) {
    cytochines2[r,c] <- cytochines[cytochines$sample == c & cytochines$molecule == r,'value']
  }
}

ot <- otu_table(ps_16S)
ot <- ot[,rownames(sample_data(ps_16S)[sample_data(ps_16S)$time =='T4',])]
rownames(ot) <- gsub('DENOVO', 'ASV_', rownames(ot))
asv <- ot[rownames(res2[res2$diff != 'NO',])]
rownames(asv) <- res2[rownames(asv),'asv_names']

quest <- data.frame(readxl::read_excel('./../questionari/valuesT4.xlsx'))
#quest_sign <- quest[,c('SF36.Role.limit.emotional','CFS', 'SF36.Phys.func','SF36.Health.change')]
rownames(quest) <- quest$sample
quest$sample <- NULL
quest$sample_id <- NULL
quest_t4 <- quest[rownames(quest) %in% rownames(sample_data(ps_16S)[sample_data(ps_16S)$time =='T4',]),]

s <- intersect(intersect(intersect(rownames(quest_t4), colnames(asv)), colnames(cytochines2)), colnames(metabolomics))

cyto_metabolites_asv <- rbind(asv[,colnames(asv) %in% s], 
                              metabolomics[,colnames(metabolomics) %in% s],cytochines2[,colnames(cytochines2) %in% s])

cyto_metabolites_asv <- t(cyto_metabolites_asv)

#rownames(df_total)[rownames(df_total) == '35'] <- 'ko00260'

quest_t4 <- quest_t4[rownames(quest_t4) %in% s,]

cyto_metabolites_asv2 <- cyto_metabolites_asv[intersect(rownames(cyto_metabolites_asv), rownames(quest_t4)),]
quest_t4 <- quest_t4[intersect(rownames(cyto_metabolites_asv), rownames(quest_t4)),]

ct1 = corr.test(cyto_metabolites_asv2,quest_t4,method="spearman", adjust="fdr")
clinic_r<-ct1$r
clinic_p<-ct1$p
clinic_p <- ifelse(clinic_p <= 0.05, '*','')

distance = dist(clinic_r, method = "euclidean")
rowclus <- hclust(dist(clinic_r, method = "euclidean"))
colclus <- hclust(dist(t(clinic_r), method = "euclidean"))

my_palette <- colorRampPalette(brewer.pal(10, "RdYlBu"))(n = 256)

pdf('../immagini_articolo/heatmap_correlation_T4_treatment_quest_vs_all.pdf', height = 12, width = 10)
pheatmap(clinic_r, cluster_rows = TRUE, cluster_cols = TRUE, color = rev(my_palette), cellwidth = 15, cellheight = 15, angle_col = 315, 
         border_color = 'white', display_numbers = clinic_p, fontsize_number = 17, number_color = 'white', fontsize_row = 15, fontsize_col = 15, )
# annotation_col  = ann, annotation_colors = ann_colors, 
dev.off()
