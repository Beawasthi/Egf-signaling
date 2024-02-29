path="/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/"
figure_path="/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/R_figures_for_paper/"

library(tidyverse)
library(VennDiagram)
library(ggfortify)
library(ggplot2)
library("RColorBrewer")
library("pheatmap")
library("grid")
library(lme4)
library(lmerTest)
library(ggsignif)
library(dendextend)
library(UpSetR)

# Functions
get_tissue_all <- function(tissue){
  phospho <- read_csv(paste0("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/final_MSstats_files/ Phospho_ProteinAdjusted_StatResults_ ", tissue, " .csv"))
  phospho <- phospho %>%
    select(Reference, PhosphoSite, paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl")) %>%
    rename(log2FC = paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           qValue = paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           significance = paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"))
  all_phospho <- phospho %>%
    mutate(protein = Reference, site = str_extract(PhosphoSite, "[[:digit:]]+")) %>%
    mutate(ProteinID = protein) %>%
    unite(protein_site, protein, site) 
  return(list(all_phospho$protein_site, all_phospho$ProteinID)) 
}

get_tissue_sig <- function(tissue){
  phospho <- read_csv(paste0("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/final_MSstats_files/ Phospho_ProteinAdjusted_StatResults_ ", tissue, " .csv"))
  phospho <- phospho %>%
    select(Reference, PhosphoSite, paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl")) %>%
    rename(log2FC = paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           qValue = paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           significance = paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"))
  sig_phospho <- phospho %>%
    mutate(protein = Reference, site = str_extract(PhosphoSite, "[[:digit:]]+")) %>%
    mutate(ProteinID = protein) %>%
    unite(protein_site, protein, site) %>%
    filter(!(significance == "n.s."))
  View(sig_phospho)
  return(list(sig_phospho$protein_site, sig_phospho$ProteinID))
}

scale_data <- function(df){
  scaled <- df %>%
    mutate_at(vars(heart_EGF_rep2_b:liver_vehctrl_rep3_b), ~.x/rowSums(select(df, heart_EGF_rep2_b:liver_vehctrl_rep3_b))*100) %>%
    mutate_at(vars(heart_EGF_rep1_a:liver_vehctrl_rep2_a), ~.x/rowSums(select(df, heart_EGF_rep1_a:liver_vehctrl_rep2_a))*100) 
  return(scaled)
}

zscorer <- function(x) {
  (x - mean(x))/(sd(x))
}

make_PCA <- function(the_data, type){
  df <- df2 <- the_data %>%
    select(heart_EGF_rep1:liver_vehctrl_rep3) %>% 
    #select(kidney_EGF_rep1:kidney_EGF_rep3, kidney_veh_ctrl_rep1:kidney_veh_ctrl_rep3) %>% 
    na.omit() %>% t() %>% data.frame()
  df2$tissue <- c(rep("heart_EGF", 3), rep("lung_EGF", 3), rep("kidney_EGF", 3), rep("liver_EGF", 3),
                  rep("heart_veh", 3), rep("lung_veh", 3), rep("kidney_veh", 3), rep("liver_veh", 3))
  #df2$tissue <- c(rep("lung_EGF", 3), rep("lung_veh_ctrl", 3))
  pca_res <- prcomp(df, scale = TRUE, center = TRUE)
  g <-  autoplot(pca_res, data = df2, colour = "tissue", frame = TRUE) + theme_classic()
  pdf(paste0(figure_path ,type, "_full_PCA.pdf"))
  plot(g)
  dev.off()
}

make_dendrogram <- function(df, type){
  df_for_dend <- df %>% 
    select(heart_EGF_rep1:liver_vehctrl_rep3) %>%
    t() %>% data.frame()
  rownames(df_for_dend) <- c("Heart EGF r1", "Heart EGF r2", "Heart EGF r3",
                             "Lung EGF r1", "Lung EGF r2", "Lung EGF r3",
                             "Kidney EGF r1", "Kidney EGF r2", "Kidney EGF r3",
                             "Liver EGF r1", "Liver EGF r2", "Liver EGF r3",
                             "Heart veh r1", "Heart veh r2", "Heart veh r3",
                             "Lung veh r1", "Lung veh r2", "Lung veh r3",
                             "Kidney veh r1", "Kidney veh r2", "Kidney veh r3",
                             "Liver veh r1", "Liver veh r2", "Liver veh r3")
  
  dend <- df_for_dend %>% 
    dist() %>% 
    hclust(method = "ward.D2") %>% 
    as.dendrogram() 
  
  EGF <- grepl("EGF", labels(dend))
  veh <- grepl("veh", labels(dend))
  
  
  pdf(paste0(figure_path, type, "_full_dendrogram.pdf"), height = 25, width = 35)
  dend %>%
    set("leaves_cex", 1) %>%
    set("leaves_pch", 15) %>%
    set("leaves_col", case_when(EGF~"red",
                                veh~"black")) %>%
    set("labels_colors", case_when(EGF~"red",
                                   veh~"black")) %>% 
    plot(horiz = TRUE)
  dev.off()
}


make_boxplot_baseline_phospho <- function(proteinsite, upper){
  filtered_sites <- mean_normalized_scaled_phospho %>%
    #mutate(sums = rowSums(select(total_for_sig_heat_map, heart_EGF_rep1:liver_vehctrl_rep3))) %>%
    select(protein_site, heart_vehctrl_rep1:heart_vehctrl_rep3,
           lung_vehctrl_rep1:lung_vehctrl_rep3,
           kidney_vehctrl_rep1:kidney_vehctrl_rep3,
           liver_vehctrl_rep1:liver_vehctrl_rep3) %>%
    filter(protein_site == proteinsite)
  filtered_sites_for_barplot <- filtered_sites %>%
    gather(key = "sample", value = "intensity", -protein_site) %>%
    mutate(tissue = c(rep("Heart", 3),
                      rep("Lung", 3),
                      rep("Kidney", 3),
                      rep("Liver", 3))) %>%
    mutate(intensity = log(intensity, 2))
  
  filtered_sites_for_barplot$tissue <- factor(filtered_sites_for_barplot$tissue,
                                              levels = unique(filtered_sites_for_barplot$tissue), ordered = TRUE)
  filtered_sites_for_barplot$sample <- factor(filtered_sites_for_barplot$sample,
                                              levels = filtered_sites_for_barplot$sample)
  
  
  g <- ggplot(filtered_sites_for_barplot,
              aes(x = factor(tissue), y = intensity,  fill = factor(tissue))) +
    geom_boxplot() +
    #geom_signif(comparisons = list(c("Heart", "Lung")), y_position = 3.0, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Heart", "Kidney")),  y_position = 3.8, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Heart", "Liver")),  y_position = 4.2, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Lung", "Kidney")),  y_position = 3.2, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Lung", "Liver")),  y_position = 4.0, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Kidney", "Liver")),  y_position = 3.6, map_signif_level = TRUE) +
    scale_fill_manual(values=c("darkred","blue", "purple", "orange")) +
    scale_y_continuous(breaks = seq(0, upper, by = 1), limits = c(0, upper)) +
    ylab("log2(intensity)") +
    labs(title=toupper(proteinsite),
         x ="Tissue", y = "log2(intensity)") +
    labs(fill = "Tissue") +
    theme_classic() + theme(plot.title=element_text(hjust=0.5, vjust = 0)) +  
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
  plot(g)
  pdf(paste0(figure_path, "full_phospho_baseline_tissue_plex_boxplot_", proteinsite, ".pdf")) 
  plot(g)
  dev.off()
}


make_boxplot_stimulated_phospho <- function(proteinsite, upper){
  filtered_sites <- mean_normalized_scaled_phospho %>%
    #mutate(sums = rowSums(select(total_for_sig_heat_map, heart_EGF_rep1:liver_vehctrl_rep3))) %>%
    select(protein_site, heart_EGF_rep1:heart_EGF_rep3, 
           lung_EGF_rep1:lung_EGF_rep3,
           kidney_EGF_rep1:kidney_EGF_rep3,
           liver_EGF_rep1:liver_EGF_rep3, ) %>%
    filter(protein_site == proteinsite)
  
  filtered_sites_for_barplot <- filtered_sites %>%
    gather(key = "sample", value = "intensity", -protein_site) %>%
    mutate(tissue = c(rep("Heart", 3),
                      rep("Lung", 3),
                      rep("Kidney", 3),
                      rep("Liver", 3))) %>%
    mutate(intensity = log(intensity, 2))
  
  filtered_sites_for_barplot$tissue <- factor(filtered_sites_for_barplot$tissue,
                                              levels = unique(filtered_sites_for_barplot$tissue), ordered = TRUE)
  filtered_sites_for_barplot$sample <- factor(filtered_sites_for_barplot$sample,
                                              levels = filtered_sites_for_barplot$sample)
  
  
  g <- ggplot(filtered_sites_for_barplot,
              aes(x = factor(tissue), y = intensity,  fill = factor(tissue))) +
    geom_boxplot() +
    #geom_signif(comparisons = list(c("Heart", "Lung")), y_position = 3.0, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Heart", "Kidney")),  y_position = 3.8, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Heart", "Liver")),  y_position = 4.2, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Lung", "Kidney")),  y_position = 3.2, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Lung", "Liver")),  y_position = 4.0, map_signif_level = TRUE) +
    #geom_signif(comparisons = list(c("Kidney", "Liver")),  y_position = 3.6, map_signif_level = TRUE) +
    scale_fill_manual(values=c("darkred","blue", "purple", "orange")) +
    scale_y_continuous(breaks = seq(0, upper, by = 1), limits = c(0, upper)) +
    ylab("log2(intensity)") +
    labs(title=toupper(proteinsite),
         x ="Tissue", y = "log2(intensity)") +
    labs(fill = "Tissue") +
    theme_classic() + theme(plot.title=element_text(hjust=0.5, vjust = 0)) +  
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
  plot(g)
  pdf(paste0(figure_path, "full_phospho_stimulated_tissue_plex_boxplot_", proteinsite, ".pdf")) 
  plot(g)
  dev.off()
}


make_total_boxplot <- function(proteinID, upper){
  filtered_total_proteins <- mean_normalized_scaled_total %>%
    #mutate(sums = rowSums(select(total_for_sig_heat_map, heart_EGF_rep1:liver_vehctrl_rep3))) %>%
    select(ProteinID, heart_EGF_rep1:heart_EGF_rep3, heart_vehctrl_rep1:heart_vehctrl_rep3,
           lung_EGF_rep1:lung_EGF_rep3, lung_vehctrl_rep1:lung_vehctrl_rep3,
           kidney_EGF_rep1:kidney_EGF_rep3, kidney_vehctrl_rep1:kidney_vehctrl_rep3,
           liver_EGF_rep1:liver_EGF_rep3, liver_vehctrl_rep1:liver_vehctrl_rep3) %>%
    filter(ProteinID == proteinID)
  print(head(filtered_total_proteins))
  
  filtered_total_proteins_for_barplot <- filtered_total_proteins %>%
    gather(key = "sample", value = "intensity", -ProteinID) %>%
    mutate(tissue = c(rep("Heart", 6),
                      rep("Lung", 6),
                      rep("Kidney", 6),
                      rep("Liver", 6))) %>%
    mutate(intensity = log(intensity, 2))
  write_csv(filtered_total_proteins_for_barplot, paste0(proteinID, "_filtered_total_proteins_for_barplot.csv"))
  filtered_total_proteins_for_barplot$tissue <- factor(filtered_total_proteins_for_barplot$tissue,
                                                       levels = unique(filtered_total_proteins_for_barplot$tissue), ordered = TRUE)
  filtered_total_proteins_for_barplot$sample <- factor(filtered_total_proteins_for_barplot$sample,
                                                       levels = filtered_total_proteins_for_barplot$sample)
  
  
  g <- ggplot(filtered_total_proteins_for_barplot,
              aes(x = factor(tissue), y = intensity,  fill = factor(tissue))) +
    geom_boxplot() +
    geom_signif(comparisons = list(c("Heart", "Lung")), y_position = 3.0, map_signif_level = TRUE) +
    geom_signif(comparisons = list(c("Heart", "Kidney")),  y_position = 3.8, map_signif_level = TRUE) +
    geom_signif(comparisons = list(c("Heart", "Liver")),  y_position = 4.2, map_signif_level = TRUE) +
    geom_signif(comparisons = list(c("Lung", "Kidney")),  y_position = 3.2, map_signif_level = TRUE) +
    geom_signif(comparisons = list(c("Lung", "Liver")),  y_position = 4.0, map_signif_level = TRUE) +
    geom_signif(comparisons = list(c("Kidney", "Liver")),  y_position = 3.6, map_signif_level = TRUE) +
    scale_fill_manual(values=c("darkred","blue", "purple", "orange")) +
    scale_y_continuous(breaks = seq(0, upper, by = 1), limits = c(0, upper)) +
    ylab("log2(intensity)") +
    labs(title=toupper(proteinID),
         x ="Tissue", y = "log2(intensity)") +
    labs(fill = "Tissue") +
    theme_classic() + theme(plot.title=element_text(hjust=0.5, vjust = 0)) +  
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
  plot(g)
  pdf(paste0(figure_path, "full_tissue_plex_boxplot_", proteinID, ".pdf")) 
  plot(g)
  dev.off()
}


calculate_total_entropy_piece <- function(x){
  (x + 1/16) /(1 + sum(x)) *log((1 + sum(x))/(x + 1/16))
} 

make_total_tissue_entropy_sig_barplots <- function(tissue){
  selected_sig_for_dist <- entropied_scale_total %>% # from entropied_scale_total, select tissue columns of interest
    select(ProteinID, paste0(tissue, "_sig"), entropy) %>%
    rename(tissue_sig = paste0(tissue, "_sig")) %>%
    filter(!tissue_sig == "n.d.") %>% # filter out all undetected
    filter(!(ProteinID %in% get_total_sig_entropy_inconsistencies$ProteinID)) # filter out all entropy inconsistencies

  selected_spec <- sig_entropy_scale_total %>% # get only tissue-specific proteins
    rename(specificity = paste0(tissue, "_spec")) %>%
    filter(specificity == "spec") %>%
    pull(., ProteinID)
  print(head(selected_spec))
  #selected_sig_for_dist <- selected_sig_for_dist %>%
  #  mutate(specificity = ifelse(ProteinID %in% selected_spec, "spec", ifelse(entropy < 2.58, "other", "non-spec"))) 
  selected_sig_for_dist <- selected_sig_for_dist %>%
    mutate(specificity = ifelse(ProteinID %in% selected_spec, "spec", "non-spec"), number = 1) 
  print(head(selected_sig_for_dist))
  g <- ggplot(data = selected_sig_for_dist, aes(fill = specificity, y = number, x = tissue_sig)) + 
    geom_bar(position = "fill", stat = "identity") + theme_bw() + 
    annotate("text", x = 2, y = 1, 
             label = nrow(selected_sig_for_dist %>% filter(tissue_sig == "sig"))) + 
    annotate("text", x = 1, y = 1, 
             label = nrow(selected_sig_for_dist %>% filter(tissue_sig == "n.s."))) + theme_bw()
  pdf(paste0(figure_path, tissue, "_total_specificity_barplots.pdf"))
  plot(g)
  dev.off()
  print(dim(selected_sig_for_dist))
  print(dim(selected_sig_for_dist %>% filter(tissue_sig == "sig")))
  print(dim(selected_sig_for_dist %>% filter(tissue_sig == "sig", specificity == "spec")))
  print(dim(selected_sig_for_dist %>% filter(tissue_sig == "n.s.")))
  print(dim(selected_sig_for_dist %>% filter(tissue_sig == "n.s.", specificity == "spec")))
  write_csv(selected_sig_for_dist %>% filter(tissue_sig == "sig", specificity == "spec"), paste0(figure_path, tissue, "_sig_Egf_tissue_specific_proteins.csv"))
}

make_linear_model_plots <- function(name){
  test_for_model <-  pre_phospho_for_model %>%
    filter(protein_site ==   name)
  print(test_for_model)
  test_for_model$tissue <- as.factor(test_for_model$tissue)
  test_for_model$treatment <- as.factor(test_for_model$treatment)
  test_for_model$treatment <- relevel(test_for_model$treatment, ref =  "vehctrl")
  test_for_model$rep <- as.factor(test_for_model$rep)
  
  test_fit <- lmerTest::lmer(phospho_abundance ~ 1  + total_abundance + treatment + tissue + tissue*treatment + (1|rep), data = test_for_model)
  print(summary(test_fit))
  test_for_model <- test_for_model %>%
    unite(tissue_treatment, tissue, treatment)
  g = plot(test_fit)
  plot(g)
  g = sjPlot::plot_model(test_fit, 
                         show.values=TRUE, show.p=TRUE) + theme_classic()
  pdf(paste(figure_path, name, "_sjplot.pdf"))
  plot(g)
  dev.off()
  g = ggplot(test_for_model) + 
    geom_boxplot( aes(x = tissue_treatment, y = phospho_abundance, fill = "green")) +
    geom_boxplot( aes(x = tissue_treatment, y = total_abundance, fill = "red")) + theme_bw()
  #plot(g)
}


# # get protein ID for significant EGF-stimulated phosphosites from MSstats
heart_all <- get_tissue_all("heart")[[1]]
lung_all <- get_tissue_all("lung")[[1]]
kidney_all <- get_tissue_all("kidney")[[1]]
liver_all <- get_tissue_all("liver")[[1]]

heart_all_proteins <- get_tissue_all("heart")[[2]]
lung_all_proteins <- get_tissue_all("lung")[[2]]
kidney_all_proteins <- get_tissue_all("kidney")[[2]]
liver_all_proteins <- get_tissue_all("liver")[[2]]

heart_sig_phospho_fold_changes <- get_tissue_sig("heart")[[1]]
lung_sig_phospho_fold_changes <- get_tissue_sig("lung")[[1]]
kidney_sig_phospho_fold_changes <- get_tissue_sig("kidney")[[1]]
liver_sig_phospho_fold_changes <- get_tissue_sig("liver")[[1]]

heart_sig_phospho_fold_changes_proteins <- get_tissue_sig("heart")[[2]]
lung_sig_phospho_fold_changes_proteins <- get_tissue_sig("lung")[[2]]
kidney_sig_phospho_fold_changes_proteins <- get_tissue_sig("kidney")[[2]]
liver_sig_phospho_fold_changes_proteins <- get_tissue_sig("liver")[[2]]



# Total part

total <- read_csv("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/2022-05-16_Bea_FinalMousePro.csv") %>%
  rename(ProteinID = `Protein Id`, 
         heart_EGF_rep2_b = "Set2~rq_126_sn sum",
         heart_EGF_rep3_b = "Set2~rq_127n_sn sum",
         lung_EGF_rep2_b = "Set2~rq_127c_sn sum",	
         lung_EGF_rep3_b = "Set2~rq_128n_sn sum",
         kidney_EGF_rep2_b = "Set2~rq_128c_sn sum",
         kidney_EGF_rep3_b = "Set2~rq_129n_sn sum",	
         liver_EGF_rep2_b = "Set2~rq_129c_sn sum",	
         liver_EGF_rep3_b = "Set2~rq_130n_sn sum",	
         heart_vehctrl_rep1_b = "Set2~rq_130c_sn sum",
         heart_vehctrl_rep3_b = "Set2~rq_131_sn sum",	
         lung_vehctrl_rep1_b = "Set2~rq_131c_sn sum",
         lung_vehctrl_rep3_b = "Set2~rq_132n_sn sum",	
         kidney_vehctrl_rep1_b = "Set2~rq_132c_sn sum",
         kidney_vehctrl_rep3_b = "Set2~rq_133n_sn sum",
         liver_vehctrl_rep1_b = "Set2~rq_133c_sn sum",
         liver_vehctrl_rep3_b = "Set2~rq_134n_sn sum",
         
         heart_EGF_rep1_a = "Set1~rq_126_sn sum",
         heart_EGF_rep2_a = "Set1~rq_127n_sn sum",
         lung_EGF_rep1_a = "Set1~rq_127c_sn sum",	
         lung_EGF_rep2_a = "Set1~rq_128n_sn sum",
         kidney_EGF_rep1_a = "Set1~rq_128c_sn sum",
         kidney_EGF_rep2_a = "Set1~rq_129n_sn sum",	
         liver_EGF_rep1_a = "Set1~rq_129c_sn sum",	
         liver_EGF_rep2_a = "Set1~rq_130n_sn sum",	
         heart_vehctrl_rep1_a = "Set1~rq_130c_sn sum",
         heart_vehctrl_rep2_a = "Set1~rq_131_sn sum",	
         lung_vehctrl_rep1_a = "Set1~rq_131c_sn sum",
         lung_vehctrl_rep2_a = "Set1~rq_132n_sn sum",	
         kidney_vehctrl_rep1_a = "Set1~rq_132c_sn sum",
         kidney_vehctrl_rep2_a = "Set1~rq_133n_sn sum",
         liver_vehctrl_rep1_a = "Set1~rq_133c_sn sum",
         liver_vehctrl_rep2_a = "Set1~rq_134n_sn sum")

scaled_total <- scale_data(total) 


# Figure 5
# Venn diagram of protein numbers in each plex
scaled_total_for_venn <- scaled_total %>%
  mutate(plex1_sums = rowSums(select(scaled_total, heart_EGF_rep1_a:liver_vehctrl_rep2_a), na.rm = TRUE),
         plex2_sums = rowSums(select(scaled_total, heart_EGF_rep2_b:liver_vehctrl_rep3_b), na.rm = TRUE))


VennDiagram::venn.diagram(x = list((scaled_total_for_venn %>% filter(plex1_sums > 0) %>% pull(ProteinID)), (scaled_total_for_venn %>% filter(plex2_sums > 0)) %>% pull(ProteinID)), category.names = c("Plex 1", "Plex 2"),
                          filename = paste0(figure_path, 'shared_total_proteins_full_tissue.png'), output = TRUE, fill = c("lightgreen", "darkgoldenrod1"))

# Calculate z-score of bridge channels for each plex and select only proteins with correlation of z(bridge) > 0.9
zscores_total <- scaled_total %>%
  select(ProteinID, "Set1 Peptides", "Set2 Peptides", heart_EGF_rep2_b, lung_EGF_rep2_b, kidney_EGF_rep2_b, liver_EGF_rep2_b,
         heart_vehctrl_rep1_b, lung_vehctrl_rep1_b, kidney_vehctrl_rep1_b, liver_vehctrl_rep1_b,
         heart_EGF_rep2_a, lung_EGF_rep2_a, kidney_EGF_rep2_a, liver_EGF_rep2_a,
         heart_vehctrl_rep1_a, lung_vehctrl_rep1_a, kidney_vehctrl_rep1_a, liver_vehctrl_rep1_a) %>%
  filter(!(rowSums(select(., heart_EGF_rep2_b:liver_vehctrl_rep1_b))==0)) %>% # filter out things that were not detected in both plexes
  filter(!(rowSums(select(., heart_EGF_rep2_a:liver_vehctrl_rep1_a))==0)) %>%
  gather(-c(ProteinID, "Set1 Peptides", "Set2 Peptides"), key = "sample", value = "scaled_value") %>% 
  mutate(plex = c(rep("a", nrow(.)/2), rep("b", nrow(.)/2))) %>%
  group_by(ProteinID, plex) %>%
  mutate(zscore = zscorer(scaled_value)) %>%
  ungroup() %>%
  na.omit() %>%
  group_by(ProteinID) %>%
  summarise(cor = tryCatch(
    {cor.test(zscore[plex == "a"], zscore[plex == "b"], 
              na.action = "na.exclude")$estimate},
    error = function(e){
      0
    }
    
  )) %>%
  filter(cor > 0.9)

# normalize "total" by bridge channel - using UNSCALED TMT values - and filter for meeting either z-score of number of peptide criteria
mean_normalized_total <- total %>%
  rowwise() %>% # calculate mean of bridges
  mutate(set2_mean = mean(c(heart_EGF_rep2_b, lung_EGF_rep2_b, kidney_EGF_rep2_b, liver_EGF_rep2_b,
                            heart_vehctrl_rep1_b, lung_vehctrl_rep1_b, kidney_vehctrl_rep1_b, liver_vehctrl_rep1_b)),
         set1_mean = mean(c(heart_EGF_rep2_a, lung_EGF_rep2_a, kidney_EGF_rep2_a, liver_EGF_rep2_a,
                            heart_vehctrl_rep1_a, lung_vehctrl_rep1_a, kidney_vehctrl_rep1_a, liver_vehctrl_rep1_a))) %>%
  ungroup() %>%
  mutate_at(vars(heart_EGF_rep2_b:liver_vehctrl_rep3_b), ~.x/set2_mean) %>% # normalize to bridges
  mutate_at(vars(heart_EGF_rep1_a:liver_vehctrl_rep2_a), ~.x/set1_mean) %>%
  rowwise() %>%
  mutate(heart_EGF_rep2 = mean(c(heart_EGF_rep2_a, heart_EGF_rep2_b)), # mean of replicate #2's
         lung_EGF_rep2 = mean(c(lung_EGF_rep2_a, lung_EGF_rep2_b)),
         kidney_EGF_rep2 = mean(c(kidney_EGF_rep2_a, kidney_EGF_rep2_b)),
         liver_EGF_rep2 = mean(c(liver_EGF_rep2_a, liver_EGF_rep2_b)),
         heart_vehctrl_rep1 = mean(c(heart_vehctrl_rep1_a, heart_vehctrl_rep1_b)),
         lung_vehctrl_rep1 = mean(c(lung_vehctrl_rep1_a, lung_vehctrl_rep1_b)),
         kidney_vehctrl_rep1 = mean(c(kidney_vehctrl_rep1_b, kidney_vehctrl_rep1_a)),
         liver_vehctrl_rep1 = mean(c(liver_vehctrl_rep1_a, liver_vehctrl_rep1_b))) %>%
  ungroup() %>%
  filter(ProteinID %in% zscores_total$ProteinID | `Set1 Peptides` > 3 | `Set2 Peptides` > 3) %>% # filter
  select(-contains(c("EGF_rep2_a", "EGF_rep2_b", "vehctrl_rep1_a", "vehctrl_rep1_b"))) %>%
  rename_with(~str_replace(., '_a', '')) %>%
  rename_with(~str_replace(., '_b', '')) %>%
  select(ProteinID, "Gene Symbol", "Set1 Peptides", "Set2 Peptides", heart_EGF_rep1, heart_EGF_rep2,
         heart_EGF_rep3, lung_EGF_rep1, lung_EGF_rep2, lung_EGF_rep3, kidney_EGF_rep1, kidney_EGF_rep2, kidney_EGF_rep3,
         liver_EGF_rep1, liver_EGF_rep2, liver_EGF_rep3, heart_vehctrl_rep1, heart_vehctrl_rep2, heart_vehctrl_rep3,
         lung_vehctrl_rep1, lung_vehctrl_rep2, lung_vehctrl_rep3, kidney_vehctrl_rep1, kidney_vehctrl_rep2, kidney_vehctrl_rep3,
         liver_vehctrl_rep1, liver_vehctrl_rep2, liver_vehctrl_rep3) %>%
  na.omit()

mean_normalized_scaled_total <- mean_normalized_total %>% # scale bridge-normalized dataset
  mutate_at(vars(heart_EGF_rep1:liver_vehctrl_rep3), ~.x/rowSums(select(mean_normalized_total, heart_EGF_rep1:liver_vehctrl_rep3))*100) 

# PCA and clustering
make_PCA(mean_normalized_scaled_total, "total")
make_dendrogram(mean_normalized_scaled_total, "total")


# compare total levels of proteins between tissues
make_total_boxplot("sp|Q01279|EGFR_MOUSE", 4.5)
make_total_boxplot("sp|P63085|MK01_MOUSE", 4.5)
make_total_boxplot("sp|Q63844|MK03_MOUSE", 4.5)
#make_total_boxplot("sp|P28028|BRAF_MOUSE", 4.5)
#make_total_boxplot("sp|Q99N57|RAF1_MOUSE", 4.5)
#make_total_boxplot("sp|P31938|MP2K1_MOUSE", 4.5)
#make_total-boxplot("sp|Q63932|MP2K2_MOUSE", 4.5)
#make_total_boxplot("sp|P31938|MP2K1_MOUSE", 4.5)
#make_total_boxplot("sp|P32883|RASK_MOUSE", 4.5)

# make heat map comparing total levels of proteins significantly affected in phosphorylation in response to Egf - this is using z-scores
all_tissue_sig_proteins <- c(heart_sig_phospho_fold_changes_proteins,
                             lung_sig_phospho_fold_changes_proteins,
                             kidney_sig_phospho_fold_changes_proteins,
                             liver_sig_phospho_fold_changes_proteins)

total_for_sig_heat_map <- mean_normalized_scaled_total %>%
  filter(ProteinID %in% all_tissue_sig_proteins) %>%
  #filter(ProteinID %in% phospho_EGF_for_sig_heat_map$ProteinID) %>%
  #mutate(sums = rowSums(select(total_for_sig_heat_map, heart_EGF_rep1:liver_vehctrl_rep3))) %>%
  select(ProteinID, heart_EGF_rep1:heart_EGF_rep3, heart_vehctrl_rep1:heart_vehctrl_rep3,
         lung_EGF_rep1:lung_EGF_rep3, lung_vehctrl_rep1:lung_vehctrl_rep3,
         kidney_EGF_rep1:kidney_EGF_rep3, kidney_vehctrl_rep1:kidney_vehctrl_rep3,
         liver_EGF_rep1:liver_EGF_rep3, liver_vehctrl_rep1:liver_vehctrl_rep3) %>%
  arrange(ProteinID)
for_total_heat_map <- data.frame(lapply(select(total_for_sig_heat_map, heart_EGF_rep1:liver_vehctrl_rep3), as.numeric))
scaled_total_for_heat_map <- data.frame(t(apply(for_total_heat_map, 1, zscorer)))
rownames(scaled_total_for_heat_map) <- make.names(total_for_sig_heat_map$ProteinID, unique = TRUE)
scaled_total_for_heat_map <- scaled_total_for_heat_map %>%
  filter(!rowSums(.) == 0) %>%
  distinct() %>%
  filter_all(all_vars(is.finite(.)))

draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 0.5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
paletteLength <- 45
myColor <- colorRampPalette(c("darkblue","white", "darkred"))(paletteLength)
#myColor <- colorRampPalette(c("darkblue","cornflowerblue", "lightblue3", "lightcyan1", "lightyellow", "lightgoldenrod1","darkgoldenrod3","darkorange", "red3"))(paletteLength)
#paletteLength <- 7
#myColor <- colorRampPalette(c("darkblue","lightblue3", "lightcyan1", "lightyellow", "lightgoldenrod1","darkorange", "red3"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
bk1 <- c(seq(-4.5, -0.2, by =0.2), -0.10001)
bk2 <- c(0.10001, seq(.2,4.5, by = 0.2))
bk <- c(bk1, bk2)

pdf(paste0(figure_path, "full_tissue_sig_phosphorylated_proteins.pdf"))
pheatmap(scaled_total_for_heat_map, cluster_cols = F, cluster_rows=T, 
         color=myColor, breaks=bk, border_color = "black")
dev.off()

# entropy


# 80% of signal in two tissues
8*((10 + (1/16))/101) *log(1/((10 + (1/16))/101)) + 
  8 * ((2.5 + (1/16))/101) * log(101/(2.5 + (1/16)))

# calculate entropy for each protein
total_entropy_pieces <- total %>% # scale "total" set again- in principle, the same as "scaled_total"
  mutate_at(vars(heart_EGF_rep2_b:liver_vehctrl_rep3_b), ~.x/rowSums(select(total, heart_EGF_rep2_b:liver_vehctrl_rep3_b))*100) %>%
  mutate_at(vars(heart_EGF_rep1_a:liver_vehctrl_rep2_a), ~.x/rowSums(select(total, heart_EGF_rep1_a:liver_vehctrl_rep2_a))*100) %>%
  select(ProteinID, `Set1 Peptides`, `Set2 Peptides`, heart_EGF_rep1_a:liver_vehctrl_rep3_b) %>%
  filter(!(rowSums(select(., heart_EGF_rep2_b:liver_vehctrl_rep3_b))==0)) %>%
  filter(!(rowSums(select(., heart_EGF_rep1_a:liver_vehctrl_rep2_a))==0)) %>%
  na.omit() %>%
  replace(is.na(.), 0) %>%
  #unite(protein_site, ProteinID, site_position) %>%
  filter(ProteinID %in% mean_normalized_scaled_total$ProteinID) %>% # these are already filtered as desired
  #filter(ProteinID %in% zscores_total$ProteinID | `Set1 Peptides` > 3 | `Set2 Peptides` > 3) %>%
  select(-c( `Set1 Peptides`, `Set2 Peptides`)) %>%
  gather(-ProteinID, key = "sample", value = "scaled_value") %>%
  mutate(plex = c(rep("a", nrow(.)/2), rep("b", nrow(.)/2))) %>%
  group_by(ProteinID, plex) %>%
  mutate(entropy_piece = calculate_total_entropy_piece(scaled_value)) %>% # calculate entropy for each protein within each plex
  #mutate(entropy_piece = sum(scaled_value))
  ungroup() %>%
  na.omit() %>%
  select(-c(scaled_value, plex)) %>%
  spread(key = sample, value = entropy_piece) %>%
  replace(is.na(.), 0) %>%
  mutate(across(heart_EGF_rep1_a:lung_vehctrl_rep3_b), is.numeric(.))

total_entropy_pieces <- total_entropy_pieces %>%
  mutate(entropy_a = rowSums(select(total_entropy_pieces,contains("_a")), na.rm = TRUE),
         entropy_b = rowSums(select(total_entropy_pieces,contains("_b")), na.rm = TRUE)) %>%
  arrange(by = entropy_a) 

# merge entropy values with signfiicant phospho fold-changes from MSstats
entropied_scale_total <- inner_join(select(total_entropy_pieces, c(ProteinID,  
                                                                   entropy_a,
                                                                   entropy_b)), scaled_total) %>%
  select(-c("Group ID", "Set1~rq_126_sn scaled":"Set2~rq_134n_sn scaled")) %>%
  mutate(entropy = rowMeans(select(., entropy_a:entropy_b))) %>%
  arrange(entropy) %>%
  mutate(heart_sig = ifelse(ProteinID %in% heart_sig_phospho_fold_changes_proteins, "sig", ifelse(ProteinID %in% heart_all_proteins, "n.s.", "n.d.")),
         lung_sig = ifelse(ProteinID  %in% lung_sig_phospho_fold_changes_proteins, "sig", ifelse(ProteinID %in% lung_all_proteins, "n.s.", "n.d.")),
         kidney_sig = ifelse(ProteinID  %in% kidney_sig_phospho_fold_changes_proteins, "sig", ifelse(ProteinID %in% kidney_all_proteins, "n.s.", "n.d.")),
         liver_sig = ifelse(ProteinID  %in% liver_sig_phospho_fold_changes_proteins, "sig", ifelse(ProteinID %in% liver_all_proteins, "n.s.", "n.d."))) %>%
  select(ProteinID, heart_sig:liver_sig,
         heart_EGF_rep1_a:liver_vehctrl_rep3_b, entropy)

# pull out only tissue-specific protein IDs
sig_total_entropy_IDs <- entropied_scale_total %>%
  filter(entropy < 2.58) %>%
  pull(ProteinID)

pdf(paste0(figure_path, "total_entropy_density_plot.pdf"))
ggplot(entropied_scale_total, aes(x = entropy)) + geom_density() + 
  theme_bw() + geom_vline(xintercept = 2.58, linetype = "dashed") 
dev.off()

# figure out which tissue the protein is enriched in
characterize_total_sig_entropy <- total %>% # first scale each plex
  mutate_at(vars(heart_EGF_rep2_b:liver_vehctrl_rep3_b), ~.x/rowSums(select(total, heart_EGF_rep2_b:liver_vehctrl_rep3_b))*100) %>%
  mutate_at(vars(heart_EGF_rep1_a:liver_vehctrl_rep2_a), ~.x/rowSums(select(total, heart_EGF_rep1_a:liver_vehctrl_rep2_a))*100) %>%
  select(ProteinID, heart_EGF_rep1_a:liver_vehctrl_rep3_b) %>%
  filter(!(rowSums(select(., heart_EGF_rep2_b:liver_vehctrl_rep3_b))==0)) %>%
  filter(!(rowSums(select(., heart_EGF_rep1_a:liver_vehctrl_rep2_a))==0)) %>%
  filter(ProteinID %in% sig_total_entropy_IDs) %>% # filter for only tissue-enriched proteins
  select(ProteinID, heart_EGF_rep1_a:liver_vehctrl_rep3_b) %>%
  gather(-ProteinID, key = sample, value = scaled_value) %>%
  separate(sample, into = c("tissue", "treatment", "rep", "plex")) %>%
  group_by(ProteinID, tissue, plex) %>% 
  summarise(mean = sum(scaled_value)) %>% # within each plex, the mean TMT value for a tissue must exceed 25 for the protein to be enriched in that tissue
  mutate(heart_specificity = ifelse(tissue == "heart" & mean > 25, "spec", "n.s."),
         lung_specificity = ifelse(tissue == "lung" & mean > 25, "spec", "n.s."),
         kidney_specificity = ifelse(tissue == "kidney" & mean > 25, "spec", "n.s."), 
         liver_specificity = ifelse(tissue == "liver" & mean > 25, "spec", "n.s.")) %>%
  arrange(ProteinID) 

get_total_sig_entropy_inconsistencies <- characterize_total_sig_entropy %>% 
  ungroup() %>%
  group_by(tissue) %>%
  ungroup() %>%
  mutate(spec = rowSums(. == "spec")) %>%
  group_by(ProteinID, tissue) %>%
  summarise(spec_sum = sum(spec)) %>%
  arrange(desc(spec_sum)) %>%
  filter(!(spec_sum == 0 | spec_sum == 2)) # anything that is only tissue-specific in a tissue in one plex is "inconsistent"

filtered_total_sig_entropies <- characterize_total_sig_entropy %>%
  ungroup() %>%
  filter(!(ProteinID %in% get_total_sig_entropy_inconsistencies$ProteinID)) %>% # filter out inconsistencies
  select(-c(plex, mean)) %>%
  distinct()  %>%
  # filter(!(rowSums(. == "n.s.") == 8)) 
  gather(-c(ProteinID, tissue), key = "name", value = "specificity") %>%
  select(-name) %>%
  filter(specificity == "spec") %>%
  pivot_wider(names_from = tissue, values_from = specificity) %>%
  replace(is.na(.), "n.s.") %>%
  rename_with(~ paste0(., "_spec"), -ProteinID)

sig_entropy_scale_total <- inner_join(filtered_total_sig_entropies, entropied_scale_total) # make entropied_scale_total include only significant proteins

sig_entropy_sig_scale_total <- sig_entropy_scale_total %>% # only tissue-specific proteins significantly changed in response to Egf
  filter(heart_sig == "sig" | kidney_sig == "sig" | lung_sig == "sig" | liver_sig == "sig") 



# percentage of tissue-enriched sites 

total_tissue_spec_df <- data.frame(cbind(c("specific", "non-specific"), 
                                         c(nrow(sig_entropy_scale_total), nrow(filter(entropied_scale_total, !(ProteinID %in% sig_entropy_scale_total$ProteinID))))))

total_tissue_spec_df$X2 = as.numeric(total_tissue_spec_df$X2)


pdf(paste0(figure_path, "total_percent_tissue_enriched.pdf"))
ggplot(total_tissue_spec_df, aes(x = X1, y = X2, fill = "darkblue")) +
  geom_bar(stat = "identity") + labs(x = "Tissue Specificity", y = "Number of Proteins") + 
  scale_fill_manual(values = "darkblue") + theme_classic() +
  annotate("text", x = 2, y = 2000, 
           label = paste0(round(100*nrow(sig_entropy_scale_total)/(nrow(entropied_scale_total)), 2), "%"))  
dev.off()

#sum(sig_entropy_scale_total$heart_spec == "spec")/nrow(sig_entropy_scale_total)
#sum(sig_entropy_scale_total$lung_spec == "spec")/nrow(sig_entropy_scale_total) 
#sum(sig_entropy_scale_total$kidney_spec == "spec")/nrow(sig_entropy_scale_total)
#sum(sig_entropy_scale_total$liver_spec == "spec")/nrow(sig_entropy_scale_total)


#sum(sig_entropy_scale_total$heart_spec == "spec")/nrow(sig_entropy_scale_total)
#sum(sig_entropy_scale_total$lung_spec == "spec")/nrow(sig_entropy_scale_total)
#sum(sig_entropy_scale_total$kidney_spec == "spec")/nrow(sig_entropy_scale_total)
#sum(sig_entropy_scale_total$liver_spec == "spec")/nrow(sig_entropy_scale_total)

#nrow(entropied_scale_total %>% filter(!heart_sig == "n.d."))
#nrow(entropied_scale_total %>% filter(!lung_sig == "n.d."))
#nrow(entropied_scale_total %>% filter(!kidney_sig == "n.d."))
#nrow(entropied_scale_total %>% filter(!liver_sig == "n.d."))

# UpSet plot

tissue_total_specificity_upset <-list(heart = sig_entropy_scale_total %>% filter(heart_spec == "spec") %>% pull(ProteinID),
                                      lung = sig_entropy_scale_total %>% filter(lung_spec == "spec") %>% pull(ProteinID),
                                      kidney = sig_entropy_scale_total %>% filter(kidney_spec == "spec") %>% pull(ProteinID),
                                      liver = sig_entropy_scale_total %>% filter(liver_spec == "spec") %>% pull(ProteinID))
pdf(paste0(figure_path, "total_tissue_specificity_upset.pdf"))
upset(fromList(tissue_total_specificity_upset), order.by = "degree",
      matrix.color = "black", sets.bar.color = "black", main.bar.color = "black")
dev.off()

# entropy barplots
make_total_tissue_entropy_sig_barplots("heart")
make_total_tissue_entropy_sig_barplots("lung")
make_total_tissue_entropy_sig_barplots("kidney")
make_total_tissue_entropy_sig_barplots("liver")


# Phospho part
phospho <- read_csv("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/2022-05-16_Bea_FinalMousePho.csv") %>%
  rename(ProteinID = `Protein Id`, 
         site_position = `Site Position`,
         heart_EGF_rep2_b = "Set2~126_sn_sum",
         heart_EGF_rep3_b = "Set2~127n_sn_sum",
         lung_EGF_rep2_b = "Set2~127c_sn_sum",	
         lung_EGF_rep3_b = "Set2~128n_sn_sum",
         kidney_EGF_rep2_b = "Set2~128c_sn_sum",
         kidney_EGF_rep3_b = "Set2~129n_sn_sum",	
         liver_EGF_rep2_b = "Set2~129c_sn_sum",	
         liver_EGF_rep3_b = "Set2~130n_sn_sum",	
         heart_vehctrl_rep1_b = "Set2~130c_sn_sum",
         heart_vehctrl_rep3_b = "Set2~131_sn_sum",	
         lung_vehctrl_rep1_b = "Set2~131c_sn_sum",
         lung_vehctrl_rep3_b = "Set2~132n_sn_sum",	
         kidney_vehctrl_rep1_b = "Set2~132c_sn_sum",
         kidney_vehctrl_rep3_b = "Set2~133n_sn_sum",
         liver_vehctrl_rep1_b = "Set2~133c_sn_sum",
         liver_vehctrl_rep3_b = "Set2~134n_sn_sum",
         
         heart_EGF_rep1_a = "Set1~126_sn_sum",
         heart_EGF_rep2_a = "Set1~127n_sn_sum",
         lung_EGF_rep1_a = "Set1~127c_sn_sum",	
         lung_EGF_rep2_a = "Set1~128n_sn_sum",
         kidney_EGF_rep1_a = "Set1~128c_sn_sum",
         kidney_EGF_rep2_a = "Set1~129n_sn_sum",	
         liver_EGF_rep1_a = "Set1~129c_sn_sum",	
         liver_EGF_rep2_a = "Set1~130n_sn_sum",	
         heart_vehctrl_rep1_a = "Set1~130c_sn_sum",
         heart_vehctrl_rep2_a = "Set1~131_sn_sum",	
         lung_vehctrl_rep1_a = "Set1~131c_sn_sum",
         lung_vehctrl_rep2_a = "Set1~132n_sn_sum",	
         kidney_vehctrl_rep1_a = "Set1~132c_sn_sum",
         kidney_vehctrl_rep2_a = "Set1~133n_sn_sum",
         liver_vehctrl_rep1_a = "Set1~133c_sn_sum",
         liver_vehctrl_rep2_a = "Set1~134n_sn_sum")

scaled_phospho <- scale_data(phospho) 


# Figure 6

scaled_phospho_for_venn <- scaled_phospho %>%
  mutate(plex1_sums = rowSums(select(scaled_phospho, heart_EGF_rep1_a:liver_vehctrl_rep2_a), na.rm = TRUE),
         plex2_sums = rowSums(select(scaled_phospho, heart_EGF_rep2_b:liver_vehctrl_rep3_b), na.rm = TRUE))


VennDiagram::venn.diagram(x = list((scaled_phospho_for_venn %>% filter(plex1_sums > 0) %>% pull(site_id)), (scaled_phospho_for_venn %>% filter(plex2_sums > 0)) %>% pull(site_id)), category.names = c("Plex 1", "Plex 2"),
                          filename = paste0(figure_path, 'shared_phospho_proteins_full_tissue.png'), output = TRUE, fill = c("lightgreen", "darkgoldenrod1"))


zscores_phospho <- scaled_phospho %>%
  #select(ProteinID, site_position, heart_EGF_rep2_b:liver_vehctrl_rep2_a, -"Set1~num_quant") %>%
  select(ProteinID, site_position, heart_EGF_rep2_b, lung_EGF_rep2_b, kidney_EGF_rep2_b, liver_EGF_rep2_b,
         heart_vehctrl_rep1_b, lung_vehctrl_rep1_b, kidney_vehctrl_rep1_b, liver_vehctrl_rep1_b,
         heart_EGF_rep2_a, lung_EGF_rep2_a, kidney_EGF_rep2_a, liver_EGF_rep2_a,
         heart_vehctrl_rep1_a, lung_vehctrl_rep1_a, kidney_vehctrl_rep1_a, liver_vehctrl_rep1_a) %>%
  filter(!(rowSums(select(., heart_EGF_rep2_b:liver_vehctrl_rep1_b))==0)) %>%
  filter(!(rowSums(select(., heart_EGF_rep2_a:liver_vehctrl_rep1_a))==0)) %>%
  unite(protein_site, ProteinID, site_position) %>%
  gather(-protein_site, key = "sample", value = "scaled_value") %>% 
  mutate(plex = c(rep("b", nrow(.)/2), rep("a", nrow(.)/2))) %>%
  group_by(protein_site, plex) %>%
  mutate(zscore = zscorer(scaled_value)) %>%
  ungroup() %>%
  na.omit() %>%
  group_by(protein_site) %>%
  summarise(cor = tryCatch(
    {cor.test(zscore[plex == "b"], zscore[plex == "a"], 
              na.action = "na.exclude")$estimate},
    error = function(e){
      0
    }
    
  )) %>%
  filter(cor > 0.9)

mean_normalized_phospho <- phospho %>%
  rowwise() %>%
  mutate(set2_mean = mean(c(heart_EGF_rep2_b, lung_EGF_rep2_b, kidney_EGF_rep2_b, liver_EGF_rep2_b,
                            heart_vehctrl_rep1_b, lung_vehctrl_rep1_b, kidney_vehctrl_rep1_b, liver_vehctrl_rep1_b)),
         set1_mean = mean(c(heart_EGF_rep2_a, lung_EGF_rep2_a, kidney_EGF_rep2_a, liver_EGF_rep2_a,
                            heart_vehctrl_rep1_a, lung_vehctrl_rep1_a, kidney_vehctrl_rep1_a, liver_vehctrl_rep1_a))) %>%
  ungroup() %>%
  mutate_at(vars(heart_EGF_rep2_b:liver_vehctrl_rep3_b), ~.x/set2_mean) %>%
  mutate_at(vars(heart_EGF_rep1_a:liver_vehctrl_rep2_a), ~.x/set1_mean) %>%
  rowwise() %>%
  mutate(heart_EGF_rep2 = mean(c(heart_EGF_rep2_a, heart_EGF_rep2_b)),
         lung_EGF_rep2 = mean(c(lung_EGF_rep2_a, lung_EGF_rep2_b)),
         kidney_EGF_rep2 = mean(c(kidney_EGF_rep2_a, kidney_EGF_rep2_b)),
         liver_EGF_rep2 = mean(c(liver_EGF_rep2_a, liver_EGF_rep2_b)),
         heart_vehctrl_rep1 = mean(c(heart_vehctrl_rep1_a, heart_vehctrl_rep1_b)),
         lung_vehctrl_rep1 = mean(c(lung_vehctrl_rep1_a, lung_vehctrl_rep1_b)),
         kidney_vehctrl_rep1 = mean(c(kidney_vehctrl_rep1_a, kidney_vehctrl_rep1_b)),
         liver_vehctrl_rep1 = mean(c(liver_vehctrl_rep1_a, liver_vehctrl_rep1_b))) %>%
  ungroup() %>%
  unite(protein_site, ProteinID, site_position) %>%
  filter(protein_site %in% zscores_phospho$protein_site) %>%
  select(-contains(c("EGF_rep2_a", "EGF_rep2_b", "vehctrl_rep1_a", "vehctrl_rep1_b"))) %>%
  rename_with(~str_replace(., '_a', '')) %>%
  rename_with(~str_replace(., '_b', '')) %>%
  select(protein_site, gene_symbol, Motif, sequence, heart_EGF_rep1, heart_EGF_rep2,
         heart_EGF_rep3, lung_EGF_rep1, lung_EGF_rep2, lung_EGF_rep3, kidney_EGF_rep1, kidney_EGF_rep2, kidney_EGF_rep3,
         liver_EGF_rep1, liver_EGF_rep2, liver_EGF_rep3, heart_vehctrl_rep1, heart_vehctrl_rep2, heart_vehctrl_rep3,
         lung_vehctrl_rep1, lung_vehctrl_rep2, lung_vehctrl_rep3, kidney_vehctrl_rep1, kidney_vehctrl_rep2, kidney_vehctrl_rep3,
         liver_vehctrl_rep1, liver_vehctrl_rep2, liver_vehctrl_rep3) %>%
  na.omit()

mean_normalized_scaled_phospho <- mean_normalized_phospho %>%
  mutate_at(vars(heart_EGF_rep1:liver_vehctrl_rep3), ~.x/rowSums(select(mean_normalized_phospho, heart_EGF_rep1:liver_vehctrl_rep3))*100) 

make_PCA(mean_normalized_scaled_phospho, "phospho")
make_dendrogram(mean_normalized_scaled_phospho, "phospho")


make_boxplot_baseline_phospho("sp|P63085|MK01_MOUSE_185", 3)
make_boxplot_baseline_phospho("sp|Q63844|MK03_MOUSE_205", 3)
make_boxplot_baseline_phospho("sp|Q01279|EGFR_MOUSE_1172", 3)
make_boxplot_baseline_phospho("sp|Q01279|EGFR_MOUSE_1166", 4)

make_boxplot_stimulated_phospho("sp|P63085|MK01_MOUSE_185", 4.5)
make_boxplot_stimulated_phospho("sp|Q01279|EGFR_MOUSE_1172", 4.5)
make_boxplot_stimulated_phospho("sp|Q63844|MK03_MOUSE_205", 4.5)
make_boxplot_stimulated_phospho("sp|Q01279|EGFR_MOUSE_1166", 4)


all_tissue_sig_phospho <- c(heart_sig_phospho_fold_changes,
                            lung_sig_phospho_fold_changes,
                            kidney_sig_phospho_fold_changes,
                            liver_sig_phospho_fold_changes)
phospho_for_sig_heat_map <- mean_normalized_scaled_phospho %>%
  filter(protein_site %in% all_tissue_sig_phospho) %>%
  separate(protein_site, into = c("ProteinID", "mouse", "site_position"), sep = "_") %>%
  unite(ProteinID, ProteinID, mouse) %>%
  select(ProteinID, heart_vehctrl_rep1:liver_vehctrl_rep3,
         heart_EGF_rep1:liver_EGF_rep3) %>%
  arrange(ProteinID) 
for_phospho_heat_map <- data.frame(lapply(select(phospho_for_sig_heat_map, heart_vehctrl_rep1:liver_EGF_rep3), as.numeric))
scaled_phospho_for_heat_map <- data.frame(t(apply(for_phospho_heat_map, 1, zscorer)))
rownames(scaled_phospho_for_heat_map) <- make.names(phospho_for_sig_heat_map$ProteinID, unique = TRUE)

scaled_phospho_for_heat_map <- scaled_phospho_for_heat_map %>%
  filter(!rowSums(.) == 0) %>%
  distinct() %>%
  filter_all(all_vars(is.finite(.)))

draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 0.5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
paletteLength <- 45
myColor <- colorRampPalette(c("darkblue","white", "darkred"))(paletteLength)
#myColor <- colorRampPalette(c("darkblue","cornflowerblue", "lightblue3", "lightcyan1", "lightyellow", "lightgoldenrod1","darkgoldenrod3","darkorange", "red3"))(paletteLength)
#paletteLength <- 7
#myColor <- colorRampPalette(c("darkblue","lightblue3", "lightcyan1", "lightyellow", "lightgoldenrod1","darkorange", "red3"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
bk1 <- c(seq(-4.5, -0.2, by =0.2), -0.10001)
bk2 <- c(0.10001, seq(.2,4.5, by = 0.2))
bk <- c(bk1, bk2)

pdf(paste0(figure_path, "full_tissue_sig_phosphosites.pdf"))
pheatmap(scaled_phospho_for_heat_map, cluster_cols = F, cluster_rows=T, 
         color=myColor, breaks=bk, border_color = "black", gaps_col = 12)
dev.off()







# Figure 7

library(lme4)
library(lmerTest)

pre_mean_normalized_total_for_model <- mean_normalized_scaled_total %>%
  select(ProteinID, heart_EGF_rep1:liver_vehctrl_rep3) %>%
  gather(-ProteinID, key = sample, value = total_abundance)

pre_phospho_for_model <- mean_normalized_scaled_phospho %>%
  separate(protein_site, into = c("ProteinID", "mouse", "site_position"), sep = "_") %>%
  unite(ProteinID, ProteinID, mouse) %>%
  select(ProteinID, site_position, heart_EGF_rep1, heart_EGF_rep2,
         heart_EGF_rep3, lung_EGF_rep1, lung_EGF_rep2, lung_EGF_rep3, kidney_EGF_rep1, kidney_EGF_rep2, kidney_EGF_rep3,
         liver_EGF_rep1, liver_EGF_rep2, liver_EGF_rep3, heart_vehctrl_rep1, heart_vehctrl_rep2, heart_vehctrl_rep3,
         lung_vehctrl_rep1, lung_vehctrl_rep2, lung_vehctrl_rep3, kidney_vehctrl_rep1, kidney_vehctrl_rep2, kidney_vehctrl_rep3,
         liver_vehctrl_rep1, liver_vehctrl_rep2, liver_vehctrl_rep3) %>%
  gather(-c(ProteinID, site_position), key = sample, value = phospho_abundance) %>%
  left_join(., pre_mean_normalized_total_for_model, by = c("ProteinID", "sample")) %>%
  unite(protein_site, ProteinID, site_position) %>%
  separate(sample, into = c("tissue", "treatment", "rep"), sep = "_") 
#filter(protein_site ==  "sp|Q63ZW7|INADL_MOUSE_1590")


zscored_proteins <- unique(pre_phospho_for_model$protein_site)

make_linear_model_plots("sp|Q01279|EGFR_MOUSE_1166")
make_linear_model_plots("sp|Q01279|EGFR_MOUSE_1172")
make_linear_model_plots("sp|P63085|MK01_MOUSE_185")

