library("RColorBrewer")
library("pheatmap")
library("grid")
library(dendextend)
library(ggfortify)
library(UpSetR)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(tidyverse)
library(clusterProfiler)

# Functions
# get significant EGF-stimulated phosphosites from MSstats
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
    unite(protein_site, protein, site) %>%
    filter(!(significance == "n.s."))
  return(sig_phospho$protein_site)
}

# get protein ID for significant EGF-stimulated phosphosites from MSstats
get_tissue_sig_proteins <- function(tissue){
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
  return(sig_phospho$ProteinID)
}

# get all phosphosites detected in a tissue
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
    unite(protein_site, protein, site) 
  return(all_phospho$protein_site)
}

# get all significant EGF : EGF + MEKi phosphosites from MSstats
get_tissue_ERK_dependent <- function(tissue){
  phospho <- read_csv(paste0("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/final_MSstats_files/ Phospho_ProteinAdjusted_StatResults_ ", tissue, " .csv"))
  phospho <- phospho %>%
    select(Reference, PhosphoSite, paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"),
           paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"),
           paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi")) %>%
    rename(log2FC = paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"),
           qValue = paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"),
           significance = paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"))
  sig_phospho <- phospho %>%
    mutate(protein = Reference, site = str_extract(PhosphoSite, "[[:digit:]]+")) %>%
    unite(protein_site, protein, site) %>%
    filter(!(significance == "n.s."))
  return(sig_phospho$protein_site)
}




# scale TMT from Core

scale_data <- function(df){
  scaled <- df  %>%
   # select(-c("heart_4.28", "lung_4.28", "kidney_4.28", "liver_4.28")) %>% # filter out unrelated sample
    mutate_at(vars(heart_4.2:heart_4.26, heart_4.30), ~.x/rowSums(select(df, heart_4.2:heart_4.26, heart_4.30))*100) %>%
    mutate_at(vars(lung_4.2:lung_4.26, lung_4.30), ~.x/rowSums(select(df, lung_4.2:lung_4.26, lung_4.30))*100) %>%
    mutate_at(vars(kidney_4.2:kidney_4.26, kidney_4.30), ~.x/rowSums(select(df, kidney_4.2:kidney_4.26, kidney_4.30))*100) %>%
    mutate_at(vars(liver_4.2:liver_4.26, liver_4.30), ~.x/rowSums(select(df, liver_4.2:liver_4.26, liver_4.30))*100) %>%
    filter(rowMeans(.[,c("heart_4.2","lung_4.2", "kidney_4.2", "liver_4.2")], na.rm=TRUE) > 0) %>% # filter out 0 values
    filter(!(grepl("##", gene_symbol)))  # filter out reverse hits
}

to_make_summary_barplots <- function(tissue){
  selected_phospho <- scaled_phospho %>%
    filter(!is.na(rowMeans(select(., paste0(tissue, "_4.2"):paste0(tissue, "_4.30"))))) %>%
    mutate(gene = gene_symbol, site = `Site Position`, protein = `Protein Id`) %>%
    unite("gene_site", gene_symbol, `Site Position`) %>%
    mutate(site_a = site) 
  selected_total <- scaled_total %>%
    filter(!is.na(rowMeans(select(., paste0(tissue, "_4.2"):paste0(tissue, "_4.30")))))
  print(c(nrow(selected_phospho), length(unique(selected_phospho$protein)),
          length(unique(selected_total$`Protein Id`))))
  return(c(nrow(selected_phospho), length(unique(selected_phospho$protein)),
           length(unique(selected_total$`Protein Id`))))
}
make_summary_barplots <- function(type){
  names <- c("Heart", "Lung", "Kidney", "Liver")
  if(type == "peptide"){
    for_summary_barplot <- data.frame(cbind(names,
                                            c(summary_barplot_info[[1]][[1]],
                                              summary_barplot_info[[2]][[1]],
                                              summary_barplot_info[[3]][[1]],
                                              summary_barplot_info[[4]][[1]])))
  }else if(type == "phospho-protein"){
    for_summary_barplot <- data.frame(cbind(names,
                                            c(summary_barplot_info[[1]][[2]],
                                              summary_barplot_info[[2]][[2]],
                                              summary_barplot_info[[3]][[2]],
                                              summary_barplot_info[[4]][[2]])))
  }else if(type == "total"){
    for_summary_barplot <- data.frame(cbind(names,
                                            c(summary_barplot_info[[1]][[3]],
                                              summary_barplot_info[[2]][[3]],
                                              summary_barplot_info[[3]][[3]],
                                              summary_barplot_info[[4]][[3]])))
  }
  colnames(for_summary_barplot) <- c("Tissue", "of_interest")
  print(for_summary_barplot)
  for_summary_barplot$Tissue <- factor(for_summary_barplot$Tissue, levels = for_summary_barplot$Tissue)
  for_summary_barplot$of_interest <- as.numeric(for_summary_barplot$of_interest)
  print(for_summary_barplot)
  g <- ggplot(for_summary_barplot, aes(x = Tissue, y = of_interest, fill = Tissue)) +
    geom_bar(stat = "identity", fill = c("darkred", "blue", "purple", "orange")) +
    theme_classic()
  
  if(type == "peptide"){
    g <- g +  labs(title = "Number of Phospho-Peptides Detected in Each Tissue") +
      ylab("Number of Phospho-Peptides") + 
      scale_y_continuous(breaks = seq(0, 12000, by = 2000), limits = c(0,12000)) +
      coord_flip()
    
  }else if(type == "phospho-protein"){
    g <- g +  labs(title = "Number of Phospho-Proteins Detected in Each Tissue") +
      ylab("Number of Phospho-Proteins") + 
      scale_y_continuous(breaks = seq(0, 5000, by = 1000), limits = c(0,5000)) + 
      coord_flip()
  }else if(type == "total"){
    g <- g +  labs(title = "Total number of Proteins Detected in Each Tissue") +
      ylab("Number of Proteins") + 
      scale_y_continuous(breaks = seq(0, 10000, by = 2000), limits = c(0,10000)) + 
      coord_flip()
  }
  plot(g)
}


coef_of_variation <- function(x){
  cov <- sd(x)/mean(x)
}
cov_for_one_treatment_group <- function(tissue, df, first, second, third){
  a <- paste0(tissue, "_", first)
  b <- paste0(tissue, "_", second)
  c <- paste0(tissue, "_", third)
  subset <- df %>% select(c(all_of(a), all_of(b), all_of(c)))
  print(head(subset))
  cov_vector <- apply(subset, 1, coef_of_variation)
  return(cov_vector)
}
make_veh_cov_plot <- function(tissue, color){
  cov_veh_control_df_tissue <- cov_veh_control_df %>%
    select(c(tissue, protein_site)) %>%
    na.omit() 
  g <- ggplot(cov_veh_control_df_tissue, aes(!! rlang::sym(tissue))) +
    geom_histogram(fill = color, bins = 50) + xlim(0,2.0) + 
    xlab("Coefficient of Variation") + ylab("Count") + theme_classic() + 
    scale_y_continuous(breaks = seq(0, 2600, by = 200), limits = c(0, 2600)) +
    geom_vline(xintercept = 0.7, linetype = 2)
  pdf(paste0(figure_path, tissue, "_vehcontrol_cov_plot.pdf"),
      height = 4, width = 4)
  plot(g)
  dev.off()
  #return(g)
}
make_egf_cov_plot <- function(tissue, color){
  cov_egf_df_tissue <- cov_egf_df %>%
    select(c(tissue, protein_site)) %>%
    na.omit() 
  g <- ggplot(cov_egf_df_tissue, aes(!! rlang::sym(tissue))) +
    geom_histogram(fill = color, bins = 50) + xlim(0,2.0) + 
    xlab("Coefficient of Variation") + ylab("Count") + theme_classic() + 
    scale_y_continuous(breaks = seq(0, 2600, by = 200), limits = c(0, 2600)) +
    geom_vline(xintercept = 0.6, linetype = 2)
  pdf(paste0(figure_path, tissue, "_egf_cov_plot.pdf"),
      height = 4, width = 4)
  plot(g)
  dev.off()
  return(g)
}

make_PCA <- function(the_data){
  df <- df2 <- the_data %>%
    t()%>%
    data.frame() %>% select(Heart.EGF.r1:Liver.EGF.r3) %>% 
    na.omit() %>% t() %>% data.frame()
  df2$tissue <- c(rep("heart_EGF", 3), rep("lung_EGF", 3), rep("kidney_EGF", 3), rep("liver_EGF", 3))
  pca_res <- prcomp(df, center = TRUE, scale = TRUE)
  g <-  autoplot(pca_res, data = df2, colour = "tissue", frame = TRUE, label = TRUE) + theme_classic()
  pdf(paste0(figure_path,"single_PCA.pdf"))
  plot(g)
  dev.off()
}

volcano_plot <- function(tissue){
  df <- read_csv(paste0("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/final_MSstats_files/ Phospho_ProteinAdjusted_StatResults_ ", tissue, " .csv"))
  
  df <- df %>%
    select(Reference, paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl")) %>%
    rename(log2FC = paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           qValue = paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           significance = paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"))
  df$color <- 0
  df$color[df$significance == "n.s."] <- 1
  g= ggplot(df, aes(x=log2FC, y=-log10(qValue), group = color)) + 
    geom_point(aes(color=as.factor(color)), show.legend = FALSE)  +  
    labs(x = "log2(fc)", y = "-log10(qvalue)", title = str_to_sentence(tissue)) +
    scale_x_continuous(limits = c(-2.5, 3.5), breaks = seq(-2.5, 3.5, by = 0.5)) + 
    scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 1)) + 
    geom_vline(xintercept = 0.58, linetype = "dashed") + 
    geom_vline(xintercept = -0.58, linetype = "dashed") + 
    geom_hline(yintercept = 1.3, linetype = "dashed")  +
    annotate("text", x = 2, y = 4, 
             label = paste0(round(100*sum(df$qValue < 0.05 & 
                                            df$log2FC > 0.58, 
                                          na.rm = TRUE)/nrow(df),1), "%")) + 
    annotate("text", x = -2, y = 4, 
             label = paste0(round(100*sum(df$qValue < 0.05 & 
                                            df$log2FC < -0.58, 
                                          na.rm = TRUE)/nrow(df),1), "%")) + 
    theme_classic() + scale_colour_manual(values = c("red", "black"))
  pdf(paste0(figure_path,tissue,"EGF_sig_volcano_plot.pdf"))
  plot(g)
  dev.off()
}

make_pairwise_scatterplot <- function(tissue1, tissue2){
  tissue_fold_change_table <- for_scatterplots %>%
    select(uniprot_site, gene_site, paste0(tissue1, "_EGF_means"), paste0(tissue2, "_EGF_means")) %>%
    na.omit()
  print(head(tissue_fold_change_table))
  colnames(tissue_fold_change_table) = c("uniprot_site", "gene_site", "X", "Y")
  g =  ggplot(tissue_fold_change_table, aes(x=X, y=Y)) + 
    geom_point()  +  geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) +
    theme_bw() +  scale_x_continuous(limits = c(-2, 4)) +
    scale_y_continuous(limits = c(-2, 4)) + 
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = 2) + labs(x =paste0("log2(",tissue1," EGF-induced FC)"),
                                                                               y = paste0("log2(",tissue2," EGF-induced FC)")) +
    annotate("text", x = 2, y = 4, size = 10,
             label = paste0(round(100*sum(tissue_fold_change_table$X > 0 & 
                                            tissue_fold_change_table$Y > 0, 
                                          na.rm = TRUE)/nrow(tissue_fold_change_table),1), "%")) + 
    annotate("text", x = -1, y = 4, size = 10,
             label = paste0(round(100*sum(tissue_fold_change_table$X < 0 & 
                                            tissue_fold_change_table$Y > 0, 
                                          na.rm = TRUE)/nrow(tissue_fold_change_table),1), "%")) + 
    annotate("text", x = -1, y = -2, size = 10,
             label = paste0(round(100*sum(tissue_fold_change_table$X < 0 & 
                                            tissue_fold_change_table$Y < 0, 
                                          na.rm = TRUE)/nrow(tissue_fold_change_table),1), "%")) + 
    annotate("text", x = 2, y = -2, size = 10,
             label = paste0(round(100*sum(tissue_fold_change_table$X > 0 & 
                                            tissue_fold_change_table$Y < 0, 
                                          na.rm = TRUE)/nrow(tissue_fold_change_table),1), "%")) +
    annotate("text", x = -1, y = 2, size = 10, label =paste0("r = ",round(cor(tissue_fold_change_table$X, tissue_fold_change_table$Y),2)))
  pdf(paste0(figure_path, tissue1, "_", tissue2, "_pairwise EGF alone scatterplot_with_labels.pdf"))
  plot(g)
  dev.off()
}

do_GO_analysis_on_EGF <- function(tissue, type){
  phospho <- read_csv(paste0("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/final_MSstats_files/ Phospho_ProteinAdjusted_StatResults_ ", tissue, " .csv"))
  phospho <- phospho %>%
    select(Reference, paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl")) %>%
    rename(log2FC = paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           qValue = paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           significance = paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"))
  
  total_all_genes <- scaled_total %>%
    mutate(protein = `Protein Id`) %>%
    select(protein, paste0(tissue, "_4.2"): paste0(tissue, "_4.30")) %>%
    na.omit()
  print(c(dim(scaled_total), dim(total_all_genes)))
  phospho_all_genes <- phospho %>%
    mutate(protein = Reference) %>%
    na.omit()
  all_genes <- full_join(select(total_all_genes, protein), select(phospho_all_genes, protein))
  print(dim(all_genes))
  
  sig_phospho_all_genes <- phospho %>%
    mutate(protein = Reference) %>%
    filter(!(significance == "n.s."))
  View(sig_phospho_all_genes)
  sig_phospho_all_genes_sep <- sig_phospho_all_genes %>%
    mutate(firstsub = gsub("sp\\|", "", protein)) %>%
    mutate(firstsub = gsub("tr\\|", "", firstsub)) %>%
    mutate(uniprot = gsub("\\|.*","",firstsub))
  print(dim(sig_phospho_all_genes_sep))
  
  all_genes_sep <- all_genes %>%
    mutate(firstsub = gsub("sp\\|", "", protein)) %>%
    mutate(firstsub = gsub("tr\\|", "", firstsub)) %>%
    mutate(uniprot = gsub("\\|.*","",firstsub))
  
  gene <- unique(sig_phospho_all_genes_sep$uniprot)
  gene_list = sort(gene, decreasing = TRUE)
  print(length(gene))
  universe <- unique(all_genes_sep$uniprot)
  print(length(universe))
  
  if(type == "MF"){
    ego <- enrichGO(gene          = gene,
                    universe = universe, 
                    OrgDb         = org.Mm.eg.db,
                    keyType = "UNIPROT",
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
  }else if(type == "BP"){
    ego <- enrichGO(gene          = gene,
                    universe = universe, 
                    OrgDb         = org.Mm.eg.db,
                    keyType = "UNIPROT",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
  }else if(type == "CC"){
    ego <- enrichGO(gene          = gene,
                    universe = universe, 
                    OrgDb         = org.Mm.eg.db,
                    keyType = "UNIPROT",
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
  }
  
  
  #barplot(ggo, showCategory=20) 
  return(ego)
  
  
}

for_kinase_enrichment_analysis <- function(tissue, df, instruction){
  
  if(instruction == "egf_alone"){
    cleaned <- df %>%
      mutate(sequence = gsub("#", "*", sequence)) %>%
      mutate(sequence = gsub("\\.", "", sequence)) %>%
      mutate(heart_EGF_means = log(heart_EGF_means, 2),
             lung_EGF_means = log(lung_EGF_means, 2),
             kidney_EGF_means = log(kidney_EGF_means, 2),
             liver_EGF_means = log(liver_EGF_means, 2)) %>%
      mutate(heart_EGF_plus_minus_MEKi_fold_change = log(heart_EGF_plus_minus_MEKi_fold_change, 2),
             lung_EGF_plus_minus_MEKi_fold_change = log(lung_EGF_plus_minus_MEKi_fold_change, 2),
             kidney_EGF_plus_minus_MEKi_fold_change = log(kidney_EGF_plus_minus_MEKi_fold_change, 2),
             liver_EGF_plus_minus_MEKi_fold_change = log(liver_EGF_plus_minus_MEKi_fold_change, 2))
    selected <- cleaned %>%
      select(c(protein_site, sequence, paste0(tissue, "_EGF_means"))) %>%
      na.omit()
    write_csv(selected, paste0(path, tissue, "_EGF_alone_fold_changes.csv"))
  }else if(instruction == "egf_plus_minus_MEKi"){
    cleaned <- df %>%
      # filter(!(grepl("S#P|T#P", sequence))) %>%
      mutate(sequence = gsub("#", "*", sequence)) %>%
      mutate(sequence = gsub("\\.", "", sequence)) %>%
      mutate(heart_EGF_means = log(heart_EGF_means, 2),
             lung_EGF_means = log(lung_EGF_means, 2),
             kidney_EGF_means = log(kidney_EGF_means, 2),
             liver_EGF_means = log(liver_EGF_means, 2)) %>%
      mutate(heart_EGF_plus_minus_MEKi_fold_change = log(heart_EGF_plus_minus_MEKi_fold_change, 2),
             lung_EGF_plus_minus_MEKi_fold_change = log(lung_EGF_plus_minus_MEKi_fold_change, 2),
             kidney_EGF_plus_minus_MEKi_fold_change = log(kidney_EGF_plus_minus_MEKi_fold_change, 2),
             liver_EGF_plus_minus_MEKi_fold_change = log(liver_EGF_plus_minus_MEKi_fold_change, 2))
    selected <- cleaned %>%
      select(c(protein_site, sequence, paste0(tissue, "_EGF_plus_minus_MEKi_fold_change"))) %>%
      na.omit()
    write_csv(selected, paste0(path, tissue, "_EGF_plus_minus_MEKi_means.csv"))
  }
  return(head(selected))
}

make_EGF_versus_plusMEKi_scatterplot <- function(tissue){
  df <- for_scatterplots %>%
    select(protein_site, paste0(tissue, "_EGF_plus_MEKi_means"), paste0(tissue, "_EGF_means")) %>%
    na.omit() %>%
    rename(EGF_plus_MEKi_means = paste0(tissue, "_EGF_plus_MEKi_means"),
           EGF_means = paste0(tissue, "_EGF_means"))
  sig_sites <- get_tissue_ERK_dependent(tissue)
  EGF_sig <- get_tissue_sig(tissue)
  all_sig <- intersect(sig_sites, EGF_sig) # all sites that are ERK-dependent and EGF-dependent
  print(c(length(sig_sites), length(EGF_sig), length(all_sig)))
  df <- df %>%
    mutate(color = ifelse(protein_site %in% all_sig, 0, 1))
  print(df$EGF_plus_MEKi_means < df$EGF_means & 
          df$color == 0)
  View(df)
  g <- ggplot(df, aes(x=EGF_plus_MEKi_means, y=EGF_means), group = color ) + 
    geom_point(aes(color=as.factor(color)), show.legend = FALSE)  +  geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + scale_color_manual(values = c("red", "black")) +
    theme_bw() + scale_x_continuous(limits = c(-2, 4)) +
    scale_y_continuous(limits = c(-2, 4)) + 
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = 2) +
    annotate("text", x = 2, y = 4, size = 10,
             label = paste0(round(100*sum(df$EGF_plus_MEKi_means < df$EGF_means & 
                                            df$color == 0, 
                                          na.rm = TRUE)/length(EGF_sig),1), "%")) +
    annotate("text", x = 2, y = -2, size = 10,
             label = paste0(round(100*sum(df$EGF_plus_MEKi_means > df$EGF_means & 
                                            df$color == 0, 
                                          na.rm = TRUE)/length(EGF_sig),1), "%"))
  
  pdf(paste0(figure_path,tissue, "_EGF_EGFplusMEKi_scatterplot.pdf"))
  plot(g) 
  dev.off()
  #print(nrow(df %>% filter(EGF_plus_MEKi_means < EGF_means & color == 0)))
  #print(nrow(df))
}

do_GO_analysis_on_EGF_MEKi <- function(tissue, sig_proteins){
  phospho <- read_csv(paste0("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/final_MSstats_files/ Phospho_ProteinAdjusted_StatResults_ ", tissue, " .csv"))
  phospho <- phospho %>%
    select(Reference, PhosphoSite, paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"),
           paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"),
           paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi")) %>%
    rename(log2FC = paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"),
           qValue = paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"),
           significance = paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"EGFMEKi"))
  
  total_all_genes <- scaled_total %>%
    mutate(protein = `Protein Id`) %>%
    select(protein, paste0(tissue, "_4.2"): paste0(tissue, "_4.30")) %>%
    na.omit()
  print(c(dim(scaled_total), dim(total_all_genes)))
  phospho_all_genes <- phospho %>%
    mutate(protein = Reference) %>%
    na.omit()
  all_genes <- full_join(select(total_all_genes, protein), select(phospho_all_genes, protein))
  print(dim(all_genes))
  
  sig_phospho_all_genes <- phospho %>%
    mutate(protein = Reference) %>%
    filter(!(significance == "n.s.")) %>%
    filter(protein %in% sig_proteins) # look at only EGF + MEKi -significant sites that are also EGF-significant
  
  sig_phospho_all_genes_sep <- sig_phospho_all_genes %>%
    mutate(firstsub = gsub("sp\\|", "", protein)) %>%
    mutate(firstsub = gsub("tr\\|", "", firstsub)) %>%
    mutate(uniprot = gsub("\\|.*","",firstsub))
  print(dim(sig_phospho_all_genes_sep))
  
  all_genes_sep <- all_genes %>%
    mutate(firstsub = gsub("sp\\|", "", protein)) %>%
    mutate(firstsub = gsub("tr\\|", "", firstsub)) %>%
    mutate(uniprot = gsub("\\|.*","",firstsub))
  
  gene <- unique(sig_phospho_all_genes_sep$uniprot)
  gene_list = sort(gene, decreasing = TRUE)
  print(length(gene))
  universe <- unique(all_genes_sep$uniprot)
  print(head(universe))
  
  ego <- enrichGO(gene          = gene,
                  universe = universe, 
                  OrgDb         = org.Mm.eg.db,
                  keyType = "UNIPROT",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  #barplot(ggo, showCategory=20) 
  return(ego)
  
  
}

get_tissue_all_MSstats <- function(tissue){
  phospho <- read_csv(paste0("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/final_MSstats_files/ Phospho_ProteinAdjusted_StatResults_ ", tissue, " .csv"))
  phospho <- phospho %>%
    select(Reference, PhosphoSite, paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl")) %>%
    rename(log2FC = paste0("log2FC_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           qValue = paste0("q.val_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"),
           significance = paste0("sig_",str_to_sentence(tissue),"EGF-",str_to_sentence(tissue),"VehCtl"))
  all_phospho <- phospho %>%
    mutate(protein = Reference, site = sub('.', '', PhosphoSite)) %>%
    unite(protein_site, protein, site) %>%
    #filter(!(significance == "n.s."))
    return(all_phospho)
}



# Choose path where files are located and where to save files to
path="/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/"
figure_path="/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/R_figures_for_paper/"

# Read in files- output from Comet-based in-house software

og_phospho <- read_csv(paste0(path, "2023-01-24_Bea_5mouseTissuesInOne_ProPho.csv")) %>%
  select(c("Protein Id", "gene_symbol", "Site Position", "sequence", "Motif", "heart_4.2":"heart_4.30", 
           "lung_4.2":"lung_4.30", "kidney_4.2":"kidney_4.30", "liver_4.2":"liver_4.30")) %>%
  slice(-1) %>%
  mutate_at(vars(heart_4.2:liver_4.30), as.numeric) 


og_total <- read_csv(paste0(path,"total 2023-01-24_Bea_5mouseTissuesInOne_ProPho.csv")) %>%
  select(c("Protein Id", "Gene Symbol", "heart_4.2":"heart_4.30", 
           "lung_4.2":"lung_4.30", "kidney_4.2":"kidney_4.30", "liver_4.2":"liver_4.30")) %>%
  slice(-1) %>%
  rename(gene_symbol = "Gene Symbol") %>%
  mutate_at(vars(heart_4.2:liver_4.30), as.numeric) 

all_tissues <- list("heart", "lung", "kidney", "liver")


# Read in corresponding MSstats significance outputs

# all proteins detected in tissues
heart_all <- get_tissue_all("heart")
lung_all <- get_tissue_all("lung")
kidney_all <- get_tissue_all("kidney")
liver_all <- get_tissue_all("liver")

# sig phosphopeptide fold-changes in response to EGF
heart_sig_phospho_fold_changes <- get_tissue_sig("heart")
lung_sig_phospho_fold_changes <- get_tissue_sig("lung")
kidney_sig_phospho_fold_changes <- get_tissue_sig("kidney")
liver_sig_phospho_fold_changes <- get_tissue_sig("liver")

# sig fold-changes between EGF and EGF + MEKi samples
heart_sig_ERK_dependent_phospho_fold_changes <- get_tissue_ERK_dependent("heart")
lung_sig_ERK_dependent_phospho_fold_changes <- get_tissue_ERK_dependent("lung")
kidney_sig_ERK_dependent_phospho_fold_changes <- get_tissue_ERK_dependent("kidney")
liver_sig_ERK_dependent_phospho_fold_changes <- get_tissue_ERK_dependent("liver")

# sig phosphoprotein fold-changes in response to EGF
heart_sig_phospho_proteins <- get_tissue_sig_proteins("heart")
lung_sig_phospho_proteins <- get_tissue_sig_proteins("lung")
kidney_sig_phospho_proteins <- get_tissue_sig_proteins("kidney")
liver_sig_phospho_proteins <- get_tissue_sig_proteins("liver")


# Scale, filter data from in-house software
scaled_phospho <- scale_data(og_phospho)
scaled_total <- scale_data(og_total)


# Figure 1

# Make barplots summarizing [phospho]peptide and [phospho]protein numbers
summary_barplot_info <- lapply(all_tissues, to_make_summary_barplots)

pdf(paste0(figure_path, "phosphoprotein_summary_barplot.pdf"))
make_summary_barplots("phospho-protein")
dev.off()
pdf(paste0(figure_path, "phosphopeptide_summary_barplot.pdf"))
make_summary_barplots("peptide")
dev.off()
pdf(paste0(figure_path, "total_protein_summary_barplot.pdf"))
make_summary_barplots("total")
dev.off()


# Make heat maps
for_pre_heat_map <- scaled_phospho %>%  
  mutate(heart_baseline_means = rowMeans(select(., heart_4.6, heart_4.16, heart_4.26), na.rm=TRUE),
         lung_baseline_means = rowMeans(select(., lung_4.6, lung_4.16, lung_4.26), na.rm=TRUE),
         kidney_baseline_means = rowMeans(select(., kidney_4.6, kidney_4.16, kidney_4.26), na.rm=TRUE),
         liver_baseline_means = rowMeans(select(., liver_4.6, liver_4.16, liver_4.26), na.rm=TRUE)) %>%
  mutate(heart_EGF_rep1 = heart_4.2/heart_baseline_means,
         heart_EGF_rep2 = heart_4.12/heart_baseline_means,
         heart_EGF_rep3 = heart_4.22/heart_baseline_means,
         lung_EGF_rep1 = lung_4.2/lung_baseline_means,
         lung_EGF_rep2 = lung_4.12/lung_baseline_means,
         lung_EGF_rep3 = lung_4.22/lung_baseline_means,
         kidney_EGF_rep1 = kidney_4.2/kidney_baseline_means,
         kidney_EGF_rep2 = kidney_4.12/kidney_baseline_means,
         kidney_EGF_rep3 = kidney_4.22/kidney_baseline_means,
         liver_EGF_rep1 = liver_4.2/liver_baseline_means,
         liver_EGF_rep2 = liver_4.12/liver_baseline_means,
         liver_EGF_rep3 = liver_4.22/liver_baseline_means,
         
         heart_EGF_plus_MEKi_rep1 = heart_4.4/heart_baseline_means,
         heart_EGF_plus_MEKi_rep2 = heart_4.14/heart_baseline_means,
         heart_EGF_plus_MEKi_rep3 = heart_4.24/heart_baseline_means,
         lung_EGF_plus_MEKi_rep1 = lung_4.4/lung_baseline_means,
         lung_EGF_plus_MEKi_rep2 = lung_4.14/lung_baseline_means,
         lung_EGF_plus_MEKi_rep3 = lung_4.24/lung_baseline_means,
         kidney_EGF_plus_MEKi_rep1 = kidney_4.4/kidney_baseline_means,
         kidney_EGF_plus_MEKi_rep2 = kidney_4.14/kidney_baseline_means,
         kidney_EGF_plus_MEKi_rep3 = kidney_4.24/kidney_baseline_means,
         liver_EGF_plus_MEKi_rep1 = liver_4.4/liver_baseline_means,
         liver_EGF_plus_MEKi_rep2 = liver_4.14/liver_baseline_means,
         liver_EGF_plus_MEKi_rep3 = liver_4.24/liver_baseline_means,
         
         heart_MEKi_rep1 = heart_4.8/heart_baseline_means,
         heart_MEKi_rep2 = heart_4.18/heart_baseline_means,
         lung_MEKi_rep1 = lung_4.8/lung_baseline_means,
         lung_MEKi_rep2 = lung_4.18/lung_baseline_means,
         kidney_MEKi_rep1 = kidney_4.8/kidney_baseline_means,
         kidney_MEKi_rep2 = kidney_4.18/kidney_baseline_means,
         liver_MEKi_rep1 = liver_4.8/liver_baseline_means,
         liver_MEKi_rep2 = liver_4.18/liver_baseline_means) %>%
  #filter(gene_symbol == "Mapk3" | gene_symbol == "Mapk1") %>%
  #filter(gene_symbol == "Tsc2" | gene_symbol == "Tsc1") %>%
  #filter(gene_symbol == "Cdkn1b" | gene_symbol == "Cdk2") %>%
  filter(gene_symbol == "Egfr") %>%
  unite("gene_site", gene_symbol, "Site Position") 

EGF_alone_heat_map <- for_pre_heat_map %>%
  select(c(gene_site, heart_EGF_rep1:heart_EGF_rep3,
           lung_EGF_rep1:lung_EGF_rep3, 
           kidney_EGF_rep1:kidney_EGF_rep3,
           liver_EGF_rep1:liver_EGF_rep3))

for_heat_map <- EGF_alone_heat_map

scaled_for_heat_map <- data.frame(log2(for_heat_map[,2:ncol(for_heat_map)]))
rownames(scaled_for_heat_map) <- make.names(for_heat_map$gene_site, unique = TRUE)

draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 0.5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
paletteLength <- 46
myColor <- colorRampPalette(c("darkblue","white", "darkred"))(paletteLength)

bk1 <- c(seq(-4.6, -0.2, by =0.2), -0.10001)
bk2 <- c(0.10001, seq(.2,4.6, by = 0.2))
bk <- c(bk1, bk2)

pdf(paste0(figure_path, "EGF_alone_EGFR_heat_map.pdf"))
pheatmap(scaled_for_heat_map, cluster_cols = F, cluster_rows=F, 
         color=myColor, breaks=bk, na_col = "grey", border_color = "black",
         gaps_col = c(3, 6, 9))
dev.off()


#get coefficient of variation across vehicle controls
cov_veh_control_list <- lapply(all_tissues, cov_for_one_treatment_group,scaled_phospho,"4.6", "4.16", "4.26")
cov_veh_control_df <- as.data.frame(do.call(cbind, cov_veh_control_list))
colnames(cov_veh_control_df) <- unlist(all_tissues)
cov_veh_control_df <- cov_veh_control_df %>%
  mutate(protein = scaled_phospho$`Protein Id`, site = scaled_phospho$`Site Position`) %>%
  unite("protein_site", protein, site)

make_veh_cov_plot("heart", "darkred")
make_veh_cov_plot("lung", "blue")
make_veh_cov_plot("kidney", "purple")
make_veh_cov_plot("liver", "orange")

cov_veh_control_df_top_varying <- cov_veh_control_df %>%
  gather(key = "tissue", value = "coefficient", -protein_site) %>%
  replace(is.na(.), 0) %>%
  filter(coefficient > 0.7)

all_stable_veh_scale_phospho <- scaled_phospho %>% 
  mutate(protein = `Protein Id`, site = `Site Position`) %>%
  unite("protein_site", protein, site) %>%
  filter(!(protein_site %in% cov_veh_control_df_top_varying$protein_site)) %>%
  select(-protein_site)

write_csv(all_stable_veh_scale_phospho, paste0(path, "all_stable_veh_scale_phospho.csv"))
          
# get coefficient of variation across EGF
cov_egf_list <- lapply(all_tissues, cov_for_one_treatment_group, all_stable_veh_scale_phospho, "4.2", "4.12", "4.22")
cov_egf_df <- as.data.frame(do.call(cbind, cov_egf_list))
colnames(cov_egf_df) <- unlist(all_tissues)
cov_egf_df <- cov_egf_df %>%
  mutate(protein = all_stable_veh_scale_phospho$`Protein Id`, site = all_stable_veh_scale_phospho$`Site Position`) %>%
  unite("protein_site", protein, site)

make_egf_cov_plot("heart", "darkred")
make_egf_cov_plot("lung", "blue")
make_egf_cov_plot("kidney", "purple")
make_egf_cov_plot("liver", "orange")

cov_egf_df_top_varying <- cov_egf_df %>%
  gather(key = "tissue", value = "coefficient", -protein_site) %>%
  replace(is.na(.), 0) %>%
  filter(coefficient > 0.6)

stable_veh_egf_scale_phospho <- all_stable_veh_scale_phospho %>% 
  mutate(protein = `Protein Id`, site = `Site Position`) %>%
  unite("protein_site", protein, site) %>%
  filter(!(protein_site %in% cov_egf_df_top_varying$protein_site)) %>%
  select(-protein_site)

write_csv(stable_veh_egf_scale_phospho, paste0(path,"stable_veh_egf_scale_phospho.csv"))


# make dendrogram
fold_change_table_for_dendrogram <- stable_veh_egf_scale_phospho %>%  
  mutate(heart_baseline_means = rowMeans(.[,c("heart_4.6", "heart_4.16", "heart_4.26")], na.rm=TRUE),
         lung_baseline_means = rowMeans(.[,c("lung_4.6", "lung_4.16", "lung_4.26")], na.rm=TRUE),
         kidney_baseline_means = rowMeans(.[,c("kidney_4.6", "kidney_4.16", "kidney_4.26")], na.rm=TRUE),
         liver_baseline_means = rowMeans(.[,c("liver_4.6", "liver_4.16", "liver_4.26")], na.rm=TRUE)) %>%
  mutate(heart_EGF_rep1 = log(heart_4.2/heart_baseline_means + 0.01,2),
         heart_EGF_rep2 = log(heart_4.12/heart_baseline_means + 0.01,2),
         heart_EGF_rep3 = log(heart_4.22/heart_baseline_means + 0.01, 2),
         lung_EGF_rep1 = log(lung_4.2/lung_baseline_means + 0.01, 2),
         lung_EGF_rep2 = log(lung_4.12/lung_baseline_means + 0.01, 2),
         lung_EGF_rep3 = log(lung_4.22/lung_baseline_means + 0.01, 2),
         kidney_EGF_rep1 = log(kidney_4.2/kidney_baseline_means + 0.01, 2),
         kidney_EGF_rep2 = log(kidney_4.12/kidney_baseline_means + 0.01,2),
         kidney_EGF_rep3 = log(kidney_4.22/kidney_baseline_means + 0.01, 2),
         liver_EGF_rep1 = log(liver_4.2/liver_baseline_means + 0.01, 2),
         liver_EGF_rep2 = log(liver_4.12/liver_baseline_means + 0.01, 2),
         liver_EGF_rep3 = log(liver_4.22/liver_baseline_means + 0.01, 2)) %>%
  #na.omit() %>%
  #replace(is.na(.), 0) %>%
  select(c("heart_EGF_rep1":"liver_EGF_rep3")) %>%
  t() %>%
  data.frame()

#dist_mat <- dist(fold_change_table, method = 'euclidean')
#hclust_avg <- hclust(dist_mat, method = 'average')


rownames(fold_change_table_for_dendrogram) <- c("Heart EGF r1", "Heart EGF r2", "Heart EGF r3",
                                                "Lung EGF r1", "Lung EGF r2", "Lung EGF r3",
                                                "Kidney EGF r1", "Kidney EGF r2", "Kidney EGF r3",
                                                "Liver EGF r1", "Liver EGF r2", "Liver EGF r3")

dend <- fold_change_table_for_dendrogram %>% 
  dist() %>% 
  hclust(method = "ward.D2") %>% 
  as.dendrogram() 
dend %>% plot(horiz = TRUE)


heart_egf <- grepl("Heart EGF r", labels(dend))
lung_egf <- grepl("Lung EGF r", labels(dend))
kidney_egf <- grepl("Kidney EGF r", labels(dend))
liver_egf <- grepl("Liver EGF r", labels(dend))


pdf(paste0(figure_path, "egf_alone_tissue_clustering.pdf"), height = 25, width = 35)


dend %>%
  set("leaves_cex", 1) %>%
  set("leaves_pch", 15) %>%
  set("leaves_col", case_when(heart_egf~"darkred",
                              lung_egf~"blue",
                              kidney_egf~"purple",
                              liver_egf~"orange")) %>%
  set("labels_colors", case_when(heart_egf~"darkred",
                                 lung_egf~"blue",
                                 kidney_egf~"purple",
                                 liver_egf~"orange")) %>%
  #set("branches_lty", 5) %>%
  plot(horiz=TRUE, axes=FALSE)

dev.off()

make_PCA(fold_change_table_for_dendrogram)

# UpSet plots
#for each tissue, make a data frame of all detected data to make summary plots
total_heart_phospho_detected <- stable_veh_egf_scale_phospho %>%
  filter(!is.na(rowMeans(select(., heart_4.2:heart_4.30)))) %>%
  mutate(gene = gene_symbol, site = `Site Position`, protein = `Protein Id`) %>%
  unite("gene_site", gene_symbol, `Site Position`) %>%
  mutate(site_a = site) %>%
  unite("protein_site", `Protein Id`, site_a)
total_lung_phospho_detected <- stable_veh_egf_scale_phospho %>%
  filter(!is.na(rowMeans(select(., lung_4.2:lung_4.30)))) %>%
  mutate(gene = gene_symbol, site = `Site Position`, protein = `Protein Id`) %>%
  unite("gene_site", gene_symbol, `Site Position`) %>%
  mutate(site_a = site) %>%
  unite("protein_site", `Protein Id`, site_a)
total_kidney_phospho_detected <- stable_veh_egf_scale_phospho %>%
  filter(!is.na(rowMeans(select(., kidney_4.2:kidney_4.30)))) %>%
  mutate(gene = gene_symbol, site = `Site Position`, protein = `Protein Id`) %>%
  unite("gene_site", gene_symbol, `Site Position`) %>%
  mutate(site_a = site) %>%
  unite("protein_site", `Protein Id`, site_a)
total_liver_phospho_detected <- stable_veh_egf_scale_phospho %>%
  filter(!is.na(rowMeans(select(., liver_4.2:liver_4.30)))) %>%
  mutate(gene = gene_symbol, site = `Site Position`, protein = `Protein Id`) %>%
  unite("gene_site", gene_symbol, `Site Position`) %>%
  mutate(site_a = site) %>%
  unite("protein_site", `Protein Id`, site_a)

for_phosphopeptide_upset <-list(Heart = total_heart_phospho_detected$protein_site, Lung = total_lung_phospho_detected$protein_site, 
                                Kidney = total_kidney_phospho_detected$protein_site, Liver = total_liver_phospho_detected$protein_site)

pdf(paste0(figure_path, "for_phosphopeptide_upset.pdf"),
    width = 12,
    height = 8)

upset(fromList(for_phosphopeptide_upset), order.by = "degree", number.angles = 90,
      text.scale = c(2.5,2,1.5,1.5, 2.5, 2.5), mainbar.y.label = "Phospho-peptide intersection",
      sets.x.label = "Phospho-peptides per tissue",  keep.order = TRUE,
      matrix.color = "black", sets.bar.color = "black", main.bar.color = "black")
dev.off()

for_phosphoprotein_upset <-list(Heart = total_heart_phospho_detected$protein, Lung = total_lung_phospho_detected$protein, 
                                Kidney = total_kidney_phospho_detected$protein, Liver = total_liver_phospho_detected$protein)
pdf(paste0(figure_path, "for_phosphoprotein_upset.pdf"),
    width = 12,
    height = 8)


upset(fromList(for_phosphoprotein_upset), order.by = "degree", number.angles = 90,
      text.scale = c(2.5,2,1.5,1.5, 2.5, 2.5), mainbar.y.label = "Phospho-protein intersection",
      sets.x.label = "Phospho-proteins per tissue",  keep.order = TRUE,
      matrix.color = "black", sets.bar.color = "black", main.bar.color = "black")
dev.off()


heart_total_detected <- scaled_total %>%
  filter(!is.na(rowMeans(select(., heart_4.2:heart_4.30))))

lung_total_detected <- scaled_total %>%
  filter(!is.na(rowMeans(select(., lung_4.2:lung_4.30))))

kidney_total_detected <- scaled_total %>%
  filter(!is.na(rowMeans(select(., kidney_4.2:kidney_4.30))))

liver_total_detected <- scaled_total %>%
  filter(!is.na(rowMeans(select(., liver_4.2:liver_4.30))))


for_total_protein_upset <-list(Heart = heart_total_detected$`Protein Id`, Lung = lung_total_detected$`Protein Id`, 
                               Kidney = kidney_total_detected$`Protein Id`, Liver = liver_total_detected$`Protein Id`)
pdf(paste0(figure_path, "for_total_protein_upset.pdf"),
    width = 12,
    height = 8)

upset(fromList(for_total_protein_upset), order.by = "degree", number.angles = 90,
      text.scale = c(2.5,2,1.5,1.5, 2.5, 2.5), mainbar.y.label = "Total protein intersection",
      sets.x.label = "Total proteins per tissue",  keep.order = TRUE,
      matrix.color = "black", sets.bar.color = "black", main.bar.color = "black")
dev.off()


# Figure 2
# Volcano plot
lapply(all_tissues, volcano_plot)

# UpSet plots to compare significant peptides between tissues
common_peptides_all <- intersect(intersect(intersect(heart_all, lung_all), kidney_all), liver_all)

# Sig EGF-altered phosphopeptides that appear in all tissues
heart_sig_common_phospho_fold_changes <- intersect(get_tissue_sig("heart"), common_peptides_all)
lung_sig_common_phospho_fold_changes <- intersect(get_tissue_sig("lung"), common_peptides_all)
kidney_sig_common_phospho_fold_changes <- intersect(get_tissue_sig("kidney"), common_peptides_all)
liver_sig_common_phospho_fold_changes <- intersect(get_tissue_sig("liver"), common_peptides_all)

sig_peptides_common <- list(Heart = heart_sig_common_phospho_fold_changes, 
                     Lung = lung_sig_common_phospho_fold_changes, 
                     Kidney = kidney_sig_common_phospho_fold_changes,
                     Liver = liver_sig_common_phospho_fold_changes)
pdf(paste0(figure_path, "upset_egf_alone_sig_common_peptides.pdf"))
upset(fromList(sig_peptides_common), order.by = "degree", number.angles = 90,
      keep.order = TRUE, nintersects = NA, text.scale = c(2,1.5,1.5,1.5, 1.5, 2))
dev.off()


sig_peptides_all <- list(Heart = heart_sig_phospho_fold_changes, 
                     Lung = lung_sig_phospho_fold_changes, 
                     Kidney = kidney_sig_phospho_fold_changes,
                     Liver = liver_sig_phospho_fold_changes)
pdf(paste0(figure_path, "upset_egf_alone_sig_all_peptides.pdf"))
upset(fromList(sig_peptides_all), order.by = "degree", number.angles = 90,
      keep.order = TRUE, nintersects = NA, text.scale = c(2,1.5,1.5,1.5, 1.5, 2))
dev.off()


# make pairwise scatterplots
for_scatterplots <- scaled_phospho %>%  
  mutate(heart_baseline_means = rowMeans(select(., heart_4.6, heart_4.16, heart_4.26), na.rm=TRUE),
         lung_baseline_means = rowMeans(select(., lung_4.6, lung_4.16, lung_4.26), na.rm=TRUE),
         kidney_baseline_means = rowMeans(select(., kidney_4.6, kidney_4.16, kidney_4.26), na.rm=TRUE),
         liver_baseline_means = rowMeans(select(., liver_4.6, liver_4.16, liver_4.26), na.rm=TRUE)) %>%
  mutate(heart_EGF_rep1 = heart_4.2/heart_baseline_means,
         heart_EGF_rep2 = heart_4.12/heart_baseline_means,
         heart_EGF_rep3 = heart_4.22/heart_baseline_means,
         lung_EGF_rep1 = lung_4.2/lung_baseline_means,
         lung_EGF_rep2 = lung_4.12/lung_baseline_means,
         lung_EGF_rep3 = lung_4.22/lung_baseline_means,
         kidney_EGF_rep1 = kidney_4.2/kidney_baseline_means,
         kidney_EGF_rep2 = kidney_4.12/kidney_baseline_means,
         kidney_EGF_rep3 = kidney_4.22/kidney_baseline_means,
         liver_EGF_rep1 = liver_4.2/liver_baseline_means,
         liver_EGF_rep2 = liver_4.12/liver_baseline_means,
         liver_EGF_rep3 = liver_4.22/liver_baseline_means,
         
         heart_EGF_plus_MEKi_rep1 = heart_4.4/heart_baseline_means,
         heart_EGF_plus_MEKi_rep2 = heart_4.14/heart_baseline_means,
         heart_EGF_plus_MEKi_rep3 = heart_4.24/heart_baseline_means,
         lung_EGF_plus_MEKi_rep1 = lung_4.4/lung_baseline_means,
         lung_EGF_plus_MEKi_rep2 = lung_4.14/lung_baseline_means,
         lung_EGF_plus_MEKi_rep3 = lung_4.24/lung_baseline_means,
         kidney_EGF_plus_MEKi_rep1 = kidney_4.4/kidney_baseline_means,
         kidney_EGF_plus_MEKi_rep2 = kidney_4.14/kidney_baseline_means,
         kidney_EGF_plus_MEKi_rep3 = kidney_4.24/kidney_baseline_means,
         liver_EGF_plus_MEKi_rep1 = liver_4.4/liver_baseline_means,
         liver_EGF_plus_MEKi_rep2 = liver_4.14/liver_baseline_means,
         liver_EGF_plus_MEKi_rep3 = liver_4.24/liver_baseline_means,
         
         heart_MEKi_rep1 = heart_4.8/heart_baseline_means,
         heart_MEKi_rep2 = heart_4.18/heart_baseline_means,
         lung_MEKi_rep1 = lung_4.8/lung_baseline_means,
         lung_MEKi_rep2 = lung_4.18/lung_baseline_means,
         kidney_MEKi_rep1 = kidney_4.8/kidney_baseline_means,
         kidney_MEKi_rep2 = kidney_4.18/kidney_baseline_means,
         liver_MEKi_rep1 = liver_4.8/liver_baseline_means,
         liver_MEKi_rep2 = liver_4.18/liver_baseline_means) %>%
  mutate(heart_EGF_means = log(rowMeans(select(., heart_EGF_rep1:heart_EGF_rep3), na.rm = TRUE),2),
         lung_EGF_means = log(rowMeans(select(., lung_EGF_rep1:lung_EGF_rep3), na.rm = TRUE),2),
         kidney_EGF_means = log(rowMeans(select(., kidney_EGF_rep1:kidney_EGF_rep3), na.rm = TRUE),2),
         liver_EGF_means = log(rowMeans(select(., liver_EGF_rep1:liver_EGF_rep3), na.rm = TRUE),2)) %>%
  mutate(heart_EGF_plus_MEKi_means = log(rowMeans(select(., heart_EGF_plus_MEKi_rep1:heart_EGF_plus_MEKi_rep3), na.rm = TRUE),2),
         lung_EGF_plus_MEKi_means = log(rowMeans(select(., lung_EGF_plus_MEKi_rep1:lung_EGF_plus_MEKi_rep3), na.rm = TRUE),2),
         kidney_EGF_plus_MEKi_means = log(rowMeans(select(., kidney_EGF_plus_MEKi_rep1:kidney_EGF_plus_MEKi_rep3), na.rm = TRUE),2),
         liver_EGF_plus_MEKi_means = log(rowMeans(select(., liver_EGF_plus_MEKi_rep1:liver_EGF_plus_MEKi_rep3), na.rm = TRUE),2)) %>%
  mutate(gene = gene_symbol, site = `Site Position`, protein = `Protein Id`) %>%
  unite("gene_site", gene_symbol, `Site Position`) %>%
  mutate(site_a = site) %>%
  unite("protein_site", `Protein Id`, site_a) %>%
  mutate(firstsub = gsub("sp\\|", "", protein)) %>%
  mutate(firstsub = gsub("tr\\|", "", firstsub)) %>%
  mutate(uniprot = gsub("\\|.*","",firstsub)) %>%
  unite(uniprot_site, uniprot, site) %>%
  select(uniprot_site, gene_site, protein_site, "heart_EGF_means":"liver_EGF_plus_MEKi_means")

make_pairwise_scatterplot("liver", "kidney")
make_pairwise_scatterplot("liver", "lung")
make_pairwise_scatterplot("liver", "heart")
make_pairwise_scatterplot("kidney", "lung")
make_pairwise_scatterplot("lung", "heart")
make_pairwise_scatterplot("kidney", "heart")


# Kinase enrichment analysis for EGF alone
# Bubble chart for kinase enrichment

for_kinase_enrichment <- stable_veh_egf_scale_phospho  %>%
  mutate(heart_baseline_means = rowMeans(select(., heart_4.6, heart_4.16, heart_4.26), na.rm=TRUE),
         lung_baseline_means = rowMeans(select(., lung_4.6, lung_4.16, lung_4.26), na.rm=TRUE),
         kidney_baseline_means = rowMeans(select(., kidney_4.6, kidney_4.16, kidney_4.26), na.rm=TRUE),
         liver_baseline_means = rowMeans(select(., liver_4.6, liver_4.16, liver_4.26), na.rm=TRUE)) %>%
  mutate(heart_EGF_rep1 = heart_4.2/heart_baseline_means,
         heart_EGF_rep2 = heart_4.12/heart_baseline_means,
         heart_EGF_rep3 = heart_4.22/heart_baseline_means,
         lung_EGF_rep1 = lung_4.2/lung_baseline_means,
         lung_EGF_rep2 = lung_4.12/lung_baseline_means,
         lung_EGF_rep3 = lung_4.22/lung_baseline_means,
         kidney_EGF_rep1 = kidney_4.2/kidney_baseline_means,
         kidney_EGF_rep2 = kidney_4.12/kidney_baseline_means,
         kidney_EGF_rep3 = kidney_4.22/kidney_baseline_means,
         liver_EGF_rep1 = liver_4.2/liver_baseline_means,
         liver_EGF_rep2 = liver_4.12/liver_baseline_means,
         liver_EGF_rep3 = liver_4.22/liver_baseline_means,
         
         heart_EGF_plus_MEKi_rep1 = heart_4.4/heart_baseline_means,
         heart_EGF_plus_MEKi_rep2 = heart_4.14/heart_baseline_means,
         heart_EGF_plus_MEKi_rep3 = heart_4.24/heart_baseline_means,
         lung_EGF_plus_MEKi_rep1 = lung_4.4/lung_baseline_means,
         lung_EGF_plus_MEKi_rep2 = lung_4.14/lung_baseline_means,
         lung_EGF_plus_MEKi_rep3 = lung_4.24/lung_baseline_means,
         kidney_EGF_plus_MEKi_rep1 = kidney_4.4/kidney_baseline_means,
         kidney_EGF_plus_MEKi_rep2 = kidney_4.14/kidney_baseline_means,
         kidney_EGF_plus_MEKi_rep3 = kidney_4.24/kidney_baseline_means,
         liver_EGF_plus_MEKi_rep1 = liver_4.4/liver_baseline_means,
         liver_EGF_plus_MEKi_rep2 = liver_4.14/liver_baseline_means,
         liver_EGF_plus_MEKi_rep3 = liver_4.24/liver_baseline_means,
         
         heart_MEKi_rep1 = heart_4.8/heart_baseline_means,
         heart_MEKi_rep2 = heart_4.18/heart_baseline_means,
         lung_MEKi_rep1 = lung_4.8/lung_baseline_means,
         lung_MEKi_rep2 = lung_4.18/lung_baseline_means,
         kidney_MEKi_rep1 = kidney_4.8/kidney_baseline_means,
         kidney_MEKi_rep2 = kidney_4.18/kidney_baseline_means,
         liver_MEKi_rep1 = liver_4.8/liver_baseline_means,
         liver_MEKi_rep2 = liver_4.18/liver_baseline_means) %>%
  
  mutate(heart_EGF_means = rowMeans(select(., heart_EGF_rep1:heart_EGF_rep3), na.rm = TRUE),
         lung_EGF_means = rowMeans(select(., lung_EGF_rep1:lung_EGF_rep3), na.rm = TRUE),
         kidney_EGF_means = rowMeans(select(., kidney_EGF_rep1:kidney_EGF_rep3), na.rm = TRUE),
         liver_EGF_means = rowMeans(select(., liver_EGF_rep1:liver_EGF_rep3), na.rm = TRUE)) %>%
  
  mutate(heart_EGF_plus_minus_MEKi_fold_change = heart_EGF_means/
           rowMeans(select(., heart_EGF_plus_MEKi_rep1:heart_EGF_plus_MEKi_rep3)),
         lung_EGF_plus_minus_MEKi_fold_change = lung_EGF_means/
           rowMeans(select(., lung_EGF_plus_MEKi_rep1:lung_EGF_plus_MEKi_rep3)),
         kidney_EGF_plus_minus_MEKi_fold_change = kidney_EGF_means/
           rowMeans(select(., kidney_EGF_plus_MEKi_rep1:kidney_EGF_plus_MEKi_rep3)),
         liver_EGF_plus_minus_MEKi_fold_change = liver_EGF_means/
           rowMeans(select(., liver_EGF_plus_MEKi_rep1:liver_EGF_plus_MEKi_rep3))) %>%
  
  mutate(heart_EGF_plus_MEKi_means = rowMeans(select(., heart_EGF_plus_MEKi_rep1:heart_EGF_plus_MEKi_rep3)),
         lung_EGF_plus_MEKi_means = rowMeans(select(., lung_EGF_plus_MEKi_rep1:lung_EGF_plus_MEKi_rep3)),
         kidney_EGF_plus_MEKi_means = rowMeans(select(., kidney_EGF_plus_MEKi_rep1:kidney_EGF_plus_MEKi_rep3)),
         liver_EGF_plus_MEKi_means =  rowMeans(select(., liver_EGF_plus_MEKi_rep1:liver_EGF_plus_MEKi_rep3))) %>%
  unite(protein_site, `Protein Id`, `Site Position`)

lapply(all_tissues, for_kinase_enrichment_analysis, for_kinase_enrichment, "egf_alone")

lung_egf_kinase_enrichment <- read_tsv(paste0(path, "final_kinase_enrichment/lung EGF alone enrichment.txt")) %>%
  filter(dominant_adjusted_p_value < 0.05) %>%
  rename(lung_enrichment_value = dominant_enrichment_value_log2) %>%
  select(kinase, lung_enrichment_value, dominant_direction) 
lung_egf_kinase_enrichment$lung_color <- 0 
lung_egf_kinase_enrichment$lung_color[lung_egf_kinase_enrichment$dominant_direction == "downregulated set"] <- 1
lung_egf_kinase_enrichment <- lung_egf_kinase_enrichment %>% select(-dominant_direction)
kidney_egf_kinase_enrichment <- read_tsv(paste0(path, "final_kinase_enrichment/kidney EGF alone enrichment.txt")) %>%
  filter(dominant_adjusted_p_value < 0.05) %>%
  rename(kidney_enrichment_value = dominant_enrichment_value_log2) %>%
  select(kinase, kidney_enrichment_value, dominant_direction) 
kidney_egf_kinase_enrichment$kidney_color <- 0 
kidney_egf_kinase_enrichment$kidney_color[kidney_egf_kinase_enrichment$dominant_direction == "downregulated set"] <- 1
kidney_egf_kinase_enrichment <- kidney_egf_kinase_enrichment %>% select(-dominant_direction)
liver_egf_kinase_enrichment <- read_tsv(paste0(path, "final_kinase_enrichment/liver EGF alone enrichment.txt")) %>%
  filter(dominant_adjusted_p_value < 0.05) %>%
  rename(liver_enrichment_value = dominant_enrichment_value_log2) %>%
  select(kinase, liver_enrichment_value, dominant_direction)
liver_egf_kinase_enrichment$liver_color <- 0 
liver_egf_kinase_enrichment$liver_color[liver_egf_kinase_enrichment$dominant_direction == "downregulated set"] <- 1
liver_egf_kinase_enrichment <- liver_egf_kinase_enrichment %>% select(-dominant_direction)

#merge all data frames in list
joined_df_kinase_enrichment <- full_join(lung_egf_kinase_enrichment, kidney_egf_kinase_enrichment) %>%
  full_join(liver_egf_kinase_enrichment) %>%
  replace(is.na(.),0)

gathered_joined_df_kinase_enrichment <- joined_df_kinase_enrichment %>%
  gather(Tissue, Value, -c(kinase, lung_color, kidney_color, liver_color)) %>%
  rename(lung_enrichment_value = lung_color,
         kidney_enrichment_value = kidney_color,
         liver_enrichment_value = liver_color) %>%
  gather(color, tissue_color, -c(Tissue, Value, kinase)) %>%
  filter(Tissue == color) %>%
  filter(!(Value == 0)) %>%
  mutate(
    Tissue = as.factor(Tissue),
    kinase= as.factor(kinase),
    color = as.factor(tissue_color)) %>%
  mutate(Value = abs(Value)) 

pdf(paste0(figure_path, "kinase_enrichment_EGF_alone_dotplot.pdf"), width = 15, height = 6)
ggplot(gathered_joined_df_kinase_enrichment, aes(kinase, Tissue)) +
  geom_point(aes(size = Value, color=as.factor(tissue_color))) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_colour_manual(values = c("black", "lightblue")) 
dev.off()




lung_MF <- do_GO_analysis_on_EGF("lung", "MF")
kidney_MF <- do_GO_analysis_on_EGF("kidney", "MF")
liver_MF <- do_GO_analysis_on_EGF("liver", "MF")
heart_MF <- do_GO_analysis_on_EGF("heart", "MF")

pdf(paste0(figure_path, "lung_MF_EGF_sig_GO_plot.pdf"))
dotplot(lung_MF, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "kidney_MF_EGF_sig_GO_plot.pdf"))
dotplot(kidney_MF, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "liver_MF_EGF_sig_GO_plot.pdf"))
dotplot(liver_MF, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "heart_MF_EGF_sig_GO_plot.pdf"))
dotplot(heart_MF, showCategory=10, font.size = 6)
dev.off()


lung_BP <- do_GO_analysis_on_EGF("lung", "BP")
kidney_BP <- do_GO_analysis_on_EGF("kidney", "BP")
liver_BP <- do_GO_analysis_on_EGF("liver", "BP")
heart_BP <- do_GO_analysis_on_EGF("heart", "BP")


pdf(paste0(figure_path, "lung_BP_EGF_sig_GO_plot.pdf"))
dotplot(lung_BP, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "kidney_BP_EGF_sig_GO_plot.pdf"))
dotplot(kidney_BP, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "liver_BP_EGF_sig_GO_plot.pdf"))
dotplot(liver_BP, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "heart_BP_EGF_sig_GO_plot.pdf"))
dotplot(heart_BP, showCategory=10, font.size = 6)
dev.off()

lung_CC <- do_GO_analysis_on_EGF("lung", "CC")
kidney_CC <- do_GO_analysis_on_EGF("kidney", "CC")
liver_CC <- do_GO_analysis_on_EGF("liver", "CC")
heart_CC <- do_GO_analysis_on_EGF("heart", "CC")

pdf(paste0(figure_path, "lung_CC_EGF_sig_GO_plot.pdf"))
dotplot(lung_CC, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "kidney_CC_EGF_sig_GO_plot.pdf"))
dotplot(kidney_CC, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "liver_CC_EGF_sig_GO_plot.pdf"))
dotplot(liver_CC, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "heart_CC_EGF_sig_GO_plot.pdf"))
dotplot(heart_CC, showCategory=10, font.size = 6)
dev.off()


# Figure 3

make_EGF_versus_plusMEKi_scatterplot("heart")
make_EGF_versus_plusMEKi_scatterplot("lung")
make_EGF_versus_plusMEKi_scatterplot("kidney")
make_EGF_versus_plusMEKi_scatterplot("liver")

# make UpSet plots for ERK-dependent sites

# intersect sig EGF phosphopeptides found in all tissues with sig EGF/EGF + MEKi fold-changes
heart_EGF_sig_common_ERK_dependent <- intersect(heart_sig_common_phospho_fold_changes,
                                                heart_sig_ERK_dependent_phospho_fold_changes)

lung_EGF_sig_common_ERK_dependent <- intersect(lung_sig_common_phospho_fold_changes,
                                               lung_sig_ERK_dependent_phospho_fold_changes)

kidney_EGF_sig_common_ERK_dependent <- intersect(kidney_sig_common_phospho_fold_changes,
                                                 kidney_sig_ERK_dependent_phospho_fold_changes)

liver_EGF_sig_common_ERK_dependent <- intersect(liver_sig_common_phospho_fold_changes,
                                                liver_sig_ERK_dependent_phospho_fold_changes)


sig_ERK_dependent_peptides <- list(Heart = heart_EGF_sig_common_ERK_dependent, 
                                   Lung = lung_EGF_sig_common_ERK_dependent, 
                                   Kidney = kidney_EGF_sig_common_ERK_dependent,
                                   Liver = liver_EGF_sig_common_ERK_dependent)

pdf(paste0(figure_path, "upset_erk_dependent_egf_alone_sig_common_peptides.pdf"))
upset(fromList(sig_ERK_dependent_peptides), order.by = "degree",
      keep.order = TRUE, nintersects = NA, text.scale = c(2,1.5,1.5,1.5, 1.5, 2))
dev.off()


# make kidney and lung pie charts for Erk dependency

lung_kidney_sig_phospho_fold_changes <- intersect(lung_sig_phospho_fold_changes,
                                                  kidney_sig_phospho_fold_changes)
lung_erk_dependent_lung_kidney_common_sig_phospho_fold_changes <- intersect(lung_kidney_sig_phospho_fold_changes,
                                                                            lung_sig_ERK_dependent_phospho_fold_changes)
kidney_erk_dependent_lung_kidney_common_sig_phospho_fold_changes <- intersect(lung_kidney_sig_phospho_fold_changes,
                                                                            kidney_sig_ERK_dependent_phospho_fold_changes)

lung_alone_sig_phospho_fold_changes <- setdiff(lung_sig_common_phospho_fold_changes,
                                               c(kidney_sig_common_phospho_fold_changes,
                                                 liver_sig_common_phospho_fold_changes,
                                                 heart_sig_common_phospho_fold_changes))
kidney_alone_sig_phospho_fold_changes <- setdiff(kidney_sig_common_phospho_fold_changes,
                                               c(lung_sig_common_phospho_fold_changes,
                                                 liver_sig_common_phospho_fold_changes,
                                                 heart_sig_common_phospho_fold_changes))

lung_erk_dependent_unique_sig_phospho_fold_changes <- intersect(lung_sig_ERK_dependent_phospho_fold_changes,
                                                                lung_alone_sig_phospho_fold_changes)
kidney_erk_dependent_unique_sig_phospho_fold_changes <- intersect(kidney_sig_ERK_dependent_phospho_fold_changes,
                                                                kidney_alone_sig_phospho_fold_changes)

labels = c("Erk-dependent", "Erk-independent")
kidney_erk_dependent_slices <- c(length(kidney_erk_dependent_lung_kidney_common_sig_phospho_fold_changes),
                                 length(setdiff(lung_kidney_sig_phospho_fold_changes, kidney_erk_dependent_lung_kidney_common_sig_phospho_fold_changes)))
kidney_erk_dependent_pct <- round(kidney_erk_dependent_slices/sum(kidney_erk_dependent_slices)*100)
lung_erk_dependent_slices <- c(length(lung_erk_dependent_lung_kidney_common_sig_phospho_fold_changes),
                                 length(setdiff(lung_kidney_sig_phospho_fold_changes, lung_erk_dependent_lung_kidney_common_sig_phospho_fold_changes)))
lung_erk_dependent_pct <- round(lung_erk_dependent_slices/sum(lung_erk_dependent_slices)*100)

pdf(paste0(figure_path, "kidney_erk_dependent_common_piechart.pdf"))
pie(kidney_erk_dependent_slices, labels = paste0(labels, " ", kidney_erk_dependent_pct, "%"))
dev.off()
pdf(paste0(figure_path, "lung_erk_dependent_common_piechart.pdf"))
pie(lung_erk_dependent_slices, labels = paste0(labels, " ", lung_erk_dependent_pct, "%"))
dev.off()

unique_kidney_erk_dependent_slices <- c(length(kidney_erk_dependent_unique_sig_phospho_fold_changes),
                                        length(setdiff(kidney_alone_sig_phospho_fold_changes, kidney_erk_dependent_unique_sig_phospho_fold_changes)))
kidney_unique_erk_dependent_pct <- round(unique_kidney_erk_dependent_slices/sum(unique_kidney_erk_dependent_slices)*100)

unique_lung_erk_dependent_slices <- c(length(lung_erk_dependent_unique_sig_phospho_fold_changes),
                               length(setdiff(lung_alone_sig_phospho_fold_changes, lung_erk_dependent_unique_sig_phospho_fold_changes)))
lung_unique_erk_dependent_pct <- round(unique_lung_erk_dependent_slices/sum(unique_lung_erk_dependent_slices)*100)

pdf(paste0(figure_path, "kidney_unique_erk_dependent_common_piechart.pdf"))
pie(unique_kidney_erk_dependent_slices, labels = paste0(labels, " ", kidney_unique_erk_dependent_pct, "%"))
dev.off()
pdf(paste0(figure_path, "lung_unique_erk_dependent_common_piechart.pdf"))
pie(unique_lung_erk_dependent_slices, labels = paste0(labels, " ", lung_unique_erk_dependent_pct, "%"))  
dev.off()

# make heat map

for_pre_heat_map <- scaled_phospho %>%  
  mutate(heart_baseline_means = rowMeans(select(., heart_4.6, heart_4.16, heart_4.26), na.rm=TRUE),
         lung_baseline_means = rowMeans(select(., lung_4.6, lung_4.16, lung_4.26), na.rm=TRUE),
         kidney_baseline_means = rowMeans(select(., kidney_4.6, kidney_4.16, kidney_4.26), na.rm=TRUE),
         liver_baseline_means = rowMeans(select(., liver_4.6, liver_4.16, liver_4.26), na.rm=TRUE)) %>%
  mutate(heart_EGF_rep1 = heart_4.2/heart_baseline_means,
         heart_EGF_rep2 = heart_4.12/heart_baseline_means,
         heart_EGF_rep3 = heart_4.22/heart_baseline_means,
         lung_EGF_rep1 = lung_4.2/lung_baseline_means,
         lung_EGF_rep2 = lung_4.12/lung_baseline_means,
         lung_EGF_rep3 = lung_4.22/lung_baseline_means,
         kidney_EGF_rep1 = kidney_4.2/kidney_baseline_means,
         kidney_EGF_rep2 = kidney_4.12/kidney_baseline_means,
         kidney_EGF_rep3 = kidney_4.22/kidney_baseline_means,
         liver_EGF_rep1 = liver_4.2/liver_baseline_means,
         liver_EGF_rep2 = liver_4.12/liver_baseline_means,
         liver_EGF_rep3 = liver_4.22/liver_baseline_means,
         
         heart_EGF_plus_MEKi_rep1 = heart_4.4/heart_baseline_means,
         heart_EGF_plus_MEKi_rep2 = heart_4.14/heart_baseline_means,
         heart_EGF_plus_MEKi_rep3 = heart_4.24/heart_baseline_means,
         lung_EGF_plus_MEKi_rep1 = lung_4.4/lung_baseline_means,
         lung_EGF_plus_MEKi_rep2 = lung_4.14/lung_baseline_means,
         lung_EGF_plus_MEKi_rep3 = lung_4.24/lung_baseline_means,
         kidney_EGF_plus_MEKi_rep1 = kidney_4.4/kidney_baseline_means,
         kidney_EGF_plus_MEKi_rep2 = kidney_4.14/kidney_baseline_means,
         kidney_EGF_plus_MEKi_rep3 = kidney_4.24/kidney_baseline_means,
         liver_EGF_plus_MEKi_rep1 = liver_4.4/liver_baseline_means,
         liver_EGF_plus_MEKi_rep2 = liver_4.14/liver_baseline_means,
         liver_EGF_plus_MEKi_rep3 = liver_4.24/liver_baseline_means,
         
         heart_MEKi_rep1 = heart_4.8/heart_baseline_means,
         heart_MEKi_rep2 = heart_4.18/heart_baseline_means,
         lung_MEKi_rep1 = lung_4.8/lung_baseline_means,
         lung_MEKi_rep2 = lung_4.18/lung_baseline_means,
         kidney_MEKi_rep1 = kidney_4.8/kidney_baseline_means,
         kidney_MEKi_rep2 = kidney_4.18/kidney_baseline_means,
         liver_MEKi_rep1 = liver_4.8/liver_baseline_means,
         liver_MEKi_rep2 = liver_4.18/liver_baseline_means) %>%
  filter(gene_symbol == "Mapk3" | gene_symbol == "Mapk1") %>%
  #filter(gene_symbol == "Tsc2" | gene_symbol == "Tsc1") %>%
  #filter(gene_symbol == "Cdkn1b" | gene_symbol == "Cdk2") %>%
  #filter(gene_symbol == "Egfr") %>%
  unite("gene_site", gene_symbol, "Site Position") 

all_conditions_heat_map <- for_pre_heat_map %>%
  select(c(gene_site, heart_EGF_rep1:heart_EGF_rep3,
           heart_EGF_plus_MEKi_rep1:heart_EGF_plus_MEKi_rep3,
           heart_MEKi_rep1:heart_MEKi_rep2, 
           lung_EGF_rep1:lung_EGF_rep3,
           lung_EGF_plus_MEKi_rep1:lung_EGF_plus_MEKi_rep3,
           lung_MEKi_rep1:lung_MEKi_rep2, 
           kidney_EGF_rep1:kidney_EGF_rep3,
           kidney_EGF_plus_MEKi_rep1:kidney_EGF_plus_MEKi_rep3, 
           kidney_MEKi_rep1:kidney_MEKi_rep2,
           liver_EGF_rep1:liver_EGF_rep3,
           liver_EGF_plus_MEKi_rep1:liver_EGF_plus_MEKi_rep3,
           liver_MEKi_rep1:liver_MEKi_rep2))


for_heat_map <- all_conditions_heat_map

scaled_for_heat_map <- data.frame(log2(for_heat_map[,2:ncol(for_heat_map)]))
rownames(scaled_for_heat_map) <- make.names(for_heat_map$gene_site, unique = TRUE)

draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 0.5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
paletteLength <- 26
myColor <- colorRampPalette(c("darkblue","white", "darkred"))(paletteLength)

bk1 <- c(seq(-2.6, -0.2, by =0.2), -0.10001)
bk2 <- c(0.10001, seq(.2,2.6, by = 0.2))
bk <- c(bk1, bk2)

pdf(paste0(figure_path, "all_conditions_ERK_heat_map.pdf"))
pheatmap(scaled_for_heat_map, cluster_cols = F, cluster_rows=F, 
         color=myColor, breaks=bk, na_col = "grey", border_color = "black",
         gaps_col = c(8, 16, 24))
dev.off()



# EGF plus MEKi kinase enrichment

lapply(all_tissues, for_kinase_enrichment_analysis, for_kinase_enrichment, "egf_plus_minus_MEKi")

heart_egf_plus_MEKi_kinase_enrichment <- read_tsv(paste0(path, "final_kinase_enrichment/heart EGF plus minus MEKi enrichment.txt")) %>%
  filter(dominant_adjusted_p_value < 0.05) %>%
  rename(heart_enrichment_value = dominant_enrichment_value_log2) %>%
  select(kinase, heart_enrichment_value, dominant_direction) 
heart_egf_plus_MEKi_kinase_enrichment$heart_color <- 0 
heart_egf_plus_MEKi_kinase_enrichment$heart_color[heart_egf_plus_MEKi_kinase_enrichment$dominant_direction == "downregulated set"] <- 1
heart_egf_plus_MEKi_kinase_enrichment <- heart_egf_plus_MEKi_kinase_enrichment %>% select(-dominant_direction)
lung_egf_plus_MEKi_kinase_enrichment <- read_tsv(paste0(path, "final_kinase_enrichment/lung EGF plus minus MEKi enrichment.txt")) %>%
  filter(dominant_adjusted_p_value < 0.05) %>%
  rename(lung_enrichment_value = dominant_enrichment_value_log2) %>%
  select(kinase, lung_enrichment_value, dominant_direction) 
lung_egf_plus_MEKi_kinase_enrichment$lung_color <- 0 
lung_egf_plus_MEKi_kinase_enrichment$lung_color[lung_egf_plus_MEKi_kinase_enrichment$dominant_direction == "downregulated set"] <- 1
lung_egf_plus_MEKi_kinase_enrichment <- lung_egf_plus_MEKi_kinase_enrichment %>% select(-dominant_direction)
kidney_egf_plus_MEKi_kinase_enrichment <- read_tsv(paste0(path, "final_kinase_enrichment/kidney EGF plus minus MEKi enrichment.txt")) %>%
  filter(dominant_adjusted_p_value < 0.05) %>%
  rename(kidney_enrichment_value = dominant_enrichment_value_log2) %>%
  select(kinase, kidney_enrichment_value, dominant_direction) 
kidney_egf_plus_MEKi_kinase_enrichment$kidney_color <- 0 
kidney_egf_plus_MEKi_kinase_enrichment$kidney_color[kidney_egf_plus_MEKi_kinase_enrichment$dominant_direction == "downregulated set"] <- 1
kidney_egf_plus_MEKi_kinase_enrichment <- kidney_egf_plus_MEKi_kinase_enrichment %>% select(-dominant_direction)
liver_egf_plus_MEKi_kinase_enrichment <- read_tsv(paste0(path, "final_kinase_enrichment/liver EGF plus minus MEKi enrichment.txt")) %>%
  filter(dominant_adjusted_p_value < 0.05) %>%
  rename(liver_enrichment_value = dominant_enrichment_value_log2) %>%
  select(kinase, liver_enrichment_value, dominant_direction)
liver_egf_plus_MEKi_kinase_enrichment$liver_color <- 0 
liver_egf_plus_MEKi_kinase_enrichment$liver_color[liver_egf_plus_MEKi_kinase_enrichment$dominant_direction == "downregulated set"] <- 1
liver_egf_plus_MEKi_kinase_enrichment <- liver_egf_plus_MEKi_kinase_enrichment %>% select(-dominant_direction)


joined_df_kinase_enrichment <- full_join(lung_egf_plus_MEKi_kinase_enrichment, kidney_egf_plus_MEKi_kinase_enrichment) %>%
  full_join(liver_egf_plus_MEKi_kinase_enrichment) %>%
  full_join(heart_egf_plus_MEKi_kinase_enrichment) %>%
  replace(is.na(.),0)

gathered_joined_df_kinase_enrichment <- joined_df_kinase_enrichment %>%
  gather(Tissue, Value, -c(kinase, heart_color, lung_color, kidney_color, liver_color)) %>%
  rename(heart_enrichment_value = heart_color, 
         lung_enrichment_value = lung_color,
         kidney_enrichment_value = kidney_color,
         liver_enrichment_value = liver_color) %>%
  gather(color, tissue_color, -c(Tissue, Value, kinase)) %>%
  filter(Tissue == color) %>%
  filter(!(Value == 0)) %>%
  mutate(
    Tissue = as.factor(Tissue),
    kinase= as.factor(kinase),
    color = as.factor(tissue_color)) %>%
  mutate(Value = abs(Value)) 

pdf(paste0(figure_path, "kinase_enrichment_EGF_plus_MEKi_dotplot.pdf"), width = 15, height = 6)
ggplot(gathered_joined_df_kinase_enrichment, aes(kinase, Tissue)) +
  geom_point(aes(size = Value, color=as.factor(tissue_color))) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_colour_manual(values = c("black", "lightblue")) 
dev.off()


# GO analysis on EGF + MEKi
#lung_egf_meki_not_egfmekisig <- do_GO_analysis_on_EGF_MEKi_not_egfmekisig("lung", lung_sig_phospho_proteins)
#kidney_egf_meki <- do_GO_analysis_on_EGF_MEKi("kidney", kidney_sig_phospho_proteins)
#liver_egf_meki_not_egfmekisig <- do_GO_analysis_on_EGF_MEKi_not_egfmekisig("liver", liver_sig_ERK_dependent_phospho_fold_changes)
#heart_egf_meki <- do_GO_analysis_on_EGF_MEKi("heart", heart_sig_phospho_proteins)

lung_egf_meki <- do_GO_analysis_on_EGF_MEKi("lung", lung_sig_phospho_proteins)
kidney_egf_meki <- do_GO_analysis_on_EGF_MEKi("kidney", kidney_sig_phospho_proteins)
liver_egf_meki <- do_GO_analysis_on_EGF_MEKi("liver", liver_sig_phospho_proteins)
heart_egf_meki <- do_GO_analysis_on_EGF_MEKi("heart", heart_sig_phospho_proteins)

pdf(paste0(figure_path, "lung_EGFMEKi_onlysig_sig_GO_plot.pdf"))
dotplot(lung_egf_meki, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "kidney_EGFMEKi_onlysig_sig_GO_plot.pdf"))
dotplot(kidney_egf_meki, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "liver_EGFMEKi_onlysig_sig_GO_plot.pdf"))
dotplot(liver_egf_meki, showCategory=10, font.size = 6)
dev.off()
pdf(paste0(figure_path, "heart_EGFMEKi_onlysig_sig_GO_plot.pdf"))
dotplot(heart_egf_meki, showCategory=10, font.size = 6)
dev.off()






# compare MSstats to Core

heart_MSstats <- get_tissue_all_MSstats("heart")
lung_MSstats <- get_tissue_all_MSstats("lung")
kidney_MSstats <- get_tissue_all_MSstats("kidney")
liver_MSstats <- get_tissue_all_MSstats("liver")

heart_Core <- read_csv(paste0(path, "heart_EGF_alone_fold_changes.csv")) 
lung_Core <- read_csv(paste0(path, "lung_EGF_alone_fold_changes.csv"))
kidney_Core <- read_csv(paste0(path, "kidney_EGF_alone_fold_changes.csv"))
liver_Core <- read_csv(paste0(path, "liver_EGF_alone_fold_changes.csv"))


heart_combined <- left_join(heart_MSstats, heart_Core) %>%
  select(log2FC, heart_EGF_means) %>%
  na.omit()
lung_combined <- left_join(lung_MSstats, lung_Core) %>%
  select(log2FC, lung_EGF_means) %>%
  na.omit()
kidney_combined <- left_join(kidney_MSstats, kidney_Core) %>%
  select(protein_site, log2FC, kidney_EGF_means) %>%
  na.omit()
liver_combined <- left_join(liver_MSstats, liver_Core) %>%
  select(protein_site, log2FC, liver_EGF_means) %>%
  na.omit()


cor(heart_combined$log2FC, heart_combined$heart_EGF_means)
cor(lung_combined$log2FC, lung_combined$lung_EGF_means)
cor(kidney_combined$log2FC, kidney_combined$kidney_EGF_means)
cor(liver_combined$log2FC, liver_combined$liver_EGF_means)

pdf(paste0(figure_path, "heart_core_MSstats.pdf"))
ggplot(heart_combined, aes(x=log2FC, y=heart_EGF_means)) + 
  geom_point(color = "darkred") + theme_bw() + 
  annotate("text", x = 0, y = 0, 
           label = round(cor(heart_combined$log2FC, heart_combined$heart_EGF_means), 2)) +
  geom_abline(slope = 1, intercept = 0)
dev.off()

pdf(paste0(figure_path, "lung_core_MSstats.pdf"))
ggplot(lung_combined, aes(x=log2FC, y=lung_EGF_means)) + 
  geom_point(color = "blue") + theme_bw() + 
  annotate("text", x = 0, y = 0, 
           label = round(cor(lung_combined$log2FC, lung_combined$lung_EGF_means), 2)) +
  geom_abline(slope = 1, intercept = 0)
dev.off()


pdf(paste0(figure_path, "kidney_core_MSstats.pdf"))
ggplot(kidney_combined, aes(x=log2FC, y=kidney_EGF_means)) + 
  geom_point(color = "purple") + theme_bw() + 
  annotate("text", x = 0, y = 0, 
           label = round(cor(kidney_combined$log2FC, kidney_combined$kidney_EGF_means), 2)) +
  geom_abline(slope = 1, intercept = 0)
dev.off()

pdf(paste0(figure_path, "liver_core_MSstats.pdf"))
ggplot(liver_combined, aes(x=log2FC, y=liver_EGF_means)) + 
  geom_point(color = "orange") + theme_bw() + 
  annotate("text", x = 0, y = 0, 
           label = round(cor(liver_combined$log2FC, liver_combined$liver_EGF_means), 2)) +
  geom_abline(slope = 1, intercept = 0)
dev.off()



