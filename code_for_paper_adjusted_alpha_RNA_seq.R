# Set path for data and for figures
path = "/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/"
figure_path = "/Users/beatrice/Desktop/2023_03_12_Haigis_lab/Coding/EGF_and_MEKi_treatment_final/R_figures_for_paper/"

library(ggplot2)
library("DESeq2")
library(AnnotationDbi)
library(org.Mm.eg.db)
library(tidyverse)
library(clusterProfiler)
library("RColorBrewer")
library(pheatmap)
library("grid")
library("msigdbr")
library(UpSetR)

# read in data
gene_count <- read_delim("/Users/beatrice/Desktop/2023_03_12_Haigis_lab/RNA_seq/R_analysis/gene_count.txt")

# Functions 

get_dds_object <- function(tissue){
  countData <- gene_count %>%
    dplyr::select(gene_id, paste0(tissue, "8") : paste0(tissue, "24")) %>% data.frame()
  
  groups <- c(rep("EGF_2h", 4), rep("veh_2h", 3), rep("EGF_6h", 4), rep("veh_6h", 3), rep("untr", 3))
  metaData <- data.frame(cbind(colnames(countData[-1]), groups))
  colnames(metaData) <- c("sample", "group")
  
  dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metaData, 
                                design=~group, tidy = TRUE)
  dds <- estimateSizeFactors(dds)
  #keep <- rowSums(counts(dds)) >= 10
  keep <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= 3
  print(sum(keep))
  dds <- dds[keep,]
  print(head(counts(dds)))
  dds <- DESeq(dds)
  return(dds)
}

get_dds_results <- function(dds_object){
  res_2h <- results(dds_object, contrast = c("group", "EGF_2h", "veh_2h"), alpha = 0.05)
  res_6h <- results(dds_object, contrast = c("group", "EGF_6h", "veh_6h"), alpha = 0.05)
  return(list(res_2h, res_6h))
}

volcano_plot <- function(tissue_res, tissue_time){
  df <- data.frame(tissue_res) # DDS result df
  df$color <- 0
  df$color[df$padj < 0.05 & abs(df$log2FoldChange)<1] <- 1
  df$color[df$padj < 0.05 & abs(df$log2FoldChange)>1] <- 2
  df$gene_id <- rownames(df)
  #gene_ids <- rownames(df[df$color == 2,])
  #gene_names <- gene_count %>% filter(gene_id %in% gene_ids) %>% select(gene_id, gene_name)
  #joined_df <- left_join(df, gene_names, by = "gene_id")
  #print(head(joined_df))
  #View(joined_df)
  #df$color[df$significance == "n.s."] <- 1
  print(head(df))
  g= ggplot(df, aes(x=log2FoldChange, y=-log10(padj), group = color)) + 
    geom_point(aes(color=as.factor(color)), show.legend = FALSE)  +  
    labs(x = "log2(fc)", y = "-log10(pvalue)") +   
    geom_vline(xintercept = 1, linetype = "dashed") + 
    geom_vline(xintercept = -1, linetype = "dashed") + 
    geom_hline(yintercept = 1.3, linetype = "dashed") + 
    annotate("text", x = 8, y = 4, 
             label = paste0("p.adj < 0.05 & \n log2(FC) > 2: n = ", sum(df$color == 2))) + 
    annotate("text", x = -8, y = 4, 
             label = paste0("p.adj < 0.05 & \n log2(FC) < 2: n = ", sum(df$color == 1))) + 
    #geom_text(aes(label=ifelse(color==2,as.character(gene_name),'')),hjust=0,vjust=0) +
    scale_colour_manual(values = c("black","blue", "red")) + xlim(c(-10, 10)) + theme_classic()
  pdf(paste0(figure_path, tissue_time, "_gene_expression_volcano_plot.pdf"))
  plot(g)
  dev.off()
}

get_sig_res <- function(tissue_res){
  df <- data.frame(tissue_res)
  sig_res <- df %>% filter(padj < 0.05 & abs(log2FoldChange)>1)
  print(head(sig_res))
  sig_res$ensembl <- rownames(sig_res)
  gene_symbols <- data.frame(mapIds(org.Mm.eg.db, keys = rownames(sig_res), keytype = "ENSEMBL", column = "SYMBOL"))
  print(head(gene_symbols))
  gene_symbols$ensembl <- rownames(gene_symbols)
  combine <- inner_join(select(sig_res, c(ensembl, log2FoldChange)), gene_symbols, by = "ensembl")
  #combine <- merge(sig_res["log2FoldChange"], gene_symbols, by = 'row.names', all = TRUE) 
  print(head(combine))
  colnames(combine) <- c("ENSEMBL_ID", "log2FoldChange", "gene_symbol")
  combine <- combine %>% arrange(log2FoldChange)
  return(combine)
}


up_do_GO_analysis_on_sig_genes_mf <- function(sig_gene_list, tissue_dds){
  up_genes <- sig_gene_list[which(sig_gene_list$log2FoldChange>1),]
  gene <- up_genes$ENSEMBL_ID
  universe <- rownames(tissue_dds) # all genes in DDS object serve as univeerse
  
  ego <- enrichGO(gene          = gene,
                  universe = universe, 
                  OrgDb         = org.Mm.eg.db,
                  keyType = "ENSEMBL",
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  #barplot(ggo, showCategory=20) 
  return(ego)
  
  
}

make_2h_vs_6h_scatterplots <- function(tissue_res, sig_2h, sig_6h, tissue){
  res2h <- data.frame(tissue_res[[1]]$log2FoldChange)
  res2h$ensembl <- rownames(tissue_res[[1]])
  colnames(res2h) <- c("log2fc_2h", "ensembl")
  res6h <- data.frame(tissue_res[[2]]$log2FoldChange)
  res6h$ensembl <- rownames(tissue_res[[2]])
  colnames(res6h) <- c("log2fc_6h", "ensembl")
  combined_df <- inner_join(res2h, res6h) %>%
    filter(ensembl %in% c(sig_2h$ENSEMBL_ID, sig_6h$ENSEMBL_ID))
  print(intersect(sig_2h$ENSEMBL_ID, sig_6h$ENSEMBL_ID))
  combined_df$color <- 0
  combined_df$color[combined_df$ensembl %in% sig_2h$ENSEMBL_ID] <-  1
  combined_df$color[combined_df$ensembl %in% sig_6h$ENSEMBL_ID] <-  2
  combined_df$color[combined_df$ensembl %in% intersect(sig_2h$ENSEMBL_ID, sig_6h$ENSEMBL_ID)] <-  3
  g =  ggplot(combined_df, aes(x=log2fc_2h, y=log2fc_6h, group = color)) + 
    geom_point(aes(color=as.factor(color)), show.legend = FALSE)  +  geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + scale_colour_manual(values = c("coral2", "deepskyblue4", "chartreuse")) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = 2) +
    scale_x_continuous(limits = c(-9, 9)) +
    scale_y_continuous(limits = c(-9, 9)) 
  pdf(paste0(figure_path, tissue, "_2h_versus_6h_DEGs.pdf"))
  plot(g)
  dev.off()
}


# make DEseq2 object on whole data set

countData <- gene_count %>%
  dplyr::select(gene_id, He8 : Liv24) %>% data.frame()
groups <- c(rep("He_EGF_2h", 4), rep("He_veh_2h", 3), rep("He_EGF_6h", 4), rep("He_veh_6h", 3), rep("He_untr", 3),
            rep("Lu_EGF_2h", 4), rep("Lu_veh_2h", 3), rep("Lu_EGF_6h", 4), rep("Lu_veh_6h", 3), rep("Lu_untr", 3),
            rep("Kid_EGF_2h", 4), rep("Kid_veh_2h", 3), rep("Kid_EGF_6h", 4), rep("Kid_veh_6h", 3), rep("Kid_untr", 3),
            rep("Liv_EGF_2h", 4), rep("Liv_veh_2h", 3), rep("Liv_EGF_6h", 4), rep("Liv_veh_6h", 3), rep("Liv_untr", 3))

metaData <- data.frame(cbind(colnames(countData[-1]), groups))
colnames(metaData) <- c("sample", "group")

all_dds <- DESeqDataSetFromMatrix(countData=countData, 
                                  colData=metaData, 
                                  design=~group, tidy = TRUE)
all_dds <- estimateSizeFactors(all_dds)
keep_all <- rowSums(counts(all_dds, normalized=TRUE) >= 10 ) >= 3
all_dds <- all_dds[keep_all,]
all_dds <- DESeq(all_dds)



# PCA on all tissues

plotPCA.adapted <-  function(object, intgroup="condition", ntop, returnData=FALSE, pcs, samples) 
  # function adapted from https://support.bioconductor.org/p/69737/, James. W. McDonald (plotPCA.mystyle)
{ 
  stopifnot(length(pcs) == 2)    ### added this to check number of PCs ####
  # calculate the variance for each gene
  rv <- rowVars(assay(object), useNames = TRUE)
  print(head(rv))
  # select the ntop genes by variance...here, we use all genes
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,samples]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[samples, intgroup, drop=FALSE])
  # add the intgroup factors together to create a new grouping factor
  print(intgroup.df)
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]][samples]
  }
  print(group)
  # assemble the data for the plot
  ########## Here we just use the pcs object passed by the end user ####
  d <- data.frame(first=pca$x[1:nrow(pca$x),pcs[1]], second=pca$x[1:nrow(pca$x),pcs[2]], group=group, intgroup.df, name=colnames(object)[samples])
  print(d)
   #d$tissue <- c(rep("Heart", 17),
  #          rep("Lung", 17),
  #         rep("Kidney", 17),
  #        rep("Liver", 17))
  d$treatment <- rep(c(rep("EGF", 4), rep("not", 3), rep("EGF", 4), rep("not", 3), rep("not", 3)), 4)
  d$timepoint <- rep(c(rep("EGF_2h", 4), rep("no_EGF", 3), rep("EGF_6h", 4), rep("no_EGF", 3), rep("no_EGF", 3)), 4)
  #d$treatment <- rep(c(rep("EGF_2h", 4), rep("not", 3), rep("not", 4), rep("not", 3), rep("not", 3)), 4)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcs]
    return(d)
  }
}

generic_plot_pca <- function(dds_object, the_pcs, the_samples, firstpc, secondpc, name){
  vsd <- vst(dds_object, blind=FALSE)
  pcaData <- plotPCA.adapted(vsd, intgroup="group", ntop = nrow(dds_object), 
                             returnData=TRUE, pcs = the_pcs, samples = the_samples)
  all_percentVar <- round(100 * attr(pcaData, "percentVar"), 2)
  g <- ggplot(pcaData, aes(first, second, color=group)) +
    geom_point(size=3) + 
    xlab(paste0(firstpc,": ",all_percentVar[1],"% variance")) +
    ylab(paste0(secondpc,": ",all_percentVar[2],"% variance")) + 
    coord_fixed() + theme_bw() + theme(aspect.ratio= 8/10)
  pdf(paste0(figure_path, name, "_", firstpc, "_", secondpc, "_PCA_plot.pdf"))
  plot(g)
  dev.off()
}

generic_plot_pca(all_dds, c(1,3), c(1:68), "PC1", "PC3", "whole_tissue_RNA")
generic_plot_pca(all_dds, c(1,2), c(1:68), "PC1", "PC2", "whole_tissue_RNA")
generic_plot_pca(all_dds, c(2,3), c(1:68), "PC2", "PC3", "whole_tissue_RNA")

# specifically for treatment
generic_plot_pca_for_treatment <- function(dds_object, the_pcs, the_samples, firstpc, secondpc, name){
  vsd <- vst(dds_object, blind=FALSE)
  pcaData <- plotPCA.adapted(vsd, intgroup="group", ntop = nrow(dds_object), 
                                  returnData=TRUE, pcs = the_pcs, samples = the_samples)
  all_percentVar <- round(100 * attr(pcaData, "percentVar"), 2)
  g <- ggplot(pcaData, aes(first, second, color=treatment)) +
    geom_point(size=3) + 
    xlab(paste0(firstpc,": ",all_percentVar[1],"% variance")) +
    ylab(paste0(secondpc,": ",all_percentVar[2],"% variance")) + 
    coord_fixed() + theme_bw() + theme(aspect.ratio= 8/10)
  pdf(paste0(figure_path, name, "_", firstpc, "_", secondpc, "_PCA_plot.pdf"))
  plot(g)
  dev.off()
}

# specifically for timepoint
generic_plot_pca_for_timepoint <- function(dds_object, the_pcs, the_samples, firstpc, secondpc, name){
  vsd <- vst(dds_object, blind=FALSE)
  pcaData <- plotPCA.adapted(vsd, intgroup="group", ntop = nrow(dds_object), 
                                  returnData=TRUE, pcs = the_pcs, samples = the_samples)
  all_percentVar <- round(100 * attr(pcaData, "percentVar"), 2)
  g <- ggplot(pcaData, aes(first, second, color=timepoint)) +
    geom_point(size=3) + 
    xlab(paste0(firstpc,": ",all_percentVar[1],"% variance")) +
    ylab(paste0(secondpc,": ",all_percentVar[2],"% variance")) + 
    coord_fixed() + theme_bw() + theme(aspect.ratio= 8/10)
  pdf(paste0(figure_path, name, "_", firstpc, "_", secondpc, "_PCA_plot.pdf"))
  plot(g)
  dev.off()
}


# PCA on all tissues for timepoint and treatment
generic_plot_pca_for_timepoint(all_dds, c(4,6), c(1:68), "PC4", "PC6", "whole_tissue_RNA")
generic_plot_pca_for_treatment(all_dds, c(7,8), c(1:68), "PC7", "PC8", "whole_tissue_RNA")

# make heatmap to correlate PCs with conditions
get_PCs_for_cor <-  function(object, intgroup="condition", ntop, returnData=FALSE, pc)
{
  print("yes")
  stopifnot(length(pc) == 1)    ### added this to check number of PCs ####
  # calculate the variance for each gene
  print("yes")
  rv <- rowVars(assay(object), useNames = TRUE)
  print("yes")
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  print("yes")
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  print(intgroup.df)
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  ########## Here we just use the pcs object passed by the end user ####
  d <- data.frame(comp = pca$x[,pc], group=group, intgroup.df, name=colnames(object))
  #d$tissue <- c(rep("Heart", 17),
  #         rep("Lung", 17),
  #         rep("Kidney", 17),
  #         rep("Liver", 17))
  #d$treatment <- rep(c(rep("EGF_2h", 4), rep("veh_2h", 3), rep("EGF_6h", 4), rep("veh_6h", 3), rep("unt", 3)), 4)
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc]
    return(pca$x[,pc])
  }
}


all_vsd <- vst(all_dds, blind = FALSE)
all_pc_list <- vector(mode = "list", length = 30)
for(i in 1:30){ # get first 30 PCs
  pc_determinant <- get_PCs_for_cor(all_vsd, intgroup="group", ntop = nrow(all_dds), 
                                    returnData=TRUE, pc = i)
  all_pc_list[[i]] <- pc_determinant
}
# make DF of all conditions of interest
df_for_all_cor <- data.frame(  heart = c(rep(1, 17), rep(0, 17), rep(0, 17), rep(0, 17)),
                           lung = c(rep(0, 17), rep(1, 17), rep(0, 17), rep(0, 17)),
                           kidney = c(rep(0, 17), rep(0, 17), rep(1, 17), rep(0, 17)),
                           liver = c(rep(0, 17), rep(0, 17), rep(0, 17), rep(1, 17)),
                           egf_2h = rep(c(rep(1, 4), rep(0, 3), rep(0, 4), rep(0, 3), rep(0, 3)), 4),
                           veh_2h = rep(c(rep(0, 4), rep(1, 3), rep(0, 4), rep(0, 3), rep(0, 3)), 4),
                           egf_6h = rep(c(rep(0, 4), rep(0, 3), rep(1, 4), rep(0, 3), rep(0, 3)), 4),
                           veh_6h = rep(c(rep(0, 4), rep(0, 3), rep(0, 4), rep(1, 3), rep(0, 3)), 4),
                           unt = rep(c(rep(0, 4), rep(0, 3), rep(0, 4), rep(0, 3), rep(1, 3)), 4),
                           egf =  rep(c(rep(1, 4), rep(0, 3), rep(1, 4), rep(0, 3), rep(0, 3)), 4),
                           non_egf =  rep(c(rep(0, 4), rep(1, 3), rep(0, 4), rep(1, 3), rep(1, 3)), 4))
list_of_all_pc_cor <- vector(mode = "list", length = 30)
for(i in 1:30){ # correlate PCs with conditions of interest
  list_of_all_pc_cor[[i]] <- cor(as.numeric(all_pc_list[[i]]), df_for_all_cor, method = "spearman")
}

list_of_all_pc_cor.z <- t(do.call(rbind, list_of_all_pc_cor))
colnames(list_of_all_pc_cor.z) <- sprintf("PC%d", 1:30)
list_of_all_pc_cor.z <- na.omit(list_of_all_pc_cor.z)

draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 0.5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
paletteLength <- 20
myColor <- colorRampPalette(c("darkblue","white", "darkred"))(paletteLength)

bk1 <- c(seq(-1.0, -0.01, by =0.1), -0.0000010001)
bk2 <- c(0.00000010001, seq(.01, 1.0, by = 0.1))
bk <- c(bk1, bk2)

pdf(paste0(figure_path, "all_tissue_pc_heatmap.pdf"))
pheatmap(list_of_all_pc_cor.z, cluster_cols = F, cluster_rows=F, border_color = "black",
         color=myColor, breaks=bk, na_col = "grey", 
         #gaps_col = c(17,34, 51) #, gaps_row = c(length(shared_cell_cycle_sig_res),
         #length(heart_alone_cycle_sig_res),
         #length(liver_alone_cycle_sig_res))
)
dev.off()



# make tissue-specific dds objects

liver_dds <- get_dds_object("Liv")
liver_res <- get_dds_results(liver_dds)

lung_dds <- get_dds_object("Lu")
lung_res <- get_dds_results(lung_dds)

kid_dds <- get_dds_object("Kid")
kid_res <- get_dds_results(kid_dds)

heart_dds <- get_dds_object("He")
heart_res <- get_dds_results(heart_dds)


# individual tissue and time-point volcano plots

volcano_plot(heart_res[[1]], "heart_2h_")
volcano_plot(lung_res[[1]], "lung_2h_")
volcano_plot(kid_res[[1]], "kidney_2h_")
volcano_plot(liver_res[[1]], "liver_2h_")

volcano_plot(heart_res[[2]], "heart_6h_")
volcano_plot(lung_res[[2]], "lung_6h_")
volcano_plot(kid_res[[2]], "kidney_6h_")
volcano_plot(liver_res[[2]], "liver_6h_")



# get significant DEG lists for each tissue

sig_liver_2h <- get_sig_res(liver_res[[1]])
sig_lung_2h <- get_sig_res(lung_res[[1]])
sig_kidney_2h <- get_sig_res(kid_res[[1]])
sig_heart_2h <- get_sig_res(heart_res[[1]])
sig_liver_6h <- get_sig_res(liver_res[[2]])
sig_lung_6h <- get_sig_res(lung_res[[2]])
sig_kidney_6h <- get_sig_res(kid_res[[2]])
sig_heart_6h <- get_sig_res(heart_res[[2]])

csv_path <- paste0(figure_path, "DEGs/")
write_csv(sig_liver_2h, paste0(csv_path, "liver_2h_DEG.csv"))
write_csv(sig_liver_6h, paste0(csv_path, "liver_6h_DEG.csv"))
write_csv(sig_lung_2h, paste0(csv_path, "lung_2h_DEG.csv"))
write_csv(sig_lung_6h, paste0(csv_path, "lung_6h_DEG.csv"))
write_csv(sig_kidney_2h, paste0(csv_path, "kidney_2h_DEG.csv"))
write_csv(sig_kidney_6h, paste0(csv_path, "kidney_6h_DEG.csv"))
write_csv(sig_heart_2h, paste0(csv_path, "heart_2h_DEG.csv"))
write_csv(sig_heart_6h, paste0(csv_path, "heart_6h_DEG.csv"))

# Look at overlaps between 2h genes and 6h DEGs for each tissue

sig_gene_upset_2h <-list(Heart_2h = sig_heart_2h$ENSEMBL_ID, 
                          Lung_2h = sig_lung_2h$ENSEMBL_ID, 
                          Kidney_2h = sig_kidney_2h$ENSEMBL_ID,
                          Liver_2h = sig_liver_2h$ENSEMBL_ID)
                          #Heart_6h = sig_heart_6h$ENSEMBL_ID, 
                          #Lung_6h = sig_lung_6h$ENSEMBL_ID, 
                          #Kidney_6h = sig_kidney_6h$ENSEMBL_ID,
                          #Liver_6h = sig_liver_6h$ENSEMBL_ID)
sig_gene_upset_6h <-list(Heart_6h = sig_heart_6h$ENSEMBL_ID, 
  Lung_6h = sig_lung_6h$ENSEMBL_ID, 
  Kidney_6h = sig_kidney_6h$ENSEMBL_ID,
  Liver_6h = sig_liver_6h$ENSEMBL_ID)

pdf(paste0(figure_path, "sig_gene_2h_upset.pdf"),
    width = 12,
    height = 8)

upset(fromList(sig_gene_upset_2h), order.by = "degree", number.angles = 90,
      text.scale = c(2.5,2,1.5,1.5, 2.5, 2.5), mainbar.y.label = "DEGs intersection",
      sets.x.label = "DEGs per tissue/time-point",  keep.order = TRUE,
      matrix.color = "black", sets.bar.color = "black", main.bar.color = "black", nsets = 10)
dev.off()

pdf(paste0(figure_path, "sig_gene_6h_upset.pdf"),
    width = 12,
    height = 8)

upset(fromList(sig_gene_upset_6h), order.by = "degree", number.angles = 90,
      text.scale = c(2.5,2,1.5,1.5, 2.5, 2.5), mainbar.y.label = "DEGs intersection",
      sets.x.label = "DEGs per tissue/time-point",  keep.order = TRUE,
      matrix.color = "black", sets.bar.color = "black", main.bar.color = "black", nsets = 10)
dev.off()


make_2h_vs_6h_scatterplots(liver_res, sig_liver_2h, sig_liver_6h, "liver")
make_2h_vs_6h_scatterplots(lung_res, sig_lung_2h, sig_lung_6h, "lung")
make_2h_vs_6h_scatterplots(kid_res, sig_kidney_2h, sig_kidney_6h, "kidney")
make_2h_vs_6h_scatterplots(heart_res, sig_heart_2h, sig_heart_6h, "heart") 


# do all GO analyses

up_liver_2h_GO_mf <- up_do_GO_analysis_on_sig_genes_mf(sig_liver_2h, liver_dds)
up_liver_6h_GO_mf <- up_do_GO_analysis_on_sig_genes_mf(sig_liver_6h, liver_dds)
up_lung_2h_GO_mf <- up_do_GO_analysis_on_sig_genes_mf(sig_lung_2h, lung_dds)
#lung_6h_GO_bp <- do_GO_analysis_on_sig_genes_bp(sig_lung_6h, lung_dds)
up_kid_2h_GO_mf <- up_do_GO_analysis_on_sig_genes_mf(sig_kidney_2h, kid_dds)
#up_kid_6h_GO_bp <- up_do_GO_analysis_on_sig_genes_bp(sig_kidney_6h, kid_dds)
up_heart_2h_GO_mf <- up_do_GO_analysis_on_sig_genes_mf(sig_heart_2h, heart_dds)
#heart_6h_GO_bp <- do_GO_analysis_on_sig_genes_bp(sig_heart_6h, heart_dds)


write_csv(up_liver_2h_GO_mf@result, paste0(csv_path, "up_liver_2h_GO_mf.csv"))
write_csv(up_lung_2h_GO_mf@result, paste0(csv_path, "up_lung_2h_GO_mf.csv"))
write_csv(up_kid_2h_GO_mf@result, paste0(csv_path, "up_kidney_2h_GO_mf.csv"))
write_csv(up_heart_2h_GO_mf@result, paste0(csv_path, "up_heart_2h_GO_mf.csv"))

# make heatplots

pdf(paste0(figure_path, "up_liver_2h_GO_mf.pdf"))
heatplot(up_liver_2h_GO_mf, showCategory=10)
dev.off()

pdf(paste0(figure_path, "up_liver_6h_GO_mf.pdf"))
heatplot(up_liver_6h_GO_mf, showCategory=10)
dev.off()

pdf(paste0(figure_path, "up_lung_2h_GO_mf.pdf"))
heatplot(up_lung_2h_GO_mf, showCategory=5)
dev.off()

pdf(paste0(figure_path, "up_lung_2h_GO_mf_20.pdf"))
heatplot(up_lung_2h_GO_mf, showCategory=20)
dev.off()

#pdf(paste0(figure_path, "lung_6h_GO_mf.pdf"))
#heatplot(lung_6h_GO_bp, showCategory=10)
#dev.off()

pdf(paste0(figure_path, "up_kid_2h_GO_mf.pdf"))
heatplot(up_kid_2h_GO_mf, showCategory=5)
dev.off()
pdf(paste0(figure_path, "up_kid_2h_GO_mf_20.pdf"))
heatplot(up_kid_2h_GO_mf, showCategory=20)
dev.off()


#pdf(paste0(figure_path, "up_kid_6h_GO_bp.pdf"))
#heatplot(up_kid_6h_GO_bp, showCategory=10)
#dev.off()

pdf(paste0(figure_path, "up_heart_2h_GO_mf.pdf"))
heatplot(up_heart_2h_GO_mf, showCategory=5)
dev.off()
pdf(paste0(figure_path, "up_heart_2h_GO_mf)20.pdf"))
heatplot(up_heart_2h_GO_mf, showCategory=20)
dev.off()


