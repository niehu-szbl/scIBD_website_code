# Functions----
getPackage <- function(pkg, check = TRUE, load = TRUE, silent = FALSE, github = NULL, bioc = NULL) {
    if (check) {
      if (!suppressMessages(suppressWarnings(require(
        pkg, character.only = TRUE, quietly = TRUE
      )))) {
        if (is.null(github) & is.null(bioc)) {
          try(install.packages(pkg), silent = TRUE)
        }
        else if(github){
          try(remotes::install_github(github))
        }else if(bioc){
          if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
          try(BiocManager::install(bioc), silent = TRUE)
        }
      }
    }
    if (load)
      suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE))
    if (load & !silent)
      message("Loaded ", pkg)
  }

# get gene expression
get_gex <- function(ad, cells = NULL,  genes = NULL ){
  # row are cells
  # columns are genes
  
  if(class(ad)[1] == "Seurat"){
    return( GetAssayData(subset(ad, cells = cells, features = genes), slot = "data", assay = "RNA")
            %>% as.data.frame() %>% t())
  }else{
    return(ad[cells, genes]$X %>% as.data.frame())
  }
}

run_gsea = function(minor_cluster_list, expressed_genes, input_geneset, markers_df){
  # minor_cluster_list: all minor clusters
  # expressed_genes: expressed genes in seurat object, genes expressed in at least 200 cells
  # input_geneset: input gene set
  # markers_df: result of FindMarkers with down-sampled data
  
  all_pval = c()
  all_ods = c()
  
  # test
  for( i in 1:length(minor_cluster_list) ){
    cs_deg = markers_df[markers_df$cluster == minor_cluster_list[i], ]$gene
    go.obj <- newGeneOverlap(cs_deg, input_geneset, genome.size=length(expressed_genes))
    go.obj <- testGeneOverlap(go.obj)
    pval = getPval(go.obj)
    ods = getOddsRatio(go.obj)
    all_pval = c(all_pval, pval)
    all_ods = c(all_ods, ods)
  }
  
  # adjust p-value
  all_padj = round(p.adjust(all_pval, "BH"), 5)
  
  # round
  all_pval = round(all_pval, 5)
  all_ods = round(all_ods, 5)
  
  # cancer scientific number
  all_pval = format(all_pval, scientific = F)
  all_padj = format(all_padj, scientific = F)
  all_ods = format(all_ods, scientific = F)
  
  # generate result
  res = cbind("Cell type" = minor_cluster_list, 
              "P-value" = all_pval, 
              "Adjusted P-value" = all_padj,
              "Odds ratio"= all_ods)
  return(res)
}

plot_heatmap.minor_cluster = function(ad, cells = NULL, genes = NULL, mj_color = NULL){
  
  
  # get gene expression
  expression = get_gex(ad, cells = cells, genes = genes)
  
  # add label
  expression$label = ad$obs[cells, ]$minor_cluster
  expression$sample = ad$obs[cells, ]$sample
  expression = reshape2::melt(expression)
  colnames(expression) <- c("minor_cluster","sample","gene","expression")
  
  # expression, long format
  avg_exp = expression %>% 
    group_by(gene, minor_cluster, sample) %>% 
    summarise(avg_exp = mean(expm1(expression))) %>% 
    group_by(gene, minor_cluster) %>%
    summarise(avg_exp = log1p(mean(avg_exp)))
  
  # expression, matrix
  mtx = reshape2::dcast(avg_exp,  gene ~ minor_cluster, value.var = "avg_exp")
  rownames(mtx) = mtx[,1]
  mtx = mtx[,-1]
  mtx = apply(mtx, MARGIN = 1, function(x) rescale(x, to=c(0,1)) ) %>% t %>% as.data.frame() %>% mutate(gene = rownames(.))
  
  # expression, long format, scaled
  avg_exp_scaled = reshape2::melt(mtx)
  # avg_exp_scaled = mtx %>% reshape2::melt
  colnames(avg_exp_scaled) <- c("gene","minor_cluster","avg_exp")
  
  major_minor_mapping = ad$obs %>%
    select(minor_cluster, major_cluster) %>% 
    unique %>% 
    mutate(group = major_cluster, p = "group")
  
  p1 = ggplot(avg_exp_scaled, aes(x = minor_cluster, y = gene)) +
    geom_tile(aes(fill = avg_exp), color = "white", lwd = 0.25, linetype = 0)+
    scale_fill_viridis_c("Expression", option = "viridis") +
    xlab("")+
    ylab("")+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 8,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 8,  color = "black", face = "italic"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12, color = "black"),
          axis.title.y = element_blank(),  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent', color='white'),
          panel.border = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"),
          legend.position = "right") + 
    geom_vline(xintercept = 0.5 + major_minor_mapping$major_cluster %>% table %>% cumsum,  size = 0.5, color = "white")
  
  p2 = ggplot(major_minor_mapping, aes(y = p, x = minor_cluster, fill = major_cluster)) + 
    geom_tile(aes(fill = major_cluster)) +
    scale_fill_manual(values = mj_color) + 
    scale_y_discrete(position = "right")+ 
    xlab(NULL)+
    ylab(NULL) + 
    theme_minimal() + 
    labs(fill = "Cell type") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  mtx = mtx[,-ncol(mtx)]
  
  phr <- hclust(dist(mtx)) %>% 
    ggtree::ggtree()
  
  p = p1 %>% insert_top(p2, height = 0.05) %>% insert_left(phr, width = 0.1)
  
  # clean
  rm(avg_exp_scaled, mtx)
  
  return(p)
}