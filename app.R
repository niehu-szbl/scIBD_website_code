# Load functions----
source("./utils.R")

# Load packages
if(T){
  packages <-
    c("shiny",
      "shinycssloaders",
      "shinythemes",
      "shinydashboard",
      "shinyWidgets",
      "shinyjs",
      "shinyBS",
      'downloadthis',
      "bsplus",
      "DT",
      "dplyr",
      "reshape2",
      "data.table",
      "ggplot2",
      "ggpubr",
      "aplot",
      "ggtree",
      "plotly",
      "ggrepel",
      "RColorBrewer",
      "stringr",
      "readxl",
      "networkD3",
      "GeneOverlap",
      "scales")
  lapply(packages, getPackage)
  
  packages <- c("SCopeLoomR", 
                "anndata",
                "SCENIC",
                "Seurat")
  lapply(packages, getPackage)
}

# Define variables and load data
if(T){
  
  # connect to data
  adata = read_h5ad("./www/scanpy/all.clean.h5ad")
  adata_obs = adata$obs
  
  # get genes
  all_genes = adata$var %>% rownames %>% sort
    
  # get major cluster and minor cluster
  if(T){
    # get major cluster, major_cluster
    major_cluster = adata_obs[,"major_cluster"] %>% levels
    
    # get minor cluster, all_cluster_list
    all_cluster_list = list()
    for(i in 1:length(major_cluster)){
      minor_cluster =  adata_obs[ adata_obs$major_cluster == major_cluster[i], ]$minor_cluster %>% 
        droplevels  %>% levels
      all_cluster_list[[i]] = minor_cluster
    }
    names(all_cluster_list) = major_cluster
  }
  
  # set colors
  if(T){
    # set color for consecutive values
    #scIBD_consecutive_color = brewer.pal.info %>% filter(category != "qual") %>% rownames
    scIBD_consecutive_color = brewer.pal.info %>% rownames
    
    # set color for major cluster
    major_cluster_color = read.table("./www/color/major_cluster_color.txt", header = F, 
                                     stringsAsFactors = F, comment.char = ">", row.names = 1)
    color_name = rownames(major_cluster_color)
    major_cluster_color = major_cluster_color[,1]
    names(major_cluster_color) = color_name
    
    # set color for minor cluster
    minor_cluster_color = read.table("./www/color/minor_cluster_color.txt", header = F, 
                                     stringsAsFactors = F, comment.char = ">", sep = "\t", row.names = 1)
    color_name = rownames(minor_cluster_color)
    minor_cluster_color = minor_cluster_color[,1]
    names(minor_cluster_color) = color_name
  }

  # get meta data, get group information
  if(T){
    # get stage, disease, tissue, study
    all_stage_list = levels(adata_obs$stage)
    all_disease_list = levels(adata_obs$disease)
    all_tissue_list = levels(adata_obs$tissue)
    all_study_list = levels(adata_obs$study)
    
    # get sample list
    all_sample_list = list()
    for(i in 1:length(all_study_list)){
      samples = adata_obs[adata_obs$study == all_study_list[i], ]$sample %>% 
        droplevels  %>% levels
      all_sample_list[[i]] = samples
    }
    names(all_sample_list) = all_study_list
    
    # get tissue.sub list
    all_tissue.sub_list = list()
    for(i in 1:length(all_tissue_list)){
      tissue.sub = adata_obs[adata_obs$tissue == all_tissue_list[i], ]$tissue.sub %>% 
        droplevels  %>% levels
      all_tissue.sub_list[[i]] = tissue.sub
    }
    names(all_tissue.sub_list) = all_tissue_list
  }
  
  # load projects information
  if(T){ 
    project_info = read_xlsx("./www/meta/projects_IBD.xlsx",sheet=1) %>% as.data.frame()
  }
  
  # load drugs, targets, clinical trails for IBD
  if(T){
    # load data, FDA approved drugs
    therapy_panel.FDA_approved_drugs = read_xlsx("./www/therapy/FDA_approved_biologics_and_novel_small_molecules_for_IBD.xlsx",sheet=1)
    therapy_panel.FDA_approved_drugs = data.frame(therapy_panel.FDA_approved_drugs, check.names = FALSE)
    rownames(therapy_panel.FDA_approved_drugs) = therapy_panel.FDA_approved_drugs[,1]
    
    # details of FDA approved drugs
    therapy_panel.FDA_approved_drugs_detail = read_xlsx("./www/therapy/Pivotal_clinical_trails_of_IBD.xlsx",sheet=1)
    therapy_panel.FDA_approved_drugs_detail = data.frame(therapy_panel.FDA_approved_drugs_detail, check.names = FALSE)
    rownames(therapy_panel.FDA_approved_drugs_detail) = paste0("more_",1:nrow(therapy_panel.FDA_approved_drugs_detail))
    
    # therapy drugs and targets for IBD
    therapy_panel.drug_and_target = readRDS("./www/therapy/EFO_0003767-known-drugs.modified.rds")
  }
  
  # load gwas study
  if(T){
    ibd_gwas_data = read_xlsx("./www/gwas/IBD_GWAS.xlsx")
    gwas_study = paste( gsub(pattern = ".*>([a-zA-Z0-9]+)<.*", "\\1", ibd_gwas_data$Study),
                        gsub(pattern = ".*>([a-zA-Z0-9]+)<.*", "\\1", ibd_gwas_data$Catalog),
                        ibd_gwas_data$Disease,
                        sep = ",")
    gwas_study_list = list()
    for(i in 1:length(gwas_study)){ gwas_study_list[[i]] = i}
    names(gwas_study_list) = gwas_study
    
    # get genes
    catalogs = gsub(pattern = ".*>([a-zA-Z0-9]+)<.*", "\\1", ibd_gwas_data$Catalog)
    gwas_gene_list = list()
    for(i in 1:length(catalogs)){
      tmp = read.table( paste0("./www/gwas/", catalogs[i], ".gene.txt"), header = F, stringsAsFactors = F, sep = "\t")
      tmp = tmp$V1
      gwas_gene_list[[i]] = tmp
    }
  }
  
  # load risk gene
  if(T){
    adult_ibd_risk_genes_data = read.table("./www/gwas/gsct.extend.txt", sep = "\t", header = T, strip.white = F)
    adult_ibd_risk_genes_data$description = str_replace(adult_ibd_risk_genes_data$description, "\\[.*", "")
    
    pediatric_ibd_risk_genes_data = read.table("./www/gwas/pibd.extend.txt", sep = "\t", header = T, strip.white = F)
    pediatric_ibd_risk_genes_data$description = str_replace(pediatric_ibd_risk_genes_data$description, "\\[.*", "")
  }
  
  # load deg data
  deg_data = readRDS("./www/deg/diff_gex.rds")
  
  # load scenic data
  if(T){
    loom_file = read.table("./www/scenic/loom.txt", header = F, sep = "\t")
    loom_data = list()
    for(i in 1:nrow(loom_file)){
      loom_data[[i]] = open_loom(loom_file[i,2])
    }
    names(loom_data) = loom_file[,1]
  }
  
  # load differential regulon table
  diff_regulon_data = readRDS("./www/scenic/diff_regulon.rds")
  
  # load tf list
  if(T){
    tf_list = readRDS("./www/scenic/tf_list.rds")
  }
  
  # load rss data
  if(T){
    rss_file = read.table("./www/scenic/rss.txt", header = F, sep = "\t")
    rss_data = list()
    for(i in 1:nrow(rss_file)){
      rss_data[[i]] = read.table(rss_file[i,2], sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    }
    names(rss_data) = rss_file[,1]
  }
  
  # load data for gsea
  if(T){
    markers_df = read.table("./www/gsea/scIBD_markers.downsampled.txt", header = T, stringsAsFactors = F, sep = "\t")
    expressed_genes = read.table("./www/gsea/scIBD.expressed_genes.txt", header = F, stringsAsFactors = F, sep = "\t")
    expressed_genes = expressed_genes$V1
  }
  
  # load sample meta data
  if(T){
    sample_meta_data = read.csv("./www/meta/ibd_meta_data.csv", header = T, stringsAsFactors = F)
  }
}

# Define UI page
# global layout: navbarPage
# lay out of IBD page: dashboardPage
if(T){
  # Shiny UI: Gene expression profile panel
  if(T){
    GEX_profile_panel = fluidRow(
      column( width = 12,
        fluidRow(
          box(title = "Options", solidHeader = T, width = 12, collapsible = F,
              status = "primary",
              column( width = 4,
                      pickerInput( inputId = "GEX_profile_panel.major_cluster", 
                                   label = "Major cluster", 
                                   choices = major_cluster,
                                   selected = major_cluster[1],
                                   multiple = FALSE,
                                   options = list(`actions-box` = TRUE, 
                                                  `live-search` = TRUE,
                                                  size = 8,
                                                  style = "background:white; color:black")),
                      pickerInput(
                        inputId = "GEX_profile_panel.Gene_list",
                        label = "Select genes",
                        choices = all_genes,
                        selected = c("TPSAB1","CPA3"),
                        multiple = TRUE,
                        options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8, 
                                       style = "background:white; color:black"))
                  ),
              column( width = 4,
                      pickerInput(
                        inputId = "GEX_profile_panel.Embedding_used",
                        label = "Embedding used",
                        choices = c("UMAP", "tSNE"),
                        selected = "UMAP",
                        multiple = FALSE,
                        options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8, 
                                       style = "background:white; color:black")),
                      pickerInput(
                        inputId = "GEX_profile_panel.Colorpanel",
                        label = "Color profile",
                        choices = scIBD_consecutive_color,
                        selected = "Reds",
                        multiple = FALSE,
                        options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8, 
                                       style = "background:white; color:black"))
                  ),
              column( width = 4,
                      sliderInput(inputId = "GEX_profile_panel.Embedding_dotsize", 
                                  label = "Dot size",
                                  value = 0.1, 
                                  min = 0, max = 1,step = 0.05),
                      actionButton(inputId = "GEX_profile_panel.Submit", 
                                   label = "GO", 
                                   width = "100px", 
                                   icon = icon("paper-plane"),
                                   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")
                  ))),
        fluidRow(
          box(title = "Dot plot of gene expression", solidHeader = T, width = 4, collapsible = T,
              status = "primary",
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_profile_panel.Scanpy_embedding.plot_exp")),
              
                dropdownButton(
                  tags$h5("Confiugre download"),
                  textInput(inputId = 'GEX_profile_panel.Scanpy_embedding.plot_exp.title', 
                            label = 'Title',
                            value = 'Dot plot of gene expression'),
                  sliderInput(inputId = 'GEX_profile_panel.Scanpy_embedding.plot_exp.width',
                              label = 'Width',
                              value = 5, min = 1, max = 10),
                  sliderInput(inputId = 'GEX_profile_panel.Scanpy_embedding.plot_exp.height',
                              label = 'Height',
                              value = 5, min = 1, max = 10),
                  selectInput(
                    inputId = "GEX_profile_panel.Scanpy_embedding.plot_exp.file_type",
                    label = "Type", 
                    choices = c("pdf","jpeg"),
                  ),
                  circle = TRUE, status = "info", icon = icon("cog"), 
                  size = 'sm', inline = TRUE, up = TRUE,
                  tooltip = tooltipOptions(title = "Click to see inputs !")
                ),
                
                downloadBttn(outputId = 'GEX_profile_panel.Scanpy_embedding.plot_exp.download_pdf',
                            label = 'Download',
                            style = 'minimal', color = 'primary', 
                            block = FALSE, size = 'sm', no_outline = TRUE)
            ),
          box(title = "Annotation of cell subsets", solidHeader = T, width = 4, collapsible = T,
              status = "primary",
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_profile_panel.Scanpy_embedding.plot_label")),
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_profile_panel.Scanpy_embedding.plot_label.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_profile_panel.Scanpy_embedding.plot_label.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_profile_panel.Scanpy_embedding.plot_label.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_profile_panel.Scanpy_embedding.plot_label.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              downloadBttn(outputId = 'GEX_profile_panel.Scanpy_embedding.plot_label.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)
              ),
          box(title="Barplot of cell numbers", solidHeader = T, width = 4, collapsible = T,
              status = "primary",
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_profile_panel.subset_cell_number.barplot")),
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_profile_panel.subset_cell_number.barplot.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_profile_panel.subset_cell_number.barplot.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_profile_panel.subset_cell_number.barplot.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_profile_panel.subset_cell_number.barplot.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              downloadBttn(outputId = 'GEX_profile_panel.subset_cell_number.barplot.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)
              )
          ),
        fluidRow(
          box(title = "Violin plot of gene expression", solidHeader = T, width = 6, collapsible = T,
              status = "primary",
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_profile_panel.Violin_plot")),
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_profile_panel.Violin_plot.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_profile_panel.Violin_plot.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_profile_panel.Violin_plot.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_profile_panel.Violin_plot.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              downloadBttn(outputId = 'GEX_profile_panel.Violin_plot.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)
          ),
          box(title = "Dot plot of gene expression", solidHeader = T, width = 6, collapsible = T,
              status = "primary",
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_profile_panel.Heatmap_plot")),
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_profile_panel.Heatmap_plot.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_profile_panel.Heatmap_plot.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_profile_panel.Heatmap_plot.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_profile_panel.Heatmap_plot.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              downloadBttn(outputId = 'GEX_profile_panel.Heatmap_plot.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)
              )
          ),
        fluidRow(
          box(title="Heatmap plot of marker genes", solidHeader = T, width = 12, collapsible = T,
              status = "primary",
              selectInput(inputId = "GEX_profile_panel.deg_topn", 
                          label = "Number of top marker genes: ", 
                          choices = 1:10, selected = 4, width = "300px"),
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_profile_panel.marker_gene.heatmap")),
              
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_profile_panel.marker_gene.heatmap.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_profile_panel.marker_gene.heatmap.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_profile_panel.marker_gene.heatmap.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_profile_panel.marker_gene.heatmap.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              downloadBttn(outputId = 'GEX_profile_panel.marker_gene.heatmap.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)
              )
            ),
      fluidRow(
        box(title = "Marker genes of each cell subtype", solidHeader = T, width = 12, collapsible = T,
            status = "primary",
            DT::DTOutput("GEX_profile_panel.marker_gene_tbl"), style = "font-size: 100%;"))
      )
    )
  }

  # Shiny UI: Regulon profile panel
  if(T){
    GRN_profile_panel = fluidRow(
      column( width = 12,
              fluidRow(
                box(title = "Options", solidHeader = T, width = 12, collapsible = T,
                    status = "primary",
                    column(width = 4,
                           pickerInput(
                             inputId = "GRN_profile_panel.major_cluster",
                             label = "Major cluster", 
                             choices = major_cluster,
                             selected = major_cluster[1],
                             multiple = FALSE,
                             options = list(`actions-box` = TRUE, 
                                            `live-search` = TRUE,
                                            size = 8,
                                            style = "background:white; color:black")),
                           pickerInput(
                               inputId = "GRN_profile_panel.TF",
                               label = "Transcript factor",
                               choices = tf_list[[ major_cluster[1] ]],
                               selected = tf_list[[ major_cluster[1] ]][1:3],
                               multiple = TRUE,
                               options = list(`actions-box` = TRUE, 
                                              `live-search` = TRUE, 
                                              size = 8, 
                                              style = "background:white; color:black"))),
                    column( width = 4,
                            pickerInput(
                              inputId = "GRN_profile_panel.Embedding_used",
                              label = "Embedding used",
                              choices = c("UMAP", "tSNE"),
                              selected = "UMAP",
                              multiple = FALSE,
                              options = list(`actions-box` = TRUE, 
                                             `live-search` = TRUE, 
                                             size = 8, 
                                             style = "background:white; color:black")),
                            pickerInput(
                              inputId = "GRN_profile_panel.Colorpanel",
                              label = "Color profile",
                              choices = scIBD_consecutive_color,
                              selected = "Reds",
                              multiple = FALSE,
                              options = list(`actions-box` = TRUE, 
                                             `live-search` = TRUE, 
                                             size = 8, 
                                             style = "background:white; color:black"))),
                    column(width = 4,
                      sliderInput(inputId = "GRN_profile_panel.Embedding_dotsize", label = "Dot size",
                                  value = 0.1, 
                                  min = 0, max = 2, step = 0.1),
                      actionButton(inputId = "GRN_profile_panel.Submit", label = "GO", width = "100px", 
                                   icon = icon("paper-plane"),
                                   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
                )),
              fluidRow(
                # display regulon activity score of given genes
                box(title = "Regulon activity", solidHeader = T, width = 4, collapsible = T,
                    status = "primary",
                    shinycssloaders::withSpinner(
                      plotOutput(outputId = "GRN_profile_panel.Scanpy_embedding.plot_activity")),
                    
                    dropdownButton(
                      tags$h5("Confiugre download"),
                      textInput(inputId = 'GRN_profile_panel.Scanpy_embedding.plot_activity.title', 
                                label = 'Title',
                                value = 'Dot plot of gene expression'),
                      sliderInput(inputId = 'GRN_profile_panel.Scanpy_embedding.plot_activity.width',
                                  label = 'Width',
                                  value = 5, min = 1, max = 10),
                      sliderInput(inputId = 'GRN_profile_panel.Scanpy_embedding.plot_activity.height',
                                  label = 'Height',
                                  value = 5, min = 1, max = 10),
                      selectInput(
                        inputId = "GRN_profile_panel.Scanpy_embedding.plot_activity.file_type",
                        label = "Type", 
                        choices = c("pdf","jpeg"),
                      ),
                      circle = TRUE, status = "info", icon = icon("cog"), 
                      size = 'sm', inline = TRUE, up = TRUE,
                      tooltip = tooltipOptions(title = "Click to see inputs !")
                    ),
                    
                    downloadBttn(outputId = 'GRN_profile_panel.Scanpy_embedding.plot_activity.download_pdf',
                                 label = 'Download',
                                 style = 'minimal', color = 'primary', 
                                 block = FALSE, size = 'sm', no_outline = TRUE)
                  ),
                
                # display cell subset annotation, scanpy embedding
                box(title = "Annotation of cell subsets", solidHeader = T, width = 4, collapsible = T,
                    status = "primary",
                    shinycssloaders::withSpinner(
                      plotOutput(outputId = "GRN_profile_panel.Scanpy_embedding.plot_label")),
                  
                    dropdownButton(
                      tags$h5("Confiugre download"),
                      textInput(inputId = 'GRN_profile_panel.Scanpy_embedding.plot_label.title', 
                                label = 'Title',
                                value = 'Dot plot of gene expression'),
                      sliderInput(inputId = 'GRN_profile_panel.Scanpy_embedding.plot_label.width',
                                  label = 'Width',
                                  value = 5, min = 1, max = 10),
                      sliderInput(inputId = 'GRN_profile_panel.Scanpy_embedding.plot_label.height',
                                  label = 'Height',
                                  value = 5, min = 1, max = 10),
                      selectInput(
                        inputId = "GRN_profile_panel.Scanpy_embedding.plot_label.file_type",
                        label = "Type", 
                        choices = c("pdf","jpeg"),
                      ),
                      circle = TRUE, status = "info", icon = icon("cog"), 
                      size = 'sm', inline = TRUE, up = TRUE,
                      tooltip = tooltipOptions(title = "Click to see inputs !")
                    ),
                    
                    downloadBttn(outputId = 'GRN_profile_panel.Scanpy_embedding.plot_label.download_pdf',
                                 label = 'Download',
                                 style = 'minimal', color = 'primary', 
                                 block = FALSE, size = 'sm', no_outline = TRUE)
                  ),
                
                # display SCENIC embedding
                box(title = "SCENIC embedding", solidHeader = T, width = 4, collapsible = T,
                    status = "primary",
                    shinycssloaders::withSpinner(
                      plotOutput(outputId = "GRN_profile_panel.Scenic_embedding.plot_label")),
                  
                    dropdownButton(
                      tags$h5("Confiugre download"),
                      textInput(inputId = 'GRN_profile_panel.Scenic_embedding.plot_label.title', 
                                label = 'Title',
                                value = 'Dot plot of gene expression'),
                      sliderInput(inputId = 'GRN_profile_panel.Scenic_embedding.plot_label.width',
                                  label = 'Width',
                                  value = 5, min = 1, max = 10),
                      sliderInput(inputId = 'GRN_profile_panel.Scenic_embedding.plot_label.height',
                                  label = 'Height',
                                  value = 5, min = 1, max = 10),
                      selectInput(
                        inputId = "GRN_profile_panel.Scenic_embedding.plot_label.file_type",
                        label = "Type", 
                        choices = c("pdf","jpeg"),
                      ),
                      circle = TRUE, status = "info", icon = icon("cog"), 
                      size = 'sm', inline = TRUE, up = TRUE,
                      tooltip = tooltipOptions(title = "Click to see inputs !")
                    ),
                    
                    downloadBttn(outputId = 'GRN_profile_panel.Scenic_embedding.plot_label.download_pdf',
                                 label = 'Download',
                                 style = 'minimal', color = 'primary', 
                                 block = FALSE, size = 'sm', no_outline = TRUE)
                  )),
              fluidRow(
                ## place violin plot and heatmap
                box(title = "Violin plot of regulon activity", solidHeader = T, width = 6, collapsible = T,
                    status = "primary",
                    shinycssloaders::withSpinner(
                      plotOutput(outputId = "GRN_profile_panel.Violin_plot")),
                  
                    dropdownButton(
                      tags$h5("Confiugre download"),
                      textInput(inputId = 'GRN_profile_panel.Violin_plot.title', 
                                label = 'Title',
                                value = 'Dot plot of gene expression'),
                      sliderInput(inputId = 'GRN_profile_panel.Violin_plot.width',
                                  label = 'Width',
                                  value = 5, min = 1, max = 10),
                      sliderInput(inputId = 'GRN_profile_panel.Violin_plot.height',
                                  label = 'Height',
                                  value = 5, min = 1, max = 10),
                      selectInput(
                        inputId = "GRN_profile_panel.Violin_plot.file_type",
                        label = "Type", 
                        choices = c("pdf","jpeg"),
                      ),
                      circle = TRUE, status = "info", icon = icon("cog"), 
                      size = 'sm', inline = TRUE, up = TRUE,
                      tooltip = tooltipOptions(title = "Click to see inputs !")
                    ),
                    
                    downloadBttn(outputId = 'GRN_profile_panel.Violin_plot.download_pdf',
                                 label = 'Download',
                                 style = 'minimal', color = 'primary', 
                                 block = FALSE, size = 'sm', no_outline = TRUE)  
                  ),
                box(title = "Heatmap plot of regulon activity", solidHeader = T, width = 6, collapsible = T,
                    status = "primary",
                    shinycssloaders::withSpinner(
                      plotOutput(outputId = "GRN_profile_panel.Heatmap_plot")),
                  
                    dropdownButton(
                      tags$h5("Confiugre download"),
                      textInput(inputId = 'GRN_profile_panel.Heatmap_plot.title', 
                                label = 'Title',
                                value = 'Dot plot of gene expression'),
                      sliderInput(inputId = 'GRN_profile_panel.Heatmap_plot.width',
                                  label = 'Width',
                                  value = 5, min = 1, max = 10),
                      sliderInput(inputId = 'GRN_profile_panel.Heatmap_plot.height',
                                  label = 'Height',
                                  value = 5, min = 1, max = 10),
                      selectInput(
                        inputId = "GRN_profile_panel.Heatmap_plot.file_type",
                        label = "Type", 
                        choices = c("pdf","jpeg"),
                      ),
                      circle = TRUE, status = "info", icon = icon("cog"), 
                      size = 'sm', inline = TRUE, up = TRUE,
                      tooltip = tooltipOptions(title = "Click to see inputs !")
                    ),
                    
                    downloadBttn(outputId = 'GRN_profile_panel.Heatmap_plot.download_pdf',
                                 label = 'Download',
                                 style = 'minimal', color = 'primary', 
                                 block = FALSE, size = 'sm', no_outline = TRUE)
                  )),
              
              fluidRow(
                box(title = "Network of regulons", solidHeader = T, width = 6, collapsible = T,
                    status = "primary",
                    shinycssloaders::withSpinner(
                      forceNetworkOutput(outputId = "GRN_profile_panel.plot_network")),
                  
                    dropdownButton(
                      tags$h5("Confiugre download"),
                      textInput(inputId = 'GRN_profile_panel.plot_network.title', 
                                label = 'Title',
                                value = ''),
                      sliderInput(inputId = 'GRN_profile_panel.plot_network.width',
                                  label = 'Width',
                                  value = 5, min = 1, max = 10),
                      sliderInput(inputId = 'GRN_profile_panel.plot_network.height',
                                  label = 'Height',
                                  value = 5, min = 1, max = 10),
                      selectInput(
                        inputId = "GRN_profile_panel.plot_network.file_type",
                        label = "Type", 
                        choices = c("html"),
                      ),
                      circle = TRUE, status = "info", icon = icon("cog"), 
                      size = 'sm', inline = TRUE, up = TRUE,
                      tooltip = tooltipOptions(title = "Click to see inputs !")
                    ),
                    
                    downloadBttn(outputId = 'GRN_profile_panel.plot_network.download_pdf',
                                 label = 'Download',
                                 style = 'minimal', color = 'primary', 
                                 block = FALSE, size = 'sm', no_outline = TRUE)  
                  ),
                box(title = "Network of regulons", solidHeader = T, width = 6, collapsible = T,
                    status = "primary",
                    # htmlOutput(outputId = "GRN_profile_panel.describe_network")
                    shinycssloaders::withSpinner(
                      uiOutput(outputId = "GRN_profile_panel.describe_network"))
                    )),
              fluidRow(
                box(title = "Regulon specificity scores (RSS) in each cell type", solidHeader = T, width = 12, collapsible = T,
                    status = "primary",
                    shinycssloaders::withSpinner(DT::DTOutput("GRN_profile_panel.rss_tbl"))))
    ))
  }
  
  # Shiny UI: Gene expression comparison panel
  if(T){
    GEX_comparison_panel = fluidRow(
      column(width = 12,
        fluidRow(
          box(title = "Select dataset", solidHeader = T, width = 12, collapsible = T,
              status = "primary",
              column(width = 2,
                     pickerInput( inputId = "GEX_comparison_panel.Gene_list",
                                  label = "Select genes ", 
                                  choices = sort(all_genes),
                                  selected = c("TPSAB1","CPA3"),
                                  multiple = TRUE,
                                  options = list(`actions-box` = TRUE, 
                                                 `live-search` = TRUE, 
                                                 size = 8, 
                                                 style = "background:white; color:black")),
                     pickerInput(
                       inputId = "GEX_comparison_panel.Embedding_used",
                       label = "Embedding used",
                       choices = c("UMAP", "tSNE"), 
                       selected = "UMAP",
                       multiple = FALSE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8, 
                                      style = "background:white; color:black")),
                     actionButton(inputId = "GEX_comparison_panel.Submit", 
                           label = "GO", 
                           width = "100px", 
                           icon = icon("paper-plane"),
                           style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")),
              column(width = 2,
                     pickerInput(
                       inputId = "GEX_comparison_panel.major_cluster",
                       label = "Major cluster",
                       choices = major_cluster,
                       selected = major_cluster,
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8, 
                                      style = "background:white; color:black")),
                     pickerInput(
                       inputId = "GEX_comparison_panel.minor_cluster",
                       label = "Minor cluster",
                       choices = all_cluster_list,
                       selected = unlist(all_cluster_list),
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8,
                                      style = "background:white; color:black")),
                     ),
              column(width = 2,
                     pickerInput(
                       inputId = "GEX_comparison_panel.tissue",
                       label = "Tissue",
                       choices = all_tissue_list,
                       selected = all_tissue_list,
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8,
                                      style = "background:white; color:black")),
                     pickerInput(
                       inputId = "GEX_comparison_panel.tissue.sub",
                       label = "Location",
                       choices = all_tissue.sub_list,
                       selected = unlist(all_tissue.sub_list),
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8,
                                      style = "background:white; color:black")),
                     ),
              column(width = 2,
                     pickerInput(
                       inputId = "GEX_comparison_panel.stage",
                       label = "Development stage",
                       choices = all_stage_list,
                       selected = all_stage_list,
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8,
                                      style = "background:white; color:black")),
                     pickerInput(
                       inputId = "GEX_comparison_panel.disease",
                       label = "Disease state",
                       choices = all_disease_list,
                       selected = all_disease_list,
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8,
                                      style = "background:white; color:black")),
                     ),
              column(width = 2,
                     pickerInput(
                       inputId = "GEX_comparison_panel.study",
                       label = "Study",
                       choices = all_study_list,
                       selected = all_study_list,
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8,
                                      style = "background:white; color:black")),
                     pickerInput(
                       inputId = "GEX_comparison_panel.sample",
                       label = "Sample",
                       choices = all_sample_list,
                       selected = unlist(all_sample_list),
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8,
                                      style = "background:white; color:black")),
                     ),
              column(width = 2,
                     pickerInput(
                       inputId = "GEX_comparison_panel.Colorpanel",
                       label = "Color profile",
                       choices = scIBD_consecutive_color,
                       selected = "Reds",
                       multiple = FALSE,
                       options = list(`actions-box` = TRUE, 
                                      `live-search` = TRUE, 
                                      size = 8,
                                      style = "background:white; color:black")),
                     sliderInput(inputId = "GEX_comparison_panel.Embedding_dotsize", 
                                 label = "Dot size",
                                 value = 0.1, 
                                 min = 0, max = 1,step = 0.05),
                     )
        )),
        fluidRow(
          # plot umap/tsne, all cells, color by average gene expression----
          box(title = "Embedding plot of gene expression", solidHeader = T, width = 4, collapsible = T,
              status = "primary",
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_comparison_panel.Scanpy_embedding.plot_exp")),
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_exp.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_exp.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_exp.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_comparison_panel.Scanpy_embedding.plot_exp.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              downloadBttn(outputId = 'GEX_comparison_panel.Scanpy_embedding.plot_exp.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)
            ),
          
          # plot umap/tsne, all cells, color by major cluster----
          box(title = "Annotation of all cells", solidHeader = T, width = 4, collapsible = T,
              status = "primary",
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_comparison_panel.Scanpy_embedding.plot_label")),
              
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_label.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_label.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_label.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_comparison_panel.Scanpy_embedding.plot_label.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              
              downloadBttn(outputId = 'GEX_comparison_panel.Scanpy_embedding.plot_label.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)  
            ),
          
          # plot umap/tsne, selected cells, color by selected group----
          box(title = "Annotation of selected cells", solidHeader = T, width = 4, collapsible = T,
              status = "primary",
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_comparison_panel.Scanpy_embedding.plot_cell")),
              
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_cell.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_cell.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_comparison_panel.Scanpy_embedding.plot_cell.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_comparison_panel.Scanpy_embedding.plot_cell.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              
              downloadBttn(outputId = 'GEX_comparison_panel.Scanpy_embedding.plot_cell.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)  
          ),
        ),
        fluidRow(
          # plot gene expression, violin plot, group by one variable ----
          box(title = "Violin plot of gene expression", solidHeader = T, width = 6, collapsible = T,
              status = "primary",
              fluidRow(
                column(width = 6,
                       selectInput(inputId = "GEX_comparison_panel.groupby", 
                                   label = "Group by", 
                                   choices = c("major_cluster","minor_cluster", "stage", "tissue", "tissue.sub", "disease",'study','sample'),
                                   selected = "major_cluster", 
                                   width = "300px"))
              ),
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_comparison_panel.Violin_plot")),
            
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_comparison_panel.Violin_plot.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_comparison_panel.Violin_plot.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_comparison_panel.Violin_plot.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_comparison_panel.Violin_plot.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              
              downloadBttn(outputId = 'GEX_comparison_panel.Violin_plot.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)
            ),
              
          # plot cell number, bar plot, group by one variable ----
          box(title = "Cells in each group", solidHeader = T, width = 6, collapsible = T,
              status = "primary",
              fluidRow(
                column(width = 6,
                       selectInput(inputId = "GEX_comparison_panel.Cell_number.Bar_plot.type", 
                                   label = "Select value to display", 
                                   choices = c("Number","Percentage"),
                                   selected = "Percentage", 
                                   width = "300px")),
              ),
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_comparison_panel.Cell_number.Bar_plot")),
            
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_comparison_panel.Cell_number.Bar_plot.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_comparison_panel.Cell_number.Bar_plot.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_comparison_panel.Cell_number.Bar_plot.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_comparison_panel.Cell_number.Bar_plot.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              
              downloadBttn(outputId = 'GEX_comparison_panel.Cell_number.Bar_plot.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)  
            ),
        ),
        fluidRow(
          # plot gene expression, violin plot, group by two variable ----
          box(title = "Violin plot of gene expression", solidHeader = T, width = 6, collapsible = T,
              status = "primary",
              fluidRow(
                column(width = 6,
                       selectInput(inputId = "GEX_comparison_panel.groupby.1", 
                                   label = "Group by.1", 
                                   #choices = c("major_cluster","minor_cluster"),
                                   choices = c("major_cluster","minor_cluster","stage","disease","tissue","tissue.sub","study","sample"),
                                   selected = "major_cluster", 
                                   width = "300px")),
                column(width = 6, 
                       selectInput(inputId = "GEX_comparison_panel.groupby.2", 
                                   label = "Group by.2", 
                                   #choices = c("stage","disease","tissue","tissue.sub","study","sample"),
                                   choices = c("major_cluster","minor_cluster","stage","disease","tissue","tissue.sub","study","sample"),
                                   selected = "disease", 
                                   width = "300px")),
              ),
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_comparison_panel.Violin_plot2")),
            
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_comparison_panel.Violin_plot2.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_comparison_panel.Violin_plot2.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_comparison_panel.Violin_plot2.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_comparison_panel.Violin_plot2.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              
              downloadBttn(outputId = 'GEX_comparison_panel.Violin_plot2.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE)  
            ),
          
          # plot cell number, bar plot, group by two variable ----
          box(title = "Cells in each group", solidHeader = T, width = 6, collapsible = T,
              status = "primary",
              fluidRow(
                column(width = 6,
                       selectInput(inputId = "GEX_comparison_panel.Bar_plot2.type", 
                                   label = "Select value to display", 
                                   choices = c("Number","Percentage"),
                                   selected = "Percentage", 
                                   width = "300px")),
                column(width = 6,
                       selectInput(inputId = "GEX_comparison_panel.Bar_plot2.faceting_group", 
                                   label = "Select faceting group", 
                                   # this need further fix
                                   choices = c("major_cluster","minor_cluster","stage","disease","tissue","tissue.sub","study","sample"),
                                   selected = "minor_cluster", 
                                   width = "300px")),
                
              ),
              shinycssloaders::withSpinner(
                plotOutput(outputId = "GEX_comparison_panel.Cell_number.Bar_plot2")),
            
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'GEX_comparison_panel.Cell_number.Bar_plot2.title', 
                          label = 'Title',
                          value = 'Dot plot of gene expression'),
                sliderInput(inputId = 'GEX_comparison_panel.Cell_number.Bar_plot2.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'GEX_comparison_panel.Cell_number.Bar_plot2.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "GEX_comparison_panel.Cell_number.Bar_plot2.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              
              downloadBttn(outputId = 'GEX_comparison_panel.Cell_number.Bar_plot2.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary',
                           block = FALSE, size = 'sm', no_outline = TRUE)
            )
        )
      ))
  }
  
  # Shiny UI: Regulon profile comparison
  if(T){
    GRN_comparison_panel = fluidRow(
      column( width = 12,
              fluidRow(
                box(title = "Select dataset", solidHeader = T, width = 12, collapsible = F,
                    status = "primary",
                    column(width = 2,
                           pickerInput(
                             inputId = "GRN_comparison_panel.major_cluster",
                             label = "Major cluster",
                             choices = major_cluster,
                             selected = major_cluster[1],
                             multiple = FALSE,
                             options = list(`actions-box` = TRUE, 
                                           `live-search` = TRUE, 
                                           size = 8, 
                                           style = "background:white; color:black")),
                          pickerInput(
                            inputId = "GRN_comparison_panel.TF",
                            label = "Transcript factor",
                            choices = NULL,
                            multiple = FALSE,
                            options = list(`actions-box` = TRUE, 
                                           `live-search` = TRUE, 
                                           size = 8, 
                                           style = "background:white; color:black"))),
                    column(width = 2,
                           pickerInput(
                             inputId = "GRN_comparison_panel.minor_cluster",
                             label = "Minor cluster",
                             choices = NULL,
                             multiple = TRUE,
                             options = list(`actions-box` = TRUE, 
                                     `live-search` = TRUE, 
                                     size = 8, 
                                     style = "background:white; color:black")),
                           pickerInput(
                             inputId = "GRN_comparison_panel.tissue",
                             label = "Tissue",
                             choices = all_tissue_list,
                             selected = all_tissue_list,
                             multiple = TRUE,
                             options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8,
                                       style = "background:white; color:black"))),
                    column(width = 2,
                           pickerInput(
                             inputId = "GRN_comparison_panel.tissue.sub",
                             label = "Location",
                             choices = all_tissue.sub_list,
                             selected = unlist(all_tissue.sub_list),
                             multiple = TRUE,
                             options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8,
                                       style = "background:white; color:black")),
                        pickerInput(
                          inputId = "GRN_comparison_panel.stage",
                          label = "Development stage",
                          choices = all_stage_list,
                          selected = all_stage_list,
                          multiple = TRUE,
                          options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8,
                                       style = "background:white; color:black"))),
                    column(width = 2,
                           pickerInput(
                             inputId = "GRN_comparison_panel.disease",
                             label = "Disease state",
                             choices = all_disease_list,
                             selected = all_disease_list,
                             multiple = TRUE,
                             options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8,
                                       style = "background:white; color:black")),
                        pickerInput(
                          inputId = "GRN_comparison_panel.study",
                          label = "Study",
                          choices = all_study_list,
                          selected = all_study_list,
                          multiple = TRUE,
                          options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8,
                                       style = "background:white; color:black"))),
                    column(width = 2,
                      pickerInput(
                        inputId = "GRN_comparison_panel.sample",
                        label = "Sample",
                        choices = all_sample_list,
                        selected = unlist(all_sample_list),
                        multiple = TRUE,
                        options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8,
                                       style = "background:white; color:black")),
                      pickerInput(
                        inputId = "GRN_comparison_panel.Colorpanel",
                        label = "Color profile",
                        choices = scIBD_consecutive_color,
                        selected = "Reds",
                        multiple = FALSE,
                        options = list(`actions-box` = TRUE, 
                                       `live-search` = TRUE, 
                                       size = 8, 
                                       style = "background:white; color:black"))
                      ),
                    column(width = 2,
                           tags$div(tags$p(tags$b("Submit"))),
                           actionButton(inputId = "GRN_comparison_panel.Submit", 
                                        label = "GO", 
                                        width = "100px", 
                                        icon = icon("paper-plane"),
                                        style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
                  )
                ),
              fluidRow(
                  box(title = "Violin plot of regulon activity", solidHeader = T, width = 6, collapsible = F,
                      status = "primary",
                      fluidRow(
                        column(width = 3, 
                               selectInput(inputId = "GRN_comparison_panel.groupby", 
                                           label = "Group by", 
                                           choices = c("stage","disease","tissue","tissue.sub","minor_cluster"),
                                           selected = "disease",
                                           width = "300px"))),
                      shinycssloaders::withSpinner(
                        plotOutput(outputId = "GRN_comparison_panel.Violin_plot")),
                      
                      dropdownButton(
                        tags$h5("Confiugre download"),
                        textInput(inputId = 'GRN_comparison_panel.Violin_plot.title', 
                                  label = 'Title',
                                  value = 'Dot plot of gene expression'),
                        sliderInput(inputId = 'GRN_comparison_panel.Violin_plot.width',
                                    label = 'Width',
                                    value = 5, min = 1, max = 10),
                        sliderInput(inputId = 'GRN_comparison_panel.Violin_plot.height',
                                    label = 'Height',
                                    value = 5, min = 1, max = 10),
                        selectInput(
                          inputId = "GRN_comparison_panel.Violin_plot.file_type",
                          label = "Type", 
                          choices = c("pdf","jpeg"),
                        ),
                        circle = TRUE, status = "info", icon = icon("cog"), 
                        size = 'sm', inline = TRUE, up = TRUE,
                        tooltip = tooltipOptions(title = "Click to see inputs !")
                      ),
                      
                      downloadBttn(outputId = 'GRN_comparison_panel.Violin_plot.download_pdf',
                                   label = 'Download',
                                   style = 'minimal', color = 'primary', 
                                   block = FALSE, size = 'sm', no_outline = TRUE)
                    ),
                  box(title = "Cells in each group", solidHeader = T, width = 6, collapsible = F,
                      status = "primary",
                      fluidRow(
                        column(width = 3, 
                               selectInput(inputId = "GRN_comparison_panel.Cell_number.Bar_plot.type", 
                                           label = "Select value to display", 
                                           choices = c("Number","Percentage"),
                                           selected = "Percentage", 
                                           width = "300px"))),
                      shinycssloaders::withSpinner(
                        plotOutput(outputId = "GRN_comparison_panel.Cell_number.Bar_plot")),
                      
                      dropdownButton(
                        tags$h5("Confiugre download"),
                        textInput(inputId = 'GRN_comparison_panel.Cell_number.Bar_plot.title', 
                                  label = 'Title',
                                  value = 'Dot plot of gene expression'),
                        sliderInput(inputId = 'GRN_comparison_panel.Cell_number.Bar_plot.width',
                                    label = 'Width',
                                    value = 5, min = 1, max = 10),
                        sliderInput(inputId = 'GRN_comparison_panel.Cell_number.Bar_plot.height',
                                    label = 'Height',
                                    value = 5, min = 1, max = 10),
                        selectInput(
                          inputId = "GRN_comparison_panel.Cell_number.Bar_plot.file_type",
                          label = "Type", 
                          choices = c("pdf","jpeg"),
                        ),
                        circle = TRUE, status = "info", icon = icon("cog"), 
                        size = 'sm', inline = TRUE, up = TRUE,
                        tooltip = tooltipOptions(title = "Click to see inputs !")
                      ),
                      
                      downloadBttn(outputId = 'GRN_comparison_panel.Cell_number.Bar_plot.download_pdf',
                                   label = 'Download',
                                   style = 'minimal', color = 'primary', 
                                   block = FALSE, size = 'sm', no_outline = TRUE)  
                    )
                ),
              fluidRow(
                box(title = "Compare regulons between healthy and UC",
                    solidHeader = T, width = 12, collapsible = F, status = "primary",
                    style = "font-size: 100%",
                    DT::DTOutput("GRN_comparison_panel.diff_regulon.healthy_UC.tbl", 
                                 width = "100%", height = "auto")
                    )
              ),
              fluidRow(
                box(title = "Compare regulons between healthy and CD",
                    solidHeader = T, width = 12, collapsible = F, status = "primary",
                    style = "font-size: 100%",
                    DT::DTOutput("GRN_comparison_panel.diff_regulon.healthy_CD.tbl",
                                 width = "100%", height = "auto"))
              )
        )
    )
  }
  
  # Shiny UI: Cellular composition
  if(T){
    cellular_composition_panel = fluidRow(
      column( width = 12,
              fluidRow(
                box(title = "Select samples", solidHeader = T, width = 12, collapsible = F,
                    status = "primary",
                    column(width = 4, 
                           pickerInput(
                             inputId = "cellular_composition_panel.study",
                             label = "Study",
                             choices = NULL,
                             multiple = TRUE,
                             options = list(`actions-box` = TRUE, 
                                            `live-search` = TRUE, 
                                            size = 8,
                                            style = "background:white; color:black")),
                           ),
                    column(width = 4,
                           pickerInput(
                             inputId = "cellular_composition_panel.sample",
                             label = "Sample",
                             choices = NULL,
                             multiple = TRUE,
                             options = list(`actions-box` = TRUE, 
                                            `live-search` = TRUE, 
                                            size = 8,
                                            style = "background:white; color:black"))
                           ),
                    column(width = 4,
                           tags$div(tags$p(tags$b("Submit"))),
                           actionButton(inputId = "cellular_composition_panel.Submit", 
                                        label = "GO", 
                                        width = "100px", 
                                        height = "40px",
                                        icon = icon("paper-plane"),
                                        style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                           ))),
              #
              fluidRow(
                  box(title = "Cell compositions of major subsets", solidHeader = T, width = 12, collapsible = TRUE,
                      status = "primary",
                      plotlyOutput(outputId = "cellular_composition_panel.output.major_barplot",
                                   width = "100%"))),
              
              fluidRow(
                box(title = "Cell compositions of sub-cell types", solidHeader = T, width = 12, collapsible = TRUE,
                    status = "primary",
                    plotlyOutput(outputId = "cellular_composition_panel.output.minor_barplot",
                                 width = "100%")))
              )
  )
}
  
  # Shiny UI: Gene set enrichment analysis
  if(T){
    gsea_panel = fluidRow(
          box(title = "Input gene set", solidHeader = T, width = 12, collapsible = F,
              status = "primary",
              tabBox( id = "Gsea_input_panel",
                      side = "left", selected = "Genes", width = 12,
                      tabPanel(title = "Genes", 
                               textAreaInput( inputId = "Gsea_panel.input.text",
                                              label = "Input a gene set:", 
                                              value = NULL) %>%
                                 shinyInput_label_embed(
                                   shiny_iconlink() %>%
                                     bs_embed_popover(
                                       title = "Note:",
                                       placement = "left",
                                       content = "Please place one gene per line"
                                     )
                                 )
                               ),
                      tabPanel(title = "Upload",
                               textAreaInput(inputId = "Gsea_panel.file.text",
                                             label = "Check genes in uploaded file:",
                                             value = NULL),
                               fileInput(
                                 inputId = "Gsea_panel.input.file",
                                 label = "Select a txt file:",
                                 multiple = FALSE ) %>%
                                 shinyInput_label_embed(
                                   shiny_iconlink() %>%
                                     bs_embed_popover(
                                       title = "Note:",
                                       placement = "left",
                                       content = "Please place one gene per line in uploaded file"
                                     )
                                 )
                               ),
                      tabPanel(title = "Select",
                               textAreaInput(inputId = "Gsea_panel.select.text",
                                             label = "Check genes in selected studies:",
                                             value = NULL),
                               checkboxGroupInput("Gsea_panel.select.predefined", 
                                                  label = "Pre-defined gene sets", 
                                                  choices = gwas_study_list,
                                                  selected = 1),
                               )
              ),
              shiny::actionButton(
                inputId = "Gsea_panel.Submit",
                label = "GO",
                width = "80px",
                icon = icon("paper-plane"),
                style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
              )),
          box(title = "Output of enrichment analysis", solidHeader = T, width = 12, collapsible = TRUE,
              status = "primary",
              div(DT::DTOutput("Gsea_panel.output.tbl"), style = "font-size: 100%;")),
          box(title = "Heatmap", solidHeader = T, width = 12, collapsible = TRUE,
              status = "primary",
              dropdownButton(
                tags$h5("Confiugre download"),
                textInput(inputId = 'Gsea_panel.output.heatmap.title', 
                          label = 'Title',
                          value = ''),
                sliderInput(inputId = 'Gsea_panel.output.heatmap.width',
                            label = 'Width',
                            value = 5, min = 1, max = 10),
                sliderInput(inputId = 'Gsea_panel.output.heatmap.height',
                            label = 'Height',
                            value = 5, min = 1, max = 10),
                selectInput(
                  inputId = "Gsea_panel.output.heatmap.file_type",
                  label = "Type", 
                  choices = c("pdf","jpeg"),
                ),
                circle = TRUE, status = "info", icon = icon("cog"), 
                size = 'sm', inline = TRUE, up = FALSE, right = TRUE,
                tooltip = tooltipOptions(title = "Click to see inputs !")
              ),
              
              downloadBttn(outputId = 'Gsea_panel.output.heatmap.download_pdf',
                           label = 'Download',
                           style = 'minimal', color = 'primary', 
                           block = FALSE, size = 'sm', no_outline = TRUE),
              shinycssloaders::withSpinner(
                plotOutput(outputId = "Gsea_panel.output.heatmap", height = "2000px", width = "100%"))
            ),
    )
  }
  
  # Shiny UI: IBD targets panel
  if(T){
    IBD_targets_panel = fluidRow(
      box(title = "FDA approved drugs for IBD", solidHeader = T, width = 12, collapsible = F,
          status = "primary",
          div(
            # FDA approved table
            DT::DTOutput("therapy_panel.FDA_approved_tbl"), style = "font-size: 90%;"),
          tags$hr(),
          tags$p("Route: iv (intravenous) represents administration within or into a vein or veins; sc (subcutaneous) represents administration beneath the skin; po (Oral) represents administration to or by way of the mouth"),
          tags$p("CD represents Crohn’s disease; UC represents ulcerative colitis; x represents that the drug can be used to treat UC or CD")
      ),
      box(title = "Therapy targets and drugs for IBD", solidHeader = T, width = 12, collapsible = F,
          status = "primary",
          div(DT::DTOutput("therapy_panel.ibd_targets_tbl"), style = "font-size: 90%;"),
          tags$hr(),
          tags$p("Drugs and targets for IBD are retrived from ", 
                 tags$a(href='https://www.opentargets.org', target = "_blank", "Open Targets database"))
          )
    )
  }
  
  # Shiny UI: IBD risk genes panel
  if(T){
    IBD_risk_genes_panel = fluidRow(
      box(title = "Major GWAS study on IBD", solidHeader = T, width = 12, collapsible = F,
          status = "primary",
          div(DT::DTOutput("ibd_gwas_tbl"), style = "font-size: 100%;")),
      box(title = "Risk genes of IBD patients", solidHeader = T, width = 12, collapsible = F,
          status = "primary",
          div(DT::DTOutput("adult_ibd_risk_genes_tbl"), style = "font-size: 100%;")),
      box(title = "Risk genes of pediatric patients", solidHeader = T, width = 12, collapsible = F,
          status = "primary",
          div(DT::DTOutput("pediatric_ibd_risk_genes_tbl"), style = "font-size: 100%;"))
    )
  }
  
  # Shiny UI: Meta data panel
  if(T){
    meta_data_panel = fluidRow(
      box(title = "Meta data of IBD datasets", solidHeader = T, width = 12, collapsible = F,
          status="primary",
          div(DT::DTOutput("meta_data_panel"), style = "font-size: 70%;")),
      box(title = "Note", solidHeader = T, width = 12, collapsible = F,
          status="primary",
          div(
            p("location1: sample location (blood, largeInt, smallInt)"),
            p("location2: sample location (e.g. colon, ileum, rectum)"),
            p("stage: development stage of the sample (e.g. adult, fetal)"),
            p("data_source: accession of the data source"),
            p("study: study name"),
            ))
    )
  }
}

# Shiny UI
ui <- navbarPage(id="nav", 
                 theme = shinytheme("flatly"),
                 collapsible = FALSE,
                 windowTitle = "scIBD",

# Create navbarPage----
              ## set style
             HTML('<a style="text-decoration:none;
                             cursor:default;
                             color:#FFFFFF;
                             font-family: Arial, Georgia, Times, serif;
                             font-weight: bold;
                             font-size: 20px;" 
                  class="active" href="#">scIBD</a>'), 
# Home page----
             tabPanel("Home",
                      tags$div( style = "margin-left: 3%; margin-right: 3%;",
                        tags$h4("Single-Cell Transcriptomic Atlas of Human Inflammatory Bowel Disease"),
                        tags$p("Inflammatory bowel disease (IBD) is a type of chronic inflammation disease whose exact etiology is still unclear. \
                               With the increasing studies of IBD by single cell RNA sequencing technique (scRNA-seq), dysregulation of immune \
                               micro-environment and pathogenesis of IBD have been uncovered successively. The enormous IBD-related scRNA-seq \
                               datasets in the past decade are calling for a burning demand to be uniformly processed and integrated for conveniently access. \
                               Here, we developed a database of Single Cell transcriptomic atlas of Inflammatory Bowel Disease (scIBD) that contains ~1.14 million \
                               single cells across multiple development stages and disease states, comprising 9 major subtypes and 96 minor subtypes. \
                               We also collected clinical trials, therapy targets, GWAS-implicated risk genes to give a quick glance of advances in the treatment of IBD. \
                               Furthermore, we developed a multi-functional and user-friendly website that provides interactive visualization for biologists \
                               to analyse the transcriptome features, gene regulatory networks and enrichment of given gene set in each cell subset.")
                      ),
                      tags$div(
                        tags$h5("Overview of scIBD datasets and annotations"), align = "middle"
                      ),
                      tags$div( style = "width: 100%; display: inline-block; ", 
                          tags$div( align = "middle",
                            tags$img(src = "./img/umap.by_major_cluster.png"),
                            tags$img(src = "./img/umap.by_disease.png"),
                            tags$img(src = "./img/umap.by_stage.png"))
                      ),
                      tags$div( style = "width: 100%; display: inline-block; ", 
                          tags$div( align = "middle",
                            tags$img(src = "./img/umap.by_tissue.png"),
                            tags$img(src = "./img/umap.by_tissue.sub.png"),
                            tags$img(src = "./img/umap.by_dataset.png"))
                      ),
                      tags$div( style = "width: 100%; display: inline-block; ", 
                        tags$div( align = "middle",
                          tags$img(src = "./img/umap.by_minor_cluster.png"))
                      ),
                      tags$div( style = "width: 100%; display: inline-block; ", 
                                align = "middle",
                        tags$p("Mailing address", tags$br(),
                               "Institute of Cancer Research,", tags$br(),
                               "Shenzhen Bay Laboratory,", tags$br(),
                               "Guangming District, Shenzhen, Guangdong,", tags$br(),
                               "P.R. China"),
                        tags$p(tags$b("Developed by LabZhangLei | © Copyright 2022")))
                      ),

# IBD page----
             tabPanel("Exploration",
                        dashboardPage(
                          dashboardHeader(disable = T),
                          dashboardSidebar(
                            sidebarMenu(
                              menuItem(text = "Gene Expression Profile", tabName = "GEX_profile"),
                              menuItem(text = "Regulon Activity Profile", tabName = "GRN_profile"),
                              menuItem(text = "Gene Expression Comparison", tabName = "GEX_comparison"),
                              menuItem(text = "Regulon Activity Comparison", tabName = "GRN_comparison"),
                              menuItem(text = "Cellular composition", tabName = "cellular_composition"),
                              menuItem(text = "Gene Enrichment Analysis", tabName = "gsea"),
                              menuItem(text = "Current Therapy Strategy", tabName = "IBD_targets"),
                              menuItem(text = "GWAS-implicated Risk Genes", tabName = "IBD_risk_genes"),
                              menuItem(text = "Meta Data Exploration", tabName = "meta_data")
                            )),
                        dashboardBody(
                          ## css
                          ## .logo, .navbar
                          tags$head(tags$style(HTML('.skin-blue .main-sidebar .sidebar .sidebar-menu .active a{background-color: #3079AE; color: #FFFFFF;}'))),
                          tags$head(tags$style(HTML('.skin-blue .main-sidebar {
                            padding-top: 10px;
                            background-color: #212F3F;
                            font-family: "Arial", "Georgia", Times, "Times New Roman", serif;
                            font-weight: bold;
                            font-size: 14px;}'))),
                          tags$head(tags$style(HTML(".content-wrapper, .right-side{
                                          background-color: #f3f6f4;}"))),
                          
                          ## content
                          tabItems(
                            tabItem(tabName = "GEX_profile", GEX_profile_panel),
                            tabItem(tabName = "GRN_profile", GRN_profile_panel),
                            tabItem(tabName = "GEX_comparison", GEX_comparison_panel),
                            tabItem(tabName = "GRN_comparison", GRN_comparison_panel),
                            tabItem(tabName = "cellular_composition", cellular_composition_panel),
                            tabItem(tabName = "gsea", gsea_panel),
                            tabItem(tabName = "IBD_risk_genes", IBD_risk_genes_panel),
                            tabItem(tabName = "IBD_targets", IBD_targets_panel),
                            tabItem(tabName = "meta_data", meta_data_panel))))
                    ),
# Documents page----
             tabPanel("Documents"),

# References page----
             tabPanel("References",
                      tags$div( style = "margin-left: 2%; margin-right: 2%",
                        tags$h4("Browse Projects"),
                        DT::DTOutput("projects_info_table")
                      )),
# Downloads page----
             tabPanel("Downloads",
                      tags$div( style = "margin-left: 2%; margin-right: 2%",
                        tags$h3("Downloads"),
                        tags$p("Data files used in this website could be accessed through the following links."),
                        tags$h4("Gene Expression Matrix"),
                        downloadLink(outputId = "download_panel.download_gex", 
                                     label = "Integrated gene expression matrix of scIBD"),
                        tags$br(),
                        downloadLink(outputId = "download_panel.download_meta", 
                                     label = "Meta data of scIBD"),
                        
                        tags$h4("Differential Expression Genes"),
                        downloadLink(outputId = "download_panel.download_deg", 
                                     label = "Differential expressed genes within each major cluster compartment, and differential expressed genes across all minor clusters"),
                        
                        tags$h4("Differential Regulons Between Disease and Healthy"),
                        downloadLink(outputId = "download_panel.download_loom", 
                                     label = "Regulons within each major cluster"),
                        tags$br(),
                        downloadLink(outputId = "download_panel.download_rss", 
                                     label = "Regulon specificity score (RSS) of each minor cluster within major cluster compartment"),
                        tags$br(),
                        downloadLink(outputId = "download_panel.download_diff_regulon",
                                     label = "Differential regulons between healthy and disease within major cluster compartment")
                      )
                    ),
# Help page----
             tabPanel("Help",
                      tags$div( style = "margin-left: 2%; margin-right: 2%",
                        tags$h3("FAQ"),
                        tags$p(tags$b("Q1: What is scIBD?")),
                        tags$p("scIBD is a database of Single Cell Transcriptomic Atlas of Human Inflammatory Bowel Disease (scIBD) that \
                               contains ~1.14 million single cells from 12 datasets across multiple development stages (including fetal, pediatric, and adult), \
                               tissue locations (includign blood, small intestine and large intestine) and different disease states (healthy, inflammed UC, inflammed CD, etc.). \
                               scIBD comprises 9 major subtypes (including Myeloid, CD4+ T cells, CD8+ T cells, ILCs, B/Plasma cells, Epithelial cells, Mesenchynal cells, Endothelial cells, and Neuronal cells), and \
                               96 minor subtypes. scIBD provides a multi-functional and user-friendly interface that provides interactive visualization for biologists to \
                               analyse the transcriptome features, gene regulatory networks and enrichment of given gene set in each cell subset."),
                        tags$p(tags$b("Q2: What are the feature functions of scIBD?")),
                        tags$p("We have integrated 13 datasets from multiple studies which investigate the pathologies of IBD, and present a comprehensive single cell transcriptomic atlas for further studying IBD. \
                               With scIBD, users are convenient to explore marker genes of each cell subtype, and compare gene expression of given genes (such as therapy targets, cytokines, IBD-GWAS related genes, or others) between healthy and disease across major clusters or minor clusters. \
                               With scIBD, users are also convenient to explore the underlying gene regulatory networks (GRNs) of each minor cluster, and compare the activities of given regulons between healthy and disease. \
                               IBD is caused by a complex interaction between genetic and environment factors (such as gut microbes). \
                               Currently, treatments for IBD including 5-ASA, antibiotics, steroids, immunosuppressants, and biologic therapies (including anti–tumor necrosis factor [TNF] antibodies, anti–α4β7 integrin antibodies, and anti–IL12/23 antibodies). \
                               scIBD also collected clinical trials, therapy targets, and GWAS-implicated risk genes to give a quick glance of advances in the treatment of IBD. \
                               With scIBD, user could capture both the enriched cell subsets and gene expression profiles of the risk genes for UC and CD or any given gene set."),
                        tags$p(tags$b("Q3: How did you perform the cell type annotations?")),
                        tags$p("Preprocessing and integration of scRNA-seq datasets were performed with Scanpy. \
                                Raw count matrix of gene expression of all samples in all datasets were merged. \
                                We used a double-MAD (median absolute deviation) method to detected the outliers of number of expressed genes in all cells with lower threshold of 1.5 and upper threshold of 10. \
                                As a result, cells with fewer than 399 genes or more than 6,849 genes, or > 25% mitochondrial UMI counts were filtered out. Samples with less than 100 cells were dropped. \
                                Then, we used a two rounds clustering strategy to integrate, cluster and annotate major clusters and minor clusters."),
                        tags$p("In the first-round clustering, immunoglobulin genes, T cell receptor (TCR) genes, and ribosome-protein-coding genes (gene symbol with string pattern “^RP[0-9]+-|[LS]”), cell cycling genes (e.g. TOP2A, MKI67), and mitochondrial genes were removed from the combined gene expression matrix. \
                               For each cell, the UMI counts of genes were divided by the total UMI count of the cell and then scaled by 1e4, and then log-transformed. \
                               Top 2,000 highly variable genes were identified across cells and samples (with highly_variable_genes function with scanpy (v1.8.2, n_top_genes=2000 and batch_key = “sampleName”). \
                               These 2,000 highly variable genes were used for downstream analysis. The expression values of each gene were scaled to unit variance, and PCA were performed. \
                               Batch correction of datasets from multiple studies which covered three developmental stages and multiple disease conditions, was performed with bbknn (v1.5.1, bath_key=\"sampleName\", neighbors_within_batch=3, metric=\"euclidean\", n_pcs =30). \
                               Dimensionality reduction (t-SNE and UMAP) and leiden clustering (resolution=0.1) was performed and major clusters were annotated based on canonical marker genes and differentially expressed genes (DEGs)."),
                        tags$p("In the second-round, integration and clustering of given major cluster were similar to first-round. We checked the gene expression of canonical marker genes, and sub clusters expressed marker genes from other major clusters were considered as doublets and dropped. 
                               Then, each major cluster was integrated and clustered again, and minor clusters were annotated based on marker genes identified from DEGs and published studies."),
                      ),
                      tags$div( style= "margin-left: 2%; margin-right: 2%;",
                        tags$h3("Contact Details"),
                        
                        tags$h4("Scientific Problems"),
                        tags$p("We welcome any suggestions regarding how to improve our database, please feel free to contact us with feedback."),
                        tags$p("Please contact Hu Nie: niehu2021@163.com"),
                        
                        tags$h4("Technical Problems"),
                        tags$p("If you have any questions about the usage of scIBD and the interpretation of results, or encounter problems when using scIBD."),
                        tags$p("Please contact Hu Nie: niehu2021@163.com"),
                        
                        tags$h4("Address"),
                        tags$p("Gaoke Innovation Center, Guangqiao Road, Guangming District, Shenzhen"),
                        tags$p("Institute of Cancer Research"),
                        tags$p("Shenzhen Bay Laboratory"),
                        tags$p("Guangdong, China"),
                        tags$p("Tel: +86 (10) 8696-7710")
                      ))
)

# shiny server----
server <- function(input, output, session) {

  # Gene expression panel
  if(T){
      # set genes
      GEX_profile_panel.genes = reactiveValues(genes = c("TPSAB1","CPA3"))

      # get genes
      observeEvent(input$GEX_profile_panel.Submit,{
        GEX_profile_panel.genes$genes = input$GEX_profile_panel.Gene_list
      })
     
      # get cells
      GEX_profile_panel.data = reactive({
        
        # get cell name
        cells = rownames( adata_obs[adata_obs$major_cluster == input$GEX_profile_panel.major_cluster,] )
        
        # return
        list(cells = cells)
      })
      
      # get umap and tsne
      if(T){
        GEX_profile_panel.scanpy_embedding = reactive({
          
          # get umap
          umap = adata_obs[GEX_profile_panel.data()$cells, ] %>% 
            select("UMAP_1","UMAP_2","minor_cluster") %>% as.data.frame()
          
          # get tsne
          tsne = adata_obs[GEX_profile_panel.data()$cells, ] %>% 
            select("TSNE_1","TSNE_2","minor_cluster") %>% as.data.frame()

          # return
          list(umap = umap, tsne = tsne)
        })
      }
      
      # get gene expression of given genes and given cells
      if(T){ 
        GEX_profile_panel.gene_data = reactive({
          
          # get gene expression
          expression = get_gex(adata, 
                               cells = GEX_profile_panel.data()$cells, 
                               genes = GEX_profile_panel.genes$genes)
          
          # get average expression of input genes
          avg_exp = rowMeans(expression) %>% as.data.frame
          
          # return
          list(avg_exp = avg_exp, expression = expression)
        })
      }
      
      # Dot plot of gene expression
      if(T){
        GEX_profile_panel.Scanpy_embedding.plot_exp.plot = reactiveValues(plot = NULL)
        output$GEX_profile_panel.Scanpy_embedding.plot_exp <- renderPlot({
          # plot setting
          mycolor <- colorRampPalette(brewer.pal(8, input$GEX_profile_panel.Colorpanel))(100)
          dot.size = input$GEX_profile_panel.Embedding_dotsize
        
          # plot
          if(input$GEX_profile_panel.Embedding_used == "UMAP"){
            # prepare data
            umap_data = cbind(GEX_profile_panel.scanpy_embedding()$umap[,c(1,2)],
                              GEX_profile_panel.gene_data()$avg_exp)
            colnames(umap_data) = c("UMAP_1","UMAP_2","expression")
            
            # plot umap
            p= ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) + 
              geom_point(aes(color = expression), size = dot.size) + 
              scale_colour_gradientn("Exp", colors = mycolor) +
              xlab("UMAP_1") + 
              ylab("UMAP_2") + 
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='transparent', color='white'),
                    panel.border = element_blank(),
                    legend.position = "right",
                    axis.line = element_line(color = "black"))
            GEX_profile_panel.Scanpy_embedding.plot_exp.plot$plot = p
            print(p)
          }else if(input$GEX_profile_panel.Embedding_used == "tSNE"){
            # prepare data
            tsne_data = cbind(GEX_profile_panel.scanpy_embedding()$tsne[,c(1,2)],
                              GEX_profile_panel.gene_data()$avg_exp)
            colnames(tsne_data) = c("tSNE_1","tSNE_2","expression")
            
            # plot tsne
            p = ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2)) + 
              geom_point(aes(color = expression), size = dot.size) + 
              scale_colour_gradientn("Exp", colors = mycolor) + 
              xlab("tSNE_1") + 
              ylab("tSNE_2") + 
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='transparent', color='white'),
                    panel.border = element_blank(),
                    legend.position = "right",
                    axis.line = element_line(color = "black"))
            GEX_profile_panel.Scanpy_embedding.plot_exp.plot$plot = p
            print(p)
          }
        })
        
        # download GEX_profile_panel.Scanpy_embedding.plot_exp
        output$GEX_profile_panel.Scanpy_embedding.plot_exp.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_profile_panel.Scanpy_embedding.plot_exp.file_type == 'pdf'){
              paste0('GEX_profile_panel.Scanpy_embedding.plot_exp_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_profile_panel.Scanpy_embedding.plot_exp.file_type == 'jpeg'){
              paste0('GEX_profile_panel.Scanpy_embedding.plot_exp_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_profile_panel.Scanpy_embedding.plot_exp.file_type == 'pdf'){
              pdf(file, width = input$GEX_profile_panel.Scanpy_embedding.plot_exp.width, height = input$GEX_profile_panel.Scanpy_embedding.plot_exp.height)
              plot(GEX_profile_panel.Scanpy_embedding.plot_exp.plot$plot)
              dev.off()
            }else if(input$GEX_profile_panel.Scanpy_embedding.plot_exp.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_profile_panel.Scanpy_embedding.plot_exp.width,
                  height = input$GEX_profile_panel.Scanpy_embedding.plot_exp.height,
                  units = 'in', res = 300)
              plot(GEX_profile_panel.Scanpy_embedding.plot_exp.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      # Violin plot of gene expression
      if(T){
        GEX_profile_panel.Violin_plot.plot = reactiveValues(plot = NULL)
        output$GEX_profile_panel.Violin_plot <- renderPlot({
          # plot setting
          mycolor <- colorRampPalette(brewer.pal(8, input$GEX_profile_panel.Colorpanel))(100)
          
          # prepare data
          annotation = GEX_profile_panel.scanpy_embedding()$umap[,3]
          expression = GEX_profile_panel.gene_data()$avg_exp
          exp_in_subtype = data.frame(cbind(annotation, expression))
          colnames(exp_in_subtype) = c("label","expression")
          exp_in_subtype$expression = as.numeric(exp_in_subtype$expression)
          
          gene_mean = exp_in_subtype %>% group_by(label) %>% summarise(avg_exp = mean(expression))
          rnames = rownames(exp_in_subtype)
          exp_in_subtype = exp_in_subtype %>% left_join(gene_mean, by = "label")
          rownames(exp_in_subtype) = rnames
          
          # violin plot
          p = ggplot( exp_in_subtype, aes(x = label, y = expression)) + 
            geom_violin(aes(fill = avg_exp), scale = "width", color = "white", kernel = "gaussian") + 
            #geom_boxplot(outlier.size = -1, width = .1, fill = "white") + 
            scale_fill_gradientn("Exp", colors = mycolor) +
            ylab("Relative gene expression") + 
            xlab("")+
            theme_bw() + 
            theme(axis.text.x = element_text(size = 10, family="Times", color = "black", angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right",
                  axis.line = element_line(color = "black"))
          GEX_profile_panel.Violin_plot.plot$plot = p
          print(p)
        })
        
        # download
        output$GEX_profile_panel.Violin_plot.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_profile_panel.Violin_plot.file_type == 'pdf'){
              paste0('GEX_profile_panel.Violin_plot_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_profile_panel.Violin_plot.file_type == 'jpeg'){
              paste0('GEX_profile_panel.Violin_plot_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_profile_panel.Violin_plot.file_type == 'pdf'){
              pdf(file, width = input$GEX_profile_panel.Violin_plot.width, height = input$GEX_profile_panel.Violin_plot.height)
              plot(GEX_profile_panel.Violin_plot.plot$plot)
              dev.off()
            }else if(input$GEX_profile_panel.Scanpy_embedding.plot_exp.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_profile_panel.Violin_plot.width,
                   height = input$GEX_profile_panel.Violin_plot.height,
                   units = 'in', res = 300)
              plot(GEX_profile_panel.Violin_plot.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      # Annotation of cell subsets
      if(T){
        GEX_profile_panel.Scanpy_embedding.plot_label.plot = reactiveValues(plot = NULL)
        output$GEX_profile_panel.Scanpy_embedding.plot_label <- renderPlot({
          # plot setting
          dot.size = input$GEX_profile_panel.Embedding_dotsize
          
          # plot cell subset annotation in Scanpy embedding
          if(input$GEX_profile_panel.Embedding_used == "UMAP"){
            
            # prepare data
            umap_data = GEX_profile_panel.scanpy_embedding()$umap
            colnames(umap_data) = c("UMAP_1","UMAP_2","label")
            
            # set color
            mycolor = minor_cluster_color[ umap_data$label %>% droplevels %>% levels ]
            
            # plot umap
            p = ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = label)) +
              geom_point(size = dot.size) +
              scale_color_manual("Label", values = mycolor) +
              guides(color = guide_legend(override.aes = list(size = 3))) + 
              xlab("UMAP_1") +
              ylab("UMAP_2") + 
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_blank(),
                    legend.position = "bottom") + 
              theme(axis.line = element_line(color = "black"))
            
            GEX_profile_panel.Scanpy_embedding.plot_label.plot$plot = p
            print(p)

          }else if(input$GEX_profile_panel.Embedding_used == "tSNE"){
            
            # prepare data
            tsne_data = GEX_profile_panel.scanpy_embedding()$tsne
            colnames(tsne_data) = c("tSNE_1","tSNE_2","label")
            
            # set color
            mycolor = minor_cluster_color[tsne_data$label %>% droplevels %>% levels] 
            
            # plot tsne
            p = ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2, color = label)) +
              geom_point(size = dot.size) + 
              scale_color_manual("Label", values = mycolor) +
              guides(color = guide_legend(override.aes = list(size = 3))) + 
              xlab("tSNE_1") +
              ylab("tSNE_2") + 
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='transparent', color='white'),
                    panel.border = element_blank(),
                    legend.position = "bottom",
                    axis.line = element_line(color = "black"))
            
            GEX_profile_panel.Scanpy_embedding.plot_label.plot$plot = p
            print(p)
          }
        })
        
        # download
        output$GEX_profile_panel.Scanpy_embedding.plot_label.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_profile_panel.Scanpy_embedding.plot_label.file_type == 'pdf'){
              paste0('GEX_profile_panel.Scanpy_embedding.plot_label_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_profile_panel.Scanpy_embedding.plot_label.file_type == 'jpeg'){
              paste0('GEX_profile_panel.Scanpy_embedding.plot_label_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_profile_panel.Scanpy_embedding.plot_label.file_type == 'pdf'){
              pdf(file, width = input$GEX_profile_panel.Scanpy_embedding.plot_label.width, height = input$GEX_profile_panel.Scanpy_embedding.plot_label.height)
              plot(GEX_profile_panel.Scanpy_embedding.plot_label.plot$plot)
              dev.off()
            }else if(input$GEX_profile_panel.Scanpy_embedding.plot_exp.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_profile_panel.Scanpy_embedding.plot_label.width,
                   height = input$GEX_profile_panel.Scanpy_embedding.plot_label.height,
                   units = 'in', res = 300)
              plot(GEX_profile_panel.Scanpy_embedding.plot_label.plot$plot)
              dev.off()
            }
          }
        )
      }
    
      # Dot plot of gene expression
      if(T){
        GEX_profile_panel.Heatmap_plot.plot = reactiveValues(plot = NULL)
        output$GEX_profile_panel.Heatmap_plot <- renderPlot({
          # set color
          mycolor <- colorRampPalette(brewer.pal(8, input$GEX_profile_panel.Colorpanel))(100)
          
          # prepare data
          annotation = droplevels(GEX_profile_panel.scanpy_embedding()$umap[,3])
          expression = GEX_profile_panel.gene_data()$expression
          
          exp_in_subtype = data.frame(cbind(annotation, expression))
          colnames(exp_in_subtype) = c("label",colnames(expression))
          exp_in_subtype$label = factor(exp_in_subtype$label, levels = levels(annotation))
          
          exp_in_subtype = melt(data = data.table(exp_in_subtype), id.vars=c("label"))
          colnames(exp_in_subtype) = c("label","gene","expression")
          exp_in_subtype$expression = as.numeric(exp_in_subtype$expression)
          
          # summarize mean expression
          group_exp = exp_in_subtype %>%
            group_by(gene, label) %>%
            summarise(avg_exp = mean(expression))
          
          # summarize fraction
          group_frac = exp_in_subtype %>%
            group_by(gene, label) %>%
            summarise(fraction = sum(expression > 0.1)/length(expression))
          
          exp_frac = merge(group_exp, group_frac, by = c("gene","label"))
          colnames(exp_frac) = c("gene","label","Exp","Fraction")
          
            # heatmap plot
          p = ggplot(exp_frac, aes(x = label, y = gene)) + 
              geom_point( aes(size = Fraction, color = Exp) ) + 
              scale_color_gradientn("Exp", colors = mycolor) +
              xlab("") +
              ylab("") + 
              theme_bw() + 
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black", angle = 90, vjust = 0.5, hjust = 1),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='transparent', color='white'),
                    legend.position = "bottom")
          GEX_profile_panel.Heatmap_plot.plot$plot = p
          print(p)
        })
        
        # download
        output$GEX_profile_panel.Heatmap_plot.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_profile_panel.Heatmap_plot.file_type == 'pdf'){
              paste0('GEX_profile_panel.Heatmap_plot_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_profile_panel.Heatmap_plot.file_type == 'jpeg'){
              paste0('GEX_profile_panel.Heatmap_plot_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_profile_panel.Heatmap_plot.file_type == 'pdf'){
              pdf(file, width = input$GEX_profile_panel.Heatmap_plot.width, height = input$GEX_profile_panel.Heatmap_plot.height)
              plot(GEX_profile_panel.Heatmap_plot.plot$plot)
              dev.off()
            }else if(input$GEX_profile_panel.Heatmap_plot.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_profile_panel.Heatmap_plot.width,
                   height = input$GEX_profile_panel.Heatmap_plot.height,
                   units = 'in', res = 300)
              plot(GEX_profile_panel.Heatmap_plot.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      # get subset number
      if(T){
        GEX_profile_panel.subset_cell_number = reactive({
          
          # count
          subset_count = adata_obs %>% 
            filter( major_cluster == input$GEX_profile_panel.major_cluster) %>% 
            select(minor_cluster) %>% 
            droplevels %>% 
            group_by(minor_cluster) %>% 
            summarise(number = n())
         
          # set colnames
          colnames(subset_count) = c("subset","number")
          
          # add color
          subset_count$mycolor = minor_cluster_color[subset_count$subset %>% levels]
          
          # return
          list(subset_count = subset_count %>% data.frame)
        })
      }
      
      # Barplot of cell numbers
      if(T){
        GEX_profile_panel.subset_cell_number.barplot.plot = reactiveValues(plot = NULL)
        output$GEX_profile_panel.subset_cell_number.barplot <- renderPlot({
          
          # get count table
          subset_cell_number = GEX_profile_panel.subset_cell_number()$subset_count
          
          # plot
          p = ggplot(data = subset_cell_number,
                 aes(x = subset, y = number, fill = subset)) +
            geom_bar(stat="identity")+
            scale_fill_manual(values = subset_cell_number$mycolor) + 
            xlab("") + 
            ylab("Number of cells") +
            theme_bw() + 
            theme(axis.text.x = element_text(size=10, family="Times", angle = 90, vjust = 0.5, hjust = 1, color = "black"),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "none") + 
            theme(axis.line = element_line(color = 'black'))
          
          GEX_profile_panel.subset_cell_number.barplot.plot$plot = p
          print(p)
        })
        
        # download
        output$GEX_profile_panel.subset_cell_number.barplot.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_profile_panel.subset_cell_number.barplot.file_type == 'pdf'){
              paste0('GEX_profile_panel.subset_cell_number.barplot_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_profile_panel.subset_cell_number.barplot.file_type == 'jpeg'){
              paste0('GEX_profile_panel.subset_cell_number.barplot_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_profile_panel.subset_cell_number.barplot.file_type == 'pdf'){
              pdf(file, width = input$GEX_profile_panel.subset_cell_number.barplot.width, height = input$GEX_profile_panel.subset_cell_number.barplot.height)
              plot(GEX_profile_panel.subset_cell_number.barplot.plot$plot)
              dev.off()
            }else if(input$GEX_profile_panel.subset_cell_number.barplot.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_profile_panel.subset_cell_number.barplot.width,
                   height = input$GEX_profile_panel.subset_cell_number.barplot.height,
                   units = 'in', res = 300)
              plot(GEX_profile_panel.subset_cell_number.barplot.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      # get top marker genes
      if(T){
        GEX_profile_panel.top_marker_gene_expression = reactive({
          # get top marker genes
          marker_gene = deg_data[[input$GEX_profile_panel.major_cluster]] %>% 
            group_by(cluster) %>% 
            top_n(n = as.numeric(input$GEX_profile_panel.deg_topn), wt = avg_log2FC)
          marker_gene = marker_gene$gene
          
          # get expression
          expression = get_gex(adata, cells = GEX_profile_panel.data()$cells, 
                               genes = marker_gene)
          
          # return
          list(expression = expression)
        })
      }

      # Heatmap plot of marker genes
      if(T){
        output$GEX_profile_panel.marker_gene.heatmap <- renderPlot({
          
          # set color
          mycolor <- colorRampPalette(brewer.pal(8, input$GEX_profile_panel.Colorpanel))(100)
 
          # get expression
          expression = GEX_profile_panel.top_marker_gene_expression()$expression
          
          # scale across cells
          expression = sapply(expression, function(x) (x-mean(x))/sd(x))

          # add cell type to expression matrix
          annotation = droplevels(GEX_profile_panel.scanpy_embedding()$umap[,3])
          exp_in_subtype = data.frame(cbind(as.character(annotation), expression))
          colnames(exp_in_subtype) = c("label",colnames(expression))
          exp_in_subtype$label = factor( exp_in_subtype$label, levels = levels(annotation) )
          
          exp_in_subtype = melt(data = data.table(exp_in_subtype), id.vars=c("label"))
          colnames(exp_in_subtype) = c("label","gene","expression")
          exp_in_subtype$expression = as.numeric(exp_in_subtype$expression)
      
          # summarize data
          group_mean = exp_in_subtype %>%
            group_by(gene, label) %>%
            summarise(avg_exp = mean(expression))
          
          # create heatmap using blue color scale
          ggplot(group_mean, aes(label, gene)) +
            geom_tile(aes(fill = avg_exp), colour = "white") +
            scale_fill_gradientn("Exp", colors = mycolor) + 
            scale_x_discrete(limits = rev(levels(annotation))) + 
            xlab("")+
            ylab("")+
            coord_flip() + 
            theme_bw() + 
            theme(axis.text.x = element_text(size = 10, family="Times", color = "black", angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  legend.position = "right")
        })
      }
      
      # display deg table
      # Marker genes of each cell subtypes
      if(T){
        output$GEX_profile_panel.marker_gene_tbl = DT::renderDataTable(
                            datatable(deg_data[[input$GEX_profile_panel.major_cluster]],
                                      filter = 'top', escape = FALSE, rownames = FALSE) %>% formatRound(3:7, 4))
      }
  }

  # Gene expression comparison
  if(T){
      
      # keep gene list in gene expression panel
      GEX_comparison_panel.genes = reactiveValues(genes = c("TPSAB1","CPA3"))
      
      # observe gene list
      observeEvent(input$GEX_comparison_panel.Submit,{
        ## parse gene list
        genes = str_split(input$GEX_comparison_panel.Gene_list,"\n")[[1]]
        genes = str_trim(genes, side = c("both")) ## trim white space
        genes = genes[genes != ""] ## remove empty elements
        genes = unique(genes) ## remove duplicated elements
        genes = toupper(genes) ## convert to uppercase
        
        # check gene list
        if( length(genes) == 0){
          message = paste0("Please input a gene name!")
          sendSweetAlert(session = session, title = message, type = "warning")
        }else{
          genes_ok = genes[genes %in% all_genes ]
          if( length(genes_ok) == 0 ){
            message = paste0( "Input genes are not found in database!")
            sendSweetAlert(session = session, title = message, type = "warning")
          }else{
            ## store gene list
            GEX_comparison_panel.genes$genes = genes_ok
          }
        }
      })
      
      # observe minor cluster and major cluster
      if(T){ 
        GEX_comparison_panel.group = reactiveValues(
          selected_major_cluster = major_cluster,
          selected_minor_cluster = unlist(all_cluster_list),
          selected_stage = unlist(all_stage_list),
          selected_disease = unlist(all_disease_list),
          selected_tissue = unlist(all_tissue_list),
          selected_tissue.sub = unlist(all_tissue.sub_list),
          selected_study = unlist(all_study_list),
          selected_sample = unlist(all_sample_list)
        )
        
        observeEvent(input$GEX_comparison_panel.Submit,{
          GEX_comparison_panel.group$selected_major_cluster = input$GEX_comparison_panel.major_cluster
          GEX_comparison_panel.group$selected_minor_cluster = input$GEX_comparison_panel.minor_cluster
          GEX_comparison_panel.group$selected_stage = input$GEX_comparison_panel.stage
          GEX_comparison_panel.group$selected_disease = input$GEX_comparison_panel.disease
          GEX_comparison_panel.group$selected_tissue = input$GEX_comparison_panel.tissue
          GEX_comparison_panel.group$selected_tissue.sub = input$GEX_comparison_panel.tissue.sub
          GEX_comparison_panel.group$selected_study = input$GEX_comparison_panel.study
          GEX_comparison_panel.group$selected_sample = input$GEX_comparison_panel.sample
        })
      }

      # get umap and tsne of all cells
      if(T){
        GEX_comparison_panel.scanpy_embedding <- reactive({
          
          # get umap
          umap = adata_obs %>% 
            select(gUMAP_1, gUMAP_2, major_cluster) %>% as.data.frame()
          
          # get tsne
          tsne = adata_obs %>%
            select(gTSNE_1, gTSNE_2, major_cluster) %>% as.data.frame()
          
          # return
          list(umap = umap, tsne = tsne)
        })
      }
      
      # select data
      # get gene expression matrix of given genes and given cells
      # get average expression
      if(T){
        GEX_comparison_panel.subsets_data = reactive({
          # subset data
          selected_genes = GEX_comparison_panel.genes$genes
        
          # subset by major cluster
          selected = adata_obs$major_cluster %in% GEX_comparison_panel.group$selected_major_cluster
          
          # subset by minor cluster
          selected = selected &
            adata_obs$minor_cluster %in% GEX_comparison_panel.group$selected_minor_cluster
          
          # subset by stage
          selected = selected &
            adata_obs$stage %in% GEX_comparison_panel.group$selected_stage
          
          # subset by disease
          selected = selected &
            adata_obs$disease %in% GEX_comparison_panel.group$selected_disease
          
          # subset by tissue
          selected = selected &
            adata_obs$tissue %in% GEX_comparison_panel.group$selected_tissue
          
          # subset by tissue
          selected = selected &
            adata_obs$tissue.sub %in% GEX_comparison_panel.group$selected_tissue.sub
          
          # subset by study
          selected = selected &
            adata_obs$study %in% GEX_comparison_panel.group$selected_study
          
          # subset by sample
          selected = selected &
            adata_obs$sample %in% GEX_comparison_panel.group$selected_sample
          
          # get selected cells
          selected_cells = rownames(adata_obs[selected, ])
          
          # get expression
          expression = get_gex(adata, 
                               cells = selected_cells, 
                               genes = selected_genes)
          
          # get subset meta data
          selected_meta_data = adata_obs[selected_cells, ]
          
          # get average expression of given input genes
          avg_exp = rowMeans(expression)
          avg_exp = as.data.frame(avg_exp)
          
          # update selected_meta_data
          selected_meta_data$major_cluster = droplevels(selected_meta_data$major_cluster)
          selected_meta_data$minor_cluster = droplevels(selected_meta_data$minor_cluster)
          selected_meta_data$stage = droplevels(selected_meta_data$stage)
          selected_meta_data$disease = droplevels(selected_meta_data$disease)
          selected_meta_data$tissue = droplevels(selected_meta_data$tissue)
          selected_meta_data$tissue.sub = droplevels(selected_meta_data$tissue.sub)
          selected_meta_data$study = droplevels(selected_meta_data$study)
          selected_meta_data$sample = droplevels(selected_meta_data$sample)
          
          # return expression
          list(avg_exp = avg_exp, 
               expression = expression,
               selected_meta_data = selected_meta_data)
        })
      }
      
      # get umap and tsne of selected cells
      if(T){
        GEX_comparison_panel.selected_scanpy_embedding <- reactive({
          
          # get umap
          umap = GEX_comparison_panel.subsets_data()$selected_meta_data %>% 
            select(gUMAP_1, gUMAP_2) %>% as.data.frame()
          
          # get tsne
          tsne = GEX_comparison_panel.subsets_data()$selected_meta_data %>%
            select(gTSNE_1, gTSNE_2) %>% as.data.frame()
          
          # return
          list(umap = umap, tsne = tsne)
        })
      }
      
      # get average expression of given genes
      # all cells
      if(T){
        GEX_comparison_panel.global_avg_exp = reactive({
          
          # get average expression of given genes
          expression = get_gex(adata, 
                               cells = adata_obs %>% rownames, 
                               genes = GEX_comparison_panel.genes$genes)
          
          # get average expression
          avg_exp = rowMeans(expression)
          
          # clean
          rm(expression)
          
          # return average expression
          list(avg_exp = avg_exp)
        })
      }
      
      # plot global scanpy embedding (left)
      # all cells
      # color by average expression of given genes
      if(T){
        GEX_comparison_panel.Scanpy_embedding.plot_exp.plot = reactiveValues(plot = NULL)
        output$GEX_comparison_panel.Scanpy_embedding.plot_exp <- renderPlot({
          
          # plot setting
          mycolor <- colorRampPalette(brewer.pal(8, input$GEX_comparison_panel.Colorpanel))(100)
          dot.size = input$GEX_comparison_panel.Embedding_dotsize
          
          # plot gene expression in scanpy embedding
          if(input$GEX_comparison_panel.Embedding_used == "UMAP"){
            
            # prepare data
            # get umap
            umap_data = cbind(GEX_comparison_panel.scanpy_embedding()$umap[, c(1,2)],
                              GEX_comparison_panel.global_avg_exp()$avg_exp )
            colnames(umap_data) = c("UMAP_1","UMAP_2","avg_exp")
            
            # plot umap
            p = ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
              geom_point(aes(color = avg_exp), size = dot.size) + 
              scale_colour_gradientn("Exp", colors = mycolor) +
              xlab("UMAP_1") +
              ylab("UMAP_2") +
              xlim( min(adata_obs$gUMAP_1), max(adata_obs$gUMAP_1) ) + 
              ylim( min(adata_obs$gUMAP_2), max(adata_obs$gUMAP_2) ) + 
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, color = "black"),
                    axis.text.y = element_text(size=10, color = "black"),
                    axis.title.x = element_text(size=12, family="Times",  color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='transparent', color='white'),
                    panel.border = element_blank(),
                    legend.position = "right") + 
              theme(axis.line = element_line(color = "black"))
            
            GEX_comparison_panel.Scanpy_embedding.plot_exp.plot$plot = p
            print(p)
            
          }else if(input$GEX_comparison_panel.Embedding_used == "tSNE"){
            
            # prepare data
            # get umap
            tsne_data = cbind(GEX_comparison_panel.scanpy_embedding()$umap[, c(1,2)],
                              GEX_comparison_panel.global_avg_exp()$avg_exp )
            colnames(tsne_data) = c("tSNE_1","tSNE_2","avg_exp")
            
            # plot tsne
            p = ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2)) + 
              geom_point(aes(color = avg_exp), size = dot.size) + 
              scale_colour_gradientn("Exp", colors = mycolor) + 
              xlab("tSNE_1") + 
              ylab("tSNE_2") + 
              xlim( min(adata_obs$gUMAP_1), max(adata_obs$gUMAP_1) ) + 
              ylim( min(adata_obs$gUMAP_2), max(adata_obs$gUMAP_2) ) + 
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, color = "black"),
                    axis.text.y = element_text(size=10, color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times",  color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='transparent', color='white'),
                    panel.border = element_blank(),
                    legend.position = "right") + 
              theme(axis.line = element_line(color = "black"))
            
            GEX_comparison_panel.Scanpy_embedding.plot_exp.plot$plot = p
            print(p)
          }
        })
        
        # download
        output$GEX_comparison_panel.Scanpy_embedding.plot_exp.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_comparison_panel.Scanpy_embedding.plot_exp.file_type == 'pdf'){
              paste0('GEX_comparison_panel.Scanpy_embedding.plot_exp_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_comparison_panel.Scanpy_embedding.plot_exp.file_type == 'jpeg'){
              paste0('GEX_comparison_panel.Scanpy_embedding.plot_exp_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_comparison_panel.Scanpy_embedding.plot_exp.file_type == 'pdf'){
              pdf(file, width = input$GEX_comparison_panel.Scanpy_embedding.plot_exp.width, height = input$GEX_comparison_panel.Scanpy_embedding.plot_exp.height)
              plot(GEX_comparison_panel.Scanpy_embedding.plot_exp.plot$plot)
              dev.off()
            }else if(input$GEX_comparison_panel.Scanpy_embedding.plot_exp.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_comparison_panel.Scanpy_embedding.plot_exp.width,
                   height = input$GEX_comparison_panel.Scanpy_embedding.plot_exp.height,
                   units = 'in', res = 300)
              plot(GEX_comparison_panel.Scanpy_embedding.plot_exp.plot$plot)
              dev.off()
            }
          }
        )
      }
    
      # plot global scanpy embedding (middle)
      # all cells
      # color by major cell type
      if(T){
        GEX_comparison_panel.Scanpy_embedding.plot_label.plot = reactiveValues(plot = NULL)
        output$GEX_comparison_panel.Scanpy_embedding.plot_label <- renderPlot({
          
          # plot setting
          dot.size = input$GEX_comparison_panel.Embedding_dotsize
          
          # plot cell subset annotation in Scanpy embedding
          if(input$GEX_comparison_panel.Embedding_used == "UMAP"){
            ## prepare data
            umap_data = GEX_comparison_panel.scanpy_embedding()$umap
            colnames(umap_data) = c("UMAP_1","UMAP_2","label")
            
            # get position of labels
            if(F){
              label_pos <- umap_data %>% 
                group_by(label) %>% 
                summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
            }
            
            # plot umap
            p = ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = label)) +
              geom_point(size = dot.size) +
              scale_color_manual("Label", values = major_cluster_color ) + 
              xlab("UMAP_1") +
              ylab("UMAP_2") + 
              xlim( min(adata_obs$gUMAP_1), max(adata_obs$gUMAP_1) ) + 
              ylim( min(adata_obs$gUMAP_2), max(adata_obs$gUMAP_2) ) + 
              guides(color = guide_legend(override.aes = list(size = 3))) + 
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_blank(),
                    legend.position = "bottom") +  # show legend on right
              theme(axis.line = element_line(color = "black"))
            
            GEX_comparison_panel.Scanpy_embedding.plot_label.plot$plot = p
            print(p)
            
          }else if(input$GEX_comparison_panel.Embedding_used == "tSNE"){
            ## prepare data
            tsne_data = GEX_comparison_panel.scanpy_embedding()$tsne
            colnames(tsne_data) = c("tSNE_1","tSNE_2","label")
            
            if(F){ # get position of labels
              label_pos <- tsne_data %>% 
                group_by(label) %>% 
                summarize(tSNE_1 = mean(tSNE_1), tSNE_2 = mean(tSNE_2))
            }
            
            # plot tsne
            ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2, color = label)) +
              geom_point(size = dot.size) + 
              scale_color_manual("Label", values = major_cluster_color ) +
              xlab("tSNE_1") +
              ylab("tSNE_2") + 
              xlim( min(adata_obs$gUMAP_1), max(adata_obs$gUMAP_1) ) + 
              ylim( min(adata_obs$gUMAP_2), max(adata_obs$gUMAP_2) ) + 
              guides(color = guide_legend(override.aes = list(size = 3))) +
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_blank(),
                    legend.position = "bottom") +  # show legend on right
              theme(axis.line = element_line(color = "black"))
            
            GEX_comparison_panel.Scanpy_embedding.plot_label.plot$plot = p
            print(p)
          }
        })
        
        # download
        output$GEX_comparison_panel.Scanpy_embedding.plot_label.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_comparison_panel.Scanpy_embedding.plot_label.file_type == 'pdf'){
              paste0('GEX_comparison_panel.Scanpy_embedding.plot_label_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_comparison_panel.Scanpy_embedding.plot_label.file_type == 'jpeg'){
              paste0('GEX_comparison_panel.Scanpy_embedding.plot_label_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_comparison_panel.Scanpy_embedding.plot_label.file_type == 'pdf'){
              pdf(file, width = input$GEX_comparison_panel.Scanpy_embedding.plot_label.width, height = input$GEX_comparison_panel.Scanpy_embedding.plot_label.height)
              plot(GEX_comparison_panel.Scanpy_embedding.plot_label.plot$plot)
              dev.off()
            }else if(input$GEX_comparison_panel.Scanpy_embedding.plot_label.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_comparison_panel.Scanpy_embedding.plot_label.width,
                   height = input$GEX_comparison_panel.Scanpy_embedding.plot_label.height,
                   units = 'in', res = 300)
              plot(GEX_comparison_panel.Scanpy_embedding.plot_label.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      # plot global scanpy embedding (left)
      # selected cells
      # colored by user defined values (major/minor/stage/dev/disease)
      if(T){
        GEX_comparison_panel.Scanpy_embedding.plot_cell.plot = reactiveValues(plot = NULL)
        output$GEX_comparison_panel.Scanpy_embedding.plot_cell <- renderPlot({
          
          # setting
          dot.size = input$GEX_comparison_panel.Embedding_dotsize
          
          # plot cell in Scanpy embedding
          if(input$GEX_comparison_panel.Embedding_used == "UMAP"){
            
            # prepare data
            # get umap data
            umap_data = GEX_comparison_panel.selected_scanpy_embedding()$umap
            # add label
            umap_data$label = GEX_comparison_panel.subsets_data()$selected_meta_data[rownames(umap_data), input$GEX_comparison_panel.groupby]
            colnames(umap_data) = c("UMAP_1","UMAP_2","label")
            
            # set color
            mycolor = ""
            if(input$GEX_comparison_panel.groupby == "major_cluster"){
              mycolor = major_cluster_color
            }else if(input$GEX_comparison_panel.groupby == "minor_cluster"){
              mycolor = minor_cluster_color[umap_data$label %>% droplevels %>% levels]
            }else{
              colourCount = umap_data$label %>% droplevels %>% levels %>% length
              getPalette = colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_comparison_panel.Colorpanel,'maxcolors'], input$GEX_comparison_panel.Colorpanel))
              mycolor = getPalette(colourCount)
            }
            
            # plot
            p = ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = label)) +
              geom_point(size = dot.size) +
              scale_color_manual("Label", values = mycolor) +
              xlab("UMAP_1") +
              ylab("UMAP_2") + 
              xlim( min(adata_obs$gUMAP_1), max(adata_obs$gUMAP_1) ) + 
              ylim( min(adata_obs$gUMAP_2), max(adata_obs$gUMAP_2) ) + 
              guides(color = guide_legend(override.aes = list(size = 3))) + 
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_blank(),
                    legend.position = "bottom") +  # show legend on right
              theme(axis.line = element_line(color = "black"))
            
            GEX_comparison_panel.Scanpy_embedding.plot_cell.plot$plot = p
            print(p)
            
          }else if(input$GEX_comparison_panel.Embedding_used == "tSNE"){
            
            # prepare data
            # get tsne data
            tsne_data = GEX_comparison_panel.selected_scanpy_embedding()$umap
            # add label
            tsne_data$label = GEX_comparison_panel.subsets_data()$selected_meta_data[rownames(tsne_data), input$GEX_comparison_panel.groupby]
            colnames(tsne_data) = c("tSNE_1","tSNE_2","label")
            
            # set color
            mycolor = ""
            if(input$GEX_comparison_panel.groupby == "major_cluster"){
              mycolor = major_cluster_color
            }else if(input$GEX_comparison_panel.groupby == "minor_cluster"){
              mycolor = minor_cluster_color[umap_data$label %>% droplevels %>% levels]
            }else{
              colourCount = umap_data$label %>% droplevels %>% levels %>% length
              getPalette = colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_comparison_panel.Colorpanel,'maxcolors'], input$GEX_comparison_panel.Colorpanel))
              mycolor = getPalette(colourCount)
            }
            
            # plot tsne
            p = ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2, color = label)) +
              geom_point(size = dot.size) + 
              scale_color_manual("Label", values = mycolor) +
              xlab("tSNE_1") +
              ylab("tSNE_2") + 
              xlim( min(adata_obs$gUMAP_1), max(adata_obs$gUMAP_1) ) + 
              ylim( min(adata_obs$gUMAP_2), max(adata_obs$gUMAP_2) ) + 
              guides(color = guide_legend(override.aes = list(size = 3))) +
              theme_bw() +
              theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                    axis.text.y = element_text(size=10, family="Times", color = "black"),
                    axis.title.x = element_text(size=12, family="Times", color = "black"),
                    axis.title.y = element_text(size=12, family="Times", color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_blank(),
                    legend.position = "bottom") +  # show legend on right
              theme(axis.line = element_line(color = "black"))
            
            GEX_comparison_panel.Scanpy_embedding.plot_cell.plot$plot = p
            print(p)
          } 
        })
        
        # download
        output$GEX_comparison_panel.Scanpy_embedding.plot_cell.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_comparison_panel.Scanpy_embedding.plot_cell.file_type == 'pdf'){
              paste0('GEX_comparison_panel.Scanpy_embedding.plot_cell_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_comparison_panel.Scanpy_embedding.plot_cell.file_type == 'jpeg'){
              paste0('GEX_comparison_panel.Scanpy_embedding.plot_cell_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_comparison_panel.Scanpy_embedding.plot_cell.file_type == 'pdf'){
              pdf(file, width = input$GEX_comparison_panel.Scanpy_embedding.plot_cell.width, height = input$GEX_comparison_panel.Scanpy_embedding.plot_cell.height)
              plot(GEX_comparison_panel.Scanpy_embedding.plot_cell.plot$plot)
              dev.off()
            }else if(input$GEX_comparison_panel.Scanpy_embedding.plot_cell.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_comparison_panel.Scanpy_embedding.plot_cell.width,
                   height = input$GEX_comparison_panel.Scanpy_embedding.plot_cell.height,
                   units = 'in', res = 300)
              plot(GEX_comparison_panel.Scanpy_embedding.plot_cell.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      # plot gene expression of selected cells
      # violin plot
      # group by major/minor/stage/disease/tissue
      if(T){
        GEX_comparison_panel.Violin_plot.plot = reactiveValues(plot = NULL)
        output$GEX_comparison_panel.Violin_plot <- renderPlot({
          
          # set color
          mycolor = colorRampPalette(brewer.pal(8, input$GEX_comparison_panel.Colorpanel))(100)
          
          # prepare data
          # get label
          label = GEX_comparison_panel.subsets_data()$selected_meta_data[, input$GEX_comparison_panel.groupby]
          
          # prepare data
          exp_in_subtype = data.frame(cbind(label, GEX_comparison_panel.subsets_data()$avg_exp))
          colnames(exp_in_subtype) = c("group","expression")
          
          # get average expression by group
          group_avg_exp = exp_in_subtype %>% group_by(group) %>% summarise(avg_exp = mean(expression))
          rnames = rownames(exp_in_subtype)
          
          exp_in_subtype = exp_in_subtype %>% left_join(group_avg_exp, by = c("group"))
          rownames(exp_in_subtype) = rnames
          
          # violin plot
          p = ggplot( exp_in_subtype, aes(x = group, y = expression)) + 
            geom_violin(aes(fill = avg_exp), scale = "width", color = "white", kernel = "gaussian") + 
            # geom_boxplot(outlier.size = -1, width = .1, fill = "white") + 
            scale_fill_gradientn("Exp", colors = mycolor) +
            ylab("Relative gene expression") + 
            xlab("")+
            theme_bw() + 
            theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size=10,  color = "black"),
                  axis.title.x = element_text(size=12, color = "black"),
                  axis.title.y = element_text(size=12,  color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") +
            theme(axis.line = element_line(color = "black"))
          
          GEX_comparison_panel.Violin_plot.plot$plot = p
          print(p)
        })
        
        # download
        output$GEX_comparison_panel.Violin_plot.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_comparison_panel.Violin_plot.file_type == 'pdf'){
              paste0('GEX_comparison_panel.Violin_plot_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_comparison_panel.Violin_plot.file_type == 'jpeg'){
              paste0('GEX_comparison_panel.Violin_plot_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_comparison_panel.Violin_plot.file_type == 'pdf'){
              pdf(file, width = input$GEX_comparison_panel.Violin_plot.width, height = input$GEX_comparison_panel.Violin_plot.height)
              plot(GEX_comparison_panel.Violin_plot.plot$plot)
              dev.off()
            }else if(input$GEX_comparison_panel.Violin_plot.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_comparison_panel.Violin_plot.width,
                   height = input$GEX_comparison_panel.Violin_plot.height,
                   units = 'in', res = 300)
              plot(GEX_comparison_panel.Violin_plot.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      # plot cell number of selected cells
      # bar plot
      # group by one variable 
      # major/minor/stage/disease/tissue
      if(T){
        GEX_comparison_panel.Cell_number.Bar_plot.plot = reactiveValues(plot = NULL)
        output$GEX_comparison_panel.Cell_number.Bar_plot <- renderPlot({

          # prepare data
          selected_cells = GEX_comparison_panel.subsets_data()$selected_meta_data %>% select(input$GEX_comparison_panel.groupby)
          colnames(selected_cells) = c("group")
          
          # count, group by group
          cell_number_count = selected_cells %>% group_by(group) %>% summarise(count = n())
          
          # calculate percentage of cells
          cell_number_pct = cell_number_count
          cell_number_pct$pct = cell_number_count$count/sum(cell_number_count$count)
          
          if(input$GEX_comparison_panel.Cell_number.Bar_plot.type == 'Number'){
            # barplot
            p = ggplot( cell_number_count, aes(x = group, y = count, fill = group)) + 
              geom_bar(stat = "identity") + 
              ylab("Number of cells") + 
              xlab("")+
              theme_bw() +
              theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                    axis.text.y = element_text(size=10,  color = "black"),
                    axis.title.x = element_text(size=12, color = "black"),
                    axis.title.y = element_text(size=12,  color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='transparent', color='white'),
                    panel.border = element_blank(),
                    legend.position = "right") +
              theme(axis.line = element_line(color = "black"),
                    legend.position = "none")
          }else if(input$GEX_comparison_panel.Cell_number.Bar_plot.type == 'Percentage'){
            # barplot
            p = ggplot( cell_number_pct, aes(x = group, y = pct, fill = group)) + 
              geom_bar(stat = "identity") +
              ylab("Percentage of cells") + 
              xlab("")+
              theme_bw() +
              theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                    axis.text.y = element_text(size=10,  color = "black"),
                    axis.title.x = element_text(size=12, color = "black"),
                    axis.title.y = element_text(size=12,  color = "black"),     
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill='transparent', color='white'),
                    panel.border = element_blank(),
                    legend.position = "right") +
              theme(axis.line = element_line(color = "black"),
                    legend.position = "none")
          }
            
          GEX_comparison_panel.Cell_number.Bar_plot.plot$plot = p
          print(p)
        })
        
        # download
        output$GEX_comparison_panel.Cell_number.Bar_plot.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_comparison_panel.Cell_number.Bar_plot.file_type == 'pdf'){
              paste0('GEX_comparison_panel.Cell_number.Bar_plot_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_comparison_panel.Cell_number.Bar_plot.file_type == 'jpeg'){
              paste0('GEX_comparison_panel.Cell_number.Bar_plot_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_comparison_panel.Cell_number.Bar_plot.file_type == 'pdf'){
              pdf(file, width = input$GEX_comparison_panel.Cell_number.Bar_plot.width, height = input$GEX_comparison_panel.Cell_number.Bar_plot.height)
              plot(GEX_comparison_panel.Cell_number.Bar_plot.plot$plot)
              dev.off()
            }else if(input$GEX_comparison_panel.Cell_number.Bar_plot.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_comparison_panel.Cell_number.Bar_plot.width,
                   height = input$GEX_comparison_panel.Cell_number.Bar_plot.height,
                   units = 'in', res = 300)
              plot(GEX_comparison_panel.Cell_number.Bar_plot.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      # plot gene expression of selected cells
      # violin plot
      # group by two variable, major/minor, stage/disease/tissue
      if(T){
        GEX_comparison_panel.Violin_plot2.plot = reactiveValues(plot = NULL)
        output$GEX_comparison_panel.Violin_plot2 <- renderPlot({
          
          # set color
          mycolor <- colorRampPalette(brewer.pal(8, input$GEX_comparison_panel.Colorpanel))(100)
          
          # prepare data
          # get label
          exp_in_subtype = data.frame( cbind(GEX_comparison_panel.subsets_data()$selected_meta_data[,c(input$GEX_comparison_panel.groupby.1, input$GEX_comparison_panel.groupby.2)],
                                       GEX_comparison_panel.subsets_data()$avg_exp))
          colnames(exp_in_subtype) = c("subset","group","expression")
          exp_in_subtype$expression = as.numeric(exp_in_subtype$expression)
          
          gene_mean = exp_in_subtype %>% group_by(subset, group) %>% summarise(avg_exp = mean(expression))
          rnames = rownames(exp_in_subtype)
          
          exp_in_subtype = exp_in_subtype %>% left_join(gene_mean, by = c("subset","group"))
          rownames(exp_in_subtype) = rnames
          
          # violin plot
          p = ggplot( exp_in_subtype, aes(x = group, y = expression)) + 
            geom_violin(aes(fill = avg_exp), scale = "width", color = "white", kernel = "gaussian") + 
            #geom_boxplot(outlier.size = -1, width = .1, fill = "white") + 
            scale_fill_gradientn("Exp", colors = mycolor) +
            ylab("Relative gene expression") + 
            xlab("")+
            facet_grid(cols = vars(subset)) + 
            theme_bw() + 
            theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size=10,  color = "black"),
                  axis.title.x = element_text(size=12, color = "black"),
                  axis.title.y = element_text(size=12,  color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") +
            theme(axis.line = element_line(color = "black"))
          
          GEX_comparison_panel.Violin_plot2.plot$plot = p
          print(p)
        })
        
        # download
        output$GEX_comparison_panel.Violin_plot2.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_comparison_panel.Violin_plot2.file_type == 'pdf'){
              paste0('GEX_comparison_panel.Violin_plot2_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_comparison_panel.Violin_plot2.file_type == 'jpeg'){
              paste0('GEX_comparison_panel.Violin_plot2_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_comparison_panel.Violin_plot2.file_type == 'pdf'){
              pdf(file, width = input$GEX_comparison_panel.Violin_plot2.width, height = input$GEX_comparison_panel.Violin_plot2.height)
              plot(GEX_comparison_panel.Violin_plot2.plot$plot)
              dev.off()
            }else if(input$GEX_comparison_panel.Violin_plot2.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_comparison_panel.Violin_plot2.width,
                   height = input$GEX_comparison_panel.Violin_plot2.height,
                   units = 'in', res = 300)
              plot(GEX_comparison_panel.Violin_plot2.plot$plot)
              dev.off()
            }
          }
        )
      }
      
      observe({
        # set faceting_group
        updateSelectInput(
          session = session,
          inputId = "GEX_comparison_panel.Bar_plot2.faceting_group",
          choices = c(input$GEX_comparison_panel.groupby.1, input$GEX_comparison_panel.groupby.2),
          selected = input$GEX_comparison_panel.groupby.1,
      )})
        
      # plot cell number of selected cells
      # bar plot
      # group by two variables
      if(T){
        GEX_comparison_panel.Cell_number.Bar_plot2.plot = reactiveValues(plot = NULL)
        output$GEX_comparison_panel.Cell_number.Bar_plot2 <- renderPlot({
          
          # prepare data
          selected_cells = GEX_comparison_panel.subsets_data()$selected_meta_data[,c(input$GEX_comparison_panel.groupby.1, input$GEX_comparison_panel.groupby.2)]
          colnames(selected_cells) = c("group1","group2")
          
          # count, group by subset and group
          cell_number_count = selected_cells %>% 
            group_by(group1, group2) %>% 
            summarise(count = n())
          
          # count, group by subset
          cell_number_total = selected_cells %>% 
            group_by(group1) %>%
            summarise(count = n())
          
          # calculate percentage of cells
          cell_number_pct = merge(cell_number_count, cell_number_total, by = "group1")
          cell_number_pct$pct = cell_number_pct$count.x / cell_number_pct$count.y
          
          # barplot
          if(input$GEX_comparison_panel.Bar_plot2.type == 'Number'){
            if(input$GEX_comparison_panel.Bar_plot2.faceting_group == input$GEX_comparison_panel.groupby.1){
              p = ggplot( cell_number_count, aes( x = group2, y = count, fill = group2)) + 
                facet_grid(cols = vars(group1)) +
                geom_bar(stat = "identity") + 
                ylab("Number of cells") + 
                xlab("")+
                theme_bw() + 
                theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                      axis.text.y = element_text(size=10,  color = "black"),
                      axis.title.x = element_text(size=12, color = "black"),
                      axis.title.y = element_text(size=12,  color = "black"),     
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill='transparent', color='white'),
                      panel.border = element_blank(),
                      legend.position = "right") +
                theme(axis.line = element_line(color = "black"),
                      legend.position = "none")
            }
            else if(input$GEX_comparison_panel.Bar_plot2.faceting_group == input$GEX_comparison_panel.groupby.2){
              p = ggplot( cell_number_count, aes( x = group1, y = count, fill = group1)) + 
                facet_grid(cols = vars(group2)) +
                geom_bar(stat = "identity") + 
                ylab("Number of cells") + 
                xlab("")+
                theme_bw() + 
                theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                      axis.text.y = element_text(size=10,  color = "black"),
                      axis.title.x = element_text(size=12, color = "black"),
                      axis.title.y = element_text(size=12,  color = "black"),     
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill='transparent', color='white'),
                      panel.border = element_blank(),
                      legend.position = "right") +
                theme(axis.line = element_line(color = "black"),
                      legend.position = "none")
            }
          }else if(input$GEX_comparison_panel.Bar_plot2.type == 'Percentage'){
            if(input$GEX_comparison_panel.Bar_plot2.faceting_group == input$GEX_comparison_panel.groupby.1){
              p = ggplot( cell_number_pct, aes( x = group2, y = pct, fill = group2)) + 
                facet_grid(cols = vars(group1)) +
                geom_bar(stat = "identity") + 
                ylab("Percentage of cells") + 
                xlab("")+
                theme_bw() + 
                theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                      axis.text.y = element_text(size=10,  color = "black"),
                      axis.title.x = element_text(size=12, color = "black"),
                      axis.title.y = element_text(size=12,  color = "black"),     
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill='transparent', color='white'),
                      panel.border = element_blank(),
                      legend.position = "right") +
                theme(axis.line = element_line(color = "black"),
                      legend.position = "none")
            }else if(input$GEX_comparison_panel.Bar_plot2.faceting_group == input$GEX_comparison_panel.groupby.2){
              p = ggplot( cell_number_pct, aes( x = group1, y = pct, fill = group1)) + 
                facet_grid(cols = vars(group2)) +
                geom_bar(stat = "identity") + 
                ylab("Percentage of cells") + 
                xlab("")+
                theme_bw() + 
                theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                      axis.text.y = element_text(size=10,  color = "black"),
                      axis.title.x = element_text(size=12, color = "black"),
                      axis.title.y = element_text(size=12,  color = "black"),     
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill='transparent', color='white'),
                      panel.border = element_blank(),
                      legend.position = "right") +
                theme(axis.line = element_line(color = "black"),
                      legend.position = "none")
            }
          }
          
          GEX_comparison_panel.Cell_number.Bar_plot2.plot$plot = p
          print(p)
        })
        
        # download 
        output$GEX_comparison_panel.Cell_number.Bar_plot2.download_pdf <- downloadHandler(
          filename = function(){
            if(input$GEX_comparison_panel.Cell_number.Bar_plot2.file_type == 'pdf'){
              paste0('GEX_comparison_panel.Cell_number.Bar_plot2_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$GEX_comparison_panel.Cell_number.Bar_plot2.file_type == 'jpeg'){
              paste0('GEX_comparison_panel.Cell_number.Bar_plot2_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$GEX_comparison_panel.Cell_number.Bar_plot2.file_type == 'pdf'){
              pdf(file, width = input$GEX_comparison_panel.Cell_number.Bar_plot2.width, height = input$GEX_comparison_panel.Cell_number.Bar_plot2.height)
              plot(GEX_comparison_panel.Cell_number.Bar_plot2.plot$plot)
              dev.off()
            }else if(input$GEX_comparison_panel.Cell_number.Bar_plot2.file_type == 'jpeg'){
              jpeg(file, width = input$GEX_comparison_panel.Cell_number.Bar_plot2.width,
                   height = input$GEX_comparison_panel.Cell_number.Bar_plot2.height,
                   units = 'in', res = 300)
              plot(GEX_comparison_panel.Cell_number.Bar_plot2.plot$plot)
              dev.off()
            }
          }
        )
      }
  }
  
  # Regulon activity panel
  if(T){
    
    # observe major cluster
    observe({
      updatePickerInput(
        session = session,
        inputId = "GRN_profile_panel.TF",
        choices = tf_list[[ input$GRN_profile_panel.major_cluster ]],
        selected = tf_list[[ input$GRN_profile_panel.major_cluster ]][1:3],
        choicesOpt = list(style = "background:white; color:black")
      )
    })
      
    # get loom file
    if(T){
      # select loom
      GRN_profile_panel.data = reactive({
        # return result
        list(scenic_loom = loom_data[[ input$GRN_profile_panel.major_cluster ]])
      })
    }
    
    # get regulon activity 
    if(T){ 
      GRN_profile_panel.scenic_activity = reactive({
        
        # get auc matrix of regulons
        regulonsAUC = get_regulons_AUC( GRN_profile_panel.data()$scenic_loom, 
                                        column.attr.name = "RegulonsAUC")
        
        # return
        list(grn_activity = getAUC(regulonsAUC) %>% t %>% as.data.frame )
      })
    }
    
    # get genes
    if(T){ 
      # use tfs shared among major clusters
      GRN_profile_panel.genes = reactiveValues(genes = tf_list[["Myeloid"]][1:3])
    }
    
    # observe submit
    if(T){
      observeEvent(input$GRN_profile_panel.Submit, {
        GRN_profile_panel.genes$genes = input$GRN_profile_panel.TF
        })
    }
    
    # get cell information
    if(T) { 
      GRN_profile_panel.cellInfo = reactive({
        # return
        list(cellInfo = get_cell_annotation(GRN_profile_panel.data()$scenic_loom))
      })
    }
    
    # get scanpy embedding
    if(T){
      GRN_profile_panel.scanpy_embedding <- reactive({
        # return umap and tsne
        list(umap = get_embeddings(GRN_profile_panel.data()$scenic_loom)$Scanpy_UMAP, 
             tsne = get_embeddings(GRN_profile_panel.data()$scenic_loom)$Scanpy_tSNE)
      })
    }
    
    # get scenic embedding
    if(T){
      GRN_profile_panel.scenic_embedding <- reactive({
        # return umap and tsne
        list(umap = get_embeddings(GRN_profile_panel.data()$scenic_loom)$SCENIC_AUC_UMAP,
             tsne = get_embeddings(GRN_profile_panel.data()$scenic_loom)$SCENIC_AUC_tSNE)
      })
    }

    # get grn network
    if(T){
      GRN_profile_panel.scenic_network <- reactive({
        
        # get network
        regulons_incidMat <- get_regulons(GRN_profile_panel.data()$scenic_loom, column.attr.name = "Regulons") 
        network <- regulonsToGeneLists(regulons_incidMat)
        
        # get network of given regulons
        regulon_name = paste0(GRN_profile_panel.genes$genes, "(+)")
        grn_network = network[regulon_name]
        names(grn_network) = GRN_profile_panel.genes$genes
        
        # get correlation between TF and downstream genes
        tf_gene_cor_table = data.frame(TF = NULL, target = NULL, correlation = NULL, p.value = NULL, expression = NULL)
        for(i in 1:length(grn_network)){
          tf = names(grn_network)[i]
          targets = grn_network[[ tf ]]
          gex = get_gex(adata, 
                  cells = adata_obs[adata_obs$major_cluster == input$GRN_profile_panel.major_cluster, ] %>% rownames, 
                  genes = c(tf, targets))

          # the first is tf, the rest is targets
          for(i in 2:length(gex)){
            tmp = gex[,c(1,i)]
            tmp = tmp[tmp[,1] > 0 & tmp[,2] > 0,]
            tmp_avg_exp = mean(tmp[,2]) %>% round(digits = 4)
            spearman_cor = cor(tmp[,1], tmp[,2], method = "spearman") %>% round(digits = 4)
            p.value = cor.test(tmp[,1], tmp[,2], method = "spearman", exact = FALSE)$p.value %>% 
              format(digits = 4, scientific = TRUE)
            tf_gene_cor_table = rbind(tf_gene_cor_table, c(tf, colnames(gex)[i], spearman_cor, p.value, tmp_avg_exp))
            colnames(tf_gene_cor_table) = c("TF","target","correlation","p-value","expression")
          }
        }
        # tf_gene_cor_table$TF = factor(tf_gene_cor_table$TF, levels = tf_gene_cor_table$TF %>% unique %>% sort)
        # clean
        rm(tf, targets, gex)
        
        # get nodes and links of the GRN networks
        src = c()
        target = c()
        for(i in 1:length(grn_network)){
          this_target = grn_network[[i]]
          this_src = names(grn_network)[i]
          src = c(src, rep(this_src, length(this_target)))
          target = c(target, this_target)
        }
        networkData <- data.frame(src, target)
        
        # record source nodes, set group =1, size = 30
        src_nodes = unique(networkData$src)
        scenic_src_nodes = data.frame( name = src_nodes,
                                       group = 1,
                                       size = 30,
                                       id = 1:length(src_nodes) - 1)
        # record target nodes, set group = 2, size = 20
        target_nodes = unique(networkData$target)
        target_nodes = target_nodes[ !(target_nodes %in% src_nodes)]
        scenic_target_nodes = data.frame( name = target_nodes,
                                          group = 2,
                                          size = 20,
                                          id = length(src_nodes) + 1:length(target_nodes) - 1)
        # combine source nodes, and target nodes
        scenic_nodes = rbind(scenic_src_nodes, scenic_target_nodes)
        
        # create links
        colnames(scenic_nodes)[1] = "src"
        tmp = merge(networkData, scenic_nodes, by = "src")
        colnames(scenic_nodes)[1] = "target"
        scenic_links = merge(tmp, scenic_nodes, by = "target")
        scenic_links = scenic_links[,c("id.x","id.y")]
        colnames(scenic_links) = c("source","target")
        scenic_links$value = 10
        names(scenic_nodes)[1] = "name"
        scenic_nodes$id = NULL
        
        # create tf-target information
        regulon_info_txt = c()
        tf_targets_info_txt = c()
        for(i in 1:length(grn_network)){
          tf_txt = paste0("<b>TF: ", names(grn_network)[i], "</b><br/>")
          targets_txt = paste(grn_network[[i]], collapse = "\t")
          targets_txt = paste0("Targets: ", targets_txt, "<br/>")
          regulon_info_txt = paste0( tf_txt, targets_txt)
          tf_targets_info_txt = c(tf_targets_info_txt, regulon_info_txt)
        }
        tf_targets_info_txt = paste(tf_targets_info_txt, collapse = "<br/>")
        
        ## return
        list(scenic_nodes = scenic_nodes, 
             scenic_links = scenic_links, 
             tf_targets_info_txt = tf_targets_info_txt,
             tf_gene_cor_table = tf_gene_cor_table)
      })
    }
    
    # get grn activity of given genes
    if(T){
      GRN_profile_panel.gene_data = reactive({
        # get activity of given regulons, a data frame
        grn_activity = GRN_profile_panel.scenic_activity()$grn_activity[, GRN_profile_panel.genes$genes, drop = FALSE]
        
        # get the average activity of regulons, a vector
        avg_grn_activity = rowMeans(grn_activity)
        
        list(grn_activity = grn_activity, avg_grn_activity = avg_grn_activity)
      })
    }
    
    # get rss data
    if(T){
      GRN_profile_panel.rss_data = reactive({
        scenic_rss_data = rss_data[[input$GRN_profile_panel.major_cluster]]
        list(scenic_rss_data = scenic_rss_data)
      })
    }
    
    # display regulon activity in Scanpy embedding
    if(T){
      GRN_profile_panel.Scanpy_embedding.plot_activity.plot = reactiveValues(plot = NULL)
      output$GRN_profile_panel.Scanpy_embedding.plot_activity <- renderPlot({
        
        ## get annotation
        cellInfo = GRN_profile_panel.cellInfo()$cellInfo
        
        if(input$GRN_profile_panel.Embedding_used == "UMAP"){
          
          # get umap
          umap_data = GRN_profile_panel.scanpy_embedding()$umap[rownames(cellInfo), ] %>% as.data.frame
          colnames(umap_data) = c("UMAP_1","UMAP_2")
          
          # add average grn activity
          umap_data = cbind(umap_data, activity = GRN_profile_panel.gene_data()$avg_grn_activity)
          
          # rename colnames
          colnames(umap_data) = c("UMAP_1","UMAP_2","activity")
          
          # plot
          mycolor = colorRampPalette(brewer.pal(8, input$GRN_profile_panel.Colorpanel))(100)
          p = ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
            geom_point(aes(color = activity), size = input$GRN_profile_panel.Embedding_dotsize) +
            scale_colour_gradientn("Activity", colors = mycolor ) +
            xlab("UMAP_1") +
            ylab("UMAP_2") +
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") + 
            theme(axis.line = element_line(color = "black"))
          
          GRN_profile_panel.Scanpy_embedding.plot_activity.plot$plot = p
          print(p)
          
        }else if(input$GRN_profile_panel.Embedding_used == "tSNE"){
          
          # get umap
          tnse_data = GRN_profile_panel.scanpy_embedding()$tsne[rownames(cellInfo), ] %>% as.data.frame
          colnames(tnse_data) = c("TSNE_1","TSNE_2")
          
          # add average grn activity
          tnse_data = cbind(tnse_data, activity = GRN_profile_panel.gene_data()$avg_grn_activity)
          
          # rename colnames
          colnames(tnse_data) = c("TSNE_1","TSNE_2","activity")
          
          # plot
          mycolor = colorRampPalette(brewer.pal(8, input$GRN_profile_panel.Colorpanel))(100)
          p = ggplot(tnse_data, aes(x = TSNE_1, y = TSNE_2)) +
            geom_point(aes(color = activity), size = input$GRN_profile_panel.Embedding_dotsize) +
            scale_colour_gradientn("Activity", colors = mycolor ) +
            xlab("TSNE_1") +
            ylab("TSNE_2") +
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") + 
            theme(axis.line = element_line(color = "black"))
          
          GRN_profile_panel.Scanpy_embedding.plot_activity.plot$plot = p
          print(p)
        }
      })
      
      # download
      output$GRN_profile_panel.Scanpy_embedding.plot_activity.download_pdf <- downloadHandler(
        filename = function(){
          if(input$GRN_profile_panel.Scanpy_embedding.plot_activity.file_type == 'pdf'){
            paste0('GRN_profile_panel.Scanpy_embedding.plot_activity_', stringi::stri_rand_strings(1, 10), '.pdf')
          }else if(input$GRN_profile_panel.Scanpy_embedding.plot_activity.file_type == 'jpeg'){
            paste0('GRN_profile_panel.Scanpy_embedding.plot_activity_', stringi::stri_rand_strings(1, 10), '.jpeg')
          }
        },
        content = function(file){
          if(input$GRN_profile_panel.Scanpy_embedding.plot_activity.file_type == 'pdf'){
            pdf(file, width = input$GRN_profile_panel.Scanpy_embedding.plot_activity.width, height = input$GRN_profile_panel.Scanpy_embedding.plot_activity.height)
            plot(GRN_profile_panel.Scanpy_embedding.plot_activity.plot$plot)
            dev.off()
          }else if(input$GRN_profile_panel.Scanpy_embedding.plot_activity.file_type == 'jpeg'){
            jpeg(file, width = input$GRN_profile_panel.Scanpy_embedding.plot_activity.width,
                 height = input$GRN_profile_panel.Scanpy_embedding.plot_activity.height,
                 units = 'in', res = 300)
            plot(GRN_profile_panel.Scanpy_embedding.plot_activity.plot$plot)
            dev.off()
          }
        }
      )
    }
    
    # display Scanpy embedding
    if(T) {
      GRN_profile_panel.Scanpy_embedding.plot_label.plot = reactiveValues(plot = NULL)
      output$GRN_profile_panel.Scanpy_embedding.plot_label <- renderPlot({
        
        ## get annotation
        cellInfo = GRN_profile_panel.cellInfo()$cellInfo
        
        ## plot annotation in umap
        if(input$GRN_profile_panel.Embedding_used == "UMAP"){
          umap_data = GRN_profile_panel.scanpy_embedding()$umap
          umap_data = cbind( umap_data[ rownames(cellInfo), ], cellInfo)
          umap_data = umap_data[,1:3]
          colnames(umap_data) = c("UMAP_1","UMAP_2","label")
          
          ## get position of labels
          label_pos <- umap_data %>%
            group_by(label) %>%
            summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
          
          # set levels
          umap_data$label = factor(umap_data$label, 
                                   levels = all_cluster_list[[input$GRN_profile_panel.major_cluster]])
          
          # plot annotation
          mycolor = minor_cluster_color[umap_data$label %>% droplevels %>% levels]
          
          # plot
          p = ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = label)) +
            geom_point(size = input$GRN_profile_panel.Embedding_dotsize) +
            scale_color_manual("Label", values = mycolor) +
            xlab("UMAP_1") +
            ylab("UMAP_2") +
            ggrepel::geom_text_repel(data = label_pos, aes(label = label), colour = "black")+ ## add label
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "none") + 
            theme(axis.line = element_line(color = "black"))
          
          GRN_profile_panel.Scanpy_embedding.plot_label.plot$plot = p
          print(p)
          
        }else if(input$GRN_profile_panel.Embedding_used == "tSNE"){
          tsne_data = GRN_profile_panel.scanpy_embedding()$tsne
          tsne_data = cbind( tsne_data[ rownames(cellInfo), ], cellInfo)
          tsne_data = tsne_data[,1:3]
          colnames(tsne_data) = c("TSNE_1","TSNE_2","label")
          
          ## get position of labels
          label_pos <- tsne_data %>%
            group_by(label) %>%
            summarize(TSNE_1 = mean(TSNE_1), TSNE_2 = mean(TSNE_2))
          
          # set levels
          tsne_data$label = factor(tsne_data$label, 
                                   levels = all_cluster_list[[input$GRN_profile_panel.major_cluster]])
          # plot annotation
          mycolor = minor_cluster_color[tsne_data$label %>% droplevels %>% levels]
          
          p = ggplot(tsne_data, aes(x = TSNE_1, y = TSNE_2, color = label)) +
            geom_point(size = input$GRN_profile_panel.Embedding_dotsize) +
            scale_color_manual("Label", values = mycolor) +
            xlab("TSNE_1") +
            ylab("TSNE_2") +
            ggrepel::geom_text_repel(data = label_pos, aes(label = label), colour = "black")+ ## add label
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "none") + 
            theme(axis.line = element_line(color = "black"))
          
          GRN_profile_panel.Scanpy_embedding.plot_label.plot$plot = p
          print(p)
        }
      })
      
      # download
      output$GRN_profile_panel.Scanpy_embedding.plot_label.download_pdf <- downloadHandler(
        filename = function(){
          if(input$GRN_profile_panel.Scanpy_embedding.plot_label.file_type == 'pdf'){
            paste0('GRN_profile_panel.Scanpy_embedding.plot_label_', stringi::stri_rand_strings(1, 10), '.pdf')
          }else if(input$GRN_profile_panel.Scanpy_embedding.plot_label.file_type == 'jpeg'){
            paste0('GRN_profile_panel.Scanpy_embedding.plot_label_', stringi::stri_rand_strings(1, 10), '.jpeg')
          }
        },
        content = function(file){
          if(input$GRN_profile_panel.Scanpy_embedding.plot_label.file_type == 'pdf'){
            pdf(file, width = input$GRN_profile_panel.Scanpy_embedding.plot_label.width, height = input$GRN_profile_panel.Scanpy_embedding.plot_label.height)
            plot(GRN_profile_panel.Scanpy_embedding.plot_label.plot$plot)
            dev.off()
          }else if(input$GRN_profile_panel.Scanpy_embedding.plot_label.file_type == 'jpeg'){
            jpeg(file, width = input$GRN_profile_panel.Scanpy_embedding.plot_label.width,
                 height = input$GRN_profile_panel.Scanpy_embedding.plot_label.height,
                 units = 'in', res = 300)
            plot(GRN_profile_panel.Scanpy_embedding.plot_label.plot$plot)
            dev.off()
          }
        }
      )
    }
    
    # display Scenic embedding
    if(T) {
      GRN_profile_panel.Scenic_embedding.plot_label.plot = reactiveValues(plot = NULL)
      output$GRN_profile_panel.Scenic_embedding.plot_label <- renderPlot({
        ## get annotation
        cellInfo = GRN_profile_panel.cellInfo()$cellInfo
        
        if(input$GRN_profile_panel.Embedding_used == "UMAP"){
          
          # get umap
          umap_data = GRN_profile_panel.scenic_embedding()$umap
          
          # add annotation
          umap_data = cbind( umap_data[ rownames(cellInfo), ], cellInfo)
          
          # rename colnames
          umap_data = umap_data[,1:3]
          colnames(umap_data) = c("UMAP_1","UMAP_2","label")
          umap_data$label = factor(umap_data$label, 
                                   levels = all_cluster_list[[input$GRN_profile_panel.major_cluster]])
          
          # set color
          mycolor = minor_cluster_color[umap_data$label %>% droplevels %>% levels] 
          
          # plot
          p = ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = label)) +
            geom_point(size = input$GRN_profile_panel.Embedding_dotsize) + 
            scale_color_manual("Label", values = mycolor) + 
            xlab("UMAP_1") +
            ylab("UMAP_2") + 
            guides(colour = guide_legend(override.aes = list(size=3))) + 
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") + 
            theme(axis.line = element_line(color = "black"))
          
          GRN_profile_panel.Scenic_embedding.plot_label.plot$plot = p
          print(p)
          
        }else if(input$GRN_profile_panel.Embedding_used == "tSNE"){
          # get tsne
          tsne_data = GRN_profile_panel.scenic_embedding()$tsne
          
          # add annotation
          tsne_data = cbind( tsne_data[ rownames(cellInfo), ], cellInfo)
          
          # rename colnames
          tsne_data = tsne_data[,1:3]
          colnames(tsne_data) = c("TSNE_1","TSNE_2","label")
          tsne_data$label = factor(tsne_data$label, 
                                   levels = all_cluster_list[[input$GRN_profile_panel.major_cluster]])
          
          # set color
          mycolor = minor_cluster_color[tsne_data$label %>% droplevels %>% levels]
          
          # plot
          p = ggplot(tsne_data, aes(x = TSNE_1, y = TSNE_2, color = label)) +
            geom_point(size = input$GRN_profile_panel.Embedding_dotsize) + 
            scale_color_manual("Label", values = mycolor) + 
            xlab("TSNE_1") +
            ylab("TSNE_2") + 
            guides(colour = guide_legend(override.aes = list(size=3))) + 
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
                  axis.text.y = element_text(size=10, family="Times", color = "black"),
                  axis.title.x = element_text(size=12, family="Times", color = "black"),
                  axis.title.y = element_text(size=12, family="Times", color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") + 
            theme(axis.line = element_line(color = "black"))
          
          GRN_profile_panel.Scenic_embedding.plot_label.plot$plot = p
          print(p)
        }
      })
      
      # download
      output$GRN_profile_panel.Scenic_embedding.plot_label.download_pdf <- downloadHandler(
        filename = function(){
          if(input$GRN_profile_panel.Scenic_embedding.plot_label.file_type == 'pdf'){
            paste0('GRN_profile_panel.Scenic_embedding.plot_label_', stringi::stri_rand_strings(1, 10), '.pdf')
          }else if(input$GRN_profile_panel.Scenic_embedding.plot_label.file_type == 'jpeg'){
            paste0('GRN_profile_panel.Scenic_embedding.plot_label_', stringi::stri_rand_strings(1, 10), '.jpeg')
          }
        },
        content = function(file){
          if(input$GRN_profile_panel.Scenic_embedding.plot_label.file_type == 'pdf'){
            pdf(file, width = input$GRN_profile_panel.Scenic_embedding.plot_label.width, height = input$GRN_profile_panel.Scenic_embedding.plot_label.height)
            plot(GRN_profile_panel.Scenic_embedding.plot_label.plot$plot)
            dev.off()
          }else if(input$GRN_profile_panel.Scenic_embedding.plot_label.file_type == 'jpeg'){
            jpeg(file, width = input$GRN_profile_panel.Scenic_embedding.plot_label.width,
                 height = input$GRN_profile_panel.Scenic_embedding.plot_label.height,
                 units = 'in', res = 300)
            plot(GRN_profile_panel.Scenic_embedding.plot_label.plot$plot)
            dev.off()
          }
        }
      )
    }
    
    # violin plot of regulon activity
    if(T){
      GRN_profile_panel.Violin_plot.plot = reactiveValues(plot = NULL)
      output$GRN_profile_panel.Violin_plot <- renderPlot({
        # prepare data
        cellInfo = GRN_profile_panel.cellInfo()$cellInfo
        activity = GRN_profile_panel.gene_data()$avg_grn_activity
        activity_in_subtype = cbind(activity, cellInfo)
        activity_in_subtype = activity_in_subtype[,c(2,1)]
        colnames(activity_in_subtype) = c("label","activity")
        activity_in_subtype$activity = as.numeric(activity_in_subtype$activity)
        
        # get the average activity of multiple cell types
        group_mean = activity_in_subtype %>% group_by(label) %>% summarise(avg_activity = mean(activity))
        rnames = rownames(activity_in_subtype)
        
        # add group_mean to activity_in_subtype
        activity_in_subtype = activity_in_subtype %>% left_join(group_mean, by = "label")
        rownames(activity_in_subtype) = rnames
        
        # cell type order
        cell_type_order = levels(droplevels(adata_obs[adata_obs$major_cluster == input$GRN_profile_panel.major_cluster,]$minor))
        
        # violin plot
        p = ggplot( activity_in_subtype, aes(x = label, y = activity)) + 
          geom_violin(aes(fill = avg_activity), scale = "width", color = "white", kernel = "gaussian") + 
          # geom_boxplot(outlier.size = -1, width = .1, fill = "white") + 
          scale_fill_gradientn("Activity", colors = colorRampPalette(brewer.pal(8, input$GRN_profile_panel.Colorpanel))(100) ) +
          scale_x_discrete( limits = cell_type_order) + 
          ylab("Relative regulon activity") + 
          theme(axis.text.x = element_text(size = 10, family="Times", color = "black", angle= 90, vjust = 0.5, hjust = 1),
                axis.text.y = element_text(size=10, family="Times", color = "black"),
                axis.title.x = element_text(size=12, family="Times", color = "black"),
                axis.title.y = element_text(size=12, family="Times", color = "black"),     
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill='transparent', color='white'),
                panel.border = element_blank(),
                legend.position = "right") + 
          theme(axis.line = element_line(color = "black"))
        
        GRN_profile_panel.Violin_plot.plot$plot = p
        print(p)
      })
      
      # download
      output$GRN_profile_panel.Violin_plot.download_pdf <- downloadHandler(
        filename = function(){
          if(input$GRN_profile_panel.Violin_plot.file_type == 'pdf'){
            paste0('GRN_profile_panel.Violin_plot_', stringi::stri_rand_strings(1, 10), '.pdf')
          }else if(input$GRN_profile_panel.Violin_plot.file_type == 'jpeg'){
            paste0('GRN_profile_panel.Violin_plot_', stringi::stri_rand_strings(1, 10), '.jpeg')
          }
        },
        content = function(file){
          if(input$GRN_profile_panel.Violin_plot.file_type == 'pdf'){
            pdf(file, width = input$GRN_profile_panel.Violin_plot.width, height = input$GRN_profile_panel.Violin_plot.height)
            plot(GRN_profile_panel.Violin_plot.plot$plot)
            dev.off()
          }else if(input$GRN_profile_panel.Violin_plot.file_type == 'jpeg'){
            jpeg(file, width = input$GRN_profile_panel.Violin_plot.width,
                 height = input$GRN_profile_panel.Violin_plot.height,
                 units = 'in', res = 300)
            plot(GRN_profile_panel.Violin_plot.plot$plot)
            dev.off()
          }
        }
      )
    }
    
    # heat map of regulon activity
    if(T){
      GRN_profile_panel.Heatmap_plot.plot = reactiveValues(plot = NULL)
      output$GRN_profile_panel.Heatmap_plot <- renderPlot({
        ## get annotation
        cellInfo = GRN_profile_panel.cellInfo()$cellInfo
        activity = GRN_profile_panel.gene_data()$grn_activity
        activity_in_subtype = cbind(activity, cellInfo)
        activity_in_subtype = activity_in_subtype[, 1:(ncol(activity)+1)]
        colnames(activity_in_subtype) = c( colnames(activity), "label")
        activity_in_subtype = reshape2::melt(data = activity_in_subtype, id.vars = c("label") )
        colnames(activity_in_subtype) = c("label","gene","activity")
        activity_in_subtype$activity = as.numeric(activity_in_subtype$activity)
        
        ## summarize data
        activity_mean = activity_in_subtype %>%
          group_by(gene, label) %>%
          summarise(avg_activity = mean(activity))
        
        # cell type order
        cell_type_order = levels(droplevels(adata_obs[adata_obs$major_cluster == input$GRN_profile_panel.major_cluster,]$minor))
        
        #create heatmap using blue color scale
        p = ggplot(activity_mean, aes(label, gene)) +
          geom_tile(aes(fill = avg_activity), colour = "white") +
          scale_fill_gradientn("Activity", colors = colorRampPalette(brewer.pal(8, input$GRN_profile_panel.Colorpanel))(100)) + 
          scale_x_discrete( limits = cell_type_order) +
          theme_bw() + 
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10))
        
        GRN_profile_panel.Heatmap_plot.plot$plot = p
        print(p)
      })
      
      # download
      output$GRN_profile_panel.Heatmap_plot.download_pdf <- downloadHandler(
        filename = function(){
          if(input$GRN_profile_panel.Heatmap_plot.file_type == 'pdf'){
            paste0('GRN_profile_panel.Heatmap_plot_', stringi::stri_rand_strings(1, 10), '.pdf')
          }else if(input$GRN_profile_panel.Heatmap_plot.file_type == 'jpeg'){
            paste0('GRN_profile_panel.Heatmap_plot_', stringi::stri_rand_strings(1, 10), '.jpeg')
          }
        },
        content = function(file){
          if(input$GRN_profile_panel.Heatmap_plot.file_type == 'pdf'){
            pdf(file, width = input$GRN_profile_panel.Heatmap_plot.width, height = input$GRN_profile_panel.Heatmap_plot.height)
            print(GRN_profile_panel.Heatmap_plot.plot$plot)
            dev.off()
          }else if(input$GRN_profile_panel.Heatmap_plot.file_type == 'jpeg'){
            jpeg(file, width = input$GRN_profile_panel.Heatmap_plot.width,
                 height = input$GRN_profile_panel.Heatmap_plot.height,
                 units = 'in', res = 300)
            print(GRN_profile_panel.Heatmap_plot.plot$plot)
            dev.off()
          }
        }
      )
    }
    
    # display GRN network data, network
    if(T){
      GRN_profile_panel.plot_network.plot = reactiveValues(plot = NULL)
      output$GRN_profile_panel.plot_network <- renderForceNetwork({
        p = forceNetwork(Links = GRN_profile_panel.scenic_network()$scenic_links, 
                     Nodes = GRN_profile_panel.scenic_network()$scenic_nodes, Source = "source",
                     Target = "target", Value = "value", NodeID = "name",
                     Group = "group", opacity = 1, zoom = TRUE, opacityNoHover = TRUE, linkWidth = 0.5)
        GRN_profile_panel.plot_network.plot$plot = p
        print(p)
      })
      
      # download
      output$GRN_profile_panel.plot_network.download_pdf <- downloadHandler(
        filename = function(){paste0('GRN_profile_panel.plot_network_', stringi::stri_rand_strings(1, 10), '.html')},
        content = function(file){saveNetwork(GRN_profile_panel.plot_network.plot$plot, file)}
      )
    }
    
    # display GRN network data, text
    if(T){
      output$GRN_profile_panel.describe_network <- renderUI({
        #HTML(GRN_profile_panel.scenic_network()$tf_targets_info_txt)
        DT::renderDataTable({
          datatable(GRN_profile_panel.scenic_network()$tf_gene_cor_table, filter = "top", escape = FALSE)
        })
      })
    }
    
    # display rss data
    if(T){
      output$GRN_profile_panel.rss_tbl <- DT::renderDataTable({
        data = GRN_profile_panel.rss_data()$scenic_rss_data
        datatable(data, options = list(scrollX=TRUE, scrollCollapse=TRUE) ) %>%
          formatRound(columns = c(1:ncol(data)), digits = 2) %>%
          formatStyle(columns = c(1:ncol(data)), 'text-align' = 'middle')})
    }
  }
  
  # Regulon activity comparison
  if(T){
    
    # set minor cluster and tfs according to selected major cluster
    observe({
      # set regulons
      updatePickerInput(
        session = session,
        inputId = "GRN_comparison_panel.TF",
        choices = tf_list[[ input$GRN_comparison_panel.major_cluster ]],
        selected = tf_list[[ input$GRN_comparison_panel.major_cluster ]][1],
        choicesOpt = list(style = "background:white; color:black")
      )
      # set minor cluster
      updatePickerInput(
        session = session,
        inputId = "GRN_comparison_panel.minor_cluster",
        choices = all_cluster_list[[input$GRN_comparison_panel.major_cluster]],
        selected = all_cluster_list[[input$GRN_comparison_panel.major_cluster]],
        choicesOpt = list(style = "background:white; color:black"))
    })
    
    # get loom file
    if(T){
      # select loom
      GRN_comparison_panel.data = reactive({
        # return result
        list(scenic_loom = loom_data[[ input$GRN_comparison_panel.major_cluster ]])
      })
    }
  
    # select data
    # get regulon activity of given tf
    if(T){
      GRN_comparison_panel.subsets_data = reactive({
        
        selected = adata_obs$major_cluster == input$GRN_comparison_panel.major_cluster
        
        # subset by minor_cluster
        selected = selected & 
          adata_obs$minor_cluster %in% input$GRN_comparison_panel.minor_cluster
        
        # subset by stage
        selected = selected & 
          adata_obs$stage %in% input$GRN_comparison_panel.stage
        
        # subset by disease
        selected = selected &
          adata_obs$disease %in% input$GRN_comparison_panel.disease
        
        # subset by tissue
        selected = selected &
          adata_obs$tissue %in% input$GRN_comparison_panel.tissue
        
        # subset by tissue
        selected = selected &
          adata_obs$tissue.sub %in% input$GRN_comparison_panel.tissue.sub
        
        # subset by study
        selected = selected &
          adata_obs$study %in% input$GRN_comparison_panel.study
        
        # subset by sample
        selected = selected &
          adata_obs$sample %in% input$GRN_comparison_panel.sample
        
        # get selected cells
        selected_cells = rownames(adata_obs[selected, ])
        
        # get subset meta data
        selected_meta_data = adata_obs[selected_cells, ]
        
        # get auc matrix of regulons
        regulonsAUC = get_regulons_AUC( GRN_comparison_panel.data()$scenic_loom, 
                                        column.attr.name = "RegulonsAUC")
        
        # get grn activity
        grn_activity = getAUC(regulonsAUC) %>% t %>% as.data.frame
        
        # get grn activity of given TF
        grn_activity = grn_activity[selected_cells, input$GRN_comparison_panel.TF]
        
        # update selected_meta_data
        selected_meta_data$minor_cluster = droplevels(selected_meta_data$minor_cluster)
        selected_meta_data$stage = droplevels(selected_meta_data$stage)
        selected_meta_data$disease = droplevels(selected_meta_data$disease)
        selected_meta_data$tissue = droplevels(selected_meta_data$tissue)
        selected_meta_data$tissue.sub = droplevels(selected_meta_data$tissue.sub)
        selected_meta_data$study = droplevels(selected_meta_data$study)
        selected_meta_data$sample = droplevels(selected_meta_data$sample)
        
        # return expression
        list(grn_activity = grn_activity, 
             selected_meta_data = selected_meta_data)
      })
    }
    
    # plot regulon activity of selected cells
    # violin plot
    # group by one variable, major/minor, stage/disease/tissue
    if(T){
      GRN_comparison_panel.Violin_plot.plot = reactiveValues(plot = NULL)
      output$GRN_comparison_panel.Violin_plot <- renderPlot({
        
        # set color
        mycolor = colorRampPalette(brewer.pal(8, input$GRN_comparison_panel.Colorpanel))(100)
        
        # prepare data
        # get label
        label = GRN_comparison_panel.subsets_data()$selected_meta_data[, input$GRN_comparison_panel.groupby]
        activity = GRN_comparison_panel.subsets_data()$grn_activity
        
        # prepare data
        activity_in_subtype = data.frame(group = label, activity = activity)
        
        # get average activity by group
        group_avg_exp = activity_in_subtype %>% group_by(group) %>% summarise(avg_activity = mean(activity))
        rnames = rownames(activity_in_subtype)
        
        activity_in_subtype = activity_in_subtype %>% left_join(group_avg_exp, by = c("group"))
        rownames(activity_in_subtype) = rnames
        
        # violin plot
        p = ggplot( activity_in_subtype, aes(x = group, y = activity)) + 
          geom_violin(aes(fill = avg_activity), scale = "width", color = "white", kernel = "gaussian") + 
          #geom_boxplot(outlier.size = -1, width = .1, fill = "white") + 
          scale_fill_gradientn("Activity", colors = mycolor) +
          ylab("Relative regulon activity") + 
          xlab("")+
          theme_bw() + 
          theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                axis.text.y = element_text(size=10,  color = "black"),
                axis.title.x = element_text(size=12, color = "black"),
                axis.title.y = element_text(size=12,  color = "black"),     
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill='transparent', color='white'),
                panel.border = element_blank(),
                legend.position = "right") +
          theme(axis.line = element_line(color = "black"))
        
        GRN_comparison_panel.Violin_plot.plot$plot = p
        print(p)
      })
      
      # download
      output$GRN_comparison_panel.Violin_plot.download_pdf <- downloadHandler(
        filename = function(){
          if(input$GRN_comparison_panel.Violin_plot.file_type == 'pdf'){
            paste0('GRN_comparison_panel.Violin_plot_', stringi::stri_rand_strings(1, 10), '.pdf')
          }else if(input$GRN_comparison_panel.Violin_plot.file_type == 'jpeg'){
            paste0('GRN_comparison_panel.Violin_plot_', stringi::stri_rand_strings(1, 10), '.jpeg')
          }
        },
        content = function(file){
          if(input$GRN_comparison_panel.Violin_plot.file_type == 'pdf'){
            pdf(file, width = input$GRN_comparison_panel.Violin_plot.width, height = input$GRN_comparison_panel.Violin_plot.height)
            plot(GRN_comparison_panel.Violin_plot.plot$plot)
            dev.off()
          }else if(input$GRN_comparison_panel.Violin_plot.file_type == 'jpeg'){
            jpeg(file, width = input$GRN_comparison_panel.Violin_plot.width,
                 height = input$GRN_comparison_panel.Violin_plot.height,
                 units = 'in', res = 300)
            plot(GRN_comparison_panel.Violin_plot.plot$plot)
            dev.off()
          }
        }
      )
    }
    
    # plot cell number of selected cells
    # bar plot
    # group by major/minor/stage/disease/tissue
    if(T){
      GRN_comparison_panel.Cell_number.Bar_plot.plot = reactiveValues(plot = NULL)
      output$GRN_comparison_panel.Cell_number.Bar_plot <- renderPlot({
        
        # sleep 2s
        Sys.sleep("2")
        
        # prepare data
        selected_cells = GRN_comparison_panel.subsets_data()$selected_meta_data %>%
          select(input$GRN_comparison_panel.groupby)
        colnames(selected_cells) = c("group")
        
        # count, group by group
        cell_number_count = selected_cells %>% group_by(group) %>% summarise(count = n())
        
        # calculate percentage of cells
        cell_number_pct = cell_number_count
        cell_number_pct$pct = cell_number_count$count/sum(cell_number_count$count)
        
        # barplot
        if(input$GRN_comparison_panel.Cell_number.Bar_plot.type == 'Number'){
          p = ggplot( cell_number_count, aes(x = group, y = count, fill = group)) + 
            geom_bar(stat = "identity") +
            ylab("Number of cells") + 
            xlab("")+
            theme_bw() +
            theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size=10,  color = "black"),
                  axis.title.x = element_text(size=12, color = "black"),
                  axis.title.y = element_text(size=12,  color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") +
            theme(axis.line = element_line(color = "black"),
                  legend.position = "none")
        }else if(input$GRN_comparison_panel.Cell_number.Bar_plot.type == 'Percentage'){
          p = ggplot( cell_number_pct, aes(x = group, y = pct, fill = group)) + 
            geom_bar(stat = "identity") +
            ylab("Percentage of cells") + 
            xlab("")+
            theme_bw() +
            theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size=10,  color = "black"),
                  axis.title.x = element_text(size=12, color = "black"),
                  axis.title.y = element_text(size=12,  color = "black"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") +
            theme(axis.line = element_line(color = "black"),
                  legend.position = "none")
        }
        
        GRN_comparison_panel.Cell_number.Bar_plot.plot$plot = p
        print(p)
      })
      
      # download
      output$GRN_comparison_panel.Cell_number.Bar_plot.download_pdf <- downloadHandler(
        filename = function(){
          if(input$GRN_comparison_panel.Cell_number.Bar_plot.file_type == 'pdf'){
            paste0('GRN_comparison_panel.Cell_number.Bar_plot_', stringi::stri_rand_strings(1, 10), '.pdf')
          }else if(input$GRN_comparison_panel.Cell_number.Bar_plot.file_type == 'jpeg'){
            paste0('GRN_comparison_panel.Cell_number.Bar_plot_', stringi::stri_rand_strings(1, 10), '.jpeg')
          }
        },
        content = function(file){
          if(input$GRN_comparison_panel.Cell_number.Bar_plot.file_type == 'pdf'){
            pdf(file, width = input$GRN_comparison_panel.Cell_number.Bar_plot.width, height = input$GRN_comparison_panel.Cell_number.Bar_plot.height)
            plot(GRN_comparison_panel.Cell_number.Bar_plot.plot$plot)
            dev.off()
          }else if(input$GRN_comparison_panel.Cell_number.Bar_plot.file_type == 'jpeg'){
            jpeg(file, width = input$GRN_comparison_panel.Cell_number.Bar_plot.width,
                 height = input$GRN_comparison_panel.Cell_number.Bar_plot.height,
                 units = 'in', res = 300)
            plot(GRN_comparison_panel.Cell_number.Bar_plot.plot$plot)
            dev.off()
          }
        }
      )
    }
    
    # diff regulon, UC vs healthy
    if(T){
      output$GRN_comparison_panel.diff_regulon.healthy_UC.tbl <- DT::renderDataTable(
        datatable(diff_regulon_data[[input$GRN_comparison_panel.major_cluster]][[1]],
                  filter = "top", escape = FALSE, rownames = FALSE) %>% formatRound(2:6, 6)
      )
    }
    
    # diff regulon, CD vs healthy
    if(T){
      output$GRN_comparison_panel.diff_regulon.healthy_CD.tbl <- DT::renderDataTable(
        datatable(diff_regulon_data[[input$GRN_comparison_panel.major_cluster]][[2]],
                  filter = "top", escape = FALSE, rownames = FALSE) %>% formatRound(2:6, 6)
        )
    }
    
  }
  
  # Cellular composition
  if(T){
    # observe samples selected
    observe({
      updatePickerInput(
        session = session,
        inputId = "cellular_composition_panel.study",
        choices = all_study_list,
        selected = all_study_list,
        choicesOpt = list(style = "background:white; color:black")
      )
      updatePickerInput(
        session = session,
        inputId = "cellular_composition_panel.sample",
        choices = all_sample_list,
        selected = unlist(all_sample_list),
        choicesOpt = list(style = "background:white; color:black")
      )
    })
    
    # update study and samples selected
    if(T){
      # default selection
      cellular_composition_panel.samples = reactiveValues(samples = head(unlist(all_sample_list), n= 50))
      cellular_composition_panel.studies = reactiveValues(studies = all_study_list)
      
      # update selection once click Submit
      observeEvent(input$cellular_composition_panel.Submit,{
        cellular_composition_panel.samples$samples = input$cellular_composition_panel.sample
        cellular_composition_panel.studies$studies = input$cellular_composition_panel.study
      })
    }
      
    # select samples
    if(T){
      cellular_composition_panel.select_meta_data = reactive({
        
        # select study
        select = adata_obs$study %in% cellular_composition_panel.studies$studies
        
        # select samples
        select = select & 
          adata_obs$sample %in% cellular_composition_panel.samples$samples
        
        select_meta_data = adata_obs[select, ]
        
        # return expression
        list(select_meta_data = select_meta_data)
      })
    }
    
    # stacked bar plot, major cluster
    if(T){
      output$cellular_composition_panel.output.major_barplot <- renderPlotly({
        
        # prepare data
        total = cellular_composition_panel.select_meta_data()$select_meta_data %>% 
          group_by(sample) %>% 
          summarise(total = n())
        major_count = cellular_composition_panel.select_meta_data()$select_meta_data %>% 
          group_by(sample, major_cluster) %>%
          summarise(count = n())
        major_pct = left_join(major_count, total, by = "sample") %>% 
          mutate(percentage = count/total) %>% as.data.frame()
        
        # plot
        p = ggplot(major_pct, aes(x = sample, y = percentage, fill = major_cluster)) + 
          geom_histogram(stat = "identity") +
          ylim(0,1) + 
          xlab("") + 
          ylab("Percentage") +
          scale_fill_manual(values = major_cluster_color) + 
          theme_classic() + 
          theme(axis.text.x =  element_text(size = 10, color = "black", angle = 90, vjust = 1, hjust = 1))
        ggplotly(p)
      })
    }
    
    # stacked bar plot, minor cluster
    if(T){
      output$cellular_composition_panel.output.minor_barplot <- renderPlotly({
        
        # prepare data
        total = cellular_composition_panel.select_meta_data()$select_meta_data %>% 
          group_by(sample) %>% 
          summarise(total = n())
        minor_count = cellular_composition_panel.select_meta_data()$select_meta_data %>% 
          group_by(sample, minor_cluster) %>% 
          summarise(count = n())
        minor_pct = left_join(minor_count, total, by = "sample") %>% 
          mutate(percentage = count/total) %>% as.data.frame()
        
       # plot
       p = ggplot(minor_pct, aes(x = sample, y = percentage, fill = minor_cluster)) + 
          geom_histogram(stat = "identity") +
          ylim(0,1) + 
          xlab("") + 
          ylab("Percentage") +
          theme_classic() + 
          theme(axis.text.x =  element_text(size = 10, color = "black", angle = 90, vjust = 1, hjust = 1))
       ggplotly(p)
      })
    }
}
  
  # Therapy panel
  if(T){
    # Therapy panel, FDA approved table
    therapy_panel.FDA_approved_tbl.buttons <- list()
    therapy_panel.FDA_approved_tbl.buttonIds <- list()
    for (r in rownames(therapy_panel.FDA_approved_drugs)) {
      id <- r
      therapy_panel.FDA_approved_tbl.buttonIds[[r]] <- id
      button <- actionButton(id, label = "View", 
                             onclick = 'Shiny.onInputChange(\"select_button\",  [this.id, Math.random()])',
                             width = "80px", height = "50px", 
                             style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")
      therapy_panel.FDA_approved_tbl.buttons[[r]] <- as.character(button)
    }
    therapy_panel.FDA_approved_drugs$`Clinical trails` = therapy_panel.FDA_approved_tbl.buttons
    output$therapy_panel.FDA_approved_tbl = DT::renderDataTable(
      therapy_panel.FDA_approved_drugs, escape = FALSE, rownames = FALSE)
      
    # more info
    therapy_panel.FDA_approved_tbl.more.buttons <- list()
    therapy_panel.FDA_approved_tbl.more.buttonIds <- list()
    for (r in 1:nrow(therapy_panel.FDA_approved_drugs_detail)) {
      id = paste0("more_",r)
      r = id
      therapy_panel.FDA_approved_tbl.more.buttonIds[[r]] <- id
      button <- actionButton(id, label = "View",
                             onclick = 'Shiny.onInputChange(\"select_button_more\",  [this.id, Math.random()])',
                             width = "80px", height = "50px", 
                             style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")
      therapy_panel.FDA_approved_tbl.more.buttons[[r]] <- as.character(button)
    }
    therapy_panel.FDA_approved_drugs_detail$More = therapy_panel.FDA_approved_tbl.more.buttons
      
    # Therapy panel, FAD approved table, clinical information tables for each drug
    observeEvent(input$select_button, {
      showModal( modalDialog({
        # select data
        detail = therapy_panel.FDA_approved_drugs_detail[therapy_panel.FDA_approved_drugs_detail[,1] == input$select_button[1],]
        detail = detail[,c("Trial acronym","Clinical trail","Stage","Year","Reference","More")]
        # display data
        DT::renderDataTable(detail, escape = FALSE, rownames = FALSE,
                            options = list(ordering=TRUE))},
        title = paste0("Clinical trails for ", input$select_button[1]),
        size = "l",
        footer = modalButton("Close")))
    })
  
    # Therapy panel, FAD approved table, clinical information tables for each drug, more information
    observeEvent(input$select_button_more, {
      showModal( modalDialog({
        # select data
        detail = data.frame(t(therapy_panel.FDA_approved_drugs_detail[input$select_button_more[1],1:(ncol(therapy_panel.FDA_approved_drugs_detail)-1)]))
        colnames(detail) = c("Value")
        DT::renderDataTable(detail, escape = FALSE, rownames = TRUE,
                            options = list(ordering=FALSE))},
        title = paste0("Detailed clinical information"), 
        size = "l",
        footer = actionButton("restoreModal",label = "Close")))
    })
    
    # restore first modal
    observeEvent(input$restoreModal, {
      showModal( modalDialog({
        # select data
        detail = therapy_panel.FDA_approved_drugs_detail[therapy_panel.FDA_approved_drugs_detail[,1] == input$select_button[1],]
        detail = detail[,c("Trial acronym","Clinical trail","Stage","Year","Reference","More")]
        # display data
        DT::renderDataTable(detail, escape = FALSE, rownames = FALSE, 
                            options = list(ordering=TRUE))},
        title = paste0("Clinical trails for ", input$select_button[1]),
        size = "l",
        footer = modalButton("Close")))
    })
    
    # drugs and targets for IBD
    if(T){
      output$therapy_panel.ibd_targets_tbl = DT::renderDataTable( 
          datatable(therapy_panel.drug_and_target, 
                    filter = 'top', escape = FALSE, rownames = FALSE))
    }
  }
  
  # Risk gene panel
  if(T){
    if(T){ # gwas study
      output$ibd_gwas_tbl = DT::renderDataTable(
        datatable(ibd_gwas_data, 
                  filter = "top", escape = FALSE, rownames = FALSE))
    }
    
    if(T){ # risk genes
      output$adult_ibd_risk_genes_tbl = DT::renderDataTable(
        datatable(adult_ibd_risk_genes_data,
                  filter = "top", escape = FALSE, rownames = FALSE))
    }
    
    if(T){ # pediatric risk genes
      output$pediatric_ibd_risk_genes_tbl = DT::renderDataTable(
        datatable(pediatric_ibd_risk_genes_data))
    }
  }

  # projects_info_table
  if(T){
    output$projects_info_table = DT::renderDataTable(
      datatable(project_info, escape = FALSE, rownames = FALSE))
  }
  
  # gene set enrichment analysis
  if(T){

    # Update Gsea_panel.file.text
    if(T){
      observe({
      if (length(input$Gsea_panel.input.file$datapath) != 0) {
        # process input file
        genes_in_file <- read.table( input$Gsea_panel.input.file$datapath, header = F, stringsAsFactors = F)
        genes_in_file = genes_in_file$V1
        
        # update text
        updateTextAreaInput(session = session,
                            inputId = "Gsea_panel.file.text",
                            value = paste(genes_in_file, collapse = "\n"))
      }})
    }
    
    # Update Gsea_panel.select.text
    if(T){
      observe({
        if(length(input$Gsea_panel.select.predefined) != 0){
          gwas_genes = c()
          for(i in as.numeric(input$Gsea_panel.select.predefined)){
            gwas_genes = c(gwas_genes, gwas_gene_list[[i]])
          }
          genes_selected = unique(gwas_genes)
          genes_selected = intersect(genes_selected, all_genes)
          
        # update text
        updateTextAreaInput(session = session,
                            inputId = "Gsea_panel.select.text",
                            value = paste(genes_selected, collapse = "\n"))
      }})
    }
    
    # get input gene set, then run gsea
    if(T){
      observeEvent(input$Gsea_panel.Submit, {
      
      genes = c()
      # get input gene set
      if(input$Gsea_input_panel == "Genes"){
        # process input gene list
        genes = input$Gsea_panel.input.text
        genes = str_split(genes,"\n")[[1]]
        genes = str_trim(genes, side = c("both"))
        genes = genes[genes != ""]
        genes = unique(genes)
        genes = toupper(genes)
        genes = intersect(genes, all_genes)
      }else if(input$Gsea_input_panel == "Upload"){
        genes = input$Gsea_panel.file.text
        genes = str_split(genes,"\n")[[1]]
        genes = str_trim(genes, side = c("both"))
        genes = genes[genes != ""]
        genes = unique(genes)
        genes = toupper(genes)
        genes = intersect(genes, all_genes)
      }else if(input$Gsea_input_panel == "Select"){
        genes = input$Gsea_panel.select.text
        genes = str_split(genes,"\n")[[1]]
        genes = str_trim(genes, side = c("both"))
        genes = genes[genes != ""]
        genes = unique(genes)
        genes = toupper(genes)
        genes = intersect(genes, all_genes)
      }
      
      # gene set enrichment analysis
      gsea.result = data.frame()
      if( length(genes) != 0){
        input_geneset = genes
        gsea.result = run_gsea(minor_cluster_list = adata_obs$minor_cluster %>% levels, 
                               expressed_genes, input_geneset, markers_df)
      }
      
      # display result
      if(T){
        output$Gsea_panel.output.tbl <- DT::renderDataTable(
          datatable(gsea.result, escape = FALSE, rownames = F)
        )
      }
      
      # display heatmap
      if(T){
        Gsea_panel.output.heatmap.plot = reactiveValues(plot = NULL)
        output$Gsea_panel.output.heatmap = renderPlot({
          p = plot_heatmap.minor_cluster(adata, cells = rownames(adata_obs), genes = genes, mj_color = major_cluster_color)
          Gsea_panel.output.heatmap.plot$plot = p
          print(p)
        }, height = length(genes)*10 + 300)
        
        # download
        output$Gsea_panel.output.heatmap.download_pdf <- downloadHandler(
          filename = function(){
            if(input$Gsea_panel.output.heatmap.file_type == 'pdf'){
              paste0('Gsea_panel.output.heatmap_', stringi::stri_rand_strings(1, 10), '.pdf')
            }else if(input$Gsea_panel.output.heatmap.file_type == 'jpeg'){
              paste0('Gsea_panel.output.heatmap_', stringi::stri_rand_strings(1, 10), '.jpeg')
            }
          },
          content = function(file){
            if(input$Gsea_panel.output.heatmap.file_type == 'pdf'){
              pdf(file, width = input$Gsea_panel.output.heatmap.width, height = input$Gsea_panel.output.heatmap.height)
              print(Gsea_panel.output.heatmap.plot$plot)
              dev.off()
            }else if(input$Gsea_panel.output.heatmap.file_type == 'jpeg'){
              jpeg(file, width = input$Gsea_panel.output.heatmap.width,
                   height = input$Gsea_panel.output.heatmap.height,
                   units = 'in', res = 300)
              print(Gsea_panel.output.heatmap.plot$plot)
              dev.off()
            }
          }
        )
       }
      })
    }
  }
  
  # Meta-data panel
  if(T){
    output$meta_data_panel = DT::renderDataTable(sample_meta_data)
  }
  
  # Download panel
  if(T){
    output$download_panel.download_gex = downloadHandler(
      filename = "scIBD.gex_matrix.h5ad",
      content = function(file){
        file.copy("./www/download/scIBD.gex_matrix.h5ad", file)}
    )
    
    output$download_panel.download_meta = downloadHandler(
      filename = "scIBD_sample_meta_data.csv",
      content = function(file){
        file.copy("./www/download/scIBD_sample_meta_data.csv", file)}
    )
    
    output$download_panel.download_deg = downloadHandler(
      filename = "scIBD.deg_by_major_cluster.zip",
      content = function(file){
        file.copy("./www/download/scIBD.deg_by_major_cluster.zip", file)}
    )
    
    output$download_panel.download_loom = downloadHandler(
      filename = "scIBD.regulon_by_major_cluster.zip",
      content = function(file){
        file.copy("./www/download/scIBD.regulon_by_major_cluster.zip", file)}
    )
    
    output$download_panel.download_rss = downloadHandler(
      filename = "scIBD.rss_by_major_cluster.zip",
      content = function(file){
        file.copy("./www/download/scIBD.rss_by_major_cluster.zip", file)}
    )
    
    output$download_panel.download_diff_regulon = downloadHandler(
      filename = "scIBD.diff_regulon.zip",
      content = function(file){
        file.copy("./www/download/scIBD.diff_regulon.zip", file)}
    )
  }
}

# Run the application
shinyApp(ui = ui, server = server)
# Done
