# Pre loading
if(T){
  
  # Load packages
  if(T){
    # install package from CRAN
    packages <- c("shiny", "shinycssloaders", "shinythemes", "shinydashboard", "shinyWidgets", "shinyjs", "shinyBS",
        'downloadthis',"bsplus","DT", "dplyr","reshape2","data.table", "stringr", "readxl", 
        "ggplot2","ggpubr", "aplot", "ggtree", "plotly","ggrepel","RColorBrewer","networkD3","scales","Seurat","anndata")
    suppressMessages(lapply(packages, getPackage))
    
    # install package from github
    suppressMessages(getPackage(pkg = 'SCopeLoomR', github = "aertslab/SCopeLoomR"))
    suppressMessages(getPackage(pkg = "SCENIC", github = "aertslab/SCENIC"))
    
    # install package from Bioconductor
    suppressMessages(getPackage(pkg = "GeneOverlap", bioc = "GeneOverlap"))
    suppressMessages(getPackage(pkg = "ggtree", bioc = "ggtree"))
  }
  
  # Define variables and load data
  if(T){
    
    # Use H5AD object
    if(T){
      # connect to data
      adata = anndata::read_h5ad("./www/scanpy/all.clean.h5ad")
      
      # get meta data
      adata_obs = adata$obs
      
      # get genes
      all_genes = adata$var %>% rownames %>% sort
    }
    
    # get major clusters and minor clusters
    if(T){
      
      # get major cluster
      major_cluster = adata_obs[,"major_cluster"] %>% levels
      
      # get minor cluster
      all_cluster_list = list()
      for(i in 1:length(major_cluster)){
        minor_cluster =  adata_obs[ adata_obs$major_cluster == major_cluster[i], ]$minor_cluster %>% 
          droplevels  %>% levels
        all_cluster_list[[i]] = minor_cluster
      }
      names(all_cluster_list) = major_cluster
    }
    
    # set color for major cluster and minor cluster
    # generate scIBD_consecutive_color, major_cluster_color, minor_cluster_color variables
    if(T){
      # set color for consecutive values
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
    
    # get meta data
    # get group information
    if(T){
      # get stage, disease, tissue, study
      all_stage_list = levels(adata_obs$stage)
      all_disease_list = levels(adata_obs$disease)
      all_tissue_list = levels(adata_obs$tissue)
      all_study_list = levels(adata_obs$study)
      
      # get sample list
      all_sample_list = list()
      adata_obs$sample = factor(adata_obs$sample)
      for(i in 1:length(all_study_list)){
        samples = adata_obs[adata_obs$study == all_study_list[i], ]$sample %>% 
          droplevels  %>% levels
        all_sample_list[[i]] = samples
      }
      names(all_sample_list) = all_study_list
      all_samples = unlist(all_sample_list)
      
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
    
    # load study track info
    if(T){
      study_tracking_info = read_xlsx("./www/meta/study_tracking_info.xlsx",sheet=1) %>% as.data.frame()
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
                          ibd_gwas_data$Disease, sep = ",")
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
    if(T){
      deg_data = readRDS("./www/deg/diff_gex.rds")
      for(cluster in names(deg_data)){
        deg_data[[cluster]]$cluster = factor(deg_data[[cluster]]$cluster, levels = all_cluster_list[[cluster]])
      }
    }
    
    # get scanpy global embedding
    if(T){
      # get umap
      global_umap = adata_obs %>% dplyr::select(gUMAP_1, gUMAP_2, major_cluster) %>% as.data.frame()
      
      # get tsne
      global_tsne = adata_obs %>% dplyr::select(gTSNE_1, gTSNE_2, major_cluster) %>% as.data.frame()
    }
    
    # load scenic data
    if(T){
      loom_file = read.table("./www/scenic/loom.txt", header = F, sep = "\t")
      loom_data = list()
      for(i in 1:nrow(loom_file)){
        loom_data[[i]] = open_loom(loom_file[i,2])
      }
      names(loom_data) = loom_file[,1]
      rm(loom_file)
    }
    
    # load regulon-related data
    if(T){
      # load differential regulon
      diff_regulon_data = readRDS("./www/scenic/diff_regulon.rds") 
      
      # load tf list
      tf_list = readRDS("./www/scenic/tf_list.rds")
      
      # load rss data
      rss_file = read.table("./www/scenic/rss.txt", header = F)
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
    
    # load markdown files
    if(T){
      manual_page_markdown = readLines("./www/document/scIBD_documentation.md", warn = FALSE)
      Case_study_markdown = readLines("./www/document/scIBD_case_study.md", warn = FALSE)
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
                                                          style = "background:white; color:black")) %>% 
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Note:",
                                      placement = "left",
                                      content = "Click submit to update results.")),
                              pickerInput(
                                inputId = "GEX_profile_panel.Gene_list",
                                label = "Select genes",
                                choices = all_genes,
                                selected = c("TPSAB1","CPA3"),
                                multiple = TRUE,
                                options = list(`actions-box` = TRUE, 
                                               `live-search` = TRUE, 
                                               size = 8, 
                                               style = "background:white; color:black")) %>% 
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Note:",
                                      placement = "left",
                                      content = "Click submit to update results."))
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
                                selected = "Blues",
                                multiple = FALSE,
                                options = list(`actions-box` = TRUE, 
                                               `live-search` = TRUE, 
                                               size = 8, 
                                               style = "background:white; color:black"))
                      ),
                      column(width = 4,
                             sliderInput(inputId = "GEX_profile_panel.Embedding_dotsize", 
                                         label = "Dot size",
                                         value = 0.1, 
                                         min = 0, max = 1,step = 0.05),
                               column(width = 8,
                                      pickerInput(inputId = "GEX_profile_panel.downsampled_cell_num", 
                                                  label = "Downsampled cells per subtype", 
                                                  choices = c(Inf, seq(200,5000,200)),
                                                  selected = 400, 
                                                  multiple = FALSE,
                                                  options = list(`actions-box` = TRUE, 
                                                                 `live-search` = TRUE, 
                                                                 size = 8, 
                                                                 style = "background:white; color:black")) %>% 
                                        shinyInput_label_embed(
                                          shiny_iconlink() %>%
                                            bs_embed_popover(
                                              title = "Note:",
                                              placement = "left",
                                              content = "Select Inf to use all cells.\nClick submit to update results."))
                               ),
                               column(width = 4,
                                      tags$div(tags$p(tags$b("Submit"))),
                                      actionButton(inputId = "GEX_profile_panel.Submit", 
                                                   label = "GO", 
                                                   width = "125px", 
                                                   icon = icon("paper-plane"),
                                                   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
                     )),
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
                                              style = "background:white; color:black")) %>% 
                               shinyInput_label_embed(
                                 shiny_iconlink() %>%
                                   bs_embed_popover(
                                     title = "Note:",
                                     placement = "left",
                                     content = "Click submit to update results.")),
                             pickerInput(
                               inputId = "GRN_profile_panel.tfs",
                               label = "Transcript factor",
                               choices = tf_list[[ major_cluster[1] ]],
                               selected = "RUNX2",
                               multiple = TRUE,
                               options = list(`actions-box` = TRUE, 
                                              `live-search` = TRUE, 
                                              size = 8, 
                                              style = "background:white; color:black")) %>% 
                               shinyInput_label_embed(
                                 shiny_iconlink() %>%
                                   bs_embed_popover(
                                     title = "Note:",
                                     placement = "left",
                                     content = "Click submit to update results."))
                             ),
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
                                selected = "Blues",
                                multiple = FALSE,
                                options = list(`actions-box` = TRUE, 
                                               `live-search` = TRUE, 
                                               size = 8, 
                                               style = "background:white; color:black"))),
                      column(width = 4,
                             sliderInput(inputId = "GRN_profile_panel.Embedding_dotsize", 
                                         label = "Dot size",
                                         value = 0.1, 
                                         min = 0, max = 1,step = 0.05),
                             column(width =8,
                                    pickerInput(inputId = "GRN_profile_panel.downsampled_cell_num", 
                                                label = "Downsampled cells per subtype", 
                                                choices = c(Inf, seq(200,5000,200)),
                                                selected = 400, 
                                                multiple = FALSE,
                                                options = list(`actions-box` = TRUE, 
                                                               `live-search` = TRUE, 
                                                               size = 8, 
                                                               style = "background:white; color:black")) %>% 
                                      shinyInput_label_embed(
                                        shiny_iconlink() %>%
                                          bs_embed_popover(
                                            title = "Note:",
                                            placement = "left",
                                            content = "Select Inf to use all cells. Click submit to update results."))
                             ),
                             column(width = 4,
                                    tags$div(tags$p(tags$b("Submit"))),
                                    actionButton(inputId = "GRN_profile_panel.Submit", 
                                                 label = "GO", 
                                                 width = "150px", 
                                                 icon = icon("paper-plane"),
                                                 style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
                      )
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
                        tooltip = tooltipOptions(title = "Click to see inputs!")
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
                        tooltip = tooltipOptions(title = "Click to see inputs!")
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
                      fluidRow(
                        column(width = 6,
                               
                               pickerInput(inputId = "GRN_profile_panel.targets_num", 
                                           label = "Top expressed targets:", 
                                           choices = c(Inf, seq(1, 50, 1)),
                                           selected = 10, 
                                           multiple = FALSE,
                                           options = list(`actions-box` = TRUE, 
                                                          `live-search` = TRUE, 
                                                          size = 8, 
                                                          style = "background:white; color:black"))
                        )),
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
                            tags$div(tags$p(tags$b("Submit"))),
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
                              selected = "Myeloid",
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
                            pickerInput(inputId = "GEX_comparison_panel.downsampled_cell_num", 
                                        label = "Downsampled cells per subtype", 
                                        choices = c(Inf, seq(200,5000,200)),
                                        selected = 400, 
                                        multiple = FALSE,
                                        options = list(`actions-box` = TRUE, 
                                                       `live-search` = TRUE, 
                                                       size = 8, 
                                                       style = "background:white; color:black")) %>% 
                              shinyInput_label_embed( shiny_iconlink() %>%
                                                        bs_embed_popover( title = "Note:", placement = "left",
                                                                          content = "Select Inf to use all cells.")),
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
                              selected = all_samples,
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
                              selected = "Blues",
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
                                          selected = "minor_cluster", 
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
                               choices = "ATF3",
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
                               selected = all_samples,
                               multiple = TRUE,
                               options = list(`actions-box` = TRUE, 
                                              `live-search` = TRUE, 
                                              size = 8,
                                              style = "background:white; color:black")),
                             pickerInput(
                               inputId = "GRN_comparison_panel.Colorpanel",
                               label = "Color profile",
                               choices = scIBD_consecutive_color,
                               selected = "Blues",
                               multiple = FALSE,
                               options = list(`actions-box` = TRUE, 
                                              `live-search` = TRUE, 
                                              size = 8, 
                                              style = "background:white; color:black"))
                      ),
                      column(width = 2,
                             pickerInput(inputId = "GRN_comparison_panel.downsampled_cell_num", 
                                         label = "Downsampled cells per subtype", 
                                         choices = c(Inf, seq(200,5000,200)),
                                         selected = 400, 
                                         multiple = FALSE,
                                         options = list(`actions-box` = TRUE, 
                                                        `live-search` = TRUE, 
                                                        size = 8, 
                                                        style = "background:white; color:black")) %>% 
                               shinyInput_label_embed( shiny_iconlink() %>%
                                                         bs_embed_popover( title = "Note:", placement = "left",
                                                                           content = "Select Inf to use all cells.")),
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
                               choices = all_study_list,
                               selected = all_study_list,
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
                               choices = all_samples,
                               selected = all_samples[1:20],
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
                      shinycssloaders::withSpinner(plotlyOutput(outputId = "cellular_composition_panel.output.major_barplot",
                                   width = "100%")))),
                
                fluidRow(
                  box(title = "Cell compositions of sub-cell types", solidHeader = T, width = 12, collapsible = TRUE,
                      status = "primary",
                      shinycssloaders::withSpinner(plotlyOutput(outputId = "cellular_composition_panel.output.minor_barplot",
                                   width = "100%"))))
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
            div(DT::DTOutput("therapy_panel.FDA_approved_tbl"), style = "font-size: 90%;"), # FDA approved table
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
    
    # Shiny UI: manual page
    if(T){
      Manual_page = fluidRow(
        tags$div( style = "margin-left: 0%; margin-right: 0%; background-color: white;", 
                  tags$div(style = "margin-left: 10%; margin-right: 20%; font-size:18px;", markdown(manual_page_markdown))
                )
      )
    }
    
    # Shiny UI: case study page
    if(T){
      Case_study_page = fluidRow(
        tags$div( style = "margin-left: 0%; margin-right: 0%; background-color: white;", 
                  tags$div(style = "margin-left: 10%; margin-right: 20%; font-size:18px;", markdown(Case_study_markdown))
        )
      )
    }
    
    # Shiny UI: reference page
    if(T){
      Reference_page = fluidRow(
        tags$div( style = "margin-left: 3%; margin-right: 3%",
                  tags$h4("Browse references"),
                  DT::DTOutput("projects_info_table")
      ))
    }
    
    # Shiny UI: study tracking
    if(T){
      Study_tracking_page = fluidRow(
        tags$div( style = "margin-left: 3%; margin-right: 3%",
                  tags$h4("Browse the recently published datasets studying IBD"),
                  DT::DTOutput("study_tracking_table")
        ))
    }
    
    # Shiny UI: download page
    if(T){
      Download_page = fluidRow(
        tags$div( style = "margin-left: 3%; margin-right: 3%",
                  tags$h3("Downloads"),
                  tags$p("Data files used in this website could be accessed through the following links."),
                  tags$h4("Gene expression matrix"),
                  tags$a(href="https://figshare.com/s/e8ebff28e5cd2a5ce6b4?file=39701551", "Integrated gene expression matrix of scIBD (Anndata object in H5AD format)", target="_blank"),
                  tags$br(),
                  tags$a(href="https://figshare.com/s/e8ebff28e5cd2a5ce6b4?file=39703201", "Integrated gene expression matrix of scIBD (Seurat object in RDS format)", target="_blank"),
                  tags$br(),
                  tags$a(href="https://figshare.com/s/e8ebff28e5cd2a5ce6b4?file=39701335", "Meta data of scIBD", target="_blank"),
                  tags$br(),
                  tags$h4("Differential expressed genes"),
                  tags$a(href="https://figshare.com/s/e8ebff28e5cd2a5ce6b4?file=39701311","Differential expressed genes across cell subtypes within each major cluster", target="_blank"),
                  tags$br(),
                  tags$h4("Differential activated regulons"),
                  tags$a(href="https://figshare.com/s/e8ebff28e5cd2a5ce6b4?file=39701332", "Differential regulons between healthy and disease in each major cluster", target="_blank")
        )
      )
    }
  }
}

# shiny ui----
ui <- navbarPage(id="nav", theme = shinytheme("flatly"), collapsible = FALSE, windowTitle = "scIBD",
          # Create navbarPage----
             HTML('<a style="text-decoration:none;
                             cursor:default;
                             color:#FFFFFF;
                             font-family: Arial, Georgia, Times, serif;
                             font-weight: bold;
                             font-size: 20px;" 
                  class="active" href="#">scIBD</a>'), 
            # Home page----
             tabPanel("Home", selected = TRUE,
                      tags$div( style = "margin-left: 3%; margin-right: 3%;",
                        tags$h4("Single-cell meta-analysis of inflammatory bowel disease with scIBD"),
                        tags$p("Inflammatory bowel disease (IBD) is a type of chronic inflammation disease whose exact etiology is still unclear. \
                               With the increasing studies of IBD by single cell RNA sequencing technique (scRNA-seq), dysregulation of immune \
                               micro-environment and pathogenesis of IBD have been uncovered successively. The enormous IBD-related scRNA-seq \
                               datasets in the past decade are calling for a burning demand to be uniformly processed and integrated for conveniently access. \
                               Here, we developed a database of Single Cell transcriptomic atlas of Inflammatory Bowel Disease (scIBD) that contains ~1.14 million \
                               single cells across multiple development stages and disease states, comprising 9 major subtypes and 101 minor subtypes. \
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

            # scIBD page----
             tabPanel("Exploration", use_bs_popover(),
                        dashboardPage(
                          dashboardHeader(disable = T),
                          dashboardSidebar(
                            sidebarMenu(
                              menuItem(text = "Gene Expression Profile", tabName = "GEX_profile", selected = NULL),
                              menuItem(text = "Regulon Activity Profile", tabName = "GRN_profile", selected = NULL),
                              menuItem(text = "Gene Expression Comparison", tabName = "GEX_comparison", selected = NULL),
                              menuItem(text = "Regulon Activity Comparison", tabName = "GRN_comparison", selected = NULL),
                              menuItem(text = "Cellular composition", tabName = "cellular_composition", selected = NULL),
                              menuItem(text = "Gene Enrichment Analysis", tabName = "gsea", selected = NULL)
                            )),
                        dashboardBody(
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
                            tabItem(tabName = "gsea", gsea_panel))))
                    ),
            # Resource page----
            tabPanel("Resouces",
                       dashboardPage(
                         dashboardHeader(disable = T),
                         dashboardSidebar(
                           sidebarMenu(
                             menuItem(text = "Current Therapy Strategy", tabName = "IBD_targets", selected = TRUE),
                             menuItem(text = "GWAS-implicated Risk Genes", tabName = "IBD_risk_genes", selected = NULL),
                             menuItem(text = "Manual", tabName = "Manual", selected = NULL),
                             menuItem(text = "Case study", tabName = "Case_study", selected = NULL),
                             menuItem(text = "Reference", tabName = "Reference", selected = NULL),
                             menuItem(text = "Download", tabName = "Download", selected = NULL),
                             menuItem(text = "Study tracking", tabName = "Study_tracking", selected = NULL)
                           )),
                         dashboardBody(
                           tags$head(tags$style(HTML('.skin-blue .main-sidebar .sidebar .sidebar-menu .active a{background-color: #3079AE; color: #FFFFFF;}'))),
                           tags$head(tags$style(HTML('.skin-blue .main-sidebar {
                                          padding-top: 10px;
                                          background-color: #212F3F;
                                          font-family: "Arial", "Georgia", Times, "Times New Roman", serif;
                                          font-weight: bold;
                                          font-size: 14px;}'))),
                           tags$head(tags$style(HTML(".content-wrapper, .right-side{
                                                        background-color: #FFFFFF;}"))),
                           
                           ## content
                           tabItems(
                             tabItem(tabName = "IBD_targets", IBD_targets_panel),
                             tabItem(tabName = "IBD_risk_genes", IBD_risk_genes_panel),
                             tabItem(tabName = "Manual", Manual_page),
                             tabItem(tabName = "Case_study", Case_study_page),
                             tabItem(tabName = "Reference", Reference_page),
                             tabItem(tabName = "Download", Download_page),
                             tabItem(tabName = "Study_tracking", Study_tracking_page)
                           )))),
            # Help page----
             tabPanel("Help",
                      tags$div( style = "margin-left: 2%; margin-right: 2%",
                        tags$h3("FAQ"),
                        tags$p(tags$b("Q1: What is scIBD?")),
                        tags$p("scIBD is a platform for single-cell meta-analysis of inflammatory bowel disease (IBD) that \
                               contains ~1.14 million single cells from 12 datasets across multiple development stages (including fetal, pediatric, and adult), \
                               tissues from multiple anatomical regions (includign blood, small intestine and large intestine, etc.) and different disease states (healthy, inflammed UC, inflammed CD, etc.). \
                               scIBD comprises 9 major subtypes (Myeloid, CD4 T cells, CD8 T cells, ILCs, B/Plasma cells, Epithelial cells, Mesenchynal cells, Endothelial cells, and Neural cells), and \
                               101 cell subtypes. scIBD provides a multi-functional and user-friendly interface that provides interactive visualization for biologists to \
                               analyse the transcriptome features, gene regulatory networks and enrichment of given gene set in each cell subset."),
                        tags$p(tags$b("Q2: What are the feature functions of scIBD?")),
                        tags$p("We have integrated 12 datasets from multiple studies which investigate the pathologies of IBD, and present a comprehensive single cell transcriptomic atlas for further studying IBD. \
                               With scIBD, users are convenient to explore signature genes of each cell subtype, and compare gene expression of given genes (such as therapy targets, cytokines, IBD-GWAS related genes, or others) between health and disease across major clusters or cell subtypes. \
                               With scIBD, users are also convenient to explore the underlying gene regulatory networks (GRNs) of each cell subtype, and compare the activities of given regulons between health and disease. \
                               IBD is caused by a complex interaction between genetic and environment factors (such as gut microbes). \
                               Currently, treatments for IBD including 5-ASA, antibiotics, steroids, immunosuppressants, and biologic therapies (including anti–tumor necrosis factor [TNF] antibodies, anti–α4β7 integrin antibodies, and anti–IL12/23 antibodies). \
                               For convenience,scIBD also collected clinical trials, therapy targets, and GWAS-implicated risk genes to give a quick glance of advances in the treatment of IBD. \
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
                        tags$p("Tel: 86-10-26849285")
                      ))
)



# shiny server----
server <- function(input, output, session) {
  
  # Gene expression panel
  if(T){
    
    # initialize
    if(T){
      GEX_profile_panel.genes = reactiveValues(genes = c("TPSAB1","CPA3"))
      GEX_profile_panel.major_cluster = reactiveValues(major_cluster = "Myeloid")
      GEX_profile_panel.downsampled_cell_num = reactiveValues(num = 400)
    }
    
    # get genes 
    observeEvent(input$GEX_profile_panel.Submit, {
      GEX_profile_panel.genes$genes = input$GEX_profile_panel.Gene_list
      GEX_profile_panel.major_cluster$major_cluster = input$GEX_profile_panel.major_cluster
      GEX_profile_panel.downsampled_cell_num$num = input$GEX_profile_panel.downsampled_cell_num
    })
    
    # get cells
    if(T){
      GEX_profile_panel.cells = reactive({
        # get cells
        if(GEX_profile_panel.downsampled_cell_num$num != Inf){
          cells = tibble::rownames_to_column(adata_obs, "cell.id") %>%
            dplyr::filter(major_cluster == GEX_profile_panel.major_cluster$major_cluster) %>%
            dplyr::group_by(minor_cluster) %>% 
            dplyr::sample_n(size = min(GEX_profile_panel.downsampled_cell_num$num %>% as.numeric(), n())) %>% 
            dplyr::ungroup() %>% dplyr::select(cell.id) %>% dplyr::pull()
        }else{
          cells = adata_obs %>% dplyr::filter(major_cluster == GEX_profile_panel.major_cluster$major_cluster) %>% rownames()
        }
        # return
        list(cells = cells)
      })
    }
    
    # get umap and tsne
    if(T){
      GEX_profile_panel.scanpy_embedding = reactive({
        
        # get umap
        umap = adata_obs[GEX_profile_panel.cells()$cells, ] %>% 
          dplyr::select("UMAP_1","UMAP_2","minor_cluster") %>% as.data.frame()
        
        # get tsne
        tsne = adata_obs[GEX_profile_panel.cells()$cells, ] %>% 
          dplyr::select("TSNE_1","TSNE_2","minor_cluster") %>% as.data.frame()
        
        # return umap and tsne
        list(umap = umap, tsne = tsne)
      })
    }
    
    # get gene expression of given genes and given cells
    if(T){
      GEX_profile_panel.gene_data = reactive({
        
        # get gene expression
        expression = get_gex(adata, cells = GEX_profile_panel.cells()$cells, genes = GEX_profile_panel.genes$genes)
        
        # get average expression of input genes
        avg_exp = rowMeans(expm1(expression)) %>% log1p %>% as.data.frame
        
        # return
        list(avg_exp = avg_exp, expression = expression)
      })
    }
    
    # Scatter plot of gene expression
    if(T){
      
      GEX_profile_panel.Scanpy_embedding.plot_exp.plot = reactiveValues(plot = NULL)
      output$GEX_profile_panel.Scanpy_embedding.plot_exp <- renderPlot({
        
        # plot setting
        dot.color <- colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_profile_panel.Colorpanel,][1] %>% unlist, input$GEX_profile_panel.Colorpanel))(100)
        dot.size = input$GEX_profile_panel.Embedding_dotsize
        
        # plot
        if(input$GEX_profile_panel.Embedding_used == "UMAP"){
          # prepare data
          umap_data = cbind(GEX_profile_panel.scanpy_embedding()$umap[,c(1,2)], GEX_profile_panel.gene_data()$avg_exp)
          colnames(umap_data) = c("UMAP_1","UMAP_2","expression")
          
          # plot umap
          p= ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) + 
            geom_point(aes(color = expression), size = dot.size) + 
            scale_colour_gradientn("Exp", colors = dot.color) +
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
          tsne_data = cbind(GEX_profile_panel.scanpy_embedding()$tsne[,c(1,2)], GEX_profile_panel.gene_data()$avg_exp)
          colnames(tsne_data) = c("tSNE_1","tSNE_2","expression")
          
          # plot tsne
          p = ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2)) + 
            geom_point(aes(color = expression), size = dot.size) + 
            scale_colour_gradientn("Exp", colors = dot.color) + 
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
            scale_color_manual("Cluster", values = mycolor) +
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
                  legend.position = "right") + 
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
            scale_color_manual("Cluster", values = mycolor) +
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
          }else if(input$GEX_profile_panel.Scanpy_embedding.plot_label.file_type == 'jpeg'){
            jpeg(file, width = input$GEX_profile_panel.Scanpy_embedding.plot_label.width,
                 height = input$GEX_profile_panel.Scanpy_embedding.plot_label.height,
                 units = 'in', res = 300)
            plot(GEX_profile_panel.Scanpy_embedding.plot_label.plot$plot)
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
        dot.color <- colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_profile_panel.Colorpanel,][1] %>% unlist, input$GEX_profile_panel.Colorpanel))(100)
        dot.size = input$GEX_profile_panel.Embedding_dotsize
        
        # prepare data
        exp_in_subtype = data.frame( label = GEX_profile_panel.scanpy_embedding()$umap[,3], expression = GEX_profile_panel.gene_data()$avg_exp[,1])
        gene_mean = exp_in_subtype %>% group_by(label) %>% summarise(avg_exp = log1p(expm1(mean(expression))))
        exp_in_subtype = dplyr::left_join(exp_in_subtype, gene_mean, by = "label")
        
        # violin plot
        p = ggplot( exp_in_subtype, aes(x = label, y = expression)) + 
          geom_violin(aes(fill = avg_exp), scale = "width", color = "white", kernel = "gaussian") + 
          #geom_boxplot(outlier.size = -1, width = .1, fill = "white") + 
          scale_fill_gradientn("Exp", colors = dot.color) +
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
          }else if(input$GEX_profile_panel.Violin_plot.file_type == 'jpeg'){
            jpeg(file, width = input$GEX_profile_panel.Violin_plot.width,
                 height = input$GEX_profile_panel.Violin_plot.height,
                 units = 'in', res = 300)
            plot(GEX_profile_panel.Violin_plot.plot$plot)
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
        mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_profile_panel.Colorpanel,][1] %>% unlist, input$GEX_profile_panel.Colorpanel))(100)
        
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
          summarise(avg_exp = log1p(mean(expm1(expression))))
        
        # summarize fraction
        group_frac = exp_in_subtype %>%
          group_by(gene, label) %>%
          summarise(fraction = sum(expression > 0)/length(expression))
        
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
          filter( major_cluster == GEX_profile_panel.major_cluster$major_cluster) %>% 
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
    
    # Bar plot of cell numbers
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
        marker_gene = deg_data[[GEX_profile_panel.major_cluster$major_cluster]] %>%
          dplyr::filter(gene %in% all_genes) %>% 
          dplyr::group_by(cluster) %>% 
          dplyr::top_n(n = as.numeric(input$GEX_profile_panel.deg_topn), wt = avg_log2FC) %>% 
          dplyr::arrange(cluster, desc(avg_log2FC))
        marker_gene = marker_gene$gene
        
        # get gene expression
        expression = get_gex(adata, cells = GEX_profile_panel.cells()$cells, genes = marker_gene)
        
        # return
        list(expression = expression)
      })
    }
    
    # Heatmap plot of marker genes
    if(T){
      
      GEX_profile_panel.marker_gene.heatmap.plot = reactiveValues(plot = NULL)
      
      output$GEX_profile_panel.marker_gene.heatmap <- renderPlot({
        # set color
        mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_profile_panel.Colorpanel,][1] %>% unlist, 
                                               input$GEX_profile_panel.Colorpanel))(100)
        
        # get expression
        expression = GEX_profile_panel.top_marker_gene_expression()$expression
        
        # add cell type to expression matrix
        annotation = droplevels(GEX_profile_panel.scanpy_embedding()$umap[,3])
        exp_in_subtype = data.frame(cbind(as.character(annotation), expression))
        colnames(exp_in_subtype) = c("label",colnames(expression))
        exp_in_subtype$label = factor( exp_in_subtype$label, levels = levels(annotation) )
        
        exp_in_subtype = melt(data = data.table(exp_in_subtype), id.vars=c("label"))
        colnames(exp_in_subtype) = c("label","gene","expression")
        exp_in_subtype$expression = as.numeric(exp_in_subtype$expression)
        
        # calculate average expression of each gene
        # group by minor cluster
        group_mean = exp_in_subtype %>%
          group_by(gene, label) %>%
          summarise(avg_exp = mean(expm1(expression)))
        
        # min-max normalization
        group_mean_dcast = reshape2::dcast(group_mean, label ~ gene)
        group_mean_dcast[,2:ncol(group_mean_dcast)] = sapply(group_mean_dcast[,2:ncol(group_mean_dcast)], 
                                                             function(x) (x-mean(x))/sd(x))
        
        # format data
        group_mean = reshape2::melt(group_mean_dcast)
        colnames(group_mean) = c("label","gene","avg_exp")
        
        # create heatmap using blue color scale
        p = ggplot(group_mean, aes(label, gene)) +
          geom_tile(aes(fill = avg_exp), colour = "white") +
          scale_fill_gradientn("Exp", colors = mycolor) + 
          scale_x_discrete(limits = rev(levels(annotation))) + 
          xlab("")+
          ylab("")+
          coord_flip() + 
          theme_bw() + 
          theme(axis.text.x = element_text(size = 10, family="Times", color = "black", angle = 90, vjust = 0.5, hjust = 1),
                axis.text.y = element_text(size=10, family="Times", color = "black"),
                axis.title.x = element_text(size=12, family="Times", face= "italic", color = "black"),
                axis.title.y = element_text(size=12, family="Times", color = "black"),     
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill='transparent', color='white'),
                legend.position = "right")
        GEX_profile_panel.marker_gene.heatmap.plot$plot = p
        print(p)
      })
      
      # download
      output$GEX_profile_panel.marker_gene.heatmap.download_pdf <- downloadHandler(
        filename = function(){
          if(input$GEX_profile_panel.marker_gene.heatmap.file_type == 'pdf'){
            paste0('GEX_profile_panel.marker_gene.heatmap_', stringi::stri_rand_strings(1, 10), '.pdf')
          }else if(input$GEX_profile_panel.marker_gene.heatmap.file_type == 'jpeg'){
            paste0('GEX_profile_panel.marker_gene.heatmap_', stringi::stri_rand_strings(1, 10), '.jpeg')
          }
        },
        content = function(file){
          if(input$GEX_profile_panel.marker_gene.heatmap.file_type == 'pdf'){
            pdf(file, width = input$GEX_profile_panel.marker_gene.heatmap.width, 
                height = input$GEX_profile_panel.marker_gene.heatmap.height)
            plot(GEX_profile_panel.marker_gene.heatmap.plot$plot)
            dev.off()
          }else if(input$GEX_profile_panel.marker_gene.heatmap.file_type == 'jpeg'){
            jpeg(file, width = input$GEX_profile_panel.marker_gene.heatmap.width,
                 height = input$GEX_profile_panel.marker_gene.heatmap.height,
                 units = 'in', res = 300)
            plot(GEX_profile_panel.marker_gene.heatmap.plot$plot)
            dev.off()
          }
        }
      )
    }
    
    # display deg table
    # Marker genes of each minor cluster
    if(T){
      output$GEX_profile_panel.marker_gene_tbl = renderDT(
        datatable(deg_data[[GEX_profile_panel.major_cluster$major_cluster]],
                  filter = 'top', escape = FALSE, rownames = FALSE, options = list(pageLength = 5)) %>% formatRound(3:7, 4))
    }
  }
  
  # Gene expression comparison
  if(T){
    
    # initialize
    if(T){
      GEX_comparison_panel.genes = reactiveValues(genes = c("CPA3","TPSAB1"))
      GEX_comparison_panel.major_cluster = reactiveValues(major_cluster = "Myeloid")
      GEX_comparison_panel.minor_cluster = reactiveValues(minor_cluster = unlist(all_cluster_list))
      GEX_comparison_panel.stage = reactiveValues(stage = unlist(all_stage_list))
      GEX_comparison_panel.disease = reactiveValues(disease = unlist(all_disease_list))
      GEX_comparison_panel.tissue = reactiveValues(tissue = unlist(all_tissue_list))
      GEX_comparison_panel.tissue.sub = reactiveValues(tissue.sub = unlist(all_tissue.sub_list))
      GEX_comparison_panel.study = reactiveValues(study = unlist(all_study_list))
      GEX_comparison_panel.sample = reactiveValues(sample = all_samples)
      GEX_comparison_panel.downsampled_cell_num = reactiveValues(num = 400)
    }
    
    # observe submit
    if(T){
      observeEvent(input$GEX_comparison_panel.Submit,{
        GEX_comparison_panel.genes$genes = input$GEX_comparison_panel.Gene_list
        GEX_comparison_panel.major_cluster$major_cluster = input$GEX_comparison_panel.major_cluster
        GEX_comparison_panel.minor_cluster$minor_cluster = input$GEX_comparison_panel.minor_cluster
        GEX_comparison_panel.stage$stage = input$GEX_comparison_panel.stage
        GEX_comparison_panel.disease$disease = input$GEX_comparison_panel.disease
        GEX_comparison_panel.tissue$tissue = input$GEX_comparison_panel.tissue
        GEX_comparison_panel.tissue.sub$tissue.sub = input$GEX_comparison_panel.tissue.sub
        GEX_comparison_panel.study$study = input$GEX_comparison_panel.study
        GEX_comparison_panel.sample$sample = input$GEX_comparison_panel.sample
        GEX_comparison_panel.downsampled_cell_num$num = input$GEX_comparison_panel.downsampled_cell_num %>% as.numeric()
      }, ignoreInit = TRUE)
    }
    
    # subset data
    if(T){
      GEX_comparison_panel.subsets_data = reactive({
        # get data
        selected.adata_obs = adata_obs %>% filter(major_cluster %in% GEX_comparison_panel.major_cluster$major_cluster &
                                                    minor_cluster %in% GEX_comparison_panel.minor_cluster$minor_cluster &
                                                    stage %in% GEX_comparison_panel.stage$stage &
                                                    disease %in% GEX_comparison_panel.disease$disease &
                                                    tissue %in% GEX_comparison_panel.tissue$tissue &
                                                    tissue.sub %in% GEX_comparison_panel.tissue.sub$tissue.sub &
                                                    study %in% GEX_comparison_panel.study$study &
                                                    sample %in% GEX_comparison_panel.sample$sample)
        selected.adata_obs$minor_cluster = droplevels(selected.adata_obs$minor_cluster)
        list(selected.adata_obs = selected.adata_obs)
      })
    }
    
    # down sample and get gene expression matrix of given genes and given cells
    # get average expression
    if(T){
      GEX_comparison_panel.subsets_data.downsample = reactive({
        
        # down sample
        if(GEX_comparison_panel.downsampled_cell_num$num != Inf){
          selected_cells = tibble::rownames_to_column(GEX_comparison_panel.subsets_data()$selected.adata_obs, "barcode") %>%
            dplyr::group_by(minor_cluster) %>% 
            dplyr::sample_n(size = min(GEX_comparison_panel.downsampled_cell_num$num %>% as.numeric(), n())) %>% 
            dplyr::ungroup() %>% dplyr::select(barcode) %>% dplyr::pull()
        }else{
          selected_cells = GEX_comparison_panel.subsets_data()$selected.adata_obs %>% rownames()
        }
        
        # get expression
        expression = get_gex(adata, cells = selected_cells, genes = GEX_comparison_panel.genes$genes)
        
        # get subset meta data
        selected_meta_data = adata_obs[selected_cells, ]
        
        # get average expression of given input genes
        avg_exp = log1p(rowMeans(expm1(expression))) %>% as.data.frame()
        
        # combine result
        result = cbind(selected_meta_data, avg_exp[,1])
        colnames(result) = c(colnames(selected_meta_data),"avg_exp")
        
        # update result
        result$major_cluster = droplevels(result$major_cluster)
        result$minor_cluster = droplevels(result$minor_cluster)
        result$stage = droplevels(result$stage)
        result$disease = droplevels(result$disease)
        result$tissue = droplevels(result$tissue)
        result$tissue.sub = droplevels(result$tissue.sub)
        result$study = droplevels(result$study)
        result$sample = droplevels(result$sample)
        
        # clean
        rm(selected_cells, expression, selected_meta_data, avg_exp)
        
        # return
        list(result = result) })
    }
    
    # subset scIBD
    if(T){
      GEX_comparison_panel.subsets_scIBD = reactive({
        if(GEX_comparison_panel.downsampled_cell_num$num != Inf){
          selected_cells = tibble::rownames_to_column(adata_obs, "barcode") %>%
            dplyr::group_by(minor_cluster) %>% 
            dplyr::sample_n(size = min(GEX_comparison_panel.downsampled_cell_num$num %>% as.numeric(), n())) %>% 
            dplyr::ungroup() %>% dplyr::select(barcode) %>% dplyr::pull()
        }else{
          selected_cells = adata_obs %>% rownames()
        }
        
        # get expression
        expression = get_gex(adata, cells = selected_cells, genes = GEX_comparison_panel.genes$genes)
        
        # get subset meta data
        selected_meta_data = adata_obs[selected_cells, ]
        
        # get average expression of given input genes
        avg_exp = log1p(rowMeans(expm1(expression))) %>% as.data.frame()
        
        # combine result
        result = cbind(selected_meta_data, avg_exp[,1])
        colnames(result) = c(colnames(selected_meta_data),"avg_exp")
        
        list(result = result)
      })
    }
    
    # plot global scanpy embedding (left)
    # down sample scIBD
    # color by average expression of given genes
    if(T){
      GEX_comparison_panel.Scanpy_embedding.plot_exp.plot = reactiveValues(plot = NULL)
      output$GEX_comparison_panel.Scanpy_embedding.plot_exp <- renderPlot({
        
        # plot setting
        dot.color <- colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_comparison_panel.Colorpanel,][1] %>% unlist, input$GEX_comparison_panel.Colorpanel))(100)
        dot.size = input$GEX_comparison_panel.Embedding_dotsize
        
        # plot gene expression in scanpy embedding
        if(input$GEX_comparison_panel.Embedding_used == "UMAP"){
          
          # plot umap
          p = ggplot(GEX_comparison_panel.subsets_scIBD()$result[,c("gUMAP_1","gUMAP_2","avg_exp")], aes(x = gUMAP_1, y = gUMAP_2)) +
            geom_point(aes(color = avg_exp), size = dot.size) + 
            scale_colour_gradientn("Exp", colors = dot.color) +
            xlab("UMAP_1") +
            ylab("UMAP_2") +
            xlim( min(adata_obs$gUMAP_1), max(adata_obs$gUMAP_1) ) + 
            ylim( min(adata_obs$gUMAP_2), max(adata_obs$gUMAP_2) ) + 
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, family = "Times",color = "black"),
                  axis.text.y = element_text(size=10, family = "Times", color = "black"),
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
          
          # plot tsne
          p = ggplot(GEX_comparison_panel.subsets_scIBD()$result[,c("gTSNE_1","gTSNE_2","avg_exp")], aes(x = gTSNE_1, y = gTSNE_2)) +
            geom_point(aes(color = avg_exp), size = dot.size) + 
            scale_colour_gradientn("Exp", colors = dot.color) + 
            xlab("tSNE_1") + 
            ylab("tSNE_2") + 
            xlim( min(adata_obs$gTSNE_1), max(adata_obs$gTSNE_2) ) + 
            ylim( min(adata_obs$gTSNE_1), max(adata_obs$gTSNE_2) ) + 
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, family = "Times",color = "black"),
                  axis.text.y = element_text(size=10, family = "Times", color = "black"),
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
    # down sample scIBD
    # color by major cell type
    if(T){
      GEX_comparison_panel.Scanpy_embedding.plot_label.plot = reactiveValues(plot = NULL)
      output$GEX_comparison_panel.Scanpy_embedding.plot_label <- renderPlot({
        
        # plot setting
        dot.size = input$GEX_comparison_panel.Embedding_dotsize
        
        # plot cell subset annotation in Scanpy embedding
        if(input$GEX_comparison_panel.Embedding_used == "UMAP"){
          
          # plot umap
          p = ggplot(GEX_comparison_panel.subsets_scIBD()$result, aes(x = gUMAP_1, y = gUMAP_2, color = major_cluster)) +
            geom_point(size = dot.size) +
            scale_color_manual("Label", values = major_cluster_color ) + 
            xlab("UMAP_1") +
            ylab("UMAP_2") + 
            xlim( min(global_umap$gUMAP_1), max(global_umap$gUMAP_1) ) + 
            ylim( min(global_umap$gUMAP_2), max(global_umap$gUMAP_2) ) + 
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
          
          # plot tsne
          p = ggplot(GEX_comparison_panel.subsets_scIBD()$result, aes(x = gTSNE_1, y = gTSNE_2, color = major_cluster)) +
            geom_point(size = dot.size) + 
            scale_color_manual("Label", values = major_cluster_color ) +
            xlab("tSNE_1") +
            ylab("tSNE_2") + 
            xlim( min(global_tsne$gTSNE_1), max(global_tsne$gTSNE_1) ) + 
            ylim( min(global_tsne$gTSNE_2), max(global_tsne$gTSNE_2) ) + 
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
        }})
      
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
    
    # plot global scanpy embedding (right)
    # selected cells
    # colored by user defined values (major/minor/stage/dev/disease)
    if(T){
      GEX_comparison_panel.Scanpy_embedding.plot_cell.plot = reactiveValues(plot = NULL)
      output$GEX_comparison_panel.Scanpy_embedding.plot_cell <- renderPlot({
        
        # plot setting
        dot.size = input$GEX_comparison_panel.Embedding_dotsize
        
        # plot cell in Scanpy embedding
        if(input$GEX_comparison_panel.Embedding_used == "UMAP"){
          
          # prepare data
          umap_data = GEX_comparison_panel.subsets_data.downsample()$result[,c("gUMAP_1","gUMAP_2",input$GEX_comparison_panel.groupby)]
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
          tsne_data = GEX_comparison_panel.subsets_data.downsample()$result[,c("gTSNE_1","gTSNE_2",input$GEX_comparison_panel.groupby)]
          colnames(tsne_data) = c("tSNE_1","tSNE_2","label")
          
          # set color
          mycolor = ""
          if(input$GEX_comparison_panel.groupby == "major_cluster"){
            mycolor = major_cluster_color
          }else if(input$GEX_comparison_panel.groupby == "minor_cluster"){
            mycolor = minor_cluster_color[tsne_data$label %>% droplevels %>% levels]
          }else{
            colourCount = tsne_data$label %>% droplevels %>% levels %>% length
            getPalette = colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_comparison_panel.Colorpanel,'maxcolors'], input$GEX_comparison_panel.Colorpanel))
            mycolor = getPalette(colourCount)
          }
          
          # plot tsne
          p = ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2, color = label)) +
            geom_point(size = dot.size) + 
            scale_color_manual("Label", values = mycolor) +
            xlab("tSNE_1") +
            ylab("tSNE_2") + 
            xlim( min(adata_obs$TSNE_1), max(adata_obs$TSNE_1) ) + 
            ylim( min(adata_obs$TSNE_2), max(adata_obs$TSNE_2) ) + 
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
        mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_comparison_panel.Colorpanel,][1] %>% unlist, input$GEX_comparison_panel.Colorpanel))(100)
        
        # prepare data
        exp_in_subtype = GEX_comparison_panel.subsets_data.downsample()$result[,c(input$GEX_comparison_panel.groupby, "avg_exp")]
        colnames(exp_in_subtype) = c("group","avg_exp")
        
        # get avg of avg by group
        group_avg_exp = exp_in_subtype %>% dplyr::group_by(group) %>% dplyr::summarise(man_of_avg = log1p(mean(expm1(avg_exp))))
        exp_in_subtype = dplyr::left_join(exp_in_subtype, group_avg_exp, by = "group")
        
        # violin plot
        p = ggplot( exp_in_subtype, aes(x = group, y = avg_exp, fill = group, color = group)) + 
          geom_violin(aes(fill = man_of_avg), scale = "width", color = "white", kernel = "gaussian") + 
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
        print(p)})
      
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
        cts = GEX_comparison_panel.subsets_data()$selected.adata_obs %>% group_by(across(input$GEX_comparison_panel.groupby)) %>%
          summarise(cts = n()) %>% as.data.frame()
        colnames(cts) = c("group","count")
        
        pct = GEX_comparison_panel.subsets_data()$selected.adata_obs %>% group_by(across(input$GEX_comparison_panel.groupby)) %>%
          summarise(pct = n()/nrow(GEX_comparison_panel.subsets_data()$selected.adata_obs)) %>% as.data.frame()
        colnames(pct) = c("group","pct")
        
        if(input$GEX_comparison_panel.Cell_number.Bar_plot.type == 'Number'){
          
          # bar plot
          p = ggplot( cts, aes(x = group, y = count, fill = group)) + 
            geom_bar(stat = "identity") + 
            ylab("Number of cells") + 
            xlab("")+
            theme_bw() +
            theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1, family = "Times"),
                  axis.text.y = element_text(size=10,  color = "black", family = "Times"),
                  axis.title.x = element_text(size=12, color = "black", family = "Times"),
                  axis.title.y = element_text(size=12,  color = "black", family = "Times"),     
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill='transparent', color='white'),
                  panel.border = element_blank(),
                  legend.position = "right") +
            theme(axis.line = element_line(color = "black"),
                  legend.position = "none")
          
        }else if(input$GEX_comparison_panel.Cell_number.Bar_plot.type == 'Percentage'){
          
          # bar plot
          p = ggplot( pct, aes(x = group, y = pct, fill = group)) + 
            geom_bar(stat = "identity") +
            ylab("Percentage of cells") + 
            xlab("")+
            theme_bw() +
            theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1, family = "Times"),
                  axis.text.y = element_text(size=10,  color = "black", family = "Times"),
                  axis.title.x = element_text(size=12, color = "black", family = "Times"),
                  axis.title.y = element_text(size=12,  color = "black", family = "Times"),     
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
        mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GEX_comparison_panel.Colorpanel,][1] %>% unlist, input$GEX_comparison_panel.Colorpanel))(100)
        
        # prepare data
        exp_in_subtype = GEX_comparison_panel.subsets_data.downsample()$result[,c(input$GEX_comparison_panel.groupby.1, input$GEX_comparison_panel.groupby.2, "avg_exp")]
        colnames(exp_in_subtype) = c("group1","group2","avg_exp")
        
        # get avg of avg by group
        group_avg_exp = exp_in_subtype %>% dplyr::group_by(group1, group2) %>% dplyr::summarise(man_of_avg = log1p(mean(expm1(avg_exp))))
        exp_in_subtype = dplyr::left_join(exp_in_subtype, group_avg_exp, by = c("group1","group2"))
        
        # violin plot
        p = ggplot( exp_in_subtype, aes(x = group2, y = avg_exp)) + 
          geom_violin(aes(fill = man_of_avg), scale = "width", color = "white", kernel = "gaussian") + 
          scale_fill_gradientn("Exp", colors = mycolor) +
          ylab("Relative gene expression") + 
          xlab("")+
          facet_grid(cols = vars(group1)) + 
          theme_bw() + 
          theme(axis.text.x = element_text(size = 10,  color = "black", angle = 90, vjust = 0.5, hjust = 1, family = "Times"),
                axis.text.y = element_text(size=10,  color = "black", family = "Times"),
                axis.title.x = element_text(size=12, color = "black", family = "Times"),
                axis.title.y = element_text(size=12,  color = "black", family = "Times"),     
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
    
    # observe groups
    if(T){
      observe({
        # set faceting_group
        updateSelectInput(
          session = session,
          inputId = "GEX_comparison_panel.Bar_plot2.faceting_group",
          choices = c(input$GEX_comparison_panel.groupby.1, input$GEX_comparison_panel.groupby.2),
          selected = input$GEX_comparison_panel.groupby.1,
        )})
    }
    
    # plot cell number of selected cells
    # bar plot
    # group by two variables
    if(T){
      GEX_comparison_panel.Cell_number.Bar_plot2.plot = reactiveValues(plot = NULL)
      output$GEX_comparison_panel.Cell_number.Bar_plot2 <- renderPlot({
        
        # prepare data
        selected_cells = GEX_comparison_panel.subsets_data()$selected.adata_obs[,c(input$GEX_comparison_panel.groupby.1, input$GEX_comparison_panel.groupby.2)]
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
    
    # initialize
    if(T){
      GRN_profile_panel.tfs = reactiveValues(tfs = c("RUNX2"))
      GRN_profile_panel.major_cluster = reactiveValues(major_cluster = c("Myeloid"))
      GRN_profile_panel.downsampled_cell_num = reactiveValues(num = 400)
    }
    
    # check genes
    observe({
      if( length(input$GRN_profile_panel.tfs) == 0){
        message = paste0("Please select a gene name!")
        sendSweetAlert(session = session, title = message, type = "warning")
      }
    })
    
    # observe submit
    if(T){
      observeEvent(input$GRN_profile_panel.Submit,{
        GRN_profile_panel.tfs$tfs = input$GRN_profile_panel.tfs
        GRN_profile_panel.major_cluster$major_cluster = input$GRN_profile_panel.major_cluster
        GRN_profile_panel.downsampled_cell_num$num = input$GRN_profile_panel.downsampled_cell_num
      }, ignoreInit = TRUE)
    }
    
    # observe major cluster
    if(T){
      observeEvent(input$GRN_profile_panel.major_cluster, {
        updatePickerInput(
          session = session,
          inputId = "GRN_profile_panel.tfs",
          choices = tf_list[[ input$GRN_profile_panel.major_cluster ]],
          selected = tf_list[[ input$GRN_profile_panel.major_cluster ]][1],
          choicesOpt = list(style = "background:white; color:black")
        )
      }, ignoreInit = TRUE)
    }
    
    # get cells
    if(T){
      GRN_profile_panel.cells = reactive({
        # get cells
        if(GRN_profile_panel.downsampled_cell_num$num == Inf){
          cells = get_cell_annotation(loom_data[[GRN_profile_panel.major_cluster$major_cluster]]) %>% rownames()
        }else{
          cells = tibble::rownames_to_column(get_cell_annotation(loom_data[[GRN_profile_panel.major_cluster$major_cluster]]), "barcode") %>%
            dplyr::group_by(Label) %>% 
            dplyr::sample_n(size = min(GRN_profile_panel.downsampled_cell_num$num %>% as.numeric(), n())) %>% 
            dplyr::ungroup() %>% dplyr::select(barcode) %>% dplyr::pull()
        }
        list(cells = cells)
      })
    }
    
    # get data
    if(T){
      GRN_profile_panel.scenic_activity = reactive({
        
        # get cells
        cells = GRN_profile_panel.cells()$cells
        
        # get meta data
        cellInfo = get_cell_annotation(loom_data[[GRN_profile_panel.major_cluster$major_cluster]])[cells, ]
        cellInfo$barcode = rownames(cellInfo)
        
        # get tfs
        genes = GRN_profile_panel.tfs$tfs
        
        # get regulon activity
        regulonsAUC = get_regulons_AUC( loom_data[[GRN_profile_panel.major_cluster$major_cluster]], column.attr.name = "RegulonsAUC")[genes, cells]
        grn_activity = getAUC(regulonsAUC) %>% t() %>% as.data.frame()
        grn_activity$barcode = rownames(grn_activity)
        grn_activity$avg_grn_activity = rowMeans(grn_activity[,genes, drop = FALSE])
        
        # combine
        result = dplyr::left_join(grn_activity, cellInfo[,c("barcode","Label")], by = "barcode") %>% as.data.frame()
        
        # clean
        rm(cells, cellInfo, genes, grn_activity)
        
        # return
        list(result = result)
      })
    }
    
    # get cell information
    if(T) { 
      GRN_profile_panel.cellInfo = reactive({
        # return
        list(cellInfo = get_cell_annotation(loom_data[[GRN_profile_panel.major_cluster$major_cluster]])[GRN_profile_panel.cells()$cells,] )
      })
    }
    
    # get scanpy embedding
    if(T){
      GRN_profile_panel.scanpy_embedding <- reactive({
        # return umap and tsne
        list(umap = get_embeddings(loom_data[[GRN_profile_panel.major_cluster$major_cluster]])$Scanpy_UMAP, 
             tsne = get_embeddings(loom_data[[GRN_profile_panel.major_cluster$major_cluster]])$Scanpy_tSNE)
      })
    }
    
    # get scenic embedding
    if(T){
      GRN_profile_panel.scenic_embedding <- reactive({
        # return umap and tsne
        list(umap = get_embeddings(loom_data[[GRN_profile_panel.major_cluster$major_cluster]])$SCENIC_AUC_UMAP,
             tsne = get_embeddings(loom_data[[GRN_profile_panel.major_cluster$major_cluster]])$SCENIC_AUC_tSNE)
      })
    }
    
    # get grn network
    if(T){
      GRN_profile_panel.scenic_network <- reactive({
        
        # get network
        regulon_name = paste0(GRN_profile_panel.tfs$tfs, "(+)")
        regulons_incidMat <- get_regulons(loom_data[[GRN_profile_panel.major_cluster$major_cluster]], column.attr.name = "Regulons")[regulon_name,,drop=FALSE] 
        grn_network <- regulonsToGeneLists(regulons_incidMat)
        names(grn_network) = GRN_profile_panel.tfs$tfs
        
        # exclude self
        for(i in 1:length(grn_network)){
          grn_network[[i]] = grn_network[[i]][ !(grn_network[[i]] %in% names(grn_network)[i]) ]
        }
        
        # get number of targets
        targets_num = c()
        grn_network.filter = list()
        if(input$GRN_profile_panel.targets_num != Inf){
          df = data.frame(total = lapply(grn_network, length) %>% unlist(), max = input$GRN_profile_panel.targets_num %>% as.integer())
          targets_num = apply(df, 1, FUN = min) %>% unlist()
        }else{
          targets_num = lapply(grn_network, length) %>% unlist()
        }
        names(targets_num) = GRN_profile_panel.tfs$tfs
        
        # get correlation between TF and downstream genes
        tf_gene_cor_table = data.frame(TF = NULL, target = NULL, correlation = NULL, p.value = NULL, expression = NULL)
        for(i in 1:length(grn_network)){
          tf = names(grn_network)[i]
          targets = grn_network[[ tf ]]
          gex = get_gex(adata, cells = GRN_profile_panel.cells()$cells, genes = c(tf, targets))
          
          # the first is tf, the rest is targets
          for(j in 2:length(gex)){
            # calculate correlation
            tmp = gex[,c(1,j)]
            tmp = tmp[tmp[,1] > 0 & tmp[,2] > 0,]
            tmp_avg_exp = log1p(mean(expm1(gex[,j]))) %>% round(digits = 4)
            if(nrow(tmp) >= 3){
              spearman_cor = cor(tmp[,1], tmp[,2], method = "spearman") %>% round(digits = 4)
              p.value = cor.test(tmp[,1], tmp[,2], method = "spearman", exact = FALSE)$p.value %>% 
                format(digits = 4, scientific = TRUE)
              tf_gene_cor_table = rbind(tf_gene_cor_table, c(tf, colnames(gex)[j], spearman_cor, p.value, tmp_avg_exp))
            }else{
              tf_gene_cor_table = rbind(tf_gene_cor_table, c(tf, colnames(gex)[j], NA, NA, tmp_avg_exp))
            }
          }
          
          # filter grn network
          targets.filter = expm1(gex[,-1]) %>% colMeans() %>% sort(decreasing = TRUE) %>% head(n=targets_num[tf]) %>% names()
          grn_network.filter[[tf]] = targets.filter
        }
        colnames(tf_gene_cor_table) = c("TF","target","correlation","p-value","expression")
        
        # clean
        rm(tf, targets, gex)
        
        # get nodes and links of the GRN networks
        grn_network = grn_network.filter
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
        
        ## return
        list(scenic_nodes = scenic_nodes, 
             scenic_links = scenic_links,
             tf_gene_cor_table = tf_gene_cor_table)
      })
    }
    
    # display regulon activity in Scanpy embedding
    # down sampled cells
    if(T){
      GRN_profile_panel.Scanpy_embedding.plot_activity.plot = reactiveValues(plot = NULL)
      output$GRN_profile_panel.Scanpy_embedding.plot_activity <- renderPlot({
        
        if(input$GRN_profile_panel.Embedding_used == "UMAP"){
          
          # prepare data
          umap_data = GRN_profile_panel.scanpy_embedding()$umap[GRN_profile_panel.scenic_activity()$result$barcode, ] %>% as.data.frame()
          umap_data = cbind(umap_data, activity = GRN_profile_panel.scenic_activity()$result$avg_grn_activity)
          colnames(umap_data) = c("UMAP_1","UMAP_2","activity")
          
          # plot tsne
          mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GRN_profile_panel.Colorpanel,][1] %>% unlist, input$GRN_profile_panel.Colorpanel))(100)
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
          
          # prepare data
          tnse_data = GRN_profile_panel.scanpy_embedding()$tsne[GRN_profile_panel.scenic_activity()$result$barcode, ] %>% as.data.frame()
          tnse_data = cbind(tnse_data, activity = GRN_profile_panel.scenic_activity()$result$avg_grn_activity)
          colnames(tnse_data) = c("TSNE_1","TSNE_2","activity")
          
          # plot umap
          mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GRN_profile_panel.Colorpanel,][1] %>% unlist, input$GRN_profile_panel.Colorpanel))(100)
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
    # down sampled cells
    if(T){
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
                                   levels = all_cluster_list[[GRN_profile_panel.major_cluster$major_cluster]])
          
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
                                   levels = all_cluster_list[[GRN_profile_panel.major_cluster$major_cluster]])
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
    # down sampled cells
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
                                   levels = all_cluster_list[[GRN_profile_panel.major_cluster$major_cluster]])
          
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
                                   levels = all_cluster_list[[GRN_profile_panel.major_cluster$major_cluster]])
          
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
        activity_in_subtype = GRN_profile_panel.scenic_activity()$result[,c("Label","avg_grn_activity")]
        colnames(activity_in_subtype) = c("label","activity")
        
        # get the average activity of multiple cell types
        group_mean = activity_in_subtype %>% group_by(label) %>% summarise(avg_activity = mean(activity))
        rnames = rownames(activity_in_subtype)
        
        # combine
        activity_in_subtype = activity_in_subtype %>% dplyr::left_join(group_mean, by = "label")
        
        # cell type order
        cell_type_order = adata_obs[adata_obs$major_cluster == GRN_profile_panel.major_cluster$major_cluster, ]$minor_cluster %>% droplevels() %>% levels
        
        # violin plot
        mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GRN_profile_panel.Colorpanel,][1] %>% unlist, input$GRN_profile_panel.Colorpanel))(100)
        p = ggplot( activity_in_subtype, aes(x = label, y = activity)) + 
          geom_violin(aes(fill = avg_activity), scale = "width", color = "white", kernel = "gaussian") + 
          scale_fill_gradientn("Activity", colors = mycolor) +
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
        
        # prepare data
        activity_in_subtype = GRN_profile_panel.scenic_activity()$result
        activity_in_subtype = activity_in_subtype[,-c(ncol(activity_in_subtype)-1, ncol(activity_in_subtype)-2)]
        colnames(activity_in_subtype)[ncol(activity_in_subtype)] = 'label'
        
        activity_in_subtype = reshape2::melt(data = activity_in_subtype, id.vars = c("label") )
        colnames(activity_in_subtype) = c("label","gene","activity")
        
        ## summarize data
        activity_mean = activity_in_subtype %>%
          dplyr::group_by(gene, label) %>%
          dplyr::summarise(avg_activity = mean(activity))
        
        # cell type order
        cell_type_order = adata_obs[adata_obs$major_cluster == GRN_profile_panel.major_cluster$major_cluster, ]$minor_cluster %>% droplevels() %>% levels
        
        # set color
        mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GRN_profile_panel.Colorpanel,][1] %>% unlist, input$GRN_profile_panel.Colorpanel))(100)
        
        # plot heatmap
        p = ggplot(activity_mean, aes(label, gene)) +
          geom_tile(aes(fill = avg_activity), colour = "white") +
          scale_fill_gradientn("Activity", colors = mycolor) + 
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
        renderDT({
          datatable(GRN_profile_panel.scenic_network()$tf_gene_cor_table, filter = "top", escape = FALSE, options = list(pageLength = 5))
        })
      })
    }
    
    # display rss data
    if(T){
      output$GRN_profile_panel.rss_tbl <- renderDT({
        data = rss_data[[GRN_profile_panel.major_cluster$major_cluster]]
        datatable(data, options = list(scrollX=TRUE, pageLength = 5) ) %>%
          formatRound(columns = c(1:ncol(data)), digits = 2) %>%
          formatStyle(columns = c(1:ncol(data)), 'text-align' = 'middle')})
    }
  }
  
  # Regulon activity comparison
  if(T){
    
    # initialize
    if(T){
      GRN_comparison_panel.major_cluster = reactiveValues(major_cluster = "Myeloid")
      GRN_comparison_panel.minor_cluster = reactiveValues(minor_cluster = all_cluster_list[["Myeloid"]])
      GRN_comparison_panel.TF = reactiveValues(tf = "ATF3")
      # select all
      GRN_comparison_panel.tissue = reactiveValues(tissue = all_tissue_list)
      GRN_comparison_panel.tissue.sub = reactiveValues(tissue.sub = unlist(all_tissue.sub_list))
      GRN_comparison_panel.stage = reactiveValues(stage = all_stage_list)
      GRN_comparison_panel.disease = reactiveValues(disease = all_disease_list)
      GRN_comparison_panel.study = reactiveValues(study = all_study_list)
      GRN_comparison_panel.sample = reactiveValues(sample = all_samples)
      GRN_comparison_panel.downsampled_cell_num = reactiveValues(num = 400)
      
    }
    
    # observe major cluster
    if(T){
      observeEvent(input$GRN_comparison_panel.major_cluster, {
        updatePickerInput(
          session = session,
          inputId = "GRN_comparison_panel.TF",
          choices = tf_list[[ input$GRN_comparison_panel.major_cluster ]],
          selected = "ATF3",
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
    }
    
    # check genes
    if(T){
      observe({
        if( length(input$GRN_comparison_panel.TF) == 0){
          message = paste0("Please select a gene name!")
          sendSweetAlert(session = session, title = message, type = "warning")
        }
      })
    }
    
    # observe submit
    if(T){
      observeEvent(input$GRN_comparison_panel.Submit,{
        GRN_comparison_panel.major_cluster$major_cluster = input$GRN_comparison_panel.major_cluster
        GRN_comparison_panel.minor_cluster$minor_cluster = input$GRN_comparison_panel.minor_cluster
        GRN_comparison_panel.TF$tf = input$GRN_comparison_panel.TF
        GRN_comparison_panel.tissue$tissue = input$GRN_comparison_panel.tissue
        GRN_comparison_panel.tissue.sub$tissue.sub = input$GRN_comparison_panel.tissue.sub
        GRN_comparison_panel.stage$stage = input$GRN_comparison_panel.stage
        GRN_comparison_panel.disease$disease = input$GRN_comparison_panel.disease
        GRN_comparison_panel.study$study = input$GRN_comparison_panel.study
        GRN_comparison_panel.sample$sample = input$GRN_comparison_panel.sample
        GRN_comparison_panel.downsampled_cell_num$num = input$GRN_comparison_panel.downsampled_cell_num
      }, ignoreInit = TRUE)
    }
    
    # select data
    if(T){
      GRN_comparison_panel.subsets_data = reactive({
        # select data
        selected.adata_obs = adata_obs %>% filter(major_cluster == GRN_comparison_panel.major_cluster$major_cluster &
                                                    minor_cluster %in% GRN_comparison_panel.minor_cluster$minor_cluster &
                                                    stage %in% GRN_comparison_panel.stage$stage &
                                                    disease %in% GRN_comparison_panel.disease$disease &
                                                    tissue %in% GRN_comparison_panel.tissue$tissue &
                                                    tissue.sub %in% GRN_comparison_panel.tissue.sub$tissue.sub &
                                                    study %in% GRN_comparison_panel.study$study &
                                                    sample %in% GRN_comparison_panel.sample$sample)
        
        # get cells
        # down sample
        selected.adata_obs$minor_cluster = droplevels(selected.adata_obs$minor_cluster)
        list(selected.adata_obs = selected.adata_obs)
      })
    }
    
    # down sample and get regulon activity of given tf
    if(T){
      GRN_comparison_panel.subsets_data.downsample = reactive({
        
        if(GRN_comparison_panel.downsampled_cell_num$num != Inf){
          selected_cells = tibble::rownames_to_column(GRN_comparison_panel.subsets_data()$selected.adata_obs, "barcode") %>%
            dplyr::filter(major_cluster == GRN_comparison_panel.major_cluster$major_cluster) %>%
            dplyr::group_by(minor_cluster) %>% 
            dplyr::sample_n(size = min(GRN_comparison_panel.downsampled_cell_num$num %>% as.numeric(), n())) %>% 
            dplyr::ungroup() %>% dplyr::select(barcode) %>% dplyr::pull()
        }else{
          selected_cells = GRN_comparison_panel.subsets_data()$selected.adata_obs %>% rownames()
        }
        
        # get subset meta data
        selected_meta_data = GRN_comparison_panel.subsets_data()$selected.adata_obs[selected_cells, ]
        
        # get regulon activity
        regulonsAUC = get_regulons_AUC( loom_data[[GRN_comparison_panel.major_cluster$major_cluster]], column.attr.name = "RegulonsAUC")[GRN_comparison_panel.TF$tf, selected_cells]
        grn_activity = getAUC(regulonsAUC) %>% t() %>% as.data.frame()
        grn_activity$barcode = rownames(grn_activity)
        
        # update selected_meta_data
        selected_meta_data$stage = droplevels(selected_meta_data$stage)
        selected_meta_data$disease = droplevels(selected_meta_data$disease)
        selected_meta_data$tissue = droplevels(selected_meta_data$tissue)
        selected_meta_data$tissue.sub = droplevels(selected_meta_data$tissue.sub)
        selected_meta_data$study = droplevels(selected_meta_data$study)
        selected_meta_data$sample = droplevels(selected_meta_data$sample)
        
        # combine
        selected_meta_data$barcode = rownames(selected_meta_data)
        result = dplyr::left_join(selected_meta_data, grn_activity, by = "barcode")
        
        # clean
        rm(selected_cells, regulonsAUC, grn_activity, selected_meta_data)
        
        # return
        list(result = result)
      })
    }
    
    # plot regulon activity of selected cells
    # violin plot
    # group by one variable, major/minor, stage/disease/tissue
    if(T){
      GRN_comparison_panel.Violin_plot.plot = reactiveValues(plot = NULL)
      output$GRN_comparison_panel.Violin_plot <- renderPlot({
        
        # set color
        mycolor <- colorRampPalette(brewer.pal(brewer.pal.info[input$GRN_comparison_panel.Colorpanel,][1] %>% unlist, input$GRN_comparison_panel.Colorpanel))(100)
        
        # prepare data
        # get label
        label = GRN_comparison_panel.subsets_data.downsample()$result[, input$GRN_comparison_panel.groupby]
        activity = GRN_comparison_panel.subsets_data.downsample()$result[,GRN_comparison_panel.TF$tf]
        
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
        
        # sleep 0.25s
        Sys.sleep("0.25")
        
        # prepare data
        selected_cells = GRN_comparison_panel.subsets_data()$selected.adata_obs %>%
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
      output$GRN_comparison_panel.diff_regulon.healthy_UC.tbl <- renderDT(
        datatable(diff_regulon_data[[input$GRN_comparison_panel.major_cluster]][[1]],
                  filter = "top", escape = FALSE, rownames = FALSE, options = list(pageLength = 5)) %>% formatRound(2:6, 6)
      )
    }
    
    # diff regulon, CD vs healthy
    if(T){
      output$GRN_comparison_panel.diff_regulon.healthy_CD.tbl <- renderDT(
        datatable(diff_regulon_data[[input$GRN_comparison_panel.major_cluster]][[2]],
                  filter = "top", escape = FALSE, rownames = FALSE, options = list(pageLength = 5)) %>% formatRound(2:6, 6)
      )
    }
    
  }
  
  # Cellular composition
  if(T){
    
    # initialize
    if(T){
      # default selection
      cellular_composition_panel.samples = reactiveValues(samples = all_samples[1:20])
      cellular_composition_panel.studies = reactiveValues(studies = all_study_list)
    }
    
    # update plot only when click submit
    if(T){
      observeEvent(input$cellular_composition_panel.Submit,{
        cellular_composition_panel.samples$samples = input$cellular_composition_panel.sample
        cellular_composition_panel.studies$studies = input$cellular_composition_panel.study
      })
    }
    
    # select samples
    if(T){
      cellular_composition_panel.select_meta_data = reactive({
        select_meta_data = adata_obs %>% filter(study %in% cellular_composition_panel.studies$studies &
                                                  sample %in% cellular_composition_panel.samples$samples)
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
  
  # Gene set enrichment analysis
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
          output$Gsea_panel.output.tbl <- renderDT(
            datatable(gsea.result, escape = FALSE, rownames = F)
          )
        }
        
        # display heat map
        if(T){
          Gsea_panel.output.heatmap.plot = reactiveValues(plot = NULL)
          
          # down sample cells
          cells = tibble::rownames_to_column(adata_obs, "cell.id") %>%
            dplyr::group_by(minor_cluster) %>% 
            dplyr::sample_n(size = min(500, n())) %>% 
            dplyr::ungroup() %>% dplyr::select(cell.id) %>% dplyr::pull()
          
          output$Gsea_panel.output.heatmap = renderPlot({
            p = plot_heatmap.minor_cluster(adata, cells = cells, genes = genes, mj_color = major_cluster_color)
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
  
  # Resource panel
  if(T){
    
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
      output$therapy_panel.FDA_approved_tbl = renderDT(
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
          renderDT(detail, escape = FALSE, rownames = FALSE,
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
          renderDT(detail, escape = FALSE, rownames = TRUE,
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
          renderDT(detail, escape = FALSE, rownames = FALSE, 
                   options = list(ordering=TRUE))},
          title = paste0("Clinical trails for ", input$select_button[1]),
          size = "l",
          footer = modalButton("Close")))
      })
      
      # Drugs and targets for IBD
      output$therapy_panel.ibd_targets_tbl = renderDT( 
        DT::datatable(therapy_panel.drug_and_target, 
                      filter = 'top', escape = FALSE, rownames = FALSE,
                      options = list(pageLength = 5)))
    }
    
    # Risk gene panel
    if(T){ 
      # from GWAS research
      output$ibd_gwas_tbl = renderDT(
        DT::datatable(ibd_gwas_data, 
                      filter = "top", escape = FALSE, rownames = FALSE,
                      options = list(pageLength = 5, deferRender = TRUE)))
      
      # adult
      output$adult_ibd_risk_genes_tbl = renderDT(
        DT::datatable(adult_ibd_risk_genes_data,
                      filter = "top", escape = FALSE, rownames = FALSE,
                      options = list(pageLength = 5, deferRender = TRUE)))
      
      # pediatric
      output$pediatric_ibd_risk_genes_tbl = renderDT(
        DT::datatable(pediatric_ibd_risk_genes_data, rownames = FALSE,
                      options = list(pageLength = 5, deferRender = TRUE)))
    }
    
    # Reference panel
    if(T){
      output$projects_info_table = renderDT(
        DT::datatable(project_info, escape = FALSE, rownames = FALSE,
                      options = list(pageLength = 5, deferRender = TRUE)))
    }
    
    # Study tracking panel
    if(T){
      output$study_tracking_table = renderDT(
        DT::datatable(study_tracking_info, escape = FALSE, rownames = FALSE,
                      options = list(pageLength = 5, deferRender = TRUE)))
    }
  }
}

# Run the application----
shinyApp(ui = ui, server = server)
# Done