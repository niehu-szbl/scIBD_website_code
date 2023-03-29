

# Manual of scIBD

Authors: Peilu Lin & Hu Nie

Date: 2023-03-29

*************************

## Outline

*************************

1. Explore gene expression for cell subtypes

2. Explore regulon activity for cell subtypes

3. Compare gene expression between healthy individuals and patients with IBD

4. Compare regulon activity between healthy individuals and patients with IBD

5. Apply gene set enrichment analysis

6. Explore clinical trials, therapy drugs/targets, and GWAS-risk genes of IBD

*************************

## 1. Explore gene expression for cell subtypes

*************************

<p align='center' width='100%'><img src='document/459d19d8be429565545e116ed49b97c3.png' alt='FIGURES' style='zoom:60%;' ></p>

### Example: Explore differential expressed genes in myeloid compartment

*************************

After selecting Myeloid as major cluster, annotations of cell subsets (left) and number of cells (right) 
in myeloid compartment will be shown.

<p align='center' width='100%'><img src='document/4724f078b87cf8d131a29082037074e9.png' alt='FIGURES' style='zoom:100%;' ></p>

<br>

<p align='center' width='100%'><img src='document/89e4d0e8c20a3a0660e47cf706ff1997.png' alt='FIGURES' style='zoom:60%;' ></p>

After choosing **Mast** cell and ranking the list according to **avg_log2FC**, differential expressed genes are listed. You can explore and compare their expression profiles in 'gene expression profile' and 'gene expression comparison'.

> cluster: cell subtypes in myeloid compartment
> 
> p_val: p value calculated by wilcoxon test
> 
> p_val_adj: adjusted p value, based on based on bonferroni correction using all genes in the dataset.
> 
> avg_log2FC: average log 2 fold change

<p align='center' width='100%'><img src='document/2f824d91447de58ba860772947ff4d6d.png' alt='FIGURES' style='zoom:60%;' ></p>

Select one gene of interest

<p align='center' width='100%'><img src='document/30e04aec7482a0a7dd54248898033030.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

Similarly, after choosing LAMP3+ DC and ranking the list according to avg_log2FC, differential expressed genes are listed. You can explore and compare their expression profiles in 'gene expression profile' and 'gene expression comparison'.

<p align='center' width='100%'><img src='document/87540f89e48689e69c62081533ce4ca9.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

Select multiple genes of interest

<p align='center' width='100%'><img src='document/3d0b13a83198f9de4d5f4877067b06e7.png' alt='FIGURES' style='zoom:60%;' ></p>

*************************

## 2. Explore regulon activity for cell subtypes

*************************

After selecting **Myeloid** as major cluster, UMAP embedding (left) and SCENIC embedding (right) with annotations of cell subsets in myeloid compartment will be shown.

<p align='center' width='100%'><img src='document/a21394e199e12a6d7eec1f39463d4aba.png' alt='FIGURES' style='zoom:60%;' ></p>

### Example: Explore regulons in myeloid compartment

*************************

<p align='center' width='100%'><img src='document/53a4b27c223158331b9ed611e860a1f5.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

**Select one regulon of interest**

<p align='center' width='100%'><img src='document/286395180c2563313267cf9a7b0c0d23.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

Activities of regulons in myeloid compartment

<p align='center' width='100%'><img src='document/5dbf897e98ec2d60bfd1b7194e816b25.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

<br>

**Select multiple regulons of interest**

<p align='center' width='100%'><img src='document/4d33e148c8ff6997d4362e57f6e08051.png' alt='FIGURES' style='zoom:60%;' ></p>

The panel 'Network of regulons' interactive

<p align='center' width='100%'><img src='document/128612587c6a60cb630db8235cd2faeb.png' alt='FIGURES' style='zoom:60%;' ></p>

*************************

## 3. Compare gene expression between healthy individuals and patients with IBD

*************************

Compared to the control panel of the first two parts, several additional choices are provided for 'Gene Expression Comparison' part, including 'Tissue', 'Developmental stage', 'Study', 'Minor cluster', 'Location', 'Disease state' and 'Sample'. 

<p align='center' width='100%'><img src='document/792071b09cf422a60b2cc61a9f78d956.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

### Example: Explore the gene expression of HLA-II molecules

*************************

<br>

**Explore the gene expression of HLA-II molecules in all major clusters**

<p align='center' width='100%'><img src='document/cc5003655c6f40bd98dd39d9e24e2b2d.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

<br>

**Explore the gene expression of HLA-II molecules in epithelial cells**

<p align='center' width='100%'><img src='document/98490bb9e794720bfc9bc285c328dd07.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

<br>

**Compare gene expressioin of MHC-II molecules between health and UC in DUOX2+ epithelial cells**

<p align='center' width='100%'><img src='document/734753714b1fa07cf2fedfc7ab1def7d.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

<br>

**Compare gene expressioin of MHC-II molecules between health and CD in enterocytes**

<p align='center' width='100%'><img src='document/ef69c693632551db281250c454c11f5a.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

*************************

## 4. Compare regulon activity between healthy individuals and patients with IBD

*************************

### Example: Explore differentially activated regulons between health and UC or CD in epithelial cells

*************************

Select **Epithelial** as major cluster to explore differentially activated regulons(in this example, **AR**) between health and UC or CD in epithelial cells.

<p align='center' width='100%'><img src='document/b32ce31dceb0c1b28f99272f7e224ddd.png' alt='FIGURES' style='zoom:60%;' ></p>

In the panel 'Compare regulons between healthy and CD (UC)', the value in the second and the third column represents the average regulon activity in inflamed tissue of CD, UC patients or healthy individuals.

<p align='center' width='100%'><img src='document/cfa24edd73f991dd671989d86c8f9778.png' alt='FIGURES' style='zoom:60%;' ></p>

### Exaple: Compare regulon activity of PITX1 between healthy individuals and patients with UC or CD in epithelial cells

*************************

<p align='center' width='100%'><img src='document/358f24267164e547f167b2e4b9f9156f.png' alt='FIGURES' style='zoom:60%;' ></p>

### Example: Compare regulon activity of PITX1 between healthy individuals and patients with UC in DUOX2+ epithelial cells in colorectum

*************************

<p align='center' width='100%'><img src='document/067602bf0bfb0739fd444f21dd66aca7.png' alt='FIGURES' style='zoom:60%;' ></p>

### Example: Compare regulon activity of PITX1 between colon and rectum in DUOX2+ epithelial cells in patients with UC

*************************

<p align='center' width='100%'><img src='document/d1ad279f04ae36995c20bcd813872250.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

*************************

## 5. Apply gene set enrichment analysis

*************************

In the 'Gene Enrichment Analysis' panel, three ways to input a gene set are provided.

You can choose one or more pre-defined risk gene sets from different studies. Here, we choose all pre-defined risk genes of UC **(1)**.After selection, gene list would be generated automatically in the box.

Similarly, you can define your interested gene set **(2)** in the same format (one gene per line) or upload a txt file **(3)** containing your gene set.

<p align='center' width='100%'><img src='document/c1b93c651be7bd3443975dd1b01631c5.png' alt='FIGURES' style='zoom:60%;' ></p>

After applying gene set enrichment analysis on the risk genes of UC, the enrichment of them in each cell types are shown.

> Odds ratio indicates how likely an outcome is to occur in one context relative to another.

<p align='center' width='100%'><img src='document/6341fc8eb5b136d813869b3402e1d376.png' alt='FIGURES' style='zoom:60%;' ></p>

#### Heatmap to show gene expression of each GWAS-risk genes in each cell subtypes

> Rows: genes 
Columns: cell subtypes
Gene expresion are scaled by row.

<p align='center' width='100%'><img src='document/30b3084c63b757200e4fabf34a7546fb.png' alt='FIGURES' style='zoom:60%;' ></p>

*************************

## 6. Explore clinical trials, therapy drugs/targets, and GWAS-risk genes of IBD

*************************

The purpose of the section '**Current Therapy Strategy**' in '**Resources**' is to provide a summary of therapy targets, drugs, and relevant clinical trials for IBD. Here, users can search for clinical trials by disease type (e.g. Crohn's disease, ulcerative colitis), therapy type (e.g. biologics, small molecules), or therapeutic target gene. The results of the search will include a list of clinical trials that match the specified criteria, along with the therapy type and target gene for each trial.

Two parts are included in this section:

1. FDA approved drugs for IBD
2. Therapy targets and drugs for IBD

<p align='center' width='100%'><img src='document/4be0751f20a3e5dc5c9a74cfdf290d42.png' alt='FIGURES' style='zoom:60%;' ></p>

<br>

**Explore clinical trials of IBD**

You can explore FDA approved drugs in the 'Current Therapy Strategy' panel.

Clinical trial information can be reached through the 'View' button.

<p align='center' width='100%'><img src='document/28059d2786b9e50c1c14d22b296ac564.png' alt='FIGURES' style='zoom:60%;' ></p>

The clinical trail information for each drug includes the clinical stage, year of publishment, and links for study record of clinical trail, PubMed page for the reference and more detailed clinical information.

<p align='center' width='100%'><img src='document/160fe157b26c3a58091ce58ae9a42c72.png' alt='FIGURES' style='zoom:60%;' ></p>

You can also explore drugs or targets under clinical trials in this panel.
Click 'View' button for clinical trail information.

<p align='center' width='100%'><img src='document/b9134b231e93cd807f32aa1b4bef7157.png' alt='FIGURES' style='zoom:60%;' ></p>

You can also explore GWAS-risk genes of IBD in '**GWAS-implicated Risk Genes**' panel.

<p align='center' width='100%'><img src='document/96d81e7e5b5901f9bbb4789c2c52f3e0.png' alt='FIGURES' style='zoom:60%;' ></p>

The risk genes of audlt IBD were retrieved from these studies listed in the table of 'Major GWAS study on IBD'

<p align='center' width='100%'><img src='document/172bd3565f0d87349c2fb18b5d0530b1.png' alt='FIGURES' style='zoom:60%;' ></p>

The risk genes of pediatric IBD were retrieved from this paper:

> B. Huang et al., Mucosal profiling of pediatric-onset colitis and IBD reveals common pathogenics and therapeutic pathways. Cell 179, 1160-1176 e1124 (2019).

<p align='center' width='100%'><img src='document/23f559f7f4eee6011d6c119689125581.png' alt='FIGURES' style='zoom:60%;' ></p>
