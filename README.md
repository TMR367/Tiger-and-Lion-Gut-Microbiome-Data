# Tiger-and-Lion-Gut-Microbiome-Data
Gut bacterial communities and predicted metabolic potential in captive Royal Bengal tigers (Panthera tigris tigris) and African lions (Panthera leo) in Bangladesh
Gut Microbiome of Captive Royal Bengal Tigers and African Lions in Bangladesh

This repository contains the R scripts and supporting analysis workflows used for the study of fecal bacterial communities in captive Royal Bengal tigers (Panthera tigris tigris) and African lions (Panthera leo) maintained in zoological and safari facilities in Bangladesh.

Study overview

The study investigated the gut microbiome of 23 captive large felids, comprising 16 Royal Bengal tigers and 7 African lions, sampled from one zoological park and two safari parks. Bacterial communities were characterized using 16S rRNA gene V3–V4 amplicon sequencing.

The analyses were designed to assess:

Taxonomic composition of the gut microbiota Alpha diversity and microbial richness Beta diversity and community structure Host species-, sex-, age-, and facility-associated variation Differentially abundant bacterial genera Core microbiome composition Phylogenetic diversity Predicted microbial metabolic potential using PICRUSt2 Main analytical workflow

Sequence processing and ASV generation were performed using QIIME 2 and DADA2. Taxonomic assignments were based on the SILVA 138.2 database. Downstream statistical analyses and visualization were performed primarily in R.

The repository includes scripts for:

Importing and organizing microbiome abundance and metadata files Taxonomic abundance summaries and visualization Alpha-diversity analyses Bray–Curtis, weighted UniFrac, and unweighted UniFrac beta-diversity analyses PCoA and NMDS ordination PERMANOVA and multivariate dispersion analyses Differential abundance analysis between zoo and safari populations Core microbiome identification and visualization Phylogenetic diversity analyses PICRUSt2-derived KEGG functional comparisons Generation of publication-quality tables and figures Statistical analysis

Depending on data distribution, comparisons were performed using Welch's t-test or the Mann–Whitney U test, with corresponding effect sizes. Multiple comparisons were controlled using the Benjamini–Hochberg false discovery rate (FDR) where appropriate.

Community-level differences were evaluated using PERMANOVA with 999 permutations, and multivariate dispersion was assessed before interpretation of PERMANOVA results.

Differentially abundant genera were evaluated using Kruskal–Wallis/LDA-based screening, with complementary FDR-adjusted statistical analysis.

Core microbiome

Core bacterial genera were defined as taxa present at a relative abundance of ≥0.1% in at least 80% of samples within the relevant group. A more stringent 90% prevalence threshold was also evaluated.

Functional prediction

Predicted microbial functional potential was generated using PICRUSt2 and summarized using KEGG orthologues and broader metabolic categories. These analyses represent predicted functional potential rather than directly measured microbial activity.

Software

Major software and R packages used include:

QIIME 2 DADA2 PICRUSt2 R vegan ape picante ggplot2 pheatmap viridis VennDiagram UpSetR

Exact software versions and analytical parameters are provided in the associated manuscript and individual scripts where applicable.

Reproducibility

Scripts are organized according to the analytical workflow used in the manuscript. Users should update file paths and input filenames according to their local directory structure before running the analyses.

Where randomization or permutation procedures are used, fixed seeds should be retained to ensure reproducibility.

Data availability

Raw 16S rRNA gene sequencing data are deposited in the NCBI Sequence Read Archive (SRA) under the accession number provided in the associated manuscript. Processed data, metadata, and analysis scripts available in this repository are intended to facilitate transparency and reproducibility of the study.

Citation

If you use these scripts or data, please cite the associated publication:

Bijoy Chowdhury, Mira Rani Baidhya, Subrata Kumar Ghosh, Mahabub Alam, Md Sarwar Uddin, Md Shahadat Hossain Suvo, Mohammad Mustafizur Rahman, Hatem Sazzat Md Zulkar Nine, Lipi Akter, Yoko Makabe-Kobayashi, Yusuke Matsumoto, Koji Hamasaki, Tofazzal Md Rakib. 2026. A conserved core and a variable fringe: gut bacterial communities and predicted metabolic potential in captive Royal Bengal tigers (Panthera tigris tigris) and African lions (Panthera leo) in Bangladesh. [Journal details to be added after publication].

Contact

For questions regarding the analysis or reuse of the scripts, please contact the corresponding author through the contact information provided in the associated publication.
