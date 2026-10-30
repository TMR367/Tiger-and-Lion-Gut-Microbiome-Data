
##  LAYOUT

base_dir <- "/Users/subratakumarghosh/Desktop/TigerMicrobiome"

## If you unzip the archive somewhere else, either edit the line above or let
## the script find itself. In RStudio the following line works without editing:
##   base_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
## and from the command line (Rscript R/Tiger_microbiome_all_figures.R) the
## script falls back to the folder it was launched from.

if (!dir.exists(base_dir)) {
  alt <- tryCatch(normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")),
                  error = function(e) NA_character_)
  if (!is.na(alt) && dir.exists(alt)) base_dir <- alt else base_dir <- getwd()
  message("base_dir not found, using: ", base_dir)
}

data_dir <- file.path(base_dir, "data")
out_dir  <- file.path(base_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(dir.exists(data_dir))
message("Reading data from : ", data_dir)
message("Writing output to : ", out_dir)

RARE_DEPTH <- 50398      # even sampling depth used throughout
set.seed(42)             # fixed so that jittered points are reproducible

## Colour scheme, kept identical across all figures ----------------------------
COL_SPECIES  <- c("Tiger" = "#508578", "Lion" = "#E69F00")
COL_SPECIES2 <- c("Royal Bengal tiger" = "#508578", "Lion" = "#E69F00")
COL_HABITAT  <- c("Zoo" = "#0072B2", "Safari" = "#D55E00")
COL_SEX      <- c("Male" = "#009E73", "Female" = "#CC79A7")
COL_AGE      <- c("Adult" = "#0072B2", "Cub" = "#E69F00")
COL_FACILITY <- c("Zoo" = "#0072B2", "Safari Park 1" = "#D55E00",
                  "Safari Park 2" = "#E1A100")
COL_KEGG <- c("Amino acid metabolism" = "#4C72B0",
              "Carbohydrate metabolism" = "#DD8452",
              "Lipid metabolism" = "#55A868",
              "Energy metabolism" = "#C44E52",
              "Vitamin biosynthesis" = "#8172B3",
              "Secondary metabolite biosynthesis" = "#937860",
              "Protein digestion" = "#DA8BC3")


## =============================================================================
## 1. PACKAGES
## =============================================================================

cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "vegan", "ape",
               "pheatmap", "viridis", "RColorBrewer", "VennDiagram", "UpSetR",
               "jsonlite", "gridExtra", "scales", "ragg")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    try(install.packages(p, repos = "https://cloud.r-project.org"), silent = TRUE)
  }
}
library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
suppressMessages({
  library(vegan); library(ape); library(pheatmap); library(viridis)
  library(RColorBrewer); library(VennDiagram); library(UpSetR); library(jsonlite)
})
## VennDiagram writes a log file per call unless this is silenced
if (requireNamespace("futile.logger", quietly = TRUE)) {
  futile.logger::flog.threshold(futile.logger::ERROR)
}

## Common theme. Arial is used throughout; on Linux install the msttcorefonts
## package or substitute "Helvetica" if Arial is unavailable.
FONT <- "Arial"
theme_pub <- theme_bw(base_size = 11, base_family = FONT) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

save_fig <- function(plot, name, width, height) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(name, ".png")), plot, width = width, height = height,
         dpi = 300, device = ragg::agg_png)
  message("  saved: ", name, ".pdf / .png")
}

kfmt <- function(x) paste0(x / 1000, "k")


## =============================================================================
## 2. LOAD DATA
## =============================================================================

message("Loading input tables ...")

meta_full <- read.csv(file.path(data_dir, "sample_metadata_table_full.csv"), check.names = FALSE)
meta <- data.frame(
  sample   = meta_full$Sample_ID,
  animal   = meta_full$Animal_ID,
  Species  = ifelse(meta_full$Species == "Royal Bengal tiger", "Tiger", "Lion"),
  SpeciesL = meta_full$Species,
  Sex      = meta_full$Sex,
  Age      = meta_full$Age,
  Facility = meta_full$Facility,
  stringsAsFactors = FALSE
)
meta$Habitat <- ifelse(meta$Facility == "Zoo", "Zoo", "Safari")
rownames(meta) <- meta$sample

alpha  <- read.csv(file.path(data_dir, "alpha_diversity_values.csv"), row.names = 1, check.names = FALSE)
depth  <- read.csv(file.path(data_dir, "depth_summary.csv"), check.names = FALSE)
dstats <- read.csv(file.path(data_dir, "denoising_stats_merged.csv"), check.names = FALSE)
rare   <- read.csv(file.path(data_dir, "rarefaction_curves.csv"), check.names = FALSE)

## metrics carried on the alpha table
ALPHA_METRICS <- c("Observed ASVs", "Chao1", "ACE", "Shannon", "Simpson", "Faith PD")
ALPHA_LABELS  <- c("Observed ASVs", "Chao1 index", "ACE index",
                   "Shannon index", "Simpson index", "Faith's PD")

read_dm <- function(f) as.dist(as.matrix(read.csv(f, row.names = 1, check.names = FALSE)))
dist_list <- list(
  "Bray-Curtis"        = read_dm(file.path(data_dir, "dm_bray.csv")),
  "Weighted UniFrac"   = read_dm(file.path(data_dir, "dm_wunifrac.csv")),
  "Unweighted UniFrac" = read_dm(file.path(data_dir, "dm_uunifrac.csv"))
)

## align metadata to the distance matrices
ord_ids <- labels(dist_list[[1]])
meta_d  <- meta[ord_ids, ]


## =============================================================================
## 3. HELPER: normality-guided two-group test
##    Shapiro-Wilk first, then Welch t test or Mann-Whitney U, with effect size.
##    This is the rule used for every two-group comparison in the manuscript.
## =============================================================================

shapiro_ok <- function(x) {
  if (length(x) < 3 || length(unique(x)) < 2) return(FALSE)
  isTRUE(tryCatch(shapiro.test(x)$p.value > 0.05, error = function(e) FALSE))
}

two_group_test <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  n1 <- length(a); n2 <- length(b)
  if (shapiro_ok(a) && shapiro_ok(b)) {
    p  <- t.test(a, b, var.equal = FALSE)$p.value
    sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2))
    es <- if (sp > 0) (mean(a) - mean(b)) / sp else NA_real_
    list(test = "Welch t-test", p = p, es_name = "d", es = es, n1 = n1, n2 = n2)
  } else {
    w  <- suppressWarnings(wilcox.test(a, b))
    es <- 1 - 2 * as.numeric(w$statistic) / (n1 * n2)   # rank-biserial
    list(test = "Mann-Whitney U", p = w$p.value, es_name = "r", es = es,
         n1 = n1, n2 = n2)
  }
}

sig_star <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}


## =============================================================================
## 4. TABLE 1  Sample metadata and sequencing summary
## =============================================================================

message("Table 1 ...")

tab1 <- data.frame(
  `Animal ID`             = meta_full$Animal_ID,
  `Sample ID`             = meta_full$Sample_ID,
  Species                 = meta_full$Species,
  Sex                     = meta_full$Sex,
  Age                     = meta_full$Age,
  Facility                = meta_full$Facility,
  `Reads obtained`        = format(meta_full$Reads_obtained, big.mark = ","),
  `Reads after filtering` = format(meta_full$Reads_after_filtering, big.mark = ","),
  `Final ASVs`            = format(meta_full$Final_ASVs, big.mark = ","),
  check.names = FALSE
)
write.csv(tab1, file.path(out_dir, "Table1_Sample_Metadata.csv"), row.names = FALSE)

th <- gridExtra::ttheme_minimal(
  core    = list(fg_params = list(fontfamily = FONT, fontsize = 9),
                 bg_params = list(fill = c("white", "#F2F5F8"), col = "#DDDDDD")),
  colhead = list(fg_params = list(fontfamily = FONT, col = "white", fontface = "bold"),
                 bg_params = list(fill = "#34495E")))
tg  <- gridExtra::tableGrob(tab1, rows = NULL, theme = th)
ttl <- grid::textGrob(
  "Table 1. Sample metadata and sequencing summary for the 23 Panthera faecal samples",
  x = 0.01, just = "left",
  gp = grid::gpar(fontsize = 12, fontface = "bold", fontfamily = FONT))
tab1_grob <- gridExtra::arrangeGrob(ttl, tg, ncol = 1, heights = c(0.07, 0.93))
ggsave(file.path(out_dir, "Table1_Sample_Metadata.pdf"), tab1_grob, width = 14.5, height = 9.2,
       device = cairo_pdf)
ggsave(file.path(out_dir, "Table1_Sample_Metadata.png"), tab1_grob, width = 14.5, height = 9.2,
       dpi = 300, device = ragg::agg_png)


## =============================================================================
## 5. FIGURE 2  Sequencing depth summary
##    (A) raw vs retained reads per sample; (B) DADA2 retention funnel
## =============================================================================

message("Figure 2 ...")

depth$Facility <- factor(depth$Facility, levels = names(COL_FACILITY))

f2a_dat <- depth %>%
  arrange(reads_obtained) %>%
  mutate(sample = factor(sample, levels = sample))

f2a <- ggplot(f2a_dat, aes(y = sample)) +
  geom_col(aes(x = reads_obtained, fill = "Reads obtained (raw)"),
           position = position_nudge(y = 0.2), width = 0.4,
           colour = "black", linewidth = 0.3) +
  geom_col(aes(x = depth, fill = Facility),
           position = position_nudge(y = -0.2), width = 0.4,
           colour = "black", linewidth = 0.3) +
  geom_vline(xintercept = RARE_DEPTH, colour = "#B00000",
             linetype = "dashed", linewidth = 0.6) +
  scale_fill_manual(values = c("Reads obtained (raw)" = "#B8C6D1", COL_FACILITY),
                    name = NULL) +
  scale_x_continuous(labels = kfmt) +
  labs(x = "Number of reads", y = NULL,
       title = "Reads obtained and retained per sample") +
  theme_pub + theme(legend.position = c(0.78, 0.16),
                    panel.grid.major.y = element_blank())

names(dstats) <- make.names(names(dstats))
stage_cols <- c("input", "filtered", "denoised", "merged", "non.chimeric")
stage_lab  <- c("Input\n(raw)", "Filtered", "Denoised", "Merged", "Non-\nchimeric")
stage_mean <- sapply(stage_cols, function(s) mean(dstats[[s]], na.rm = TRUE))
funnel <- data.frame(
  stage = factor(stage_lab, levels = stage_lab),
  mean  = stage_mean,
  pct   = stage_mean / stage_mean[1] * 100)

f2b <- ggplot(funnel, aes(stage, mean, fill = stage)) +
  geom_col(colour = "black", width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.0fk\n(%.0f%%)", mean / 1000, pct)),
            vjust = -0.2, size = 3, family = FONT) +
  scale_fill_manual(values = c("#95A5A6", "#7FB3D5", "#5499C7", "#5DADE2", "#48C9B0")) +
  scale_y_continuous(labels = kfmt, expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Mean reads per sample",
       title = "Read retention through DADA2") +
  theme_pub

fig2 <- (f2a | f2b) + plot_layout(widths = c(1.7, 1)) +
  plot_annotation(title = "Sequencing depth summary",
    theme = theme(plot.title = element_text(face = "bold", size = 15,
                                            hjust = 0.5, family = FONT)))
save_fig(fig2, "Figure2_Sequencing_Depth", 14.5, 6.6)

## Rarefaction curves (supporting the depth assessment) ------------------------
rare$SpeciesL <- ifelse(rare$Species == "Tiger", "Royal Bengal tiger", "Lion")

rar_a <- ggplot(rare, aes(depth, observed, group = sample, colour = SpeciesL)) +
  geom_line(linewidth = 0.5, alpha = 0.75) +
  geom_vline(xintercept = RARE_DEPTH, colour = "#B00000",
             linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = COL_SPECIES2, name = NULL) +
  scale_x_continuous(labels = kfmt) +
  labs(x = "Sequencing depth (reads)", y = "Observed ASVs",
       title = "Per-sample rarefaction (coloured by species)") +
  theme_pub + theme(legend.position = c(0.78, 0.15))

rar_grp <- rare %>%
  group_by(Facility, depth) %>%
  summarise(mean = mean(observed), sd = sd(observed), n = dplyr::n(),
            .groups = "drop") %>%
  filter(n >= 2)
rar_grp$Facility <- factor(rar_grp$Facility, levels = names(COL_FACILITY))

rar_b <- ggplot(rar_grp, aes(depth, mean, colour = Facility, fill = Facility)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = RARE_DEPTH, colour = "#B00000",
             linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = COL_FACILITY) +
  scale_fill_manual(values = COL_FACILITY) +
  scale_x_continuous(labels = kfmt) +
  labs(x = "Sequencing depth (reads)", y = "Observed ASVs",
       title = "Group mean \u00b1 SD (by facility)") +
  theme_pub + theme(legend.position = c(0.78, 0.16))

fig_rar <- (rar_a | rar_b) +
  plot_annotation(title = "Rarefaction curves of observed ASV richness",
    theme = theme(plot.title = element_text(face = "bold", size = 15,
                                            hjust = 0.5, family = FONT)))
save_fig(fig_rar, "Figure2_Rarefaction_Curves", 14, 6.4)


## =============================================================================
## 6. FIGURE 3  Community composition at class (A) and genus (B) level
##    Taxa reaching >= 2% in at least one sample are shown; the rest are Others.
## =============================================================================

message("Figure 3 ...")

comp_panel <- function(file, level_label, ncol_legend, italic_legend) {
  m <- read.csv(file, row.names = 1, check.names = FALSE)
  d <- as.data.frame(t(m))
  d$sample <- rownames(d)
  long <- pivot_longer(d, -sample, names_to = "Taxon", values_to = "abundance")
  ## keep the column order of the matrix (already grouped by species/habitat)
  long$sample <- factor(long$sample, levels = colnames(m))
  ## Others always last, remaining taxa by decreasing mean abundance
  taxa <- rownames(m)
  taxa <- c(setdiff(taxa, "Others"), "Others")
  long$Taxon <- factor(long$Taxon, levels = rev(taxa))

  pal <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B","#E377C2",
           "#17BECF","#BCBD22","#393B79","#E6550D","#31A354","#756BB1","#636363",
           "#3182BD","#FD8D3C","#74C476","#9E9AC8","#969696","#6BAED6","#FDAE6B",
           "#A1D99B","#BCBDDC","#D9D9D9","#AD494A","#8CA252","#B5CF6B","#E7CB94",
           "#843C39","#7B4173","#A55194","#CE6DBD","#DE9ED6","#00868B","#8B0A50")
  cols <- setNames(rep(pal, length.out = length(taxa)), taxa)
  cols["Others"] <- "#C0C0C0"

  ggplot(long, aes(sample, abundance, fill = Taxon)) +
    geom_col(width = 0.82, colour = "black", linewidth = 0.25) +
    scale_fill_manual(values = cols, breaks = taxa,
                      guide = guide_legend(ncol = ncol_legend, reverse = FALSE)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
    labs(x = "Sample", y = "Relative abundance (%)", fill = level_label) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.text = element_text(size = 8,
            face = if (italic_legend) "italic" else "plain"),
          legend.key.size = unit(0.4, "cm"),
          panel.grid = element_blank())
}

f3a <- comp_panel(file.path(data_dir, "relabund_class_2pct.csv"), "Class", 1, FALSE) +
  ggtitle("(A) Class level")
f3b <- comp_panel(file.path(data_dir, "relabund_genus_2pct.csv"), "Genus", 2, TRUE) +
  ggtitle("(B) Genus level")

save_fig(f3a, "Figure3A_Composition_Class", 9.5, 5.6)
save_fig(f3b, "Figure3B_Composition_Genus", 11.0, 6.6)


## =============================================================================
## 7. FIGURES 4 and S1  Alpha diversity
##    Figure 4  : tiger vs lion, zoo vs safari
##    Figure S1 : male vs female, adult vs cub (exploratory, only 2 cubs)
## =============================================================================

message("Figures 4 and S1 ...")

alpha$SpeciesL <- ifelse(alpha$Species == "Tiger", "Royal Bengal tiger", "Lion")
## the alpha table carries habitat as ZooSafari (Zoo / Safari) and sex as Gender;
## alias them so that either naming works
if (!"ZooSafari" %in% names(alpha) && "Habitat" %in% names(alpha)) {
  alpha$ZooSafari <- ifelse(alpha$Habitat == "Captive", "Zoo", "Safari")
}
if (!"Gender" %in% names(alpha) && "Sex" %in% names(alpha)) alpha$Gender <- alpha$Sex

alpha_figure <- function(group_col, levels_in, labels_in, cols, title_txt) {
  panels <- list()
  for (i in seq_along(ALPHA_METRICS)) {
    m  <- ALPHA_METRICS[i]
    df <- data.frame(value = alpha[[m]], grp = alpha[[group_col]])
    df <- df[df$grp %in% levels_in & !is.na(df$value), ]
    df$grp <- factor(df$grp, levels = levels_in)

    a <- df$value[df$grp == levels_in[1]]
    b <- df$value[df$grp == levels_in[2]]
    st <- two_group_test(a, b)

    xlabs <- setNames(
      paste0(labels_in, "\n(n=", c(length(a), length(b)), ")"), levels_in)

    rng  <- range(df$value); span <- diff(rng); if (span == 0) span <- 1
    ybar <- rng[2] + 0.10 * span

    panels[[i]] <- ggplot(df, aes(grp, value, fill = grp)) +
      geom_boxplot(width = 0.55, alpha = 0.85, outlier.shape = NA,
                   colour = "black", linewidth = 0.5) +
      geom_jitter(aes(colour = grp), width = 0.12, height = 0,
                  size = 1.7, shape = 21, colour = "black", stroke = 0.4) +
      annotate("segment", x = 1, xend = 2, y = ybar, yend = ybar, linewidth = 0.4) +
      annotate("text", x = 1.5, y = ybar, label = sig_star(st$p),
               vjust = -0.3, size = 4, fontface = "bold", family = FONT) +
      annotate("label", x = -Inf, y = Inf,
               label = sprintf("%s\np = %.3f (%s)\n%s = %.2f",
                               st$test, st$p, sig_star(st$p), st$es_name, st$es),
               hjust = -0.05, vjust = 1.15, size = 2.5, label.size = 0.2,
               fill = "white", alpha = 0.9, lineheight = 0.9, family = FONT) +
      scale_fill_manual(values = cols, guide = "none") +
      scale_x_discrete(labels = xlabs) +
      scale_y_continuous(expand = expansion(mult = c(0.06, 0.24))) +
      labs(x = NULL, y = ALPHA_LABELS[i]) +
      theme_grey(base_size = 11, base_family = FONT) +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank())
  }
  wrap_plots(panels, nrow = 2) +
    plot_annotation(title = title_txt,
      theme = theme(plot.title = element_text(face = "bold", size = 13,
                                              hjust = 0.5, family = FONT)))
}

fig4a <- alpha_figure("SpeciesL",
                      c("Royal Bengal tiger", "Lion"),
                      c("Royal Bengal\ntiger", "Lion"),
                      COL_SPECIES2,
                      "Alpha diversity: Royal Bengal tiger vs African lion")
save_fig(fig4a, "Figure4_AlphaDiversity_Species", 10.5, 7.2)

fig4b <- alpha_figure("ZooSafari", c("Zoo", "Safari"), c("Zoo", "Safari"),
                      COL_HABITAT, "Alpha diversity: zoo vs safari park")
save_fig(fig4b, "Figure4_AlphaDiversity_Habitat", 10.5, 7.2)

figS1a <- alpha_figure("Gender", c("Male", "Female"), c("Male", "Female"),
                       COL_SEX, "Alpha diversity: male vs female")
save_fig(figS1a, "FigureS1_AlphaDiversity_Sex", 10.5, 7.2)

figS1b <- alpha_figure("Age", c("Adult", "Cub"), c("Adult", "Cub"),
                       COL_AGE,
                       "Alpha diversity: adult vs cub (exploratory, n = 2 cubs)")
save_fig(figS1b, "FigureS1_AlphaDiversity_Age", 10.5, 7.2)


## =============================================================================
## 8. FIGURES 5 and S2  Beta diversity
##    PCoA and NMDS on three distances, with PERMANOVA on each panel.
## =============================================================================

message("Figures 5 and S2 ...")

## PERMANOVA for every distance and factor
perm_tab <- do.call(rbind, lapply(names(dist_list), function(dn) {
  do.call(rbind, lapply(c("Species", "Habitat", "Sex", "Age"), function(f) {
    set.seed(42)
    a <- adonis2(dist_list[[dn]] ~ meta_d[[f]], permutations = 999)
    data.frame(Distance = dn, Factor = f, R2 = a$R2[1], F = a$F[1],
               p = a$`Pr(>F)`[1])
  }))
}))
write.csv(perm_tab, file.path(out_dir, "PERMANOVA_per_factor.csv"), row.names = FALSE)

## combined model: does habitat survive adjustment for species and sex?
perm_comb <- do.call(rbind, lapply(names(dist_list), function(dn) {
  set.seed(42)
  a <- adonis2(dist_list[[dn]] ~ Species + Habitat + Sex, data = meta_d,
               permutations = 999, by = "terms")
  data.frame(Distance = dn, Term = rownames(a), Df = a$Df,
             R2 = round(a$R2, 4), F = round(a$F, 3), p = a$`Pr(>F)`)
}))
write.csv(perm_comb, file.path(out_dir, "PERMANOVA_combined_model.csv"), row.names = FALSE)


## ---------------------------------------------------------------------------
## PERMDISP: homogeneity of multivariate dispersion.
##
## A significant PERMANOVA can come from a shift in group centroid or from one
## group simply being more variable than the other, and the test statistic on
## its own does not separate the two. Every grouping is therefore also run
## through the distance-based dispersion test of Anderson (2006), so that the
## facility result can be read as a difference in composition rather than in
## spread.
## ---------------------------------------------------------------------------

message("PERMDISP ...")

disp_tab <- do.call(rbind, lapply(names(dist_list), function(dn) {
  do.call(rbind, lapply(c("Species", "Habitat", "Sex"), function(f) {
    g  <- factor(meta_d[[f]])
    bd <- vegan::betadisper(dist_list[[dn]], g, type = "median")
    set.seed(42)
    pt <- vegan::permutest(bd, permutations = 999)
    data.frame(Distance   = dn,
               Factor     = f,
               PERMDISP_F = round(pt$tab$F[1], 3),
               PERMDISP_p = pt$tab$`Pr(>F)`[1])
  }))
}))
write.csv(disp_tab, file.path(out_dir, "PERMDISP_per_factor.csv"), row.names = FALSE)
print(disp_tab)


## ---------------------------------------------------------------------------
## Relatedness sensitivity analysis.
##
## Six of the eight tigers sampled at the zoological park belong to a single
## breeding lineage (see Table 1). Because the headline result contrasts the
## zoological park with the safari parks, the facility comparison is repeated
## with that lineage reduced to the one founder male already in the dataset.
## That leaves 17 animals with no known first-degree relatives in either group.
##
## If kinship were producing the facility signal, removing it should remove the
## signal. It does not.
## ---------------------------------------------------------------------------

message("Relatedness sensitivity analysis ...")

RELATED_DROP <- c("CZ-06", "CZ-09", "CZ-14", "CZ-15", "CZ-16", "CZ-20")

kin_tab <- do.call(rbind, lapply(names(dist_list), function(dn) {
  m    <- as.matrix(dist_list[[dn]])
  keep <- setdiff(rownames(m), RELATED_DROP)
  d    <- as.dist(m[keep, keep])
  md   <- meta_d[keep, ]
  set.seed(42)
  a <- adonis2(d ~ Habitat, data = md, permutations = 999)
  data.frame(Distance = dn,
             n        = length(keep),
             n_zoo    = sum(md$Habitat == "Zoo"),
             n_safari = sum(md$Habitat == "Safari"),
             R2       = round(a$R2[1], 4),
             F        = round(a$F[1], 3),
             p        = a$`Pr(>F)`[1])
}))
write.csv(kin_tab, file.path(out_dir, "PERMANOVA_relatedness_sensitivity.csv"),
          row.names = FALSE)
print(kin_tab)


## ordinations, computed once per distance
ords <- lapply(dist_list, function(d) {
  pc  <- ape::pcoa(d)
  pcs <- as.data.frame(pc$vectors[, 1:2]); colnames(pcs) <- c("Axis1", "Axis2")
  pct <- 100 * pc$values$Relative_eig[1:2]
  set.seed(42)
  nm  <- vegan::metaMDS(d, k = 2, trymax = 100, trace = 0)
  nmd <- as.data.frame(vegan::scores(nm, display = "sites"))
  colnames(nmd) <- c("NMDS1", "NMDS2")
  list(pcoa = pcs, pct = pct, nmds = nmd, stress = nm$stress)
})

beta_panel <- function(dn, method, factor_name, cols) {
  o  <- ords[[dn]]
  pr <- perm_tab[perm_tab$Distance == dn & perm_tab$Factor == factor_name, ]
  if (method == "PCoA") {
    df <- data.frame(x = o$pcoa$Axis1, y = o$pcoa$Axis2, Group = meta_d[[factor_name]])
    xl <- sprintf("PCo1 (%.1f%%)", o$pct[1]); yl <- sprintf("PCo2 (%.1f%%)", o$pct[2])
    extra <- NULL
  } else {
    df <- data.frame(x = o$nmds$NMDS1, y = o$nmds$NMDS2, Group = meta_d[[factor_name]])
    xl <- "NMDS1"; yl <- "NMDS2"
    extra <- sprintf("2D stress = %.3f", o$stress)
  }
  df$Group <- factor(df$Group, levels = names(cols))
  lab <- sprintf("PERMANOVA (%s)\nR\u00b2 = %.3f\npseudo-F = %.2f\np = %.3f (%s)",
                 factor_name, pr$R2, pr$F, pr$p, sig_star(pr$p))

  g <- ggplot(df, aes(x, y, colour = Group, fill = Group)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    stat_ellipse(aes(fill = Group), geom = "polygon", type = "norm",
                 level = 0.95, alpha = 0.12, colour = NA, show.legend = FALSE) +
    stat_ellipse(type = "norm", level = 0.95, linewidth = 0.6, show.legend = FALSE) +
    geom_point(size = 2.6, shape = 21, colour = "black", stroke = 0.4) +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    annotate("label", x = -Inf, y = Inf, label = lab, hjust = -0.03, vjust = 1.03,
             size = 2.5, label.size = 0.2, fill = "white", alpha = 0.9,
             lineheight = 0.95, family = FONT) +
    labs(x = xl, y = yl, title = paste0(dn, ": ", method)) +
    theme_pub + theme(plot.title = element_text(size = 10.5, hjust = 0.5),
                      legend.position = "top",
                      legend.title = element_text(face = "bold"))
  if (!is.null(extra)) {
    g <- g + annotate("label", x = Inf, y = -Inf, label = extra,
                      hjust = 1.03, vjust = -0.3, size = 2.5, fontface = "italic",
                      label.size = 0.15, fill = "white", alpha = 0.85, family = FONT)
  }
  g
}

beta_figure <- function(factor_name, cols, title_txt, caveat = NULL) {
  plots <- list()
  for (dn in names(dist_list)) {
    for (m in c("PCoA", "NMDS")) {
      plots[[length(plots) + 1]] <- beta_panel(dn, m, factor_name, cols)
    }
  }
  wrap_plots(plots, ncol = 2, byrow = TRUE, guides = "collect") +
    plot_annotation(title = title_txt, caption = caveat,
      theme = theme(
        plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5,
                                  family = FONT),
        plot.caption = element_text(colour = "red3", face = "italic",
                                    hjust = 0.5, size = 9, family = FONT))) &
    theme(legend.position = "top")
}

fig5a <- beta_figure("Species", COL_SPECIES,
                     "Beta diversity: Royal Bengal tiger vs African lion")
save_fig(fig5a, "Figure5_BetaDiversity_Species", 9.5, 11.5)

fig5b <- beta_figure("Habitat", COL_HABITAT,
                     "Beta diversity: zoo vs safari park")
save_fig(fig5b, "Figure5_BetaDiversity_Habitat", 9.5, 11.5)

figS2a <- beta_figure("Sex", COL_SEX, "Beta diversity: male vs female")
save_fig(figS2a, "FigureS2_BetaDiversity_Sex", 9.5, 11.5)

figS2b <- beta_figure("Age", COL_AGE, "Beta diversity: adult vs cub",
                      caveat = "Cub n = 2 (both captive tigers); exploratory and confounded.")
save_fig(figS2b, "FigureS2_BetaDiversity_Age", 9.5, 11.5)


## =============================================================================
## 9. TABLE 2  Twenty most abundant genera
## =============================================================================

message("Table 2 ...")

t20 <- read.csv(file.path(data_dir, "top20_taxa_table.csv"), check.names = FALSE)
ms  <- function(m, s) sprintf("%.2f \u00b1 %.2f", m, s)
tab2 <- data.frame(
  Rank  = seq_len(nrow(t20)),
  Genus = t20$Genus,
  Overall = ms(t20$Overall_mean, t20$Overall_SD),
  Tiger   = ms(t20$Tiger_mean,   t20$Tiger_SD),
  Lion    = ms(t20$Lion_mean,    t20$Lion_SD),
  Prevalence = sprintf("%.0f", t20$`Prevalence_%`),
  check.names = FALSE)
colnames(tab2) <- c("Rank", "Genus", "Overall (Mean \u00b1 SD %)",
                    "Tiger (Mean \u00b1 SD %)", "Lion (Mean \u00b1 SD %)",
                    "Prevalence (%)")
write.csv(tab2, file.path(out_dir, "Table2_Top20_Genera.csv"), row.names = FALSE)

tg2  <- gridExtra::tableGrob(tab2, rows = NULL, theme = th)
ttl2 <- grid::textGrob(
  "Table 2. Twenty most abundant bacterial genera (relative abundance, mean \u00b1 SD)",
  gp = grid::gpar(fontsize = 12.5, fontface = "bold", fontfamily = FONT))
tab2_grob <- gridExtra::arrangeGrob(ttl2, tg2, ncol = 1, heights = c(0.08, 0.92))
ggsave(file.path(out_dir, "Table2_Top20_Genera.pdf"), tab2_grob, width = 11, height = 8.2,
       device = cairo_pdf)
ggsave(file.path(out_dir, "Table2_Top20_Genera.png"), tab2_grob, width = 11, height = 8.2,
       dpi = 300, device = ragg::agg_png)


## =============================================================================
## 10. FIGURE 6  LEfSe discriminant genera (zoo vs safari)
##     Selection: Kruskal-Wallis p < 0.05 and |LDA| > 2 (uncorrected, per the
##     standard LEfSe procedure). None survived FDR correction; see Results.
## =============================================================================

message("Figure 6 ...")

lef <- read.csv(file.path(data_dir, "lefse_ZooVsSafari.csv"), check.names = FALSE)
lef <- lef[order(lef$LDA), ]
lef$Genus <- factor(lef$Genus, levels = lef$Genus)

fig6 <- ggplot(lef, aes(Genus, LDA, fill = Enriched)) +
  geom_col(colour = "black", width = 0.7, linewidth = 0.6) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = COL_HABITAT, name = "Enriched in") +
  scale_y_continuous(limits = c(-max(abs(lef$LDA)) - 1.4,
                                 max(abs(lef$LDA)) + 1.4)) +
  labs(x = NULL, y = expression(LDA~score~(log[10])),
       title = "Differentially abundant genera between zoo and safari animals") +
  theme_pub +
  theme(axis.text.y = element_text(face = "italic", size = 11),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank(),
        legend.position = c(0.12, 0.15))
save_fig(fig6, "Figure6_LEfSe_ZooVsSafari", 10, 5.2)


## =============================================================================
## 11. FIGURE 7  Core microbiome (four-set Venn plus UpSet)
##     Core = >= 0.1% relative abundance in >= 80% of a group's samples.
## =============================================================================

message("Figure 7 ...")

core80 <- jsonlite::fromJSON(file.path(data_dir, "core80_sets.json"))
oc     <- read.csv(file.path(data_dir, "core_taxa_overall.csv"), check.names = FALSE)

## SILVA writes multi-word genus names with underscores; make them readable
tidy_name <- function(x) gsub("_", " ", x)
core80 <- lapply(core80, tidy_name)
oc$Genus <- tidy_name(oc$Genus)

venn_grob <- venn.diagram(
  x = list(Tiger = core80$Tiger, Lion = core80$Lion,
           Zoo = core80$Zoo, Safari = core80$Safari),
  filename = NULL,
  fill = c("#508578", "#E69F00", "#0072B2", "#D55E00"), alpha = 0.45,
  col = "black", cex = 1.1, fontfamily = FONT,
  cat.fontfamily = FONT, cat.cex = 1.3, cat.fontface = "bold", margin = 0.10)

p_venn <- patchwork::wrap_elements(grid::grobTree(venn_grob)) +
  ggtitle("Core genera shared across groups\n(\u2265 0.1% in \u2265 80% of a group's samples)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11,
                                  family = FONT))

oc$Genus <- factor(oc$Genus, levels = oc$Genus[order(oc[["Mean_relabund_%"]])])
p_bar <- ggplot(oc, aes(Genus, `Mean_relabund_%`, fill = `Core_90%`)) +
  geom_col(colour = "black", width = 0.7) + coord_flip() +
  geom_text(aes(label = paste0(`Prevalence_%`, "%")), hjust = -0.2,
            size = 3, family = FONT) +
  scale_fill_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "#A6BDDB"),
                    labels = c(`TRUE` = "Core at \u2265 90%",
                               `FALSE` = "Core at \u2265 80% only"), name = NULL) +
  labs(x = NULL, y = "Mean relative abundance (%)",
       title = "Overall core genera (\u2265 80% of all samples)") +
  expand_limits(y = max(oc$`Mean_relabund_%`) * 1.15) +
  theme_pub + theme(axis.text.y = element_text(face = "italic"),
                    legend.position = c(0.75, 0.15))

fig7 <- (p_venn | p_bar) + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(title = "Core microbiome of the Panthera gut",
    theme = theme(plot.title = element_text(face = "bold", size = 13.5,
                                            hjust = 0.5, family = FONT)))
save_fig(fig7, "Figure7_CoreMicrobiome_Venn", 13, 5.8)

## UpSet panel (drawn with base graphics, so saved via device calls)
core_sets <- list(Tiger = core80$Tiger, Lion = core80$Lion,
                  Zoo = core80$Zoo, Safari = core80$Safari)
cairo_pdf(file.path(out_dir, "Figure7_CoreMicrobiome_UpSet.pdf"), width = 10, height = 6.5)
print(upset(fromList(core_sets), order.by = "freq", main.bar.color = "#34495E",
            sets.bar.color = "#7F8C8D", text.scale = 1.3,
            mainbar.y.label = "Intersection size", sets.x.label = "Set size"))
dev.off()
ragg::agg_png(file.path(out_dir, "Figure7_CoreMicrobiome_UpSet.png"), width = 10, height = 6.5,
              units = "in", res = 300)
print(upset(fromList(core_sets), order.by = "freq", main.bar.color = "#34495E",
            sets.bar.color = "#7F8C8D", text.scale = 1.3,
            mainbar.y.label = "Intersection size", sets.x.label = "Set size"))
dev.off()
message("  saved: Figure7_CoreMicrobiome_UpSet.pdf / .png")


## =============================================================================
## 12. FIGURE 8  Phylogeny and Faith's phylogenetic diversity
##     (A) circular cladogram of all ASVs, tips coloured by phylum
##     (B) Faith's PD across the four groupings
##     Note: panel A is drawn as a cladogram so that 3,370 tips remain legible;
##     branch lengths are therefore not to scale. Faith's PD in panel B is
##     computed from the full branch-length tree, not from the drawing.
## =============================================================================

message("Figure 8 ...")

tree <- ape::read.tree(file.path(data_dir, "tree.nwk"))

## phylum of each tip, taken from the genus-level composition run upstream.
## If a per-ASV taxonomy file is available, substitute it here.
phylum_file <- file.path(data_dir, "asv_phylum.csv")      # optional: columns asv_id, phylum
if (file.exists(phylum_file)) {
  ph <- read.csv(phylum_file, stringsAsFactors = FALSE)
  tip_phylum <- setNames(ph$phylum, ph$asv_id)[tree$tip.label]
  tip_phylum[is.na(tip_phylum)] <- "Unassigned"
} else {
  tip_phylum <- rep("Unassigned", length(tree$tip.label))
  names(tip_phylum) <- tree$tip.label
  message("  note: asv_phylum.csv not found, tips drawn in a single colour")
}

top_phyla <- c("Firmicutes", "Bacillota", "Actinobacteriota", "Actinomycetota",
               "Proteobacteria", "Pseudomonadota", "Bacteroidota",
               "Fusobacteriota", "Chloroflexi", "Patescibacteria",
               "Cyanobacteria", "Verrucomicrobiota")
ph_pal <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
            "#CCB974", "#64B5CD", "#E377C2", "#7F7F7F")
## colour the nine most ASV-rich phyla, group the remainder as Other
ph_counts <- sort(table(tip_phylum[tip_phylum != "Unassigned"]), decreasing = TRUE)
uph <- names(ph_counts)[seq_len(min(9, length(ph_counts)))]
ph_cols <- setNames(rep(ph_pal, length.out = length(uph)), uph)
tip_cols <- ifelse(tip_phylum %in% names(ph_cols),
                   ph_cols[tip_phylum], "#D9D9D9")
legend_lab  <- if (length(uph)) paste0(uph, " (", ph_counts[uph], ")") else "Unassigned"
legend_fill <- if (length(uph)) unname(ph_cols[uph]) else "#D9D9D9"

cairo_pdf(file.path(out_dir, "Figure8A_Phylogeny.pdf"), width = 9, height = 9)
plot.phylo(tree, type = "fan", show.tip.label = FALSE,
           edge.width = 0.3, edge.color = "grey55", no.margin = TRUE)
tiplabels(pch = 15, col = tip_cols, cex = 0.25)
legend("topleft", legend = legend_lab, fill = legend_fill, bty = "n",
       cex = 0.8, title = "Phylum (no. of ASVs)")
dev.off()
ragg::agg_png(file.path(out_dir, "Figure8A_Phylogeny.png"), width = 9, height = 9, units = "in",
              res = 300)
plot.phylo(tree, type = "fan", show.tip.label = FALSE,
           edge.width = 0.3, edge.color = "grey55", no.margin = TRUE)
tiplabels(pch = 15, col = tip_cols, cex = 0.25)
legend("topleft", legend = legend_lab, fill = legend_fill, bty = "n",
       cex = 0.8, title = "Phylum (no. of ASVs)")
dev.off()
message("  saved: Figure8A_Phylogeny.pdf / .png")

## Panel B: Faith's PD across all four groupings, on one axis
pd_groups <- list(
  list(col = "SpeciesL", lv = c("Royal Bengal tiger", "Lion"),
       lb = c("Bengal\ntiger", "Lion"), cols = COL_SPECIES2),
  list(col = "ZooSafari", lv = c("Zoo", "Safari"),
       lb = c("Zoo", "Safari"), cols = COL_HABITAT),
  list(col = "Gender", lv = c("Male", "Female"),
       lb = c("Male", "Female"), cols = COL_SEX),
  list(col = "Age", lv = c("Adult", "Cub"),
       lb = c("Adult", "Cub"), cols = COL_AGE))

pd_rows <- list(); pd_ann <- list(); xpos <- 0
for (g in pd_groups) {
  for (j in seq_along(g$lv)) {
    v <- alpha$`Faith PD`[alpha[[g$col]] == g$lv[j]]
    v <- v[!is.na(v)]
    xpos <- xpos + 1
    pd_rows[[length(pd_rows) + 1]] <- data.frame(
      x = xpos, value = v,
      lab = paste0(g$lb[j], "\n(n=", length(v), ")"),
      grp = g$lv[j], fill = unname(g$cols[g$lv[j]]),
      stringsAsFactors = FALSE)
  }
  a <- alpha$`Faith PD`[alpha[[g$col]] == g$lv[1]]
  b <- alpha$`Faith PD`[alpha[[g$col]] == g$lv[2]]
  st <- two_group_test(a, b)
  pd_ann[[length(pd_ann) + 1]] <- data.frame(
    x = xpos - 0.5, p = st$p, star = sig_star(st$p),
    lab = sprintf("%s, p = %.2f", ifelse(st$test == "Welch t-test", "t", "MWU"), st$p))
  xpos <- xpos + 0.6      # gap between comparisons
}
pd_df  <- do.call(rbind, pd_rows)
pd_ann <- do.call(rbind, pd_ann)
pd_lab <- pd_df %>% group_by(x) %>% summarise(lab = first(lab), .groups = "drop")

fig8b <- ggplot(pd_df, aes(factor(x), value, fill = I(fill))) +
  geom_boxplot(width = 0.62, alpha = 0.8, outlier.shape = NA, colour = "black") +
  geom_jitter(width = 0.11, height = 0, size = 1.6, shape = 21,
              colour = "black", stroke = 0.4) +
  geom_text(data = pd_ann, aes(x = factor(round(x)), y = max(pd_df$value) + 1.4,
                               label = star),
            inherit.aes = FALSE, fontface = "bold", size = 4, family = FONT) +
  geom_text(data = pd_ann, aes(x = factor(round(x)), y = max(pd_df$value) + 2.4,
                               label = lab),
            inherit.aes = FALSE, size = 2.8, colour = "grey35", family = FONT) +
  scale_x_discrete(labels = setNames(pd_lab$lab, factor(pd_lab$x))) +
  labs(x = NULL, y = "Faith's phylogenetic diversity",
       title = "Faith's PD across groups") +
  theme_pub + theme(plot.title = element_text(hjust = 0.5))
save_fig(fig8b, "Figure8B_FaithPD", 11, 6)


## =============================================================================
## 13. FIGURES 9 and 10 and TABLE 3  Predicted metabolic potential (PICRUSt2)
## =============================================================================

message("Figures 9, 10 and Table 3 ...")

pw    <- read.csv(file.path(data_dir, "picrust_pathway_relabund.csv"), row.names = 1, check.names = FALSE)
pinfo <- read.csv(file.path(data_dir, "picrust_pathway_info.csv"), row.names = 1, check.names = FALSE)
pstat <- read.csv(file.path(data_dir, "picrust_pathway_stats_ZooVsSafari.csv"), check.names = FALSE)
cstat <- read.csv(file.path(data_dir, "picrust_category_stats_ZooVsSafari.csv"), check.names = FALSE)
qc    <- read.csv(file.path(data_dir, "picrust_nsti_qc.csv"), row.names = 1, check.names = FALSE)

pw <- pw[, rownames(meta)[rownames(meta) %in% colnames(pw)], drop = FALSE]

## Figure 9: heatmap of pathways with mean >= 0.05%, plus protein digestion
keep <- union(rownames(pw)[rowMeans(pw) >= 0.05], "ko04974")
keep <- intersect(keep, rownames(pinfo))
mat  <- as.matrix(pw[keep, ])
rownames(mat) <- pinfo[keep, "description"]

ann_col <- data.frame(`Host species` = meta[colnames(mat), "SpeciesL"],
                      Habitat = meta[colnames(mat), "Habitat"],
                      row.names = colnames(mat), check.names = FALSE)
ann_row <- data.frame(`KEGG category` = pinfo[keep, "category"],
                      row.names = rownames(mat), check.names = FALSE)
ann_colors <- list(`Host species` = COL_SPECIES2, Habitat = COL_HABITAT,
                   `KEGG category` = COL_KEGG)

for (dev_ext in c("pdf", "png")) {
  pheatmap(mat, scale = "row",
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks = seq(-2.5, 2.5, length.out = 101),
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "average",
           annotation_col = ann_col, annotation_row = ann_row,
           annotation_colors = ann_colors,
           fontsize_row = 8, fontsize_col = 8.5, border_color = "white",
           main = "Predicted KEGG pathways (PICRUSt2) across samples",
           filename = file.path(out_dir, paste0("Figure9_PICRUSt2_Heatmap.", dev_ext)),
           width = 13.5, height = 15)
}
message("  saved: Figure9_PICRUSt2_Heatmap.pdf / .png")

## Figure 10: bubble plot with fold-change panel
CAT_ORDER <- c("Amino acid metabolism", "Carbohydrate metabolism",
               "Energy metabolism", "Vitamin biosynthesis", "Lipid metabolism",
               "Secondary metabolite biosynthesis", "Protein digestion")
sel <- union(keep, pstat$feature[pstat$p_value < 0.05])
d10 <- pstat %>%
  filter(feature %in% sel) %>%
  mutate(overall = (Zoo_mean + Safari_mean) / 2,
         category = factor(category, levels = CAT_ORDER)) %>%
  arrange(category, overall) %>%
  mutate(description = factor(description, levels = unique(description)))

d10_long <- d10 %>%
  select(description, category, Zoo = Zoo_mean, Safari = Safari_mean) %>%
  pivot_longer(c(Zoo, Safari), names_to = "Group", values_to = "abundance") %>%
  mutate(Group = factor(Group, levels = c("Zoo", "Safari")))

p10a <- ggplot(d10_long, aes(Group, description, size = abundance, fill = category)) +
  geom_point(shape = 21, colour = "black", stroke = 0.4, alpha = 0.85) +
  scale_size_area(max_size = 11, name = "Relative\nabundance (%)",
                  breaks = c(0.05, 0.5, 2, 5)) +
  scale_fill_manual(values = COL_KEGG, guide = "none") +
  labs(x = NULL, y = NULL, title = "Mean relative abundance") +
  theme_pub +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom")

p10b <- ggplot(d10, aes(log2FC, description)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_segment(aes(x = 0, xend = log2FC, yend = description), colour = "grey75") +
  geom_point(aes(fill = category), shape = 21, size = 3, colour = "black",
             stroke = 0.4) +
  geom_text(data = subset(d10, p_value < 0.05), aes(label = "*"),
            hjust = -0.6, vjust = 0.75, size = 6, colour = "#B00000",
            family = FONT) +
  scale_fill_manual(values = COL_KEGG, name = "KEGG category") +
  labs(x = expression(log[2]~fold~change~(Zoo/Safari)), y = NULL,
       title = "Difference between groups") +
  theme_pub +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom")

fig10 <- (p10a | p10b) + plot_layout(widths = c(1, 0.62)) +
  plot_annotation(
    title = "Predicted KEGG pathway abundance (PICRUSt2): zoo versus safari",
    caption = "* nominal p < 0.05; no pathway passed FDR correction",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5, family = FONT),
      plot.caption = element_text(face = "italic", size = 9, hjust = 1, family = FONT)))
save_fig(fig10, "Figure10_PICRUSt2_Bubble", 13.2, 15.5)

## Table 3: KEGG level-2 category comparison
tab3 <- data.frame(
  category   = cstat$feature,
  zoo        = ms(cstat$Zoo_mean, cstat$Zoo_SD),
  safari     = ms(cstat$Safari_mean, cstat$Safari_SD),
  log2fc     = round(cstat$log2FC, 3),
  test       = cstat$test,
  p          = round(cstat$p_value, 4),
  fdr        = round(cstat$FDR, 4),
  effect     = sprintf("%s = %.2f", sub(" .*", "", cstat$effect_name),
                       cstat$effect_size),
  higher_in  = cstat$Enriched,
  check.names = FALSE, stringsAsFactors = FALSE)
colnames(tab3) <- c("KEGG category", "Zoo (mean \u00b1 SD %)",
                    "Safari (mean \u00b1 SD %)", "log2 FC", "Test",
                    "p", "FDR", "Effect size", "Higher in")
write.csv(tab3, file.path(out_dir, "Table3_PICRUSt2_KEGG_Category.csv"), row.names = FALSE)

tg3  <- gridExtra::tableGrob(tab3, rows = NULL, theme = th)
ttl3 <- grid::textGrob(
  "Table 3. Predicted functional capacity (PICRUSt2, KEGG level 2), zoo versus safari",
  gp = grid::gpar(fontsize = 12.5, fontface = "bold", fontfamily = FONT))
tab3_grob <- gridExtra::arrangeGrob(ttl3, tg3, ncol = 1, heights = c(0.14, 0.86))
ggsave(file.path(out_dir, "Table3_PICRUSt2_KEGG_Category.pdf"), tab3_grob, width = 14.5, height = 3.6,
       device = cairo_pdf)
ggsave(file.path(out_dir, "Table3_PICRUSt2_KEGG_Category.png"), tab3_grob, width = 14.5, height = 3.6,
       dpi = 300, device = ragg::agg_png)

cat("\nNSTI weighted mean =", round(mean(qc$weighted_mean_NSTI), 4),
    " range", round(min(qc$weighted_mean_NSTI), 4), "to",
    round(max(qc$weighted_mean_NSTI), 4), "\n")


## =============================================================================
## 14. SESSION INFORMATION
## =============================================================================

cat("\nAll figures and tables written to:", normalizePath(base_dir), "\n\n")
sessionInfo()

###############################################################################
## End of script
###############################################################################
