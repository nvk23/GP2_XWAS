# make a manhattan plot (combined male + female) for the X chromosome
#
# Same logic as miami.R, with the single "combined" result instead of the
# male/female pair. Set REGION_SET below to switch between the two analyses:
#   "non_PAR" -> non-PAR region of chrX (single region)
#   "PAR"     -> PAR1 + PAR2 combined side by side in one plot

if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(reshape2)) install.packages('reshape2')
if (!require(ggrepel)) install.packages('ggrepel')

library(tidyverse)
library(data.table)
library(reshape2)
library(ggrepel)

### --- CONFIG ---
REGION_SET <- "non_PAR"          # "non_PAR" or "PAR"

### Working directory. Everything below is relative to it, so this is the only
### path that changes between machines.
wd        <- "/YOUR/WORKING/DIR"

base_dir  <- file.path(wd, "X_Chrom/locuszoom/filtered")
fig_dir   <- file.path(wd, "X_Chrom/locuszoom/figures")
out_dir   <- file.path(fig_dir, REGION_SET)

### Nearest-gene key. Relative to the repository root, so run this script from
### there; set an absolute path if you keep the key somewhere else.
gene_key_file <- "plots/ref/nearest_genes.csv"

### Minimum -log10(P) kept for plotting. 0 keeps every variant (these are
### single-chromosome files, so this is cheap). Use 3.114074 to thin the cloud.
log10P_floor <- 0

### Blank band between adjacent regions on the x axis, as a fraction of the
### total plotted width. Separates PAR1 from PAR2; ignored for non-PAR, which
### has only one region. 0 butts the regions together.
x_gap_frac <- 0.03

### Window (bp) used to pick one lead variant per locus for labelling.
### 5e5 keeps DHRSX and XG in PAR1 apart (they sit 552 kb from each other).
hit_window <- 5e5

### Saved figure size (inches), the same for both region sets.
fig_width  <- 12
fig_height <- 6

### --- READ INPUT FILES ---
### Each read_gwas() call returns CHR / BP / SNP / P plus a REGION label that
### replaces CHR as the grouping variable along the x axis.
read_gwas <- function(path, region){
  g <- fread(path)
  data.table(
    CHR    = g$CHR,
    BP     = as.numeric(g$POS),
    SNP    = g$rsid,
    P      = as.numeric(g$Pvalue),
    REGION = region
  )
}

if (REGION_SET == "non_PAR") {

  in_dir <- file.path(base_dir, "non_PAR/MAF_0.01")
  gwas_c <- read_gwas(file.path(in_dir, "GP2_FIGS_EUR_AJ_non_PAR_combined_filtered.tsv.gz"), "non-PAR")

} else if (REGION_SET == "PAR") {

  in_dir <- file.path(base_dir, "PAR/MAF_0.01")
  gwas_c <- rbind(
    read_gwas(file.path(in_dir, "PAR1/GP2_EUR_AJ_PAR1_combined_filtered.tsv.gz"), "PAR1"),
    read_gwas(file.path(in_dir, "PAR2/GP2_EUR_AJ_PAR2_combined_filtered.tsv.gz"), "PAR2")
  )

} else {
  stop("REGION_SET must be 'non_PAR' or 'PAR'")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

### --- NEAREST-GENE KEY ---
### rsID -> nearest gene + previously-identified/novel status, used to label the
### hits. Edit nearest_genes.csv to add or correct an entry; a hit that is not
### in the key keeps its rsID and is treated as novel.
gene_key <- fread(gene_key_file)
gene_key <- unique(gene_key[, .(SNP, GENE, STATUS)])

### --- PROCESS THE GWAS FILE ---
process_gwas <- function(gwas) {
  gwas$log10Praw <- -log10(gwas$P)
  gwas$log10P    <- ifelse(gwas$log10Praw > 40, 40, gwas$log10Praw)
  gwas$Plevel    <- NA
  gwas$Plevel[gwas$P < 5E-08] <- "possible"
  gwas$Plevel[gwas$P < 5E-09] <- "likely"
  gwas_filtered <- subset(gwas, log10P > log10P_floor)
  return(gwas_filtered)
}

gwas_c_filt <- process_gwas(gwas_c)

### --- BUILD THE HITS TABLE ---
### Genome-wide significant variants taken from summary stats
### thinned to one lead variant per hit_window so the labels stay readable.
### GENE and STATUS both come from the nearest-gene key above.
get_hits <- function(gwas_filt, window = hit_window){
  sig  <- gwas_filt[gwas_filt$P < 5E-08, ]
  sig  <- sig[order(sig$P), ]
  keep <- sig[0, ]
  for (i in seq_len(nrow(sig))) {
    if (!any(keep$REGION == sig$REGION[i] & abs(keep$BP - sig$BP[i]) < window)) {
      keep <- rbind(keep, sig[i, ])
    }
  }
  i    <- match(keep$SNP, gene_key$SNP)
  gene <- gene_key$GENE[i]
  stat <- gene_key$STATUS[i]
  data.table(SNP    = keep$SNP,
             GENE   = ifelse(is.na(gene), keep$SNP, gene),
             STATUS = ifelse(is.na(stat), "novel", stat))
}

hits_combined <- get_hits(gwas_c_filt)

### Convert to a character vector
combined_highlight_snps <- hits_combined$SNP

### --- CALCULATE REGION OFFSETS (before filtering) ---
### Everything sits on chrX, so REGION takes the place of CHR here: this is what
### lays PAR1 and PAR2 next to each other instead of leaving a 153 Mb gap.
chr_offsets <- gwas_c %>%
  select(REGION, BP) %>%
  group_by(REGION) %>%
  summarize(chr_len = as.numeric(max(BP, na.rm=TRUE) - min(BP, na.rm=TRUE))) %>%
  arrange(REGION)

### Each region after the first is pushed along by its own width plus x_gap
x_gap <- x_gap_frac * sum(chr_offsets$chr_len)

chr_offsets <- chr_offsets %>%
  mutate(tot = lag(cumsum(chr_len + x_gap), default = 0))

### Offsets are relative to the start of each region so the gap before PAR1 and
### between PAR1 and PAR2 is dropped
region_starts <- gwas_c %>%
  select(REGION, BP) %>%
  group_by(REGION) %>%
  summarize(region_start = as.numeric(min(BP, na.rm=TRUE)))

chr_offsets <- left_join(chr_offsets, region_starts, by="REGION")

### --- MODIFIED prep_plot_data to use shared offsets ---
prep_plot_data <- function(gwas_filt, hits, chr_offsets){

  g <- merge(gwas_filt, hits, by="SNP", all.x=TRUE)
  g <- g[order(REGION, BP), ]

  # Use the SHARED region offsets
  g2 <- g %>%
    left_join(chr_offsets, by="REGION") %>%
    mutate(BPcum = as.numeric(BP - region_start + tot)) %>%
    arrange(REGION, BP) %>%
    mutate(
      is_highlight = ifelse(SNP %in% hits$SNP, "yes", "no"),
      is_annotate  = ifelse(log10P > 7.3, "yes", "no")
    )

  return(g2)
}

plot_c <- prep_plot_data(gwas_c_filt, hits_combined, chr_offsets)

### Annotate highlight status
plot_c$is_highlight_c <- ifelse(plot_c$SNP %in% combined_highlight_snps, "yes", "no")

### Label the lead variant of each genome-wide significant locus
plot_c$is_annotate <- ifelse(plot_c$is_highlight_c=="yes" & plot_c$log10P > 7.3, "yes", "no")

plot_c$BPcum <- as.numeric(plot_c$BPcum)

axisdf <- plot_c %>%
  select(REGION, BPcum) %>%
  group_by(REGION) %>%
  summarize(center = mean(range(BPcum, na.rm=TRUE))) %>%
  arrange(REGION)

### --- X AXIS ---
### Both region sets use real chrX positions in Mb. PAR1 and PAR2 are drawn back
### to back, so the breaks are built one region at a time: the tick labels simply
### jump from the end of PAR1 to the start of PAR2.
make_breaks <- function(chr_offsets, n_total = 10){
  b <- lapply(seq_len(nrow(chr_offsets)), function(i){
    r    <- chr_offsets[i, ]
    n    <- max(2, round(n_total * r$chr_len / sum(chr_offsets$chr_len)))
    gpos <- pretty(c(r$region_start, r$region_start + r$chr_len), n = n)
    gpos <- gpos[gpos >= r$region_start & gpos <= r$region_start + r$chr_len]
    data.frame(pos   = gpos - r$region_start + r$tot,
               label = as.character(round(gpos / 1e6, 1)))
  })
  do.call(rbind, b)
}

x_breaks <- make_breaks(chr_offsets)
x_scale  <- scale_x_continuous(breaks = x_breaks$pos, labels = x_breaks$label)
x_title  <- "Chromosome X position (Mb)"

### PAR1 / PAR2 are marked inside the panel rather than on the axis: a divider
### where one region ends and the next begins, plus a label over each region
if (nrow(axisdf) > 1) {
  region_sep <- list(
    geom_vline(xintercept = chr_offsets$tot[-1] - x_gap/2, color="black", linewidth=0.3),
    annotate("text", x = axisdf$center, y = Inf, label = axisdf$REGION,
             vjust = 1.4, size = 3.5, fontface = "bold")
  )
} else {
  region_sep <- NULL
}

### --- Y AXIS ---
### Keep gene labels off the significance line: ggrepel is told to place every
### label box above it.
sig_line    <- 7.30103              ### -log10(5E-08)
label_pad   <- 0.2
label_floor <- sig_line + label_pad

### Gridlines are placed by hand so that none of them runs past the top of the
### data, and the panel is then given a blank band above the highest point for
### the PAR1 / PAR2 labels to sit in. Left to itself the scale adds a break
### above the data and the region labels land on top of it.
y_max <- max(c(plot_c$log10P, sig_line), na.rm=TRUE)
y_brk <- pretty(c(0, y_max), n = 4)
y_brk <- y_brk[y_brk >= 0 & y_brk <= y_max]
y_mnr <- seq(diff(y_brk)[1]/2, y_max, by = diff(y_brk)[1])

y_scale <- scale_y_continuous(
  breaks       = y_brk,
  minor_breaks = y_mnr,
  expand       = expansion(mult = c(0.05, 0.16))
)

### --- CREATE MANHATTAN PLOT ---
manhattan <- ggplot() +

  region_sep +

  ### Combined (male + female) GWAS
  geom_point(data = plot_c,
             aes(x=BPcum, y=log10P, color=as.factor(REGION)),
             alpha=0.8, size=1.3) +

  ### Genome-wide significance (5E-08)
  geom_hline(yintercept = sig_line, linetype="dashed", color="black", linewidth=0.3) +

  ### Highlight significant loci
  geom_point(data=subset(plot_c, is_highlight_c=="yes" & Plevel=="likely"),
             aes(x=BPcum, y=log10P), color="red", size=2) +
  geom_point(data=subset(plot_c, is_highlight_c=="yes" & Plevel=="possible"),
             aes(x=BPcum, y=log10P), color="orange", size=2) +

  ### Variant label annotation
  geom_label_repel(
    data=subset(plot_c, is_annotate=="yes"),
    aes(x=BPcum, y=log10P, label=GENE, fill=factor(STATUS)),
    ylim = c(label_floor, Inf),
    color="black", size=3, max.overlaps=40
  ) +

  ### the 0.5 alpha lives in the fill, so the box stays translucent while the
  ### label text and border stay solid black
  scale_fill_manual(values=c("novel"              = alpha("cyan", 0.5),
                             "previously-identified" = alpha("aquamarine", 0.5))) +
  scale_color_manual(values = rep(c("darkgrey","lightgrey"), 22)) +

  x_scale +
  y_scale +
  ylab("-log10(P)") +
  xlab(x_title) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(manhattan)

### --- SAVE ---
ggsave(file.path(out_dir, paste0("ManhattanPlot_", REGION_SET, ".png")), manhattan,
       width=fig_width, height=fig_height, dpi=600)
