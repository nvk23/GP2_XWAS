library(locuscomparer)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(LDlinkR)

# free token: https://ldlink.nih.gov/?tab=apiaccess
# Export it in your shell before running, e.g. export LDLINK_TOKEN='your-token'
ldlink_token <- Sys.getenv("LDLINK_TOKEN")

# To run, from the repository root: Rscript plots/locusCompare_setAxes.R

### Working directory. Everything below is relative to it, so this is the only
### path that changes between machines.
wd <- "/YOUR/WORKING/DIR"

base_dir <- file.path(wd, "X_Chrom/locuszoom/filtered")

# ---- locus to plot (hg38; set these) ----
chr      <- "X"
bp       <- "19707642"
lead_snp <- "rs73459702"
# LD anchor. NULL = use lead_snp. Some variants (triallelic sites, and anything
# 1000G phase 3 never called) are in no LD panel at all; point this at a nearby
# variant that is, and the lead stays labelled while r^2 is measured to the proxy.
ld_snp   <- NULL
window   <- 500000   # +/- bp around `bp`

r2_legend       <- TRUE    # draw the r^2 color key on the scatter panel
legend_position <- "topleft" # "bottomright", "topright" or "topleft"

# ---- inputs ----
### Non-PAR
# male_file <- file.path(base_dir, "non_PAR/GP2_FIGS_EUR_AJ_non_PAR_male_filtered.tsv.gz")
# female_file <- file.path(base_dir, "non_PAR/GP2_FIGS_EUR_AJ_non_PAR_female_filtered.tsv.gz")

male_file <- file.path(base_dir, "non_PAR/GP2_EUR_AJ_non_PAR_male_filtered.tsv.gz")
female_file <- file.path(base_dir, "non_PAR/GP2_EUR_AJ_non_PAR_female_filtered.tsv.gz")

### PAR
# male_file <- file.path(base_dir, "PAR/GP2_EUR_AJ_PAR_male_filtered.tsv.gz")
# female_file <- file.path(base_dir, "PAR/GP2_EUR_AJ_PAR_female_filtered.tsv.gz")

out_dir   <- file.path(wd, "X_Chrom/locuszoom/locus_compare/output")

read_metal <- function(in_fn, marker_col = "rsid", pval_col = "pval") {
  if (is.character(in_fn)) {
    d <- read.table(in_fn, header = TRUE, stringsAsFactors = FALSE)
    colnames(d)[colnames(d) == marker_col] <- "rsid"
    colnames(d)[colnames(d) == pval_col]   <- "pval"
  } else if (is.data.frame(in_fn)) {
    d <- in_fn
  } else {
    stop('in_fn must be a string or data.frame')
  }
  d$logp <- -log10(d$pval)
  d[, c("rsid", "pval", "logp")]
}

get_position <- function(x, genome = c("hg19","hg38")) {
  genome <- match.arg(genome)
  stopifnot("rsid" %in% colnames(x))
  stopifnot(exists("pos_lookup"))
  merge(x, pos_lookup, by = "rsid")
}

# LD sources. LDlink's 1000G panel is autosome-only, so any chrX lead SNP comes
# back as "not in 1000G reference panel"; Ensembl REST does cover chrX, with a
# 500 kb max window (so +/-250 kb, narrower than `window`). Each returns NULL
# rather than stopping, so the next source gets a turn.
ld_from_ldlink <- function(snp, population) {
  if (!nzchar(ldlink_token)) {
    warning("No LDlink token. Get one at https://ldlink.nih.gov/?tab=apiaccess and export LDLINK_TOKEN=...")
    return(NULL)
  }
  # NIH serves the idp.nih.gov certificate on ldlink.nih.gov, so the hostname
  # check fails. Skip that check only; ssl_verifypeer stays on, so the chain is
  # still validated against a trusted CA.
  res <- try(httr::with_config(
    httr::config(ssl_verifyhost = 0),
    LDlinkR::LDproxy(
      snp          = snp,
      pop          = population,
      r2d          = "r2",
      token        = ldlink_token,
      genome_build = "grch38",
      win_size     = window
    )
  ), silent = TRUE)

  if (inherits(res, "try-error") || !"RS_Number" %in% colnames(res)) return(NULL)
  res <- res[res$RS_Number != "." & !is.na(res$R2), ]
  if (nrow(res) == 0) return(NULL)
  data.frame(
    SNP_A = snp,
    SNP_B = as.character(res$RS_Number),
    R2    = as.numeric(res$R2),
    stringsAsFactors = FALSE
  )
}

ld_from_ensembl <- function(snp, population) {
  res <- try(httr::GET(
    sprintf("https://rest.ensembl.org/ld/human/%s/1000GENOMES:phase_3:%s?window_size=500",
            snp, population),
    httr::accept_json()
  ), silent = TRUE)

  if (inherits(res, "try-error") || httr::status_code(res) != 200) return(NULL)
  d <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))
  if (!is.data.frame(d) || nrow(d) == 0) return(NULL)
  data.frame(
    SNP_A = snp,
    SNP_B = as.character(d$variation2),
    R2    = as.numeric(d$r2),
    stringsAsFactors = FALSE
  )
}

retrieve_LD <- function(chr, snp, population, anchor = snp) {
  ld <- ld_from_ldlink(anchor, population)

  if (is.null(ld)) {
    message(sprintf("No LDlink LD for %s (its panel is autosome-only); trying Ensembl.", anchor))
    ld <- ld_from_ensembl(anchor, population)
  }
  if (is.null(ld)) {
    warning(sprintf("No LD for %s from LDlink or Ensembl; drawing the plot without r^2 shading.", anchor))
    ld <- data.frame(SNP_A = character(), SNP_B = character(), R2 = numeric(),
                     stringsAsFactors = FALSE)
  }
  # r^2 is to `anchor`, but file it under `snp` so assign_color() still finds it
  # when the anchor is a proxy standing in for a lead that no panel carries.
  if (nrow(ld) > 0) ld$SNP_A <- snp
  ld
}

get_lead_snp <- function(merged, snp = NULL) {
  if (is.null(snp)) {
    snp <- merged[which.min(merged$pval1 + merged$pval2), "rsid"]
  } else if (!snp %in% merged$rsid) {
    stop(sprintf("%s not found in intersection of in_fn1 and in_fn2.", snp))
  }
  as.character(snp)
}

assign_color <- function(rsid, snp, ld) {
  ld <- ld[ld$SNP_A == snp, ]
  ld$color <- as.character(cut(
    ld$R2,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c("blue4", "skyblue", "darkgreen", "orange", "red"),
    include.lowest = TRUE
  ))
  
  color <- data.frame(rsid, stringsAsFactors = FALSE)
  color <- merge(color, ld[, c("SNP_B", "color")],
                 by.x = "rsid", by.y = "SNP_B", all.x = TRUE)
  color[is.na(color$color), "color"] <- "blue4"
  
  if (snp %in% color$rsid) {
    color$color[color$rsid == snp] <- "purple"
  } else {
    color <- rbind(color, data.frame(rsid = snp, color = "purple"))
  }
  
  res <- color$color
  names(res) <- color$rsid
  res
}

add_label <- function(merged, snp) {
  merged$label <- ifelse(merged$rsid %in% snp, merged$rsid, "")
  merged
}

make_scatterplot <- function(merged, title1, title2, color, shape, size) {
  # title2 (FEMALE) on x, title1 (MALE) on y
  ggplot(merged, aes(logp2, logp1)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(
      data = merged[merged$label != "", ],
      aes(logp2, logp1, fill = rsid, size = rsid, shape = rsid)
    ) +
    xlab(bquote(.(title2) ~ -log[10] * "(P)")) +
    ylab(bquote(.(title1) ~ -log[10] * "(P)")) +
    scale_fill_manual(values = color, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf) +
    theme_classic()
}

make_locuszoom <- function(metal, title, chr, color, shape, size) {
  ggplot(metal, aes(x = pos, logp)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(
      data = metal[metal$label != "", ],
      aes(x = pos, logp, fill = rsid, size = rsid, shape = rsid)
    ) +
    scale_fill_manual(values = color, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x / 1e6)) +
    ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf) +
    xlab(paste0("chr", chr, " (Mb)")) +
    ylab(bquote(.(title) ~ -log[10] * "(P)")) +
    theme_classic() +
    theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"))
}

add_r2_legend <- function(p, position = c("bottomright","topright","topleft")) {
  position <- match.arg(position)
  legend_box <- if (position == "bottomright") {
    data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
  } else if (position == "topright") {
    data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
  } else {
    data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
  }

  cowplot::ggdraw(p) +
    geom_rect(
      data = legend_box,
      aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
      color = "black",
      fill  = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))
    ) +
    cowplot::draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
    cowplot::draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
    cowplot::draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
    cowplot::draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
    cowplot::draw_label(parse(text = "r^2"),
                        x = legend_box$x[1] + 0.025,          # centred on the column
                        y = legend_box$y[1] + 0.05 + 0.015,   # clear of the top swatch
                        hjust = 0.5, vjust = 0, size = 10)
}

make_combined_plot <- function(merged, title1, title2, ld, chr, snp = NULL) {
  snp   <- get_lead_snp(merged, snp)
  color <- assign_color(merged$rsid, snp, ld)
  
  shape <- ifelse(merged$rsid == snp, 23, 21)
  names(shape) <- merged$rsid
  
  size <- ifelse(merged$rsid == snp, 3, 2)
  names(size) <- merged$rsid
  
  merged <- add_label(merged, snp)
  
  max_logp <- max(c(merged$logp1, merged$logp2), na.rm = TRUE)
  if (!is.finite(max_logp) || max_logp <= 0) max_logp <- 1
  max_logp <- ceiling(max_logp)
  
  p1 <- make_scatterplot(merged, title1, title2, color, shape, size) +
    scale_x_continuous(limits = c(0, max_logp)) +
    scale_y_continuous(limits = c(0, max_logp))

  # must come after the scales: ggdraw() output can no longer take scale_*()
  if (isTRUE(r2_legend)) p1 <- add_r2_legend(p1, legend_position)
  
  metal1 <- merged[, c("rsid", "logp1", "chr", "pos", "label")]
  colnames(metal1)[colnames(metal1) == "logp1"] <- "logp"
  p2 <- make_locuszoom(metal1, title1, chr, color, shape, size) +
    scale_y_continuous(limits = c(0, max_logp))
  
  metal2 <- merged[, c("rsid", "logp2", "chr", "pos", "label")]
  colnames(metal2)[colnames(metal2) == "logp2"] <- "logp"
  p3 <- make_locuszoom(metal2, title2, chr, color, shape, size) +
    scale_y_continuous(limits = c(0, max_logp))
  
  p2 <- p2 + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank())
  
  right <- cowplot::plot_grid(p2, p3, align = "v", nrow = 2,
                              rel_heights = c(0.8, 1))
  
  cowplot::plot_grid(
    p1, right,
    align      = "h",
    rel_widths = c(1.6, 2)
  )
}

locuscompare <- function(in_fn1, in_fn2,
                         title1 = "eQTL", title2 = "GWAS",
                         snp = NULL, ld_snp = NULL, population = "EUR",
                         genome = c("hg19","hg38")) {
  
  d1 <- read_metal(in_fn1)
  d2 <- read_metal(in_fn2)
  
  merged <- merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
  genome <- match.arg(genome)
  merged <- get_position(merged, genome)
  
  chr <- unique(merged$chr)
  if (length(chr) != 1) stop("There must be one and only one chromosome.")
  
  snp <- get_lead_snp(merged, snp)
  ld  <- retrieve_LD(chr, snp, population,
                     anchor = if (is.null(ld_snp)) snp else ld_snp)
  make_combined_plot(merged, title1, title2, ld, chr, snp)
}


start <- max(0, as.numeric(bp) - window)
end   <- as.numeric(bp) + window

region_cmd <- function(fn) {
  sprintf("gzip -dc %s | awk -F'\\t' 'NR==1 || ($1==\"%s\" && $2>=%.0f && $2<=%.0f)'",
          shQuote(fn), chr, start, end)
}

gwas_fn <- fread(cmd = region_cmd(male_file))
setnames(gwas_fn, 1, "chrom")
gwas_fn <- gwas_fn[rsid != "." & !is.na(neg_log_pvalue)]
gwas_fn[, pval := 10^(-neg_log_pvalue)]
gwas_fn <- unique(gwas_fn, by = "rsid")
pos_male <- gwas_fn[, .(rsid, chr = as.character(chrom), pos = as.numeric(pos))]
gwas_fn <- as.data.frame(gwas_fn[, .(rsid, pval)])

eqtl_fn <- fread(cmd = region_cmd(female_file))
setnames(eqtl_fn, 1, "chrom")
eqtl_fn <- eqtl_fn[rsid != "." & !is.na(neg_log_pvalue)]
eqtl_fn[, pval := 10^(-neg_log_pvalue)]
eqtl_fn <- unique(eqtl_fn, by = "rsid")
pos_female <- eqtl_fn[, .(rsid, chr = as.character(chrom), pos = as.numeric(pos))]
eqtl_fn <- as.data.frame(eqtl_fn[, .(rsid, pval)])

pos_lookup <- unique(rbind(pos_male, pos_female), by = "rsid")
pos_lookup <- as.data.frame(pos_lookup)

p <- locuscompare(
  in_fn1 = gwas_fn,
  in_fn2 = eqtl_fn,
  title1 = "MALE",
  title2 = "FEMALE",
  snp    = lead_snp,
  ld_snp = ld_snp,
  population = "EUR",
  genome     = "hg38"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
base_name <- paste0("locuscompare_chr", chr, "_", bp, "_", lead_snp)

ggsave(
  filename = file.path(out_dir, paste0(base_name, ".png")),
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)
