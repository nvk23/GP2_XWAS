"""chrX plotting helpers for the GP2 XWAS pipeline.

Matplotlib Manhattan and Miami plots for chrX, drawn straight from PLINK2 ``.glm``
output.

Chromosome coding
-----------------
``load_plink_glm`` maps chrX to 23 and PAR (plink's 'XY', or plink2's PAR1/PAR2
from ``--split-par``) to 25, so the plots can colour the two regions differently.

Naming
------
``*_with_PAR``              - non-PAR and PAR drawn together, two colours
``*_PAR_only``              - the pseudoautosomal regions on their own, PAR1 and PAR2
                              in separate panels
``plot_mqq_mb`` / ``miami`` - Manhattan (+ QQ) and mirrored two-panel Miami drawn from a
                              GWASlab frame, with per-locus lead labels

Usage
-----
From a notebook at the repository root::

    from plots.plotting import load_plink_glm, miami_plot_with_PAR

Or as a CLI, from the repository root::

    python plots/plotting.py --male-nonpar male.glm --female-nonpar female.glm --output miami.png
"""

import argparse
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, MaxNLocator
from scipy.stats import chi2


__all__ = [
    # loading
    "load_plink_glm", "filter_chr", "label_par_region", "split_par",
    # plots
    "miami_plot_with_PAR", "manhattan_plot_with_PAR",
    "manhattan_plot_PAR_only", "manhattan_plot_PAR_with_highlights",
    # GWASlab-frame plots
    "plot_mqq_mb", "miami", "get_leads", "nice_ymax",
    "make_xmap", "add_axis_breaks", "draw_manhattan", "draw_qq",
]


CHR_MAP = {"X": 23, "Y": 24, "XY": 25, "MT": 26,
           "chrX": 23, "chrY": 24, "chrXY": 25, "chrM": 26, "chrMT": 26,
           "PAR1": 25, "PAR2": 25}

# Shared palette - every plot in this module defaults to these.
NONPAR_COLOR = "#76b7e5"   # non-PAR (chrX) points
PAR_COLOR = "#4f77be"      # PAR (chrXY) points
SIG_COLOR = "crimson"      # genome-wide significant points

NONPAR_CHR = 23   # chrX outside the pseudoautosomal regions
PAR_CHR = 25      # PAR1 + PAR2, i.e. plink's 'XY'. plink2 --split-par writes
                  # PAR1/PAR2 instead, so CHR_MAP folds both onto 25.

PAR1_END = 2_781_479       # GRCh38 chrX PAR1 ends here
PAR2_START = 155_701_383   # GRCh38 chrX PAR2 starts here

# --- Loading PLINK2 .glm results ---------------------------------------------

def load_plink_glm(path, keep_test="ADD"):
    """Load PLINK2 glm (logistic/logistic.hybrid/linear) and standardize columns."""
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"#CHROM":"CHR","CHROM":"CHR","POS":"BP","ID":"SNP"})
    if "TEST" in df.columns and keep_test is not None:
        df = df[df["TEST"].eq(keep_test)].copy()
    chr_str = df["CHR"].astype(str)
    chr_num = chr_str.map(CHR_MAP).fillna(chr_str)
    df["CHR"] = pd.to_numeric(chr_num, errors="coerce")
    for c in ["BP","P"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["CHR","BP","P"]).copy()
    df["CHR"] = df["CHR"].astype(int)
    if "SNP" not in df.columns:
        df["SNP"] = df.apply(lambda r: f"{int(r['CHR'])}:{int(r['BP'])}", axis=1)
    return df

def filter_chr(df, chr_code, region=None):
    """Keep rows on `chr_code` (an int, or an iterable of ints).

    If the filter empties a frame that had rows, the chromosome coding isn't what the
    caller assumed — e.g. a PAR frame (CHR 25) handed to a non-PAR slot. That would
    otherwise plot as a silently blank panel, so say so.
    """
    if df is None or len(df)==0:
        return df
    codes = {chr_code} if np.isscalar(chr_code) else set(chr_code)
    out = df[df["CHR"].isin(codes)].copy()
    if len(out) == 0:
        where = f' for the {region} panel' if region else ''
        print(f'filter_chr: no rows on CHR {sorted(codes)}{where}; '
              f'frame has CHR {sorted(df["CHR"].unique())} — nothing to plot')
    return out

def label_par_region(df, par1_end=PAR1_END, par2_start=PAR2_START):
    """Return a 'PAR1'/'PAR2' Series for a PAR frame (both regions share chr 25).

    Uses the variant ID when it names the region (plink2 --split-par writes
    chrPAR1:… / chrPAR2:…), else falls back to GRCh38 coordinates.
    """
    region = pd.Series(pd.NA, index=df.index, dtype="object")
    if "SNP" in df.columns:
        ids = df["SNP"].astype(str).str.upper()
        region[ids.str.startswith("CHRPAR1") | ids.str.startswith("PAR1")] = "PAR1"
        region[ids.str.startswith("CHRPAR2") | ids.str.startswith("PAR2")] = "PAR2"
    unlabeled = region.isna()
    if unlabeled.any():
        bp = df.loc[unlabeled, "BP"]
        region[unlabeled] = np.where(bp <= par1_end, "PAR1",
                                     np.where(bp >= par2_start, "PAR2", "PAR1"))
    return region

def split_par(nonpar, par, nonpar_chr=NONPAR_CHR, par_chr=PAR_CHR):
    """Split a single frame that carries both non-PAR and PAR rows into (non-PAR, PAR).

    Lets a whole-chrX result be passed as one frame and still get the two-color treatment.
    Frames that are already separated are returned untouched.
    """
    if par is not None or nonpar is None or len(nonpar) == 0 or "CHR" not in nonpar:
        return nonpar, par
    chroms = set(nonpar["CHR"].unique())
    if {nonpar_chr, par_chr} <= chroms:
        return (nonpar[nonpar["CHR"] == nonpar_chr].copy(),
                nonpar[nonpar["CHR"] == par_chr].copy())
    return nonpar, par


# --- PAR-aware plots (non-PAR + PAR, two colours) ----------------------------

def miami_plot_with_PAR(
    male_nonpar, male_par,
    fem_nonpar, fem_par,
    title="Miami Plot chrX (PAR + non-PAR)",
    genomewide=5e-8, suggestive=1e-6,
    nonpar_color=NONPAR_COLOR,
    par_color=PAR_COLOR,
    point_size=6,
    figsize=(12,6),
    annotate_top_n=0,
    nonpar_chr=NONPAR_CHR,
    par_chr=PAR_CHR
):
    # prep helper
    def prep(df, color):
        if df is None or len(df)==0:
            return pd.DataFrame(columns=["BP","P","neglog10p","color","SNP"])
        out = df.copy()
        out["neglog10p"] = -np.log10(out["P"].clip(lower=1e-300))
        out["color"] = color
        return out

    # chr codes: non-PAR on X (23), PAR on XY (25); a frame holding both is split
    male_nonpar, male_par = split_par(male_nonpar, male_par, nonpar_chr, par_chr)
    fem_nonpar, fem_par = split_par(fem_nonpar, fem_par, nonpar_chr, par_chr)

    male_nonpar = prep(filter_chr(male_nonpar, nonpar_chr, 'male non-PAR'), nonpar_color)
    male_par    = prep(filter_chr(male_par,    par_chr,    'male PAR'),     par_color)
    fem_nonpar  = prep(filter_chr(fem_nonpar,  nonpar_chr, 'female non-PAR'), nonpar_color)
    fem_par     = prep(filter_chr(fem_par,     par_chr,    'female PAR'),     par_color)

    top_df = pd.concat([male_nonpar, male_par], ignore_index=True)
    bot_df = pd.concat([fem_nonpar, fem_par], ignore_index=True)

    y_gw  = -np.log10(genomewide) if genomewide else None
    y_sug = -np.log10(suggestive) if suggestive else None

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=figsize, sharex=True,
        gridspec_kw=dict(hspace=0.1, height_ratios=[1,1])
    )
    fig.suptitle(title, fontsize=16, y=0.98)

    # ---- TOP (males) ----
    if len(top_df):
        ax_top.scatter(top_df["BP"], top_df["neglog10p"],
                       c=top_df["color"], s=point_size, linewidths=0, alpha=0.85)
        if genomewide:
            sig = top_df["P"] <= genomewide
            if sig.any():
                ax_top.scatter(top_df.loc[sig,"BP"], top_df.loc[sig,"neglog10p"],
                               c=SIG_COLOR, s=point_size+4, linewidths=0, alpha=0.95)
    if y_sug: ax_top.axhline(y_sug, ls="--", lw=1, c="gray", alpha=0.8)
    if y_gw:  ax_top.axhline(y_gw,  ls="--", lw=1.2, c="gray", alpha=0.9)
    ax_top.set_ylabel("Males\n" + r"$-\log_{10}(P)$")
    ax_top.set_ylabel(r"$\mathbf{Males}$" + "\n" + r"$-\log_{10}(P)$")
    ax_top.yaxis.set_major_locator(MaxNLocator(nbins=6))

    # ---- BOTTOM (females, mirrored) ----
    if len(bot_df):
        ax_bot.scatter(bot_df["BP"], -bot_df["neglog10p"],
                       c=bot_df["color"], s=point_size, linewidths=0, alpha=0.85)
        if genomewide:
            sig = bot_df["P"] <= genomewide
            if sig.any():
                ax_bot.scatter(bot_df.loc[sig,"BP"], -bot_df.loc[sig,"neglog10p"],
                               c=SIG_COLOR, s=point_size+4, linewidths=0, alpha=0.95)
    if y_sug: ax_bot.axhline(-y_sug, ls="--", lw=1, c="gray", alpha=0.8)
    if y_gw:  ax_bot.axhline(-y_gw,  ls="--", lw=1.2, c="gray", alpha=0.9)
    ax_bot.axhline(0, color="black", lw=1)
    # ax_bot.set_ylabel("Females\n" + r"$-\log_{10}(P)$")
    ax_bot.set_ylabel(r"$\mathbf{Females}$" + "\n" + r"$-\log_{10}(P)$")


    # X axis in Mb
    ax_bot.set_xlabel("Position (Mb)")
    ax_bot.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x/1e6)}"))

    # symmetric y-lims
    ymax = max(
        (top_df["neglog10p"].max() if len(top_df) else 0),
        (bot_df["neglog10p"].max() if len(bot_df) else 0),
        y_gw or 0, y_sug or 0
    ) + 0.5
    ymax = max(1.0, ymax)
    ax_top.set_ylim(0,  ymax)
    ax_bot.set_ylim(-ymax, 0)

    # optional annotations
    def annotate(ax, dframe, sign):
        if annotate_top_n > 0 and len(dframe):
            top = dframe.nsmallest(annotate_top_n, "P")
            for _, r in top.iterrows():
                ax.annotate(
                    r["SNP"], xy=(r["BP"], sign * r["neglog10p"]),
                    xytext=(0, 6 if sign>0 else -10),
                    textcoords="offset points",
                    ha="center", va="bottom" if sign>0 else "top",
                    fontsize=8, color="dimgray"
                )
    annotate(ax_top, top_df, +1)
    annotate(ax_bot, bot_df, -1)

    # clean style
    for ax in (ax_top, ax_bot):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.4)
        ax.set_facecolor("white")

    plt.tight_layout(rect=[0,0,1,0.96])
    return fig

def manhattan_plot_with_PAR(
    merged_nonpar,
    merged_par,
    title="chrX Manhattan Plot (PAR + non-PAR)",
    genomewide=5e-8,
    suggestive=1e-6,
    nonpar_color=NONPAR_COLOR,
    par_color=PAR_COLOR,
    point_size=6,
    figsize=(12, 6),
    annotate_top_n=0,
    nonpar_chr=NONPAR_CHR,
    par_chr=PAR_CHR,
):
    # prep helper
    def prep(df, color):
        if df is None or len(df) == 0:
            return pd.DataFrame(columns=["BP", "P", "neglog10p", "color", "SNP"])
        out = df.copy()
        out["neglog10p"] = -np.log10(out["P"].clip(lower=1e-300))
        out["color"] = color
        return out

    # a single frame holding both non-PAR and PAR rows is split into the two colors
    merged_nonpar, merged_par = split_par(merged_nonpar, merged_par, nonpar_chr, par_chr)
    merged_nonpar = prep(filter_chr(merged_nonpar, nonpar_chr, 'non-PAR'), nonpar_color)
    merged_par = prep(filter_chr(merged_par, par_chr, 'PAR'), par_color)
    df = pd.concat([merged_nonpar, merged_par], ignore_index=True)

    y_gw = -np.log10(genomewide) if genomewide else None
    y_sug = -np.log10(suggestive) if suggestive else None

    fig, ax = plt.subplots(figsize=figsize)
    fig.suptitle(title, fontsize=16, y=0.98)

    if len(df):
        ax.scatter(df["BP"], df["neglog10p"], c=df["color"], s=point_size, linewidths=0, alpha=0.85)
        if genomewide:
            sig = df["P"] <= genomewide
            if sig.any():
                ax.scatter(df.loc[sig, "BP"], df.loc[sig, "neglog10p"],
                           c=SIG_COLOR, s=point_size + 4, linewidths=0, alpha=0.95)

    if y_sug:
        ax.axhline(y_sug, ls="--", lw=1, c="gray", alpha=0.8)
    if y_gw:
        ax.axhline(y_gw, ls="--", lw=1.2, c="gray", alpha=0.9)

    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.set_xlabel("Position (Mb)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x/1e6)}"))

    ymax = max((df["neglog10p"].max() if len(df) else 0), y_gw or 0, y_sug or 0) + 0.5
    ax.set_ylim(0, max(1.0, ymax))

    if annotate_top_n > 0 and len(df):
        top = df.nsmallest(annotate_top_n, "P")
        for _, r in top.iterrows():
            ax.annotate(
                r["SNP"],
                xy=(r["BP"], r["neglog10p"]),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
                color="dimgray",
            )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.4)
    ax.set_facecolor("white")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig

def manhattan_plot_PAR_only(
    merged_par,
    title="chrX Manhattan Plot (PAR only)",
    genomewide=5e-8,
    suggestive=1e-6,
    par_color=PAR_COLOR,
    point_size=6,
    figsize=(12, 6),
    annotate_top_n=0,
    par_chr=PAR_CHR,
    split_regions=True,
):
    """Manhattan plot of the pseudoautosomal regions only.

    PAR1 and PAR2 both live on chr 25 but sit ~153 Mb apart on chrX, so a single
    x-axis squeezes each into a thin spike with a huge empty gap between them —
    which is what a whole-PAR frame (all of chr 25) looks like. By default
    each region gets its own panel, scaled to its own span, sharing one y-axis.

    Parameters
    ----------
    merged_par : pandas.DataFrame
        PAR results from `load_plink_glm`. Needs 'BP', 'P' and (for
        annotation) 'SNP'.
    split_regions : bool
        True (default) draws PAR1 and PAR2 as separate panels when both are
        present. False keeps the old single-axis behaviour.
    """
    # ---- Prep ----
    def prep(df, color):
        if df is None or len(df) == 0:
            return pd.DataFrame(columns=["BP", "P", "neglog10p", "color", "SNP", "PAR"])
        out = df.copy()
        out["neglog10p"] = -np.log10(out["P"].clip(lower=1e-300))
        out["color"] = color
        out["PAR"] = label_par_region(out)
        return out

    # Only PAR (chr25 in your encoding)
    df = prep(filter_chr(merged_par, par_chr, 'PAR'), par_color)

    y_gw = -np.log10(genomewide) if genomewide else None
    y_sug = -np.log10(suggestive) if suggestive else None

    regions = [r for r in ("PAR1", "PAR2") if (df["PAR"] == r).any()] if len(df) else []
    if not split_regions or len(regions) < 2:
        regions = [None]   # one panel holding everything

    # panel widths track each region's span so bp-per-inch stays comparable, but the
    # short region (PAR2 is ~1/8 of PAR1) gets a floor so its axis stays readable
    spans = []
    for r in regions:
        sub = df if r is None else df[df["PAR"] == r]
        spans.append(max(sub["BP"].max() - sub["BP"].min(), 1) if len(sub) else 1)
    widths = [max(s, 0.25 * max(spans)) for s in spans]

    fig, axes = plt.subplots(
        1, len(regions), figsize=figsize, sharey=True, layout="constrained",
        gridspec_kw=dict(width_ratios=widths, wspace=0.06)
    )
    axes = np.atleast_1d(axes)
    fig.suptitle(title, fontsize=16, y=0.98)

    for ax, region in zip(axes, regions):
        sub = df if region is None else df[df["PAR"] == region]

        if len(sub):
            ax.scatter(sub["BP"], sub["neglog10p"], c=sub["color"],
                       s=point_size, linewidths=0, alpha=0.85)

            if genomewide:
                sig = sub["P"] <= genomewide
                if sig.any():
                    ax.scatter(sub.loc[sig, "BP"], sub.loc[sig, "neglog10p"],
                               c=SIG_COLOR, s=point_size + 4, linewidths=0, alpha=0.95)

        # ---- Threshold lines ----
        if y_sug:
            ax.axhline(y_sug, ls="--", lw=1, c="gray", alpha=0.8)
        if y_gw:
            ax.axhline(y_gw, ls="--", lw=1.2, c="gray", alpha=0.9)

        # ---- Axes ----
        # kb reads well for PAR1 (a few Mb from the start of chrX); PAR2 sits at
        # ~155.7 Mb, where kb would be six digits, so that panel switches to Mb
        as_mb = len(sub) and sub["BP"].max() >= 10e6
        unit = "Mb" if as_mb else "kb"
        ax.set_xlabel(f"Position ({unit})" if region is None else f"{region} position ({unit})")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4 if len(regions) > 1 else 6))
        ax.xaxis.set_major_formatter(FuncFormatter(
            (lambda x, _: f"{x/1e6:.2f}") if as_mb else (lambda x, _: f"{int(x/1e3)}")
        ))
        if len(sub):
            pad = max(spans[list(axes).index(ax)] * 0.03, 1)
            ax.set_xlim(sub["BP"].min() - pad, sub["BP"].max() + pad)

        # ---- Optional annotations (top N overall, drawn on its own panel) ----
        if annotate_top_n > 0 and len(sub):
            top = df.nsmallest(annotate_top_n, "P")
            top = top if region is None else top[top["PAR"] == region]
            for _, r in top.iterrows():
                ax.annotate(
                    r["SNP"],
                    xy=(r["BP"], r["neglog10p"]),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color="dimgray",
                )

        # ---- Style ----
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.4)
        ax.set_facecolor("white")

    # extra headroom when annotating, so labels don't run into the title
    headroom = 1.2 if annotate_top_n > 0 else 0.5
    ymax = max((df["neglog10p"].max() if len(df) else 0), y_gw or 0, y_sug or 0) + headroom
    axes[0].set_ylim(0, max(1.0, ymax))
    axes[0].set_ylabel(r"$-\log_{10}(P)$")
    axes[0].yaxis.set_major_locator(MaxNLocator(nbins=6))

    # constrained layout handles the suptitle, so no tight_layout call here
    return fig

def manhattan_plot_PAR_with_highlights(
    merged_par,
    green_variant,              # e.g. "rs123"
    purple_variant,             # e.g. "rs456"
    title="chrX PAR Manhattan Plot",
    genomewide=5e-8,
    suggestive=1e-6,
    par_color=PAR_COLOR,
    point_size=6,
    figsize=(12, 6),
    annotate_top_n=0,
    par_chr=PAR_CHR,
):
    def prep(df, color):
        if df is None or len(df) == 0:
            return pd.DataFrame(columns=["BP", "P", "neglog10p", "color", "SNP"])
        out = df.copy()
        out["neglog10p"] = -np.log10(out["P"].clip(lower=1e-300))
        out["color"] = color
        return out

    df = prep(filter_chr(merged_par, par_chr, 'PAR'), par_color)

    y_gw = -np.log10(genomewide) if genomewide else None
    y_sug = -np.log10(suggestive) if suggestive else None

    fig, ax = plt.subplots(figsize=figsize)
    fig.suptitle(title, fontsize=16, y=0.98)

    if len(df):
        ax.scatter(df["BP"], df["neglog10p"], c=df["color"],
                   s=point_size, linewidths=0, alpha=0.85)

        if genomewide:
            sig = df["P"] <= genomewide
            if sig.any():
                ax.scatter(df.loc[sig, "BP"], df.loc[sig, "neglog10p"],
                           c=SIG_COLOR, s=point_size + 4, linewidths=0, alpha=0.95)

    # ⭐ Highlight variants
    for variant, color in [(green_variant, "limegreen"), (purple_variant, "purple")]:
        hit = df[df["SNP"] == variant]
        if len(hit):
            ax.scatter(hit["BP"], hit["neglog10p"],
                       marker="*", s=200, c=color, edgecolors="black",
                       linewidths=0.6, zorder=10)

    if y_sug:
        ax.axhline(y_sug, ls="--", lw=1, c="gray", alpha=0.8)
    if y_gw:
        ax.axhline(y_gw, ls="--", lw=1.2, c="gray", alpha=0.9)

    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.set_xlabel("Position (kb)")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x/1e3)}"))

    ymax = max((df["neglog10p"].max() if len(df) else 0), y_gw or 0, y_sug or 0) + 0.5
    ax.set_ylim(0, max(1.0, ymax))

    if annotate_top_n > 0 and len(df):
        top = df.nsmallest(annotate_top_n, "P")
        for _, r in top.iterrows():
            ax.annotate(
                r["SNP"],
                xy=(r["BP"], r["neglog10p"]),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
                color="dimgray",
            )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.4)
    ax.set_facecolor("white")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


# --- Manhattan / QQ / Miami from a GWASlab frame ----------------------------
# These take a GWASlab-standardised frame (CHR / POS / P / rsID columns) rather
# than a PLINK .glm one, label one lead variant per locus, and draw the QQ panel
# alongside. make_xmap/add_axis_breaks squeeze the empty stretch between PAR1 and
# PAR2 out of the x-axis.

BLUE, RED = "#597FBD", "#CB132D"   # GWASLab default colors


def get_leads(d, sig, win):
    """Most significant variant per locus, loci separated by > win bp."""
    s = d[d.P < sig].sort_values("P")
    keep = []
    for _, r in s.iterrows():
        if all(abs(r.POS - k.POS) > win for k in keep):
            keep.append(r)
    return pd.DataFrame(keep, columns=d.columns)


def nice_ymax(*ps, sig=5e-8):
    m = max([-np.log10(p.min()) for p in ps] + [-np.log10(sig)])
    return max(10, int(np.ceil(m * 1.4 / 2) * 2))   # headroom for labels


def make_xmap(pos, split_mb=10, gap_frac=0.04):
    """Map bp -> plot x (Mb), closing any gap > split_mb. Returns mapper, ticks, tick labels, break points."""
    pos = np.sort(np.asarray(pos)) / 1e6
    cut = np.where(np.diff(pos) > split_mb)[0]
    starts = np.r_[pos[0], pos[cut + 1]]
    ends = np.r_[pos[cut], pos[-1]]
    gap = gap_frac * (ends - starts).sum()                    # small visual gap between segments
    offsets = np.r_[0, np.cumsum((ends - starts)[:-1] + gap)] - starts

    def xmap(bp):
        mb = np.asarray(bp, float) / 1e6
        i = np.clip(np.searchsorted(starts, mb, side="right") - 1, 0, len(starts) - 1)
        return mb + offsets[i]

    ticks, labels = [], []
    total = (ends - starts).sum()
    for s, e, o in zip(starts, ends, offsets):
        nb = max(2, int(6 * (e - s) / total) + 1)             # more ticks on wider segments
        for t in MaxNLocator(nbins=nb).tick_values(s, e):
            if s <= t <= e:
                ticks.append(t + o)
                labels.append(f"{t:g}")                        # real Mb position
    breaks = [e + o + gap / 2 for e, o in zip(ends[:-1], offsets[:-1])]
    return xmap, ticks, labels, breaks


def add_axis_breaks(ax, breaks, y=0, size=0.012):
    """Draw // across the spine at axes height y (0 = bottom, 1 = top)."""
    for b in breaks:
        xa = ax.transLimits.transform((b, 0))[0]
        ax.plot([xa - 0.004, xa + 0.004], [y, y], transform=ax.transAxes,
                c="white", lw=3, clip_on=False, zorder=4)     # cut the spine
        for dx in (-0.004, 0.004):
            ax.plot([xa + dx - size / 2, xa + dx + size / 2], [y - 0.02, y + 0.02],
                    transform=ax.transAxes, c="black", lw=1, clip_on=False, zorder=5)


def draw_manhattan(ax, d, leads, sig, win, ymax, flip, label_gap=0.035, hl_leads=None, xmap=None):
    xmap = xmap or (lambda bp: np.asarray(bp, float) / 1e6)  # default: plain Mb
    y = -np.log10(d.P.values)
    x = xmap(d.POS.values)

    # colour red: every variant within the window of a highlighted lead
    near = np.zeros(len(d), dtype=bool)
    hl = leads if hl_leads is None else hl_leads
    for pos in hl.POS:
        near |= np.abs(d.POS.values - pos) <= win            # window stays in real bp
    size = np.clip(1 + y * 1.2, 1, 12)

    ax.scatter(x[~near], y[~near], s=size[~near], c=BLUE, lw=0, rasterized=True)
    ax.scatter(x[near], y[near], s=size[near] * 1.5, c=RED, lw=0, rasterized=True, zorder=3)
    ax.axhline(-np.log10(sig), ls=(0, (5, 3)), lw=2, c="grey", zorder=2)

    pad = 0.02 * (x.max() - x.min())
    xlo, xhi = x.min() - pad, x.max() + pad
    base_y = ymax * (0.72 if flip else 0.88)          # keep bottom labels clear of the Mb axis

    L = leads.sort_values("POS")
    if len(L):
        px = xmap(L.POS.values)
        py = -np.log10(L.P.values)

        # spread label positions so neighbours are at least min_gap apart
        min_gap = label_gap * (xhi - xlo)
        lx = px.copy()
        for i in range(1, len(lx)):
            lx[i] = max(lx[i], lx[i - 1] + min_gap)
        over = lx[-1] - (xhi - 0.02 * (xhi - xlo))
        if over > 0:                                  # pushed past the right edge: shift back
            lx -= over
            for i in range(len(lx) - 2, -1, -1):
                lx[i] = min(lx[i], lx[i + 1] - min_gap)

        lab_y = max(base_y, py.max() + ymax * 0.12)   # one label row, clear of the tallest peak
        bend_y = lab_y - ymax * 0.06

        for (_, l), x_l, x_p, y_p in zip(L.iterrows(), lx, px, py):
            lab = l["rsID"] if pd.notna(l.get("rsID")) else f"X:{int(l.POS)}"
            # one continuous connector: stub under the label, angled leg, arrow tip at the peak
            tip_y = y_p + ymax * 0.015
            ax.plot([x_l, x_l, x_p], [lab_y, bend_y, tip_y], c="#D9D9D9", lw=1.2,
                    solid_joinstyle="round", solid_capstyle="butt", zorder=1)
            ax.annotate("", xy=(x_p, tip_y),
                        xytext=(x_p + (x_l - x_p) * 0.02, tip_y + (bend_y - tip_y) * 0.02),
                        arrowprops=dict(arrowstyle="-|>", color="#D9D9D9", lw=1.2,
                                        shrinkA=0, shrinkB=0, mutation_scale=12), zorder=1)
            ax.text(x_l, lab_y + ymax * 0.005 * (-1 if flip else 1), lab,
                    rotation=-40 if flip else 40, rotation_mode="anchor",
                    ha="left", va="top" if flip else "bottom", fontsize=11)

    ax.set_xlim(xlo, xhi)
    ax.set_ylim((ymax, 0) if flip else (0, ymax))
    ax.spines[["top", "right"]].set_visible(False)
    if flip:
        ax.spines["top"].set_visible(True)            # baseline at 0 for the mirrored half


def draw_qq(ax, d, ymax):
    p = np.sort(d.P.values)
    n = len(p)
    obs = -np.log10(p)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    ax.scatter(exp, obs, s=np.clip(1 + obs * 1.5, 1, 14), c=BLUE, lw=0, rasterized=True)
    ax.plot([0, exp.max()], [0, exp.max()], ls="--", c="grey", lw=1.5)
    lam = chi2.isf(np.median(p), 1) / chi2.ppf(0.5, 1)
    ax.text(0.1, 0.97, rf"$\lambda_{{\mathrm{{GC}}}}$ = {lam:.4f}",
            transform=ax.transAxes, va="top", fontsize=11)
    ax.set_ylim(0, ymax)                              # never flipped
    ax.spines[["top", "right"]].set_visible(False)


def plot_mqq_mb(d, sig=5e-8, nominal=None, win_kb=500, out="mqq.png", ymax=None, squeeze_gaps=False):
    """Manhattan + QQ for one sumstats frame.

    nominal: second, looser threshold -- everything past it is labelled, only
             genome-wide hits are coloured red. None draws genome-wide only.
    squeeze_gaps: close gaps > 10 Mb on the x-axis (PAR1 vs PAR2).
    """
    win = win_kb * 1000
    peak = max(-np.log10(d.P.min()), -np.log10(sig))
    ymax = ymax or max(10, int(np.ceil(peak / 0.78)))   # just enough room for the labels
    fig, (ax_m, ax_q) = plt.subplots(1, 2, figsize=(16, 5.2),
                                     gridspec_kw={"wspace": 0.14, "width_ratios": [3.2, 1]})

    gw_leads  = get_leads(d, sig, win)                                  # genome-wide loci -> red
    lab_leads = get_leads(d, nominal, win) if nominal else gw_leads     # labelled

    xmap = ticks = tick_labels = breaks = None
    if squeeze_gaps:
        xmap, ticks, tick_labels, breaks = make_xmap(d.POS.values)

    draw_manhattan(ax_m, d, lab_leads, sig, win, ymax, flip=False,
                   hl_leads=gw_leads, xmap=xmap)
    if nominal:
        ax_m.axhline(-np.log10(nominal), ls=":", lw=1.8, c="grey", zorder=2)   # nominal line
    if squeeze_gaps:
        ax_m.set_xticks(ticks)
        ax_m.set_xticklabels(tick_labels)
        add_axis_breaks(ax_m, breaks)
    ax_m.set_ylabel(r"$-\log_{10}(\mathrm{P})$", fontsize=12)
    ax_m.set_xlabel("Chromosome X position (Mb)", fontsize=12)

    draw_qq(ax_q, d, ymax)
    ax_q.set_ylabel(r"Observed $-\log_{10}(\mathrm{P})$", fontsize=12)
    ax_q.set_xlabel(r"Expected $-\log_{10}(\mathrm{P})$", fontsize=12)
    ax_q.set_xlim(0, None)

    for ax in (ax_m, ax_q):
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))   # whole-number ticks

    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.show()
    return fig


def miami(top, bot, titles, sig=5e-8, nominal=None, win_kb=500, out="miami.png", qq=True,
          squeeze_gaps=False):
    """Two Manhattans mirrored around a shared x-axis, optionally with QQ panels."""
    win = win_kb * 1000
    ymax = nice_ymax(top.P, bot.P, sig=sig)           # same scale on both halves

    xmap = ticks = tick_labels = breaks = None
    if squeeze_gaps:
        # one shared x mapping so both halves line up
        xmap, ticks, tick_labels, breaks = make_xmap(np.r_[top.POS.values, bot.POS.values])

    ncol = 2 if qq else 1
    fig, axes = plt.subplots(2, ncol, figsize=(16 if qq else 12, 9.5), squeeze=False,
                             gridspec_kw={"hspace": 0.06, "wspace": 0.18,
                                          "width_ratios": [3.2, 1] if qq else [1]})
    for row, (d, flip) in enumerate([(top, False), (bot, True)]):
        gw_leads  = get_leads(d, sig, win)                                   # red
        lab_leads = get_leads(d, nominal, win) if nominal else gw_leads      # labelled
        ax = axes[row, 0]
        draw_manhattan(ax, d, lab_leads, sig, win, ymax, flip, hl_leads=gw_leads, xmap=xmap)
        if nominal:
            ax.axhline(-np.log10(nominal), ls=":", lw=1.8, c="grey", zorder=2)
        ax.set_ylabel(titles[row] + "\n" + r"$-\log_{10}(\mathrm{P})$", fontsize=12)
        if qq:
            draw_qq(axes[row, 1], d, ymax)
            axes[row, 1].set_ylabel(titles[row] + "\n" + r"Observed $-\log_{10}(\mathrm{P})$", fontsize=12)

    # shared x range and ticks; labels only under the bottom panel
    xlim = (min(a.get_xlim()[0] for a in axes[:, 0]), max(a.get_xlim()[1] for a in axes[:, 0]))
    for ax in axes[:, 0]:
        ax.set_xlim(xlim)
        if squeeze_gaps:
            ax.set_xticks(ticks)
    axes[0, 0].tick_params(axis="x", labelbottom=False)
    if squeeze_gaps:
        axes[1, 0].set_xticklabels(tick_labels)
        add_axis_breaks(axes[0, 0], breaks, y=0)   # top panel baseline
        add_axis_breaks(axes[1, 0], breaks, y=1)   # mirrored panel baseline
        add_axis_breaks(axes[1, 0], breaks, y=0)   # bottom Mb axis
    else:
        axes[1, 0].set_xticks(axes[0, 0].get_xticks())
        axes[1, 0].set_xlim(xlim)
    axes[1, 0].set_xlabel("Chromosome X position (Mb)", fontsize=12)

    if qq:
        axes[1, 1].set_xlim(axes[0, 1].get_xlim())
        axes[0, 1].tick_params(axis="x", labelbottom=False)
        axes[1, 1].set_xlabel(r"Expected $-\log_{10}(\mathrm{P})$", fontsize=12)

    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.show()
    return fig


# --- Command-line interface --------------------------------------------------

def _optional_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    lowered = str(value).strip().lower()
    if lowered in {"none", "null", "na"}:
        return None
    return float(value)

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a Miami plot for chrX PAR and non-PAR results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--male-nonpar", type=Path, help="PLINK2 glm results for males (non-PAR, chrX).")
    parser.add_argument("--male-par", type=Path, help="PLINK2 glm results for males (PAR, chrXY).")
    parser.add_argument("--female-nonpar", type=Path, help="PLINK2 glm results for females (non-PAR, chrX).")
    parser.add_argument("--female-par", type=Path, help="PLINK2 glm results for females (PAR, chrXY).")
    parser.add_argument("--keep-test", default="ADD", help="Value of the TEST column to keep (set to 'NONE' to keep all).")
    parser.add_argument("--title", default="Miami Plot chrX (PAR + non-PAR)", help="Plot title.")
    parser.add_argument("--genomewide", type=_optional_float, default=5e-8, help="Genome-wide significance P-value threshold.")
    parser.add_argument("--suggestive", type=_optional_float, default=1e-6, help="Suggestive P-value threshold.")
    parser.add_argument("--point-size", type=float, default=6.0, help="Marker size for scatter points.")
    parser.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=(12.0, 6.0),
        help="Figure size in inches.",
    )
    parser.add_argument("--annotate-top-n", type=int, default=0, help="Annotate the top N SNPs per panel.")
    parser.add_argument("--output", type=Path, help="Path to save the figure. If omitted, the plot is shown interactively.")
    parsed = parser.parse_args(argv)

    if not any([parsed.male_nonpar, parsed.male_par, parsed.female_nonpar, parsed.female_par]):
        parser.error("Provide at least one PLINK result file.")
    if isinstance(parsed.genomewide, float) and parsed.genomewide <= 0:
        parser.error("--genomewide must be > 0 or 'none'.")
    if isinstance(parsed.suggestive, float) and parsed.suggestive <= 0:
        parser.error("--suggestive must be > 0 or 'none'.")
    parsed.figsize = (parsed.figsize[0], parsed.figsize[1])
    if parsed.keep_test and parsed.keep_test.strip().upper() == "NONE":
        parsed.keep_test = None
    return parsed

def _load_optional_plink(path: Optional[Path], keep_test: Optional[str]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    return load_plink_glm(path, keep_test=keep_test)

def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)

    male_nonpar = _load_optional_plink(args.male_nonpar, args.keep_test)
    male_par = _load_optional_plink(args.male_par, args.keep_test)
    fem_nonpar = _load_optional_plink(args.female_nonpar, args.keep_test)
    fem_par = _load_optional_plink(args.female_par, args.keep_test)

    fig = miami_plot_with_PAR(
        male_nonpar=male_nonpar,
        male_par=male_par,
        fem_nonpar=fem_nonpar,
        fem_par=fem_par,
        title=args.title,
        genomewide=args.genomewide,
        suggestive=args.suggestive,
        point_size=args.point_size,
        figsize=args.figsize,
        annotate_top_n=args.annotate_top_n,
    )

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
