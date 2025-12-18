import argparse
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, MaxNLocator


CHR_MAP = {"X":23,"Y":24,"XY":25,"MT":26,"chrX":23,"chrY":24,"chrXY":25,"chrM":26,"chrMT":26, 'PAR1': 25, 'PAR2': 25}

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

def filter_chr(df, chr_code):
    if df is None or len(df)==0:
        return df
    return df[df["CHR"] == chr_code].copy()

def miami_plot_with_PAR(
    male_nonpar, male_par,
    fem_nonpar, fem_par,
    title="Miami Plot chrX (PAR + non-PAR)",
    genomewide=5e-8, suggestive=1e-6,
    nonpar_color="steelblue",
    par_color="#4f77be",
    point_size=6,
    figsize=(12,6),
    annotate_top_n=0
):
    # prep helper
    def prep(df, color):
        if df is None or len(df)==0:
            return pd.DataFrame(columns=["BP","P","neglog10p","color","SNP"])
        out = df.copy()
        out["neglog10p"] = -np.log10(out["P"].clip(lower=1e-300))
        out["color"] = color
        return out

    # chr codes: non-PAR on X (23), PAR on XY (25) after --split-par b38
    male_nonpar = prep(filter_chr(male_nonpar, 23), nonpar_color)
    male_par    = prep(filter_chr(male_par,    25), par_color)
    fem_nonpar  = prep(filter_chr(fem_nonpar,  23), nonpar_color)
    fem_par     = prep(filter_chr(fem_par,     25), par_color)

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
                               c="crimson", s=point_size+4, linewidths=0, alpha=0.95)
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
                               c="crimson", s=point_size+4, linewidths=0, alpha=0.95)
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
    nonpar_color="steelblue",
    par_color="#4f77be",
    point_size=6,
    figsize=(12, 6),
    annotate_top_n=0,
):
    # prep helper
    def prep(df, color):
        if df is None or len(df) == 0:
            return pd.DataFrame(columns=["BP", "P", "neglog10p", "color", "SNP"])
        out = df.copy()
        out["neglog10p"] = -np.log10(out["P"].clip(lower=1e-300))
        out["color"] = color
        return out

    merged_nonpar = prep(filter_chr(merged_nonpar, 23), nonpar_color)
    merged_par = prep(filter_chr(merged_par, 25), par_color)
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
                           c="crimson", s=point_size + 4, linewidths=0, alpha=0.95)

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
