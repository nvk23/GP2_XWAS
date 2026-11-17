"""
METAL meta-analysis of PLINK2 logistic GWAS results (*.glm.logistic.hybrid).
Based on run_metal.nf from https://github.com/GP2code/progwas-pipeline.

Per input file:
  * keep TEST == "ADD" (drops covariate rows if --hide-covar wasn't used)
  * drop rows with ERRCODE != "." and rows with missing / non-finite OR, SE, P
  * BETA = ln(OR), SE = LOG(OR)_SE
  * A2 = the other of REF/ALT (A1==REF -> ALT, A1==ALT -> REF), else OMITTED
  * write a slim, gzipped METAL-ready table to <outdir>/prepared/

Usage in Jupyter, from the repository root (either):
    %run src/run_metal_plink2.py      # then call run_metal(...) in the next cell
or
    from src.run_metal_plink2 import run_metal
    run_metal(metal_input="/data/gwas/*/*.glm.logistic.hybrid", metal_outdir="metal_results")

METAL runs inside the GP2 pipeline image by default
(ghcr.io/gp2code/progwas-pipeline:1.0.0, same as conf/profiles/standard.config),
whose METAL build supports ANALYZE RANDOM.
  * docker                 -> docker run --platform linux/amd64 -v <outdir>:<outdir> ...

Set container=None to use a local METAL binary instead (metal_bin= / $METAL_BIN).
"""

# %% Imports
import glob
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", force=True)
log = logging.getLogger("metal")

PLINK_COLS = ["#CHROM", "POS", "ID", "REF", "ALT", "A1", "OMITTED",
              "A1_FREQ", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "P", "ERRCODE"]

DEFAULT_CONTAINER = "ghcr.io/gp2code/progwas-pipeline:1.0.0"

# Columns a fixed-effects (+ HETEROGENEITY) METAL run produces; anything beyond these
# in an ANALYZE RANDOM run means the random-effects model was actually computed.
FIXED_TBL_COLS = {"MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq",
                  "Effect", "StdErr", "P-value", "Direction", "TotalSampleSize",
                  "HetISq", "HetChiSq", "HetDf", "HetPVal"}

PLINK_SUFFIX_RE =r"\.glm\.logistic(\.hybrid)?(\.gz|\.zst)?$"


# %% Helpers
def resolve_inputs(metal_input):
    """Glob pattern(s): a string (comma-separated allowed) or a list."""
    patterns = metal_input if isinstance(metal_input, (list, tuple)) else str(metal_input).split(",")
    files = []
    for pat in (p.strip() for p in patterns if str(p).strip()):
        hits = sorted(glob.glob(os.path.expanduser(pat)))
        if not hits:
            raise FileNotFoundError(f"No files matched: {pat}")
        files.extend(hits)
    # de-duplicate, keep order
    return list(dict.fromkeys(os.path.abspath(f) for f in files))


def match_pheno(filename, pheno_names):
    """PLINK2 names outputs <out>.<pheno>.glm.logistic.hybrid; pick the longest matching phenotype."""
    base = re.sub(PLINK_SUFFIX_RE, "", os.path.basename(filename))
    hits = [p for p in pheno_names if base.endswith(f".{p}") or base == p]
    return max(hits, key=len) if hits else None


def group_by_pheno(files, metal_pheno_name):
    names = [p.strip() for p in (metal_pheno_name or "").split(",") if p.strip()] \
        if isinstance(metal_pheno_name, (str, type(None))) else list(metal_pheno_name)
    if not names:
        return {"": files}
    groups = {}
    for f in files:
        pheno = match_pheno(f, names)
        if pheno is None:
            raise ValueError(f"File does not match any metal_pheno_name ({', '.join(names)}): {f}")
        groups.setdefault(pheno, []).append(f)
    return groups


def prepare_file(path, out_path, test="ADD", marker_id="ID"):
    """Convert one PLINK2 logistic file to a METAL-ready TSV.gz. Returns (stats, marker->chr/pos table)."""
    df = pd.read_csv(
        path, sep="\t", usecols=lambda c: c in PLINK_COLS,
        dtype=str, keep_default_na=False, compression="infer",
    )
    missing = [c for c in PLINK_COLS if c not in df.columns and c != "OMITTED"]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    n_raw = len(df)

    df = df[df["TEST"] == test]
    n_test = len(df)
    df = df[df["ERRCODE"] == "."]
    n_err = n_test - len(df)

    or_ = pd.to_numeric(df["OR"], errors="coerce")
    se = pd.to_numeric(df["LOG(OR)_SE"], errors="coerce")
    p = pd.to_numeric(df["P"], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.log(or_)
    ok = np.isfinite(beta) & np.isfinite(se) & (se > 0) & p.between(0, 1)
    n_bad = int((~ok).sum())
    df, beta, se, p = df[ok], beta[ok], se[ok], p[ok]

    # A2 from REF/ALT (as in run_metal.nf), fall back to OMITTED for multi-allelic / odd cases
    a1_ref = df["A1"] == df["REF"]
    a1_alt = df["A1"] == df["ALT"]
    a2 = np.where(a1_ref, df["ALT"], np.where(a1_alt, df["REF"], "."))
    if "OMITTED" in df.columns:
        a2 = np.where(a2 == ".", df["OMITTED"], a2)
    a2 = pd.Series(a2, index=df.index)
    has_a2 = (a2 != ".") & (a2 != "")
    n_no_a2 = int((~has_a2).sum())

    if marker_id == "ID":
        marker = df["ID"]
    elif marker_id == "chrpos":
        # chr:pos:alleles in sorted order -> matches across cohorts regardless of REF/ALT orientation
        lo = np.minimum(df["A1"].str.upper(), a2.str.upper())
        hi = np.maximum(df["A1"].str.upper(), a2.str.upper())
        chrom = df["#CHROM"].str.replace(r"^chr", "", regex=True)
        marker = chrom + ":" + df["POS"] + ":" + lo + ":" + hi
    else:
        raise ValueError("marker_id must be 'ID' or 'chrpos'")

    out = pd.DataFrame({
        "MARKER": marker, "CHR": df["#CHROM"], "POS": df["POS"],
        "REF": df["REF"], "ALT": df["ALT"],
        "A1": df["A1"], "A2": a2, "A1_FREQ": df["A1_FREQ"],
        "BETA": beta, "SE": se, "P": p, "OBS_CT": df["OBS_CT"],
    })[has_a2.values]

    n_dup = int(out["MARKER"].duplicated().sum())
    if n_dup:
        log.warning("  %d duplicated MARKER values in %s (METAL keeps the first). "
                    "Consider marker_id='chrpos' or unique IDs.", n_dup, os.path.basename(path))
    if (out["MARKER"] == ".").any():
        log.warning("  %d rows have ID '.'; use marker_id='chrpos'.", int((out["MARKER"] == ".").sum()))

    out.to_csv(out_path, sep="\t", index=False, compression="gzip", float_format="%.6g")

    stats = dict(file=os.path.basename(path), rows_in=n_raw, rows_test=n_test,
                 dropped_errcode=n_err, dropped_na=n_bad, dropped_no_a2=n_no_a2,
                 A1_EQ_REF=int(a1_ref[ok].sum()), A1_EQ_ALT=int(a1_alt[ok].sum()),
                 rows_out=len(out))
    return stats, out[["MARKER", "CHR", "POS", "REF", "ALT"]]


def build_script(prefix, prepared_files, analyze="RANDOM"):
    lines = [
        "SCHEME STDERR",
        "GENOMICCONTROL OFF",
        "AVERAGEFREQ ON",
        "MINMAXFREQ ON",
        "CUSTOMVARIABLE TotalSampleSize",
        "LABEL TotalSampleSize AS OBS_CT",
        "MARKER MARKER",
        "FREQ A1_FREQ",
        "ALLELE A1 A2",
        "EFFECT BETA",
        "STDERR SE",
        "PVAL P",
        "WEIGHT OBS_CT",
        "",
    ]
    lines += [f"PROCESS {f}" for f in prepared_files]
    lines += ["", f"OUTFILE {prefix} .TBL", f"ANALYZE {analyze}".strip(), "QUIT", ""]
    return "\n".join(lines)


def find_metal(metal_bin=None):
    cand = metal_bin or os.environ.get("METAL_BIN") or shutil.which("metal")
    if not cand or not (os.path.isfile(cand) and os.access(cand, os.X_OK)):
        raise FileNotFoundError(
            "METAL binary not found. Pass metal_bin='/path/to/metal' or set os.environ['METAL_BIN'] "
            "(on Biowulf: `module show metal` to find the path).")
    return os.path.abspath(cand)


def _docker_usable(docker):
    try:
        return subprocess.run([docker, "info"], capture_output=True, timeout=20).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def find_engine(container_engine="auto", engine_bin=None):
    """Return (engine_name, engine_path). engine_name is 'docker' or 'singularity'-style."""
    if engine_bin:
        name = os.path.basename(engine_bin)
        kind = "docker" if "docker" in name or "podman" in name else "singularity"
        return kind, engine_bin
    order = {"auto": ["docker", "apptainer", "singularity"],
             "docker": ["docker"], "podman": ["podman"],
             "apptainer": ["apptainer"], "singularity": ["singularity"]}[container_engine]
    for exe in order:
        path = shutil.which(exe)
        if not path:
            continue
        if exe in ("docker", "podman"):
            if _docker_usable(path):
                return "docker", path
            if container_engine != "auto":
                raise RuntimeError(f"{exe} found but not usable (daemon not running / no permission)")
            continue
        return "singularity", path
    raise FileNotFoundError(
        f"No container engine found for container_engine='{container_engine}'. Install/enable Docker, "
        "or pass engine_bin='/path/to/apptainer' (Biowulf: `module show apptainer`), "
        "or set container=None and metal_bin= to use a local METAL.")


def pull_container(sif_path, container=DEFAULT_CONTAINER, container_engine="auto", engine_bin=None):
    """Build a local .sif from the Docker image (Apptainer/Singularity). No-op if it already exists."""
    sif_path = Path(sif_path).expanduser().resolve()
    if sif_path.exists():
        log.info("Container already present: %s", sif_path)
        return str(sif_path)
    if engine_bin:
        exe = engine_bin
    elif container_engine == "auto":
        exe = shutil.which("apptainer") or shutil.which("singularity")
    else:
        exe = shutil.which(container_engine)
    if not exe or "docker" in os.path.basename(exe) or "podman" in os.path.basename(exe):
        raise RuntimeError("pull_container builds a .sif and needs apptainer/singularity (pass engine_bin= "
                           "if it isn't on PATH). With Docker, run_metal() pulls the image automatically.")
    sif_path.parent.mkdir(parents=True, exist_ok=True)
    log.info("Pulling docker://%s -> %s", container, sif_path)
    subprocess.run([exe, "pull", str(sif_path), f"docker://{container}"], check=True)
    return str(sif_path)


def metal_command(script_name, workdir, container=DEFAULT_CONTAINER, container_engine="auto",
                  engine_bin=None, metal_bin=None, container_metal="metal"):
    """Build the argv that runs METAL on script_name inside workdir (host path)."""
    if not container:
        return [find_metal(metal_bin), script_name]
    kind, exe = find_engine(container_engine, engine_bin)
    workdir = str(workdir)
    if kind == "docker":
        return [exe, "run", "--rm", "--platform", "linux/amd64",
                "-u", f"{os.getuid()}:{os.getgid()}",
                "-v", f"{workdir}:{workdir}", "-w", workdir,
                "--entrypoint", container_metal, container, script_name]
    image = container if container.endswith(".sif") or os.path.exists(container) \
        else f"docker://{container}"
    return [exe, "exec", "--cleanenv", "-B", workdir, "--pwd", workdir,
            image, container_metal, script_name]


def check_random_effects(tbl_path, analyze):
    """Warn if RANDOM was requested but the METAL build only wrote fixed-effects columns."""
    if "RANDOM" not in (analyze or "").upper():
        return []
    with open(tbl_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    extra = [c for c in header if c not in FIXED_TBL_COLS]
    if extra:
        log.info("Random-effects columns present: %s", ", ".join(extra))
    else:
        log.warning("ANALYZE RANDOM requested but %s has no random-effects columns -- the METAL "
                    "build that ran ignores RANDOM. Check that container= points at the random-effects "
                    "image and container_metal is the right binary inside it.", tbl_path)
    return extra


def orient_alleles(res, effect_allele="ALT"):
    """
    METAL always writes Allele1/Allele2 in alphabetical order (lowercase). Re-orient each row so
    Allele1 is the REF or ALT allele from the input (first file the marker appeared in), flipping
    the sign of every Effect* column, Freq1/MinFreq/MaxFreq and the Direction string to match.
    Rows where neither allele matches (e.g. multi-allelic ALT) keep METAL's order, ALIGNED=False.
    """
    res["Allele1"] = res["Allele1"].str.upper()
    res["Allele2"] = res["Allele2"].str.upper()
    if effect_allele == "metal":
        return res
    target = res[effect_allele].astype(str).str.upper()
    flip = (res["Allele2"] == target).to_numpy()
    res["ALIGNED"] = flip | (res["Allele1"] == target).to_numpy()

    for c in [c for c in res.columns if c.startswith("Effect")]:
        res.loc[flip, c] = -res.loc[flip, c]
    if "Freq1" in res:
        res.loc[flip, "Freq1"] = 1 - res.loc[flip, "Freq1"]
    if {"MinFreq", "MaxFreq"} <= set(res.columns):
        lo, hi = res.loc[flip, "MinFreq"].copy(), res.loc[flip, "MaxFreq"].copy()
        res.loc[flip, "MinFreq"], res.loc[flip, "MaxFreq"] = 1 - hi, 1 - lo
    if "Direction" in res:
        res.loc[flip, "Direction"] = res.loc[flip, "Direction"].str.translate(str.maketrans("+-", "-+"))
    a1 = res["Allele1"].copy()
    res.loc[flip, "Allele1"] = res.loc[flip, "Allele2"]
    res.loc[flip, "Allele2"] = a1[flip]
    n_bad = int((~res["ALIGNED"]).sum())
    log.info("Oriented to %s allele: %d flipped, %d could not be aligned (kept METAL order)",
             effect_allele, int(flip.sum()), n_bad)
    return res


def annotate_results(tbl_path, positions, out_path, effect_allele="ALT"):
    """Add CHR/POS/REF/ALT, orient alleles, add OR columns; sort by position."""
    res = pd.read_csv(tbl_path, sep="\t")
    pos = positions.drop_duplicates("MARKER")
    res = res.merge(pos, left_on="MarkerName", right_on="MARKER", how="left").drop(columns="MARKER")
    res = orient_alleles(res, effect_allele)
    res["OR"] = np.exp(res["Effect"])
    if "EffectARE" in res:
        res["OR_ARE"] = np.exp(res["EffectARE"])
    res["CHR_num"] = pd.to_numeric(res["CHR"].astype(str).str.replace("chr", "")
                                   .replace({"X": "23", "Y": "24", "XY": "25", "MT": "26", "M": "26"}),
                                   errors="coerce")
    res["POS"] = pd.to_numeric(res["POS"], errors="coerce")
    res = res.sort_values(["CHR_num", "POS"]).drop(columns="CHR_num")
    first = ["CHR", "POS", "MarkerName", "REF", "ALT", "Allele1", "Allele2"]
    cols = first + [c for c in res.columns if c not in first]
    res[cols].to_csv(out_path, sep="\t", index=False, compression="gzip")
    return res[cols]


# %% Main entry point
def run_metal(metal_input, metal_outdir="metal_results", metal_prefix="GWAS_META",
              metal_pheno_name=None, test="ADD", marker_id="ID", analyze="RANDOM",
              container=DEFAULT_CONTAINER, container_engine="auto", engine_bin=None,
              container_metal="metal", metal_bin=None, annotate=True, work_dir=None,
              effect_allele="ALT"):
    """
    metal_input       glob pattern, comma-separated patterns, or list of files
    metal_outdir      output directory
    metal_prefix      output prefix (phenotype appended when metal_pheno_name is set)
    metal_pheno_name  comma-separated phenotype names; files grouped by '.<pheno>.glm.logistic.hybrid'
                      and one METAL run per phenotype
    test              PLINK2 TEST row to meta-analyse (default 'ADD')
    marker_id         'ID' (use the ID column as-is) or 'chrpos' (build chr:pos:allele1:allele2)
    analyze           'RANDOM' (fixed + random effects, as in the .nf) or '' for fixed only
    container         Docker image ref or local .sif (default: GP2 progwas-pipeline image with
                      random-effects METAL). None -> run a local METAL binary (metal_bin).
    container_engine  'auto' (docker, then apptainer, then singularity), 'docker', 'podman',
                      'apptainer' or 'singularity'
    engine_bin        explicit path to the engine executable (useful in Jupyter on Biowulf)
    container_metal   METAL command inside the container (default 'metal' on the image PATH)
    metal_bin         path to a local METAL executable (only used when container=None)
    annotate          write <prefix>_annotated.tsv.gz with CHR/POS/OR added
    effect_allele     allele reported as Allele1 (effect allele) in the annotated file:
                      'ALT' (default), 'REF', or 'metal' (METAL's alphabetical order)
    work_dir          local directory where METAL actually runs (the only path mounted into the
                      container); results are copied to metal_outdir afterwards. Use this when
                      metal_outdir is on a mount the Docker daemon can't see (e.g. a gcsfuse bucket).
                      Default: run directly in metal_outdir.

    Returns {pheno: {"tbl": path, "info": path, "log": path, "stats": DataFrame, "results": DataFrame|None}}
    """
    outdir = Path(metal_outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    wd = Path(work_dir).expanduser().resolve() if work_dir else outdir
    (wd / "prepared").mkdir(parents=True, exist_ok=True)
    # Resolve the engine/binary up front so a missing engine fails before the (slow) prep step.
    metal_command("x", wd, container, container_engine, engine_bin, metal_bin, container_metal)
    log.info("METAL via %s", f"container {container}" if container else "local binary")

    files = resolve_inputs(metal_input)
    groups = group_by_pheno(files, metal_pheno_name)
    log.info("Found %d input file(s) in %d group(s)", len(files), len(groups))

    results = {}
    for pheno, group_files in groups.items():
        prefix = f"{metal_prefix}_{pheno}" if pheno else metal_prefix
        log.info("=== %s: %d file(s) ===", prefix, len(group_files))

        prepared, all_stats, positions = [], [], []
        for i, f in enumerate(group_files, 1):
            base = re.sub(r"\.(gz|zst)$", "", os.path.basename(f))
            out_path = wd / "prepared" / f"{prefix}_input{i:02d}_{base}.metal.tsv.gz"
            stats, pos = prepare_file(f, out_path, test=test, marker_id=marker_id)
            log.info("[A2_BUILD] file=%s A1_EQ_REF=%d A1_EQ_ALT=%d kept=%d (errcode=%d, NA=%d, noA2=%d)",
                     stats["file"], stats["A1_EQ_REF"], stats["A1_EQ_ALT"], stats["rows_out"],
                     stats["dropped_errcode"], stats["dropped_na"], stats["dropped_no_a2"])
            prepared.append(str(out_path))
            all_stats.append(stats)
            positions.append(pos)

        stats_df = pd.DataFrame(all_stats)
        log.info("[A2_BUILD] total A1_EQ_REF=%d A1_EQ_ALT=%d",
                 stats_df["A1_EQ_REF"].sum(), stats_df["A1_EQ_ALT"].sum())
        stats_df.to_csv(outdir / f"{prefix}_input_summary.tsv", sep="\t", index=False)

        script_path = wd / f"{prefix}_metal_script.txt"
        script_path.write_text(build_script(prefix, prepared, analyze))

        log_path = wd / f"{prefix}_command.log"
        cmd = metal_command(script_path.name, wd, container, container_engine,
                            engine_bin, metal_bin, container_metal)
        log.info("Running: %s", " ".join(cmd))
        with open(log_path, "w") as lf:
            lf.write("# " + " ".join(cmd) + "\n")
            lf.flush()
            proc = subprocess.run(cmd, cwd=wd, stdout=lf, stderr=subprocess.STDOUT, text=True)
        if proc.returncode != 0:
            tail = "".join(open(log_path).readlines()[-5:])
            hint = ""
            if proc.returncode == 125 and "mount source path" in tail:
                hint = ("\nDocker can't mount this directory (its daemon doesn't see the same filesystem, "
                        "e.g. a FUSE/bucket mount). Pass work_dir= on a local disk, e.g. "
                        "work_dir='~/metal_tmp'; results are copied to metal_outdir afterwards.")
            raise RuntimeError(f"METAL failed for {prefix} (exit {proc.returncode}); see {log_path}\n"
                               f"{tail}{hint}")

        if wd != outdir:
            (outdir / "prepared").mkdir(exist_ok=True)
            for f in list(wd.glob(f"{prefix}*")) + list((wd / "prepared").glob(f"{prefix}_*")):
                shutil.copy2(f, outdir / f.relative_to(wd))
            log.info("Copied %s outputs from %s to %s", prefix, wd, outdir)
            log_path = outdir / log_path.name

        tbl = outdir / f"{prefix}1.TBL"  # METAL appends '1' to the OUTFILE prefix
        if not tbl.exists():
            raise RuntimeError(f"Expected {tbl} not found; see {log_path}")
        log.info("METAL done: %s", tbl)
        check_random_effects(tbl, analyze)

        res_df = None
        if annotate:
            ann = outdir / f"{prefix}_annotated.tsv.gz"
            res_df = annotate_results(tbl, pd.concat(positions), ann, effect_allele)
            log.info("Annotated: %s (%d markers)", ann, len(res_df))

        results[pheno or prefix] = dict(tbl=str(tbl), info=str(tbl) + ".info", log=str(log_path),
                                        stats=stats_df, results=res_df)
    return results


# %% Example (edit and run this cell)
if __name__ == "__main__" and "get_ipython" not in globals():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metal_input", required=True)
    ap.add_argument("--metal_outdir", default="metal_results")
    ap.add_argument("--metal_prefix", default="GWAS_META")
    ap.add_argument("--metal_pheno_name", default=None)
    ap.add_argument("--test", default="ADD")
    ap.add_argument("--marker_id", default="ID", choices=["ID", "chrpos"])
    ap.add_argument("--analyze", default="RANDOM")
    ap.add_argument("--container", default=DEFAULT_CONTAINER,
                    help="Docker image or .sif; pass 'none' to use a local METAL")
    ap.add_argument("--container_engine", default="auto",
                    choices=["auto", "docker", "podman", "apptainer", "singularity"])
    ap.add_argument("--engine_bin", default=None)
    ap.add_argument("--container_metal", default="metal")
    ap.add_argument("--metal_bin", default=None)
    ap.add_argument("--effect_allele", default="ALT", choices=["ALT", "REF", "metal"])
    ap.add_argument("--work_dir", default=None)
    a = vars(ap.parse_args())
    if str(a["container"]).lower() == "none":
        a["container"] = None
    run_metal(**a)

# Jupyter example (Biowulf, Apptainer): pull the image once, then run
# sif = pull_container("/data/kuznetsovn2/containers/progwas-pipeline_1.0.0.sif")
# res = run_metal(
#     metal_input="/data/kuznetsovn2/gwas/*/*.glm.logistic.hybrid",
#     metal_outdir="/data/kuznetsovn2/metal_results",
#     metal_prefix="PD_META",
#     metal_pheno_name=None,          # e.g. "PHENO1,PHENO2" for one run per phenotype
#     marker_id="ID",                 # "chrpos" if IDs differ between cohorts
#     container=sif,                  # or DEFAULT_CONTAINER with Docker
#     # engine_bin="/path/to/apptainer",  # if apptainer isn't on the kernel's PATH
# )
# res["PD_META"]["results"].head()
