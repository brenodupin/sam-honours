#!/usr/bin/env python
# coding: utf-8

# In[1]:


import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import polars as pl
from Bio import SeqIO
from ete4 import Tree

def log(msg: str) -> None:
    """Print a timestamped message."""
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}")


# In[2]:


CWD = Path().resolve()
BASE = CWD.parent

log(f"Current working directory: {CWD}")
log(f"Folder where groups (data) folder should be: {BASE}")

groups_map = {
    "plants_plt": ["solanaceae"],
}

RUN_NAME = "first_pass"  # output goes to <subgroup>/cds_tree/iqtree/<RUN_NAME>/
MAFFT_RUN_NAME = "first_pass"  # alignments from <subgroup>/cds_tree/mafft/<MAFFT_RUN_NAME>/

# outgroup per subgroup: the tree is rooted on these taxa after IQ-TREE (they must form a clade).
# Not passed to IQ-TREE with -o, which can crash ModelFinder when an outgroup taxon lacks some genes.
OUTGROUPS = {
    "solanaceae": ["NC_040178.1", "NC_050782.1", "NC_050169.1"],  # Petunia x2, Calibrachoa
}

IQTREE = "iqtree3"
THREADS = "AUTO"  # IQ-TREE picks the best number of threads ...
THREADS_MAX = 16  # ... up to this many (only used with "AUTO"; or set THREADS to a number)
MODEL = "MFP+MERGE"  # ModelFinder per partition, merging similar partitions
RCLUSTER = 10  # % of partition pairs examined while merging (keeps MERGE fast)
UFBOOT = 1000  # ultrafast bootstrap replicates
ALRT = 1000  # SH-aLRT replicates
PARTITION_BY_CODON = False  # one partition per codon position; only for codon-aware alignments
SEED = 12345  # fixed seed, so a rerun gives the same tree
EXTRA_ARGS: list[str] = []  # anything else to pass to IQ-TREE
OVERWRITE = False  # True reruns IQ-TREE from scratch and allows changing the settings of an existing run

# flagging thresholds, nothing is removed
LOW_UFBOOT = 95
LOW_ALRT = 80
LONG_BRANCH_RATIO = 5  # terminal branch this many times longer than its sister clade's average depth

if shutil.which(IQTREE) is None:
    raise FileNotFoundError(f"'{IQTREE}' not found in PATH")
iqtree_version = subprocess.run([IQTREE, "--version"], capture_output=True, text=True).stdout.splitlines()[0].strip()
log(f"IQ-TREE: {iqtree_version}")


# In[3]:


def check_run_settings(run_json: Path, settings: dict) -> None:
    """Stop if RUN_NAME already exists with different settings (unless OVERWRITE)."""
    if not run_json.exists():
        return
    old = json.loads(run_json.read_text())
    changed = {k: (old.get(k), v) for k, v in settings.items() if old.get(k) != v}
    if changed and not OVERWRITE:
        details = "; ".join(f"{k}: {a!r} -> {b!r}" for k, (a, b) in changed.items())
        raise RuntimeError(
            f"{run_json.parent.name} already exists with different settings ({details}). "
            "Use a new RUN_NAME, or set OVERWRITE = True."
        )


def charset_name(label: str) -> str:
    """NEXUS-safe partition name (no hyphens)."""
    return label.replace("-", "_")


def write_partitions(path: Path, aln_files: dict[str, Path], by_codon: bool) -> None:
    """NEXUS partition file pointing at each gene alignment by absolute path."""
    lines = ["#nexus", "begin sets;"]
    for label, aln in aln_files.items():
        name = charset_name(label)
        if by_codon:
            for pos in (1, 2, 3):
                lines.append(f"    charset {name}_pos{pos} = {aln.resolve()}: {pos}-.\\3;")
        else:
            lines.append(f"    charset {name} = {aln.resolve()}: *;")
    lines.append("end;")
    path.write_text("\n".join(lines) + "\n")


def run_iqtree(cmd: list[str], log_file: Path) -> None:
    """Run IQ-TREE; full output goes to a log, key lines are echoed here."""
    keys = ("Best-fit", "BEST SCORE", "Log-likelihood of the tree", "Total wall-clock",
            "ERROR", "WARNING", "Analysis results written", "Merging", "Optimal log-likelihood")
    with open(log_file, "a") as log:
        log.write(f"\n# {datetime.now():%Y-%m-%d %H:%M:%S} | {' '.join(cmd)}\n")
        log.flush()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            log.write(line)
            if line.startswith(keys) or any(k in line for k in ("ERROR", "WARNING")):
                print(f"    {line.rstrip()}")
        proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"IQ-TREE failed (exit {proc.returncode}), see {log_file}")


def split_support(name: str | None) -> tuple[float | None, float | None]:
    """IQ-TREE labels internal nodes as 'SH-aLRT/UFBoot'."""
    if not name or "/" not in name:
        return None, None
    alrt, ufb = name.split("/")[:2]
    try:
        return float(alrt), float(ufb)
    except ValueError:
        return None, None


def root_on(tree: Tree, outgroup: list[str]) -> bool:
    """Root the tree on the outgroup clade; returns False if the outgroup is not monophyletic."""
    leaves = {leaf.name for leaf in tree.leaves()}
    # move the root away from the outgroup first, so their common ancestor is a proper node
    ingroup_leaf = next(n for n in sorted(leaves) if n not in outgroup)
    tree.set_outgroup(tree[ingroup_leaf])
    mrca = tree.common_ancestor(outgroup) if len(outgroup) > 1 else tree[outgroup[0]]
    if {leaf.name for leaf in mrca.leaves()} != set(outgroup):
        return False
    tree.set_outgroup(mrca)
    return True


def root_split(tree: Tree) -> frozenset[frozenset[str]]:
    """The two leaf sets on either side of the root."""
    return frozenset(frozenset(leaf.name for leaf in child.leaves()) for child in tree.children)


# In[ ]:


for g, subgs in groups_map.items():
    for subg in subgs:
        subg_dir = BASE / g / subg
        tree_dir = subg_dir / "cds_tree"
        mafft_dir = tree_dir / "mafft" / MAFFT_RUN_NAME
        run_dir = tree_dir / "iqtree" / RUN_NAME
        run_dir.mkdir(parents=True, exist_ok=True)
        run_json = run_dir / "run.json"
        prefix = run_dir / "iqtree"
        log(f"Processing {g}/{subg} | run: {RUN_NAME} | alignments: mafft/{MAFFT_RUN_NAME}")

        # --- the alignment run: must be complete; its mode decides what partitioning is allowed
        mafft_json = mafft_dir / "run.json"
        if not mafft_json.exists():
            raise FileNotFoundError(f"{mafft_json} not found; run 03_mafft with RUN_NAME = {MAFFT_RUN_NAME!r}")
        mafft_info = json.loads(mafft_json.read_text())
        if mafft_info.get("status") != "complete":
            raise RuntimeError(f"mafft/{MAFFT_RUN_NAME} is not complete (status: {mafft_info.get('status')})")
        if PARTITION_BY_CODON and mafft_info.get("mode") != "codon":
            raise ValueError(
                f"PARTITION_BY_CODON needs codon-aware alignments, mafft/{MAFFT_RUN_NAME} is mode={mafft_info.get('mode')!r}"
            )

        # --- alignments: exactly the genes listed in the MAFFT run, each one rectangular
        genes = mafft_info["genes"]
        aln_files, taxa_per_gene, aln_problems = {}, {}, []
        for label in genes:
            aln = mafft_dir / "aln" / f"{label}.fasta"
            if not aln.exists():
                aln_problems.append(f"{label}: missing {aln.name}")
                continue
            records = list(SeqIO.parse(aln, "fasta"))
            if len({len(r.seq) for r in records}) != 1:
                aln_problems.append(f"{label}: sequences of different lengths (not an alignment)")
                continue
            aln_files[label] = aln
            taxa_per_gene[label] = {r.id for r in records}
        if aln_problems:
            for p in aln_problems:
                log(f"    {p}")
            raise RuntimeError(f"Problems with alignments in mafft/{MAFFT_RUN_NAME}")
        all_taxa = set().union(*taxa_per_gene.values())

        outgroup = OUTGROUPS.get(subg, [])
        missing_out = [t for t in outgroup if t not in all_taxa]
        if not outgroup or missing_out:
            raise ValueError(f"OUTGROUPS[{subg!r}] is empty or has taxa not in the alignments: {missing_out}")
        log(f"  {len(aln_files)} genes | {len(all_taxa)} taxa | outgroup: {outgroup}")

        # --- settings that define this run; an existing run with other settings is protected
        settings = {
            "mafft_run_name": MAFFT_RUN_NAME,
            "model": MODEL,
            "rcluster": RCLUSTER,
            "ufboot": UFBOOT,
            "alrt": ALRT,
            "partition_by_codon": PARTITION_BY_CODON,
            "outgroup": outgroup,
            "seed": SEED,
            "extra_args": EXTRA_ARGS,
        }
        check_run_settings(run_json, settings)

        partition_file = run_dir / "partitions.nex"
        cmd = [
            IQTREE,
            "-p", str(partition_file),
            "-m", MODEL,
            "-B", str(UFBOOT),
            "--alrt", str(ALRT),
            "-T", str(THREADS),
            "--seed", str(SEED),
            "--prefix", str(prefix),
        ]
        if THREADS == "AUTO":
            cmd += ["--threads-max", str(THREADS_MAX)]
        if "MERGE" in MODEL:
            cmd += ["-rcluster", str(RCLUSTER)]
        cmd += EXTRA_ARGS

        old_info = json.loads(run_json.read_text()) if run_json.exists() else {}
        done = old_info.get("status") == "complete" and Path(f"{prefix}.treefile").exists()

        # --- run IQ-TREE (skipped when this run already finished)
        if done and not OVERWRITE:
            log(f"  IQ-TREE already finished for {RUN_NAME}, skipping (OVERWRITE = True to redo)")
            run_info = old_info
        else:
            if OVERWRITE:
                cmd.append("-redo")
            write_partitions(partition_file, aln_files, PARTITION_BY_CODON)
            run_info = {
                "run_name": RUN_NAME,
                "group": g,
                "subgroup": subg,
                "status": "running",
                "date_started": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
                "date_finished": None,
                **settings,
                "mafft_run_dir": str(mafft_dir.relative_to(subg_dir)),
                "mafft_run_dir_abs": str(mafft_dir.resolve()),
                "mafft_mode": mafft_info.get("mode"),
                "mafft_args": mafft_info.get("mafft_args"),
                "mafft_version": mafft_info.get("mafft_version"),
                "iqtree": IQTREE,
                "iqtree_path": shutil.which(IQTREE),
                "iqtree_version": iqtree_version,
                "threads": THREADS,
                "threads_max": THREADS_MAX,
                "command": " ".join(cmd),
                "partition_file": str(partition_file.relative_to(subg_dir)),
                "partition_file_abs": str(partition_file.resolve()),
                "n_genes": len(aln_files),
                "n_taxa": len(all_taxa),
            }
            run_json.write_text(json.dumps(run_info, indent=2))
            log(f"  running: {' '.join(cmd)}")
            run_iqtree(cmd, run_dir / "iqtree_stdout.log")

        # --- root, write the final trees
        ml_tree = Tree(Path(f"{prefix}.treefile").read_text(), parser=1)
        tree = ml_tree.copy()
        outgroup_ok = root_on(tree, outgroup)
        if not outgroup_ok:
            log(f"  WARNING: outgroup {outgroup} is not monophyletic in the ML tree; tree left as IQ-TREE wrote it")

        midpoint = ml_tree.copy()
        midpoint.set_outgroup(midpoint.get_midpoint_outgroup())
        midpoint_agrees = outgroup_ok and root_split(midpoint) == root_split(tree)

        tree.write(outfile=str(run_dir / "tree_support.nwk"), parser=1)  # with SH-aLRT/UFBoot labels
        tree.write(outfile=str(run_dir / "tree.nwk"), parser=5)  # branch lengths only, for PRANK

        # --- checks
        leaves = list(tree.leaves())
        leaf_names = {leaf.name for leaf in leaves}
        internal = [n for n in tree.traverse() if not n.is_leaf and not n.is_root]
        supports = [split_support(n.name) for n in internal]
        alrt_vals = [a for a, _ in supports if a is not None]
        ufb_vals = [u for _, u in supports if u is not None]
        bifurcating = all(len(n.children) == 2 for n in tree.traverse() if not n.is_leaf)

        rows = []
        for leaf in leaves:
            parent = leaf.up
            sisters = [c for c in parent.children if c is not leaf]
            sister_tips = [t for s in sisters for t in s.leaves()]
            sister_avg = (
                sum(parent.get_distance(parent, t) for t in sister_tips) / len(sister_tips) if sister_tips else None
            )
            rows.append(
                {
                    "AN": leaf.name,
                    "terminal_branch": leaf.dist,
                    "sister_avg_depth": sister_avg,
                    "ratio_to_sister": (leaf.dist / sister_avg) if sister_avg else None,
                    "n_genes": sum(leaf.name in taxa for taxa in taxa_per_gene.values()),
                    "in_outgroup": leaf.name in outgroup,
                }
            )
        tree_qc = pl.DataFrame(rows)
        # long branch: longer than LONG_BRANCH_RATIO x the sister clade's average depth, with the
        # tree's median terminal branch as a floor (sisters that are identical have depth ~0)
        floor = float(tree_qc["terminal_branch"].median())
        tree_qc = tree_qc.with_columns(
            (
                pl.col("terminal_branch")
                > LONG_BRANCH_RATIO * pl.max_horizontal(pl.col("sister_avg_depth").fill_null(0.0), pl.lit(floor))
            ).alias("long_branch")
        ).sort("terminal_branch", descending=True)
        tree_qc.write_csv(run_dir / "tree_qc.tsv", separator="\t")

        n_zero = tree_qc.filter(pl.col("terminal_branch") < 1e-5).height
        long_tips = tree_qc.filter("long_branch")
        max_depth = max(tree.get_distance(tree, leaf) for leaf in leaves)

        run_info.update(
            status="complete",
            date_finished=run_info.get("date_finished") or f"{datetime.now():%Y-%m-%d %H:%M:%S}",
            outgroup_monophyletic=outgroup_ok,
            midpoint_root_agrees=midpoint_agrees,
            tree=str((run_dir / "tree.nwk").relative_to(subg_dir)),
            tree_abs=str((run_dir / "tree.nwk").resolve()),
        )
        run_json.write_text(json.dumps(run_info, indent=2))

        # --- summary
        log(f"  trees -> {run_dir.relative_to(subg_dir)}/tree.nwk (PRANK), tree_support.nwk (with support)")
        log(f"  taxa in tree: {len(leaf_names)} / {len(all_taxa)} in alignments | bifurcating: {bifurcating}")
        if leaf_names != all_taxa:
            log(f"  WARNING: missing from tree: {sorted(all_taxa - leaf_names)} | extra: {sorted(leaf_names - all_taxa)}")
        log(f"  outgroup monophyletic: {outgroup_ok} | midpoint root agrees: {midpoint_agrees}")
        log(f"  max root-to-tip: {max_depth:.4f} subst/site | tree length: {sum(n.dist or 0 for n in tree.traverse()):.4f}")
        log(
            f"  internal nodes: {len(internal)} | UFBoot < {LOW_UFBOOT}: {sum(u < LOW_UFBOOT for u in ufb_vals)}"
            f" | SH-aLRT < {LOW_ALRT}: {sum(a < LOW_ALRT for a in alrt_vals)}"
        )
        log(f"  zero-length terminal branches (< 1e-5, identical or nearly identical sequences): {n_zero}")
        log(
            f"  long terminal branches (> {LONG_BRANCH_RATIO}x max(sister depth, median terminal {floor:.5f})): {long_tips.height}"
        )
        if long_tips.height:
            log(long_tips.select("AN", "terminal_branch", "sister_avg_depth", "ratio_to_sister", "n_genes"))

log(f"Done with all groups/subgroups for run {RUN_NAME}")