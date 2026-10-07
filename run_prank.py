#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import shutil

import polars as pl


# In[ ]:


CWD = Path().resolve()
BASE = CWD

print(f"Folder where groups (data) folder should be: {BASE}")


# In[ ]:


groups_map = {
    "plants_plt": ["solanaceae"],
    "plants_mit": ["solanaceae"],
}

RUN_NAME = "prank_2026-09-30"
RE_RUN = False
WORKERS = 32  # PRANK is single-threaded: one job per core
PRANK_BIN = str(Path("~/tools/prank-msa/src/prank").expanduser())

if not Path(PRANK_BIN).exists():
    raise FileNotFoundError(
        f"PRANK binary not found at {PRANK_BIN}. Please check the path."
    )


# In[ ]:


def build_cmd(in_fasta, tree_file, out_prefix):
    return [
        PRANK_BIN,
        f"-d={in_fasta}",
        f"-t={tree_file}",
        f"-o={out_prefix}",
        "-DNA",
        "+F",  # insertions are never matched later (recommended for DNA)
        "-once",  # use the given guide tree, don't re-estimate it
        "-shortnames",  # names are cut at the first space, so they match the tree tips
    ]


def run_prank(gp, gp_dir, tree_dir, run_dir, logs_dir, force_run):
    in_fasta = gp_dir / f"gp_{gp}.fasta"
    tree_file = tree_dir / f"gp_{gp}.nwk"
    out_prefix = run_dir / f"gp_{gp}"
    raw_out = run_dir / f"gp_{gp}.best.fas"  # PRANK appends .best.fas to the prefix
    out_fasta = run_dir / f"gp_{gp}_aligned.fasta"
    log_file = logs_dir / f"gp_{gp}_prank.log"

    if not in_fasta.exists():
        return gp, "missing input"
    if not tree_file.exists():
        return gp, "missing tree"
    if out_fasta.exists() and not force_run:
        return gp, "skipped"

    # PRANK creates its tmpdirprankmsa* folders in the current working directory,
    # so each job runs inside its own scratch folder (all other paths are absolute)
    tmp_dir = run_dir / "tmp" / gp
    tmp_dir.mkdir(parents=True, exist_ok=True)

    cmd = build_cmd(in_fasta, tree_file, out_prefix)
    res = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=tmp_dir
    )
    log_file.write_text(res.stdout)

    if res.returncode != 0 or not raw_out.exists():
        (logs_dir / f"gp_{gp}_prank_error.log").write_text(res.stdout)
        return gp, f"error (exit {res.returncode})"  # tmp_dir kept for debugging

    raw_out.rename(out_fasta)
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return gp, "ok"


# In[ ]:


for g, subgs in groups_map.items():
    for subg in subgs:
        subg_dir = BASE / g / subg
        print(f"Processing {g}/{subg}")

        gp_file = subg_dir / "gene_pairs_filtered.tsv"
        if not gp_file.exists():
            print(f"  no {gp_file.name}, skipping")
            continue

        gp_dir = subg_dir / "gene_pair_fastas"
        tree_dir = subg_dir / "gene_pair_trees"
        run_dir = subg_dir / "alignments" / RUN_NAME
        logs_dir = run_dir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)

        template = build_cmd("gp_{gp}.fasta", "gp_{gp}.nwk", "gp_{gp}")
        (run_dir / "params.txt").write_text(" ".join(template) + "\n")

        gp_df = pl.read_csv(gp_file, separator="\t")
        pairs = gp_df.sort(pl.col("igr_count") * pl.col("mean_len"), descending=True)[
            "gene_pair"
        ].to_list()

        with ThreadPoolExecutor(max_workers=WORKERS) as executor:
            tasks = [
                executor.submit(
                    run_prank, gp, gp_dir, tree_dir, run_dir, logs_dir, RE_RUN
                )
                for gp in pairs
            ]

            for future in as_completed(tasks):
                gp, status = future.result()
                if status.startswith("error"):
                    print(f"  [ERROR] {gp}: {status}")
                else:
                    print(f"  {gp}: {status}")

print("All done.")
