# sam-honours

Phylogeny-guided alignment of organelle **intergenic regions (IGRs)**. For each subgroup of
genomes (e.g. Solanaceae plastomes, `plants_plt/solanaceae`), IGRs are grouped by the genes that flank them
(**gene pairs**, e.g. `RBCL_ACCD_++`), and each gene pair is aligned with PRANK using a guide
tree derived from a CDS species tree. Alignment metrics are then computed per gene pair and
per sequence.

## Pipelines

The work is split into two sets of notebooks:

| Folder | What it does | How often |
|---|---|---|
| [`notebooks/`](notebooks) | The main analysis: get data, build gene pairs, align with PRANK, compute metrics | Every time the analysis changes |
| [`guide_trees/`](guide_trees) | Builds the species tree and, from it, one PRANK guide tree per gene pair | Once per subgroup (see its [README](guide_trees/README.md)) |

Run order for a new subgroup:

```
notebooks/01_get_data → notebooks/02_gene_pairs
  → guide_trees/01 … 07   (species tree, guide trees, scale fit; see guide_trees/README.md)
  → notebooks/03_prank → notebooks/04_aln_analysis
```

| Notebook | Input | Output (in `<group>/<subgroup>/`) |
|---|---|---|
| `01_get_data` | the atlas (`source_path`), NCBI Taxonomy | `<subgroup>.tsv`, `<subgroup>_summary_atlas.tsv`, `source/` |
| `02_gene_pairs` | subgroup atlas, `source/igr_by_genome/` | `gene_pairs/*.tsv`, `gene_pairs/igr_by_pair/` |
| `03_prank` | `gene_pairs/igr_by_pair/`, `guide_trees/<run>/` | `alignments/<prank run>/` |
| `04_aln_analysis` | an alignment run | `alignments/<run>/metrics_by_pair.tsv`, `metrics_by_seq.tsv` |

## Setup

**Python** ≥ 3.10 with `polars`, `numpy`, `biopython`, `ete4`, `pyyaml`, `python-dotenv` and `gdt`.

**Tools** (paths or commands set in `config.yaml` → `tools:`):
[PRANK](https://github.com/ariloytynoja/prank-msa),
[Clustal Omega](http://www.clustal.org/omega/),
[MAFFT](https://mafft.cbrc.jp/alignment/software/),
[IQ-TREE 3](https://iqtree.github.io/).

**NCBI** (used by `01_get_data` for lineages): a `.env` file at the repo root, not committed:

```
ENTREZ_EMAIL=you@example.org
ENTREZ_API_KEY=...        # optional, raises the request limit
```

**Input data**: the atlas folder set by `source_path` in `config.yaml`, laid out as

```
<source_path>/<group>/
  <group>.tsv  <group>_summary_atlas.tsv  <group>.gdt  tree.nwk
  <AN>/<AN>.fasta  <AN>/<AN>.gff3  <AN>/<AN>_intergenic.fasta
```

## Configuration

Settings live in two places; notebooks read them through `project.py`.

**`config.yaml`** (repo root): shared by everything.

| Key | Meaning |
|---|---|
| `source_path` | the atlas folder, relative to the repo root |
| `tools` | command or path for each external tool |
| `groups` | the subgroups that exist, `{group: [subgroup, ...]}` |
| `defaults` | default settings, overridable per subgroup (e.g. the gene-pair filter thresholds) |
| `runs` | default run name for each step, overridable per subgroup |

**`<group>/<subgroup>/<subgroup>.yaml`**: one per subgroup, next to its TSV. Anything set here
wins over `config.yaml` (`defaults:` keys at the top level, run names under `runs:`).

| Key | Used by | Meaning |
|---|---|---|
| `lineage_keywords` | `01_get_data` | keep genomes whose NCBI lineage contains any of these taxa |
| `min_genomes`, `min_genome_frac`, `min_mean_len` | `02_gene_pairs` | keep a gene pair if `mean_len > min_mean_len` and it is in ≥ `max(min_genomes, min_genome_frac × n)` genomes |
| `outgroup` | `guide_trees/04` | ANs the species tree is rooted on (must form a clade) |
| `guide_tree_scale`, `guide_tree_scale_source` | `guide_trees/05` | branch-length scale for the guide trees (set by hand, see `guide_trees/README.md`) |
| `exclude_pairs` | `03_prank` | gene pairs left out of the PRANK run |
| `runs` | all | per-subgroup run names |

Settings that belong to one run (thresholds, tool flags) stay as constants at the top of each
notebook and are recorded in that run's `run.json`.

### `project.py`

```python
import sys; sys.path.insert(0, "..")       # notebooks live one level down
from project import get_info, get_all, tool, log

info = get_info("plants_plt", "solanaceae")
info["paths"]["igr_by_pair"]   # every folder and file of the subgroup
info["runs"]["prank"]          # run names (config.yaml, overridden by <subgroup>.yaml)
info["config"]["outgroup"]     # defaults + <subgroup>.yaml

get_all()                       # {group: {subgroup: info}} for every subgroup in config.yaml
```

Each notebook starts with `ALL = get_all()`; replace it with a hand-picked dict to run only some
subgroups.

## Runs

Every step that runs a tool writes into a named run folder (`alignments/<run>/`,
`species_tree/mafft/<run>/`, ...) with a `run.json` recording its settings, inputs, tool
versions, status and the runs it was built on. This makes it possible to keep several versions
side by side (e.g. PRANK with different flags) and to trace any output back to its inputs.

- To try new settings, use a **new run name** (in `<subgroup>.yaml` → `runs:`).
- Re-running an existing run with **different settings** stops with an error, unless the
  notebook's `OVERWRITE = True`.
- Each step **refuses to start** unless the run it depends on has `status: complete`.
- `notebooks/03_prank` refuses guide trees whose scale was never fitted, unless `ALLOW_UNSCALED = True`.

## Adding a subgroup

1. Add it under `groups:` in `config.yaml`.
2. Create `<group>/<subgroup>/<subgroup>.yaml` with at least `lineage_keywords` (and
   `qc_excluded_ANs.tsv`, one `AN` column, if genomes must be excluded).
3. Run `01_get_data` and `02_gene_pairs`.
4. Add `outgroup` to the YAML and build the guide trees (`guide_trees/README.md`).
5. Run `03_prank` and `04_aln_analysis`.

## Layout

```
config.yaml             project settings
project.py              helpers used by every notebook
notebooks/              main pipeline
guide_trees/            notebooks for the species tree and guide trees
<group>/
  <group>.tsv, <group>_summary_atlas.tsv, tree.nwk     copied from the atlas
  <subgroup>/
    <subgroup>.yaml, <subgroup>.tsv, <subgroup>_summary_atlas.tsv, qc_excluded_ANs.tsv
    source/             genomes/, gff3/, igr_by_genome/   (copies from the atlas)
    gene_pairs/         gene_pairs*.tsv, igr_by_pair/
    species_tree/       cds/<run>/, mafft/<run>/, iqtree/<run>/
    guide_trees/<run>/  one guide tree per gene pair, scale/<clustalo run>/
    alignments/<run>/   one alignment per gene pair, run.json, metrics
```

**In git:** configuration, notebooks, every `run.json`, QC and summary tables, and the final
species trees. **Not in git** (re-created by the notebooks): copies from the atlas, FASTAs,
alignments, per-gene-pair trees, logs and per-sequence tables. See `.gitignore`.