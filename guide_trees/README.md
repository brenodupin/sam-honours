# guide_trees

Builds the guide trees PRANK needs: one tree per gene pair, with the gene pair's IGR IDs as
tips and branch lengths in **IGR** substitutions per site.

They come from a species tree inferred from plastid CDS, pruned to the genomes in each gene
pair and scaled up, because IGRs evolve faster than CDS. The scale is fitted on tree-free
(Clustal Omega) IGR alignments, so it does not depend on the trees it calibrates.

The species tree does not change once built, so these notebooks normally run **once per
subgroup**; the main analysis in [`notebooks/`](../notebooks) can then be rerun freely.

## Steps

All outputs are in `<group>/<subgroup>/`. Run names come from `runs:` in `config.yaml` /
`<subgroup>.yaml`.

| Notebook | Input | Output | Run name |
|---|---|---|---|
| `01_cds_tree` | `source/gff3/`, `<group>.gdict` | `species_tree/cds/<run>/gff3_raw.tsv`, `gff3_cds.tsv` | `cds` |
| `02_cds_fasta` | the above, `source/genomes/` | `species_tree/cds/<run>/fasta/`, `label_qc.tsv` | `cds` |
| `03_mafft` | a complete cds run | `species_tree/mafft/<run>/aln/`, `aln_qc.tsv` | `mafft` |
| `04_iq_tree` | a complete mafft run, `outgroup` | `species_tree/iqtree/<run>/tree.nwk`, `tree_support.nwk` | `iqtree` |
| `05_igr_tree` | species tree, `gene_pairs/`, `guide_tree_scale` | `guide_trees/<run>/gp_<pair>.nwk` | `guide_trees` |
| `06_clustalo` | `gene_pairs/igr_by_pair/` | `alignments/<run>/` | `clustalo` |
| `07_refine_scale` | guide trees + Clustal Omega alignments | `guide_trees/<run>/scale/<clustalo run>/` | (both) |

1. **CDS table**: every CDS in the GFF3s is traced up its `Parent` chain to its gene and
   labelled with the group's gdict (`LABEL_SYNONYMS` merges plastid label variants).
2. **CDS FASTAs**: each CDS is built from the genome (circular genomes and trans-spliced
   rps12 handled), one copy kept per genome and gene (the one closest to the gene's typical
   length), out-of-frame and short records dropped, and genes present in < `COVERAGE_MIN` of
   genomes left out.
3. **MAFFT**: one alignment per gene, with per-gene and per-sequence QC.
4. **IQ-TREE**: partitioned ML tree (ModelFinder + merging, UFBoot and SH-aLRT), rooted on the
   subgroup's `outgroup` afterwards (not with `-o`, which can crash ModelFinder when an
   outgroup taxon lacks some genes).
5. **Guide trees**: species tree × `guide_tree_scale`, pruned per gene pair, tips renamed to IGR
   IDs (genomes with several IGRs in a gene pair, e.g. IR copies, get sister tips).
6. **Clustal Omega**: tree-free IGR alignments, used only by step 7. Written next to the PRANK
   runs so `notebooks/04_aln_analysis` can also use them as a baseline.
7. **Scale fit**: for each gene pair, IQ-TREE fits one factor on its fixed guide tree; the
   weighted median (by variable sites) is the suggested `guide_tree_scale`.

## The two passes

The scale is set **by hand** in `<subgroup>.yaml`, so the guide trees are built twice:

```yaml
# pass 1: unscaled
guide_tree_scale: 1
guide_tree_scale_source: ""
runs:
  guide_trees: cds_first_pass
```

Run `05`, `06`, `07`. Step 7 prints the lines for pass 2:

```yaml
# pass 2: fitted scale
guide_tree_scale: 14.632
guide_tree_scale_source: guide_trees/cds_first_pass/scale/clustalo_2026-10-08 (weighted median)
runs:
  guide_trees: cds_refined
```

Run `05` again. `notebooks/03_prank` uses `runs.guide_trees` and refuses trees that were never
scaled (scale 1 with no `scale_source`).

## Run chain

Run order:

```
01_cds_tree → 02_cds_fasta → 03_mafft → 04_iq_tree → 05_igr_tree → 06_clustalo → 07_refine_scale → 05_igr_tree (scale from 07, set by hand)
```

Every step writes a `run.json` and refuses to start unless the runs it reads from are complete
(06 only needs `gene_pairs/igr_by_pair/`; 07 reads both 05 and 06). The final guide-tree run is
then read by `notebooks/03_prank`.

- 01 and 02 share one cds run: 01 sets `status: gff3_done`, 02 sets `complete`. Re-running 01
  resets it, so `03_mafft` refuses to run until 02 has run again.
- Each `run.json` records the run names it was built on, so any guide tree (and the PRANK
  alignment that used it) traces back to the CDS settings.
- Settings stay as constants in each notebook and are recorded in `run.json`; changing them
  for an existing run name stops with an error unless `OVERWRITE = True`.

## Notes

- `LABEL_SYNONYMS` in `01_cds_tree` are plastid labels; a mitochondrial group needs its own list.
- IQ-TREE results are only reproducible with the same `--seed` **and** thread count; set
  `THREADS` to a number rather than `"AUTO"` (which picks a count from a timing test).
- The 3 genomes rooted as outgroup for Solanaceae are two *Petunia* and one *Calibrachoa*;
  `04_iq_tree` reports whether they form a clade and whether midpoint rooting agrees.