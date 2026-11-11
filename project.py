"""Shared helpers for the notebooks: subgroup info, logging, run.json checks.

    from project import get_info, get_all, tool, log, check_run_settings, update_run_json

    info = get_info("plants_plt", "solanaceae")      # one subgroup
    everything = get_all()                            # {group: {subgroup: info}} from config.yaml

Notebooks in guide_trees/ need the project root on the path first:

    import sys; sys.path.insert(0, "..")
"""

import json
from datetime import datetime
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent
CONFIG = yaml.safe_load((ROOT / "config.yaml").read_text())
SOURCE_PATH = (ROOT / CONFIG["source_path"]).resolve()


def tool(name: str) -> str:
    """Tool command from config.yaml, with '~' expanded."""
    return str(Path(CONFIG["tools"][name]).expanduser())


def get_info(group: str, subgroup: str) -> dict:
    """Everything about one subgroup: its paths, run names and <subgroup>.yaml settings.

    config: config.yaml `defaults:`, overridden by the subgroup's <subgroup>.yaml.
    runs: config.yaml `runs:`, overridden by the subgroup's own `runs:`.
    """
    gd = ROOT / group
    d = gd / subgroup
    yaml_file = d / f"{subgroup}.yaml"
    own = (yaml.safe_load(yaml_file.read_text()) or {}) if yaml_file.exists() else {}
    config = {
        **CONFIG.get("defaults", {}),
        **own,
    }  # subgroup values win over config.yaml defaults
    runs = {**CONFIG.get("runs", {}), **(own.get("runs") or {})}
    return {
        "group": group,
        "subgroup": subgroup,
        "runs": runs,
        "config": config,  # config.yaml defaults + <subgroup>.yaml
        "paths": {
            # --- group level (shared by every subgroup of the group)
            "atlas_source": SOURCE_PATH / group,  # the atlas folder 01 copies from
            "group_dir": gd,
            "group_tsv": gd / f"{group}.tsv",
            "group_atlas": gd / f"{group}_summary_atlas.tsv",
            "group_gdict": gd / f"{group}.gdict",
            "group_tree": gd / "tree.nwk",
            "group_lineage": gd / f"{group}_lineage.tsv",
            # --- subgroup level
            "dir": d,
            "yaml": yaml_file,
            "tsv": d / f"{subgroup}.tsv",
            "atlas": d / f"{subgroup}_summary_atlas.tsv",
            "qc_excluded": d / "qc_excluded_ANs.tsv",
            # 01: verbatim copies from the atlas, re-creatable
            "source": d / "source",
            "genomes": d / "source" / "genomes",
            "gff3": d / "source" / "gff3",
            "igr_by_genome": d / "source" / "igr_by_genome",
            # 02
            "gene_pairs_dir": d / "gene_pairs",
            "gene_pairs": d / "gene_pairs" / "gene_pairs.tsv",
            "gene_pairs_filtered": d / "gene_pairs" / "gene_pairs_filtered.tsv",
            "gene_pairs_flipped": d / "gene_pairs" / "gene_pairs_flipped.tsv",
            "igr_by_pair": d / "gene_pairs" / "igr_by_pair",
            # guide_trees/
            "species_tree": d / "species_tree",
            "guide_trees": d / "guide_trees",
            # 03a/03b/04
            "alignments": d / "alignments",
        },
    }


def get_all(groups: list[str] | None = None) -> dict:
    """{group: {subgroup: get_info(...)}} for every subgroup listed in config.yaml
    (or only the groups named in `groups`)."""
    listed = CONFIG["groups"]
    unknown = set(groups or []) - set(listed)
    if unknown:
        raise KeyError(f"Groups not in config.yaml: {sorted(unknown)}")
    return {
        g: {s: get_info(g, s) for s in subs}
        for g, subs in listed.items()
        if groups is None or g in groups
    }


def excluded_pairs(info: dict, listed=None) -> set[str]:
    """Gene pairs to leave out, from `exclude_pairs` in <subgroup>.yaml.

    Every step after 02_gene_pairs (05, 06, 07, 03, 04) skips these. If `listed` (the gene
    pairs of the step) is given, names not in it are reported, to catch typos.
    """
    exclude = {str(gp) for gp in (info["config"].get("exclude_pairs") or [])}
    if listed is not None:
        unknown = exclude - {str(gp) for gp in listed}
        if unknown:
            log(f"  note: exclude_pairs has gene pairs not in this step: {sorted(unknown)}")
    return exclude


def log(msg: str) -> None:
    """Print a timestamped message."""
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}")


def update_run_json(run_json: Path, **fields) -> dict:
    """Merge `fields` into run.json (created if missing) and return the full content."""
    info = json.loads(run_json.read_text()) if run_json.exists() else {}
    info.update(fields)
    run_json.parent.mkdir(parents=True, exist_ok=True)
    run_json.write_text(json.dumps(info, indent=2))
    return info


def check_run_settings(run_json: Path, settings: dict, overwrite: bool) -> None:
    """Stop if a run folder already exists with different settings (unless overwrite)."""
    if not run_json.exists():
        return
    old = json.loads(run_json.read_text())
    # keys the run never recorded are new settings for it (e.g. a later step), not changes
    changed = {k: (old[k], v) for k, v in settings.items() if k in old and old[k] != v}
    if changed and not overwrite:
        details = "; ".join(f"{k}: {a!r} -> {b!r}" for k, (a, b) in changed.items())
        raise RuntimeError(
            f"{run_json.parent} already exists with different settings ({details}). "
            "Use a new run name, or set OVERWRITE = True."
        )
