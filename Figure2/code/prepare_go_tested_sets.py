#!/usr/bin/env python3
"""Prepare canonical GO sets for DEG or diurnal-adjusted response enrichment.

Direct Sorghum bicolor GO annotations are propagated through ``is_a`` and
``part_of`` relationships using the matching ontology release. Terms are
restricted to 20--3,000 genes in the selected tested-gene universe. An eligible
ancestor is removed when an eligible descendant contains at least 90% of the
ancestor's members.

Propagation and the ancestor traversal used by the redundancy rule both come
from GOATOOLS and are driven by the same ``PROPAGATED_RELATIONSHIPS`` set, so
the two steps cannot disagree about what counts as an ancestor.

The default ``deseq2`` mode supplies conventional DEG enrichment.
``residual-response`` rebuilds memberships in the rhythm-quantifiable /
DESeq2-tested universe used for the early/late diurnal-adjusted response calls.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

from goatools.anno.update_association import update_association
from goatools.godag.go_tasks import get_go2ancestors
from goatools.obo_parser import GODag

from _datasets import gene_level_id, read_residual_universe, read_tested_universe
from _tables import read_table, write_table


MIN_TERM_SIZE = 20
MAX_TERM_SIZE = 3_000
REDUNDANT_ANCESTOR_THRESHOLD = 0.90
PROPAGATED_RELATIONSHIPS = {"part_of"}
ONTOLOGY_ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
NAMESPACE_ABBREVIATIONS = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}


def repository_paths(mode: str) -> dict[str, Path]:
    analysis_dir = Path(__file__).resolve().parents[1]
    residual = mode == "residual-response"
    return {
        "universe": analysis_dir / "data/generatable" / (
            "diurnal_adjusted_response_categories.csv.gz" if residual
            else "deseq2_timepoint_contrasts.csv.gz"
        ),
        "associations": analysis_dir
        / "data/provided/sorghum_v5.1_go_associations.tsv.gz",
        "obo": analysis_dir / "data/externalsourcerequired/go_2026-06-19.obo",
        "output": analysis_dir / "data/generatable" / (
            "go_tested_sets_residual_response.tsv.gz" if residual
            else "go_tested_sets.tsv.gz"
        ),
    }


def read_go_associations(path: Path, dag: GODag) -> dict[str, set[str]]:
    frame = read_table(path, separator="\t")
    if list(frame.columns) != ["target_gene_id", "go_term_ids"]:
        raise ValueError("GO association columns must be ['target_gene_id', 'go_term_ids']")
    if frame["target_gene_id"].duplicated().any():
        duplicated = frame.loc[frame["target_gene_id"].duplicated(), "target_gene_id"].iloc[0]
        raise ValueError(f"Duplicate GO association row for {duplicated}")
    associations: dict[str, set[str]] = {}
    for gene_id, go_field in zip(frame["target_gene_id"], frame["go_term_ids"]):
        go_ids = set(go_field.split(";"))
        unknown = go_ids.difference(dag)
        if unknown:
            raise ValueError(f"Unknown GO IDs for {gene_id}: {sorted(unknown)}")
        obsolete = sorted(go_id for go_id in go_ids if dag[go_id].is_obsolete)
        if obsolete:
            raise ValueError(f"Obsolete GO IDs for {gene_id}: {obsolete}")
        associations[gene_id] = go_ids
    return associations


def prepare_sets(
    associations: dict[str, set[str]], universe: set[str], dag: GODag
) -> tuple[dict[str, set[str]], set[str]]:
    update_association(
        associations,
        dag,
        relationships=PROPAGATED_RELATIONSHIPS,
        prt=None,
    )
    genes_by_term: dict[str, set[str]] = defaultdict(set)
    for gene_id in universe:
        for go_id in associations.get(gene_id, set()):
            if go_id in dag and not dag[go_id].is_obsolete:
                genes_by_term[go_id].add(gene_id)

    candidates = {
        go_id: genes
        for go_id, genes in genes_by_term.items()
        if go_id not in ONTOLOGY_ROOTS
        and MIN_TERM_SIZE <= len(genes) <= MAX_TERM_SIZE
    }
    ancestors_by_term = get_go2ancestors(
        {dag[go_id] for go_id in candidates}, PROPAGATED_RELATIONSHIPS, prt=None
    )
    removed_ancestors: set[str] = set()
    for child_id, child_genes in candidates.items():
        for ancestor_id in ancestors_by_term.get(child_id, set()):
            if ancestor_id not in candidates:
                continue
            ancestor_genes = candidates[ancestor_id]
            if not child_genes.issubset(ancestor_genes):
                raise ValueError(
                    f"Propagated child {child_id} is not a subset of ancestor {ancestor_id}"
                )
            if len(child_genes) / len(ancestor_genes) >= REDUNDANT_ANCESTOR_THRESHOLD:
                removed_ancestors.add(ancestor_id)

    tested_terms = set(candidates).difference(removed_ancestors)
    return {go_id: candidates[go_id] for go_id in tested_terms}, removed_ancestors


def build_output_frame(sets: dict[str, set[str]], dag: GODag) -> pd.DataFrame:
    rows = []
    for go_id in sorted(sets):
        genes = sorted(sets[go_id])
        rows.append(
            {
                "go_id": go_id,
                "namespace": NAMESPACE_ABBREVIATIONS[dag[go_id].namespace],
                "go_name": dag[go_id].name,
                "population_count": len(genes),
                "gene_ids": ";".join(genes),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("deseq2", "residual-response"),
        default="deseq2",
        help="tested-gene universe and output family (default: %(default)s)",
    )
    args = parser.parse_args()
    paths = repository_paths(args.mode)
    missing = [str(path) for name, path in paths.items() if name != "output" and not path.exists()]
    if missing:
        raise SystemExit("Missing required input files: " + ", ".join(missing))

    universe = (
        read_residual_universe(paths["universe"])
        if args.mode == "residual-response"
        else read_tested_universe(paths["universe"])
    )
    dag = GODag(
        str(paths["obo"]),
        optional_attrs={"relationship"},
        load_obsolete=False,
        prt=None,
    )
    associations = read_go_associations(paths["associations"], dag)
    tested_sets, removed_ancestors = prepare_sets(associations, universe, dag)
    write_table(paths["output"], build_output_frame(tested_sets, dag))

    namespace_counts: dict[str, int] = defaultdict(int)
    for go_id in tested_sets:
        namespace_counts[NAMESPACE_ABBREVIATIONS[dag[go_id].namespace]] += 1
    print(f"Tested-gene universe: {len(universe):,}")
    print(
        "GO terms retained: "
        f"{len(tested_sets):,} "
        f"(BP {namespace_counts['BP']:,}; MF {namespace_counts['MF']:,}; "
        f"CC {namespace_counts['CC']:,})"
    )
    print(f"Redundant ancestor terms excluded: {len(removed_ancestors):,}")
    print(f"Wrote {paths['output']}")


if __name__ == "__main__":
    main()
