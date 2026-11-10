"""Shared readers for the generated tables the enrichment pipeline consumes.

Both the GO and regulon workflows read the same DESeq2 contrast table and the
same diurnal-adjusted response table, and both depend on the same identifier
convention. Those contracts live here so a validation tightened for one
workflow cannot leave the other blind.

The tested-gene universe is whatever the supplied table contains: its size is a
consequence of the expression threshold and the DESeq2 tested set, derived in
``scripts/analysis_universe.R``, so it is reported rather than asserted against
a literal. The structural invariants that would actually indicate breakage --
one shared universe across every contrast and both durations, no duplicated
genes -- are enforced.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from _tables import read_table


DIRECTION_ORDER = {"Up": 0, "Down": 1}
TREATMENT_ORDER = {"DAWN": 0, "DUSK": 1}
RESPONSE_DURATIONS = {"3", "27"}
NON_CATEGORIES = {"not_called", "not_evaluated"}


def gene_level_id(identifier: str) -> str:
    """Strip a trailing transcript suffix (``Sobic.001G000100.1`` -> gene)."""
    stem, separator, suffix = identifier.rpartition(".")
    return stem if separator and suffix.isdigit() else identifier


def enrichment_status(fold_enrichment: float) -> str:
    if fold_enrichment > 1:
        return "Enriched"
    if fold_enrichment < 1:
        return "Depleted"
    return "Equal"


def read_deg_table(path: Path) -> tuple[set[str], dict[tuple[str, int, str], set[str]]]:
    """Return the tested-gene universe and the up/down foreground per contrast."""
    frame = read_table(path)
    missing = {"gene_id", "treatment", "time_h", "de_status"}.difference(frame.columns)
    if missing:
        raise ValueError(f"DEG table is missing columns: {sorted(missing)}")
    unexpected = set(frame["de_status"]).difference({"Up", "Down", "Not DEG"})
    if unexpected:
        raise ValueError(f"Unexpected DEG status: {sorted(unexpected)}")
    frame = frame.assign(gene=frame["gene_id"].map(gene_level_id))
    grouped = frame.groupby(["treatment", "time_h"], sort=False)["gene"]
    if grouped.ngroups != 14:
        raise ValueError(f"Expected 14 contrasts, found {grouped.ngroups}")
    if (grouped.nunique() != grouped.size()).any():
        raise ValueError("A gene appears more than once within a contrast")
    universes = grouped.apply(frozenset)
    if universes.nunique() != 1:
        raise ValueError("Every contrast must contain the same tested-gene universe")

    contrast_genes: dict[tuple[str, int, str], set[str]] = {}
    for treatment, time_h in universes.index:
        for direction in DIRECTION_ORDER:
            contrast_genes[(treatment, int(time_h), direction)] = set()
    called = frame[frame["de_status"] != "Not DEG"]
    for (treatment, time_h, direction), block in called.groupby(
        ["treatment", "time_h", "de_status"], sort=False
    ):
        contrast_genes[(treatment, int(time_h), direction)] = set(block["gene"])
    return set(universes.iloc[0]), contrast_genes


def read_response_table(path: Path) -> tuple[set[str], dict[tuple[int, str], set[str]]]:
    """Return the response universe and one foreground per called category."""
    frame = read_table(path)
    required = {
        "gene_id", "hours_after_onset", "called_category",
        "collection_phase_specific_category",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Response table is missing columns: {sorted(missing)}")
    frame = frame.assign(gene=frame["gene_id"].map(gene_level_id))
    grouped = frame.groupby("hours_after_onset", sort=True)["gene"]
    if (grouped.nunique() != grouped.size()).any():
        raise ValueError("A gene appears more than once at one duration")
    universes = grouped.apply(frozenset)
    if set(universes.index) != RESPONSE_DURATIONS:
        raise ValueError(
            f"Expected response durations {sorted(RESPONSE_DURATIONS)}; "
            f"found {sorted(universes.index)}"
        )
    if universes.nunique() != 1:
        raise ValueError("Response durations do not share one universe")

    foregrounds: dict[tuple[int, str], set[str]] = defaultdict(set)
    called = frame[frame["called_category"] != "not_called"]
    for (time, category), block in called.groupby(
        ["hours_after_onset", "called_category"], sort=False
    ):
        foregrounds[(int(time), category)].update(block["gene"])
    phase = frame[~frame["collection_phase_specific_category"].isin(NON_CATEGORIES)]
    if not phase.empty and set(phase["hours_after_onset"]) != {"27"}:
        raise ValueError("Collection-phase-specific calls are only defined at 27 h")
    for (time, category), block in phase.groupby(
        ["hours_after_onset", "collection_phase_specific_category"], sort=False
    ):
        foregrounds[(int(time), category)].update(block["gene"])
    return set(universes.iloc[0]), dict(foregrounds)


def read_tested_universe(path: Path) -> set[str]:
    return read_deg_table(path)[0]


def read_residual_universe(path: Path) -> set[str]:
    return read_response_table(path)[0]
