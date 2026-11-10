#!/usr/bin/env python3
"""Build every panel's plotted data from data/source/ and check it against the
source tables. Writes data/panel_<letter>/. Plotting scripts read only these."""
import numpy as np
import pandas as pd
from scipy.optimize import brentq, curve_fit

from figure1_style import DATA, SOURCE

WORKBOOK = SOURCE / "Ion_leakage_relative_injury_summary_v2.xlsx"
DESEQ = SOURCE / "deseq2_timepoint_contrasts.csv.gz"
OVERLAP = SOURCE / "dawndusk_overlap_summary.csv"
GO = SOURCE / "go_enrichment_timepoint_contrasts.tsv.gz"

FREEZING_TABS = {"A": "D1 Sorg dawn ctl v cold1", "B": "D1 Sorg dusk ctl v cold1"}
TEMP, INJURY = "Temperature (°C)", "Relative injury (%)"
DAWN_TIMES = [3, 7, 11, 15, 19, 23, 27]
DUSK_TIMES = [15, 19, 23, 27, 31, 35, 39]
STATES = ["Down", "Not DEG", "Up"]
DUSK_OFFSET = 12
HOURS = (11, 15)
LIM = 8
THRESH = -np.log10(0.05)

# Panels I and J highlight three groups. Every other BP term is plotted in gray.
LIGHT_HARVESTING = "photosynthesis, light harvesting"
PHOTOSYNTHESIS = ["photosynthesis", "photosynthesis, light reaction",
                  "photosynthetic electron transport chain",
                  "photosynthetic electron transport in photosystem I", "photosystem II assembly",
                  "regulation of photosynthesis", "regulation of photosynthesis, light reaction",
                  "photosynthesis, dark reaction", "carbon fixation"]


def is_ribosomal(name):
    """Named ribosome or rRNA terms, plus ribonucleoprotein complex biogenesis."""
    lower = name.lower()
    if lower.startswith("nonribosomal"):
        return False
    return "ribosom" in lower or "rrna" in lower or lower == "ribonucleoprotein complex biogenesis"


def plot_group(name):
    if name == LIGHT_HARVESTING:
        return "Light harvesting"
    if name in PHOTOSYNTHESIS:
        return "Photosynthesis"
    return "Ribosomal" if is_ribosomal(name) else "other BP term"


def out_dir(panel):
    path = DATA / f"panel_{panel}"
    path.mkdir(parents=True, exist_ok=True)
    return path


# ---------------------------------------------------------------- A, B -----
def logistic(x, top, bottom, mid, k):
    return bottom + (top - bottom) / (1 + np.exp(k * (x - mid)))


def prepare_freezing():
    raw = pd.read_excel(WORKBOOK, sheet_name="Raw data")
    stats = pd.read_excel(WORKBOOK, sheet_name="Statistics")
    params = pd.read_excel(WORKBOOK, sheet_name="Set parameters & LT50")
    lt50_col = [c for c in params.columns if c.startswith("LT50")][0]
    pooled_cols = [c for c in stats.columns if c.startswith("Pooled-curve")]
    for panel, tab in FREEZING_TABS.items():
        row = stats[stats["Comparison tab"] == tab].iloc[0]
        groups = [row["Reference group"], row["Compared group"]]
        expected = [float(row[c]) for c in pooled_cols]
        table = pd.read_excel(WORKBOOK, sheet_name=tab, header=24, nrows=10)
        points, curves, lt50 = [], [], []
        for i, group in enumerate(groups):
            d = raw[(raw["Group"] == group) & raw[INJURY].notna()]
            p, _ = curve_fit(logistic, d[TEMP], d[INJURY], p0=[100, 5, -2.2, 5], maxfev=20000)
            t50 = brentq(lambda x: logistic(x, *p) - 50, -5, 0)
            assert abs(t50 - expected[i]) < 0.002, (group, t50, expected[i])
            curves.append(dict(zip(["group", "top", "bottom", "midpoint", "k", "temp_50pct_injury"],
                                   [group, *p, t50])))
            s = (d.assign(temperature=d[TEMP].replace({-2.25: -2.2, -2.75: -2.7}))
                 .groupby("temperature")[INJURY].agg(mean="mean", sem="sem", n="size")
                 .reset_index().sort_values("temperature", ascending=False))
            shown = table.set_index(table.columns[0])[f"{group}: mean (%)"]
            assert np.allclose(s["mean"], shown.reindex(s["temperature"]), atol=0.006), group
            s.insert(0, "group", group)
            points.append(s)
            v = params[(params["Group"] == group) & params[lt50_col].notna()]
            lt50.append(pd.DataFrame({"group": group, "set": v["Set"].astype(int), "lt50_C": v[lt50_col]}))
        out = out_dir(panel)
        pd.concat(points).to_csv(out / "points_mean_sem.csv", index=False)
        pd.DataFrame(curves).to_csv(out / "curve_fit_parameters.csv", index=False)
        pd.concat(lt50).to_csv(out / "inset_lt50_per_set.csv", index=False)
        row.to_frame().T.to_csv(out / "inset_statistics.csv", index=False)
        raw[raw["Group"].isin(groups)].to_csv(out / "raw_samples.csv", index=False)


# ---------------------------------------------------------------- C, D -----
def prepare_alluvial(deg):
    n_all = deg["gene_id"].nunique()
    for panel, treatment, times in (("C", "DAWN", DAWN_TIMES), ("D", "DUSK", DUSK_TIMES)):
        d = deg[deg["treatment"] == treatment]
        wide = d.pivot(index="gene_id", columns="time_h", values="de_status")[times]
        assert len(wide) == n_all and not wide.isna().any().any()
        # Every tested gene is kept, including genes that are Not DEG at every time.
        strata = []
        for t in times:
            for state in STATES:
                strata.append({"treatment": treatment, "time_h": t, "state": state,
                               "genes": int((wide[t] == state).sum())})
        flows = []
        for left, right in zip(times, times[1:]):
            for dest in STATES:
                for source in STATES:
                    flows.append({"treatment": treatment, "from_time_h": left, "to_time_h": right,
                                  "from_state": source, "to_state": dest,
                                  "genes": int(((wide[left] == source) & (wide[right] == dest)).sum())})
        strata, flows = pd.DataFrame(strata), pd.DataFrame(flows)
        assert (strata.groupby("time_h")["genes"].sum() == n_all).all()
        assert (flows.groupby("from_time_h")["genes"].sum() == n_all).all()
        out = out_dir(panel)
        wide.rename(columns=lambda t: f"T{t}").reset_index().to_csv(out / "gene_states_all_tested_genes.csv", index=False)
        strata.to_csv(out / "strata_counts.csv", index=False)
        flows.to_csv(out / "flow_counts.csv", index=False)
        never = int((wide == "Not DEG").all(axis=1).sum())
        print(f"{panel} {treatment}: {n_all} genes, {never} Not DEG at every time point")


# ---------------------------------------------------------------- E, F -----
def prepare_lines(deg):
    o = pd.read_csv(OVERLAP)
    o = o[o["alignment"] == "matched_duration"].sort_values(["direction", "x_value"])
    counts = (deg[deg["de_status"] != "Not DEG"]
              .groupby(["treatment", "time_h", "de_status"]).size())
    for _, r in o.iterrows():
        assert counts["DAWN", r["dawn_time_h"], r["direction"]] == r["dawn_deg"]
        assert counts["DUSK", r["dusk_time_h"], r["direction"]] == r["dusk_deg"]
    o[["hours_after_onset", "dawn_time_h", "dusk_time_h", "direction", "dawn_deg", "dusk_deg"]].to_csv(
        out_dir("E") / "deg_counts_matched_duration.csv", index=False)
    o[["hours_after_onset", "dawn_time_h", "dusk_time_h", "direction", "tested_genes", "dawn_deg",
       "dusk_deg", "observed_overlap", "dodecile_expected_overlap",
       "dodecile_log2_observed_expected", "dodecile_bh_adjusted_p"]].to_csv(
        out_dir("F") / "overlap_observed_expected.csv", index=False)


# ---------------------------------------------------------------- G, H -----
def prepare_scatter(deg):
    expected = {11: {"DAWN down / DUSK up": 1723, "DAWN up / DUSK down": 1689, "both up": 653, "both down": 137},
                15: {"DAWN down / DUSK up": 1971, "DAWN up / DUSK down": 1817, "both up": 768, "both down": 491}}
    for panel, h in (("G", 11), ("H", 15)):
        a = deg[(deg["treatment"] == "DAWN") & (deg["time_h"] == h)].set_index("gene_id")
        b = deg[(deg["treatment"] == "DUSK") & (deg["time_h"] == h + DUSK_OFFSET)].set_index("gene_id")
        p = pd.DataFrame({"DAWN_log2FC": a["log2_fold_change"], "DUSK_log2FC": b["log2_fold_change"],
                          "DAWN_padj": a["adjusted_p_value"], "DUSK_padj": b["adjusted_p_value"],
                          "DAWN_status": a["de_status"], "DUSK_status": b["de_status"]}).dropna(
                              subset=["DAWN_log2FC", "DUSK_log2FC"])
        s, t = p["DAWN_status"], p["DUSK_status"]
        cat = np.select([(s == "Up") & (t == "Down"), (s == "Down") & (t == "Up"), (s == "Up") & (t == "Up"),
                         (s == "Down") & (t == "Down"), (s != "Not DEG") ^ (t != "Not DEG")],
                        ["DAWN up / DUSK down", "DAWN down / DUSK up", "both up", "both down",
                         "DEG in one treatment"], "not DEG")
        p["category"] = cat
        n = p["category"].value_counts()
        for name, value in expected[h].items():
            assert n[name] == value, (h, name, n[name])
        out = out_dir(panel)
        p.reset_index().to_csv(out / "gene_log2fc_and_category.csv", index=False)
        n.rename_axis("category").rename("genes").reset_index().to_csv(out / "category_counts.csv", index=False)


# ---------------------------------------------------------------- I, J -----
def prepare_go():
    go = pd.read_csv(GO, sep="\t")
    bp = go[go["namespace"] == "BP"]
    enriched = bp[bp["enrichment_status"] == "Enriched"].assign(
        score=lambda d: -np.log10(d["bh_adjusted_p_value"]))
    ids = bp["go_id"].unique()
    names = bp.drop_duplicates("go_id").set_index("go_id")["go_name"]

    def scores(treatment, time_h):
        d = enriched[(enriched["treatment"] == treatment) & (enriched["time_h"] == time_h)]
        up = d[d["direction"] == "Up"].set_index("go_id")["score"].reindex(ids).fillna(0)
        down = d[d["direction"] == "Down"].set_index("go_id")["score"].reindex(ids).fillna(0)
        return pd.Series(np.where(up >= down, up, -down), index=ids)

    def sig(treatment, time_h, direction):
        d = enriched[(enriched["treatment"] == treatment) & (enriched["time_h"] == time_h)
                     & (enriched["direction"] == direction) & enriched["significant"]]
        return set(d["go_id"])

    def opposite(h, dawn_dir, dusk_dir):
        return sig("DAWN", h, dawn_dir) & sig("DUSK", h + DUSK_OFFSET, dusk_dir)

    both = {"DAWN up / DUSK down": opposite(11, "Up", "Down") & opposite(15, "Up", "Down"),
            "DAWN down / DUSK up": opposite(11, "Down", "Up") & opposite(15, "Down", "Up")}
    assert {k: len(v) for k, v in both.items()} == {"DAWN up / DUSK down": 27, "DAWN down / DUSK up": 22}

    group = names.reindex(ids).map(plot_group)
    assert group.value_counts().to_dict().get("Light harvesting") == 1
    assert (group == "Photosynthesis").sum() == len(PHOTOSYNTHESIS)
    assert (group == "Ribosomal").sum() == 19
    assert len(ids) == 1558

    # One shared axis limit keeps panels I and J directly comparable.
    all_scores = {h: (scores("DAWN", h), scores("DUSK", h + DUSK_OFFSET)) for h in HOURS}
    limit = 1.1 * max(np.abs(v).max() for pair in all_scores.values() for v in pair)
    for panel, h, known in (("I", 11, (91, 9)), ("J", 15, (78, 0))):
        table = pd.DataFrame({"go_id": ids, "go_name": names.reindex(ids).values,
                              "DAWN_signed_log10_padj": all_scores[h][0].values,
                              "DUSK_signed_log10_padj": all_scores[h][1].values,
                              "plot_group": group.values})
        # significant = BH-adjusted P < 0.05 in that DEG set
        table["significant_DAWN"] = table["DAWN_signed_log10_padj"].abs() > THRESH
        table["significant_DUSK"] = table["DUSK_signed_log10_padj"].abs() > THRESH
        table["orientation_at_both_times"] = [
            "DAWN up / DUSK down" if g in both["DAWN up / DUSK down"]
            else "DAWN down / DUSK up" if g in both["DAWN down / DUSK up"] else "" for g in ids]
        opp = opposite(h, "Up", "Down") | opposite(h, "Down", "Up")
        same = opposite(h, "Up", "Up") | opposite(h, "Down", "Down")
        assert (len(opp), len(same)) == known, (h, len(opp), len(same))
        out = out_dir(panel)
        table.to_csv(out / "go_term_scores.csv", index=False)
        table[table["plot_group"] != "other BP term"].sort_values(["plot_group", "go_name"]).to_csv(
            out / "highlighted_terms.csv", index=False)
        pd.DataFrame([{"hours_after_onset": h, "opposite_direction_terms": len(opp),
                       "same_direction_terms": len(same), "axis_limit": limit}]).to_csv(
            out / "term_counts_and_axis_limit.csv", index=False)


def prepare_selected_go():
    """Panels I and J: the three chosen GO terms, 11 and 15 h after onset, up and down DEGs."""
    terms = ["ribosome biogenesis", "photosynthesis, light harvesting", "cellular response to cold"]
    go = pd.read_csv(GO, sep="\t")
    go = go[(go["namespace"] == "BP") & go["go_name"].isin(terms)]
    parts = []
    for h in HOURS:
        for treatment, time_h in (("DAWN", h), ("DUSK", h + DUSK_OFFSET)):
            d = go[(go["treatment"] == treatment) & (go["time_h"] == time_h)].copy()
            d.insert(0, "hours_after_onset", h)
            parts.append(d)
    d = pd.concat(parts)
    assert len(d) == 3 * 2 * 2 * 2, "expected one row per term, time, treatment and direction"
    d["percent_of_term_genes"] = 100 * d["observed_overlap"] / d["go_population_count"]
    d["percent_expected"] = 100 * d["expected_overlap"] / d["go_population_count"]
    p = d["bh_adjusted_p_value"]
    d["stars"] = np.select([p < 0.001, p < 0.01, p < 0.05], ["***", "**", "*"], default="")
    d.loc[d["enrichment_status"] != "Enriched", "stars"] = ""  # stars mark enrichment only
    assert d["stars"].ne("").sum() == 6
    cols = ["hours_after_onset", "go_id", "go_name", "treatment", "time_h", "direction", "go_population_count",
            "observed_overlap", "expected_overlap", "percent_of_term_genes", "percent_expected",
            "fold_enrichment", "enrichment_status", "bh_adjusted_p_value", "stars"]
    d[cols].sort_values(["hours_after_onset", "go_name", "treatment", "direction"]).to_csv(
        out_dir("I_J") / "selected_go_terms.csv", index=False)


def main():
    prepare_freezing()
    deg = pd.read_csv(DESEQ)
    assert not deg.duplicated(["gene_id", "treatment", "time_h"]).any()
    prepare_alluvial(deg)
    prepare_lines(deg)
    prepare_scatter(deg)
    prepare_go()
    prepare_selected_go()
    print("All panel data written and checked against the source tables.")


if __name__ == "__main__":
    main()
