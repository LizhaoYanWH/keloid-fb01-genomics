"""Figure 3: existing full ablations and fixed normal-repair readouts.

Editorial rendering only. No statistical test, score fitting, member selection,
clustering, embedding, or new biological replication is performed here.
Expected ctx contract: .data (source-key -> DataFrame), .panel(ax, label,
title), and .save(fig, stem, caption, source_keys).
"""
from __future__ import annotations

import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


SOURCE_KEYS = [
    "ablation", "time", "time_members", "time_deltas", "marker_qualification"
]
PROGRAM_25 = "fixed_shared_nonGS_FB01_25"
PROGRAM_1000 = "original_MAGMA1000_Z_RNA"
VARIANTS = [
    "remove_ECM",
    "remove_largest_RNA_blind_chromosome",
    "remove_all22_published_window_union",
]
DAYS = [0, 1, 7, 30]
CONTRASTS = [
    "D1_minus_D0", "D7_minus_D0", "D30_minus_D0", "D30_minus_D7"
]
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
GRAY = "#666666"
DONOR_COLORS = [BLUE, ORANGE, GREEN]
DONOR_MARKERS = ["o", "s", "^"]


def _symmetric_limit(values, step):
    """A disclosed display limit enclosing all finite source values."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    assert finite.size, "No finite data for the display"
    maximum = float(np.max(np.abs(finite)))
    return max(step, math.ceil(maximum / step) * step)


def render_fig3(ctx):
    data = ctx.data
    missing = set(SOURCE_KEYS).difference(data)
    assert not missing, f"Missing frozen sources: {sorted(missing)}"
    ablation = data["ablation"].copy()
    time = data["time"].copy()
    members = data["time_members"].copy()
    deltas = data["time_deltas"].copy()
    qualifications = data["marker_qualification"].copy()
    assert len(ablation) == 45 and set(ablation["variant"]) == set(VARIANTS)
    assert len(time) == 24 and len(deltas) == 24
    assert set(time["program"]) == {PROGRAM_25, PROGRAM_1000}
    assert set(deltas["program"]) == {PROGRAM_25, PROGRAM_1000}
    assert not time.duplicated(["program", "patient_id", "day"]).any()
    assert not deltas.duplicated(["program", "patient_id", "contrast"]).any()
    time["day"] = time["day"].astype(int)
    members["day"] = members["day"].astype(int)
    donors = sorted(time["patient_id"].unique())
    assert len(donors) == 3
    for program in [PROGRAM_25, PROGRAM_1000]:
        selected = time.loc[time["program"].eq(program)]
        assert len(selected) == 12
        for donor in donors:
            assert set(selected.loc[selected["patient_id"].eq(donor), "day"]) == set(DAYS)
    coverage_25 = set(time.loc[time["program"].eq(PROGRAM_25), "assayed_members"].astype(int))
    coverage_1000 = set(time.loc[time["program"].eq(PROGRAM_1000), "assayed_members"].astype(int))
    assert coverage_25 == {25} and coverage_1000 == {983}
    assert set(time["HGNC_background_genes"].astype(int)) == {26158}

    fig = plt.figure(figsize=(170 / 25.4, 185 / 25.4))
    fig.subplots_adjust(left=0.105, right=0.910, top=0.955, bottom=0.075)
    grid = fig.add_gridspec(
        3, 2, height_ratios=[0.90, 0.82, 1.20], hspace=0.63, wspace=0.58
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    top_images = []

    # A/B: all 45 state-by-exclusion entries. The existing BH45 is read, not refit.
    for ax, letter, cohort, state_count in [
        (ax_a, "a", "GSE220300", 9), (ax_b, "b", "GSE335482", 6)
    ]:
        selected = ablation.loc[ablation["cohort"].eq(cohort)]
        tab = selected.pivot(
            index="state", columns="variant", values="combined_B30_C15_BH45_q"
        ).sort_index().reindex(columns=VARIANTS)
        assert tab.shape == (state_count, 3)
        values = tab.to_numpy(dtype=float)
        assert np.all((values[np.isfinite(values)] > 0) & (values[np.isfinite(values)] <= 1))
        transformed = -np.log10(values)
        top_images.append(ax.imshow(
            transformed, aspect="auto", interpolation="nearest", vmin=0, vmax=2.6,
            cmap=plt.get_cmap("Blues").with_extremes(bad="#BDBDBD")
        ))
        ax.set(
            yticks=np.arange(state_count), yticklabels=list(tab.index),
            xticks=np.arange(3), xticklabels=["−ECM20", "−chr1 (92)", "−22 windows (61)"]
        )
        ax.tick_params(axis="x", labelrotation=22, labelsize=6)
        ax.tick_params(axis="y", labelsize=6.3, length=0)
        for iy in range(state_count):
            for ix in range(3):
                q = values[iy, ix]
                if np.isfinite(q) and q < 0.05:
                    ax.text(ix, iy, "•", ha="center", va="center", fontsize=8,
                            color="white" if q < 0.01 else "black")
                elif not np.isfinite(q):
                    ax.text(ix, iy, "NA", ha="center", va="center", fontsize=5.5)
        ctx.panel(ax, letter, cohort + ": full ablation family")
    colorbar_top = fig.colorbar(top_images[0], ax=[ax_a, ax_b], fraction=0.025, pad=0.025)
    colorbar_top.set_label("−log₁₀(BH45 q)", fontsize=6.2)
    colorbar_top.ax.tick_params(labelsize=6)
    colorbar_top.set_ticks([0, 1, 2])

    # C/D: all three donor profiles at all four categorical visits; no smoothing.
    for col, letter, program, title, coverage in [
        (0, "c", PROGRAM_25, "Normal repair: shared 25", "25/25 measured"),
        (1, "d", PROGRAM_1000, "Normal repair: genetic-set RNA", "983/1,000 measured"),
    ]:
        ax = fig.add_subplot(grid[1, col])
        for j, donor in enumerate(donors):
            selected = time.loc[
                time["program"].eq(program) & time["patient_id"].eq(donor)
            ].set_index("day").reindex(DAYS)
            ax.plot(
                np.arange(4), selected["score"].to_numpy(dtype=float),
                color=DONOR_COLORS[j], marker=DONOR_MARKERS[j], markersize=3.3,
                linewidth=1, label=donor.split("::")[-1]
            )
        ax.set(
            xticks=np.arange(4), xticklabels=["D0", "D1", "D7", "D30"],
            xlabel="Visit (categorical)", ylabel="Within-assay mean rank"
        )
        ax.tick_params(labelsize=6.2)
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.4)
        ax.margins(y=0.12)
        ax.legend(frameon=False, fontsize=5.7, handlelength=1.5, loc="upper left")
        ax.text(0.98, 0.03, coverage, ha="right", va="bottom", transform=ax.transAxes,
                color=GRAY, fontsize=5.8)
        ctx.panel(ax, letter, title)

    # E: complete fixed 25 x 12 profiles. The only transform is an explicitly
    # labelled same-gene, same-donor D0 subtraction of already computed ranks.
    selected = members.loc[members["program"].eq(PROGRAM_25)].copy()
    assert len(selected) == 25 * 12
    assert not selected.duplicated(["ensembl_gene_id", "profile_id"]).any()
    genes = sorted(selected["ensembl_gene_id"].unique())
    assert len(genes) == 25
    symbol_rows = qualifications.loc[
        qualifications["Ensembl"].isin(genes), ["Ensembl", "HGNC_symbol"]
    ].drop_duplicates()
    assert set(symbol_rows["Ensembl"]) == set(genes)
    assert not symbol_rows["Ensembl"].duplicated().any(), "Ambiguous frozen gene label"
    assert not symbol_rows["HGNC_symbol"].duplicated().any(), "Nonunique frozen symbol"
    symbols = dict(zip(symbol_rows["Ensembl"], symbol_rows["HGNC_symbol"]))
    profiles = []
    column_labels = []
    for donor in donors:
        for day in DAYS:
            profile = selected.loc[
                selected["patient_id"].eq(donor) & selected["day"].eq(day), "profile_id"
            ].unique()
            assert len(profile) == 1
            profiles.append(profile[0])
            column_labels.append("D" + str(day))
    rank_tab = selected.pivot(
        index="ensembl_gene_id", columns="profile_id", values="rank_fraction"
    ).reindex(index=genes, columns=profiles)
    raw_rank = rank_tab.to_numpy(dtype=float)
    assert raw_rank.shape == (25, 12) and np.isfinite(raw_rank).all()
    assert ((raw_rank >= 0) & (raw_rank <= 1)).all()
    difference = raw_rank - np.repeat(raw_rank[:, [0, 4, 8]], 4, axis=1)
    gene_limit = _symmetric_limit(difference, 0.05)
    ax_e = fig.add_subplot(grid[2, 0])
    im_e = ax_e.imshow(
        difference, interpolation="nearest", aspect="auto",
        cmap=plt.get_cmap("RdBu_r").with_extremes(bad="#BDBDBD"),
        norm=TwoSlopeNorm(vmin=-gene_limit, vcenter=0, vmax=gene_limit)
    )
    ax_e.set(
        yticks=np.arange(25), yticklabels=[symbols[g] for g in genes],
        xticks=np.arange(12), xticklabels=column_labels
    )
    ax_e.tick_params(axis="y", labelsize=5.6, length=0)
    ax_e.tick_params(axis="x", labelrotation=90, labelsize=5.6, length=0, pad=1.5)
    for boundary in [3.5, 7.5]:
        ax_e.axvline(boundary, color="white", linewidth=0.7)
    for j, donor in enumerate(donors):
        ax_e.text((4 * j + 2) / 12, -0.16, donor.split("::")[-1],
                  transform=ax_e.transAxes, ha="center", va="top", fontsize=5.8)
    cbar_e = fig.colorbar(im_e, ax=ax_e, fraction=0.035, pad=0.025)
    cbar_e.set_label("Gene Δrank vs donor D0", fontsize=5.8)
    cbar_e.ax.tick_params(labelsize=5.7)
    ctx.panel(ax_e, "e", "All 25 fixed genes across 12 profiles")

    # F: every existing donor/program contrast. No P, interval or FDR is derived.
    row_keys = [(p, donor) for p in [PROGRAM_25, PROGRAM_1000] for donor in donors]
    tab = deltas.pivot(index=["program", "patient_id"], columns="contrast", values="delta")
    tab = tab.reindex(index=row_keys, columns=CONTRASTS)
    vals = tab.to_numpy(dtype=float)
    assert vals.shape == (6, 4) and np.isfinite(vals).all()
    delta_limit = _symmetric_limit(vals, 0.01)
    ax_f = fig.add_subplot(grid[2, 1])
    im_f = ax_f.imshow(
        vals, interpolation="nearest", aspect="auto",
        cmap=plt.get_cmap("RdBu_r").with_extremes(bad="#BDBDBD"),
        norm=TwoSlopeNorm(vmin=-delta_limit, vcenter=0, vmax=delta_limit)
    )
    ax_f.set(
        yticks=np.arange(6),
        yticklabels=[("25" if p == PROGRAM_25 else "1,000") + " | " + d.split("::")[-1].replace("PWH", "") for p, d in row_keys],
        xticks=np.arange(4), xticklabels=["D1−D0", "D7−D0", "D30−D0", "D30−D7"],
        ylabel="Program | donor"
    )
    ax_f.tick_params(axis="y", labelsize=6, length=0)
    ax_f.tick_params(axis="x", labelrotation=35, labelsize=5.8, length=0)
    ax_f.axhline(2.5, color="white", linewidth=0.9)
    for iy in range(6):
        for ix in range(4):
            value = vals[iy, ix]
            ax_f.text(ix, iy, f"{value:+.4f}", ha="center", va="center", fontsize=5.6,
                      color="white" if abs(value) > 0.65 * delta_limit else "black")
    cbar_f = fig.colorbar(im_f, ax=ax_f, fraction=0.035, pad=0.025)
    cbar_f.set_label("Program rank change", fontsize=5.8)
    cbar_f.ax.tick_params(labelsize=5.7)
    ctx.panel(ax_f, "f", "All 24 descriptive donor contrasts")

    caption = (
        "Figure 3. Complete sensitivity analyses and shared physiological repair. "
        "(a,b) All forty-five state-by-exclusion associations retain the existing joint BH45 family. "
        "Exclusions separately remove twenty direct ECM members, ninety-two chromosome-one members, "
        "or sixty-one members within twenty-two published windows. Dots mark q<0.05; these are "
        "matched-gene-set associations, not patient tests or retained fractions of genetic effect. "
        "(c,d) The fixed shared twenty-five-gene program and original MAGMA1000-Z RNA proxy across "
        "all four categorical visits in three normally healing donors. Coverage is 25/25 and "
        "983/1,000 against the original 26,158-gene HGNC background. Lines connect donor profiles, "
        "not tracked individual cells. (e) All fixed twenty-five members across the twelve profiles, "
        "in fixed gene order and donor/visit order. Color is existing rank_fraction minus the same "
        "gene's D0 rank in the same donor, an explicitly visual subtraction, not a new differential-"
        "expression analysis; symbols come from the frozen marker qualification table. (f) All "
        "twenty-four original donor/program contrast values; labels round only for display. "
        "D0 is reused, not a new donor cohort. No new P, FDR, genotype change, pathological persistence "
        "or state-abundance inference is supplied by these expression descriptions."
    )
    assert len(caption.split()) <= 300
    ctx.save(fig, "Fig3_Complete_ablations_and_shared_repair", caption, SOURCE_KEYS)
    return {"panels": 6, "source_keys": SOURCE_KEYS, "scientific_recomputation": False}
