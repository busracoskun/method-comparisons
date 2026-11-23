"""Figures 4, 7, S1, S2, and S3: patient-level UMAP embeddings.

These need the authorized private coordinate archive (see docs/reproducibility.md); the
coordinates are only read and plotted, never exported.
"""

import numpy as np
import pandas as pd
from aggregate import read
from figure_style import (
    plt,
    WIDTH,
    METHODS,
    DIAGNOSIS,
    AVAIL,
    MOD_NAMES,
    availability_label,
    marker_handles,
    save,
)


def draw(ax, data, scope, field):
    # Scatter each label group with its color and marker
    mapping = DIAGNOSIS if field == "diagnosis" else AVAIL
    assert data[field].notna().all()
    assert set(data[field]) <= set(mapping), (field, set(data[field]) - set(mapping))
    for key, (color, marker) in mapping.items():
        q = data[data[field].eq(key)]
        if len(q):
            ax.scatter(
                q.umap1, q.umap2, color=color, marker=marker, s=7.5, alpha=0.7, linewidths=0,
                rasterized=True
            )

    # Shared limits from the full scope, equal aspect, no ticks
    xy = scope[["umap1", "umap2"]].to_numpy()
    assert np.isfinite(xy).all() and len(data)
    low, high = xy.min(axis=0), xy.max(axis=0)
    low, high = low - 0.055 * (high - low), high + 0.055 * (high - low)
    ax.set(xlim=(low[0], high[0]), ylim=(low[1], high[1]), xticks=[], yticks=[])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)


def legend(ax, mapping, data, field, title):
    ax.axis("off")
    keys = [x for x in mapping if x in set(data[field])]
    ax.legend(
        handles=marker_handles(
            mapping, keys, availability_label if field == "availability" else str
        ),
        title=title, title_fontsize=7.5, loc="center left", borderaxespad=0, handletextpad=0.25,
        labelspacing=0.56
    )


def mark_panels(fig, axes):
    # Equal-aspect axes shrink by different amounts, so panel letters are anchored to the
    # subplot slots to keep them clear of the centered titles
    for i, ax in enumerate(axes):
        slot = ax.get_position(original=True)
        fig.text(
            slot.x0, slot.y1 + 0.017, chr(65 + i), fontsize=10, fontweight="bold", ha="left",
            va="bottom"
        )


def observed(path, project, supplemental=False):
    """Figure 4 (or S3 with supplemental=True): observed-data UMAPs."""

    # Method/strategy combinations and their coordinate source
    archive = project / "results/historical_manuscript_figures/tables"
    datasets = ["ADNI", "TCGA"]
    combos = (
        [("gpcca", "NA"), ("mcia", "CC_SF"), ("gfa", "MI_SF")]
        if supplemental
        else [("mofa", "NA"), ("tica", "MI_SF")]
    )
    height = len(combos)
    if supplemental:
        coordinates = read(project / "supplementary/tables/S3_observed_UMAP.csv")
    else:
        manifest = read(archive / "observed_UMAP_manifest.csv")

    # One row per combination, diagnosis and availability side by side
    fig = plt.figure(figsize=(WIDTH, 9.1 if supplemental else 6.5))
    grid = fig.add_gridspec(2 * height, 3, width_ratios=[1, 1, 0.89], hspace=0.52, wspace=0.10)
    plots = []
    for i, ds in enumerate(datasets):
        groups = []

        for j, (model, strategy) in enumerate(combos):
            if supplemental:
                q = coordinates[
                    (coordinates.dataset == ds.lower())
                    & (coordinates.model == model)
                    & (coordinates.strategy == strategy)
                ]
                expected = {"ADNI": (889, 450), "TCGA": (841, 465)}[ds][strategy == "CC_SF"]
                assert len(q) == expected
            else:
                q = read(archive / f"UMAP_{ds}_{model}_{strategy}.csv").rename(
                    columns={
                        "UMAP1": "umap1", "UMAP2": "umap2", "missingness_pattern": "availability"
                    }
                )
                row = manifest[
                    (manifest.dataset == ds)
                    & (manifest.model == model)
                    & (manifest["style"] == strategy)
                ].iloc[0]
                assert len(q) == row.n

            assert len(q) and q.patient_id.is_unique
            groups.append(q)

            for k, field in enumerate(["diagnosis", "availability"]):
                if k == 1 and strategy == "CC_SF":
                    continue
                ax = fig.add_subplot(grid[i * height + j, k])
                draw(ax, q, q, field)
                plots.append(ax)
                ax.set_title(
                    f"{ds} · {METHODS[model]} · {strategy}" if k == 0 else "Availability", pad=6
                )

        # Legends for this dataset
        all_data = pd.concat(groups)
        keygrid = grid[i * height : (i + 1) * height, 2].subgridspec(
            2, 1, height_ratios=[1, 2], hspace=0
        )
        legend(fig.add_subplot(keygrid[0]), DIAGNOSIS, all_data, "diagnosis", "Diagnosis")
        legend(
            fig.add_subplot(keygrid[1]), AVAIL, all_data, "availability", "Modality availability"
        )

    assert len(plots) == (10 if supplemental else 8)
    fig.supxlabel("UMAP 1", fontsize=8, y=0.017)
    fig.supylabel("UMAP 2", fontsize=8, x=0.014)
    fig.subplots_adjust(left=0.065, right=0.975, top=0.94, bottom=0.065)
    mark_panels(fig, plots)
    save(fig, path)


def recovery(path, project, dataset=None):
    """Figure 7 (or S1/S2 for one dataset): tICA, MCIA, and GFA recovery UMAPs."""

    # Panels to draw: MCIA/GFA for one dataset, or tICA for both
    if dataset:
        filename = f"S{1 if dataset == 'ADNI' else 2}_{dataset}_semisynthetic_UMAP.csv"
        data = read(project / "supplementary/tables" / filename)
        mods = ["gene", "meta"] if dataset == "ADNI" else ["methylation", "rppa"]
        specs = [(dataset, m, mod) for m in ["mcia", "gfa"] for mod in mods]
    else:
        archive = project / "results/historical_manuscript_figures"
        specs = [
            ("ADNI", "tica", "gene"), ("ADNI", "tica", "meta"),
            ("TCGA", "tica", "methylation"), ("TCGA", "tica", "rppa"),
        ]

    # Reference, 20%, and 70% embeddings for each panel row
    fig, axes = plt.subplots(4, 3, figsize=(WIDTH, 7.1))
    for i, (ds, model, mod) in enumerate(specs):
        if dataset:
            q = data[(data.model == model) & (data.modality == mod)]
            frames = [q[q.panel.eq(j)] for j in range(3)]
            assert all(len(z) == {"ADNI": 450, "TCGA": 465}[dataset] for z in frames)
        else:
            clinical = pd.read_csv(
                archive / "cache" / f"{ds}_clinical.csv", dtype={"patient_id": str}
            )
            clinical.diagnosis = clinical.diagnosis.str.upper()
            manifest = read(archive / "tables" / f"{ds}_{mod}_visualization_manifest.csv")
            assert set(manifest.model) == {"tica"}

            frames = []
            for condition in ["reference", "20", "70"]:
                q = pd.read_csv(
                    archive / "cache" / f"{ds}_{mod}_{condition}_UMAP.csv",
                    dtype={"patient_id": str}, float_precision="round_trip"
                )
                q = q.rename(columns={"UMAP1": "umap1", "UMAP2": "umap2"})
                q = q.merge(
                    clinical[["patient_id", "diagnosis"]], on="patient_id", validate="one_to_one"
                )
                assert len(q) == manifest.n.iloc[0]
                frames.append(q)

        assert all(
            set(q.patient_id) == set(frames[0].patient_id) and q.patient_id.is_unique
            for q in frames
        )

        scope = pd.concat(frames)
        for j, q in enumerate(frames):
            ax = axes[i, j]
            draw(ax, q, scope, "diagnosis")

    # Legend and row/column labels
    handles = []
    for ds, diagnoses in [("ADNI", ["CN", "MCI", "AD"]), ("TCGA", ["GBM", "LGG"])]:
        if dataset in (None, ds):
            handles.extend(
                marker_handles(DIAGNOSIS, diagnoses, lambda k: f"{ds} {k}" if not dataset else k)
            )
    fig.legend(handles=handles, loc="upper center", ncol=len(handles))
    fig.supxlabel("Reference UMAP 1", fontsize=8, y=0.015)
    fig.supylabel("Reference UMAP 2", fontsize=8, x=0.015)
    fig.subplots_adjust(left=0.166, right=0.98, top=0.866, bottom=0.055, hspace=0.25, wspace=0.13)
    for j, summary in enumerate(["Complete reference", "20%", "70%"]):
        slot = axes[0, j].get_position(original=True)
        fig.text((slot.x0 + slot.x1) / 2, 0.917, summary, ha="center", va="center", fontsize=9)
    for i, (ds, model, mod) in enumerate(specs):
        slot = axes[i, 0].get_position(original=True)
        fig.text(
            0.148, (slot.y0 + slot.y1) / 2,
            f"{ds if not dataset else METHODS[model]}\n{MOD_NAMES[mod]}", fontsize=7, ha="right",
            va="center"
        )

    mark_panels(fig, axes.flat)
    save(fig, path)
