"""Render the three observed-data aggregate figures."""
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[2]
from aggregate import fitted, method_axis, read
from figure_style import plt, METHODS, STRATEGY, WIDTH, marker_handles, panel, save

DATA = ROOT / "results/observed_data"
COHORT = {"ADNI": ("#0072B2", "o"), "TCGA": ("#B64A00", "s")}
CLINICAL_METRICS = ["rand_diag", "R2", "rand_miss"]

def clinical_data():
    folder = DATA
    shown = read(folder/"clinical_alignment_plotted_summary.csv").drop(columns="panel")
    bootstrap = read(folder/"clinical_alignment_bootstrap_summary.csv")
    eligible = ~(
        bootstrap.dataset.eq("TCGA") & bootstrap.metric.eq("rand_diag") & bootstrap["style"].eq("CC_SF")
    )
    keys = ["dataset", "metric", "model", "style"]
    orig = bootstrap[eligible & bootstrap.metric.isin(CLINICAL_METRICS)]
    pd.testing.assert_frame_equal(shown.sort_values(keys).reset_index(drop=True),
                                  orig.sort_values(keys).reset_index(drop=True),
                                  check_exact=False, rtol=1e-13, atol=1e-15)
    ref = read(folder/"diagnosis_missingness_references.csv").set_index("dataset").RI
    return shown, ref

def clinical_limits(data, references):
    limits = {}
    for metric, group in data.groupby("metric"):
        values = group[["low", "mean", "high"]].to_numpy().flatten()
        if metric == "rand_miss":
            values = np.r_[values, references]
        lo, hi = np.nanmin(values), np.nanmax(values)
        margin = .045*(hi-lo)
        limits[metric] = (max(0, lo-margin), min(1, hi+margin))
    return limits

def clinical_figure(data, references, limits, dataset, path):
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 4.95))
    clinical = ["Diagnosis", "MMSE"] if dataset == "ADNI" else ["Disease group", "Karnofsky score"]
    plotted = 0
    for idx, (ax, metric, title) in enumerate(zip(axes, CLINICAL_METRICS, [*clinical, "Modality availability"])):
        q = data[data.metric.eq(metric)].set_index(["model", "style"])
        for y, model in enumerate(METHODS):
            if y % 2 == 0:
                ax.axhspan(y-.5, y+.5, color="#f1f3f4", lw=0, zorder=-1)
            for j, (strategy, style) in enumerate(STRATEGY.items()):
                key = (model, strategy)
                if key not in q.index or pd.isna(q.loc[key, "mean"]):
                    continue
                r = q.loc[key]
                fitted(ax, r["mean"], y + (j-1)*.265, r.low, r.high, style)
                plotted += 1
        if metric == "rand_miss":
            ax.axvline(references.loc[dataset.lower()], color=".3", lw=.8, ls="--", zorder=0)
        method_axis(ax, list(METHODS), idx == 0)
        ax.set(xlim=limits[metric], xlabel="R²" if metric == "R2" else "Rand index")
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.set_title(title, pad=10)
        panel(ax, chr(65+idx))
    assert plotted == data["mean"].notna().sum()
    handles = marker_handles(STRATEGY, STRATEGY)
    handles += [Line2D([], [], color=".3", lw=.8, ls="--", label="Diagnosis–availability Rand index")]
    fig.suptitle(dataset, fontsize=11, y=.985)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .943), ncol=4)
    fig.subplots_adjust(left=.11, right=.978, top=.814, bottom=.103, wspace=.18)
    save(fig, path)

