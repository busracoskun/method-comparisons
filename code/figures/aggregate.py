"""Render all aggregate figures directly from the frozen manuscript tables."""
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, FuncFormatter
from figure_style import (
    plt, WIDTH, METHODS, STRATEGY, AVAIL, HATCH, MOD_STYLE, MOD_NAMES, MODS,
    MECHANISMS, METRICS, availability_label, marker_handles, panel, save,
)

DATASETS = list(MODS)
RATES = [.2, .7]
STATUS = {"unsupported": "U", "not attempted": "–", "attempted but failed": "F"}
MECH_STYLE = dict(zip(MECHANISMS, [("#0072B2", "o"), ("#B64A00", "^"), ("#00805D", "s")]))

def read(path):
    return pd.read_csv(path, keep_default_na=False, na_values=["", "NaN", "nan"],
                       float_precision="round_trip")

def method_axis(ax, methods, labels=True):
    ax.set_yticks(range(len(methods)), [METHODS[m] if labels else "" for m in methods])
    ax.set_ylim(len(methods)-.5, -.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)

def fitted(ax, x, y, lo, hi, style):
    color, marker = style
    if np.isfinite([lo, hi]).all():
        ax.hlines(y, lo, hi, lw=.85, color=color)
    ax.scatter(x, y, s=12, color=color, marker=marker, linewidths=.15, zorder=3)

def missingness(d, path, *, hatch=True):
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.75), sharey=True)
    for ax, ds, groups, letter in zip(axes, DATASETS, [["CN", "MCI", "AD"], ["GBM", "LGG"]], "AB"):
        q = d[d.dataset.eq(ds)]
        totals = q.groupby("diagnosis").total.first().reindex(groups)
        assert q.groupby("diagnosis").n.sum().reindex(groups).equals(totals)
        bottom = np.zeros(len(groups))
        for k, (color, _) in AVAIL.items():
            g = q[q.pattern.eq(k)]
            if g.empty:
                continue
            heights = g.set_index("diagnosis").percent.reindex(groups, fill_value=0).to_numpy()
            ax.bar(groups, heights, bottom=bottom, color=color, edgecolor="white", linewidth=.5,
                   width=.68, hatch=HATCH[k] if hatch else None, label=availability_label(k))
            for i, value in enumerate(heights):
                if value >= 9:
                    ax.text(i, bottom[i]+value/2, f"{value:.1f}%", color="white", ha="center",
                            va="center", bbox={"facecolor": color, "edgecolor": "none", "pad": .7})
            bottom += heights
        np.testing.assert_allclose(bottom, 100)
        ax.set(ylim=(0, 100), title=ds)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_xticks(range(len(groups)), [f"{g}\n(n = {totals[g]})" for g in groups])
        ax.legend(loc="upper left", bbox_to_anchor=(-.01, -.22), borderaxespad=0,
                  ncol=1 if ds == "ADNI" else 2, columnspacing=.8, handlelength=1.5)
        panel(ax, letter)
    axes[0].set_ylabel("Patients (%)")
    fig.subplots_adjust(left=.08, right=.98, bottom=.34, top=.90, wspace=.14)
    save(fig, path)

def clinical(d, reference, path, by_strategy=False, reliability=None):
    ari = reliability is not None
    first, last = ("ari_diag", "ari_miss") if ari else ("rand_diag", "rand_miss")
    order = ([("ADNI", first), ("TCGA", first), ("ADNI", last), ("TCGA", last)] if ari else
             [(ds, metric) for ds in DATASETS for metric in [first, "R2", last]])
    methods = list(METHODS)
    assert set(d.model) <= set(methods) and set(d["style"]) <= set(STRATEGY)
    assert not d.duplicated(["dataset", "metric", "model", "style"]).any()
    if not ari:
        assert not ((d.dataset == "TCGA") & (d.metric == first) & (d["style"] == "CC_SF")).any()
        assert d.loc[(d.metric == last) & (d["style"] == "CC_SF"), "mean"].isna().all()
    lo, hi = ("p025", "p975") if ari else ("low", "high")
    bounds = {}
    for metric, group in d.groupby("metric"):
        v = group[[lo, hi, "mean"]].to_numpy()
        left, right = np.nanmin(v), np.nanmax(v)
        if metric == last and not ari:
            left, right = min(left, min(reference)), max(right, max(reference))
        margin = .06 * (right-left)
        bounds[metric] = (min(0, left-margin) if ari else max(0, left-margin),
                          max(0, right+margin) if ari else min(1, right+margin))
    if by_strategy:
        return clinical_facets(d, reference, path, order, methods, bounds)
    cols = 2 if ari else 3
    fig, axes = plt.subplots(2, cols, figsize=(WIDTH, 8.0))
    plotted = 0
    for idx, (ax, (ds, metric)) in enumerate(zip(axes.flat, order)):
        q = d[(d.dataset == ds) & (d.metric == metric)]
        text = ("MMSE R²" if ds == "ADNI" else "Karnofsky R²") if metric == "R2" else (
            "Availability" if metric == last else "Diagnosis" if ds == "ADNI" else "Disease group")
        ax.set_title(f"{ds} · {text}", pad=8)
        panel(ax, chr(65+idx))
        strategies = list(STRATEGY)
        blocked = metric == last or (not ari and ds == "TCGA" and metric == first)
        ax.set_xlim(bounds[metric])
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.set_xlabel("Adjusted Rand index" if ari else "R²" if metric == "R2" else "Rand index")
        method_axis(ax, methods, idx % cols == 0)
        if metric == last and not ari:
            ax.axvline(reference.loc[ds], color=".35", lw=.75, ls="--", zorder=0)
        if ari:
            ax.axvline(0, color=".8", lw=.65, zorder=0)
        for i, model in enumerate(methods):
            if i % 2 == 0:
                ax.axhspan(i-.5, i+.5, color="#f3f5f6", lw=0, zorder=-1)
            for j, strategy in enumerate(strategies):
                if blocked and strategy == "CC_SF":
                    continue
                g = q[(q.model == model) & (q["style"] == strategy)]
                y = i + (j-1)*.25
                if len(g) and pd.notna(g.iloc[0]["mean"]):
                    r = g.iloc[0]
                    fitted(ax, r["mean"], y, r[lo], r[hi], STRATEGY[strategy])
                    plotted += 1
                elif ari:
                    r = reliability.loc[(ds, model, strategy)]
                    symbol = STATUS.get(r.state, "–" if not r.attempted_fits else "M")
                    ax.text(1.02+j*.055, i, symbol, transform=ax.get_yaxis_transform(),
                            color=STRATEGY[strategy][0], va="center", ha="left", fontsize=6)
        if ari:
            for j, strategy in enumerate(strategies):
                if not (blocked and strategy == "CC_SF"):
                    ax.text(1.025+j*.055, 1.006, strategy.split("_")[0], fontsize=5,
                            transform=ax.transAxes, color=STRATEGY[strategy][0])
    expected = d["mean"].notna()
    if ari:
        expected &= ~((d.metric == last) & (d["style"] == "CC_SF"))
    assert plotted == expected.sum(), (path, plotted, expected.sum())
    handles = marker_handles(STRATEGY, STRATEGY)
    if not ari:
        handles.append(Line2D([], [], color=".35", ls="--", lw=.75, label="Diagnosis–availability Rand index"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .995), ncol=len(handles))
    fig.subplots_adjust(left=.105, right=.85 if ari else .98, bottom=.054, top=.924,
                        wspace=.47 if ari else .24, hspace=.28)
    save(fig, path)

def clinical_facets(d, reference, path, order, methods, bounds):
    """Alternate Figure 2 with one column per strategy and one row per outcome."""
    fig, axes = plt.subplots(6, 3, figsize=(WIDTH, 15.1))
    plotted = 0
    for i, (ds, metric) in enumerate(order):
        text = {"rand_diag": "Diagnosis" if ds == "ADNI" else "Disease group",
                "R2": "MMSE R²" if ds == "ADNI" else "Karnofsky R²",
                "rand_miss": "Modality availability"}[metric]
        for j, (strategy, style) in enumerate(STRATEGY.items()):
            ax = axes[i,j]
            q = d[d.dataset.eq(ds) & d.metric.eq(metric) & d["style"].eq(strategy)].set_index("model")
            ax.set_xlim(bounds[metric])
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.set_xlabel("R²" if metric == "R2" else "Rand index", labelpad=1)
            method_axis(ax, methods, j == 0)
            ax.tick_params(axis="y", labelsize=6)
            if j == 0:
                panel(ax, chr(65+i))
            if j == 1:
                ax.set_title(f"{ds} · {text}", pad=6)
            if metric == "rand_miss" and strategy != "CC_SF":
                ax.axvline(reference.loc[ds], color=".35", lw=.7, ls="--")
            for y, m in enumerate(methods):
                if y % 2 == 0:
                    ax.axhspan(y-.5, y+.5, color="#f3f5f6", lw=0, zorder=-1)
                if m in q.index and pd.notna(q.loc[m,"mean"]):
                    r = q.loc[m]
                    fitted(ax, r["mean"], y, r.low, r.high, style)
                    plotted += 1
    assert plotted == d["mean"].notna().sum()
    handles = marker_handles(STRATEGY, STRATEGY)
    handles.append(Line2D([], [], color=".35", ls="--", lw=.75, label="Diagnosis–availability Rand index"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5,.998), ncol=4)
    fig.subplots_adjust(left=.105,right=.98,top=.956,bottom=.026,hspace=.45,wspace=.16)
    save(fig,path)

def metric_scale(metric):
    if metric == "wgene_rankcor":
        return plt.get_cmap("RdBu"), TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    return plt.get_cmap("cividis"), Normalize(0, 1)

def heat(ax, matrix, metric, digits=2, colorbar=None):
    cmap, norm = metric_scale(metric)
    cmap = cmap.with_extremes(bad="#e8e9eb")
    values = np.asarray(matrix, dtype=float)
    finite = values[np.isfinite(values)]
    if len(finite):
        assert finite.min() >= norm.vmin and finite.max() <= norm.vmax, (metric, finite.min(), finite.max())
    mesh = ax.pcolormesh(np.arange(values.shape[1]+1)-.5, np.arange(values.shape[0]+1)-.5,
                         values, cmap=cmap, norm=norm, edgecolors="white", linewidth=.25)
    for (i, j), val in np.ndenumerate(values):
        rgba = cmap(norm(val)) if np.isfinite(val) else cmap.get_bad()
        light = sum(c*w for c, w in zip(rgba, [.2126, .7152, .0722]))
        ax.text(j, i, f"{val:.{digits}f}" if np.isfinite(val) else "–",
                ha="center", va="center", color="#161616" if light > .5 else "white", fontsize=6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if colorbar is not None:
        cb = ax.figure.colorbar(mesh, cax=colorbar, ticks=[norm.vmin, (norm.vmin+norm.vmax)/2, norm.vmax])
        cb.outline.set_visible(False)
    return mesh

def recovery(d, dataset, path, alternative=False):
    d = d[d.dataset.eq(dataset) & d.metric.isin(["rv", "recon_cor"]) & d.n.gt(0)]
    methods = [m for m in METHODS if m in set(d.model)]
    mods = MODS[dataset]
    assert set(d.model) == set(methods)
    assert set(d.dropped_modality) <= set(mods)
    assert set(d.mechanism) <= set(MECHANISMS)
    assert set(d.target_rate) <= set(RATES)
    fig = plt.figure(figsize=(WIDTH, 7.1 if alternative else 6.8))
    cols = len(mods) if alternative else 2
    spec = [1]*cols + ([] if alternative else [.046])
    grid = fig.add_gridspec(2, len(spec), width_ratios=spec, hspace=.38, wspace=.17 if alternative else .12)
    plotted = 0
    for i, metric in enumerate(["rv", "recon_cor"]):
        keys = pd.MultiIndex.from_product([mods, MECHANISMS])
        for j in range(cols):
            ax = fig.add_subplot(grid[i, j])
            q = d[d.metric.eq(metric)]
            if alternative:
                mod = mods[j]
                q = q[q.dropped_modality.eq(mod)]
                for y, m in enumerate(methods):
                    for k, mechanism in enumerate(MECHANISMS):
                        pair = q[(q.model == m) & (q.mechanism == mechanism)].sort_values("target_rate")
                        offset = (k-1)*.23
                        xx = pair["mean"].to_numpy()
                        yy = y + offset + np.where(pair.target_rate.eq(.2), -.065, .065)
                        color, marker = MECH_STYLE[mechanism]
                        if len(pair) == 2:
                            ax.plot(xx, yy, color=color, lw=.65, alpha=.75)
                        for x, yp, rate in zip(xx, yy, pair.target_rate):
                            ax.scatter(x, yp, marker=marker, s=10, linewidths=.7, edgecolors=color,
                                       facecolors=color if rate == .2 else "white", zorder=3)
                            plotted += 1
                ax.set(xlim=(-.02, 1.02), title=MOD_NAMES[mod], xlabel="Mean recovery")
                ax.set_xticks([0, .5, 1])
                if j == 0:
                    panel(ax, chr(65+i))
            else:
                q = q[q.target_rate.eq(RATES[j])]
                table = q.pivot(index="model", columns=["dropped_modality", "mechanism"], values="mean")
                matrix = table.reindex(index=methods, columns=keys)
                plotted += matrix.notna().to_numpy().sum()
                heat(ax, matrix, metric, colorbar=fig.add_subplot(grid[i, -1]) if j == 1 else None)
                for k, mod in enumerate(mods):
                    ax.text(k*3+1, 1.025, MOD_NAMES[mod], transform=ax.get_xaxis_transform(),
                            ha="center", va="bottom", fontsize=7)
                    if k:
                        ax.axvline(k*3-.5, color="white", lw=2)
                ax.set_xticks(range(len(keys)), MECHANISMS*len(mods))
                ax.tick_params(axis="x", labelsize=5.5)
                ax.set_title(f"{RATES[j]:.0%}", y=1.085)
                if j == 0:
                    panel(ax, chr(65+i))
            method_axis(ax, methods, j == 0)
        fig.text(.116, .944-i*.474, METRICS[metric], fontsize=9)
    assert plotted == d["mean"].notna().sum(), (path, plotted, len(d))
    if alternative:
        handles = marker_handles(MECH_STYLE, MECHANISMS)
        handles += [Line2D([], [], marker="o", ls="", color=".25", markerfacecolor=c, label=l)
                    for c,l in [(".25","20%"), ("white","70%")]]
        fig.legend(handles=handles, ncol=5, loc="lower center")
    fig.subplots_adjust(left=.116, right=.936, bottom=.085 if alternative else .065, top=.867)
    save(fig, path)

def feature(d, path):
    fig, axes = plt.subplots(4, 2, figsize=(WIDTH, 9.0))
    for row, (ds, metric) in enumerate((ds, m) for ds in DATASETS for m in ["wgene_rankcor", "wgene_jaccard"]):
        q = d[d.dataset.eq(ds.lower()) & d.metric.eq(metric)]
        methods = [m for m in METHODS if m in set(q.model)]
        for col, rate in enumerate(RATES):
            ax = axes[row,col]
            matrix = q[q.rate.eq(rate)].pivot(index="model", columns="modality", values="value")
            matrix = matrix.reindex(index=methods, columns=MODS[ds])
            heat(ax, matrix, metric)
            method_axis(ax, methods, col == 0)
            ax.set_xticks(range(len(MODS[ds])), [MOD_NAMES[m] for m in MODS[ds]])
            ax.set_title(f"{rate:.0%}")
            panel(ax, chr(65+row*2+col))
        axes[row,0].text(0, 1.13, f"{ds} · {METRICS[metric]}", transform=axes[row,0].transAxes, fontsize=8.5)
    fig.subplots_adjust(left=.11, right=.98, top=.94, bottom=.037, hspace=.43, wspace=.09)
    save(fig, path)

def reliability(d, path):
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 6.0))
    for ax, ds, letter in zip(axes, DATASETS, "AB"):
        q = d[d.dataset.eq(ds.lower())].set_index(["model", "style"])
        for i, model in enumerate(METHODS):
            for j, strategy in enumerate(STRATEGY):
                r = q.loc[(model, strategy)]
                if r.attempted_fits:
                    np.testing.assert_allclose(r.success_rate, r.successful_fits/r.attempted_fits)
                    text = f"{int(r.successful_fits)}/{int(r.attempted_fits)}"
                    color = "#f7f8fa" if r.success_rate == 1 else "#f7e5ae" if r.successful_fits else "#eab2a2"
                    if 0 < r.success_rate < 1:
                        ax.barh(i+.31, .84*r.success_rate, left=j-.42, height=.05, color="#69571e")
                else:
                    text = STATUS.get(r.state, "–")
                    color = "#eceff0"
                ax.add_patch(plt.Rectangle((j-.48,i-.46), .96,.92, facecolor=color, lw=0, zorder=-1))
                ax.text(j, i-.03, text, va="center", ha="center", fontsize=6.8)
        method_axis(ax, list(METHODS), ds == "ADNI")
        ax.set(xlim=(-.5, 2.5), title=ds)
        ax.set_xticks(range(3), list(STRATEGY))
        ax.xaxis.tick_top()
        ax.tick_params(length=0)
        ax.spines["bottom"].set_visible(False)
        panel(ax, letter)
    fig.legend(handles=[plt.Rectangle((0,0),1,1, color=c, label=l) for c,l in [
        ("#f7e5ae","Partial success"), ("#eab2a2","All failed"), ("#eceff0","U Unsupported · – No record")
    ]], loc="lower center", ncol=3)
    fig.subplots_adjust(left=.106, right=.98, top=.89, bottom=.07, wspace=.10)
    save(fig, path)

def patient_feature(d, path):
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 5.9), sharex=True, sharey=True)
    for i, ds in enumerate(DATASETS):
        for j, rate in enumerate(RATES):
            ax = axes[i,j]
            q = d[d.dataset.eq(ds.lower()) & d.rate.eq(rate)]
            assert set(q.modality) <= set(MODS[ds])
            for m in MODS[ds]:
                g = q[q.modality.eq(m)]
                ax.scatter(g.rv, g.jaccard, color=MOD_STYLE[m][0], marker=MOD_STYLE[m][1], s=18,
                           alpha=.78, edgecolors="white", linewidths=.25)
            ax.set(xlim=(-.02,1.02), ylim=(-.02,1.02), title=f"{ds} · {rate:.0%}")
            panel(ax, chr(65+i*2+j))
            if j == 0:
                ax.set_ylabel("Top-50 Jaccard overlap")
            if i == 1:
                ax.set_xlabel("Latent-score RV")
        axes[i,1].legend(handles=marker_handles(MOD_STYLE, MODS[ds], MOD_NAMES.get), loc="upper left")
    fig.subplots_adjust(left=.095, right=.98, bottom=.09, top=.92, hspace=.28, wspace=.13)
    save(fig, path)

def sensitivity(d, path):
    fig, axes = plt.subplots(4, 2, figsize=(WIDTH, 8.65))
    for i, m in enumerate(METRICS):
        for j, ds in enumerate(DATASETS):
            ax = axes[i,j]
            q = d[d.dataset.eq(ds) & d.metric.eq(m)]
            keys = pd.MultiIndex.from_product([MECHANISMS, RATES])
            matrix = q.pivot(index="dropped_modality", columns=["mechanism","target_rate"], values="value")
            heat(ax, matrix.reindex(index=MODS[ds], columns=keys), m, digits=3)
            ax.set_ylim(len(MODS[ds])-.5, -.5)
            ax.set_yticks(range(len(MODS[ds])), [MOD_NAMES[x] for x in MODS[ds]])
            ax.set_xticks(range(6), ["20%", "70%"]*3)
            for k, mech in enumerate(MECHANISMS):
                ax.text(k*2+.5, 1.02, mech, transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
            ax.set_title(f"{ds} · {METRICS[m]}", y=1.11, fontsize=8)
            panel(ax, chr(65+i*2+j))
    fig.text(.5,.983,"intNMF · native nonnegative preprocessing", ha="center", va="top", fontsize=9)
    fig.subplots_adjust(left=.165, right=.98, bottom=.03, top=.90, wspace=.43, hspace=.58)
    save(fig, path)

def runtime(d, path, observed=False):
    keys = list(STRATEGY) if observed else RATES
    fig, axes = plt.subplots(2, 1 if observed else 2, figsize=(WIDTH, 7.9 if observed else 6.7),
                             squeeze=False, sharex=True)
    values = d[["q25_seconds", "q75_seconds", "median_seconds"]].to_numpy()
    limits = (10**np.floor(np.log10(np.nanmin(values))), 10**np.ceil(np.log10(np.nanmax(values))))
    count = 0
    for i, ds in enumerate(DATASETS):
        group = d[d.dataset.eq(ds.lower())]
        methods = [m for m in METHODS if m in set(group.model)]
        for j, ax in enumerate(axes[i]):
            q = group if observed else group[np.isclose(group.target_rate, keys[j])]
            level = "strategy" if observed else "dropped_modality"
            series = keys if observed else MODS[ds]
            mapping = STRATEGY if observed else MOD_STYLE
            for idx, model in enumerate(methods):
                if idx % 2 == 0:
                    ax.axhspan(idx-.5, idx+.5, color="#f3f5f6", lw=0, zorder=-1)
                for k, item in enumerate(series):
                    g = q[q.model.eq(model) & q[level].eq(item)]
                    y = idx + (k-(len(series)-1)/2)*.25
                    if len(g) and pd.notna(g.iloc[0].median_seconds):
                        r = g.iloc[0]
                        fitted(ax, r.median_seconds, y, r.q25_seconds, r.q75_seconds, mapping[item])
                        count += 1
                    elif observed:
                        r = g.iloc[0]
                        symbol = {"not represented in saved workflow":"U", "no successful fits":"F",
                                  "timing missing":"M", "available":"M"}[r.coverage]
                        ax.text(1.025+k*.045, idx, symbol, transform=ax.get_yaxis_transform(),
                                ha="center", va="center", color=mapping[item][0], fontsize=6)
            if observed:
                for k, item in enumerate(series):
                    ax.text(1.025+k*.045, 1.008, item.split("_")[0], transform=ax.transAxes,
                            ha="center", color=mapping[item][0], fontsize=6)
            method_axis(ax, methods, j == 0)
            ax.set(xscale="log", xlim=limits, title=ds if observed else f"{ds} · {keys[j]:.0%}")
            tick_powers = np.arange(round(np.log10(limits[0])), round(np.log10(limits[1]))+1, 2)
            ax.set_xticks(10.**tick_powers)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
            ax.tick_params(axis="x", which="minor", length=0)
            ax.grid(axis="x", color=".86", lw=.5)
            ax.set_axisbelow(True)
            ax.tick_params(axis="x", labelbottom=True)
            panel(ax, chr(65+i*len(axes[i])+j))
            if not observed and j == 1:
                ax.legend(handles=marker_handles(mapping, series, MOD_NAMES.get), loc="lower right")
    assert count == d.median_seconds.notna().sum()
    fig.supxlabel("Elapsed seconds per successful fit" if observed else "Elapsed seconds per successful refit",
                  fontsize=8)
    if observed:
        fig.legend(handles=marker_handles(STRATEGY, STRATEGY), loc="upper center", ncol=3)
    fig.subplots_adjust(left=.11, right=.865 if observed else .955, top=.92 if observed else .945,
                        bottom=.07, wspace=.16, hspace=.25)
    save(fig, path)
