"""Shared publication style and encodings for all manuscript panels."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "legend.frameon": False, "axes.linewidth": .6,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": "#23272c", "text.color": "#23272c",
    "pdf.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white",
    "savefig.facecolor": "white",
})
WIDTH = 7.15
DPI = 400
METHODS = dict(zip(
    "mofa gfa icluster fabia gpcca rgcca sgcca block_pls o2pls cia mcia jive gipca bidifac bsf slide intnmf cmf tica nemo".split(),
    "MOFA|GFA|iCluster|FABIA|GPCCA|RGCCA|sGCCA|block PLS|O2PLS|CIA|MCIA|JIVE|GIPCA|BIDIFAC|BSF|SLIDE|intNMF|CMF|tICA|NEMO".split("|"),
))
STRATEGY = {"CC_SF": ("#0072B2", "s"), "MI_SF": ("#B64A00", "^"), "NA": ("#00805D", "o")}
DIAGNOSIS = {"CN": ("#4B479F", "o"), "MCI": ("#AC7100", "^"), "AD": ("#B03978", "s"),
             "GBM": ("#B03978", "s"), "LGG": ("#4B479F", "o")}
MOD_NAMES = {"gene": "Gene expression", "meta": "Metabolomics", "mrna": "mRNA",
             "methylation": "Methylation", "rppa": "RPPA"}
MODS = {"ADNI": ["gene", "meta"], "TCGA": ["mrna", "methylation", "rppa"]}
MOD_STYLE = {"gene": ("#0072B2", "o"), "meta": ("#B64A00", "^"),
             "mrna": ("#0072B2", "o"), "methylation": ("#B64A00", "^"), "rppa": ("#00805D", "s")}
AVAIL = {"CC": ("#585858", "o"), **{"missing_"+k: v for k, v in MOD_STYLE.items()},
         "missing_mrna+methylation": ("#BA4B91", "D"),
         "missing_methylation+rppa": ("#997000", "v"),
         "missing_mrna+rppa": ("#604297", "P")}
HATCH = {k: h for k, h in zip(AVAIL, ["", "/", "\\", "/", "\\", ".", "x", "-", "+"])}
MECHANISMS = ["MCAR", "MAR", "MNAR"]
METRICS = {"rv": "Latent-score RV", "recon_cor": "Reconstruction correlation",
           "wgene_rankcor": "Loading-rank Spearman correlation", "wgene_jaccard": "Top-50 Jaccard overlap"}

def availability_label(key):
    if key == "CC":
        return "Complete"
    return "Missing " + " + ".join(MOD_NAMES[m] for m in key.removeprefix("missing_").split("+"))

def marker_handles(mapping, keys, label=str):
    return [Line2D([], [], marker=mapping[k][1], color=mapping[k][0], linestyle="none",
                   markersize=4.5, label=label(k)) for k in keys]

def panel(ax, letter):
    ax.text(0, 1.025, letter, transform=ax.transAxes, ha="left", va="bottom",
            fontweight="bold", fontsize=10)

def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=DPI)
    plt.close(fig)
