"""Draw S5 from the original semi-synthetic scalar exports."""
from io import BytesIO

import fitz
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator

from figure_style import plt

MECHANISMS = ["MCAR", "MAR", "MNAR"]
EXAMPLES = [
    ('tcga','mofa','methylation',.2,'rv','MOFA · 20% methylation'),
    ('tcga','mofa','methylation',.7,'rv','MOFA · 70% methylation'),
    ('tcga','icluster','methylation',.2,'rv','iCluster · 20% methylation'),
    ('tcga','icluster','methylation',.7,'rv','iCluster · 70% methylation'),
    ('tcga','mofa','methylation',.7,'recon_cor','MOFA · 70% methylation'),
    ('adni','cia','meta',.7,'recon_cor','CIA · 70% metabolomics'),
]

def render_mechanisms(results, path):
    records = pd.concat(
        [pd.read_csv(results / f"{dataset}_semi_synthetic_results.csv").assign(dataset=dataset)
         for dataset in ("adni", "tcga")],
        ignore_index=True,
    )
    for metric in ("rv", "recon_cor"):
        records[metric] = pd.to_numeric(records[metric], errors="coerce").replace([np.inf, -np.inf], np.nan)
    with plt.rc_context(rc=plt.rcParamsDefault):
        draw_mechanisms(records, path)


def draw_mechanisms(records, path):
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.spines.left':False,'axes.spines.bottom':False,'text.color':'#202b36',
                         'axes.labelcolor':'#202b36','xtick.color':'#202b36','ytick.color':'#455260'})
    fig,axs=plt.subplots(3,2,figsize=(8.6,10.2),sharey=True)
    rng=np.random.default_rng(19)
    for ax,(ds,model,mod,rate,metric,title) in zip(axs.flat,EXAMPLES, strict=True):
        subset=records.loc[
            records.dataset.eq(ds) & records.model.eq(model) & records.dropped_modality.eq(mod)
            & records.target_rate.eq(rate) & records.status.eq("success")
        ]
        labels=[]
        for i,(mech,color) in enumerate(zip(MECHANISMS,['#0072B2','#A45A00','#8A438D'])):
            x=subset.loc[subset.mechanism.eq(mech),metric].dropna()
            if len(x) != 100:
                raise ValueError(f"{ds} {model} {mod} {rate} {mech} has {len(x)} successful finite results")
            labels.append(f'{mech}\nn={len(x)}')
            ax.scatter(i+rng.uniform(-.21,.21,len(x)),x,s=11,c=color,alpha=.36,linewidths=0,zorder=2,rasterized=True)
            ax.errorbar(i,x.mean(),yerr=x.std(),fmt='D',color='#19232d',markerfacecolor='white',markersize=5,
                        markeredgewidth=1.15,elinewidth=1.4,capsize=4,zorder=3)
        ax.set_title(f'{ds.upper()}  /  {title}',loc='left',fontsize=10.5,fontweight='bold',pad=9)
        ax.set_xticks(range(3),labels,fontsize=9)
        ax.set(ylim=(0,1.04),xlim=(-.48,2.48))
        ax.yaxis.set_major_locator(MultipleLocator(.2))
        ax.grid(axis='y',alpha=.24)
        ax.tick_params(length=0)
    for ax,label in zip(axs[:,0],['Latent-score RV','Latent-score RV','Reconstruction correlation']):
        ax.set_ylabel(label,fontsize=10)
    fig.subplots_adjust(top=.958,bottom=.065,left=.09,right=.982,hspace=.57,wspace=.18)

    # The manuscript PNG was rasterized from the original one-page PDF.
    pdf = BytesIO()
    fig.savefig(pdf, format="pdf", dpi=240)
    plt.close(fig)
    path.parent.mkdir(parents=True, exist_ok=True)
    with fitz.open(stream=pdf.getvalue(), filetype="pdf") as doc:
        doc[0].get_pixmap(dpi=300).save(path.with_suffix(".png"))
