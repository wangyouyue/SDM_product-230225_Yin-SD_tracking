"""Matplotlib style helpers for publication-quality GMD2026 figures."""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
OKABE_ITO = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000"]
def apply_style():
    """Apply a compact, color-blind-safe Matplotlib style."""
    mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.labelsize":8,"axes.titlesize":9,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,"axes.linewidth":0.8,"xtick.direction":"in","ytick.direction":"in","lines.linewidth":1.0,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","savefig.bbox":"tight"})
    plt.rcParams["axes.prop_cycle"] = mpl.cycler(color=OKABE_ITO)
def save_figure(fig, output_base, dpi=600):
    """Save PDF, SVG, and high-resolution PNG versions of a figure."""
    base = Path(output_base); base.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".svg", ".png"):
        fig.savefig(base.with_suffix(suffix), dpi=dpi if suffix == ".png" else None)
