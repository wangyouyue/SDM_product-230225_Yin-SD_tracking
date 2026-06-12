"""Plot compact wallclock and output-size benchmark summaries."""
import argparse, csv
from pathlib import Path
import matplotlib.pyplot as plt
from nature_style import apply_style, save_figure
def main():
    """Read a summary CSV and write PDF/SVG/PNG plots."""
    ap=argparse.ArgumentParser(); ap.add_argument("csv_path", type=Path); ap.add_argument("--output", type=Path, default=Path("benchmark_summary")); args=ap.parse_args()
    rows=list(csv.DictReader(args.csv_path.open())); names=[r["case_name"] for r in rows]; wall=[float(r.get("wallclock_s") or 0) for r in rows]; out=[float(r.get("total_output_bytes") or 0)/1e6 for r in rows]
    apply_style(); fig,ax=plt.subplots(1,2,figsize=(6.6,2.4)); ax[0].bar(names,wall); ax[0].set_ylabel("wallclock (s)"); ax[1].bar(names,out); ax[1].set_ylabel("output (MB)")
    for a in ax: a.tick_params(axis="x",rotation=35); a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
    fig.tight_layout(); save_figure(fig,args.output)
if __name__ == "__main__": main()
