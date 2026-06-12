"""File-size collectors for GMD2026 post-processing."""
from pathlib import Path
def collect_output_inventory(case_dir):
    """Return output byte counts and file counts for one model case."""
    pats = ["SD_selected_NetCDF_*", "SD_all_NetCDF_*", "SD_coal_output_NetCDF_*", "history*.nc"]
    files = [p for pat in pats for p in Path(case_dir).glob(pat) if p.is_file()]
    total = sum(p.stat().st_size for p in files); count = len(files)
    return {"total_output_bytes":total,"number_of_output_files":count,"mean_file_size_bytes":total/count if count else 0.0,"sd_selected_output_bytes":sum(p.stat().st_size for p in Path(case_dir).glob("SD_selected_NetCDF_*") if p.is_file()),"sd_all_output_bytes":sum(p.stat().st_size for p in Path(case_dir).glob("SD_all_NetCDF_*") if p.is_file()),"coalescence_log_bytes":sum(p.stat().st_size for p in Path(case_dir).glob("SD_coal_output_NetCDF_*") if p.is_file()),"tpht_id_bytes":sum(p.stat().st_size for p in Path(case_dir).glob("**/*.ids") if p.is_file())}
