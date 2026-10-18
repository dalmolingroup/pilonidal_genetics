#!/usr/bin/env python3
"""
Step 7: Inspect and score GEO Hair Follicle Mesenchymal scRNA-seq (GSE287304: Dermal Papilla,
Dermal Sheath, and Interfollicular Dermal Fibroblasts) and Hidradenitis Suppurativa / Scalp
scRNA-seq to dissect which exact hair-follicle mesenchymal sub-compartment expresses the
PSD Mesenchymal/Wnt risk module (TWIST1, HDAC9, DKK1, TBX15, WARS2, PDGFA, WNT9A).
"""
import urllib.request
import json
import re
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"


def check_geo_supp_files(gse_list):
    for gse in gse_list:
        prefix = gse[:-3] + "nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/"
        try:
            html = urllib.request.urlopen(url, timeout=15).read().decode("utf-8", errors="ignore")
            files = re.findall(r'href="([^"?/]+)"[^>]*>[^<]+</a>\s+(\d+-\d+-\d+\s+\d+:\d+)\s+(\S+)', html)
            print(f"{gse} ({url}): {files}")
        except Exception as e:
            print(f"{gse}: {e}")


if __name__ == "__main__":
    check_geo_supp_files(["GSE287304", "GSE274880", "GSE158955", "GSE212447", "GSE273109"])
