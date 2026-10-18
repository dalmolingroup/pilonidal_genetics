#!/usr/bin/env python3
"""
Step 4: Query GTEx Portal API v2 for median expression and single-tissue eQTLs
of PSD risk genes in Skin (Suprapubic & Lower Leg), Cultured Fibroblasts, and
Subcutaneous Adipose tissue.
"""
import json
import urllib.request
import urllib.parse
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

TARGET_TISSUES = [
    "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg",
    "Cells_Cultured_fibroblasts",
    "Adipose_Subcutaneous",
    "Whole_Blood",
]


def fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.loads(r.read().decode("utf-8"))


def main():
    curated = pd.read_csv(RESULTS_DIR / "curated_psd_risk_genes.tsv", sep="\t")
    extra_exome = pd.DataFrame([
        {"gene_symbol": "MON1B", "ensembl_id": "ENSG00000103119", "locus": "16q23.1"},
        {"gene_symbol": "IFNA6", "ensembl_id": "ENSG00000120235", "locus": "9p21.3"},
        {"gene_symbol": "WNT9A", "ensembl_id": "ENSG00000143816", "locus": "1q42.13"},
        {"gene_symbol": "MCEMP1", "ensembl_id": "ENSG00000183019", "locus": "19p13.2"},
        {"gene_symbol": "CHST5", "ensembl_id": "ENSG00000135702", "locus": "16q23.1"},
    ])
    genes_df = pd.concat([curated[["gene_symbol", "ensembl_id", "locus"]], extra_exome], ignore_index=True)

    print("[1/2] Fetching GTEx v8 median gene expression across Skin, Fibroblasts, and Subcutaneous Adipose...")
    expr_rows = []
    eqtl_rows = []

    for _, row in genes_df.iterrows():
        sym = row["gene_symbol"]
        ensg = row["ensembl_id"]
        # Resolve versioned gencodeId via GTEx reference gene endpoint
        try:
            g_url = f"https://gtexportal.org/api/v2/reference/gene?geneId={urllib.parse.quote(sym)}&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38"
            g_data = fetch_json(g_url).get("data", [])
            if not g_data:
                continue
            gencode_id = g_data[0]["gencodeId"]
        except Exception as e:
            print(f"  Warning resolving {sym}: {e}")
            continue

        # Median gene expression
        try:
            e_url = f"https://gtexportal.org/api/v2/expression/medianGeneExpression?gencodeId={urllib.parse.quote(gencode_id)}&datasetId=gtex_v8"
            e_data = fetch_json(e_url).get("data", [])
            for item in e_data:
                t_id = item.get("tissueSiteDetailId")
                if t_id in TARGET_TISSUES:
                    expr_rows.append({
                        "gene_symbol": sym,
                        "gencode_id": gencode_id,
                        "tissue": t_id,
                        "median_tpm": item.get("median", 0.0),
                    })
        except Exception as e:
            print(f"  Warning expression {sym}: {e}")

        # Single-tissue eQTLs in Skin / Fibroblasts / Adipose
        for t_id in TARGET_TISSUES[:4]:
            try:
                q_url = f"https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId={urllib.parse.quote(gencode_id)}&tissueSiteDetailId={t_id}&datasetId=gtex_v8"
                q_data = fetch_json(q_url).get("data", [])
                if q_data:
                    q_sorted = sorted(q_data, key=lambda x: float(x.get("pValue", 1.0)))
                    best = q_sorted[0]
                    eqtl_rows.append({
                        "gene_symbol": sym,
                        "gencode_id": gencode_id,
                        "tissue": t_id,
                        "n_significant_eqtls": len(q_data),
                        "lead_snp": best.get("snpId"),
                        "variant_id": best.get("variantId"),
                        "nes": best.get("nes"),
                        "p_value": best.get("pValue"),
                    })
            except Exception:
                pass

    if expr_rows:
        df_e = pd.DataFrame(expr_rows)
        df_piv = df_e.pivot(index="gene_symbol", columns="tissue", values="median_tpm").fillna(0.0)
        out_e = RESULTS_DIR / "gtex_skin_adipose_median_tpm.tsv"
        df_piv.to_csv(out_e, sep="\t")
        print(f"  Saved GTEx median TPM matrix -> {out_e}")
        print(df_piv.round(2))

    if eqtl_rows:
        df_q = pd.DataFrame(eqtl_rows)
        out_q = RESULTS_DIR / "gtex_skin_adipose_lead_eqtls.tsv"
        df_q.to_csv(out_q, sep="\t", index=False)
        print(f"  Saved {len(df_q)} GTEx tissue eQTL summaries -> {out_q}")


if __name__ == "__main__":
    main()
