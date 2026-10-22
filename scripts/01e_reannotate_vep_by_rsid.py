#!/usr/bin/env python3
"""
Step 01e (lightweight): re-annotate fine-mapped credible-set variants with Ensembl VEP by rsID, so VEP uses the true
GRCh38 reference allele (01c passed alphabetically-sorted alleles, which mislabels REF for ~half the variants).
Adds CADD-free consequence, amino-acid change, regulatory feature and JASPAR motif disruption.
Input : results/psd_finemap_credible_sets.tsv   Output: results/psd_finemap_vep_by_rsid.tsv
"""
import json, time, urllib.request
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
cs = pd.read_csv(BASE / "results" / "psd_finemap_credible_sets.tsv", sep="\t")
cs["pip_max"] = cs[["abf_pip", "susie_pip"]].max(axis=1)
sel = cs[(cs.abf_cs95 == True) | (cs.susie_cs.notna()) | (cs.pip_max > 0.02)]
sel = sel[sel.rsid.astype(str).str.startswith("rs")].sort_values("pip_max", ascending=False).drop_duplicates("rsid")
sel = sel.groupby("locus").head(40)
ids = sel.rsid.tolist()
print(f"{len(ids)} rsIDs to annotate", flush=True)
rows = []
for i in range(0, len(ids), 150):
    payload = json.dumps({"ids": ids[i:i + 150]}).encode()
    for attempt in range(4):
        try:
            req = urllib.request.Request("https://rest.ensembl.org/vep/human/id?regulatory=1&canonical=1&symbol=1&mane=1",
                                         data=payload, headers={"Content-Type": "application/json", "Accept": "application/json"})
            res = json.loads(urllib.request.urlopen(req, timeout=180).read()); break
        except Exception as e:
            print("  retry", attempt, e, flush=True); time.sleep(5 * (attempt + 1)); res = []
    for v in res:
        tcs = v.get("transcript_consequences", []) or []
        canon = [t for t in tcs if t.get("canonical")] or tcs
        aa = ";".join(sorted({f"{t.get('gene_symbol')}:{t.get('amino_acids')}{t.get('protein_start','')}" for t in canon if t.get("amino_acids")}))
        motifs = v.get("motif_feature_consequences", []) or []
        rows.append({"rsid": v.get("id"), "allele_string": v.get("allele_string"), "most_severe_consequence": v.get("most_severe_consequence"),
                     "genes_canonical": ",".join(dict.fromkeys(t.get("gene_symbol", "") for t in canon if t.get("gene_symbol"))),
                     "aa_change": aa, "sift": ";".join(sorted({str(t.get('sift_prediction')) for t in canon if t.get('sift_prediction')})),
                     "polyphen": ";".join(sorted({str(t.get('polyphen_prediction')) for t in canon if t.get('polyphen_prediction')})),
                     "regulatory_features": ";".join(sorted({f"{r.get('biotype')}" for r in (v.get("regulatory_feature_consequences") or [])})),
                     "tf_motifs": ";".join(f"{(m.get('transcription_factors') or [m.get('motif_name','')])[0]}({m.get('motif_score_change')}{',HI' if m.get('high_inf_pos')=='Y' else ''})" for m in motifs[:6])})
    time.sleep(1)
ann = pd.DataFrame(rows).drop_duplicates("rsid")
out = sel[["locus", "meta", "chr", "pos", "rsid", "a1", "a2", "eaf", "p", "abf_pip", "susie_pip", "susie_cs", "pip_max"]].merge(ann, on="rsid", how="left")
out.to_csv(BASE / "results" / "psd_finemap_vep_by_rsid.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(out[out.pip_max >= 0.2][["locus", "rsid", "allele_string", "pip_max", "most_severe_consequence", "genes_canonical", "aa_change", "sift", "polyphen", "regulatory_features", "tf_motifs"]].to_string(index=False))
print("\ncoding variants in credible sets:")
print(out[out.aa_change.fillna("") != ""][["locus", "rsid", "pip_max", "most_severe_consequence", "aa_change", "sift", "polyphen"]].to_string(index=False))
