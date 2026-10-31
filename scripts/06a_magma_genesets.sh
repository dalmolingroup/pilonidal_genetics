#!/bin/bash
# Step 06a: MAGMA v1.10 competitive gene-set analyses on MHC-free gene scores.
# (1) MHC genes (GRCh37 chr6:28.4-33.5 Mb, NCBI37.3 gene loc) are dropped from the step-03a SNP->gene annotation and the gene
#     analysis is rerun for every trait, so MHC genes are absent from the gene universe (not only from the sets);
# (2) gene sets: MSigDB 2025.1 GO biological process + Reactome (Entrez), 10-500 genes; a priori follicle/keratinocyte/WNT sets flagged;
# (3) robustness for PSD: top sets rerun without their strongest locus (top gene +/-1 Mb) and without genome-wide significant genes.
# Inputs : data/magma/in/<trait>.txt and data/magma/ref/* (step 03a)
# Outputs: results/magma_genesets_nomhc.tsv, results/magma_genesets_psd_robustness.tsv
set -euo pipefail
BASE=$(cd "$(dirname "$0")/.." && pwd)
PY=${PY:-python3}
M=$BASE/data/magma; R=$M/ref; O=$M/out_nomhc; G=$M/genesets; mkdir -p $O $G
MAGMA=$R/magma

# ---- (1) MHC-free annotation + gene analysis ----
$PY - <<PYEOF
import pandas as pd
R = "$R"
loc = pd.read_csv(f"{R}/NCBI37.3.gene.loc", sep=r"\s+", header=None, names=["GENE","CHR","START","STOP","STRAND","SYMBOL"])
mhc = set(loc[(loc.CHR.astype(str) == "6") & (loc.STOP > 28_400_000) & (loc.START < 33_500_000)].GENE.astype(str))
n = 0
with open(f"{R}/annot.genes.annot") as fi, open(f"{R}/annot_noMHC.genes.annot", "w") as fo:
    for line in fi:
        if not line.startswith("#") and line.split("\t", 1)[0] in mhc: n += 1; continue
        fo.write(line)
print(f"annotation: removed {n} MHC genes")
PYEOF
run_gene() { t=$1; [ -s $O/$t.genes.raw ] && return; $MAGMA --bfile $R/g1000_eur --pval $M/in/$t.txt use=SNP,P ncol=N \
    --gene-annot $R/annot_noMHC.genes.annot --out $O/$t > $O/$t.stdout 2>&1 || echo "MAGMA gene analysis failed $t"; }
export -f run_gene; export M R O MAGMA
TRAITS=$(ls $M/in/*.txt | xargs -n1 basename | sed 's/\.txt$//')
echo "$TRAITS" | xargs -P 7 -I{} bash -c 'run_gene {}'

# ---- (2) gene sets ----
V=2025.1.Hs; U=https://data.broadinstitute.org/gsea-msigdb/msigdb/release/$V
for c in c5.go.bp c2.cp.reactome; do [ -s $G/$c.gmt ] || curl -sSL --retry 5 -o $G/$c.gmt $U/$c.v$V.entrez.gmt; done
$PY - <<PYEOF
G = "$G"
n = 0
with open(f"{G}/msigdb_gobp_reactome.sets", "w") as fo:
    for c in ("c5.go.bp", "c2.cp.reactome"):
        for line in open(f"{G}/{c}.gmt"):
            f = line.rstrip("\n").split("\t")
            if 10 <= len(f) - 2 <= 500: fo.write(" ".join([f[0]] + f[2:]) + "\n"); n += 1
print(f"{n} gene sets (10-500 genes)")
PYEOF

run_sets() { t=$1
  [ -s $O/$t.gsa.out ] || $MAGMA --gene-results $O/$t.genes.raw --set-annot $G/msigdb_gobp_reactome.sets --out $O/$t > $O/$t.sets.stdout 2>&1 || echo "set analysis failed $t"; }
export -f run_sets; export G
echo "$TRAITS" | xargs -P 7 -I{} bash -c 'run_sets {}'

# ---- summaries + (3) PSD robustness sets ----
$PY - <<PYEOF
import pandas as pd, numpy as np, glob, os
from statsmodels.stats.multitest import multipletests
B, R, O, G = "$BASE", "$R", "$O", "$G"
APRIORI = ["GOBP_MOLTING_CYCLE", "GOBP_MOLTING_CYCLE_PROCESS", "GOBP_KERATINOCYTE_DIFFERENTIATION",   # MSigDB 2025.1 has no HAIR_CYCLE/HAIR_FOLLICLE_DEVELOPMENT; molting cycle = GO hair cycle "GOBP_KERATINIZATION",
           "GOBP_EPIDERMIS_DEVELOPMENT", "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY", "GOBP_WNT_SIGNALING_PATHWAY",
           "REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE", "REACTOME_KERATINIZATION", "REACTOME_SIGNALING_BY_WNT"]
def read_gsa(f):
    d = pd.read_csv(f, sep=r"\s+", comment="#")
    if "FULL_NAME" in d: d["VARIABLE"] = d.FULL_NAME
    return d.drop(columns=[c for c in ("FULL_NAME", "TYPE") if c in d])
sets = []
for f in sorted(glob.glob(f"{O}/*.gsa.out")):
    t = os.path.basename(f).replace(".gsa.out", "")
    if t.endswith("_robust") or t.endswith("_gtex"): continue
    d = read_gsa(f)
    d.insert(0, "trait", t); d["collection"] = np.where(d.VARIABLE.str.startswith("GOBP_"), "GO_BP", "Reactome")
    d["FDR"] = multipletests(d.P, method="fdr_bh")[1]; d["P_bonf"] = np.minimum(d.P * len(d), 1)
    d["a_priori"] = d.VARIABLE.isin(APRIORI)
    d["P_bonf_apriori"] = np.where(d.a_priori, np.minimum(d.P * d.a_priori.sum(), 1), np.nan); sets.append(d)
sets = pd.concat(sets).rename(columns={"VARIABLE": "gene_set"})
sets.to_csv(f"{B}/results/magma_genesets_nomhc.tsv", sep="\t", index=False)
for t in sorted(sets.trait.unique()):
    s = sets[sets.trait == t].sort_values("P")
    print(f"\n== {t}: {len(s)} sets, FDR<0.05 {(s.FDR < 0.05).sum()}, Bonferroni {(s.P_bonf < 0.05).sum()}")
    print(s.head(8)[["gene_set", "NGENES", "BETA", "P", "FDR"]].to_string(index=False))
print("\n== a priori sets (Bonferroni over", len(APRIORI), ")")
print(sets[sets.a_priori].pivot(index="gene_set", columns="trait", values="P").map(lambda p: f"{p:.1e}").to_string())

# robustness: PSD top sets (FDR<0.05 or top 10, plus a priori) without the strongest locus / without GWS genes
loc = pd.read_csv(f"{R}/NCBI37.3.gene.loc", sep=r"\s+", header=None, names=["GENE","CHR","START","STOP","STRAND","SYMBOL"]).set_index("GENE")
members = {}
for line in open(f"{G}/msigdb_gobp_reactome.sets"):
    f = line.split(); members[f[0]] = [int(x) for x in f[1:]]
with open(f"{G}/psd_robust.sets", "w") as fo, open(f"{G}/psd_robust_drivers.tsv", "w") as fd:
    fd.write("trait\tgene_set\ttop_gene\ttop_gene_p\tn_dropped_locus\tn_dropped_gws\n")
    for t in ("PSD_META_EUR", "PSD_META_ALL"):
        gp = pd.read_csv(f"{O}/{t}.genes.out", sep=r"\s+").set_index("GENE")
        s = sets[sets.trait == t].sort_values("P")
        for name in sorted(set(s[s.FDR < 0.05].gene_set) | set(s.head(10).gene_set) | set(APRIORI)):
            m = [g for g in members.get(name, []) if g in gp.index]
            if not m: continue
            top = gp.loc[m, "P"].idxmin(); c, p0 = loc.loc[top, "CHR"], loc.loc[top, "START"]
            near = [g for g in m if loc.loc[g, "CHR"] == c and abs(loc.loc[g, "START"] - p0) < 1e6]
            gws = [g for g in m if gp.loc[g, "P"] < 2.5e-6]
            fo.write(" ".join([f"{t}__noTopLocus__{name}"] + [str(g) for g in m if g not in near]) + "\n")
            fo.write(" ".join([f"{t}__noGWS__{name}"] + [str(g) for g in m if g not in gws]) + "\n")
            fd.write(f"{t}\t{name}\t{loc.loc[top, 'SYMBOL']}\t{gp.loc[top, 'P']:.2e}\t{len(near)}\t{len(gws)}\n")
PYEOF
for t in PSD_META_EUR PSD_META_ALL; do
  grep "^${t}__" $G/psd_robust.sets > $G/${t}_robust.sets
  $MAGMA --gene-results $O/$t.genes.raw --set-annot $G/${t}_robust.sets --out $O/${t}_robust > $O/$t.robust.stdout 2>&1 || echo "robustness failed $t"
done
$PY - <<PYEOF
import pandas as pd
B, O, G = "$BASE", "$O", "$G"
full = pd.read_csv(f"{B}/results/magma_genesets_nomhc.tsv", sep="\t")
drv = pd.read_csv(f"{G}/psd_robust_drivers.tsv", sep="\t")
out = []
for t in ("PSD_META_EUR", "PSD_META_ALL"):
    d = pd.read_csv(f"{O}/{t}_robust.gsa.out", sep=r"\s+", comment="#")
    if "FULL_NAME" in d: d["VARIABLE"] = d.FULL_NAME
    d[["trait", "variant", "gene_set"]] = d.VARIABLE.str.split("__", n=2, expand=True)
    w = d.pivot(index=["trait", "gene_set"], columns="variant", values="P").add_prefix("P_").reset_index()
    out.append(w)
r = pd.concat(out).merge(full[["trait", "gene_set", "NGENES", "BETA", "P", "FDR", "a_priori"]], on=["trait", "gene_set"]).merge(drv, on=["trait", "gene_set"])
r = r.sort_values(["trait", "P"])
r.to_csv(f"{B}/results/magma_genesets_psd_robustness.tsv", sep="\t", index=False)
print("\n== PSD robustness (P full / without top locus / without GWS genes)")
print(r[["trait", "gene_set", "NGENES", "P", "FDR", "P_noTopLocus", "P_noGWS", "top_gene", "n_dropped_locus"]].to_string(index=False))
PYEOF
echo MAGMA_GENESETS_DONE
