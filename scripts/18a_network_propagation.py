#!/usr/bin/env python3
"""Step 18a: network propagation of PSD gene scores, following Visona et al. 2024 (Brief Bioinform 25:bbae014) and their code
(github.com/gvisona/GWAS_NetworkPropagation @ 796ba7a, rwr_functions.py / constants.py / ensemble notebook).

- Gene scores: MHC-free MAGMA gene p-values (step 16a), score = max(1e-16, -log10 p) (their load_pegasus_results).
- Networks: their processed NCBI edge lists (BioPlex3, HumanNet v3, PCNet, ProteomeHD, STRING), undirected networkx graphs.
- Propagation: networkx personalized PageRank with the gene scores as personalization (seed 0 for genes without a score),
  over their ALPHAS grid (0, 0.001-0.01, 0.02-0.9; 100 values).
- Ensemble: average of rank / n_nodes across networks (genes present in >= 2 networks), and average of score / (0.01 * sum of seeds).
Outputs: data/netprop/<meta>_<network>.npz (genes x alphas scores), results/netprop_psd_ensemble_top.tsv (top 200 genes per alpha),
         results/netprop_psd_gene_summary.tsv (per-gene ensemble ranks at selected alphas and top-100 frequency across the grid)
"""
import json, os, sys, urllib.request
from multiprocessing import Pool
from pathlib import Path
import numpy as np, pandas as pd
from networkx import pagerank, from_pandas_edgelist

BASE = Path(__file__).resolve().parents[1]
NP = BASE / "data/netprop"; (NP / "networks").mkdir(parents=True, exist_ok=True)
REPO = "https://raw.githubusercontent.com/gvisona/GWAS_NetworkPropagation/796ba7af2158e89fb0aaedf254142a96d59306e6/processed_data/networks"
NETWORKS = {"BioPlex3": "BioPlex3_shared", "HumanNet": "HumanNetV3", "PCNet": "PCNet", "ProteomeHD": "ProteomeHD", "STRING": "STRING"}
ALPHAS = np.concatenate(([0.0], np.linspace(0.001, 0.01, 10), np.linspace(0.02, 0.9, 89)))   # constants.py
METAS = ["PSD_META_EUR", "PSD_META_ALL"]
REPORT_ALPHAS = [0.1, 0.3, 0.5, 0.7, 0.9]

def edges(net):
    f = NP / "networks" / f"{net}.csv"
    if not f.exists():
        urllib.request.urlretrieve(f"{REPO}/{NETWORKS[net]}/edges_list_ncbi.csv", f)
    return f

def scores(meta):
    d = pd.read_csv(BASE / f"data/magma/out_nomhc/{meta}.genes.out", sep=r"\s+")
    return dict(zip(d.GENE.astype(str), np.maximum(1e-16, -np.log10(d.P))))

def run_network(net):
    df = pd.read_csv(edges(net), dtype={"node1": str, "node2": str})[["node1", "node2"]]
    g = from_pandas_edgelist(df, source="node1", target="node2")                     # load_graph_nx
    nodes = list(g.nodes)
    for meta in METAS:
        out = NP / f"{meta}_{net}.npz"
        if out.exists(): continue
        s = scores(meta)
        seeds = {n: s.get(n, 0) for n in nodes}                                       # init_rwr_scores_nx
        mat = np.zeros((len(nodes), len(ALPHAS)), dtype=np.float64)
        for j, a in enumerate(ALPHAS):
            r = pagerank(g, alpha=a, personalization=seeds)                           # perform_rwr_nx
            mat[:, j] = [r[n] for n in nodes]
        np.savez_compressed(out, nodes=np.array(nodes), alphas=ALPHAS, scores=mat, seeds=np.array([seeds[n] for n in nodes]))
        print(f"{meta} {net}: {len(nodes)} nodes, {g.number_of_edges()} edges, {sum(v > 0 for v in seeds.values())} seeded", flush=True)
    return net

if __name__ == "__main__":
    with Pool(int(os.environ.get("NCPU", 5))) as p:
        p.map(run_network, list(NETWORKS))

    loc = pd.read_csv(BASE / "data/magma/ref/NCBI37.3.gene.loc", sep=r"\s+", header=None, names=["GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"])
    sym = dict(zip(loc.GENE.astype(str), loc.SYMBOL))
    l2g = pd.read_csv(BASE / "results/psd_locus_to_gene.tsv", sep="\t")
    locus_of = {}
    for r in l2g.itertuples():
        for gname in str(r.genes_250kb).split(",") + [str(r.nearest_gene)] + str(r.prioritised_gene).split("/"):
            locus_of.setdefault(gname, r.locus)
    tops, summ = [], []
    for meta in METAS:
        magma = pd.read_csv(BASE / f"data/magma/out_nomhc/{meta}.genes.out", sep=r"\s+")
        mp = dict(zip(magma.GENE.astype(str), magma.P))
        ranks, nsc = {}, {}
        for net in NETWORKS:
            z = np.load(NP / f"{meta}_{net}.npz")
            nodes, mat, seeds = z["nodes"], z["scores"], z["seeds"]
            order = np.argsort(-mat, axis=0, kind="stable")                           # rank per alpha (1 = best)
            rk = np.empty_like(order); rk[order, np.arange(mat.shape[1])] = np.arange(1, len(nodes) + 1)[:, None]
            ranks[net] = pd.DataFrame(rk / len(nodes), index=nodes)                   # rank / n_nodes
            nsc[net] = pd.DataFrame(mat / (seeds.sum() * 0.01), index=nodes)         # score / (0.01 * sum of seeds)
        R = pd.concat(ranks, axis=1); S = pd.concat(nsc, axis=1)
        n_in = R.xs(0, axis=1, level=1).notna().sum(axis=1)
        keep = n_in[n_in >= 2].index
        avg_rank = pd.concat({a: R.xs(j, axis=1, level=1).loc[keep].mean(axis=1) for j, a in enumerate(ALPHAS)}, axis=1)
        avg_score = pd.concat({a: S.xs(j, axis=1, level=1).loc[keep].mean(axis=1) for j, a in enumerate(ALPHAS)}, axis=1)
        top100 = (avg_rank.rank(axis=0, method="min") <= 100)
        for j, a in enumerate(ALPHAS):
            t = avg_rank[a].nsmallest(200)
            tops.append(pd.DataFrame({"meta": meta, "alpha": round(a, 3), "ensemble_rank": np.arange(1, len(t) + 1), "gene_id": t.index,
                                      "symbol": [sym.get(x, "-") for x in t.index], "avg_rank_norm": t.values,
                                      "avg_score_norm": avg_score.loc[t.index, a].values, "magma_p": [mp.get(x, np.nan) for x in t.index],
                                      "n_networks": n_in.loc[t.index].values, "psd_locus": [locus_of.get(sym.get(x, "-"), "") for x in t.index]}))
        s = pd.DataFrame({"meta": meta, "gene_id": keep, "symbol": [sym.get(x, "-") for x in keep], "magma_p": [mp.get(x, np.nan) for x in keep],
                          "n_networks": n_in.loc[keep].values, "psd_locus": [locus_of.get(sym.get(x, "-"), "") for x in keep],
                          "top100_frac_alpha_grid": top100.loc[keep].mean(axis=1).values})
        for a in REPORT_ALPHAS:
            j = int(np.argmin(np.abs(ALPHAS - a)))
            s[f"ens_rank_a{a}"] = avg_rank[ALPHAS[j]].rank(method="min").loc[keep].values
        summ.append(s)
    tops = pd.concat(tops); summ = pd.concat(summ).sort_values(["meta", "ens_rank_a0.5"])
    tops.to_csv(BASE / "results/netprop_psd_ensemble_top.tsv", sep="\t", index=False)
    summ.to_csv(BASE / "results/netprop_psd_gene_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    for meta in METAS:
        for a in [0.5, 0.9]:
            t = tops[(tops.meta == meta) & (np.isclose(tops.alpha, a))].head(40)
            print(f"\n== {meta} alpha {a}: top 40 ensemble genes (magma_p; PSD locus)")
            print(t[["ensemble_rank", "symbol", "magma_p", "n_networks", "psd_locus"]].to_string(index=False))
        s = summ[summ.meta == meta]
        new = s[(s.psd_locus == "") & (s.magma_p > 1e-3)].sort_values("ens_rank_a0.5").head(30)
        print(f"\n== {meta}: top genes outside PSD loci with MAGMA p > 1e-3 (ensemble rank, alpha 0.5)")
        print(new[["symbol", "magma_p", "n_networks", "top100_frac_alpha_grid"] + [f"ens_rank_a{a}" for a in REPORT_ALPHAS]].to_string(index=False))
    print("NETPROP_DONE")
