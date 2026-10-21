#!/usr/bin/env python3
"""Step 18b: go/no-go validation of the step-18a network propagation, replacing Visona et al.'s GWAS Catalog benchmark.

Pre-set design (fixed before looking at any result):
- Chromosome split: seed the propagation with MHC-free MAGMA scores of autosomal genes on odd chromosomes only (all other genes seed 0),
  and evaluate the ensemble ranking of genes on even chromosomes; and vice versa. Same networks, propagation and ensemble as 18a.
- Targets on the held-out chromosomes: genes with MAGMA p < 1e-3 (same meta) plus the prioritised genes of the 25 PSD loci
  (results/psd_locus_to_gene.tsv). Universe: held-out autosomal genes in the ensemble (present in >= 2 networks).
- Metric: average precision (AP) of the targets in the held-out ranking; statistic = mean AP over the two folds.
- Propagation strength: alpha is chosen by the cross-check (argmax of the statistic over the 18a grid, alpha = 0 excluded).
- Connectivity baseline: genes ranked by degree (ensemble of degree rank / n_nodes across networks, same rule as the propagation ensemble).
- Null: 1,000 degree-matched permutations of the seed scores (training-chromosome genes with a MAGMA score, shuffled within
  20 quantile bins of ensemble degree rank; same shuffle in all networks). The null statistic is also maximised over alpha,
  so the permutation p accounts for choosing alpha.
- Rule (primary meta PSD_META_EUR; ALL reported alongside): GO only if the statistic beats the degree baseline AND permutation p < 0.05.
  Otherwise report a null supplementary result and stop.
Deviation from 18a: per-network ranks use average ranks for ties (18a used ordinal), so genes tied at zero score are not ordered arbitrarily.
PageRank is a vectorised re-implementation of networkx's _pagerank_scipy (same matrix, tolerance and update); it is checked against the
18a output before use.
Outputs: results/netprop_validation_curve.tsv (per meta/fold/alpha: observed AP, null mean and 95th percentile),
         results/netprop_validation_summary.tsv (selected alpha, statistic, degree baseline, permutation p, verdict)
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
from multiprocessing import Pool
from pathlib import Path
import numpy as np, pandas as pd
from networkx import from_pandas_edgelist, to_scipy_sparse_array, pagerank
from scipy.stats import rankdata

BASE = Path(__file__).resolve().parents[1]
NP = BASE / "data/netprop"
NETWORKS = ["BioPlex3", "HumanNet", "PCNet", "ProteomeHD", "STRING"]
ALPHAS = np.concatenate(([0.0], np.linspace(0.001, 0.01, 10), np.linspace(0.02, 0.9, 89)))   # same grid as 18a
METAS = ["PSD_META_EUR", "PSD_META_ALL"]
FOLDS = {"odd_to_even": 1, "even_to_odd": 0}   # value = parity of the seeded chromosomes
NPERM, NBINS, P_TARGET, SEED = 1000, 20, 1e-3, 18

loc = pd.read_csv(BASE / "data/magma/ref/NCBI37.3.gene.loc", sep=r"\s+", header=None, names=["GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"])
loc["GENE"] = loc.GENE.astype(str)
CHR = {g: int(c) for g, c in zip(loc.GENE, loc.CHR) if str(c).isdigit()}            # autosomes only
SYM2ID = dict(zip(loc.SYMBOL, loc.GENE))

def ppr(AT, P, alpha, tol=1e-6, max_iter=100):
    """networkx _pagerank_scipy, one column per personalization vector (columns of P sum to 1, no dangling nodes)."""
    N = P.shape[0]
    X = np.full_like(P, 1.0 / N)
    for _ in range(max_iter):
        Xn = alpha * (AT @ X) + (1 - alpha) * P
        done = np.abs(Xn - X).sum(axis=0).max() < N * tol
        X = Xn
        if done: return X
    print(f"WARNING: not converged at alpha {alpha}", flush=True)
    return X

# graphs: transposed row-stochastic adjacency, exactly as networkx builds it
GRAPHS = {}
for net in NETWORKS:
    df = pd.read_csv(NP / "networks" / f"{net}.csv", dtype={"node1": str, "node2": str})[["node1", "node2"]]
    g = from_pandas_edgelist(df, source="node1", target="node2")
    nodes = np.array(list(g.nodes))
    A = to_scipy_sparse_array(g, nodelist=list(nodes), dtype=float)
    S = np.asarray(A.sum(axis=1)).ravel()
    assert (S > 0).all()
    A = (A.multiply(1.0 / S[:, None])).tocsr()
    GRAPHS[net] = dict(nodes=nodes, AT=A.T.tocsr(), deg=S)
    del g                                                                              # keep forked workers light

# ensemble universe (>= 2 networks) and degree rank ensemble
n_in = pd.Series(np.concatenate([GRAPHS[n]["nodes"] for n in NETWORKS])).value_counts()
ENS = set(n_in[n_in >= 2].index)
deg_rank = pd.concat({n: pd.Series(rankdata(-GRAPHS[n]["deg"]) / len(GRAPHS[n]["nodes"]), index=GRAPHS[n]["nodes"]) for n in NETWORKS}, axis=1)
deg_ens = deg_rank.mean(axis=1)                                                        # mean over the networks a gene is in

l2g = pd.read_csv(BASE / "results/psd_locus_to_gene.tsv", sep="\t")
LOCUS_GENES = {SYM2ID[s] for p in l2g.prioritised_gene for s in (x.replace(" (lncRNA)", "").strip() for x in p.split("/")) if s in SYM2ID}

def magma(meta):
    d = pd.read_csv(BASE / f"data/magma/out_nomhc/{meta}.genes.out", sep=r"\s+")
    return pd.Series(d.P.values, index=d.GENE.astype(str))

# per (meta, fold): gene-level seed matrix (observed in column 0, degree-matched permutations after), held-out universe and targets
SETUP = {}
rng = np.random.default_rng(SEED)
for meta in METAS:
    p = magma(meta)
    score = np.maximum(1e-16, -np.log10(p))
    for fold, par in FOLDS.items():
        train = [g for g in score.index if g in CHR and CHR[g] % 2 == par and g in deg_ens.index]
        bins = pd.qcut(deg_ens.loc[train].rank(method="first"), NBINS, labels=False).values
        base = score.loc[train].values
        mat = np.empty((len(train), NPERM + 1)); mat[:, 0] = base
        for k in range(1, NPERM + 1):
            col = base.copy()
            for b in range(NBINS):
                idx = np.where(bins == b)[0]; col[idx] = base[rng.permutation(idx)]
            mat[:, k] = col
        held = np.array(sorted(g for g in ENS if g in CHR and CHR[g] % 2 != par))
        y = np.array([(p.get(g, 1.0) < P_TARGET) or (g in LOCUS_GENES) for g in held])
        SETUP[(meta, fold)] = dict(row={g: i for i, g in enumerate(train)}, mat=mat, held=held, y=y)

def average_precision(rank_mat, y):
    """AP per column; lower rank value = better."""
    order = np.argsort(rank_mat, axis=0, kind="stable")
    ys = y[order]
    hits = np.cumsum(ys, axis=0)
    prec = hits / np.arange(1, len(y) + 1)[:, None]
    return (prec * ys).sum(axis=0) / y.sum()

def held_ensemble(meta, fold, per_net):
    """Average rank / n_nodes over the networks each held-out gene is in."""
    st = SETUP[(meta, fold)]
    hidx = {g: i for i, g in enumerate(st["held"])}
    acc = None; cnt = np.zeros(len(st["held"]))
    for net, R in per_net.items():
        nodes = GRAPHS[net]["nodes"]
        sel = np.array([i for i, n in enumerate(nodes) if n in hidx]); tgt = np.array([hidx[n] for n in nodes[sel]])
        if acc is None: acc = np.zeros((len(st["held"]), R.shape[1]))
        acc[tgt] += R[sel]; cnt[tgt] += 1
    return acc / cnt[:, None]

def task(args):
    meta, fold, j = args
    st = SETUP[(meta, fold)]
    per_net = {}
    for net in NETWORKS:
        G = GRAPHS[net]
        P = np.zeros((len(G["nodes"]), NPERM + 1))
        rows = [(i, st["row"][n]) for i, n in enumerate(G["nodes"]) if n in st["row"]]
        P[[a for a, _ in rows]] = st["mat"][[b for _, b in rows]]
        P /= P.sum(axis=0)
        X = ppr(G["AT"], P, ALPHAS[j])
        per_net[net] = (rankdata(-X, axis=0) / len(G["nodes"])).astype(np.float32)
    return meta, fold, j, average_precision(held_ensemble(meta, fold, per_net), st["y"])

if __name__ == "__main__":
    # check the vectorised PageRank against networkx and the 18a output
    G = GRAPHS["STRING"]; s = magma("PSD_META_EUR"); sc = np.maximum(1e-16, -np.log10(s))
    seeds = np.array([sc.get(n, 0.0) for n in G["nodes"]])
    z = np.load(NP / "PSD_META_EUR_STRING.npz")
    df = pd.read_csv(NP / "networks/STRING.csv", dtype={"node1": str, "node2": str})
    g = from_pandas_edgelist(df, source="node1", target="node2")
    assert (z["nodes"] == G["nodes"]).all()
    for a in (0.1, 0.5, 0.9):
        mine = ppr(G["AT"], (seeds / seeds.sum())[:, None], a)[:, 0]
        nx_ = pagerank(g, alpha=a, personalization=dict(zip(G["nodes"], seeds)))
        ref = z["scores"][:, int(np.argmin(np.abs(ALPHAS - a)))]
        print(f"PPR check alpha {a}: max |mine - networkx| {np.abs(mine - np.array([nx_[n] for n in G['nodes']])).max():.2e}, "
              f"max |mine - 18a| {np.abs(mine - ref).max():.2e}", flush=True)
    for (meta, fold), st in SETUP.items():
        print(f"{meta} {fold}: {len(st['row'])} seeded genes, {len(st['held'])} held-out genes, {st['y'].sum()} targets "
              f"({sum(g in LOCUS_GENES for g in st['held'][st['y']])} locus genes)", flush=True)

    jobs = [(m, f, j) for m in METAS for f in FOLDS for j in range(1, len(ALPHAS))]
    jobs.sort(key=lambda t: -ALPHAS[t[2]])                                            # slow (high alpha) tasks first
    AP = {}
    with Pool(int(os.environ.get("NCPU", 20))) as pool:
        for n, (m, f, j, ap) in enumerate(pool.imap_unordered(task, jobs), 1):
            AP[(m, f, j)] = ap
            if n % 20 == 0: print(f"{n}/{len(jobs)} tasks done", flush=True)

    curve, summ = [], []
    for meta in METAS:
        deg_ap = {}
        for fold in FOLDS:
            st = SETUP[(meta, fold)]
            deg_ap[fold] = average_precision(deg_ens.loc[st["held"]].values[:, None], st["y"])[0]
            for j in range(1, len(ALPHAS)):
                ap = AP[(meta, fold, j)]
                curve.append(dict(meta=meta, fold=fold, alpha=round(ALPHAS[j], 3), ap_obs=ap[0], ap_null_mean=ap[1:].mean(),
                                  ap_null_q95=np.quantile(ap[1:], 0.95), p_perm_this_alpha=(1 + (ap[1:] >= ap[0]).sum()) / (NPERM + 1),
                                  ap_degree=deg_ap[fold], prevalence=st["y"].mean()))
        T = np.mean([np.vstack([AP[(meta, f, j)] for j in range(1, len(ALPHAS))]) for f in FOLDS], axis=0)   # alphas x (1 + NPERM)
        best = int(np.argmax(T[:, 0])); t_obs = T[best, 0]; t_null = T[:, 1:].max(axis=0)
        p_perm = (1 + (t_null >= t_obs).sum()) / (NPERM + 1); deg = np.mean(list(deg_ap.values()))
        verdict = "GO" if (t_obs > deg and p_perm < 0.05) else "NO-GO"
        summ.append(dict(meta=meta, selected_alpha=round(ALPHAS[best + 1], 3), mean_ap=t_obs,
                         **{f"ap_{f}": AP[(meta, f, best + 1)][0] for f in FOLDS}, mean_ap_degree=deg,
                         **{f"ap_degree_{f}": deg_ap[f] for f in FOLDS}, null_max_mean=t_null.mean(), null_max_q95=np.quantile(t_null, 0.95),
                         p_perm_alpha_adjusted=p_perm, mean_prevalence=np.mean([SETUP[(meta, f)]["y"].mean() for f in FOLDS]),
                         n_targets=sum(SETUP[(meta, f)]["y"].sum() for f in FOLDS), verdict=verdict))
    curve = pd.DataFrame(curve); summ = pd.DataFrame(summ)
    curve.to_csv(BASE / "results/netprop_validation_curve.tsv", sep="\t", index=False)
    summ.to_csv(BASE / "results/netprop_validation_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(summ.T.to_string())
    c = curve[curve.alpha.isin([0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9])]
    print(c.drop(columns=["prevalence"]).to_string(index=False))
    print("NETPROP_VALIDATION_DONE")
