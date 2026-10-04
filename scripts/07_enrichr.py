#!/usr/bin/env python3
# =============================================================================
# 07_enrichr.py — GO / pathway over-representation μέσω Enrichr REST API.
# Χωρίς clusterProfiler/gseapy: μόνο requests + matplotlib (ήδη στο base).
# Είσοδος: λίστες συμβόλων γονιδίων (UP/DOWN/SIG) από το DESeq2.
# Έξοδος: πίνακες CSV + bar plots ανά library.
# =============================================================================
import os, sys, time, json
import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = "/home/haris/Desktop/Bioinformatics_Project"
TBL  = os.path.join(PROJ, "results/DESeq2/tables")
OUT  = os.path.join(PROJ, "results/DESeq2/enrichment")
os.makedirs(OUT, exist_ok=True)

ADDLIST = "https://maayanlab.cloud/Enrichr/addList"
ENRICH  = "https://maayanlab.cloud/Enrichr/enrich"
LIBS = ["GO_Biological_Process_2023","GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023","KEGG_2021_Human","Reactome_2022",
        "MSigDB_Hallmark_2020"]
COLS = ["rank","term","pvalue","zscore","combined_score","genes",
        "adj_pvalue","old_pvalue","old_adj_pvalue"]

def submit(genes, desc):
    r = requests.post(ADDLIST, files={"list": (None,"\n".join(genes)),
                                      "description": (None, desc)}, timeout=60)
    r.raise_for_status()
    return r.json()["userListId"]

def enrich(uid, lib):
    r = requests.get(ENRICH, params={"userListId": uid, "backgroundType": lib}, timeout=60)
    r.raise_for_status()
    data = r.json().get(lib, [])
    df = pd.DataFrame(data, columns=COLS)
    if not df.empty:
        df["n_overlap"] = df["genes"].apply(len)
        df["genes"] = df["genes"].apply(lambda g: ";".join(g))
    return df

def barplot(df, title, path, n=12):
    d = df.sort_values("adj_pvalue").head(n).iloc[::-1]
    if d.empty: return
    import numpy as np
    y = range(len(d)); lp = -np.log10(d["adj_pvalue"].clip(lower=1e-300))
    plt.figure(figsize=(9, max(3, 0.5*len(d))))
    plt.barh(list(y), lp, color="#2471a3")
    plt.yticks(list(y), [t[:60] for t in d["term"]], fontsize=8)
    plt.xlabel("-log10(adjusted p-value)"); plt.title(title, fontsize=10)
    plt.axvline(-1*__import__("math").log10(0.05), ls="--", c="grey", lw=1)
    plt.tight_layout(); plt.savefig(path, dpi=150); plt.close()

def read_list(fn):
    p = os.path.join(TBL, fn)
    if not os.path.exists(p): return []
    return [x.strip() for x in open(p) if x.strip()]

sets = {"UP": read_list("genes_UP_symbols.txt"),
        "DOWN": read_list("genes_DOWN_symbols.txt"),
        "SIG": read_list("genes_SIG_symbols.txt")}

summary = []
for setname, genes in sets.items():
    if len(genes) < 5:
        print(f"[skip] {setname}: {len(genes)} genes (<5)"); continue
    print(f"[{setname}] {len(genes)} genes -> Enrichr")
    uid = submit(genes, f"{setname}_Crohn_vs_Control"); time.sleep(1)
    for lib in LIBS:
        try:
            df = enrich(uid, lib); time.sleep(0.6)
        except Exception as e:
            print(f"   ! {lib}: {e}"); continue
        if df.empty:
            print(f"   {lib}: no terms"); continue
        csv = os.path.join(OUT, f"{setname}__{lib}.csv")
        df.to_csv(csv, index=False)
        sig = df[df["adj_pvalue"] < 0.05]
        print(f"   {lib}: {len(sig)} terms padj<0.05 (of {len(df)})")
        summary.append([setname, lib, len(sig), len(df)])
        if len(sig) > 0:
            barplot(sig, f"{setname} — {lib}",
                    os.path.join(OUT, f"{setname}__{lib}.png"))

pd.DataFrame(summary, columns=["gene_set","library","n_sig_padj0.05","n_total"])\
  .to_csv(os.path.join(OUT,"enrichment_summary.csv"), index=False)
print("ENRICHMENT DONE ->", OUT)
