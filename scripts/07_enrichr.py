#!/usr/bin/env python3
"""GO and pathway over-representation analysis with the Enrichr API.

Input:      the gene-symbol lists written by 06_deseq2.R (UP, DOWN, SIG).
Background: every gene tested by DESeq2 (after pre-filtering), sent through
            Enrichr's background-enrichment API. Only genes expressed in this
            tissue can be drawn, which avoids the inflated significance of a
            whole-genome background.
Output:     one CSV per gene set and library, plus a bar plot when at least one
            term has adjusted p < 0.05, and enrichment_summary.csv.
"""

import math
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # render without a display

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402

PROJ = Path(__file__).resolve().parents[1]
TABLES = PROJ / "results/DESeq2/tables"
OUT = PROJ / "results/DESeq2/enrichment"
API = "https://maayanlab.cloud/speedrichr/api"
LIBRARIES = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
]
COLUMNS = ["rank", "term", "pvalue", "odds_ratio", "combined_score", "genes", "adj_pvalue"]


def post(endpoint: str, **kwargs) -> dict:
    response = requests.post(f"{API}/{endpoint}", timeout=120, **kwargs)
    response.raise_for_status()
    return response.json()


def add_list(genes: list[str], description: str) -> int:
    files = {"list": (None, "\n".join(genes)), "description": (None, description)}
    return post("addList", files=files)["userListId"]


def add_background(genes: list[str]) -> str:
    return post("addbackground", data={"background": "\n".join(genes)})["backgroundid"]


def enrich(list_id: int, background_id: str, library: str) -> pd.DataFrame:
    data = {"userListId": list_id, "backgroundid": background_id, "backgroundType": library}
    rows = post("backgroundenrich", data=data).get(library, [])
    table = pd.DataFrame([row[: len(COLUMNS)] for row in rows], columns=COLUMNS)
    if not table.empty:
        table["n_overlap"] = table["genes"].apply(len)
        table["genes"] = table["genes"].apply(";".join)
    return table


def barplot(table: pd.DataFrame, title: str, path: Path, n: int = 12) -> None:
    top = table.sort_values("adj_pvalue").head(n).iloc[::-1]
    y = range(len(top))
    plt.figure(figsize=(9, max(3, 0.5 * len(top))))
    plt.barh(list(y), -np.log10(top["adj_pvalue"].clip(lower=1e-300)), color="#2471a3")
    plt.yticks(list(y), [term[:60] for term in top["term"]], fontsize=8)
    plt.xlabel("-log10(adjusted p-value)")
    plt.title(title, fontsize=10)
    plt.axvline(-math.log10(0.05), ls="--", c="grey", lw=1)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


def read_list(name: str) -> list[str]:
    path = TABLES / name
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    for old in [*OUT.glob("*__*.csv"), *OUT.glob("*__*.png")]:
        old.unlink()  # avoid stale plots from earlier runs

    tested = pd.read_csv(TABLES / "DE_results_all.csv")["symbol"].dropna()
    background = sorted(set(tested))
    gene_sets = {
        "UP": read_list("genes_UP_symbols.txt"),
        "DOWN": read_list("genes_DOWN_symbols.txt"),
        "SIG": read_list("genes_SIG_symbols.txt"),
    }

    try:
        background_id = add_background(background)
        print(f"Background: {len(background)} tested genes")
        summary = []
        for set_name, genes in gene_sets.items():
            if len(genes) < 5:
                print(f"[skip] {set_name}: {len(genes)} genes (< 5)")
                continue
            print(f"[{set_name}] {len(genes)} genes -> Enrichr")
            list_id = add_list(genes, f"{set_name}_Crohn_vs_Control")
            for library in LIBRARIES:
                table = enrich(list_id, background_id, library)
                time.sleep(0.5)  # be gentle with the public API
                if table.empty:
                    print(f"   {library}: no terms")
                    continue
                table.to_csv(OUT / f"{set_name}__{library}.csv", index=False)
                significant = table[table["adj_pvalue"] < 0.05]
                print(f"   {library}: {len(significant)} terms padj<0.05 (of {len(table)})")
                summary.append([set_name, library, len(significant), len(table)])
                if not significant.empty:
                    barplot(significant, f"{set_name} — {library}", OUT / f"{set_name}__{library}.png")
    except requests.RequestException as exc:
        print(f"ERROR: Enrichr request failed: {exc}", file=sys.stderr)
        return 1

    pd.DataFrame(summary, columns=["gene_set", "library", "n_sig_padj0.05", "n_total"]).to_csv(
        OUT / "enrichment_summary.csv", index=False
    )
    print("ENRICHMENT DONE ->", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
