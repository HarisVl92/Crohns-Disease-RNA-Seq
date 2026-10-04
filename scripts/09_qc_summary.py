#!/usr/bin/env python3
"""Per-sample QC summary of the whole pipeline.

Collects the read numbers from the Trim Galore JSON reports, the HISAT2 logs,
the featureCounts summary and the DESeq2 size factors into one table:
results/qc_summary.tsv. The same table is printed as Markdown for the README.
"""

import json
import re
from pathlib import Path

import pandas as pd

PROJ = Path(__file__).resolve().parents[1]


def trimming(run: str) -> dict:
    report = PROJ / f"data/trimmed/{run}.fastq.gz_trimming_report.json"
    reads = json.loads(report.read_text())["read_processing"]
    return {"raw_reads": reads["total_reads"], "trimmed_reads": reads["reads_written"]}


def alignment(run: str) -> dict:
    text = (PROJ / f"logs/align/{run}.hisat2.log").read_text()
    overall = re.search(r"([\d.]+)% overall alignment rate", text)
    multi = re.search(r"\(([\d.]+)%\) aligned >1 times", text)
    return {"aligned_pct": float(overall.group(1)), "multimapped_pct": float(multi.group(1))}


def assigned_reads() -> pd.Series:
    summary = pd.read_csv(PROJ / "results/counts/counts.txt.summary", sep="\t", index_col=0)
    summary.columns = [Path(c).name.removesuffix(".sorted.bam") for c in summary.columns]
    return summary.loc["Assigned"]


def to_markdown(table: pd.DataFrame) -> str:
    lines = ["| " + " | ".join(table.columns) + " |", "|" + " --- |" * len(table.columns)]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in table.itertuples(index=False)]
    return "\n".join(lines)


def main() -> None:
    manifest = pd.read_csv(PROJ / "data/samples_manifest.tsv", sep="\t", comment="#")
    assigned = assigned_reads()
    size_factors = pd.read_csv(PROJ / "results/DESeq2/tables/size_factors.csv", index_col="sample")

    rows = []
    for sample in manifest.itertuples(index=False):
        run = sample.run_accession
        reads = trimming(run) | alignment(run)
        rows.append(
            {
                "run": run,
                "condition": sample.condition,
                "sex": sample.sex,
                "raw_reads_M": round(reads["raw_reads"] / 1e6, 2),
                "after_trimming_pct": round(100 * reads["trimmed_reads"] / reads["raw_reads"], 1),
                "aligned_pct": reads["aligned_pct"],
                "multimapped_pct": reads["multimapped_pct"],
                "assigned_M": round(assigned[run] / 1e6, 2),
                "assigned_pct_of_trimmed": round(100 * assigned[run] / reads["trimmed_reads"], 1),
                "size_factor": round(size_factors.loc[run, "sizeFactor"], 3),
            }
        )

    table = pd.DataFrame(rows)
    table.to_csv(PROJ / "results/qc_summary.tsv", sep="\t", index=False)
    print(to_markdown(table))


if __name__ == "__main__":
    main()
