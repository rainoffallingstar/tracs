"""Generate the committed HOMER annotatePeaks.pl fixtures.

These are trimmed but structurally real: the same columns annotatePeaks.pl
writes by default, with the annotation suffixes it appends.

Coverage, chosen so the port's edge cases are all exercised:

  - `promoter-TSS`, `exon`, ... carry a `(NM_...)` suffix, which the R code
    strips by splitting on " (".
  - `Intergenic` and `NA` carry no suffix at all, so a parser that assumes a
    suffix is always present would corrupt them.
  - `3' UTR` / `5' UTR` are renamed to `3pUTR` / `5pUTR`.
  - `NA` is a literal annotation, not a missing value.
  - `ncRNA` is absent from the palette, so it is dropped from the plot while
    still counting toward `sum(N)`.
  - `non-coding` appears in H3K27ac but not H3K4me3, exercising the
    `dcast(fill = 0)` path where a kept row has no count for one sample.
  - The two samples have different peak totals, so `npeaks` and every fraction
    differ between them.
"""

HEADER = [
    "PeakID", "Chr", "Start", "End", "Strand", "Peak Score",
    "Focus Ratio/Region Size", "Annotation", "Detailed Annotation",
    "Distance to TSS", "Nearest PromoterID", "Entrez ID", "Nearest Unigene",
    "Nearest Refseq", "Nearest Ensembl", "Gene Name", "Gene Alias",
    "Gene Description", "Gene Type",
]

# Which annotations HOMER annotates with a nearest-feature suffix.
HAS_SUFFIX = {
    "promoter-TSS": True,
    "exon": True,
    "intron": True,
    "Intergenic": False,
    "TTS": True,
    "5' UTR": True,
    "3' UTR": True,
    "non-coding": True,
    "NA": False,
    "ncRNA": True,
}

GENE_TYPES = {
    "promoter-TSS": "protein-coding",
    "exon": "protein-coding",
    "intron": "protein-coding",
    "TTS": "protein-coding",
    "5' UTR": "protein-coding",
    "3' UTR": "protein-coding",
    "non-coding": "ncRNA",
    "ncRNA": "ncRNA",
    "Intergenic": "",
    "NA": "",
}

COUNTS = {
    "H3K27ac": {
        "promoter-TSS": 6, "exon": 4, "intron": 3, "Intergenic": 4,
        "TTS": 1, "5' UTR": 1, "3' UTR": 1, "non-coding": 1, "NA": 1,
        "ncRNA": 2,
    },
    "H3K4me3": {
        "promoter-TSS": 2, "exon": 5, "intron": 2, "Intergenic": 1,
        "TTS": 1, "5' UTR": 1, "3' UTR": 1, "non-coding": 0, "NA": 0,
        "ncRNA": 1,
    },
}


def write_sample(path, counts):
    rows = []
    peak = 0
    for category, count in counts.items():
        for _ in range(count):
            peak += 1
            needs_gene = HAS_SUFFIX[category]
            if category == "NA":
                annotation = "NA"
            elif needs_gene:
                annotation = f"{category} (NM_{peak:06d})"
            else:
                annotation = category
            start = 1000 + peak * 250
            rows.append([
                f"Peak{peak}",
                "chr1",
                str(start),
                str(start + 200),
                "+" if peak % 2 else "-",
                f"{10 + peak:.4f}",
                "0.000000",
                annotation,
                f"{category} (detailed)",
                str(-500 + peak * 10),
                f"NM_{peak:06d}" if needs_gene else "",
                str(1000 + peak) if needs_gene else "",
                "",
                "",
                f"ENSG{peak:011d}" if needs_gene else "",
                f"GENE{peak:03d}" if needs_gene else "",
                "",
                "test gene",
                GENE_TYPES[category],
            ])
    with open(path, "w") as handle:
        handle.write("\t".join(HEADER) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    return len(rows)


if __name__ == "__main__":
    import os
    import sys

    out_dir = sys.argv[1]
    os.makedirs(out_dir, exist_ok=True)
    for name, counts in COUNTS.items():
        path = os.path.join(out_dir, f"{name}.homer.txt")
        total = write_sample(path, counts)
        print(f"{path}: {total} peaks")
