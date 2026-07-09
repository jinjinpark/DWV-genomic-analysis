#!/usr/bin/env python3
"""
dnds_analysis.py

Pairwise dN/dS analysis (Nei-Gojobori 1986 method, Jukes-Cantor corrected)
for each Reconstructed Ancestral Sequence (RAS) lineage versus the most
recent common ancestor (MRCA), computed separately for each protein-coding
region of the DWV-A polyprotein.

The codon-aligned nucleotide sequences (one sequence per RAS lineage plus the
MRCA reference) are read directly from a FASTA file; no pooled read-depth data
is required. For every lineage-vs-MRCA pair, synonymous (S) and nonsynonymous
(N) sites and differences are counted per codon, proportions are corrected with
the Jukes-Cantor formula, and omega (dN/dS) is reported per protein region.

Usage
-----
    python3 dnds_analysis.py \
        --fasta RAS_nucleotide_sequence.fas \
        --reference Node167 \
        --outdir results

Inputs
------
--fasta       Codon-aligned nucleotide FASTA. All sequences must share the same
              length and reading frame (first position = codon position 1).
--reference   Header name of the MRCA reference sequence (default: Node167).
--outdir      Output directory (default: dnds_results).

Outputs
-------
<outdir>/dnds_per_lineage.csv   Per-lineage, per-protein dN, dS, dN/dS.
<outdir>/dnds_summary.csv       Per-protein mean dN, dS, dN/dS across lineages.

Requirements
------------
    Python 3.7+
    biopython

Method reference
----------------
Nei M, Gojobori T (1986). Simple methods for estimating the numbers of
synonymous and nonsynonymous nucleotide substitutions. Mol Biol Evol 3:418-426.
"""

import argparse
import csv
import math
import os
from itertools import product, permutations

from Bio.Seq import Seq


# --- DWV-A polyprotein protein-coding regions (nt coordinates, 1-based inclusive) ---
# Derived from confirmed amino-acid coordinates of the 2893-aa polyprotein.
PROTEIN_REGIONS = [
    ("Leader_protein", 1, 633),
    ("VP2", 634, 1392),
    ("VP4", 1393, 1455),
    ("VP3", 1456, 2703),
    ("VP1", 2704, 3477),
    ("2A-like_protein", 3478, 3861),
    ("Helicase", 3862, 5280),
    ("cleavage_site_Hel-3C", 5281, 6351),
    ("3C-like_protease", 6352, 7179),
    ("RdRp", 7180, 8679),
    ("Full_polyprotein", 1, 8679),
]

# Standard genetic code translation table for all 64 codons.
_BASES = "TCAG"
_CODON_TABLE = {
    "".join(c): str(Seq("".join(c)).translate())
    for c in product(_BASES, repeat=3)
}


def read_fasta(path):
    """Read a FASTA file into an ordered dict {name: uppercase_sequence}."""
    seqs = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            elif name is not None:
                seqs[name] += line.strip().upper()
    return seqs


def codon_sites(codon):
    """Return (synonymous_sites, nonsynonymous_sites) for a codon (Nei-Gojobori)."""
    if codon not in _CODON_TABLE or "N" in codon or "-" in codon:
        return None
    aa = _CODON_TABLE[codon]
    syn = 0.0
    for pos in range(3):
        for nt in "ACGT":
            if nt == codon[pos]:
                continue
            mutant = codon[:pos] + nt + codon[pos + 1:]
            if _CODON_TABLE[mutant] == aa:
                syn += 1.0 / 3.0
    return syn, 3.0 - syn


def count_diffs(c1, c2):
    """Count synonymous (Sd) and nonsynonymous (Nd) differences between two codons.

    When two or three positions differ, all mutational paths are averaged.
    Paths passing through a stop codon are discarded.
    """
    if c1 == c2:
        return 0.0, 0.0
    diff_pos = [i for i in range(3) if c1[i] != c2[i]]
    if not diff_pos:
        return 0.0, 0.0
    syn_total = nonsyn_total = 0.0
    n_paths = 0
    for perm in permutations(diff_pos):
        cur = c1
        s = n = 0
        valid = True
        for pos in perm:
            nxt = cur[:pos] + c2[pos] + cur[pos + 1:]
            if _CODON_TABLE.get(cur) == "*" or _CODON_TABLE.get(nxt) == "*":
                valid = False
                break
            if _CODON_TABLE[cur] == _CODON_TABLE[nxt]:
                s += 1
            else:
                n += 1
            cur = nxt
        if valid:
            syn_total += s
            nonsyn_total += n
            n_paths += 1
    if n_paths == 0:
        return 0.0, float(len(diff_pos))
    return syn_total / n_paths, nonsyn_total / n_paths


def jukes_cantor(p):
    """Jukes-Cantor correction for an observed proportion of differences."""
    if p is None or p >= 0.75:
        return None
    val = 1.0 - (4.0 / 3.0) * p
    if val <= 0:
        return None
    return -0.75 * math.log(val)


def pairwise_dnds(seq1, seq2, start, end):
    """Compute dN and dS between two codon-aligned sequences over [start, end]."""
    s_sites = n_sites = 0.0
    sd = nd = 0.0
    sub1 = seq1[start - 1:end]
    sub2 = seq2[start - 1:end]
    for i in range(0, len(sub1) - 2, 3):
        c1, c2 = sub1[i:i + 3], sub2[i:i + 3]
        r1, r2 = codon_sites(c1), codon_sites(c2)
        if r1 is None or r2 is None:
            continue
        s_sites += (r1[0] + r2[0]) / 2.0
        n_sites += (r1[1] + r2[1]) / 2.0
        d_s, d_n = count_diffs(c1, c2)
        sd += d_s
        nd += d_n
    if s_sites == 0 or n_sites == 0:
        return None
    p_s, p_n = sd / s_sites, nd / n_sites
    d_s, d_n = jukes_cantor(p_s), jukes_cantor(p_n)
    omega = (d_n / d_s) if (d_n is not None and d_s not in (None, 0)) else None
    return {
        "N_sites": n_sites, "S_sites": s_sites,
        "Nd": nd, "Sd": sd,
        "dN": d_n, "dS": d_s, "omega": omega,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Pairwise dN/dS (Nei-Gojobori) per protein for RAS lineages vs MRCA."
    )
    parser.add_argument("--fasta", required=True,
                        help="Codon-aligned nucleotide FASTA file.")
    parser.add_argument("--reference", default="Node167",
                        help="Header name of the MRCA reference (default: Node167).")
    parser.add_argument("--exclude", default="NC_004830.2",
                        help="Comma-separated sequence names to exclude from the "
                             "lineage set (default: NC_004830.2). Use '' to keep all.")
    parser.add_argument("--outdir", default="dnds_results",
                        help="Output directory (default: dnds_results).")
    args = parser.parse_args()

    seqs = read_fasta(args.fasta)
    if args.reference not in seqs:
        raise SystemExit(f"Reference '{args.reference}' not found in {args.fasta}.")

    excluded = {n.strip() for n in args.exclude.split(",") if n.strip()}
    ref_seq = seqs[args.reference]
    lineages = [n for n in seqs if n != args.reference and n not in excluded]

    os.makedirs(args.outdir, exist_ok=True)
    per_lineage_path = os.path.join(args.outdir, "dnds_per_lineage.csv")
    summary_path = os.path.join(args.outdir, "dnds_summary.csv")

    # Per-lineage results.
    with open(per_lineage_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["protein", "lineage", "N_sites", "S_sites",
                         "Nd", "Sd", "dN", "dS", "dN/dS"])
        summary = {}
        for pname, start, end in PROTEIN_REGIONS:
            summary[pname] = {"dN": [], "dS": [], "omega": []}
            for lin in lineages:
                r = pairwise_dnds(seqs[lin], ref_seq, start, end)
                if r is None:
                    writer.writerow([pname, lin, "", "", "", "", "", "", ""])
                    continue
                writer.writerow([
                    pname, lin,
                    f"{r['N_sites']:.2f}", f"{r['S_sites']:.2f}",
                    f"{r['Nd']:.2f}", f"{r['Sd']:.2f}",
                    "" if r["dN"] is None else f"{r['dN']:.4f}",
                    "" if r["dS"] is None else f"{r['dS']:.4f}",
                    "" if r["omega"] is None else f"{r['omega']:.4f}",
                ])
                if r["dN"] is not None:
                    summary[pname]["dN"].append(r["dN"])
                if r["dS"] is not None:
                    summary[pname]["dS"].append(r["dS"])
                if r["omega"] is not None:
                    summary[pname]["omega"].append(r["omega"])

    # Per-protein summary (mean across lineages).
    def mean(xs):
        return sum(xs) / len(xs) if xs else None

    with open(summary_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["protein", "mean_dN", "mean_dS",
                         "mean_dN/dS", "n_lineages_with_omega"])
        for pname, _, _ in PROTEIN_REGIONS:
            m_dn = mean(summary[pname]["dN"])
            m_ds = mean(summary[pname]["dS"])
            m_om = mean(summary[pname]["omega"])
            writer.writerow([
                pname,
                "" if m_dn is None else f"{m_dn:.4f}",
                "" if m_ds is None else f"{m_ds:.4f}",
                "" if m_om is None else f"{m_om:.4f}",
                len(summary[pname]["omega"]),
            ])

    print(f"Wrote {per_lineage_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
