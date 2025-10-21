"""Shared helper utilities for TL1A in-silico analyses."""

from __future__ import annotations

import math
import os
import random
import re
from math import comb, exp
from typing import Dict, Iterable, List, Sequence, Tuple

AA_GROUPS = {
    "aromatic": list("FWY"),
    "polar": list("STNQ"),
    "acidic": list("DE"),
    "basic": list("KRH"),
    "aliphatic": list("AVILM"),
    "small": list("AGST"),
    "special": list("CP"),
}

AA_TO_GROUP = {aa: group for group, aas in AA_GROUPS.items() for aa in aas}


def avoid_glyco(seq: str) -> bool:
    """Return True when no NXS/T glycosylation motif exists."""
    return re.search(r"N[^P][ST]", seq) is None


def conservative_mutation(residue: str) -> str:
    """Pick a conservative substitution within the same physicochemical group."""
    group = AA_TO_GROUP.get(residue)
    if not group:
        return residue
    choices = [aa for aa in AA_GROUPS[group] if aa != residue]
    if not choices:
        return residue
    return random.choice(choices)


def mutate_cdr_sequence(seq: str, num_mutations: int, rng: random.Random | None = None) -> str:
    """Apply conservative mutations while screening out glycosylation motifs."""
    if num_mutations <= 0:
        return seq
    prng = rng or random
    idxs = list(range(len(seq)))
    prng.shuffle(idxs)
    seq_list = list(seq)
    changed = 0
    for idx in idxs:
        original = seq_list[idx]
        new = conservative_mutation(original)
        if new == original:
            continue
        seq_list[idx] = new
        candidate = "".join(seq_list)
        if avoid_glyco(candidate):
            changed += 1
        else:
            seq_list[idx] = original
        if changed >= num_mutations:
            break
    return "".join(seq_list)


def rebuild_chain(full_seq: str, cdr_map: Dict[str, str], new_cdrs: Dict[str, str], order: Sequence[str]) -> str:
    """Rebuild VH or VL sequences by substituting CDR segments in order."""
    updated = full_seq
    for key in order:
        original = cdr_map[key]
        replacement = new_cdrs[key]
        pos = updated.find(original)
        if pos >= 0:
            updated = updated[:pos] + replacement + updated[pos + len(original) :]
    return updated


def qc(seq: str) -> List[str]:
    """Return QC issues for a sequence (non-AA, whitespace, glyco motif)."""
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    issues: List[str] = []
    if any(ch not in valid for ch in seq):
        issues.append("non-AA")
    if re.search(r"\s", seq):
        issues.append("whitespace")
    if not avoid_glyco(seq):
        issues.append("NXS/T glyco motif")
    return issues


BJL = {"Cterm": 3.55, "Nterm": 7.50, "C": 9.0, "D": 4.05, "E": 4.45, "H": 5.98, "K": 10.0, "R": 12.0, "Y": 10.0}


def net_charge(seq: str, ph: float) -> float:
    """Compute net charge for a peptide at a given pH."""
    n_term = 1 / (1 + 10 ** (ph - BJL["Nterm"]))
    c_term = -1 / (1 + 10 ** (BJL["Cterm"] - ph))
    charge = n_term + c_term
    for aa in seq:
        if aa in "KRH":
            charge += 1 / (1 + 10 ** (ph - BJL[aa]))
        if aa in "DECY":
            charge -= 1 / (1 + 10 ** (BJL[aa] - ph))
    return charge


def calc_pI(seq: str, lo: float = 2.0, hi: float = 12.5) -> float:
    """Bisection search to approximate the isoelectric point."""
    a, b = lo, hi
    for _ in range(60):
        midpoint = (a + b) / 2
        fm = net_charge(seq, midpoint)
        if abs(fm) < 1e-4:
            return round(midpoint, 3)
        fa, fb = net_charge(seq, a), net_charge(seq, b)
        a, b = (midpoint, b) if fa * fm > 0 else (a, midpoint)
    return round((a + b) / 2, 3)


def hydrophobic_pct(seq: str) -> float:
    """Return hydrophobic residue percentage."""
    hydrophobic = sum(aa in "AFILMVWY" for aa in seq)
    return round(100 * hydrophobic / len(seq), 1) if seq else 0.0


def liabilities(seq: str) -> Dict[str, int]:
    """Count selected sequence liabilities."""
    return {"NG": seq.count("NG"), "DG": seq.count("DG"), "Met": seq.count("M"), "Trp": seq.count("W")}


def dar_stats(k_accessible: int, equivalents: int, eff: float = 0.45) -> Tuple[float, float, float]:
    """Return DAR probabilities and expectation for a given equivalent load."""
    paccept = 1 - exp(-eff * equivalents / max(1, k_accessible))

    def pmf(k: int) -> float:
        return comb(k_accessible, k) * (paccept**k) * ((1 - paccept) ** (k_accessible - k))

    p12 = sum(pmf(k) for k in (1, 2))
    p_ge4 = sum(pmf(k) for k in range(4, k_accessible + 1))
    expected = sum(k * pmf(k) for k in range(k_accessible + 1))
    return p12, p_ge4, expected


def TBR(Bmax_nM: float, Kd_nM: float, alpha: float = 0.4) -> float:
    """Simple binding potential model → Tissue-to-blood ratio."""
    binding_potential = Bmax_nM / max(1e-12, Kd_nM)
    return 1 + alpha * binding_potential


KD_SCALE = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}


def window_proxy(seq: str, window: int = 9) -> float:
    """Return the maximum hydropathy+charge aggregation proxy across a sliding window."""
    best = -1e9
    for idx in range(len(seq) - window + 1):
        win = seq[idx : idx + window]
        hyd = sum(KD_SCALE[aa] for aa in win) / window
        pos = sum(aa in "KRH" for aa in win) / window
        score = hyd + 0.5 * pos
        best = max(best, score)
    return round(best, 2)


def chem_motifs(seq: str) -> Dict[str, int]:
    """Return motif counts that correlate with charge heterogeneity."""
    return {motif: seq.count(motif) for motif in ("NS", "DS", "DP", "PR", "KK")}


def immu_proxy(seq: str, window: int = 15) -> int:
    """Count 15-mers enriched for aromatic anchors (rough immunogenicity proxy)."""
    anchors = "FYW"
    hits = 0
    for idx in range(len(seq) - window + 1):
        segment = seq[idx : idx + window]
        hits += sum(res in anchors for res in segment) >= 3
    return hits


def _enrich(segment: str) -> float:
    return sum(ch in "YSDNR" for ch in segment) / max(1, len(segment))


def paratope_score(H1: str, H2: str, H3: str, L1: str, L3: str) -> float:
    """Heuristic paratope enrichment score."""
    score = (
        0.5 * _enrich(H3)
        + 0.2 * _enrich(H2)
        + 0.15 * _enrich(L3)
        + 0.1 * _enrich(H1)
        + 0.05 * _enrich(L1)
    )
    return round(score, 3)


def dr3_adj(H3: str) -> float:
    """Approximate DR3 compatibility based on acidic enrichment and Lys penalty."""
    acid = sum(ch in "DEY" for ch in H3) / max(1, len(H3))
    penalty = 0.2 if "K" in H3 else 0.0
    return round(max(0.0, min(1.0, acid - penalty)), 3)


def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file if present."""
    if not os.path.exists(path):
        return {}
    seqs: Dict[str, str] = {}
    name: str | None = None
    buf: List[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def kmer_set(seq: str, k: int = 6) -> set[str]:
    """Return k-mer set for quick overlap checks."""
    return {seq[idx : idx + k] for idx in range(len(seq) - k + 1)}


def paratope_concat(heavy: Dict[str, str], light: Dict[str, str]) -> str:
    """Concatenate paratope features for cross-reactivity heuristics."""
    return heavy["H3"] + heavy["H2"] + light["L3"] + heavy["H1"] + light["L1"]


def free_fraction(Kd_nM: float, soluble_nM: float) -> float:
    """Soluble sink model: fraction of tracer remaining free."""
    return 1.0 - (soluble_nM / (Kd_nM + soluble_nM))
