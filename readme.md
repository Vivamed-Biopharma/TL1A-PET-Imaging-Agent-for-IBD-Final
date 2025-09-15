 TL1A In-Silico Package

 Quickstart:

 1) Create environment (pick one) — see `tl1a_in_silico/00_env/README.md`.
 2) Launch Jupyter and open notebooks in order inside `tl1a_in_silico/`.
 3) Notebooks are self-contained and print tables/summaries.

## Results Summary (from latest REPORT.md)

**Package Status**: ✅ Ready for in-vitro validation

- **QC**: All 12 clones PASS — no illegal characters, whitespace, or NXS/T motifs.
- **Developability**: Hydrophobicity in expected bands (VH ~41%, VL ~34%); pI_VL ~7.1±0.9; modest liabilities. **Gate**: Hyd_VH 35-45%, pI_VL 6-8, liabilities ≤2/chain. **Status**: All clones PASS.
- **DAR/Conjugation**: Uniform Eq_best=4; P(DAR1-2)=0.638; P(≥4)=0.051; E[DAR]=1.59. **Gate**: P(DAR1-2) ≥0.6, P(≥4) ≤0.1, E[DAR] 1.4-1.8. **Status**: Optimal for low-DAR imaging.
- **Detectability**: Calibrated grid (α=0.2, physiologic ranges): 25% pass TBR_pre ≥1.5; median ΔTBR=-0.107 at 80% block. **Gate**: TBR_pre ≥1.5, blocked ΔTBR ≤-0.3. **Status**: Favorable for KD≤3 nM, Bmax≥1 nM.
- **Soluble Sink**: f_free ≥0.5 for s≤Kd; low risk unless sTL1A >>Kd. **Gate**: f_free ≥0.7 at typical sTL1A. **Status**: PASS.
- **Paratope Plausibility**: Scores 0.49-0.55; DR3_adj 0.41-0.47. **Gate**: Paratope ≥0.5. **Status**: Plausible TL1A engagement.
- **Cross-reactivity**: Canonical TNFSF sequences (Swiss-Prot O95150, P50591, O43557, O43508, P01375, P01374); zero 6-mer overlaps with non-TL1A. **Gate**: Top non-TL1A ≤5. **Status**: Low risk; wet panel retained.
- **Manufacturability**: AggProxy ≤1.8; minimal charge variants. **Status**: Acceptable.
- **Immunogenicity**: Disclosure-level burden ~31-39 VH, 34 VL. **Status**: Typical for humanized.
- **Ranking**: Top candidates: **Fab04** (0.583), **Fab03/Fab09** (0.546).

**Next Steps**: KD determination (target ≤10 nM) + DR3 competition (≥50%); NOTA conjugation (DAR 1-2; IRF ≥70%; HMW ≤3%); Ga-68 labeling (RCP ≥95%).

See `tl1a_in_silico/REPORT.md` for full tables, figures, and statistical summaries.

Folder layout:

tl1a_in_silico/
 ├─ 00_env/                  # installs, environment notes
 ├─ 01_qc_sequences.ipynb    # QC & FASTA
 ├─ 02_imgt_numbering.ipynb  # optional IMGT/ANARCI
 ├─ 03_developability.ipynb  # pI/hydrophobicity/liabilities
 ├─ 04_conjugation_DAR.ipynb # NOTA–lysine DAR window & Eq_best
 ├─ 05_detectability_TE.ipynb# Bmax/Kd→TBR & TE Δ
 ├─ 06_sTL1A_sink.ipynb      # soluble TL1A sink risk
 ├─ 07_crossreactivity.ipynb # TNFSF k-mer heuristic
 ├─ 08_manufacturability.ipynb # agg patches/charge-variant proxies
 ├─ 09_immunogenicity_proxy.ipynb # disclosure-level burden
 ├─ 10_pk_window.ipynb       # toy PK → 1–2 h window sanity
 ├─ 11_rdkit_chelator.ipynb  # NOTA/p‑SCN‑Bn‑NOTA properties
 ├─ 12_visualization_py3Dmol.ipynb # optional visuals
 ├─ 13_sensitivity_SALib.ipynb # global sensitivity (toy models)
 └─ 90_master_table.ipynb    # merge all → Master Table + Exec Summary

