 TL1A In-Silico Package

 Quickstart:

 1) Create environment (pick one) — see `tl1a_in_silico/00_env/README.md`.
 2) Launch Jupyter and open notebooks in order inside `tl1a_in_silico/`.
 3) Notebooks are self-contained and print tables/summaries.

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

